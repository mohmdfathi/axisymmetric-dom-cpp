#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

/*
  Simple axisymmetric discrete-ordinates radiation solver.

  The organization intentionally resembles the original Fortran program used 
  in my previous publication in Heat and Mass transfer journal.

  M. Fathi, R. Hosseini, and M. Rahmani,
  “Calculations of non-gray gas radiative heat transfer by coupling the discrete
   ordinates method with the Leckner model in 3D rectangular enclosures,”
   *Heat and Mass Transfer* (2016).

*/

constexpr int mmax = 10, lmax = 10;
int imax = 50, jmax = 50;
int ni = imax - 1, nj = jmax - 1;
constexpr double pi = 3.14159265358979323846;
constexpr double sigma = 5.670374419e-8;

std::vector<double> r, z;
double mu[mmax][lmax] = {}, xi[mmax][lmax] = {};
double weight[mmax][lmax] = {}, alpha[mmax][lmax + 1] = {};
int numberOfDirections[mmax] = {};

std::vector<std::vector<double>> kappa, Ib, Qr;
std::vector<std::vector<std::vector<double>>> Ih;
std::vector<std::vector<std::vector<std::vector<double>>>> Ip;

void resizeGrid(int cellsR, int cellsZ)
{
    ni = cellsR;
    nj = cellsZ;
    imax = ni + 1;
    jmax = nj + 1;

    r.assign(imax, 0.0);
    z.assign(jmax, 0.0);
    kappa.assign(ni, std::vector<double>(nj));
    Ib.assign(ni, std::vector<double>(nj));
    Qr.assign(ni, std::vector<double>(nj));
    Ih.assign(mmax, std::vector<std::vector<double>>(
        ni, std::vector<double>(nj)));
    Ip.assign(mmax, std::vector<std::vector<std::vector<double>>>(
        lmax, std::vector<std::vector<double>>(
            ni, std::vector<double>(nj))));
}


void initializeDefaultProblem()
{
    resizeGrid(49, 49);
    for (int i = 0; i < imax; ++i)
        r[i] = static_cast<double>(i) / static_cast<double>(ni);

    for (int j = 0; j < jmax; ++j)
        z[j] = static_cast<double>(j) / static_cast<double>(nj);
    for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
            kappa[i][j] = 0.1;
            double T2 = 1000.0 * 1000.0;
            Ib[i][j] = sigma * T2 * T2 / pi;
        }
    }
}

bool loadProblem(const char* filename)
{
    std::ifstream input(filename);
    if (!input)
        return false;

    int inputNi, inputNj;
    if (!(input >> inputNi >> inputNj) || inputNi < 1 || inputNj < 1) {
        std::cerr << "Error: invalid grid dimensions in " << filename << "\n";
        return false;
    }

    resizeGrid(inputNi, inputNj);

    for (double& value : r)
        if (!(input >> value)) return false;
    for (double& value : z)
        if (!(input >> value)) return false;

    for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
            double temperature;
            if (!(input >> temperature >> kappa[i][j]) ||
                temperature < 0.0 || kappa[i][j] < 0.0)
                return false;
            Ib[i][j] = sigma * std::pow(temperature, 4) / pi;
        }
    }
    return true;
}

bool loadQuadrature(const char* filename)
{
    std::ifstream input(filename);
    if (!input) {
        std::cerr << "Error: cannot open " << filename << "\n";
        return false;
    }

    int m, l;

    // S10.txt contains: m  l  mu  xi  doubled_weight
    while (input >> m >> l) {
        if (m < 0 || m >= mmax || l < 0 || l >= lmax ||
            !(input >> mu[m][l] >> xi[m][l] >> weight[m][l]) ||
            weight[m][l] <= 0.0) {
            std::cerr << "Error: invalid entry in " << filename << "\n";
            return false;
        }
        numberOfDirections[m] = l + 1;
    }

    for (int mm = 0; mm < mmax; ++mm) {
        for (int ll = 0; ll < numberOfDirections[mm]; ++ll)
            alpha[mm][ll + 1] = alpha[mm][ll] + weight[mm][ll] * mu[mm][ll];
    }
    return true;
}

void setSweep(double direction, int numberOfCells,
              int& first, int& last, int& step)
{
    if (direction > 0.0) {
        first = 0; last = numberOfCells; step = 1;
    } else {
        first = numberOfCells - 1; last = -1; step = -1;
    }
}

double radialUpwind(int m, int l, int i, int j)
{
    if (mu[m][l] > 0.0) {
        if (i > 0)
            return Ip[m][l][i - 1][j];

        // Reflection at the axis r=0.
        return Ip[m][numberOfDirections[m] - 1 - l][i][j];
    }

    if (i < ni - 1)
        return Ip[m][l][i + 1][j];

    return 0.0; // cold black outer wall
}

double axialUpwind(int m, int l, int i, int j)
{
    if (xi[m][l] > 0.0) {
        if (j > 0)
            return Ip[m][l][i][j - 1];

        return 0.0; // cold black bottom wall
    }

    if (j < nj - 1)
        return Ip[m][l][i][j + 1];

    return 0.0; // cold black top wall
}

// Compute the intensity at the first angular boundary, eta=0.
void solveAngularBoundary()
{
    for (int m = 0; m < mmax; ++m) {
        double xiH = xi[m][0];
        double muH = -std::sqrt(1.0 - xiH * xiH);
        int iFirst, iLast, iStep;
        int jFirst, jLast, jStep;

        setSweep(muH, ni, iFirst, iLast, iStep);
        setSweep(xiH, nj, jFirst, jLast, jStep);

        // j is the last array index, so keeping it in the inner loop gives
        // contiguous memory access in C++.
        for (int i = iFirst; i != iLast; i += iStep) {
            for (int j = jFirst; j != jLast; j += jStep) {
                // muHalf is negative, so the radial sweep moves inward.
                double rin = r[i + 1];
                double rout = r[i];
                double rave = 0.5 * (rin + rout);
                double dr = rout - rin;

                double zin  = xiH > 0.0 ? z[j]     : z[j + 1];
                double zout = xiH > 0.0 ? z[j + 1] : z[j];
                double dz = zout - zin;

                double Ir = 0.0;
                if (i < ni - 1)
                    Ir = Ih[m][i + 1][j];

                double Iz = 0.0;
                if (xiH > 0.0 && j > 0)
                    Iz = Ih[m][i][j - 1];
                if (xiH < 0.0 && j < nj - 1)
                    Iz = Ih[m][i][j + 1];

                double numerator =
                    kappa[i][j] * Ib[i][j]
                    + muH * Ir * rin / (rave * dr)
                    + xiH * Iz / dz;

                double denominator =
                    muH * rout / (rave * dr)
                    + xiH / dz
                    - muH / rave
                    + kappa[i][j];

                Ih[m][i][j] = numerator / denominator;
            }
        }
    }
}

// Sweep the spatial mesh for one discrete direction.
void solveDirection(int m, int l)
{
    int iFirst, iLast, iStep;
    int jFirst, jLast, jStep;

    setSweep(mu[m][l], ni, iFirst, iLast, iStep);
    setSweep(xi[m][l], nj, jFirst, jLast, jStep);

    // j is the last array index, so keeping it in the inner loop gives
    // contiguous memory access in C++.
    for (int i = iFirst; i != iLast; i += iStep) {
        for (int j = jFirst; j != jLast; j += jStep) {
            double rin  = mu[m][l] > 0.0 ? r[i]     : r[i + 1];
            double rout = mu[m][l] > 0.0 ? r[i + 1] : r[i];
            double zin  = xi[m][l] > 0.0 ? z[j]     : z[j + 1];
            double zout = xi[m][l] > 0.0 ? z[j + 1] : z[j];

            double rave = 0.5 * (rin + rout);
            double dr = rout - rin;
            double dz = zout - zin;

            double Ir = radialUpwind(m, l, i, j);
            double Iz = axialUpwind(m, l, i, j);
            double angularIn = Ih[m][i][j];

            double alphaIn = alpha[m][l];
            double alphaOut = alpha[m][l + 1];

            double numerator =
                kappa[i][j] * Ib[i][j]
                + mu[m][l] * Ir * rin / (rave * dr)
                - (alphaIn + alphaOut) * angularIn
                    / (rave * weight[m][l])
                + xi[m][l] * Iz / dz;

            double denominator =
                mu[m][l] * rout / (rave * dr)
                - 2.0 * alphaOut / (rave * weight[m][l])
                + xi[m][l] / dz
                + kappa[i][j];

            Ip[m][l][i][j] = numerator / denominator;

            // Diamond angular differencing:
            // I_l = (I_(l-1/2) + I_(l+1/2))/2.
            Ih[m][i][j] = 2.0 * Ip[m][l][i][j] - angularIn;
        }
    }
}

void radiativeHeatLoss()
{
    solveAngularBoundary();

    for (int m = 0; m < mmax; ++m)
        for (int l = 0; l < numberOfDirections[m]; ++l)
            solveDirection(m, l);

    for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
            double incidentRadiation = 0.0;

            for (int m = 0; m < mmax; ++m)
                for (int l = 0; l < numberOfDirections[m]; ++l)
                    incidentRadiation += weight[m][l] * Ip[m][l][i][j];

            // Negative Qr means radiative cooling of the material.
            Qr[i][j] =
                kappa[i][j] * (incidentRadiation - 4.0 * pi * Ib[i][j]);
        }
    }
}

bool writeResults(const char* filename)
{
    std::ofstream output(filename);

    if (!output) {
        std::cerr << "Error: " << filename << " could not be opened.\n";
        return false;
    }

    output << "VARIABLES=\"r\",\"z\",\"Qr\"\n";
    output << "ZONE I=" << ni << ", J=" << nj << ", F=POINT\n";
    output << std::scientific << std::setprecision(10);

    for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
            output
                << 0.5 * (r[i] + r[i + 1]) << " "
                << 0.5 * (z[j] + z[j + 1]) << " "
                << Qr[i][j] << "\n";
        }
    }
    return true;
}

int main(int argc, char* argv[])
{
    const char* problemFile = argc > 1 ? argv[1] : "problem_input.dat";
    const char* quadratureFile = argc > 2 ? argv[2] : "S10.txt";
    const char* outputFile = argc > 3 ? argv[3] : "QrFV.plt";

    if (loadProblem(problemFile))
        std::cout << "Loaded " << problemFile << "\n";
    else {
        if (argc > 1) {
            std::cerr << "Error: could not load " << problemFile << "\n";
            return 1;
        }
        std::cout << "Using the built-in uniform problem\n";
        initializeDefaultProblem();
    }
    if (!loadQuadrature(quadratureFile))
        return 1;

    radiativeHeatLoss();
    if (!writeResults(outputFile))
        return 1;

    double sumQr = 0.0;
    double minimumQr = Qr[0][0];
    double maximumQr = Qr[0][0];

    for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
            sumQr += Qr[i][j];

            if (Qr[i][j] < minimumQr)
                minimumQr = Qr[i][j];

            if (Qr[i][j] > maximumQr)
                maximumQr = Qr[i][j];
        }
    }

    std::cout << "Solved " << ni * nj << " cells\n";
    std::cout << "Qr sum   = " << sumQr << " W/m^3 (unweighted cell sum)\n";
    std::cout << "Qr range = [" << minimumQr << ", " << maximumQr << "] W/m^3\n";
    std::cout << "Wrote " << outputFile << "\n";

    return 0;
}
