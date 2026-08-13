# Numerical model

The code solves the steady radiative transfer equation for a gray,
absorbing-emitting, non-scattering medium in axisymmetric `(r, z)` coordinates.

## Governing equation

Along a ray, the radiative transfer equation is

$$
\frac{dI}{ds}=\kappa(I_b-I), \qquad I_b=\frac{\sigma T^4}{\pi},
$$

where $I$ is radiative intensity, $\kappa$ is the absorption coefficient, and
$I_b$ is blackbody intensity. In axisymmetric coordinates, the conservative
form includes angular redistribution:

$$
\frac{\mu}{r}\frac{\partial(rI)}{\partial r}
-\frac{1}{r}\frac{\partial(\eta I)}{\partial\psi}
+\xi\frac{\partial I}{\partial z}
=\kappa(I_b-I).
$$

Here, $\mu$ and $\xi$ are the radial and axial direction cosines. The extra
angular term accounts for the rotation of the cylindrical basis.

## Angular discretization

`S10.txt` stores each ordinate as

```text
level  direction  mu  xi  weight
```

Directions with the same axial cosine form one level. The redistribution
coefficient is constructed recursively:

$$
\alpha_{m,0}=0, \qquad
\alpha_{m,l+1}=\alpha_{m,l}+w_{m,l}\mu_{m,l}.
$$

`loadQuadrature()` reads the ordinates and evaluates this recursion.
`solveAngularBoundary()` computes the starting angular intensity for every
level. `solveDirection()` then advances through that level using diamond
differencing in angle.

## Spatial discretization

The domain is divided into finite volumes. For each ordinate, `setSweep()`
chooses increasing or decreasing radial and axial traversal from the signs of
$\mu$ and $\xi$. `radialUpwind()` and `axialUpwind()` provide the incoming
intensity. This is the STEP, or first-order upwind, scheme.

The boundary conditions are mirror symmetry at $r=0$ and cold black boundaries
at the outer radius and both axial ends.

## Radiative source term

After all directions are swept, the incident radiation and volumetric source
term are

$$
G=\sum_{m,l}w_{m,l}I_{m,l}, \qquad
Q_r=\kappa(G-4\pi I_b).
$$

Negative $Q_r$ denotes net radiative cooling. `radiativeHeatLoss()` performs
the angular sweeps and evaluates $Q_r$; `writeResults()` writes cell-center
coordinates and $Q_r$ to the output file.

## Input and output

The input file contains:

1. number of radial and axial cells;
2. radial edge coordinates;
3. axial edge coordinates;
4. one `(temperature, kappa)` pair per cell, ordered by axial row and then
   radial column.

The output is a Tecplot-style text file with `r`, `z`, and `Qr` columns.

## Reference

M. Fathi, R. Hosseini, and M. Rahmani,
“Calculations of non-gray gas radiative heat transfer by coupling the discrete
ordinates method with the Leckner model in 3D rectangular enclosures,”
*Heat and Mass Transfer* (2016).
[doi:10.1007/s00231-015-1748-3](https://doi.org/10.1007/s00231-015-1748-3)

The cited work presents a 3D Cartesian, non-gray formulation. This repository
uses the same DOM and STEP foundations for a gray, axisymmetric example.
