# VFEM2D
Verified Finite Element Method (2D)

VFEM2D aims to provide rigorous matrix assembly for FEM computations.
Currently, this library provides matrix computations for the following FEMs:
- Lagrange FEM of arbitrary order
- Raviart–Thomas FEM of arbitrary order

It also provides algorithms to obtain rigorous eigenvalue bounds for the Laplace operator and the Steklov operator. Liu's projection-based eigenvalue estimation algorithm [3,4], together with linear conforming FEM, is used to obtain guaranteed lower bounds for eigenvalues. The Lehmann–Goerisch method (Chap. 5 of [1]) is used to obtain high-precision eigenvalue bounds with higher-degree Lagrange finite elements.

Revision history:
- This MATLAB library was originally developed by Chun'guang You and Xuefeng Liu in 2016 for rigorous eigenvalue estimation of Steklov problems [2].
- 2024: Xuefeng Liu — revision for Dirichlet eigenvalue computation.
- 2025-11-21: Xuefeng Liu — sorted version for release.

References
1. Xuefeng Liu, Guaranteed Computational Methods for Self-Adjoint Differential Eigenvalue Problems, 2024, SpringerBriefs in Mathematics, Springer Singapore.
2. Chun'guang You, Hehu Xie, and Xuefeng Liu, Guaranteed Eigenvalue Bounds for the Steklov Eigenvalue Problem, SIAM Journal on Numerical Analysis, 57(3), 1395–1410, 2019. DOI:10.1137/18M1189592.
3. Xuefeng Liu, A framework of verified eigenvalue bounds for self-adjoint differential operators, Applied Mathematics and Computation, 267, pp. 341–355, 2015. DOI:10.1016/j.amc.2015.03.048.
4. Xuefeng Liu and Shin'ichi Oishi, Verified eigenvalue evaluation for the Laplacian over polygonal domains of arbitrary shape, SIAM Journal on Numerical Analysis, 51(3), 1634–1654, 2013. DOI:10.1137/120878446.