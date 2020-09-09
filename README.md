Navier-Stokes
===============================================

Multi dimensional Navier-Stokes Equation solver using pseudo-spectral algorithm.


Contributors
------------

Developers:

- [Rupak Mukherjee](mailto:rupakmukherjee06@gmail.com): architecture and data structures, pushers, overall maintenance
- [Sayan Adhikari](mailto:sayanadhikari207@gmail.com): maintenance

Installation
------------
#### Common Prerequisites
1. gfortran compiler for fortran95 (for serial codes)
2. FFTW library (fftw3)
3. git

Directory Details
-----------------
## CPU1D (Serial CPU Code)
- 1D Solver for Burgers' Equation. Burgers' equation or Bateman–Burgers equation is a fundamental partial differential equation occurring in various areas of applied mathematics, such as fluid mechanics, nonlinear acoustics, gas dynamics, and traffic flow.
- Burger_Turbulence_AB.f95: Using Adam Bashforth Solver
- Burgers_Turbulence_RK4.f95: Using RK-4 Solver
## CPU2D (CPU Code)
### Hydro2D (Serial CPU Code)
A Serial Benchmarked Two Dimensional Incompressible Viscous Fluid code, using Pseudo-Spectral Method for Spatial Discritisation and Adams-Bashforth Technique for time evolution.
### MHD2D (Parallel CPU Code)
- compressible_mhd2d.f95: An OPENMP Parallel Benchmarked Compressible Viscous Neutral Fluid code, using Pseudo-Spectral Method with Multiple Time Solvers.
- mhd_parallel.f95: An OPENMP Parallel Benchmarked Two Dimensional Compressible Viscous Resistive MHD code, using Pseudo-Spectral Method with Multiple Time Solvers.
### Screened (Serial CPU Code)
- kukharkin.f95: (To be updated)
- kaladze2008_Eq_57and58.f95: (To be updated)
## CPU3D (Parallel CPU Code)
- mhd_3d.f95: An OPENMP Parallel Benchmarked Three Dimensional Compressible Viscous Resistive MHD code, using Pseudo-Spectral Method with Multiple Time Solvers.
## GPU3D (Parallel GPU Code)
### GMHD3D
#### large
- 512 x 512 x 512 grid simulation 
#### small
- 128 x 128 x 128 grid simulation 
