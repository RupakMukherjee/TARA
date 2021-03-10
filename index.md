The TARA simulation architechture is a multi-dimensional pseudo-spectral solver for weakly compressible and incompressible magneto-hydrodynamic flows. TARA is flexible for adding higher order fluid-moment equations with minimal computational overhead. This framework runs efficiently on GPU architechture. In addition, the performance scales efficiently under MPI on massively parallel shared- or distributed-memory computers.

The TARA simulation framework have been used for many different applications in astrophysical studies as well as terrestrial laboratory plasma simulations.

# Contributors

- [Rupak Mukherjee](mailto:rupakmukherjee06@gmail.com): architecture and data structures, pushers, overall maintenance
- [Sayan Adhikari](mailto:sayanadhikari207@gmail.com): visualization toolkit and maintenance
- [Shishir Biswas](mailto:shishirbeafriend@gmail.com): benchmarking, visualization toolkit and maintenance

# Installation

## Common Prerequisites for CPU version
1. [GNU Compiler (higher than version 4.0.0)](https://gcc.gnu.org/)
2. [OpenMP](https://www.openmp.org/) and [MPI](https://www.open-mpi.org/) architechture
3. [FFTW library (higher than Version 3.3.3)](http://www.fftw.org/)
4. [git](https://git-scm.com/)

## Common Prerequisites for single-GPU version
1. [CUDA Toolkit (higher than CUDA 9.1)](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html)
2. [PGI Compilers & Tools (higher than 18.1)](https://www.pgroup.com/support/new_rel_80.htm) 
3. [cuFFT library](https://developer.nvidia.com/cufft)
4. [git](https://git-scm.com/)

## Common Prerequisites for multi-GPU version
1. [CUDA Toolkit (higher than CUDA 9.1)](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html)
2. [PGI Compilers & Tools (higher than 18.1)](https://www.pgroup.com/support/new_rel_80.htm) 
3. [AccFFT library](http://accfft.org/about/)
4. [PnetCDF library (higher than 1.11.0)](https://parallel-netcdf.github.io/)
5. [git](https://git-scm.com/)

# Directory Details

## [CPU1D (Serial CPU Code)](cpu1d.md)
- 1D Solver for Burgers' Equation. Burgers' equation or Bateman–Burgers equation is a fundamental partial differential equation occurring in various areas of applied mathematics, such as fluid mechanics, nonlinear acoustics, gas dynamics, and traffic flow.
```markdown
gfortran -I/usr/include -L/usr/include/lib Burger_Turbulence_AB.f95 -lfftw3 -lm
```
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

## GMHD3D

### large
- 512 x 512 x 512 grid simulation 

### small
- 128 x 128 x 128 grid simulation 



## Welcome to GitHub Pages

You can use the [editor on GitHub](https://github.com/RupakMukherjee/TARA/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/RupakMukherjee/TARA/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
