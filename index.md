##TARA

Multi dimensional Navier-Stokes equation solver using pseudo-spectral and WENO algorithms for spatial discretization and explicit time iterators.


Contributors
------------

Developers:

- [Rupak Mukherjee](mailto:rupakmukherjee06@gmail.com): architecture and data structures, pushers, overall maintenance
- [Sayan Adhikari](mailto:sayanadhikari207@gmail.com): visualization toolkit and maintenance

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
