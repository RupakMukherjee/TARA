1d codes -- all serial

Burger_Turbulence_AB.f95: 1D Burgers equation + PS scheme + time solver AB

Burger_Turbulence_RK4.f95: 1D Burgers equation + PS scheme + time solver RK4

inhom_wave_eq.f95: 1D inhomogeneous wave equation + PS scheme + time solver RK4 + Second order derivative in time

weno_advection.f95: 1D advection problem using WENO scheme + time solver TVD-RK

weno_burgers.f95: 1D Burgers problem using WENO scheme + time solver TVD-RK

weno_advection.cpp: same as weno_advection.f95 but in C++ 

movie.scr.py: python script to generate gnu commands for plotting data from  weno_advection.f95 and weno_burgers.f95

movie_plot_c++.py: python script for plotting data from weno_advection.cpp

movie_plot_fortran.py: python script for plotting data from weno_advection.f95
