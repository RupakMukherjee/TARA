## CPU1D (Serial CPU Code)
- 1D Solver for Burgers' Equation. Burgers' equation or Bateman–Burgers equation is a fundamental partial differential equation occurring in various areas of applied mathematics, such as fluid mechanics, nonlinear acoustics, gas dynamics, and traffic flow.
```markdown
gfortran -I/usr/include -L/usr/include/lib Burger_Turbulence_AB.f95 -lfftw3 -lm
```
- Burger_Turbulence_AB.f95: Using Adam Bashforth Solver
- Burgers_Turbulence_RK4.f95: Using RK-4 Solver
