# User Manual

TARA user manual is coming soon! Stay tuned!
```
! This is an OPENMP Parallel Benchmarked Compressible Viscous Neutral Fluid code, 
! using Pseudo-Spectral Method with Multiple Time Solvers.
! Here Adams-Bashforth Subroutine is Parallalised. Other Time Solvers are NOT Parallalised.
!
! Wave-form /sim exp(i(kx-wt)). Hence derivative w.r.t. x gives ik only (and NOT -ik).
! 
! To Run this code in Desktop: gfortran -fopenmp -I/usr/local/include -L/usr/local/lib velocity_OMP.f95 -lfftw3 -lm
!
!------------------------------------------------------------------------------------------------------------------------
!
!                                         BASIC EQUATIONS 
!
!               \frac{\partial n}{\partial t} + \vec{\nabla} \cdot (n \vec{v}) = 0
!
!          \frac{\partial \vec{v}}{\partial t} + ( \vec{v} \cdot \vec{\nabla} ) \vec{v} 
!               = \frac{\mu}{n} {\nabla}^2 \vec{v} - \frac{C_s^2}{n} \vec{\nabla} n 
!            
!------------------------------------------------------------------------------------------------------------------------
!
!                                           PSEUDOCODE 
!
!                             Input-- n(x,y), ux(x,y), uy(x,y)  
!
!                                             FFTW
!             [At Time = t]     nk(kx,ky), ukx(kx,ky), uky(kx,ky)   
!
!                                           Evaluate
!                                       i*kx*nk, i*ky*nk    
!                            i*kx*ukx, i*ky*ukx, i*kx*uky, i*ky*uky
!                        (kx^2)*ukx, (ky^2)*ukx, (kx^2)*uky, (ky^2)*uky        
!                  
!                                            IFFTW
!                                          n, ux, uy
!                                        dn/dx, dn/dy
!                              dux/dx, dux/dy, duy/dx, duy/dy
!                    d^2(ux)/dx^2, d^2(ux)/dy^2, d^2(uy)/dx^2, d^2(uy)/dy^2  
!                                                    
!                                           Evaluate
!                              Px = (1/n)*dn/dx, Py = (1/n)*dn/dy    
!                   Viscxx = (1/n)*d^2(ux)/dx^2, Viscxy = (1/n)*d^2(ux)/dy^2
!                   Viscyx = (1/n)*d^2(uy)/dx^2, Viscyy = (1/n)*d^2(uy)/dy^2      
! 
!                                           Evaluate
!                       NLn = ux*(dn/dx) + uy*(dn/dy) + n*(dux/dx + duy/dy)                      -- Density Equation           
!                               NLx = ux*(dux/dx) + uy*(dux/dy)                 
!                               NLy = ux*(duy/dx) + uy*(duy/dy)                
!
!                                            FFTW
!                                          Pkx, Pky
!                             Visckxx, Visckxy, Visckyx, Visckyy        
!                                      NLkn, NLkx, NLky
!
!                                           Evaluate   
!                     NLkx = NLkx + mux*Visckxx + muy*Visckxy + (CS^2)*Pkx      -- X - Component of Momentum Equation
!                     NLky = NLky + mux*Visckyx + muy*Visckyy + (CS^2)*Pky      -- Y - Component of Momentum Equation              
! 
!                                           Evaluate
!                                      d(nk)/dt = - NLkn
!                                     d(ukx)/dt = - NLkx
!                                     d(uky)/dt = - NLky
!
!             Time Solver -- Adams-Bashforth/Euler/Predictor-Corrector/Runge-Kutta(4)
!                                       
!        [At Time = t + dt]     nk(kx,ky), ukx(kx,ky), uky(kx,ky)    
!
!------------------------------------------------------------------------------------------------------------------------ 
```
