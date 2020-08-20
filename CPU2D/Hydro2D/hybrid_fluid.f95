PROGRAM fluid_fdps

  IMPLICIT NONE

  ! Define Grid Size.
integer ( kind = 4 ), parameter :: Nx = 64
integer ( kind = 4 ), parameter :: Ny = 64

real ( kind = 8 ), parameter :: pi = 3.14159265358979323846d0

integer ( kind = 4 ) i, j, t
double precision :: Lx, Ly, dx, dy, time, time_min, time_max, dt
real ( kind = 8 ) Energy
real ( kind = 8 ) x(Nx), y(Ny), omega(Nx,Ny),  domega_dx(Nx,Ny), domega_dy(Nx,Ny)

common/comm/Lx, Ly, dx, dy, dt

open(unit=40,file='Energy.dat',status='unknown')

!===================== USER INPUTS ============================================		

! System Length.
Lx = 2.0d0*pi
Ly = 2.0d0*pi

dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)

! Total Runtime and Time-Stepping.
time_min = 0.00d0
time_max = 1.00d0 
dt = 0.0010d0    

! Grid Generation.
do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0+real(j-1)*dy
    omega(i,j) = dsin(x(i))*dcos(y(j)) 
  end do ! j
end do ! i

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do time = time_min,time_max,dt

t = nint(time/dt) - int(time_min/dt)

  call weno_reconstruction (Nx,Ny,omega,domega_dx)
  call spectral_derivative(Nx,Ny,pi,omega,domega_dy)

  call ssp_rk3(Nx,Ny,omega,domega_dx,domega_dy)

Energy = 0.0d0 

do i = 1, Nx
  do j = 1, Ny
    Energy = Energy + (omega(i,j)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny))
      if (mod(dfloat(t),100.0) == 0.0) then
      write(t+100,*) x(i), y(j), omega(i,j)
      endif
  end do ! j
end do ! i

if (t /= 0 .and. mod(dfloat(t),100.0) == 0.0) then
  close (t+100)
endif
      
write(40,*) time, Energy
  
enddo ! time

END PROGRAM fluid_fdps

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine weno_reconstruction (Nx, Ny, omega, domega_dx)
implicit none
integer ( kind = 4 ) i, j
integer ( kind = 4 ) Nx, Ny
double precision :: Lx, Ly, dx, dy, dt
double precision :: epslion, gamma_1R, gamma_2R, gamma_3R
double precision :: omega, omega_back, omegahR
double precision :: omegahR_1, omegahR_2, omegahR_3
double precision :: beta_1R, beta_2R, beta_3R
double precision :: w1R_bar, w2R_bar, w3R_bar
double precision :: wR_d, w1R, w2R, w3R, domega_dx
dimension omega(Nx,Ny), omega_back(-1:Nx+2,Ny), omegahR(0:Nx,Ny)
dimension omegahR_1(Nx,Ny), omegahR_2(Nx,Ny), omegahR_3(Nx,Ny)
dimension beta_1R(Nx,Ny), beta_2R(Nx,Ny), beta_3R(Nx,Ny)
dimension w1R_bar(Nx,Ny), w2R_bar(Nx,Ny), w3R_bar(Nx,Ny)
dimension wR_d(Nx,Ny), w1R(Nx,Ny), w2R(Nx,Ny), w3R(Nx,Ny), domega_dx(Nx,Ny)

common/comm/Lx, Ly, dx, dy, dt

epslion = 1.0d0/10.0d0**6
gamma_1R = 1.0d0/10.0d0
gamma_2R = 3.0d0/5.0d0
gamma_3R = 3.0d0/10.0d0

do i = 1, Nx
  do j = 1, Ny
    omega_back(i,j) = omega(i,j)
  enddo
enddo

do j = 1, Ny
omega_back(-1,j) = omega_back(Nx-2,j)
omega_back(0,j) = omega_back(Nx-1,j)
omega_back(Nx+1,j) = omega_back(2,j)
omega_back(Nx+2,j) = omega_back(3,j)
enddo

do i = 1, Nx
  do j = 1, Ny
    omegahR_1(i,j) = (1.0d0/3.0d0)*omega_back(i-2,j) - (7.0d0/6.0d0)*omega_back(i-1,j) + (11.0d0/6.0d0)*omega_back(i,j)
    omegahR_2(i,j) = (-1.0d0/6.0d0)*omega_back(i-1,j) + (5.0d0/6.0d0)*omega_back(i,j) + (1.0d0/3.0d0)*omega_back(i+1,j)
    omegahR_3(i,j) = (1.0d0/3.0d0)*omega_back(i,j) + (5.0d0/6.0d0)*omega_back(i+1,j) - (1.0d0/6.0d0)*omega_back(i+2,j)

    beta_1R(i,j) = (13.0d0/12.0d0)*(omega_back(i-2,j) - 2.0d0*omega_back(i-1,j) + omega_back(i,j))**2 &
                +(1.0d0/4.0d0)*(omega_back(i-2,j) - 4.0d0*omega_back(i-1,j) + 3.0d0*omega_back(i,j))**2
    beta_2R(i,j) = (13.0d0/12.0d0)*(omega_back(i-1,j) - 2.0d0*omega_back(i,j) + omega_back(i+1,j))**2 &
                +(1.0d0/4.0d0)*(omega_back(i-1,j) - omega_back(i+1,j) )**2
    beta_3R(i,j) = (13.0d0/12.0d0)*(omega_back(i,j) - 2.0d0*omega_back(i+1,j) + omega_back(i+2,j))**2 &
                +(1.0d0/4.0d0)*(3.0d0*omega_back(i,j) - 4.0d0*omega_back(i+1,j) + omega_back(i+2,j))**2

    w1R_bar(i,j) = (gamma_1R)/(epslion+beta_1R(i,j))**2
    w2R_bar(i,j) = (gamma_2R)/(epslion+beta_2R(i,j))**2
    w3R_bar(i,j) = (gamma_3R)/(epslion+beta_3R(i,j))**2

    wR_d(i,j) = w1R_bar(i,j) + w2R_bar(i,j) + w3R_bar(i,j)
      
    w1R(i,j) = w1R_bar(i,j)/wR_d(i,j)
    w2R(i,j) = w2R_bar(i,j)/wR_d(i,j)
    w3R(i,j) = w3R_bar(i,j)/wR_d(i,j)
    
    omegahR(i,j) = w1R(i,j)*omegahR_1(i,j) + w2R(i,j)*omegahR_2(i,j) + w3R(i,j)*omegahR_3(i,j)
  enddo
enddo

do j = 1, Ny
  omegahR(0,j) = omegahR(Nx-1,j)
enddo

do i = 1, Nx
  do j = 1,Ny
    domega_dx(i,j) = - (omegahR(i,j) - omegahR(i-1,j)) / dx
  enddo
enddo

end subroutine weno_reconstruction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine spectral_derivative(Nx, Ny, pi, omega, domega_dy)
implicit none
include "fftw3.f"
integer ( kind = 4 ) i, j, is, js
integer ( kind = 4 ) Nx, Ny
real ( kind = 8 ) Lx, Ly, dx, dy, dt, pi, ky
real ( kind = 8 ) omega(Nx,Ny), domega_dy(Nx,Ny)
complex ( kind = 8 ) omegak(Nx/2+1,Ny), i_ky_omegak(Nx/2+1,Ny)
integer ( kind = 8 ) plan_forward, plan_backward ! FFTW Variables.

common/comm/Lx, Ly, dx, dy, dt

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, omega, omegak, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

! Evaluation of Derivatives in Spectral Space.
do is = 1,Nx/2+1
  do js = 1,Ny/2
    ky = 2.0d0*pi*dfloat(js-1)/Ly
    i_ky_omegak(is,js) = (0.0d0,1.0d0)*ky*omegak(is,js)
  end do
  do js = Ny/2+1,Ny
    ky = 2.0d0*pi*dfloat((js-1)-Ny)/Ly
    i_ky_omegak(is,js) = (0.0d0,1.0d0)*ky*omegak(is,js)
  end do ! j
end do ! i

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_omegak, domega_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

do i = 1,Nx
  do j = 1,Ny
    domega_dy(i,j) = domega_dy(i,j)/(dfloat(Nx)*dfloat(Ny))
  end do ! j
end do ! i

end subroutine spectral_derivative

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ssp_rk3(Nx,Ny,omega,domega_dx,domega_dy) 
implicit none
integer ( kind = 4 ) i, j
integer ( kind = 4 ) Nx, Ny
real ( kind = 8 ) Lx, Ly, dx, dy, dt, pi
real ( kind = 8 ) omega(Nx,Ny), domega_dy(Nx,Ny)
double precision :: omega_dum, omega_1, omega_2, domega_dx, RHS
dimension omega_dum(Nx,Ny), omega_1(Nx,Ny), omega_2(Nx,Ny), domega_dx(Nx,Ny), RHS(Nx,Ny)

common/comm/Lx, Ly, dx, dy, dt

  do i = 1, Nx
    do j = 1,Ny
    RHS = domega_dx(i,j) + domega_dy(i,j)
    enddo
  enddo

  do i = 1, Nx
    do j = 1, Ny
    omega_dum(i,j) = omega(i,j)
    omega_1(i,j) = omega(i,j) + dt * RHS(i,j)
    enddo
  enddo

  call weno_reconstruction (Nx,Ny,omega_1,domega_dx)
  call spectral_derivative(Nx,Ny,pi,omega_1,domega_dy)

  do i = 1, Nx
    do j = 1,Ny
    RHS = domega_dx(i,j) + domega_dy(i,j)
    enddo
  enddo

  do i = 1, Nx
    do j = 1, Ny
    omega_2(i,j) = 3.0d0*omega_dum(i,j)/4.0d0 + omega_1(i,j)/4.0d0 + dt * RHS(i,j)/4.0d0
    enddo
  enddo

  call weno_reconstruction (Nx,Ny,omega_2,domega_dx)
  call spectral_derivative(Nx,Ny,pi,omega_2,domega_dy)

  do i = 1, Nx
    do j = 1,Ny
    RHS = domega_dx(i,j) + domega_dy(i,j)
    enddo
  enddo

  do i = 1, Nx
    do j = 1, Ny
    omega(i,j) = omega_dum(i,j)/3.0d0 + 2.0d0*omega_2(i,j)/3.0d0 + dt * 2.0d0 * RHS(i,j)/3.0d0
    enddo
  enddo

end subroutine ssp_rk3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

