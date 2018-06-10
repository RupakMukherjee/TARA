Program rupak
!
! This is a Serial Benchmarked Two Dimensional Incompressible Viscous Fluid code, 
! using Pseudo-Spectral Method for Spatial Discritisation and Adams-Bashforth Technique for time evolution.
!
! Wave-form /sim exp(i(kx-wt)). Hence derivative w.r.t. x gives ik only (and NOT -ik).
!
!___________________________________________________________________________________________________________________________________________
!
! Equations Solved in this code are in Vorticity (\omega) - Stream Function (\psi) Formalism.
!
!\begin{eqnarray}
!&& \frac{\partial \vec{\omega}}{\partial t} + \vec{u} \cdot \vec{\nabla} \vec{\omega} = \nu \nabla^2 \vec{\omega} \\
!&& \vec{\omega} = \vec{\nabla} \times \vec{u} \\
!&& \nabla^2 \psi = - \omega 
!\end{eqnarray}
!
!___________________________________________________________________________________________________________________________________________
!
implicit none

! Define Grid Size.
integer ( kind = 4 ), parameter :: Nx = 128
integer ( kind = 4 ), parameter :: Ny = 128
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi = 3.14159265358979323846d0

include "fftw3.f"

integer ( kind = 4 ) i,j,t
real ( kind = 8 ) Lx,Ly,dx,dy,kx,ky,time,time_min,time_max,dt,Rex,Rey,nux,nuy,W0,m,d
real ( kind = 8 ) Energy,y_Energy
real ( kind = 8 ) x(Nx),y(Ny),psi(Nx,Ny),omega(Nx,Ny),omega_dum(Nx,Ny),ux(Nx,Ny),uy(Nx,Ny),ux_dum(Nx,Ny),uy_dum(Nx,Ny)
complex ( kind = 8 ) omegak(Nh,Ny),omegak_dum(Nh,Ny),omegak_new(Nh,Ny),dt_omegak_old(Nh,Ny),dt_omegak_new(Nh,Ny)
complex ( kind = 8 ) psik(Nh,Ny),psik_dum(Nh,Ny),ukx(Nh,Ny),uky(Nh,Ny),ukx_dum(Nh,Ny),uky_dum(Nh,Ny),Ek(Nh,Ny)
integer ( kind = 8 ) plan_forward,plan_backward ! FFTW Variables.

open(unit=10,file='Initial_Condition_Omega.dat',status='unknown')
open(unit=15,file='Initial_Condition_Velocity.dat',status='unknown')
open(unit=20,file='Initial_Condition_Velocity_Reproduced.dat',status='unknown')
open(unit=40,file='Energy.dat',status='unknown')

!===================== USER INPUTS ============================================		

! System Length.
Lx = 2.0*pi
Ly = 2.0*pi

dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)

! Total Runtime and Time-Stepping.
time_min = 0.00d0
time_max = 50.0d0 
dt = 0.0010d0    

! Reynold's Number.
Rex = 10000.00
Rey = 10000.00

nux = 1.0d0/Rex
nuy = 1.0d0/Rey

! Initial Parameters.
W0 = 2.0d0
m = 3.0d0
d = 3.0d0*pi/128.0d0

! Grid Generation.
do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0+real(j-1)*dy
    ux(i,j) = 0.0d0
    uy(i,j) = 0.0d0
    ! Initial Vorticity Profile.
    omega(i,j) = W0/dcosh((y(j)+0.50d0*pi)/d)**2.0-W0/dcosh((y(j)-0.50d0*pi)/d)**2.0
    omega(i,j) = omega(i,j) - W0/dcosh((y(j)+1.50d0*pi)/d)**2.0+W0/dcosh((y(j)-1.50d0*pi)/d)**2.0
    omega(i,j) = omega(i,j)+0.01*dcos(m*x(i))
    ! Keep backup for FFTW.
    omega_dum(i,j) = omega(i,j)
    write(15,*) x(i),y(j),omega(i,j),ux(i,j),uy(i,j)
  end do ! j
end do ! i

 close (15)
!===================== INITIAL TIME DATA ===============================================

! Move to Spectral Space.
  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, omega_dum, omegak, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

! Evaluation of Initial Strema-Function and Velocity Profile.
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik(i,j) = (0.0d0,0.0d0)
      ukx(i,j) = (0.0d0,0.0d0)
      uky(i,j) = (0.0d0,0.0d0)
      else
      psik(i,j) = omegak(i,j)/( kx*kx + ky*ky ) 
      ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j) 
      uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j) 
      endif
    psik_dum(i,j) = psik(i,j)  
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  end do ! j
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    psik(i,j) = omegak(i,j)/( kx*kx + ky*ky ) 
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j) 
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j) 
    psik_dum(i,j) = psik(i,j)  
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  end do ! j
end do ! i

! Move to Real Space.
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, psik_dum, psi, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

! FFTW Normalisations.
do i = 1, Nx
  do j = 1, Ny
    psi(i,j) = psi(i,j)/(dfloat(Nx)*dfloat(Ny))
    omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
    write(10,*) x(i),y(j),omega(i,j),psi(i,j)
    write(20,*) x(i),y(j),omega(i,j),ux(i,j),uy(i,j)
  end do ! j
end do ! i

 close (10)
 close (20)
!======================= MAIN PROGRAM =====================================================

do time = time_min,time_max,dt

t = nint(time/dt) - int(time_min/dt)

!====================================================================================
  call derive (Nx,Ny,Nh,nux,nuy,time,omegak,omegak_new,dt_omegak_old,dt_omegak_new) 
  call ab (Nx,Ny,Nh,nux,nuy,time,dt,omegak,omegak_new,dt_omegak_old,dt_omegak_new)
!====================================================================================

! Reset the Values.
do i = 1,Nx/2+1
  do j = 1,Ny
    dt_omegak_old(i,j) = dt_omegak_new(i,j)
    omegak(i,j) = omegak_new(i,j)
  end do ! j
end do ! i

! Evaluation of Stream-Function and Velocity Profile.
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik(i,j) = (0.0d0,0.0d0)
      ukx(i,j) = (0.0d0,0.0d0)
      uky(i,j) = (0.0d0,0.0d0)
      else
      psik(i,j) = omegak(i,j)/( kx*kx + ky*ky ) 
      ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j) 
      uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j) 
      endif
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  end do ! j
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    psik(i,j) = omegak(i,j)/( kx*kx + ky*ky ) 
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j) 
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j) 
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  end do ! j
end do ! i

! Move to Real Space.
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)
  
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)
  
Energy = 0.0d0 
y_Energy = 0.0d0

! FFTW Normalisation and Data Printing.   
do i = 1,Nx
  do j = 1,Ny
    omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
    Energy = Energy + (ux(i,j)**2 + uy(i,j)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny))
    y_Energy = y_Energy + (uy(i,j)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny))
      if (mod(dfloat(t),100.0) == 0.0) then
      write(t+100,*) x(i),y(j),omega(i,j),ux(i,j),uy(i,j)
      endif
  end do ! j
end do ! i

if (mod(dfloat(t),100.0) == 0.0) then
  close (t+100)
endif
      
write(40,*) time,Energy,y_Energy
  
end do ! time

contains

!===================================================================================
!================== SUBROUTINE NONLINEAR DERIVATIVE ================================
!===================================================================================

subroutine derive(Nx,Ny,Nh,nux,nuy,time,omegak,omegak_new,dt_omegak_old,dt_omegak_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh
real ( kind = 8 ) time,nux,nuy,kx,ky
real ( kind = 8 ) ux(Nx,Ny),uy(Nx,Ny),domega_dx(Nx,Ny),domega_dy(Nx,Ny)
real ( kind = 8 ) ux_domega_dx(Nx,Ny),uy_domega_dy(Nx,Ny)
complex ( kind = 8 ) ukx(Nh,Ny),ukx_dum(Nh,Ny),uky(Nh,Ny),uky_dum(Nh,Ny)
complex ( kind = 8 ) psik(Nh,Ny),omegak(Nh,Ny),i_kx_omegak(Nh,Ny),i_ky_omegak(Nh,Ny),omegak_new(Nh,Ny)
complex ( kind = 8 ) NLkx(Nh,Ny),NLky(Nh,Ny),NLk(Nh,Ny),dt_omegak_old(Nh,Ny),dt_omegak_new(Nh,Ny)

! Keep backup for FFTW.
do i = 1,Nx/2+1
  do j = 1,Ny
    omegak_dum(i,j) = omegak(i,j) 
  end do ! j
end do ! i
 
! Evaluation of Strema-Function and Velocity Profile from Vorticity Data.
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik(i,j) = (0.0d0,0.0d0)
      ukx(i,j) = (0.0d0,0.0d0)
      uky(i,j) = (0.0d0,0.0d0)
      else
      psik(i,j) = omegak_dum(i,j)/( kx*kx + ky*ky ) 
      ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j) 
      uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j) 
      endif
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  end do
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    psik(i,j) = omegak_dum(i,j)/( kx*kx + ky*ky ) 
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j) 
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j) 
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  end do ! j 
end do ! i

! Evaluation of Derivatives in Spectral Space.
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    i_kx_omegak(i,j) = (0.0d0,1.0d0)*kx*omegak_dum(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0)*ky*omegak_dum(i,j)
  end do
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    i_kx_omegak(i,j) = (0.0d0,1.0d0)*kx*omegak_dum(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0)*ky*omegak_dum(i,j)
  end do ! j
end do ! i

! Move to Real Space.
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_omegak, domega_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_omegak, domega_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

! FFTW Normalisations and Evaluation of Non-Linear Terms.  
do i = 1,Nx
  do j = 1,Ny
    ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_dx(i,j) = domega_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_dy(i,j) = domega_dy(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux_domega_dx(i,j) = ux(i,j)*domega_dx(i,j)
    uy_domega_dy(i,j) = uy(i,j)*domega_dy(i,j)
  end do ! j
end do ! i

! Move to Spectral Space.
  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_domega_dx, NLkx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_domega_dy, NLky, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)
  
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      ! De - Aliazing Technique - 2/3 Truncation...
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0) then    
      NLkx(i,j) = 0.0d0
      NLky(i,j) = 0.0d0
      endif 
  end do ! j
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      ! De - Aliazing Technique - 2/3 Truncation...
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0) then
      NLkx(i,j) = 0.0d0
      NLky(i,j) = 0.0d0
      endif 
  end do ! j
end do ! i 

! Evaluation of Nonlinear Term.  
do i = 1,Nx/2+1
  do j = 1,Ny/2  
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    NLk(i,j) = ( NLkx(i,j) + NLky(i,j) )
    Nlk(i,j) = NLk(i,j) + (nux*kx*kx + nuy*ky*ky)*omegak(i,j)
  end do ! j
  do j = Ny/2+1,Ny  
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    NLk(i,j) = ( NLkx(i,j) + NLky(i,j) )
    Nlk(i,j) = NLk(i,j) + (nux*kx*kx + nuy*ky*ky)*omegak(i,j)
  end do ! j
end do ! i

! Preparing for Time Evolution.
do i = 1,Nx/2+1
  do j = 1,Ny  
    dt_omegak_new(i,j) = - NLk(i,j)
  end do ! j
end do ! i 

return

end subroutine derive

!===================================================================================
!=================== SUBROUTINE ADAMS BASHFORTH ====================================
!===================================================================================

subroutine ab(Nx,Ny,Nh,nux,nuy,time,dt,omegak,omegak_new,dt_omegak_old,dt_omegak_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh
real ( kind = 8 ) time,dt,nux,nuy
complex ( kind = 8 ) omegak(Nh,Ny),omegak_new(Nh,Ny)
complex ( kind = 8 ) dt_omegak_old(Nh,Ny)
complex ( kind = 8 ) dt_omegak_new(Nh,Ny)

! Adams-Bashforth Method for Time Evolution.
do i = 1,Nh
  do j = 1,Ny
    omegak_new(i,j) = omegak(i,j) + ( (3.0d0/2.0d0)*dt_omegak_new(i,j) - (1.0d0/2.0d0)*dt_omegak_old(i,j) )*dt
  end do ! j
end do ! i

return

end subroutine ab

!====================================================================================

end program rupak  

