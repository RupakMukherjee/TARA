Program rupak
!
! This is a Serial Two Dimensional Two Fluid Incompressible Viscous code, 
! using Pseudo-Spectral Method for Spatial Discritisation and Adams-Bashforth Technique for time evolution.
!
! Wave-form /sim exp(i(kx-wt)). Hence derivative w.r.t. x gives ik only (and NOT -ik).
!
!___________________________________________________________________________________________________________________________________________
!
! Equations Solved in this code are in Vorticity (\omega) - Stream Function (\psi) Formalism.
!
!\begin{eqnarray}
!&& \frac{\partial \vec{\omega_1}}{\partial t} + \vec{u_1} \cdot \vec{\nabla} \vec{\omega_1} = \nu_1 \nabla^2 \vec{\omega_1} \\
!&& \frac{\partial \vec{\omega_2}}{\partial t} + \vec{u_2} \cdot \vec{\nabla} \vec{\omega_2} = \nu_2 \nabla^2 \vec{\omega_2} \\
!&& \vec{\omega_1} = \vec{\nabla} \times \vec{u_1} \\
!&& \vec{\omega_2} = \vec{\nabla} \times \vec{u_2} \\
!&& \nabla^2 \psi_1 = - \omega_1 
!&& \nabla^2 \psi_2 = - \omega_2 
!\end{eqnarray}
!
!___________________________________________________________________________________________________________________________________________
!
implicit none

! Define Grid Size.
integer ( kind = 4 ), parameter :: Nx = 64
integer ( kind = 4 ), parameter :: Ny = 64
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi = 3.14159265358979323846d0

include "fftw3.f"

integer ( kind = 4 ) i,j,t
real ( kind = 8 ) Lx,Ly,dx,dy,kx,ky,time,time_min,time_max,dt,Rex_1,Rey_1,nux_1,nuy_1,W0_1,m_1,d_1
real ( kind = 8 ) Rex_2,Rey_2,nux_2,nuy_2,W0_2,m_2,d_2
real ( kind = 8 ) Energy_1,y_Energy_1,Energy_2,y_Energy_2
real ( kind = 8 ) x(Nx),y(Ny),psi_1(Nx,Ny),omega_1(Nx,Ny),omega_1_dum(Nx,Ny),ux_1(Nx,Ny),uy_1(Nx,Ny),ux_1_dum(Nx,Ny),uy_1_dum(Nx,Ny)
real ( kind = 8 ) psi_2(Nx,Ny),omega_2(Nx,Ny),omega_2_dum(Nx,Ny),ux_2(Nx,Ny),uy_2(Nx,Ny),ux_2_dum(Nx,Ny),uy_2_dum(Nx,Ny)
complex ( kind = 8 ) omegak_1(Nh,Ny),omegak_1_dum(Nh,Ny),omegak_1_new(Nh,Ny),dt_omegak_1_old(Nh,Ny),dt_omegak_1_new(Nh,Ny)
complex ( kind = 8 ) omegak_2(Nh,Ny),omegak_2_dum(Nh,Ny),omegak_2_new(Nh,Ny),dt_omegak_2_old(Nh,Ny),dt_omegak_2_new(Nh,Ny)
complex ( kind = 8 ) psik_1(Nh,Ny),psik_1_dum(Nh,Ny),ukx_1(Nh,Ny),uky_1(Nh,Ny),ukx_1_dum(Nh,Ny),uky_1_dum(Nh,Ny),Ek_1(Nh,Ny)
complex ( kind = 8 ) psik_2(Nh,Ny),psik_2_dum(Nh,Ny),ukx_2(Nh,Ny),uky_2(Nh,Ny),ukx_2_dum(Nh,Ny),uky_2_dum(Nh,Ny),Ek_2(Nh,Ny)
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
Rex_1 = 10000.00
Rey_1 = 10000.00
Rex_2 = 10000.00
Rey_2 = 10000.00

nux_1 = 1.0d0/Rex_1
nuy_1 = 1.0d0/Rey_1
nux_2 = 1.0d0/Rex_2
nuy_2 = 1.0d0/Rey_2

! Initial Parameters.
W0_1 = 4.0d0
W0_2 = 0.50d0
m_1 = 3.0d0
m_2 = 1.0d0
d_1 = 1.0d0*pi/32.0d0
d_2 = 3.0d0*pi/64.0d0

! Grid Generation.
do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0+real(j-1)*dy
    ux_1(i,j) = 0.0d0
    uy_1(i,j) = 0.0d0
    ux_2(i,j) = 0.0d0
    uy_2(i,j) = 0.0d0
    ! Initial Vorticity Profile.
    omega_1(i,j) = W0_1/dcosh((y(j)+0.50d0*pi)/d_1)**2.0-W0_1/dcosh((y(j)-0.50d0*pi)/d_1)**2.0
    omega_1(i,j) = omega_1(i,j) - W0_1/dcosh((y(j)+1.50d0*pi)/d_1)**2.0+W0_1/dcosh((y(j)-1.50d0*pi)/d_1)**2.0
    omega_1(i,j) = omega_1(i,j)+0.01*dcos(m_1*x(i))
    omega_2(i,j) = dsin(m_2*x(i))*dcos(m_2*y(j))!W0_2/dcosh((y(j)+0.50d0*pi)/d_2)**2.0-W0_2/dcosh((y(j)-0.50d0*pi)/d_2)**2.0
    !omega_2(i,j) = !omega_2(i,j) - W0_2/dcosh((y(j)+1.50d0*pi)/d_2)**2.0+W0_2/dcosh((y(j)-1.50d0*pi)/d_2)**2.0
    !omega_2(i,j) = !omega_2(i,j)+0.01*dcos(m_2*x(i))
    ! Keep backup for FFTW.
    omega_1_dum(i,j) = omega_1(i,j)
    omega_2_dum(i,j) = omega_2(i,j)
    write(15,*) x(i),y(j),omega_1(i,j),ux_1(i,j),uy_1(i,j),omega_2(i,j),ux_2(i,j),uy_2(i,j)
  end do ! j
end do ! i

 close (15)
!===================== INITIAL TIME DATA ===============================================

! Move to Spectral Space.
  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, omega_1_dum, omegak_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, omega_2_dum, omegak_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

! Evaluation of Initial Strema-Function and Velocity Profile.
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik_1(i,j) = (0.0d0,0.0d0)
      ukx_1(i,j) = (0.0d0,0.0d0)
      uky_1(i,j) = (0.0d0,0.0d0)
      psik_2(i,j) = (0.0d0,0.0d0)
      ukx_2(i,j) = (0.0d0,0.0d0)
      uky_2(i,j) = (0.0d0,0.0d0)
      else
      psik_1(i,j) = omegak_1(i,j)/( kx*kx + ky*ky ) 
      ukx_1(i,j) = + (0.0d0,1.0d0) * ky * psik_1(i,j) 
      uky_1(i,j) = - (0.0d0,1.0d0) * kx * psik_1(i,j) 
      psik_2(i,j) = omegak_2(i,j)/( kx*kx + ky*ky ) 
      ukx_2(i,j) = + (0.0d0,1.0d0) * ky * psik_2(i,j) 
      uky_2(i,j) = - (0.0d0,1.0d0) * kx * psik_2(i,j) 
      endif
    psik_1_dum(i,j) = psik_1(i,j)  
    omegak_1_dum(i,j) = omegak_1(i,j)
    ukx_1_dum(i,j) = ukx_1(i,j)
    uky_1_dum(i,j) = uky_1(i,j)
    psik_2_dum(i,j) = psik_2(i,j)  
    omegak_2_dum(i,j) = omegak_2(i,j)
    ukx_2_dum(i,j) = ukx_2(i,j)
    uky_2_dum(i,j) = uky_2(i,j)  
  end do ! j
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    psik_1(i,j) = omegak_1(i,j)/( kx*kx + ky*ky ) 
    ukx_1(i,j) = + (0.0d0,1.0d0) * ky * psik_1(i,j) 
    uky_1(i,j) = - (0.0d0,1.0d0) * kx * psik_1(i,j) 
    psik_1_dum(i,j) = psik_1(i,j)  
    omegak_1_dum(i,j) = omegak_1(i,j)
    ukx_1_dum(i,j) = ukx_1(i,j)
    uky_1_dum(i,j) = uky_1(i,j)
    psik_2(i,j) = omegak_2(i,j)/( kx*kx + ky*ky ) 
    ukx_2(i,j) = + (0.0d0,1.0d0) * ky * psik_2(i,j) 
    uky_2(i,j) = - (0.0d0,1.0d0) * kx * psik_2(i,j) 
    psik_2_dum(i,j) = psik_2(i,j)  
    omegak_2_dum(i,j) = omegak_2(i,j)
    ukx_2_dum(i,j) = ukx_2(i,j)
    uky_2_dum(i,j) = uky_2(i,j)
  end do ! j
end do ! i

! Move to Real Space.
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, psik_1_dum, psi_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_1_dum, omega_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_1_dum, ux_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_1_dum, uy_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, psik_2_dum, psi_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_2_dum, omega_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_2_dum, ux_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_2_dum, uy_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

! FFTW Normalisations.
do i = 1, Nx
  do j = 1, Ny
    psi_1(i,j) = psi_1(i,j)/(dfloat(Nx)*dfloat(Ny))
    omega_1(i,j) = omega_1(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux_1(i,j) = ux_1(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy_1(i,j) = uy_1(i,j)/(dfloat(Nx)*dfloat(Ny))
    psi_2(i,j) = psi_2(i,j)/(dfloat(Nx)*dfloat(Ny))
    omega_2(i,j) = omega_2(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux_2(i,j) = ux_2(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy_2(i,j) = uy_2(i,j)/(dfloat(Nx)*dfloat(Ny))
    write(10,*) x(i),y(j),omega_1(i,j),psi_1(i,j),omega_2(i,j),psi_2(i,j)
    write(20,*) x(i),y(j),omega_1(i,j),ux_1(i,j),uy_1(i,j),omega_2(i,j),ux_2(i,j),uy_2(i,j)
  end do ! j
end do ! i

 close (10)
 close (20)
!======================= MAIN PROGRAM =====================================================

do time = time_min,time_max,dt

t = nint(time/dt) - int(time_min/dt)

!====================================================================================
  call derive (Nx,Ny,Nh,nux_1,nuy_1,nux_2,nuy_2,time,omegak_1,omegak_1_new,dt_omegak_1_old,dt_omegak_1_new, &
                      omegak_2,omegak_2_new,dt_omegak_2_old,dt_omegak_2_new) 
  call ab (Nx,Ny,Nh,nux_1,nuy_1,nux_2,nuy_2,time,dt,omegak_1,omegak_1_new,dt_omegak_1_old,dt_omegak_1_new, &
               omegak_2,omegak_2_new,dt_omegak_2_old,dt_omegak_2_new)
!====================================================================================

! Reset the Values.
do i = 1,Nx/2+1
  do j = 1,Ny
    dt_omegak_1_old(i,j) = dt_omegak_1_new(i,j)
    omegak_1(i,j) = omegak_1_new(i,j)
    dt_omegak_2_old(i,j) = dt_omegak_2_new(i,j)
    omegak_2(i,j) = omegak_2_new(i,j)
  end do ! j
end do ! i

! Evaluation of Stream-Function and Velocity Profile.
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik_1(i,j) = (0.0d0,0.0d0)
      ukx_1(i,j) = (0.0d0,0.0d0)
      uky_1(i,j) = (0.0d0,0.0d0)
      psik_2(i,j) = (0.0d0,0.0d0)
      ukx_2(i,j) = (0.0d0,0.0d0)
      uky_2(i,j) = (0.0d0,0.0d0)
      else
      psik_1(i,j) = omegak_1(i,j)/( kx*kx + ky*ky ) 
      ukx_1(i,j) = + (0.0d0,1.0d0) * ky * psik_1(i,j) 
      uky_1(i,j) = - (0.0d0,1.0d0) * kx * psik_1(i,j) 
      psik_2(i,j) = omegak_2(i,j)/( kx*kx + ky*ky ) 
      ukx_2(i,j) = + (0.0d0,1.0d0) * ky * psik_2(i,j) 
      uky_2(i,j) = - (0.0d0,1.0d0) * kx * psik_2(i,j) 
      endif
    omegak_1_dum(i,j) = omegak_1(i,j)
    ukx_1_dum(i,j) = ukx_1(i,j)
    uky_1_dum(i,j) = uky_1(i,j)
    omegak_2_dum(i,j) = omegak_2(i,j)
    ukx_2_dum(i,j) = ukx_2(i,j)
    uky_2_dum(i,j) = uky_2(i,j)
  end do ! j
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    psik_1(i,j) = omegak_1(i,j)/( kx*kx + ky*ky ) 
    ukx_1(i,j) = + (0.0d0,1.0d0) * ky * psik_1(i,j) 
    uky_1(i,j) = - (0.0d0,1.0d0) * kx * psik_1(i,j) 
    omegak_1_dum(i,j) = omegak_1(i,j)
    ukx_1_dum(i,j) = ukx_1(i,j)
    uky_1_dum(i,j) = uky_1(i,j)
    psik_2(i,j) = omegak_2(i,j)/( kx*kx + ky*ky ) 
    ukx_2(i,j) = + (0.0d0,1.0d0) * ky * psik_2(i,j) 
    uky_2(i,j) = - (0.0d0,1.0d0) * kx * psik_2(i,j) 
    omegak_2_dum(i,j) = omegak_2(i,j)
    ukx_2_dum(i,j) = ukx_2(i,j)
    uky_2_dum(i,j) = uky_2(i,j)
  end do ! j
end do ! i

! Move to Real Space.
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_1_dum, omega_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)
  
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_1_dum, ux_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_1_dum, uy_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)
  
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_2_dum, omega_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)
  
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_2_dum, ux_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_2_dum, uy_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

Energy_1 = 0.0d0
Energy_2 = 0.0d0
y_Energy_1 = 0.0d0
y_Energy_2 = 0.0d0

! FFTW Normalisation and Data Printing.   
do i = 1,Nx
  do j = 1,Ny
    omega_1(i,j) = omega_1(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux_1(i,j) = ux_1(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy_1(i,j) = uy_1(i,j)/(dfloat(Nx)*dfloat(Ny))
    Energy_1 = Energy_1 + (ux_1(i,j)**2 + uy_1(i,j)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny))
    y_Energy_1 = y_Energy_1 + (uy_1(i,j)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny))
    omega_2(i,j) = omega_2(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux_2(i,j) = ux_2(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy_2(i,j) = uy_2(i,j)/(dfloat(Nx)*dfloat(Ny))
    Energy_2 = Energy_2 + (ux_2(i,j)**2 + uy_2(i,j)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny))
    y_Energy_2 = y_Energy_2 + (uy_2(i,j)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny))
      if (mod(dfloat(t),100.0) == 0.0) then
      write(t+100,*) x(i),y(j),omega_1(i,j),ux_1(i,j),uy_1(i,j),omega_2(i,j),ux_2(i,j),uy_2(i,j)
      endif
  end do ! j
end do ! i

if (mod(dfloat(t),100.0) == 0.0) then
  close (t+100)
endif
      
write(40,*) time,Energy_1,y_Energy_1,Energy_2,y_Energy_2
  
end do ! time

contains

!===================================================================================
!================== SUBROUTINE NONLINEAR DERIVATIVE ================================
!===================================================================================

subroutine derive(Nx,Ny,Nh,nux_1,nuy_1,nux_2,nuy_2,time,omegak_1,omegak_1_new,dt_omegak_1_old,dt_omegak_1_new, &
                                omegak_2,omegak_2_new,dt_omegak_2_old,dt_omegak_2_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh
real ( kind = 8 ) time,nux_1,nuy_1,nux_2,nuy_2,kx,ky
real ( kind = 8 ) ux_1(Nx,Ny),uy_1(Nx,Ny),domega_1_dx(Nx,Ny),domega_1_dy(Nx,Ny)
real ( kind = 8 ) ux_1_domega_1_dx(Nx,Ny),uy_1_domega_1_dy(Nx,Ny)
complex ( kind = 8 ) ukx_1(Nh,Ny),ukx_1_dum(Nh,Ny),uky_1(Nh,Ny),uky_1_dum(Nh,Ny)
complex ( kind = 8 ) psik_1(Nh,Ny),omegak_1(Nh,Ny),i_kx_omegak_1(Nh,Ny),i_ky_omegak_1(Nh,Ny),omegak_1_new(Nh,Ny)
complex ( kind = 8 ) NLkx_1(Nh,Ny),NLky_1(Nh,Ny),NLk_1(Nh,Ny),dt_omegak_1_old(Nh,Ny),dt_omegak_1_new(Nh,Ny)
real ( kind = 8 ) ux_2(Nx,Ny),uy_2(Nx,Ny),domega_2_dx(Nx,Ny),domega_2_dy(Nx,Ny)
real ( kind = 8 ) ux_2_domega_2_dx(Nx,Ny),uy_2_domega_2_dy(Nx,Ny)
complex ( kind = 8 ) ukx_2(Nh,Ny),ukx_2_dum(Nh,Ny),uky_2(Nh,Ny),uky_2_dum(Nh,Ny)
complex ( kind = 8 ) psik_2(Nh,Ny),omegak_2(Nh,Ny),i_kx_omegak_2(Nh,Ny),i_ky_omegak_2(Nh,Ny),omegak_2_new(Nh,Ny)
complex ( kind = 8 ) NLkx_2(Nh,Ny),NLky_2(Nh,Ny),NLk_2(Nh,Ny),dt_omegak_2_old(Nh,Ny),dt_omegak_2_new(Nh,Ny)

! Keep backup for FFTW.
do i = 1,Nx/2+1
  do j = 1,Ny
    omegak_1_dum(i,j) = omegak_1(i,j) 
    omegak_2_dum(i,j) = omegak_2(i,j) 
  end do ! j
end do ! i
 
! Evaluation of Strema-Function and Velocity Profile from Vorticity Data.
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik_1(i,j) = (0.0d0,0.0d0)
      ukx_1(i,j) = (0.0d0,0.0d0)
      uky_1(i,j) = (0.0d0,0.0d0)
      psik_2(i,j) = (0.0d0,0.0d0)
      ukx_2(i,j) = (0.0d0,0.0d0)
      uky_2(i,j) = (0.0d0,0.0d0)
      else
      psik_1(i,j) = omegak_1_dum(i,j)/( kx*kx + ky*ky ) 
      ukx_1(i,j) = + (0.0d0,1.0d0) * ky * psik_1(i,j) 
      uky_1(i,j) = - (0.0d0,1.0d0) * kx * psik_1(i,j) 
      psik_2(i,j) = omegak_2_dum(i,j)/( kx*kx + ky*ky ) 
      ukx_2(i,j) = + (0.0d0,1.0d0) * ky * psik_2(i,j) 
      uky_2(i,j) = - (0.0d0,1.0d0) * kx * psik_2(i,j) 
      endif
    ukx_1_dum(i,j) = ukx_1(i,j)
    uky_1_dum(i,j) = uky_1(i,j)
    ukx_2_dum(i,j) = ukx_2(i,j)
    uky_2_dum(i,j) = uky_2(i,j)
  end do
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    psik_1(i,j) = omegak_1_dum(i,j)/( kx*kx + ky*ky ) 
    ukx_1(i,j) = + (0.0d0,1.0d0) * ky * psik_1(i,j) 
    uky_1(i,j) = - (0.0d0,1.0d0) * kx * psik_1(i,j) 
    ukx_1_dum(i,j) = ukx_1(i,j)
    uky_1_dum(i,j) = uky_1(i,j)
    psik_2(i,j) = omegak_2_dum(i,j)/( kx*kx + ky*ky ) 
    ukx_2(i,j) = + (0.0d0,1.0d0) * ky * psik_2(i,j) 
    uky_2(i,j) = - (0.0d0,1.0d0) * kx * psik_2(i,j) 
    ukx_2_dum(i,j) = ukx_2(i,j)
    uky_2_dum(i,j) = uky_2(i,j)
  end do ! j 
end do ! i

! Evaluation of Derivatives in Spectral Space.
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    i_kx_omegak_1(i,j) = (0.0d0,1.0d0)*kx*omegak_1_dum(i,j)
    i_ky_omegak_1(i,j) = (0.0d0,1.0d0)*ky*omegak_1_dum(i,j)
    i_kx_omegak_2(i,j) = (0.0d0,1.0d0)*kx*omegak_2_dum(i,j)
    i_ky_omegak_2(i,j) = (0.0d0,1.0d0)*ky*omegak_2_dum(i,j)
  end do
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    i_kx_omegak_1(i,j) = (0.0d0,1.0d0)*kx*omegak_1_dum(i,j)
    i_ky_omegak_1(i,j) = (0.0d0,1.0d0)*ky*omegak_1_dum(i,j)
    i_kx_omegak_2(i,j) = (0.0d0,1.0d0)*kx*omegak_2_dum(i,j)
    i_ky_omegak_2(i,j) = (0.0d0,1.0d0)*ky*omegak_2_dum(i,j)
  end do ! j
end do ! i

! Move to Real Space.
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_1_dum, ux_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_1_dum, uy_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_omegak_1, domega_1_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_omegak_1, domega_1_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_2_dum, ux_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_2_dum, uy_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_omegak_2, domega_2_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_omegak_2, domega_2_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

! FFTW Normalisations and Evaluation of Non-Linear Terms.  
do i = 1,Nx
  do j = 1,Ny
    ux_1(i,j) = ux_1(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy_1(i,j) = uy_1(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_1_dx(i,j) = domega_1_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_1_dy(i,j) = domega_1_dy(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux_1_domega_1_dx(i,j) = ux_1(i,j)*domega_1_dx(i,j)
    uy_1_domega_1_dy(i,j) = uy_1(i,j)*domega_1_dy(i,j)

    ux_2(i,j) = ux_2(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy_2(i,j) = uy_2(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_2_dx(i,j) = domega_2_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_2_dy(i,j) = domega_2_dy(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux_2_domega_2_dx(i,j) = ux_2(i,j)*domega_2_dx(i,j)
    uy_2_domega_2_dy(i,j) = uy_2(i,j)*domega_2_dy(i,j)
  end do ! j
end do ! i

! Move to Spectral Space.
  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_1_domega_1_dx, NLkx_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_1_domega_1_dy, NLky_1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_2_domega_2_dx, NLkx_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_2_domega_2_dy, NLky_2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)
  
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      ! De - Aliazing Technique - 2/3 Truncation...
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0) then    
      NLkx_1(i,j) = 0.0d0
      NLky_1(i,j) = 0.0d0
      NLkx_2(i,j) = 0.0d0
      NLky_2(i,j) = 0.0d0
      endif 
  end do ! j
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      ! De - Aliazing Technique - 2/3 Truncation...
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0) then
      NLkx_1(i,j) = 0.0d0
      NLky_1(i,j) = 0.0d0
      NLkx_2(i,j) = 0.0d0
      NLky_2(i,j) = 0.0d0
      endif 
  end do ! j
end do ! i 

! Evaluation of Nonlinear Term.  
do i = 1,Nx/2+1
  do j = 1,Ny/2  
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    NLk_1(i,j) = ( NLkx_1(i,j) + NLky_1(i,j) ) + nux_1*(omegak_1(i,j) - omegak_2(i,j))
    Nlk_1(i,j) = NLk_1(i,j) + (nux_1*kx*kx + nuy_1*ky*ky)*omegak_1(i,j)
    NLk_2(i,j) = ( NLkx_2(i,j) + NLky_2(i,j) ) + nux_1*(omegak_2(i,j) - omegak_1(i,j))
    Nlk_2(i,j) = NLk_2(i,j) + (nux_2*kx*kx + nuy_2*ky*ky)*omegak_2(i,j)
  end do ! j
  do j = Ny/2+1,Ny  
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    NLk_1(i,j) = ( NLkx_1(i,j) + NLky_1(i,j) ) + nux_1*(omegak_1(i,j) - omegak_2(i,j))
    Nlk_1(i,j) = NLk_1(i,j) + (nux_1*kx*kx + nuy_1*ky*ky)*omegak_1(i,j)
    NLk_2(i,j) = ( NLkx_2(i,j) + NLky_2(i,j) ) + nux_1*(omegak_2(i,j) - omegak_1(i,j))
    Nlk_2(i,j) = NLk_2(i,j) + (nux_2*kx*kx + nuy_2*ky*ky)*omegak_2(i,j)
  end do ! j
end do ! i

! Preparing for Time Evolution.
do i = 1,Nx/2+1
  do j = 1,Ny  
    dt_omegak_1_new(i,j) = - NLk_1(i,j)
    dt_omegak_2_new(i,j) = - NLk_2(i,j)
  end do ! j
end do ! i 

return

end subroutine derive

!===================================================================================
!=================== SUBROUTINE ADAMS BASHFORTH ====================================
!===================================================================================

subroutine ab(Nx,Ny,Nh,nux_1,nuy_1,nux_2,nuy_2,time,dt,omegak_1,omegak_1_new,dt_omegak_1_old,dt_omegak_1_new, &
                         omegak_2,omegak_2_new,dt_omegak_2_old,dt_omegak_2_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh
real ( kind = 8 ) time,dt,nux_1,nuy_1,nux_2,nuy_2
complex ( kind = 8 ) omegak_1(Nh,Ny),omegak_1_new(Nh,Ny)
complex ( kind = 8 ) dt_omegak_1_old(Nh,Ny)
complex ( kind = 8 ) dt_omegak_1_new(Nh,Ny)
complex ( kind = 8 ) omegak_2(Nh,Ny),omegak_2_new(Nh,Ny)
complex ( kind = 8 ) dt_omegak_2_old(Nh,Ny)
complex ( kind = 8 ) dt_omegak_2_new(Nh,Ny)

! Adams-Bashforth Method for Time Evolution.
do i = 1,Nh
  do j = 1,Ny
    omegak_1_new(i,j) = omegak_1(i,j) + ( (3.0d0/2.0d0)*dt_omegak_1_new(i,j) - (1.0d0/2.0d0)*dt_omegak_1_old(i,j) )*dt
    omegak_2_new(i,j) = omegak_2(i,j) + ( (3.0d0/2.0d0)*dt_omegak_2_new(i,j) - (1.0d0/2.0d0)*dt_omegak_2_old(i,j) )*dt
  end do ! j
end do ! i

return

end subroutine ab

!====================================================================================

end program rupak  


