!#############################################################################################################
!This is a fft based Serial Incompressible Hydrodynamic (2D) solver.
!Developer: Dr.Rupak Mukherjee (IPR Phd 2019) 
!GPU Accleration is given by : Mr. Shishir Biswas (IPR Phd)
!Benchmarking Problem: KH Instability by: PG Drazin. Discontinuous velocity profiles for the orr-sommerfeld equation. Journal of Fluid Mechanics, 10(4):571–583, 1961.
!Date: 13/07/2020
!#############################################################################################################

!########################################## Module for invoking cufft ########################################

module cufft2D
    integer, public :: CUFFT_FORWARD = -1
    integer, public :: CUFFT_INVERSE = 1
    integer, public :: CUFFT_R2C = Z'2a' ! Real to Complex (interleaved)
    integer, public :: CUFFT_C2R = Z'2c' ! Complex (interleaved) to Real
    integer, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved
    integer, public :: CUFFT_D2Z = Z'6a' ! Double to Double-Complex
    integer, public :: CUFFT_Z2D = Z'6c' ! Double-Complex to Double
    integer, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex

    interface
        subroutine createCUFFTPlan2D(plan, nx, ny, planType, stream) bind (C, name = 'createCUFFTPlan2D')
            use iso_c_binding
            use openacc
            implicit none
            type (c_ptr) :: plan
            integer (c_int), value :: nx, ny, planType
            integer (acc_handle_kind), value :: stream
        end subroutine createCUFFTPlan2D
    end interface

    interface
        subroutine executeCUFFT2D(plan, iData, oData, planType) bind (C, name = 'executeCUFFT2D')
            use iso_c_binding
            use openacc
            implicit none
            type (c_ptr) ::  plan
            type (c_ptr), value :: iData
            type (c_ptr), value :: oData
            integer (c_int), value :: planType
        end subroutine executeCUFFT2D
    end interface

    interface
        subroutine destroyCUFFTPlan2D(plan) bind (C, name = 'destroyCUFFTPlan2D')
            use iso_c_binding
            use openacc
            implicit none
            type (c_ptr) :: plan
        end subroutine destroyCUFFTPlan2D
    end interface

end module cufft2D
!####################################################################################################################


Program PSHD2


use omp_lib
use cufft2D
use openacc
implicit none

integer ( kind = 4 ), parameter :: Nx = 256
integer ( kind = 4 ), parameter :: Ny = 256
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

!include "fftw3.f"

integer ( kind = 4 ) i,j,t
real ( kind = 8 ) Lx,Ly,dx,dy,kx,ky,time,time_min,time_max,dt,nux,nuy,C0,C1,C2,C3,W0,m,d,G,Gb
real ( kind = 8 ) Energy,Energy0,y_Energy,y_Energy0
real ( kind = 8 ) x(Nx),y(Ny),psi(Nx,Ny),omega(Nx,Ny),omega_dum(Nx,Ny),ux(Nx,Ny),uy(Nx,Ny),ux_dum(Nx,Ny),uy_dum(Nx,Ny)
complex ( kind = 8 ) omegak(Nh,Ny),omegak_dum(Nh,Ny),omegak_new(Nh,Ny),dt_omegak_old(Nh,Ny),dt_omegak_new(Nh,Ny)
complex ( kind = 8 ) psik(Nh,Ny),psik_dum(Nh,Ny),ukx(Nh,Ny),uky(Nh,Ny),ukx_dum(Nh,Ny),uky_dum(Nh,Ny),Ek(Nh,Ny)
integer ( kind = 8 ) plan_forward,plan_backward

!cufft
integer (acc_handle_kind) :: stream
type (c_ptr) :: fftPlanD2ZMain, fftPlanZ2DMain

real ( kind = 8 ) t1,t2
call acc_init(acc_device_nvidia)

open(unit=5,file='System_information.dat',status='unknown')
open(unit=10,file='Initial_Condition_Omega.dat',status='unknown')
!open(unit=15,file='Initial_Condition_Velocity.dat',status='unknown')
open(unit=20,file='Initial_Condition_Velocity_REP.dat',status='unknown')
!open(unit=25,file='OMEGA',status='old')
!open(unit=35,file='Energy_Spectra.dat',status='unknown')
open(unit=40,file='Energy.dat',status='unknown')
!open(unit=50,file='OMEGA',status='old')

!===================== USER INPUTS ============================================		

Lx = 2.0*pi
Ly = 2.0*pi
dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)

time_min = 0.00d0
time_max = 200.0d0 
dt = 0.0010d0    

nux = 0.00010d0
nuy = 0.00010d0

W0 = 2.0d0
m = 3.0d0
d = 3.0d0*pi/128.0d0

do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0+real(j-1)*dy
    ux(i,j) = 0.0d0
    uy(i,j) = 0.0d0
    !read(25,*) x(i),y(j),omega(i,j)
    !omega(i,j) = 2.0d0*sin(x(i))*cos(y(j))
    !read(50,*) x(i),y(j),omega(i,j)
    omega(i,j) = W0/dcosh((y(j)+0.50d0*pi)/d)**2.0-W0/dcosh((y(j)-0.50d0*pi)/d)**2.0
    omega(i,j) = omega(i,j) - W0/dcosh((y(j)+1.50d0*pi)/d)**2.0+W0/dcosh((y(j)-1.50d0*pi)/d)**2.0
    omega(i,j) = omega(i,j)+0.01*dcos(m*x(i))
    omega_dum(i,j) = omega(i,j)
!    write(15,*) x(i),y(j),omega(i,j),ux(i,j),uy(i,j)
  end do
end do

!===================== INITIAL TIME DATA ===============================================

  call acc_set_device_num(0, 0)
  stream = acc_get_cuda_stream(acc_async_sync)
  call createCUFFTPlan2D(fftPlanD2ZMain,  Ny, Nx, CUFFT_D2Z, stream) !indices are swapped as fortran is column major and c is row major
  call createCUFFTPlan2D(fftPlanZ2DMain,  Ny, Nx, CUFFT_Z2D, stream) !indices are swapped as fortran is column major and c is row major

!  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, omega_dum, omegak, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_forward)

!$acc data copyin (x,y,omega_dum) create (psik_dum,psi, omega,Energy0,y_Energy0,C0,C1,C2,C3) &
!$acc & copyout(omegak,omegak_dum,psik,ukx,uky,ukx_dum,uky_dum,ux,uy) 

!$acc host_data use_device(omega_dum, omegak)
   call executeCUFFT2D(fftPlanD2ZMain, C_LOC(omega_dum), C_LOC(omegak), CUFFT_D2Z)
!$acc end host_data

!$acc parallel loop collapse(2)
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
!    write(25,*) i-1,j-1,abs(psik(i,j)),abs(omegak(i,j))
  enddo
enddo
!$acc end parallel


!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
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
!    write(25,*) i-1,j-1,abs(psik(i,j)),abs(omegak(i,j))
  enddo
enddo
!$acc end parallel

! call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, psik_dum, psi, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)

!  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)

!  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)

!  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
  
!  call dfftw_destroy_plan_ (plan_forward)
!  call dfftw_destroy_plan_ (plan_backward)

!$acc host_data use_device(psik_dum, omegak_dum, ukx_dum, uky_dum, psi, omega, ux, uy)  
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(psik_dum), C_LOC(psi), CUFFT_Z2D)
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(omegak_dum), C_LOC(omega), CUFFT_Z2D)
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(ukx_dum), C_LOC(ux), CUFFT_Z2D)
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(uky_dum), C_LOC(uy), CUFFT_Z2D)
!$acc end host_data

Energy0 = 0.0d0
y_Energy0 = 0.0d0
 C0 = 0.0d0
 C1 = 0.0d0
 C2 = 0.0d0
 C3 = 0.0d0
!$acc parallel loop collapse(2) reduction(+:Energy0,y_Energy0,C0,C1,C2,C3)
do i = 1, Nx
  do j = 1, Ny
  psi(i,j) = psi(i,j)/(dfloat(Nx)*dfloat(Ny))
  omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
  ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
  uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
  Energy0 = Energy0 + (ux(i,j)**2 + uy(i,j)**2)
  y_Energy0 = y_Energy0 + (uy(i,j)**2)
  C0 = C0 + omega(i,j)**0.0d0
  C1 = C1 + omega(i,j)**1.0d0
  C2 = C2 + omega(i,j)**2.0d0
  C3 = C3 + omega(i,j)**3.0d0
  
  end do
end do
!$acc end parallel
!$acc update host(omega,psi,ux,uy)
do i = 1,Nx
  do j = 1,Ny
    write(10,9) x(i),y(j),omega(i,j),psi(i,j)
    write(20,9) x(i),y(j),omega(i,j),ux(i,j),uy(i,j)
  enddo 
enddo
9 FORMAT( 12(2x,f22.16) ) 
!$acc end data
!!$acc end data
!write(40,*) 0, Energy0,y_Energy0, C0,C1,C2,C3
!======================= MAIN PROGRAM =====================================================
call cpu_time ( t1 )
do time = time_min,time_max,dt

t = nint(time/dt) - int(time_min/dt)

if (mod(t,100) == 0) then
print*, t/100
endif

!====================================================================================
  call derive (Nx,Ny,Nh,nux,nuy,time,omegak,omegak_new,dt_omegak_old,dt_omegak_new) 
!  call rk4 (Nx,Ny,Nh,nux,nuy,time,dt,omegak,omegak_new,dt_omegak_new) 
  call ab (Nx,Ny,Nh,nux,nuy,time,dt,omegak,omegak_new,dt_omegak_old,dt_omegak_new)
!====================================================================================

!$acc data copyin (omegak_new) create (psik,ukx,uky,ukx_dum,uky_dum,omegak_dum) &
!$acc & copyout(omega,ux,uy) 
!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
  do j = 1,Ny
    dt_omegak_old(i,j) = dt_omegak_new(i,j)
    omegak(i,j) = omegak_new(i,j)
  enddo
enddo
!$acc end parallel

!$acc parallel loop collapse(2)
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
!    write(25,*) i-1,j-1,abs(psik(i,j)),abs(omegak(i,j))
  enddo
enddo
!$acc end parallel

!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    psik(i,j) = omegak(i,j)/( kx*kx + ky*ky ) 
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j) 
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j) 
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
!    write(25,*) i-1,j-1,abs(psik(i,j)),abs(omegak(i,j))
  enddo
enddo
!$acc end parallel
!  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
  
!  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)

!  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
  
!  call dfftw_destroy_plan_ (plan_backward)
!$acc host_data use_device (omegak_dum, ukx_dum, uky_dum, omega, ux, uy)   
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(omegak_dum), C_LOC(omega), CUFFT_Z2D)
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(ukx_dum), C_LOC(ux), CUFFT_Z2D)
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(uky_dum), C_LOC(uy), CUFFT_Z2D)
!$acc end host_data
!if (t == 10000) then 
!do i = 1,Nx/2+1
!  do j = 1,Ny
!  Ek(i,j) = Ek(i,j) + sqrt(abs(ukx(i,j))**2 + abs(uky(i,j))**2)
!  write(35,*) i,j,sqrt(float((i-1)**2)+float((j-1)**2)),abs(Ek(i,j))
!  enddo
!enddo  
!endif 
!$acc end data


Energy = 0.0d0 
y_Energy = 0.0d0
 C0 = 0.0d0
 C1 = 0.0d0
 C2 = 0.0d0
 C3 = 0.0d0
!$acc parallel loop collapse(2) reduction(+:Energy,y_Energy,C0,C1,C2,C3)    
do i = 1,Nx
  do j = 1,Ny
  omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
  ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
  uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
  Energy = Energy + (ux(i,j)**2 + uy(i,j)**2)
  y_Energy = y_Energy + (uy(i,j)**2)
  C0 = C0 + omega(i,j)**0.0d0
  C1 = C1 + omega(i,j)**1.0d0
  C2 = C2 + omega(i,j)**2.0d0
  C3 = C3 + omega(i,j)**3.0d0
!    if (t == 10000 .or. t == 0) then
!    if (mod(dfloat(t),100.0) == 0.0) then
!    write(t+100,*) x(i),y(j),ux(i,j),uy(i,j),omega(i,j)
!    endif
  enddo
enddo  
!$acc end parallel

!!$acc update host(omega,ux,uy)
!do i = 1,Nx
!  do j = 1,Ny
!    if (mod(dfloat(t),100.0) == 0.0) then
!    write(t+100,7) x(i),y(j),ux(i,j),uy(i,j),omega(i,j)
!    endif
! enddo
!enddo
! 7 FORMAT( 12(2x,f22.16) ) 
 
write(40,*) time,Energy,y_Energy/(2.0d0*dfloat(Nx)*dfloat(Ny))!/y_Energy0!, C0,C1,C2,C3 

!call destroyCUFFTPlan2D(fftPlanD2ZMain)
!call destroyCUFFTPlan2D(fftPlanZ2DMain) 

enddo ! time
call cpu_time ( t2 )
write(5,*) "Time taken for the run = ",t2 - t1
 close(5)
 close(40)
 call destroyCUFFTPlan2D(fftPlanD2ZMain)
 call destroyCUFFTPlan2D(fftPlanZ2DMain)
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

integer (acc_handle_kind) :: stream
type (c_ptr) :: fftPlanD2Z, fftPlanZ2D

!  call acc_set_device_num(0, 0)
!  stream = acc_get_cuda_stream(acc_async_sync)
!  call createCUFFTPlan2D(fftPlanD2Z,  Ny, Nx, CUFFT_D2Z, stream) !indices are swapped as fortran is column major and c is row major
!  call createCUFFTPlan2D(fftPlanZ2D,  Ny, Nx, CUFFT_Z2D, stream) !indices are swapped as fortran is column major and c is row major
  
!$acc data copyin (omegak_dum) create (i_kx_omegak,i_ky_omegak,ux,uy) &
!$acc & create(domega_dx,domega_dy,ux_domega_dx,uy_domega_dy,NLkx,NLky,NLk,Nlk) &
!$acc & copyout(dt_omegak_new,psik,ukx,uky,ukx_dum,uky_dum,omegak) 
  
!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
  do j = 1,Ny/2
    omegak_dum(i,j) = omegak(i,j) 
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      ! De - Aliazing Technique...
      if (sqrt(kx*kx + ky*ky) .gt. (float(Nx+Ny)/2.0)/3.0 + 1) then    
      !omegak_dum(i,j) = (0.0d0,0.0d0) 
      endif 
  enddo
enddo 
!$acc end parallel

!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
  do j = Ny/2+1,Ny
    omegak_dum(i,j) = omegak(i,j) 
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      ! De - Aliazing Technique...
      if (sqrt(kx*kx + ky*ky) .gt. (float(Nx+Ny)/2.0)/3.0 + 1) then    
      !omegak_dum(i,j) = (0.0d0,0.0d0) 
      endif 
  enddo
enddo 
!$acc end parallel



!$acc parallel loop collapse(2)
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
!    write(25,*) i-1,j-1,abs(psik(i,j)),abs(omegak(i,j))
  enddo
enddo
!$acc end parallel

!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    psik(i,j) = omegak_dum(i,j)/( kx*kx + ky*ky ) 
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j) 
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j) 
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
!    write(25,*) i-1,j-1,abs(psik(i,j)),abs(omegak(i,j))
  enddo
enddo
!$acc end parallel


!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    i_kx_omegak(i,j) = (0.0d0,1.0d0)*kx*omegak_dum(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0)*ky*omegak_dum(i,j)
  enddo
enddo
!$acc end parallel

!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    i_kx_omegak(i,j) = (0.0d0,1.0d0)*kx*omegak_dum(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0)*ky*omegak_dum(i,j)
  enddo
enddo
!$acc end parallel

 ! call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
 ! call dfftw_execute_ (plan_backward)

 ! call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
 ! call dfftw_execute_ (plan_backward)

! call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_omegak, domega_dx, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)

!  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_omegak, domega_dy, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
 
!$acc host_data use_device (ukx_dum, uky_dum, i_kx_omegak, i_ky_omegak, ux, uy, domega_dx, domega_dy )     
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(ukx_dum), C_LOC(ux), CUFFT_Z2D)
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(uky_dum), C_LOC(uy), CUFFT_Z2D)
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(i_kx_omegak), C_LOC(domega_dx), CUFFT_Z2D)
  call executeCUFFT2D(fftPlanZ2DMain, C_LOC(i_ky_omegak), C_LOC(domega_dy), CUFFT_Z2D)
!$acc end host_data  

!$acc parallel loop collapse(2) 
do i = 1,Nx
  do j = 1,Ny
    ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_dx(i,j) = domega_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_dy(i,j) = domega_dy(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux_domega_dx(i,j) = ux(i,j)*domega_dx(i,j)
    uy_domega_dy(i,j) = uy(i,j)*domega_dy(i,j)
  enddo
enddo
!$acc end parallel

 ! call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_domega_dx, NLkx, FFTW_ESTIMATE)
 ! call dfftw_execute_ (plan_forward)

 ! call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_domega_dy, NLky, FFTW_ESTIMATE)
 ! call dfftw_execute_ (plan_forward)
  
 ! call dfftw_destroy_plan_ (plan_forward)
 ! call dfftw_destroy_plan_ (plan_backward)

!$acc host_data use_device( ux_domega_dx, uy_domega_dy, NLkx, NLky )  
  call executeCUFFT2D(fftPlanD2ZMain, C_LOC(ux_domega_dx), C_LOC(NLkx), CUFFT_D2Z)
  call executeCUFFT2D(fftPlanD2ZMain, C_LOC(uy_domega_dy), C_LOC(NLky), CUFFT_D2Z)
!$acc end host_data

!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      ! De - Aliazing Technique...
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 0) then    
      NLkx(i,j) = 0.0d0
      NLky(i,j) = 0.0d0
      endif 
  enddo
enddo
!$acc end parallel

!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      ! De - Aliazing Technique...
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 0) then
      NLkx(i,j) = 0.0d0
      NLky(i,j) = 0.0d0
      endif 
  enddo
enddo  
!$acc end parallel


!$acc parallel loop collapse(2) 
do i = 1,Nx/2+1
  do j = 1,Ny/2  
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    NLk(i,j) = ( NLkx(i,j) + NLky(i,j) )
    Nlk(i,j) = NLk(i,j) + (nux*kx*kx + nuy*ky*ky)*omegak(i,j)
  enddo
enddo 
!$acc end parallel

!$acc parallel loop collapse(2)
do i = 1,Nx/2+1 
  do j = Ny/2+1,Ny  
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    NLk(i,j) = ( NLkx(i,j) + NLky(i,j) )
    Nlk(i,j) = NLk(i,j) + (nux*kx*kx + nuy*ky*ky)*omegak(i,j)
  enddo
enddo   
!$acc end parallel

!$acc parallel loop collapse(2)
do i = 1,Nx/2+1
  do j = 1,Ny  
    dt_omegak_new(i,j) = -NLk(i,j)
  enddo
enddo   
!$acc end parallel
!$acc end data
return

end subroutine derive

!===================================================================================
!=================== SUBROUTINE RUNGE - KUTTA 4 ====================================
!===================================================================================

subroutine rk4(Nx,Ny,Nh,nux,nuy,time,dt,omegak,omegak_new,dt_omegak)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh
real ( kind = 8 ) time,dt,nux,nuy
complex ( kind = 8 ) komega1(Nh,Ny),komega2(Nh,Ny),komega3(Nh,Ny),komega4(Nh,Ny)
complex ( kind = 8 ) omegak(Nh,Ny),omegak_new(Nh,Ny),dum_omegak(Nh,Ny),dt_omegak(Nh,Ny)
complex ( kind = 8 ) dt_omegak_old(Nh,Ny)
complex ( kind = 8 ) dt_omegak_new(Nh,Ny)

do i = 1,Nh
  do j = 1,Ny
    komega1(i,j) = dt_omegak(i,j)
    dum_omegak(i,j) = omegak(i,j) + komega1(i,j)*dt/2.0
  end do
end do

  call derive(Nx,Ny,Nh,nux,nuy,time+dt/2.0,dum_omegak,omegak_new,dt_omegak_old,dt_omegak_new)

do i = 1,Nh
  do j = 1,Ny
    komega2(i,j) = dt_omegak(i,j)
    dum_omegak(i,j) = omegak(i,j) + komega2(i,j)*dt/2.0
  end do
end do

  call derive(Nx,Ny,Nh,nux,nuy,time+dt/2.0,dum_omegak,omegak_new,dt_omegak_old,dt_omegak_new)

do i = 1,Nh
  do j = 1,Ny
    komega3(i,j) = dt_omegak(i,j)
    dum_omegak(i,j) = omegak(i,j) + komega3(i,j)*dt/2.0
  end do
end do

  call derive(Nx,Ny,Nh,nux,nuy,time+dt,dum_omegak,omegak_new,dt_omegak_old,dt_omegak_new)

do i = 1,Nh
  do j = 1,Ny
    komega4(i,j) = dt_omegak(i,j)
    omegak_new(i,j) = omegak(i,j) + dt/6.0d0*(komega1(i,j) + 2.0d0*komega2(i,j) + 2.0d0*komega3(i,j) + komega4(i,j))
!    omegak_new(i,j) = omegak(i,j) + komega1(i,j)*dt     ! EULER SOLVER
  end do
end do

return

end subroutine rk4

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
!$acc data copyin (omegak,dt_omegak_new,dt_omegak_old) &
!$acc & copyout(omegak_new)
!$acc parallel loop collapse(2)
do i = 1,Nh
  do j = 1,Ny
    omegak_new(i,j) = omegak(i,j) + ( (3.0d0/2.0d0)*dt_omegak_new(i,j) - (1.0d0/2.0d0)*dt_omegak_old(i,j) )*dt
  end do
end do
!$acc end parallel
!$acc end data
return

end subroutine ab

!====================================================================================

end program PSHD2

