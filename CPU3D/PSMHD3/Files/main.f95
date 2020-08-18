Program PSMHD3
 
use omp_lib
implicit none

! Define Grid Size.
integer ( kind = 4 ), parameter :: Nx = 16
integer ( kind = 4 ), parameter :: Ny = 16
integer ( kind = 4 ), parameter :: Nz = 16
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

include "fftw3.f"

integer ( kind = 4 ) G,i,j,k,t
real ( kind = 8 ) Lx,Ly,Lz,dx,dy,dz,kx,ky,kz,time,time_min,time_max,dt,G1,G2,G3
real ( kind = 8 ) spheat,rho0,ms,U0,mu,kf,mode,a,sigma,uy0,M,CS,Re,Rm,PM,P0,eta,mu_0,MA,VA,B0,ba
real ( kind = 8 ) Pressure,Energy,Energy0,y_Energy,y_Energy0,B_Field,C0,C1,C2,C3,div_B_tot

real ( kind = 8 ) x(Nx),y(Ny),z(Nz),rho(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
real ( kind = 8 ) rho_ux(Nx,Ny,Nz),rho_uy(Nx,Ny,Nz),rho_uz(Nx,Ny,Nz)
real ( kind = 8 ) omega_x(Nx,Ny,Nz),omega_y(Nx,Ny,Nz),omega_z(Nx,Ny,Nz),omega2(Nx,Ny,Nz)
real ( kind = 8 ) P(Nx,Ny,Nz),E(Nx,Ny,Nz),Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz),div_B(Nx,Ny,Nz)
real ( kind = 8 ) rho_dum(Nx,Ny,Nz),ux_dum(Nx,Ny,Nz),uy_dum(Nx,Ny,Nz),uz_dum(Nx,Ny,Nz)
real ( kind = 8 ) rho_ux_dum(Nx,Ny,Nz),rho_uy_dum(Nx,Ny,Nz),rho_uz_dum(Nx,Ny,Nz)
real ( kind = 8 ) P_dum(Nx,Ny,Nz),E_dum(Nx,Ny,Nz),Bx_dum(Nx,Ny,Nz),By_dum(Nx,Ny,Nz),Bz_dum(Nx,Ny,Nz)

complex ( kind = 8 ) rho_k(Nh,Ny,Nz),ux_k(Nh,Ny,Nz),uy_k(Nh,Ny,Nz),uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_ux_k(Nh,Ny,Nz),rho_uy_k(Nh,Ny,Nz),rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) omega_x_k(Nh,Ny,Nz),omega_y_k(Nh,Ny,Nz),omega_z_k(Nh,Ny,Nz)
complex ( kind = 8 ) P_k(Nh,Ny,Nz),E_k(Nh,Ny,Nz),Ek(Nh,Ny,Nz),Bk(Nh,Ny,Nz)
complex ( kind = 8 ) Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz),div_B_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_k_dum(Nh,Ny,Nz),ux_k_dum(Nh,Ny,Nz),uy_k_dum(Nh,Ny,Nz),uz_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) omega_x_k_dum(Nh,Ny,Nz),omega_y_k_dum(Nh,Ny,Nz),omega_z_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) rho_ux_k_dum(Nh,Ny,Nz),rho_uy_k_dum(Nh,Ny,Nz),rho_uz_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) E_k_dum(Nh,Ny,Nz),Bx_k_dum(Nh,Ny,Nz),By_k_dum(Nh,Ny,Nz),Bz_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) rho_k_new(Nh,Ny,Nz),rho_ux_k_new(Nh,Ny,Nz),rho_uy_k_new(Nh,Ny,Nz),rho_uz_k_new(Nh,Ny,Nz)
complex ( kind = 8 ) E_k_new(Nh,Ny,Nz),Bx_k_new(Nh,Ny,Nz),By_k_new(Nh,Ny,Nz),Bz_k_new(Nh,Ny,Nz)

complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny,Nz),d_rho_ux_k_dt_old(Nh,Ny,Nz),d_rho_uy_k_dt_old(Nh,Ny,Nz),d_rho_uz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_old(Nh,Ny,Nz),d_Bx_k_dt_old(Nh,Ny,Nz),d_By_k_dt_old(Nh,Ny,Nz),d_Bz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny,Nz),d_rho_ux_k_dt_new(Nh,Ny,Nz),d_rho_uy_k_dt_new(Nh,Ny,Nz),d_rho_uz_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_new(Nh,Ny,Nz),d_Bx_k_dt_new(Nh,Ny,Nz),d_By_k_dt_new(Nh,Ny,Nz),d_Bz_k_dt_new(Nh,Ny,Nz)

integer ( kind = 8 ) plan_forward,plan_backward ! FFTW

integer ( kind = 4 ) thread_num,num_thread,proc_num ! OMP
real ( kind = 8 ) threads,t1,t2

common/comm/time_max,Lx,Ly,Lz,dx,dy,dz,spheat,mu,ms,CS,mu_0,eta,dt,kf

integer,parameter :: seed = 99999999
call srand(seed)

!===================== FILENAMES ==============================================	

open(unit=1,file='INPUT.DATA',status='old')
open(unit=5,file='System_information.dat',status='unknown')
open(unit=15,file='Initial_Grid_Data.dat',status='unknown')
!open(unit=20,file='INPUT.dat',status='old')  ! This input file is the file generated from Vorticity code. 
!open(unit=25,file='Initial_Grid_Data_Reproduced.dat',status='unknown')
!open(unit=30,file='Initial_Energy_Spectra.dat',status='unknown')
open(unit=35,file='Energy_Spectra.dat',status='unknown')
open(unit=40,file='Energy.dat',status='unknown')

!===================== USER INPUTS ============================================		

10 FORMAT(100X, I1)
20 FORMAT(36X, F10.3)

! Define Number of Threads.
proc_num = omp_get_num_procs()
read (1,10) G
read (1,20) threads
thread_num = int(threads)
call omp_set_num_threads (thread_num)

write(5,*) "Number of Processors=", proc_num, "Number of Threads=", thread_num  

! System Size.
read (1,10) G
read (1,10) G
read (1,20) Lx
read (1,20) Ly
read (1,20) Lz
Lx = Lx*pi; Ly = Ly*pi; Lz = Lz*pi;

! Grid Resolution.
dx = Lx/dfloat(Nx); dy = Ly/dfloat(Ny); dz = Lz/dfloat(Nz)

! Runtime Details and Time Step Width.
read (1,10) G
read (1,10) G
read (1,20) time_min 
read (1,20) time_max
read (1,20) dt   

! Ratio of Specific Heats/ Adiabetic Index/ Heat Capacity Ratio/ Poisson Constant.
read (1,10) G
read (1,20) spheat 

! Co-efficient of Viscosity.
read (1,10) G
read (1,20) mu

! Co-efficient of Resistivity.
read (1,10) G
read (1,20) eta

! Magnetic Permiability.
read (1,10) G
read (1,20) mu_0

! Mass of the Fluid Element.
read (1,10) G
read (1,20) ms

! Background Density.
read (1,10) G
read (1,20) rho0

! Initial Pressure.
read (1,10) G
read (1,20) P0

! Mach Number.
!M = 0.50d0
read (1,10) G
read (1,20) M

! Alfven Mach Number.
read (1,10) G
read (1,20) MA

! Forcing Length Scale.
read (1,10) G
read (1,20) kf

! Maximum Velocity.
!U0 = 0.14732247502913884d0
!U0 = 0.6450d0 

! Sound Speed.
 CS = dsqrt(spheat*P0/rho0)
 U0 = M * CS

! Sound Speed.
! CS = U0/M

! Alfven Speed.
VA = U0/MA

! Initial Magnetic Field.
B0 = VA*dsqrt(mu_0*rho0) 

write(5,*) "Sound Speed=", CS, "Initial Velocity=", U0
write(5,*) "Alfven Speed=", VA, "Initial Magnetic Field=", B0

! Raynold's Number.
Re = ms*rho0*U0*Lx/mu; write(5,*) "Re =", Re
Rm = U0*Lx/eta; write(5,*) "Magnetic Raynold's Number =", Rm
PM = mu/eta; write(5,*) "Prandtl Number =", PM

! Grid Generation.
do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0d0+real(j-1)*dy
    do k = 1, Nz
      z(k)=0.0d0+real(k-1)*dz
      ! Initial Density Distribution.
      rho(i,j,k) = rho0
      ! Initial Velocity Distribution.
      ux(i,j,k) = U0 * rand() ! * (dtanh((y(j)-Ly/3.0d0)/a) - dtanh((y(j)-2.0d0*Ly/3.0d0)/a) - 1.0d0) 
      !ux(i,j,k) = U0*tanh((y(j)-Ly/2.0d0)/a)
      uy(i,j,k) = U0 * rand() !0.0d0
      ! Initial Velocity Perturbation.
      !uy(i,j,k) = uy(i,j,k) + uy0*dsin(mode*x(i))*dexp(-((y(j)-Ly/3.0d0)**2)/(sigma**2))
      !uy(i,j,k) = uy(i,j,k) + uy0*dsin(mode*x(i))*dexp(-((y(j)-2.0d0*Ly/3.0d0)**2)/(sigma**2))
      !uy(i,j,k) = uy(i,j,k) + uy0*dsin(mode*x(i))*exp(-((y(j)-Ly/2.0d0)**2)/(sigma**2))
      uz(i,j,k) = U0 * rand() !0.0d0
      ! Initial Velocity Distribution.
      !read(20,*) G1,G2,G3,ux(i,j,k),uy(i,j,k),uz(i,j,k)  ! You may define your own function here.
      ! Initial Pressure Distribution.
      P(i,j,k) = P0
      ! Initial Magnetic Field Profile.
      Bx(i,j,k) = B0 * rand() ! * dtanh((y(j)-Ly/3.0d0)/ba)
      By(i,j,k) = B0 * rand() !0.0d0
      Bz(i,j,k) = B0 * rand() !0.0d0
      ! Initial Energy Distribution.
      E(i,j,k) = P(i,j,k)/(spheat -1.0d0) + 0.50d0*rho(i,j,k) * (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)!&
      ! + 0.50d0*(Bx(i,j)**2 + By(i,j)**2 + Bz(i,j)**2)
      ! Initial Combined Variables.
      rho_ux(i,j,k) = rho(i,j,k) * ux(i,j,k)
      rho_uy(i,j,k) = rho(i,j,k) * uy(i,j,k)
      rho_uz(i,j,k) = rho(i,j,k) * uz(i,j,k)
      ! Keep Backup of the Arrays for FFTW. 
      rho_dum(i,j,k) = rho(i,j,k)
      ux_dum(i,j,k) = ux(i,j,k)
      uy_dum(i,j,k) = uy(i,j,k)
      uz_dum(i,j,k) = uz(i,j,k)
      rho_ux_dum(i,j,k) = rho_ux(i,j,k)
      rho_uy_dum(i,j,k) = rho_uy(i,j,k)
      rho_uz_dum(i,j,k) = rho_uz(i,j,k)
      P_dum(i,j,k) = P(i,j,k)
      E_dum(i,j,k) = E(i,j,k) 
      Bx_dum(i,j,k) = Bx(i,j,k)
      By_dum(i,j,k) = By(i,j,k)
      Bz_dum(i,j,k) = Bz(i,j,k)
      ! Write Initial Density and Velocity Distribution in File.
      write(15,*) x(i),y(j),z(k),rho(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k),P(i,j,k),E(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    enddo
  end do
end do

 close(15)

!===================== INITIAL TIME DATA ===============================================

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_dum, rho_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_ux_dum, rho_ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_uy_dum, rho_uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_uz_dum, rho_uz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, P_dum, P_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, E_dum, E_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Bx_dum, Bx_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, By_dum, By_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Bz_dum, Bz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

! Evaluate Initial Vorticity Spectra.
do i = 1, Nx/2+1
  do j = 1, Ny/2
    do k = 1, Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1, Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo 
 
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_x_k_dum, omega_x, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_y_k_dum, omega_y, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_z_k_dum, omega_z, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_destroy_plan_ (plan_forward)
  call dfftw_destroy_plan_ (plan_backward)

! Evaluate Initial Energy Spectra.
!do i = 1, Nx/2+1
  !do j = 1, Ny
    !do k = 1, Nz
      !write(30,*) i-1,j-1,k-1,abs(nk(i,j,k)),abs(ux_k(i,j,k)),abs(uy_k(i,j,k)),abs(uz_k(i,j,k)) 
    !enddo
  !end do
!end do

!Energy0 = 0.0d0; y_Energy0 = 0.0d0
 
do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      ! FFTW Normalisation.
       omega_x(i,j,k) = omega_x(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
       omega_y(i,j,k) = omega_y(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
       omega_z(i,j,k) = omega_z(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Initial Kinetic Energy.
      !Energy0 = Energy0 + (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)
      !y_Energy0 = y_Energy0 + (uy(i,j,k)**2)
      ! Store Grid Data Reproduced from FFTW.
      !write(25,*) x(i),y(j),z(k),omega_x(i,j,k),omega_y(i,j,k),omega_z(i,j,k), &
      !             rho(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k), &
      !             P(i,j,k),E(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    end do
  end do
enddo

! Write the Initial Energy in File.
!write(40,*) 0.0d0, Energy0,y_Energy0

!======================= MAIN PROGRAM =====================================================
t1 = omp_get_wtime()

! Time Loop Starts.
do time = time_min,time_max,dt

t = nint(time/dt) - dint(time_min/dt)

! Online Check of the Progress.
!if (mod(t,int(1.0/dt)) == 0) then
  !print*, t/int(1.0/dt)
!endif

!======================================================================================

! Non-Linear Term Evaluator.
  call DERIVE (Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
               d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
               d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
               d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)
  
! Time Solvers.
! Adams-Bashforth
  call ADAMS_BASHFORTH (Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
           rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new, &
           d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
           d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
           d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)
  
!======================================================================================

!$OMP PARALLEL SHARED(d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old),&
!$OMP & SHARED(d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old),&
!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new),& 
!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new),& 
!$OMP & SHARED(rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k),&
!$OMP & SHARED(rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new),&
!$OMP & SHARED(rho_k_dum,rho_ux_k_dum,rho_uy_k_dum,rho_uz_k_dum,E_k_dum,Bx_k_dum,By_k_dum,Bz_k_dum) PRIVATE(i,j,k)
!$OMP DO 

do i = 1,Nx/2+1
  do j = 1,Ny
    do k = 1,Nz
      ! Set the Variables in Proper Format for Next Time Iteration.
      d_rho_k_dt_old(i,j,k) = d_rho_k_dt_new(i,j,k)
    
      d_rho_ux_k_dt_old(i,j,k) = d_rho_ux_k_dt_new(i,j,k)
      d_rho_uy_k_dt_old(i,j,k) = d_rho_uy_k_dt_new(i,j,k)
      d_rho_uz_k_dt_old(i,j,k) = d_rho_uz_k_dt_new(i,j,k)
    
      d_E_k_dt_old(i,j,k) = d_E_k_dt_new(i,j,k)
    
      d_Bx_k_dt_old(i,j,k) = d_Bx_k_dt_new(i,j,k)
      d_By_k_dt_old(i,j,k) = d_By_k_dt_new(i,j,k)
      d_Bz_k_dt_old(i,j,k) = d_Bz_k_dt_new(i,j,k)
    
      rho_k(i,j,k) = rho_k_new(i,j,k)
    
      rho_ux_k(i,j,k) = rho_ux_k_new(i,j,k)
      rho_uy_k(i,j,k) = rho_uy_k_new(i,j,k)
      rho_uz_k(i,j,k) = rho_uz_k_new(i,j,k)
  
      E_k(i,j,k) = E_k_new(i,j,k)
    
      Bx_k(i,j,k) = Bx_k_new(i,j,k)
      By_k(i,j,k) = By_k_new(i,j,k)
      Bz_k(i,j,k) = Bz_k_new(i,j,k)
    
      ! Keep Backup of the Arrays for FFTW.
      rho_k_dum(i,j,k) = rho_k(i,j,k)   
    
      rho_ux_k_dum(i,j,k) = rho_ux_k(i,j,k)
      rho_uy_k_dum(i,j,k) = rho_uy_k(i,j,k)
      rho_uz_k_dum(i,j,k) = rho_uz_k(i,j,k)

      E_k_dum(i,j,k) = E_k(i,j,k)
    
      Bx_k_dum(i,j,k) = Bx_k(i,j,k)
      By_k_dum(i,j,k) = By_k(i,j,k)
      Bz_k_dum(i,j,k) = Bz_k(i,j,k)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL 

! Evaluate Divergence of B in Spectral Space.

!$OMP PARALLEL SHARED(Lx,Ly,Lz,Bx_k,By_k,Bz_k,div_B_k) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO
 
do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k) 
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL 

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_k_dum, rho, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_ux_k_dum, rho_ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uy_k_dum, rho_uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uz_k_dum, rho_uz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, E_k_dum, E, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bx_k_dum, Bx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, By_k_dum, By, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bz_k_dum, Bz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, div_B_k, div_B, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
!$OMP PARALLEL SHARED(spheat,rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz,div_B), &
!$OMP & SHARED(ux,uy,uz,ux_dum,uy_dum,uz_dum,P) PRIVATE(i,j,k)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      ! FFTW Normalisation.
      rho(i,j,k) = rho(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
   
      rho_ux(i,j,k) = rho_ux(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      rho_uy(i,j,k) = rho_uy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz)) 
      rho_uz(i,j,k) = rho_uz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
   
      E(i,j,k) = E(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
      Bx(i,j,k) = Bx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      By(i,j,k) = By(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Bz(i,j,k) = Bz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

      ! Evaluate Divergence of B.
      div_B(i,j,k) = div_B(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
      ! Evaluate Velocity in Real Space.
      ux(i,j,k) = rho_ux(i,j,k)/rho(i,j,k)
      uy(i,j,k) = rho_uy(i,j,k)/rho(i,j,k)
      uz(i,j,k) = rho_uz(i,j,k)/rho(i,j,k)
  
      ! Evaluate Square of Magnetic Field.
      B2(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)
  
      ! Keep Backup of the Arrays for FFTW.
      ux_dum(i,j,k) = ux(i,j,k)
      uy_dum(i,j,k) = uy(i,j,k)
      uz_dum(i,j,k) = uz(i,j,k)
  
      ! Evaluate Pressure
      P(i,j,k) = ( spheat - 1.0d0 ) * ( E(i,j,k) - 0.50d0 * &
                 ( rho_ux(i,j,k)*ux(i,j,k)+rho_uy(i,j,k)*uy(i,j,k)+rho_uz(i,j,k)*uz(i,j,k) ) )! - B2(i,j) ) )
    enddo
  enddo
enddo   

!$OMP END DO
!$OMP END PARALLEL 

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

! Evaluate Vorticity in Spectral Space.

!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,omega_x_k,omega_y_k,omega_z_k),&
!$OMP & SHARED(omega_x_k_dum,omega_y_k_dum,omega_z_k_dum) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO 

do i = 1, Nx/2+1
  do j = 1, Ny/2
    do k = 1, Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
    do k = Nz/2+1, Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo 
  enddo
  do j = Ny/2+1, Ny
    do k = 1, Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
    do k = Nz/2+1, Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo 

!$OMP END DO
!$OMP END PARALLEL
 
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_x_k_dum, omega_x, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_y_k_dum, omega_y, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_z_k_dum, omega_z, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_destroy_plan_ (plan_forward)
  call dfftw_destroy_plan_ (plan_backward)

! FFTW Normalisation and omega^2 Evaluation.

!$OMP PARALLEL SHARED(omega_x,omega_y,omega_z,omega2) PRIVATE(i,j,k)
!$OMP DO 

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      omega_x(i,j,k) = omega_x(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega_y(i,j,k) = omega_y(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega_z(i,j,k) = omega_z(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega2(i,j,k) = omega_x(i,j,k)**2 + omega_y(i,j,k)**2 + omega_z(i,j,k)**2
    enddo
  enddo
enddo  

!$OMP END DO
!$OMP END PARALLEL

! Evaluate Energy Spectra at a Specified Time.
if (t >= int((time_max-time_max/10.0d0)/dt)) then

!!! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file. !!!
! !!$OMP PARALLEL SHARED(ux_k,uy_k,uz_k) PRIVATE(i,j,k) 
! !!$OMP DO REDUCTION (+:Ek) 

  do i = 1,Nx/2+1
    do j = 1,Ny
      do k = 1,Nz
        Ek(i,j,k) = Ek(i,j,k) + sqrt(abs(ux_k(i,j,k))**2 + abs(uy_k(i,j,k))**2 + abs(uz_k(i,j,k))**2)&
                    /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
        Bk(i,j,k) = Bk(i,j,k) + sqrt(abs(Bx_k(i,j,k))**2 + abs(By_k(i,j,k))**2 + abs(Bz_k(i,j,k))**2)&
                    /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      enddo
    enddo
  enddo  

! !!$OMP END PARALLEL

endif 
 
Pressure = 0.0d0; Energy = 0.0d0; y_Energy = 0.0d0; B_Field = 0.0d0
 C0 = 0.0d0; C1 = 0.0d0; C2 = 0.0d0; C3 = 0.0d0; div_B_Tot = 0.0d0
   
! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file.
! !$OMP PARALLEL SHARED(mu_0,t,dt,x,y,z,omega_x,omega_y,omega_z,omega2,rho,ux,uy,uz) PRIVATE(i,j,k)
! !$OMP DO REDUCTION (+:Pressure,Energy,y_Energy,B_Field,C0,C1,C2,C3,div_B_Tot) 

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      ! Evaluate Pressure
      Pressure = Pressure + P(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Energy
      !Energy = Energy + E(i,j,k)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Energy = Energy + (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Growth Rate.
      y_Energy = y_Energy + rho(i,j,k)*(uy(i,j,k)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Magnetic Field.
      B_Field = B_Field + (Bx(i,j,k)**2+By(i,j,k)**2+Bz(i,j,k)**2)/(2.0d0*mu_0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Casimirs.
      C0 = C0 + dsqrt(omega2(i,j,k))**0.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C1 = C1 + dsqrt(omega2(i,j,k))**1.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C2 = C2 + dsqrt(omega2(i,j,k))**2.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C3 = C3 + dsqrt(omega2(i,j,k))**3.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Check for Div B = 0
      div_B_Tot = div_B_Tot + div_B(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Write Grid Data in State Files.  
        !if (mod(float(t),100.0) == 0.0) then
        !write(t+100,*) x(i),y(j),z(k),omega_x(i,j,k),omega_y(i,j,k),omega_z(i,j,k), &
        !                rho(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k), &
        !                P(i,j,k),E(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
        !endif
    enddo
  enddo
enddo  

! !$OMP END PARALLEL

write(40,*) time,Pressure,Energy,y_Energy/2.0d0,B_Field,C0,C1,C2,C3,div_B_Tot
  
! close(t+100)
 
enddo ! time

t2 = omp_get_wtime()

write(5,*) "Time taken for the run =",(t2 - t1)/(60.0d0*60.0d0),"Hours"

do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    kz = 2.0d0*pi*float(k-1)/Lz
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
      Ek(i,j,k) = 0.0d0
      Bk(i,j,k) = 0.0d0
      endif
    write(35,*) i,j,k,sqrt(float((i-1)**2)+float((j-1)**2)+float((k-1)**2)),&
              abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*float(i-1)/Lx
    ky = 2.0d0*pi*float(j-1)/Ly
    kz = 2.0d0*pi*float((k-1)-Nz)/Lz
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
      Ek(i,j,k) = 0.0d0
      Bk(i,j,k) = 0.0d0
      endif
    write(35,*) i,j,k-Nz,sqrt(float((i-1)**2)+float((j-1)**2)+float(((k-1)-Nz)**2)),&
              abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
  enddo  
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    kz = 2.0d0*pi*float(k-1)/Lz
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
      Ek(i,j,k) = 0.0d0
      Bk(i,j,k) = 0.0d0
      endif
    write(35,*) i,j-Ny,k,sqrt(float((i-1)**2)+float(((j-1)-Ny)**2)+float((k-1)**2)),&
              abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*float(i-1)/Lx
    ky = 2.0d0*pi*float((j-1)-Ny)/Ly
    kz = 2.0d0*pi*float((k-1)-Nz)/Lz
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
      Ek(i,j,k) = 0.0d0
      Bk(i,j,k) = 0.0d0
      endif
    write(35,*) i,j-Ny,k-Nz,sqrt(float((i-1)**2)+float(((j-1)-Ny)**2)+float(((k-1)-Nz)**2)),&
              abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
  enddo  
enddo 
 
 close(5)
 close(35)
 close(40)

!====================================================================================

end program PSMHD3  

