Program MHD3D
use, intrinsic :: iso_c_binding
implicit none
include 'mpif.h'

! Define Grid Size.
integer (C_INTPTR_T), parameter :: Nx = 8
integer (C_INTPTR_T), parameter :: Ny = 8
integer (C_INTPTR_T), parameter :: Nz = 8
integer (C_INTPTR_T), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

include "fftw3-mpi.f03"

integer(C_INTPTR_T) :: i,j,k,l,m,n,o,t,P_s,P_max,q,q_max
real(C_DOUBLE) :: Lx,Ly,Lz,dx,dy,dz,h,dr,delta,kx,ky,kz,time,time_min,time_max,dt,G1,G2,G3,A,B,C,tmp
real(C_DOUBLE) :: spheat,rho0,U0,mu,kf,mode,sigma,uy0,MS,CS,Re,Rm,PM,P0,eta,mu_0,MA,VA,B0,ba,t1,t2
real(C_DOUBLE) :: Pressure,Energy,Energy0,T_Energy,y_Energy0,B_Field,C0,C1,C2,C3,div_B_tot,Hf,Hm,HT,HuB,HAw,Rayleigh

real(C_DOUBLE) :: x(Nx),y(Ny),z(Nz)
real(C_DOUBLE), dimension(Nx,Ny,Nz) :: r,rho,ux,uy,uz,u2,rho_ux,rho_uy,rho_uz
real(C_DOUBLE), dimension(Nx,Ny,Nz) :: jx,jy,jz,Ax,Ay,Az,j2,A2,omega_x,omega_y,omega_z,omega2
real(C_DOUBLE), dimension(Nx,Ny,Nz) :: P,E,Bx,By,Bz,B2,div_B,rho_dum,ux_dum,uy_dum,uz_dum
real(C_DOUBLE), dimension(Nx,Ny,Nz) :: rho_ux_dum,rho_uy_dum,rho_uz_dum,P_dum,E_dum,Bx_dum,By_dum,Bz_dum

real(C_DOUBLE), dimension(Nx,Ny,Nz) :: Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3
real(C_DOUBLE), dimension(Nx,Ny,Nz) :: d_ux_dx,d_uy_dy,d_uz_dz,Fx,Fy,Fz,Energy_x,Energy_y,Energy_z,E_Visc
real(C_DOUBLE), dimension(Nx,Ny,Nz) :: Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2
real(C_DOUBLE), dimension(Nx,Ny,Nz) :: d_Bx_dy,d_By_dx,d_Bx_dz,d_By_dz,d_Bz_dx,d_Bz_dy,curl_x_B,curl_y_B,curl_z_B

complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: rho_k,ux_k,uy_k,uz_k,rho_ux_k,rho_uy_k,rho_uz_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: omega_x_k,omega_y_k,omega_z_k,jx_k,jy_k,jz_k,Ax_k,Ay_k,Az_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: P_k,E_k,Ek,Bk,Bx_k,By_k,Bz_k,div_B_k,rho_k_dum,ux_k_dum,uy_k_dum,uz_k_dum
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: omega_x_k_dum,omega_y_k_dum,omega_z_k_dum
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: rho_ux_k_dum,rho_uy_k_dum,rho_uz_k_dum
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: E_k_dum,Bx_k_dum,By_k_dum,Bz_k_dum
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: E_k_new,Bx_k_new,By_k_new,Bz_k_new

complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new

complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: nu,Fx_k,Fy_k,Fz_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: Mom_x_1_k,Mom_x_2_k,Mom_x_3_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: Mom_y_1_k,Mom_y_2_k,Mom_y_3_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: Mom_z_1_k,Mom_z_2_k,Mom_z_3_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: i_kx_ux_k,i_ky_uy_k,i_kz_uz_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: Energy_x_k,Energy_y_k,Energy_z_k,E_Visc_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: Mag_x_1_k,Mag_x_2_k,Mag_y_1_k,Mag_y_2_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: Mag_z_1_k,Mag_z_2_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: i_ky_Bx_k,i_kx_By_k,i_kx_Bz_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: i_kz_Bx_k,i_ky_Bz_k,i_kz_By_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k
complex(C_DOUBLE_COMPLEX), dimension(Nh,Ny,Nz) :: kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k

integer :: ierr, myid, nproc, seed, root
type(C_PTR) :: planf, planb, cdatar, cdatac
integer(C_INTPTR_T) :: alloc_local, local_Nz, local_k_offset
real(C_DOUBLE), pointer :: idata(:,:,:)
complex(C_DOUBLE_complex), pointer :: odata(:,:,:)

!===================== FILENAMES ==============================================	

open(unit=5,file='System_information.dat',status='unknown')
open(unit=15,file='Initial_Grid_Data.dat',status='unknown')
!open(unit=20,file='INPUT.dat',status='old')  ! This input file is the file generated from Vorticity code. 
!open(unit=25,file='Initial_Grid_Data_Reproduced.dat',status='unknown')
!open(unit=30,file='Initial_Energy_Spectra.dat',status='unknown')
open(unit=35,file='Energy_Spectra.dat',status='unknown')
open(unit=40,file='Energy.dat',status='unknown')
open(unit=50,file='Structure_Function.dat',status='unknown')
open(unit=60,file='State.dat',status='unknown')

!===================== USER INPUTS ============================================		

  call mpi_init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

root = 0

! System Size.
Lx = 2.0d0*pi; Ly = 2.0d0*pi; Lz = 2.0d0*pi 

! Grid Resolution.
dx = Lx/dfloat(Nx); dy = Ly/dfloat(Ny); dz = Lz/dfloat(Nz)

! Runtime Details and Time Step Width.
time_min = 0.0d0
time_max = 100.00d0
dt = 0.00010d0

! Ratio of Specific Heats/ Adiabetic Index/ Heat Capacity Ratio/ Poisson Constant.
spheat = 1.0d0

! Co-efficient of Viscosity.
Re = 1000.0d0

! Co-efficient of Resistivity.
Rm = 1000.0d0

! Background Density.
rho0 = 1.0d0

! Initial Pressure.
P0 = 1.0d0

! Mach Number.
MS = 0.10d0

! Alfven Mach Number.
MA = 1.0d0

! Maximum Velocity.
U0 = 0.10d0

! Sound Speed.
 CS = 1.0d0/MS

! Alfven Speed.
VA = U0/MA

! Initial Magnetic Field.
B0 = VA*dsqrt(rho0) 

! Forcing Length Scale.
kf = 1.0d0

 A = 1.0d0
 B = 1.0d0
 C = 1.0d0

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
      ux(i,j,k) = U0 * ( A*dsin(kf*z(k)) + C*dcos(kf*y(j)) )
      uy(i,j,k) = U0 * ( B*dsin(kf*x(i)) + A*dcos(kf*z(k)) )
      uz(i,j,k) = U0 * ( C*dsin(kf*y(j)) + B*dcos(kf*x(i)) )
      ! Initial Pressure Distribution.
      P(i,j,k) = P0
      ! Initial Magnetic Field Profile.
      Bx(i,j,k) = B0!*(A*dsin(kf*z(k)) + C*dcos(kf*y(j)))
      By(i,j,k) = B0!*(B*dsin(kf*x(i)) + A*dcos(kf*z(k)))
      Bz(i,j,k) = B0!*(C*dsin(kf*y(j)) + B*dcos(kf*x(i)))
      ! Initial Energy Distribution.
      E(i,j,k) = P(i,j,k) + 0.50d0*rho(i,j,k) * (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)&
                 + 0.50d0*(Bx(i,j,k)**2 + By(i,j,k)**2 + Bz(i,j,k)**2)
      ! Evaluate Forcing.
      Fx(i,j,k) = rho(i,j,k) * ( A*dsin(kf*z(k)) + C*dcos(kf*y(j)) )
      Fy(i,j,k) = rho(i,j,k) * ( B*dsin(kf*x(i)) + A*dcos(kf*z(k)) )
      Fz(i,j,k) = rho(i,j,k) * ( C*dsin(kf*y(j)) + B*dcos(kf*x(i)) )
    enddo
  end do
end do

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
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
      write(15,*) ux(i,j,k),uy(i,j,k),uz(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    enddo
  end do
end do

 close(15)

!===================== FFTW MPI PLAN CREATION ==========================================

  call fftw_mpi_init()
  
  ! get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_3d(Nz, Ny, Nh, MPI_COMM_WORLD, local_Nz, local_k_offset)
  cdatar = fftw_alloc_real(2 * alloc_local)
  cdatac = fftw_alloc_complex(alloc_local)
  
  call c_f_pointer(cdatar, idata, [2*Nh,Ny,local_Nz])
  call c_f_pointer(cdatac, odata, [Nh,Ny,local_Nz])

  ! create MPI plan for out-of-place DFT (note dimension reversal)
  planf = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, idata, odata, MPI_COMM_WORLD, FFTW_MEASURE)
  planb = fftw_mpi_plan_dft_c2r_3d(Nz, Ny, Nx, odata, idata, MPI_COMM_WORLD, FFTW_MEASURE)
                                                                 
!===================== INITIAL TIME DATA ===============================================

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = rho_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      rho_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = ux_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      ux_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = uy_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      uy_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = uz_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      uz_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = rho_ux_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      rho_ux_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = rho_uy_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      rho_uy_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = rho_uz_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      rho_uz_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = P_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      P_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = E_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      E_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Bx_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Bx_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = By_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      By_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Bz_dum(i, j, k + local_k_offset)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Bz_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

! Evaluate Initial Vorticity Spectra.
do i = 1, Nh
  do j = 1, Ny/2
    do k = 1, local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
        omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
        omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
        omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
        omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
        endif
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1, local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
        omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
        omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
        omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
        omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
        endif
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo 
 
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = omega_x_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      omega_x(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = omega_y_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      omega_y(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = omega_z_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      omega_z(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

!======================= MAIN PROGRAM =====================================================

t1 = MPI_Wtime() 

! Time Loop Starts.
do time = time_min,time_max,dt

t = nint(time/dt) - dint(time_min/dt)

! Online Check of the Progress.
!if (mod(t,int(1.0/dt)) == 0) then
  !print*, t/int(1.0/dt)
!endif

! Keep Backup of the Arrays for FFTW.

do i = 1,Nh
  do j = 1,Ny
    do k = 1,local_Nz
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

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = rho_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      rho(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = rho_ux_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      rho_ux(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = rho_uy_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      rho_uy(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = rho_uz_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      rho_uz(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do


  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = E_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      E(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = Bx_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      Bx(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = By_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      By(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = Bz_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      Bz(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

do i = 1,Nx
  do j = 1,Ny
    do k = 1,local_Nz
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
    enddo
  enddo
enddo   

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = ux_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      ux_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = uy_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      uy_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = uz_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      uz_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Fx(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Fx_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Fy(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Fy_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Fz(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Fz_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

! Evaluate derivatives of Velocity and Magnetic Field.

do i = 1,Nh
  do j = 1,Ny/2
    do k = 1,local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
        i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
        i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)
      
        i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
        i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
        i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
        i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
        i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
        i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
        i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
        i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)
      
        i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
        i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
        i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
        i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
        i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
        i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
        endif
    enddo 
  enddo
  do j = Ny/2+1,Ny
    do k = 1,local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
        i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
        i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)
      
        i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
        i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
        i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
        i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
        i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
        i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
        i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
        i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)
      
        i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
        i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
        i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
        i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
        i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
        i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
        endif
    enddo
  enddo   
enddo

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = i_kx_ux_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      d_ux_dx(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = i_ky_uy_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      d_uy_dy(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = i_kz_uz_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      d_uz_dz(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = i_ky_Bx_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      d_Bx_dy(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = i_kz_Bx_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      d_Bx_dz(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = i_kx_By_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      d_By_dx(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = i_kz_By_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      d_By_dz(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = i_kx_Bz_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      d_Bz_dx(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = i_ky_Bz_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      d_Bz_dy(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
do i = 1,Nx
  do j = 1,Ny
    do k = 1,local_Nz
    ! Evaluate Curl of Magnetic Field.
    curl_x_B(i,j,k) = d_Bz_dy(i,j,k) - d_By_dz(i,j,k)
    curl_y_B(i,j,k) = d_Bz_dx(i,j,k) - d_Bx_dz(i,j,k)
    curl_z_B(i,j,k) = d_By_dx(i,j,k) - d_Bx_dy(i,j,k)
    
    ! Evaluate Pressure
    P(i,j,k) = CS*CS*rho(i,j,k)!( spheat - 1.0d0 ) * ( E(i,j,k) &
               !- 0.50d0 * ( rho_ux(i,j,k)*ux(i,j,k)+rho_uy(i,j,k)*uy(i,j,k)+rho_uz(i,j,k)*uz(i,j,k) &
               !- B2(i,j,k) ) )
  
    ! Evaluate LHS of Momentum Equation.
    Mom_x_1(i,j,k) = rho_ux(i,j,k)*ux(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 - Bx(i,j,k)*Bx(i,j,k)
    Mom_x_2(i,j,k) = rho_ux(i,j,k)*uy(i,j,k) - Bx(i,j,k)*By(i,j,k)
    Mom_x_3(i,j,k) = rho_ux(i,j,k)*uz(i,j,k) - Bx(i,j,k)*Bz(i,j,k)
  
    Mom_y_1(i,j,k) = rho_ux(i,j,k)*uy(i,j,k) - Bx(i,j,k)*By(i,j,k) 
    Mom_y_2(i,j,k) = rho_uy(i,j,k)*uy(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 - By(i,j,k)*By(i,j,k)
    Mom_y_3(i,j,k) = rho_uy(i,j,k)*uz(i,j,k) - By(i,j,k)*Bz(i,j,k)
  
    Mom_z_1(i,j,k) = rho_uz(i,j,k)*ux(i,j,k) - Bz(i,j,k)*Bx(i,j,k) 
    Mom_z_2(i,j,k) = rho_uz(i,j,k)*uy(i,j,k) - Bz(i,j,k)*By(i,j,k) 
    Mom_z_3(i,j,k) = rho_uz(i,j,k)*uz(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 - Bz(i,j,k)*Bz(i,j,k)
  
    ! Evaluate LHS of Energy Equation.
    Energy_x(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * ux(i,j,k) &
                      - ux(i,j,k)*Bx(i,j,k)*Bx(i,j,k) - uy(i,j,k)*Bx(i,j,k)*By(i,j,k) - uz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) &
                      - eta * ( By(i,j,k) * curl_z_B(i,j,k) + Bz(i,j,k) * curl_y_B(i,j,k)) 
                      
    Energy_y(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * uy(i,j,k) &
                      - ux(i,j,k)*Bx(i,j,k)*By(i,j,k) - uy(i,j,k)*By(i,j,k)*By(i,j,k) - uz(i,j,k)*By(i,j,k)*Bz(i,j,k) &
                      + eta * ( Bx(i,j,k) * curl_z_B(i,j,k) - Bz(i,j,k) * curl_x_B(i,j,k))
                      
    Energy_z(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * uz(i,j,k) &
                      - ux(i,j,k)*Bz(i,j,k)*Bx(i,j,k) - uy(i,j,k)*Bz(i,j,k)*By(i,j,k) - uz(i,j,k)*Bz(i,j,k)*Bz(i,j,k) & 
                      + eta * ( Bx(i,j,k) * curl_y_B(i,j,k) + By(i,j,k) * curl_x_B(i,j,k))
                      
    ! Evaluate RHS of Energy Equation.
    E_Visc(i,j,k) = ( d_ux_dx(i,j,k) + d_uy_dy(i,j,k) + d_uz_dz(i,j,k) )**2
  
    ! Evaluate LHS of Magnetic Field Equation.
    Mag_x_1(i,j,k) = ux(i,j,k)*By(i,j,k) - uy(i,j,k)*Bx(i,j,k) 
    Mag_x_2(i,j,k) = ux(i,j,k)*Bz(i,j,k) - uz(i,j,k)*Bx(i,j,k)
  
    Mag_y_1(i,j,k) = ux(i,j,k)*By(i,j,k) - uy(i,j,k)*Bx(i,j,k)
    Mag_y_2(i,j,k) = uy(i,j,k)*Bz(i,j,k) - uz(i,j,k)*By(i,j,k)
  
    Mag_z_1(i,j,k) = ux(i,j,k)*Bz(i,j,k) - uz(i,j,k)*Bx(i,j,k) 
    Mag_z_2(i,j,k) = uy(i,j,k)*Bz(i,j,k) - uz(i,j,k)*By(i,j,k)
    enddo
  enddo
enddo   

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mom_x_1(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mom_x_1_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mom_x_2(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mom_x_2_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mom_x_3(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mom_x_3_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mom_y_1(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mom_y_1_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mom_y_2(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mom_y_2_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mom_y_3(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mom_y_3_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mom_z_1(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mom_z_1_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mom_z_2(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mom_z_2_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mom_z_3(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mom_z_3_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Energy_x(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Energy_x_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Energy_y(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Energy_y_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Energy_z(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Energy_z_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = E_Visc(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      E_Visc_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mag_x_1(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mag_x_1_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mag_x_2(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mag_x_2_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mag_y_1(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mag_y_1_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mag_y_2(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mag_y_2_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mag_z_1(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mag_z_1_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = Mag_z_2(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      Mag_z_2_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

! Evaluate the Derivatives in Spectral Space.

do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
        i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k) 
        i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k) 
        
        i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)    
        i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
        i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)   
        i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
        i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
        i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k) 
        i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)    
        i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)    
        i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)
     
        kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
        ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
        kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
        kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
        ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k) 
        kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
        kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
        ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
        kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

        i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
        i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
        i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)
    
        i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
        i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
        i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
        i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
        i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
        i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)
    
        kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
        ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
        kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
        kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
        ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k) 
        kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
        kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
        ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
        kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
        i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k) 
        i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k) 
        
        i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)    
        i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
        i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)   
        i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
        i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
        i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k) 
        i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)    
        i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)    
        i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)
     
        kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
        ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
        kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
        kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
        ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k) 
        kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
        kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
        ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
        kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

        i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
        i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
        i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)
    
        i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
        i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
        i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
        i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
        i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
        i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)
    
        kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
        ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
        kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
        kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
        ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k) 
        kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
        kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
        ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
        kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
        endif
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1,local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
        i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k) 
        i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k) 
        
        i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)    
        i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
        i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)   
        i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
        i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
        i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k) 
        i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)    
        i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)    
        i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)
     
        kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
        ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
        kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
        kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
        ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k) 
        kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
        kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
        ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
        kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

        i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
        i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
        i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)
    
        i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
        i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
        i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
        i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
        i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
        i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)
    
        kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
        ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
        kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
        kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
        ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k) 
        kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
        kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
        ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
        kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
        i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k) 
        i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k) 
        
        i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)    
        i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
        i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)   
        i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
        i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
        i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k) 
        i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)    
        i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)    
        i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)
     
        kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
        ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
        kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
        kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
        ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k) 
        kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
        kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
        ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
        kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

        i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
        i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
        i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)
    
        i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
        i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
        i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
        i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
        i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
        i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)
    
        kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
        ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
        kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
        kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
        ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k) 
        kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
        kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
        ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
        kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
        endif

    enddo
  enddo
enddo

! De - Aliazing Technique With 2/3 Rule for All the Non-Linear Terms.  

do i = 1,Nh
  do j = 1,Ny/2
    do k = 1,local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        endif
        if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k(i,j,k) = 0.0d0
        i_ky_rho_uy_k(i,j,k) = 0.0d0
        i_kz_rho_uz_k(i,j,k) = 0.0d0
        i_kx_Mom_x_1_k(i,j,k) = 0.0d0
        i_ky_Mom_x_2_k(i,j,k) = 0.0d0
        i_kz_Mom_x_3_k(i,j,k) = 0.0d0
        i_kx_Mom_y_1_k(i,j,k) = 0.0d0
        i_ky_Mom_y_2_k(i,j,k) = 0.0d0
        i_kz_Mom_y_3_k(i,j,k) = 0.0d0
        i_kx_Mom_z_1_k(i,j,k) = 0.0d0
        i_ky_Mom_z_2_k(i,j,k) = 0.0d0
        i_kz_Mom_z_3_k(i,j,k) = 0.0d0
        kx2_ux_k(i,j,k) = 0.0d0
        ky2_ux_k(i,j,k) = 0.0d0
        kz2_ux_k(i,j,k) = 0.0d0
        kx2_uy_k(i,j,k) = 0.0d0
        ky2_uy_k(i,j,k) = 0.0d0
        kz2_uy_k(i,j,k) = 0.0d0
        kx2_uz_k(i,j,k) = 0.0d0
        ky2_uz_k(i,j,k) = 0.0d0
        kz2_uz_k(i,j,k) = 0.0d0
        i_kx_Energy_x_k(i,j,k) = 0.0d0
        i_ky_Energy_y_k(i,j,k) = 0.0d0
        i_kz_Energy_z_k(i,j,k) = 0.0d0
        E_Visc_k(i,j,k) = 0.0d0
        i_ky_Mag_x_1_k(i,j,k) = 0.0d0
        i_kz_Mag_x_2_k(i,j,k) = 0.0d0
        i_kx_Mag_y_1_k(i,j,k) = 0.0d0
        i_kz_Mag_y_2_k(i,j,k) = 0.0d0
        i_kx_Mag_z_1_k(i,j,k) = 0.0d0
        i_ky_Mag_z_2_k(i,j,k) = 0.0d0
        kx2_Bx_k(i,j,k) = 0.0d0
        ky2_Bx_k(i,j,k) = 0.0d0
        kz2_Bx_k(i,j,k) = 0.0d0
        kx2_By_k(i,j,k) = 0.0d0
        ky2_By_k(i,j,k) = 0.0d0
        kz2_By_k(i,j,k) = 0.0d0
        kx2_Bz_k(i,j,k) = 0.0d0
        ky2_Bz_k(i,j,k) = 0.0d0
        kz2_Bz_k(i,j,k) = 0.0d0
        nu(i,j,k) = 0.0d0
        endif
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1,local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        endif
        if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k(i,j,k) = 0.0d0
        i_ky_rho_uy_k(i,j,k) = 0.0d0
        i_kz_rho_uz_k(i,j,k) = 0.0d0
        i_kx_Mom_x_1_k(i,j,k) = 0.0d0
        i_ky_Mom_x_2_k(i,j,k) = 0.0d0
        i_kz_Mom_x_3_k(i,j,k) = 0.0d0
        i_kx_Mom_y_1_k(i,j,k) = 0.0d0
        i_ky_Mom_y_2_k(i,j,k) = 0.0d0
        i_kz_Mom_y_3_k(i,j,k) = 0.0d0
        i_kx_Mom_z_1_k(i,j,k) = 0.0d0
        i_ky_Mom_z_2_k(i,j,k) = 0.0d0
        i_kz_Mom_z_3_k(i,j,k) = 0.0d0
        kx2_ux_k(i,j,k) = 0.0d0
        ky2_ux_k(i,j,k) = 0.0d0
        kz2_ux_k(i,j,k) = 0.0d0
        kx2_uy_k(i,j,k) = 0.0d0
        ky2_uy_k(i,j,k) = 0.0d0
        kz2_uy_k(i,j,k) = 0.0d0
        kx2_uz_k(i,j,k) = 0.0d0
        ky2_uz_k(i,j,k) = 0.0d0
        kz2_uz_k(i,j,k) = 0.0d0
        i_kx_Energy_x_k(i,j,k) = 0.0d0
        i_ky_Energy_y_k(i,j,k) = 0.0d0
        i_kz_Energy_z_k(i,j,k) = 0.0d0
        E_Visc_k(i,j,k) = 0.0d0
        i_ky_Mag_x_1_k(i,j,k) = 0.0d0
        i_kz_Mag_x_2_k(i,j,k) = 0.0d0
        i_kx_Mag_y_1_k(i,j,k) = 0.0d0
        i_kz_Mag_y_2_k(i,j,k) = 0.0d0
        i_kx_Mag_z_1_k(i,j,k) = 0.0d0
        i_ky_Mag_z_2_k(i,j,k) = 0.0d0
        kx2_Bx_k(i,j,k) = 0.0d0
        ky2_Bx_k(i,j,k) = 0.0d0
        kz2_Bx_k(i,j,k) = 0.0d0
        kx2_By_k(i,j,k) = 0.0d0
        ky2_By_k(i,j,k) = 0.0d0
        kz2_By_k(i,j,k) = 0.0d0
        kx2_Bz_k(i,j,k) = 0.0d0
        ky2_Bz_k(i,j,k) = 0.0d0
        kz2_Bz_k(i,j,k) = 0.0d0
        nu(i,j,k) = 0.0d0
        endif
    enddo
  enddo
enddo  

do i = 1,Nh
  do j = 1,Ny
    do k = 1,local_Nz
      ! Density Equation.
      d_rho_k_dt_new(i,j,k) = - ( i_kx_rho_ux_k(i,j,k) + i_ky_rho_uy_k(i,j,k) + i_kz_rho_uz_k(i,j,k) )
    
      ! Momentum Equation.
      d_rho_ux_k_dt_new(i,j,k) = - ( i_kx_Mom_x_1_k(i,j,k) + i_ky_Mom_x_2_k(i,j,k) + i_kz_Mom_x_3_k(i,j,k) ) &
                                 - ( kx2_ux_k(i,j,k) + ky2_ux_k(i,j,k) + kz2_ux_k(i,j,k) ) / Re !+ Fx_k(i,j,k) 
    
      d_rho_uy_k_dt_new(i,j,k) = - ( i_kx_Mom_y_1_k(i,j,k) + i_ky_Mom_y_2_k(i,j,k) + i_kz_Mom_y_3_k(i,j,k) ) &
                                 - ( kx2_uy_k(i,j,k) + ky2_uy_k(i,j,k) + kz2_uy_k(i,j,k) ) / Re !+ Fy_k(i,j,k)
    
      d_rho_uz_k_dt_new(i,j,k) = - ( i_kx_Mom_z_1_k(i,j,k) + i_ky_Mom_z_2_k(i,j,k) + i_kz_Mom_z_3_k(i,j,k) ) &
                                 - ( kx2_uz_k(i,j,k) + ky2_uz_k(i,j,k) + kz2_uz_k(i,j,k) ) / Re !+ Fz_k(i,j,k) 
    
      ! Energy Equation.
      d_E_k_dt_new(i,j,k) = - ( i_kx_Energy_x_k(i,j,k) + i_ky_Energy_y_k(i,j,k) + i_kz_Energy_z_k(i,j,k) ) &
                            + mu * E_Visc_k(i,j,k)
    
      ! Magnetic Field Equation.
      d_Bx_k_dt_new(i,j,k) = + ( i_ky_Mag_x_1_k(i,j,k) + i_kz_Mag_x_2_k(i,j,k) ) &
                             - ( kx2_Bx_k(i,j,k) + ky2_Bx_k(i,j,k) + kz2_Bx_k(i,j,k) ) / Rm
                           
      d_By_k_dt_new(i,j,k) = - ( i_kx_Mag_y_1_k(i,j,k) - i_kz_Mag_y_2_k(i,j,k) ) &
                             - ( kx2_By_k(i,j,k) + ky2_By_k(i,j,k) + kz2_By_k(i,j,k) ) / Rm
                           
      d_Bz_k_dt_new(i,j,k) = - ( i_kx_Mag_z_1_k(i,j,k) + i_ky_Mag_z_2_k(i,j,k) ) &                       
                             - ( kx2_Bz_k(i,j,k) + ky2_Bz_k(i,j,k) + kz2_Bz_k(i,j,k) ) / Rm
    enddo                       
  enddo
enddo  

do i = 1,Nh
  do j = 1,Ny
    do k = 1,local_Nz
      ! Density Equation Evolution.
      rho_k_new(i,j,k) = rho_k(i,j,k) + ( (3.0d0/2.0d0)*d_rho_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_k_dt_old(i,j,k) )*dt
      
      ! Momentum Equation Evolution.
      rho_ux_k_new(i,j,k) = rho_ux_k(i,j,k) + ( (3.0d0/2.0d0)*d_rho_ux_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_ux_k_dt_old(i,j,k) )*dt
      rho_uy_k_new(i,j,k) = rho_uy_k(i,j,k) + ( (3.0d0/2.0d0)*d_rho_uy_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_uy_k_dt_old(i,j,k) )*dt
      rho_uz_k_new(i,j,k) = rho_uz_k(i,j,k) + ( (3.0d0/2.0d0)*d_rho_uz_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_uz_k_dt_old(i,j,k) )*dt
      
      ! Energy Equation Evolution.
      E_k_new(i,j,k) = E_k(i,j,k) !+ ( (3.0d0/2.0d0)*d_E_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_E_k_dt_old(i,j,k) )*dt 
      
      ! Energy Equation Evolution.
      Bx_k_new(i,j,k) = Bx_k(i,j,k) + ( (3.0d0/2.0d0)*d_Bx_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_Bx_k_dt_old(i,j,k) )*dt
      By_k_new(i,j,k) = By_k(i,j,k) + ( (3.0d0/2.0d0)*d_By_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_By_k_dt_old(i,j,k) )*dt
      Bz_k_new(i,j,k) = Bz_k(i,j,k) + ( (3.0d0/2.0d0)*d_Bz_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_Bz_k_dt_old(i,j,k) )*dt
    enddo
  end do
end do

do i = 1,Nh
  do j = 1,Ny
    do k = 1,local_Nz
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

do i = 1,Nh
  do j = 1,Ny/2
    do k = 1,local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
        jx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k) - (0.0d0,1.0d0)*kz*By_k(i,j,k)
        jy_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k) - (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
        jz_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k) - (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
        jx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k) - (0.0d0,1.0d0)*kz*By_k(i,j,k)
        jy_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k) - (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
        jz_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k) - (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
        endif
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1,local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k) 
        jx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k) - (0.0d0,1.0d0)*kz*By_k(i,j,k)
        jy_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k) - (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
        jz_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k) - (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k) 
        jx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k) - (0.0d0,1.0d0)*kz*By_k(i,j,k)
        jy_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k) - (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
        jz_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k) - (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
        endif
    enddo
  enddo
enddo

do i = 1,Nh
  do j = 1,Ny/2
    do k = 1,local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (kx .eq. 0.0d0 .and. ky .eq. 0.0d0 .and. kz .eq. 0.0d0) then
        Ax_k(i,j,k) = jx_k(i,j,k)
        Ay_k(i,j,k) = jy_k(i,j,k)
        Az_k(i,j,k) = jz_k(i,j,k)
        elseif (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        Ax_k(i,j,k) = jx_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        Ay_k(i,j,k) = jy_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        Az_k(i,j,k) = jz_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        Ax_k(i,j,k) = jx_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        Ay_k(i,j,k) = jy_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        Az_k(i,j,k) = jz_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        endif
    enddo  
  enddo
  do j = Ny/2+1,Ny
    do k = 1,local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        Ax_k(i,j,k) = jx_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        Ay_k(i,j,k) = jy_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        Az_k(i,j,k) = jz_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        Ax_k(i,j,k) = jx_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        Ay_k(i,j,k) = jy_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        Az_k(i,j,k) = jz_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        endif
    enddo  
  enddo
enddo

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = jx_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      jx(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = jy_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      jy(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = jz_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      jz(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = Ax_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      Ax(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = Ay_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      Ay(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = Az_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      Az(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = rho_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      rho(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = rho_ux_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      rho_ux(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = rho_uy_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      rho_uy(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = rho_uz_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      rho_uz(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = E_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      E(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = Bx_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      Bx(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = By_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      By(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = Bz_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      Bz(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do
  
  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = div_B_k(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      div_B(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

do i = 1,Nx
  do j = 1,Ny
    do k = 1,local_Nz
      ! Evaluate Velocity in Real Space.
      ux(i,j,k) = rho_ux(i,j,k)/rho(i,j,k)
      uy(i,j,k) = rho_uy(i,j,k)/rho(i,j,k)
      uz(i,j,k) = rho_uz(i,j,k)/rho(i,j,k)
  
      ! Evaluate Square of Velocity and Magnetic Field.
      u2(i,j,k) = ux(i,j,k)*ux(i,j,k) + uy(i,j,k)*uy(i,j,k) + uz(i,j,k)*uz(i,j,k)
      B2(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)
      j2(i,j,k) = jx(i,j,k)*jx(i,j,k) + jy(i,j,k)*jy(i,j,k) + jz(i,j,k)*jz(i,j,k)
      A2(i,j,k) = Ax(i,j,k)*Ax(i,j,k) + Ay(i,j,k)*Ay(i,j,k) + Az(i,j,k)*Az(i,j,k)
  
      ! Keep Backup of the Arrays for FFTW.
      rho_dum(i,j,k) = rho(i,j,k)
      ux_dum(i,j,k) = ux(i,j,k)
      uy_dum(i,j,k) = uy(i,j,k)
      uz_dum(i,j,k) = uz(i,j,k)
  
      ! Evaluate Pressure
      P(i,j,k) = CS*CS*rho(i,j,k)!( spheat - 1.0d0 ) * ( E(i,j,k) - 0.50d0 * &
                 !( rho_ux(i,j,k)*ux(i,j,k)+rho_uy(i,j,k)*uy(i,j,k)+rho_uz(i,j,k)*uz(i,j,k) - B2(i,j,k) ) )
    enddo
  enddo
enddo   

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = ux_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      ux_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = uy_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      uy_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      idata(i, j, k) = uz_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      uz_k(i,j,k) = odata(i,j,k)
      end do
    end do
  end do

! Evaluate Vorticity in Spectral Space.

do i = 1, Nh
  do j = 1, Ny/2
    do k = 1, local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
        omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
        omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
        omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
        omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
        endif
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
  do j = Ny/2+1, Ny
    do k = 1, local_Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        if (k+local_k_offset-1 <= Nz/2) then
        kz = 2.0d0*pi*float(k+local_k_offset-1)/Lz
        omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
        omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
        omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
        else 
        kz = 2.0d0*pi*float(k+local_k_offset-1-Nz)/Lz
        omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
        omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
        omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
        endif
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo 

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = omega_x_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      omega_x(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = omega_y_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      omega_y(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

  do i = 1, Nh
    do j = 1, Ny
      do k = 1, local_Nz
      odata(i, j, k) = omega_z_k_dum(i, j, k)
      end do
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, Ny
      do k = 1, local_Nz
      omega_z(i,j,k) = idata(i,j,k)/dfloat(Nx*Ny*Nz)
      end do
    end do
  end do

! omega^2 Evaluation.

do i = 1,Nx
  do j = 1,Ny
    do k = 1,local_Nz
      omega2(i,j,k) = omega_x(i,j,k)**2 + omega_y(i,j,k)**2 + omega_z(i,j,k)**2
    enddo
  enddo
enddo  

Pressure = 0.0d0; Energy = 0.0d0; T_Energy = 0.0d0; B_Field = 0.0d0
 C0 = 0.0d0; C1 = 0.0d0; C2 = 0.0d0; C3 = 0.0d0; div_B_Tot = 0.0d0
 Hf = 0.0d0; Hm = 0.0d0; HT = 0.0d0; HuB = 0.0d0; HAw = 0.0d0; Rayleigh = 0.0d0
   
do i = 1,Nx
  do j = 1,Ny
    do k = 1,local_Nz
      ! Evaluate Pressure
      Pressure = Pressure + P(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Energy
      Energy = Energy + u2(i,j,k)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Growth Rate.
      !T_Energy = T_Energy + E(i,j,k)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Magnetic Field.
      B_Field = B_Field + B2(i,j,k)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Casimirs.
      C0 = C0 + dsqrt(omega2(i,j,k))**0.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C1 = C1 + dsqrt(omega2(i,j,k))**1.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C2 = C2 + dsqrt(omega2(i,j,k))**2.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C3 = C3 + dsqrt(omega2(i,j,k))**3.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Helicities.
      Hf = Hf + ( ux(i,j,k)*omega_x(i,j,k) + uy(i,j,k)*omega_y(i,j,k) + uz(i,j,k)*omega_z(i,j,k) ) &
                /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Hm = Hm + ( Ax(i,j,k)*Bx(i,j,k) + Ay(i,j,k)*By(i,j,k) + Az(i,j,k)*Bz(i,j,k) )/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))  
      HT = HT + ( (ux(i,j,k)+Ax(i,j,k))*(omega_x(i,j,k)+Bx(i,j,k)) + (uy(i,j,k)+Ay(i,j,k))*(omega_y(i,j,k)+By(i,j,k)) &
                + (uz(i,j,k)+Az(i,j,k))*(omega_z(i,j,k)+Bz(i,j,k)) )/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz)) 
      HuB = HuB + ( ux(i,j,k)*Bx(i,j,k) + uy(i,j,k)*By(i,j,k) + uz(i,j,k)*Bz(i,j,k) )/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      HAw = HAw + ( Ax(i,j,k)*omega_x(i,j,k) + Ay(i,j,k)*omega_y(i,j,k) + Az(i,j,k)*omega_z(i,j,k) ) &
                /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))            
      ! Check for Div B = 0
      div_B_Tot = div_B_Tot + div_B(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Rayleigh Quotient.
      Rayleigh = Rayleigh + (omega2(i,j,k)+0.50d0*j2(i,j,k))/( (u2(i,j,k)+0.50d0*B2(i,j,k)) * dfloat(Nx)*dfloat(Ny)*dfloat(Nz) )
      ! Write Grid Data in State Files.  
      !  if (mod(float(t),1000.0) == 0.0) then
      !  write(t+110,*) ux(i,j,k),uy(i,j,k),uz(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
      !  endif
      !  if (mod(float(t),1000.0) == 0.0 .and. i == 1) then
      !  write(t+120,*) y(j),z(k),u2(i,j,k),B2(i,j,k),omega2(i,j,k),j2(i,j,k),A2(i,j,k)
      !  elseif (mod(float(t),1000.0) == 0.0 .and. j == 1) then
      !  write(t+130,*) x(i),z(k),u2(i,j,k),B2(i,j,k),omega2(i,j,k),j2(i,j,k),A2(i,j,k)
      !  elseif (mod(float(t),1000.0) == 0.0 .and. k == 1) then
      !  write(t+140,*) x(i),y(j),u2(i,j,k),B2(i,j,k),omega2(i,j,k),j2(i,j,k),A2(i,j,k)
      !  endif
    enddo
  enddo
enddo  

  CALL MPI_REDUCE(Pressure, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  Pressure = tmp

  CALL MPI_REDUCE(Energy, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  Energy = tmp
  
  CALL MPI_REDUCE(B_Field, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  B_Field = tmp
  
  CALL MPI_REDUCE(C0, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  C0 = tmp
  
  CALL MPI_REDUCE(C1, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  C1 = tmp
  
  CALL MPI_REDUCE(C2, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  C2 = tmp
  
  CALL MPI_REDUCE(C3, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  C3 = tmp

  CALL MPI_REDUCE(Hf, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  Hf = tmp
  
  CALL MPI_REDUCE(Hm, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  Hm = tmp

  CALL MPI_REDUCE(HT, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  HT = tmp

  CALL MPI_REDUCE(HuB, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  HuB = tmp

  CALL MPI_REDUCE(HAw, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  HAw = tmp

  CALL MPI_REDUCE(Rayleigh, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  Rayleigh = tmp
  
  CALL MPI_REDUCE(div_B_Tot, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  div_B_Tot = tmp
 
if (myid == root) then  
  if (mod(float(t),1000.0) == 0.0) then
  write(40,*) time,Pressure,Energy,B_Field,C1,C2,C3,Hf,Hm,HT,HuB,HAw,Rayleigh,div_B_Tot
  call flush(40) 
  endif
endif

!if (mod(float(t),1000.0) == 0.0) then
!  close(t+100)
!  close(t+110)
!  close(t+120)
!  close(t+130)
!  close(t+140)
!endif
 
enddo ! time

t2 = MPI_Wtime()

write(5,*) "Elapsed time is", t2 - t1

! close(5)
! close(40)
! close(60)

  ! deallocate and destroy plans  
  call fftw_destroy_plan(planf)
  call fftw_destroy_plan(planb)
  call fftw_mpi_cleanup()
  call fftw_free(cdatar)
  call fftw_free(cdatac)

  ! FINALIZE MPI

  CALL MPI_FINALIZE(IERR)
  
end program MHD3D  

