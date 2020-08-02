!Module for invoking cufft

module cufft3D
    integer, public :: CUFFT_FORWARD = -1
    integer, public :: CUFFT_INVERSE = 1
    integer, public :: CUFFT_R2C = Z'2a' ! Real to Complex (interleaved)
    integer, public :: CUFFT_C2R = Z'2c' ! Complex (interleaved) to Real
    integer, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved
    integer, public :: CUFFT_D2Z = Z'6a' ! Double to Double-Complex
    integer, public :: CUFFT_Z2D = Z'6c' ! Double-Complex to Double
    integer, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex

    interface
        subroutine createCUFFTPlan3D(plan, nx, ny, nz, planType, stream) bind (C, name = 'createCUFFTPlan3D')
            use iso_c_binding
            use openacc
            implicit none
            type (c_ptr) :: plan
            integer (c_int), value :: nx, ny, nz, planType
            integer (acc_handle_kind), value :: stream
        end subroutine createCUFFTPlan3D
    end interface

    interface
        subroutine executeCUFFT3D(plan, iData, oData, planType) bind (C, name = 'executeCUFFT3D')
            use iso_c_binding
            use openacc
            implicit none
            type (c_ptr) ::  plan
            type (c_ptr), value :: iData
            type (c_ptr), value :: oData
            integer (c_int), value :: planType
        end subroutine executeCUFFT3D
    end interface

    interface
        subroutine destroyCUFFTPlan3D(plan) bind (C, name = 'destroyCUFFTPlan3D')
            use iso_c_binding
            use openacc
            implicit none
            type (c_ptr) :: plan
        end subroutine destroyCUFFTPlan3D
    end interface

end module cufft3D



Program PSMHD3
!
! Rupak, IPR, August 6, 2017
!
! This is an OPENMP Parallel Benchmarked Three Dimensional Compressible Viscous Resistive MHD code, 
! using Pseudo-Spectral Method with Multiple Time Solvers.
! Here Adams-Bashforth Subroutine is Parallalised. Other Time Solvers are NOT Parallalised.
!
! Wave-form /sim exp(i(kx-wt)). Hence derivative w.r.t. x gives ik only (and NOT -ik).
!
!___________________________________________________________________________________________________________________________________________
!
! To Run this code in udbhav.ipr.res.in: 
!
!#!/bin/bash
!#Give your job name
!#BSUB -J MHD2D
!#Give required number of processor
!#BSUB -n 15
!#Give your email id
!##BSUB -u <rupak@ipr.res.in>
!#Give your current working directory
!##BSUB -cwd ~/<your home area>
!#Redirect error to some file and %J is Job ID
!#BSUB -e err.%J
!#Redirect output to some file and %J is Job ID
!#BSUB -o out.%J
!#Define number of processor per host user wants to execute a job
!#BSUB -R "span[hosts=1]"
!gfortran -fopenmp -I/home/rupak/bin/include 3dmhd.f95 -lfftw3 -lm
!./a.out
!
! Terminal Command: bsub < runjob.sh
!___________________________________________________________________________________________________________________________________________
!
! To Run this code in uday.ipr.res.in:
!
!#!/bin/bash
!# declare a name for this job to be sample_job
!#PBS -N MHD2D
!# request the parallel queue for this job
!##PBS -q parallel
!# request a total of 25 processors for this job (1 node and 25 processors per node)
!#PBS -l nodes=1:ppn=25
!# Request to run for specified hrs:min:secs
!##PBS -l walltime=500:00:00
!# combine PBS standard output and error files
!#PBS -j oe
!# mail is sent to you when the job starts and when it terminates or aborts
!##PBS -m bea
!# specify your email address
!#PBS -M rupak@gmail.com
!#change to the directory where you submitted the job
!cd $PBS_O_WORKDIR
!#include the full path to the name of your OMP program
!gfortran -fopenmp -I/soft/fftw-3.3.3/include -L/soft/fftw-3.3.3/lib 3dmhd.f95 -lfftw3 -lm
!./a.out
!
! Terminal Command: qsub jobfile.sh
!___________________________________________________________________________________________________________________________________________
!
! To Run this code in Desktop: 
!
!gfortran -fopenmp -I/usr/local/include 2dmhd.f95 -lfftw3 -lm
!./a.out
!___________________________________________________________________________________________________________________________________________
!
!\documentclass[12pt]{article} 
!\usepackage{a4wide,color}
!\begin{document}
!
!MHD Equations in Conservative form:\\
!
!\begin{eqnarray*}
!&& \frac{\partial \rho}{\partial t} + \vec{\nabla} \cdot \left(\rho \vec{u}\right) = 0\\
!&& \frac{\partial (\rho \vec{u})}{\partial t} + \vec{\nabla} \cdot \left[ \rho \vec{u} \vec{u} + \left(P + \frac{B^2}{2}\right){\bf{I}} - \vec{B}\vec{B} \right] = \mu \nabla^2 \vec{u}\\
!&& \frac{\partial E}{\partial t} + \vec{\nabla} \cdot \left[ \left( E + P + \frac{B^2}{2} \right)\vec{u} - \vec{u}\cdot\left( \vec{B} \vec{B} \right)  - \eta \vec{B} \times \\
!\left(\vec{\nabla} \times \vec{B} \right) \right] = \mu \left(\vec{\nabla} \cdot \vec{u} \right)^2\\
!&& \frac{\partial \vec{B}}{\partial t} + \vec{\nabla} \cdot \left( \vec{u} \vec{B} - \vec{B} \vec{u}\right) = \eta \nabla^2 \vec{B}\\
!\end{eqnarray*}
!
!In Two Dimensions, the MHD Equations in Conservative form becomes:\\
!
!\begin{eqnarray*}
!&& \frac{\partial \rho}{\partial t} + \frac{\partial}{\partial x} (\rho u_x) + \frac{\partial}{\partial y} (\rho u_y) = 0\\
!&& ~\\
!&& ~\\
!&& \frac{\partial \left( \rho u_x \right)}{\partial t} + \frac{\partial}{\partial x} \left[ \rho u_x u_x + P + \frac{B_x^2+B_y^2}{2} - B_x B_x \right] + \frac{\partial}{\partial y} \\
!\left[ \rho u_x u_y - B_x B_y \right] = \mu \left( \frac{\partial^2 u_x}{\partial x^2} + \frac{\partial^2 u_x}{\partial y^2} \right)^2\\
!&& \frac{\partial \left( \rho u_y \right)}{\partial t} + \frac{\partial}{\partial x} \left[ \rho u_x u_y - B_x B_y \right] + \frac{\partial}{\partial y} \left[ \rho u_y u_y + P + \\
!\frac{B_x^2+B_y^2}{2} - B_y B_y \right] = \mu \left( \frac{\partial^2 u_y}{\partial x^2} + \frac{\partial^2 u_y}{\partial y^2} \right)^2\\
!&& ~\\
!&& ~\\
!&& \frac{\partial E}{\partial t} + \frac{\partial}{\partial x} \left[ \left( E + P + \frac{B_x^2+B_y^2}{2} \right)u_x - u_x B_x B_x - u_y B_x B_y - \eta B_y \left( \frac{\partial B_y}{\partial x} - \frac{\partial B_x}{\partial y} \right) \right] \\
!&& ~~~~~ + \frac{\partial}{\partial y} \left[ \left( E + P + \frac{B_x^2+B_y^2}{2} \right)u_y - u_x B_x B_y - u_y B_y B_y + \eta B_x \left( \frac{\partial B_y}{\partial x} - \\
!\frac{\partial B_x}{\partial y} \right) \right] = \mu \left( \frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} \right)^2\\
!&& ~\\
!&& ~\\
!&& \frac{\partial B_x}{\partial t} + \frac{\partial}{\partial y} \left( u_y B_x - B_y u_x \right) = \eta \left( \frac{\partial^2 B_x}{\partial x^2} + \frac{\partial^2 B_x}{\partial y^2} \right)\\
!&& \frac{\partial B_y}{\partial t} + \frac{\partial}{\partial x} \left( u_x B_y - B_x u_y \right) = \eta \left( \frac{\partial^2 B_y}{\partial x^2} + \frac{\partial^2 B_y}{\partial y^2} \right)\\
!\end{eqnarray*}
!
!\newpage
!
!Pseudocode for 2DMHD Code.
!
!\begin{eqnarray*}
!&& [At~time~ (t)] \\
!&& \\
!\rho, \rho u_x, \rho u_y, E, B_x, B_y ~~~~~ \Leftarrow ~~~ && \hat{\rho}, \widehat{\rho u_x}, \widehat{\rho u_y}, \hat{E}, \hat{B}_x, \hat{B}_y \longrightarrow \hat{\omega} = i k_x \hat{u}_y - i k_y \hat{u}_x\\
!\Downarrow ~~~~~~~~~~~~~~~~~~~~~~&& ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \Downarrow\\
!\vec{\nabla}\cdot \vec{B} \longleftarrow TE,E_y,P, B^2 \longleftarrow u_x, u_y, P &&  \Downarrow IFFT ~~~~~~~~~~~~~~~~~~~~~~~ \omega \rightarrow C^0, C^1, C^2, C^3  \\
!&& \\
!&& \rho, \rho u_x, \rho u_y, E, B_x, B_y \\
!&& \downarrow\\
!&& u_x = \frac{\rho u_x}{\rho}, u_y = \frac{\rho u_y}{\rho}\\
!&& B^2 = B_x^2 + B_y^2\\
!&& \\
!&& \Downarrow FFT\\
!&& \\
!&& \hat{u}_x, \widehat{u}_y\\
!&& \downarrow\\
!&& i k_x \hat{u}_x, i k_y \hat{u}_y\\
!&& i k_y \hat{B}_x, i k_x \hat{B}_y\\
!&& \\
!&& \Downarrow IFFT\\
!&& \\
!&& \frac{du_x}{dx}, \frac{du_y}{dy}\\
!&& \frac{dB_x}{dy}, \frac{dB_y}{dx}\\
!&& \downarrow\\
!&& CurlB = \frac{dB_y}{dx} - \frac{dB_x}{dy}\\
!&& P = \left(\gamma-1\right)\left[E - \frac{1}{2}\left(\rho u_x u_x + \rho u_y u_y + B^2\right)\right]\\
!&& \downarrow\\
!\end{eqnarray*}
!\newpage
!\begin{eqnarray*}
!&& \downarrow\\
!&& Mom_x^1 = \rho u_x u_x + P + \frac{B^2}{2} - B_x B_x\\
!&& Mom_x^2 = \rho u_x u_y - B_x B_y\\
!&& Mom_y^1 = \rho u_x u_y - B_x B_y\\
!&& Mom_y^2 = \rho u_y u_y + P + \frac{B^2}{2} - B_y B_y\\
!&& \downarrow\\
!&& Energy_x = \left(E+P+\frac{B^2}{2}\right)u_x - u_x B_x B_x - u_y B_x B_y - \eta ~ B_y ~ CurlB\\
!&& Energy_y = \left(E+P+\frac{B^2}{2}\right)u_y - u_x B_x B_y - u_y B_y B_y + \eta ~ B_x ~ CurlB\\
!&& \downarrow\\
!&& E_{Visc} = \left(\frac{du_x}{dx}+\frac{du_y}{dy}\right)^2\\
!&& \downarrow\\
!&& Mag_x = u_y B_x - B_y u_x\\
!&& Mag_y = u_x B_y - B_x u_y\\
!&& \\
!&& \Downarrow FFT\\
!&& \\
!&& i k_x \widehat{\rho u_x}, i k_y \widehat{\rho u_y}\\
!&& i k_x \widehat{Mom}_x^1, i k_y \widehat{Mom}_x^2, i k_x \widehat{Mom}_y^1, i k_y \widehat{Mom}_y^2\\
!&& k_x^2 \hat{u}_x, k_y^2 \hat{u}_x, k_x^2 \hat{u}_y, k_y^2 \hat{u}_y\\
!&& i k_x \widehat{Energy}_x, i k_y \widehat{Energy}_y\\
!&& i k_y \widehat{Mag}_x, i k_x \widehat{Mag}_y\\
!&& k_x^2 \hat{B}_x, k_y^2 \hat{B}_x, k_x^2 \hat{B}_y, k_y^2 \hat{B}_y\\
!&& \downarrow\\
!\end{eqnarray*}
!\newpage
!\begin{eqnarray*}
!&& \downarrow\\
!&& \frac{d\hat{\rho}}{dt} = - \left(i k_x \widehat{\rho u_x} + i k_y \widehat{\rho u_y}\right)\\
!&& \\
!&& \frac{d\widehat{\rho u_x}}{dt} = - \left(i k_x \widehat{Mom}_x^1 + i k_y \widehat{Mom}_x^2 \right) - \mu \left(k_x^2 \hat{u}_x + k_y^2 \hat{u}_x \right)\\
!&& \frac{d\widehat{\rho u_y}}{dt} = - \left(i k_x \widehat{Mom}_y^1 + i k_y \widehat{Mom}_y^2 \right) - \mu \left(k_x^2 \hat{u}_y + k_y^2 \hat{u}_y \right)\\
!&& \\
!&& \frac{d\hat{E}}{dt} = - \left(i k_x \widehat{Energy}_x + i k_y \widehat{Energy}_y\right) + \mu \widehat{E_{Visc}}\\
!&& \\
!&& \frac{d\hat{B}_x}{dt} = - i k_y \widehat{Mag}_x - \eta \left(k_x^2 \hat{B}_x + k_y^2 \hat{B}_x\right)\\
!&& \frac{d\hat{B}_y}{dt} = - i k_x \widehat{Mag}_y - \eta \left(k_x^2 \hat{B}_y + k_y^2 \hat{B}_y\right)\\
!&& \downarrow\\
!&& Adams ~ Bashforth\\
!&& \downarrow\\
!&& \hat{\rho}, \widehat{\rho u_x}, \widehat{\rho u_y}, \hat{E}, \hat{B}_x, \hat{B}_y\\
!&& \\
!&& [At~time~ (t+dt)]\\
!\end{eqnarray*}
!
!\end{document}
!___________________________________________________________________________________________________________________________________________
!
! 
use omp_lib
use cufft3D
use openacc
implicit none
double precision rand !added by Naga to avoid the following compilation error:PGF90-S-0038-Symbol, rand, has not been explicitly declared (PSMHD3.f95)
! Define Grid Size.
integer ( kind = 4 ), parameter :: Np = 10000
integer ( kind = 4 ), parameter :: Nx = 256
integer ( kind = 4 ), parameter :: Ny = 256 
integer ( kind = 4 ), parameter :: Nz = 256 
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

integer ( kind = 4 ) i,j,k,t,m,n,q
real ( kind = 8 ) Lx,Ly,Lz,dx,dy,dz,kx,ky,kz,time,time_min,time_max,dt,G1,G2,G3,h,sx,sy,sz,tx,ty,tz
real ( kind = 8 ) spheat,rho0,U0,mu,kf,mode,sigma,uy0,MS,CS,Re,Rm,PM,P0,eta,mu_0,MA,VA,B0,ba,A,B,C
real ( kind = 8 ) Pressure,Energy,Energy0,T_Energy,y_Energy0,B_Field,C0,C1,C2,C3,div_B_tot,Hf,Hm,HT,HuB,HAw,Rayleigh

!real ( kind = 8 ) x(Nx),y(Ny),z(Nz),rho(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
real ( kind = 8 ), dimension (:), allocatable :: x,y,z
real (kind = 8), dimension (:,:,:), allocatable :: rho,ux,uy,uz

!real ( kind = 8 ) rho_ux(Nx,Ny,Nz),rho_uy(Nx,Ny,Nz),rho_uz(Nx,Ny,Nz)
real ( kind = 8 ), dimension (:,:,:), allocatable :: rho_ux,rho_uy,rho_uz

!real ( kind = 8 ) omega_x(Nx,Ny,Nz),omega_y(Nx,Ny,Nz),omega_z(Nx,Ny,Nz),omega2(Nx,Ny,Nz)
real ( kind = 8 ), dimension (:,:,:), allocatable :: omega_x,omega_y,omega_z,omega2

!real ( kind = 8 ) P(Nx,Ny,Nz),E(Nx,Ny,Nz),Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz),div_B(Nx,Ny,Nz)
real ( kind = 8 ), dimension (:,:,:), allocatable :: P,E,Bx,By,Bz,B2,div_B,jx,jy,jz,Ax,Ay,Az,j2,A2,u2

real ( kind = 8 ), dimension (:), allocatable ::  xp,yp,zp,uxp,uyp,uzp,Bpx,Bpy,Bpz,vxp,vyp,vzp,vxd,vyd,vzd

!Backup variables for FFTW have been commented
!real ( kind = 8 ) rho_dum(Nx,Ny,Nz),ux_dum(Nx,Ny,Nz),uy_dum(Nx,Ny,Nz),uz_dum(Nx,Ny,Nz)
!real ( kind = 8 ) rho_ux_dum(Nx,Ny,Nz),rho_uy_dum(Nx,Ny,Nz),rho_uz_dum(Nx,Ny,Nz)
!real ( kind = 8 ) P_dum(Nx,Ny,Nz),E_dum(Nx,Ny,Nz),Bx_dum(Nx,Ny,Nz),By_dum(Nx,Ny,Nz),Bz_dum(Nx,Ny,Nz)

!complex ( kind = 8 ) rho_k(Nh,Ny,Nz),ux_k(Nh,Ny,Nz),uy_k(Nh,Ny,Nz),uz_k(Nh,Ny,Nz)
!complex ( kind = 8 ), dimension (:,:,:), allocatable :: rho_k,ux_k,uy_k,uz_k

!complex ( kind = 8 ) rho_ux_k(Nh,Ny,Nz),rho_uy_k(Nh,Ny,Nz),rho_uz_k(Nh,Ny,Nz)
!complex ( kind = 8 ), dimension (:,:,:), allocatable :: rho_ux_k,rho_uy_k,rho_uz_k

!complex ( kind = 8 ) omega_x_k(Nh,Ny,Nz),omega_y_k(Nh,Ny,Nz),omega_z_k(Nh,Ny,Nz)
!complex ( kind = 8 ), dimension (:,:,:), allocatable :: omega_x_k,omega_y_k,omega_z_k

!complex ( kind = 8 ) P_k(Nh,Ny,Nz),E_k(Nh,Ny,Nz),Ek(Nh,Ny,Nz),Bk(Nh,Ny,Nz)
!complex ( kind = 8 ), dimension (:,:,:), allocatable :: P_k,E_k,Ek,Bk
complex ( kind = 8 ), dimension (:,:,:), allocatable :: Ek,Bk

!complex ( kind = 8 ) Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz),div_B_k(Nh,Ny,Nz)
!complex ( kind = 8 ), dimension (:,:,:), allocatable :: Bx_k,By_k,Bz_k,div_B_k

!complex ( kind = 8 ) rho_k_dum(Nh,Ny,Nz),ux_k_dum(Nh,Ny,Nz),uy_k_dum(Nh,Ny,Nz),uz_k_dum(Nh,Ny,Nz)
!complex ( kind = 8 ) omega_x_k_dum(Nh,Ny,Nz),omega_y_k_dum(Nh,Ny,Nz),omega_z_k_dum(Nh,Ny,Nz)
!complex ( kind = 8 ) rho_ux_k_dum(Nh,Ny,Nz),rho_uy_k_dum(Nh,Ny,Nz),rho_uz_k_dum(Nh,Ny,Nz)
!complex ( kind = 8 ) E_k_dum(Nh,Ny,Nz),Bx_k_dum(Nh,Ny,Nz),By_k_dum(Nh,Ny,Nz),Bz_k_dum(Nh,Ny,Nz)

!complex ( kind = 8 ) rho_k_new(Nh,Ny,Nz),rho_ux_k_new(Nh,Ny,Nz),rho_uy_k_new(Nh,Ny,Nz),rho_uz_k_new(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new

!complex ( kind = 8 ) E_k_new(Nh,Ny,Nz),Bx_k_new(Nh,Ny,Nz),By_k_new(Nh,Ny,Nz),Bz_k_new(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: E_k_new,Bx_k_new,By_k_new,Bz_k_new

!complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny,Nz),d_rho_ux_k_dt_old(Nh,Ny,Nz),d_rho_uy_k_dt_old(Nh,Ny,Nz),d_rho_uz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old

!complex ( kind = 8 ) d_E_k_dt_old(Nh,Ny,Nz),d_Bx_k_dt_old(Nh,Ny,Nz),d_By_k_dt_old(Nh,Ny,Nz),d_Bz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old

!complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny,Nz),d_rho_ux_k_dt_new(Nh,Ny,Nz),d_rho_uy_k_dt_new(Nh,Ny,Nz),d_rho_uz_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new

!complex ( kind = 8 ) d_E_k_dt_new(Nh,Ny,Nz),d_Bx_k_dt_new(Nh,Ny,Nz),d_By_k_dt_new(Nh,Ny,Nz),d_Bz_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new

!Variables used in derive
real ( kind = 8 ), dimension(:,:,:), allocatable :: Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3

real ( kind = 8 ), dimension(:,:,:), allocatable :: d_ux_dx,d_uy_dy,d_uz_dz,Fx,Fy,Fz
real ( kind = 8 ), dimension(:,:,:), allocatable :: Energy_x,Energy_y,Energy_z,E_Visc
real ( kind = 8 ), dimension(:,:,:), allocatable :: Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2
real ( kind = 8 ), dimension(:,:,:), allocatable :: d_Bx_dy,d_By_dx,d_Bx_dz
real ( kind = 8 ), dimension(:,:,:), allocatable :: d_By_dz,d_Bz_dx,d_Bz_dy
real ( kind = 8 ) :: curl_x_B,curl_y_B,curl_z_B
complex ( kind = 8 ) :: i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k
complex ( kind = 8 ) :: i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k
complex ( kind = 8 ) :: i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k
complex ( kind = 8 ) :: i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k
complex ( kind = 8 ) :: kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k
complex ( kind = 8 ) :: i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k
complex ( kind = 8 ) :: kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k

!End variables used in derive


!cufft
integer (acc_handle_kind) :: stream
type (c_ptr) :: fftPlanD2ZMain, fftPlanZ2DMain

integer ( kind = 4 ) thread_num,num_thread,proc_num ! OMP
real ( kind = 8 ) t1,t2

common/comm/time_max,Lx,Ly,Lz,dx,dy,dz,spheat,mu,ms,CS,mu_0,eta,dt,kf

integer,parameter :: seed = 99999999
call srand(seed)

!Do all the allocations
allocate(x(Nx))
allocate(y(Ny))
allocate(z(Nz))

allocate(xp(Np))
allocate(yp(Np))
allocate(zp(Np))
allocate(uxp(Np))
allocate(uyp(Np))
allocate(uzp(Np))
allocate(Bpx(Np))
allocate(Bpy(Np))
allocate(Bpz(Np))
allocate(vxp(Np))
allocate(vyp(Np))
allocate(vzp(Np))
allocate(vxd(Np))
allocate(vyd(Np))
allocate(vzd(Np))

allocate(rho(Nx+2, Ny, Nz))
allocate(ux(Nx+2, Ny, Nz))
allocate(uy(Nx+2, Ny, Nz))
allocate(uz(Nx+2, Ny, Nz))

allocate(rho_ux(Nx+2, Ny, Nz))
allocate(rho_uy(Nx+2, Ny, Nz))
allocate(rho_uz(Nx+2, Ny, Nz))

allocate(omega_x(Nx+2, Ny, Nz))
allocate(omega_y(Nx+2, Ny, Nz))
allocate(omega_z(Nx+2, Ny, Nz))
allocate(omega2(Nx+2, Ny, Nz))

allocate(P(Nx+2, Ny, Nz))
allocate(E(Nx+2, Ny, Nz))
allocate(Bx(Nx+2, Ny, Nz))
allocate(By(Nx+2, Ny, Nz))
allocate(Bz(Nx+2, Ny, Nz))
allocate(B2(Nx+2, Ny, Nz))
allocate(div_B(Nx+2, Ny, Nz))
allocate(jx(Nx+2, Ny, Nz))
allocate(jy(Nx+2, Ny, Nz))
allocate(jz(Nx+2, Ny, Nz))
allocate(j2(Nx+2, Ny, Nz))
allocate(Ax(Nx+2, Ny, Nz))
allocate(Ay(Nx+2, Ny, Nz))
allocate(Az(Nx+2, Ny, Nz))
allocate(A2(Nx+2, Ny, Nz))
allocate(u2(Nx+2, Ny, Nz))


!allocate(rho_k(Nh, Ny, Nz))
!allocate(ux_k(Nh, Ny, Nz))
!allocate(uy_k(Nh, Ny, Nz))
!allocate(uz_k(Nh, Ny, Nz))

!allocate(rho_ux_k(Nh, Ny, Nz))
!allocate(rho_uy_k(Nh, Ny, Nz))
!allocate(rho_uz_k(Nh, Ny, Nz))

!allocate(omega_x_k(Nh, Ny, Nz))
!allocate(omega_y_k(Nh, Ny, Nz))
!allocate(omega_z_k(Nh, Ny, Nz))

!allocate(P_k(Nh, Ny, Nz))
!allocate(E_k(Nh, Ny, Nz))
allocate(Ek(Nh, Ny, Nz))
allocate(Bk(Nh, Ny, Nz))

!allocate(Bx_k(Nh, Ny, Nz))
!allocate(By_k(Nh, Ny, Nz))
!allocate(Bz_k(Nh, Ny, Nz))
!allocate(div_B_k(Nh, Ny, Nz))

allocate(rho_k_new(Nh, Ny, Nz))
allocate(rho_ux_k_new(Nh, Ny, Nz))
allocate(rho_uy_k_new(Nh, Ny, Nz))
allocate(rho_uz_k_new(Nh, Ny, Nz))

allocate(E_k_new(Nh, Ny, Nz))
allocate(Bx_k_new(Nh, Ny, Nz))
allocate(By_k_new(Nh, Ny, Nz))
allocate(Bz_k_new(Nh, Ny, Nz))

allocate(d_rho_k_dt_old(Nh, Ny, Nz))
allocate(d_rho_ux_k_dt_old(Nh, Ny, Nz))
allocate(d_rho_uy_k_dt_old(Nh, Ny, Nz))
allocate(d_rho_uz_k_dt_old(Nh, Ny, Nz))

allocate(d_E_k_dt_old(Nh, Ny, Nz))
allocate(d_Bx_k_dt_old(Nh, Ny, Nz))
allocate(d_By_k_dt_old(Nh, Ny, Nz))
allocate(d_Bz_k_dt_old(Nh, Ny, Nz))

allocate(d_rho_k_dt_new(Nh, Ny, Nz))
allocate(d_rho_ux_k_dt_new(Nh, Ny, Nz))
allocate(d_rho_uy_k_dt_new(Nh, Ny, Nz))
allocate(d_rho_uz_k_dt_new(Nh, Ny, Nz))

allocate(d_E_k_dt_new(Nh, Ny, Nz))
allocate(d_Bx_k_dt_new(Nh, Ny, Nz))
allocate(d_By_k_dt_new(Nh, Ny, Nz))
allocate(d_Bz_k_dt_new(Nh, Ny, Nz))

!Variables used in derive
allocate(Mom_x_1(Nx+2,Ny,Nz))
allocate(Mom_x_2(Nx+2,Ny,Nz))
allocate(Mom_x_3(Nx+2,Ny,Nz))
allocate(Mom_y_1(Nx+2,Ny,Nz))
allocate(Mom_y_2(Nx+2,Ny,Nz))
allocate(Mom_y_3(Nx+2,Ny,Nz))
allocate(Mom_z_1(Nx+2,Ny,Nz))
allocate(Mom_z_2(Nx+2,Ny,Nz))
allocate(Mom_z_3(Nx+2,Ny,Nz))

allocate(d_ux_dx(Nx+2,Ny,Nz))
allocate(d_uy_dy(Nx+2,Ny,Nz))
allocate(d_uz_dz(Nx+2,Ny,Nz))
allocate(Fx(Nx+2,Ny,Nz))
allocate(Fy(Nx+2,Ny,Nz))
allocate(Fz(Nx+2,Ny,Nz))

allocate(Energy_x(Nx+2,Ny,Nz))
allocate(Energy_y(Nx+2,Ny,Nz))
allocate(Energy_z(Nx+2,Ny,Nz))
allocate(E_Visc(Nx+2,Ny,Nz))

allocate(Mag_x_1(Nx+2,Ny,Nz))
allocate(Mag_x_2(Nx+2,Ny,Nz))
allocate(Mag_y_1(Nx+2,Ny,Nz))
allocate(Mag_y_2(Nx+2,Ny,Nz))
allocate(Mag_z_1(Nx+2,Ny,Nz))
allocate(Mag_z_2(Nx+2,Ny,Nz))

allocate(d_Bx_dy(Nx+2,Ny,Nz))
allocate(d_By_dx(Nx+2,Ny,Nz))
allocate(d_Bx_dz(Nx+2,Ny,Nz))
allocate(d_By_dz(Nx+2,Ny,Nz))
allocate(d_Bz_dx(Nx+2,Ny,Nz))
allocate(d_Bz_dy(Nx+2,Ny,Nz))

!allocate(curl_x_B(Nx+2,Ny,Nz))
!allocate(curl_y_B(Nx+2,Ny,Nz))
!allocate(curl_z_B(Nx+2,Ny,Nz))

!allocate(i_kx_rho_ux_k(Nh,Ny,Nz))
!allocate(i_ky_rho_uy_k(Nh,Ny,Nz))
!allocate(i_kz_rho_uz_k(Nh,Ny,Nz))

!allocate(i_kx_Mom_x_1_k(Nh,Ny,Nz))
!allocate(i_ky_Mom_x_2_k(Nh,Ny,Nz))
!allocate(i_kz_Mom_x_3_k(Nh,Ny,Nz))
!allocate(i_kx_Mom_y_1_k(Nh,Ny,Nz))
!allocate(i_ky_Mom_y_2_k(Nh,Ny,Nz))
!allocate(i_kz_Mom_y_3_k(Nh,Ny,Nz))
!allocate(i_kx_Mom_z_1_k(Nh,Ny,Nz))
!allocate(i_ky_Mom_z_2_k(Nh,Ny,Nz))
!allocate(i_kz_Mom_z_3_k(Nh,Ny,Nz))

!allocate(kx2_ux_k(Nh,Ny,Nz))
!allocate(ky2_ux_k(Nh,Ny,Nz))
!allocate(kz2_ux_k(Nh,Ny,Nz))
!allocate(kx2_uy_k(Nh,Ny,Nz))
!allocate(ky2_uy_k(Nh,Ny,Nz))
!allocate(kz2_uy_k(Nh,Ny,Nz))
!allocate(kx2_uz_k(Nh,Ny,Nz))
!allocate(ky2_uz_k(Nh,Ny,Nz))
!allocate(kz2_uz_k(Nh,Ny,Nz))

!allocate(i_kx_Energy_x_k(Nh,Ny,Nz))
!allocate(i_ky_Energy_y_k(Nh,Ny,Nz))
!allocate(i_kz_Energy_z_k(Nh,Ny,Nz))
!allocate(i_ky_Mag_x_1_k(Nh,Ny,Nz))
!allocate(i_kz_Mag_x_2_k(Nh,Ny,Nz))
!allocate(i_kx_Mag_y_1_k(Nh,Ny,Nz))
!allocate(i_kz_Mag_y_2_k(Nh,Ny,Nz))
!allocate(i_kx_Mag_z_1_k(Nh,Ny,Nz))
!allocate(i_ky_Mag_z_2_k(Nh,Ny,Nz))
!
!allocate(kx2_Bx_k(Nh,Ny,Nz))
!allocate(ky2_Bx_k(Nh,Ny,Nz))
!allocate(kz2_Bx_k(Nh,Ny,Nz))
!allocate(kx2_By_k(Nh,Ny,Nz))
!allocate(ky2_By_k(Nh,Ny,Nz))
!allocate(kz2_By_k(Nh,Ny,Nz))
!allocate(kx2_Bz_k(Nh,Ny,Nz))
!allocate(ky2_Bz_k(Nh,Ny,Nz))
!allocate(kz2_Bz_k(Nh,Ny,Nz))

!End variables used in derive


!===================== FILENAMES ==============================================	

open(unit=5,file='System_information.dat',status='unknown')
open(unit=100,file='omegaval-inplace.dat',status='unknown')
open(unit=15,file='Initial_Grid_Data.dat',status='unknown')
!open(unit=20,file='INPUT.dat',status='old')  ! This input file is the file generated from Vorticity code. 
!open(unit=25,file='Initial_Grid_Data_Reproduced.dat',status='unknown')
!open(unit=30,file='Initial_Energy_Spectra.dat',status='unknown')
open(unit=35,file='Energy_Spectra.dat',status='unknown')
open(unit=40,file='Energy.dat',status='unknown')


!===================== USER INPUTS ============================================		

! Define Number of Threads.
!proc_num = omp_get_num_procs()
!thread_num = 25
!call omp_set_num_threads (thread_num)

!write(5,*) "Number of Processors=", proc_num, "Number of Threads=", thread_num  

! System Size.
Lx = 2.0d0*pi; Ly = 2.0d0*pi; Lz = 2.0d0*pi 
!Lx = 1.0d0; Ly = 1.0d0; Lz = 1.0d0

! Grid Resolution.
dx = Lx/dfloat(Nx); dy = Ly/dfloat(Ny); dz = Lz/dfloat(Nz)

! Runtime Details and Time Step Width.
time_min = 0.0d0
!time_max = 100.0d0
!time_max = 1.0d0
time_max = 0.00001d0
dt = 0.000010d0    

! Ratio of Specific Heats/ Adiabetic Index/ Heat Capacity Ratio/ Poisson Constant.
spheat = 5.0d0/3.0d0

! Co-efficient of Viscosity.
mu = 1000.0d0

! Co-efficient of Resistivity.
eta = 1000.0d0

! Magnetic Permeability.
mu_0 = 1.0d0

! Background Density.
rho0 = 1.0d0

! Initial Pressure.
P0 = 1.0d0

! Mach Number.
MS = 0.10d0

! Alfven Mach Number.
MA = 1.0d0

! Maximum Velocity.
!U0 = 0.14732247502913884d0
U0 = 0.10d0 

! Sound Speed.
! CS = dsqrt(spheat*P0/rho0)
! U0 = M * CS

! Sound Speed.
 CS = 1.0d0/MS

! Alfven Speed.
VA = U0/MA

! Initial Magnetic Field.
B0 = VA*dsqrt(mu_0*rho0) 

write(5,*) "Sound Speed=", CS, "Initial Velocity=", U0
write(5,*) "Alfven Speed=", VA, "Initial Magnetic Field=", B0

! Forcing Length Scale.
kf = 1.0d0

 A = 1.0d0
 B = 1.0d0
 C = 1.0d0

! Raynold's Number.
Re = ms*rho0*U0*Lx/mu; write(5,*) "Raynold's Number =", Re
Rm = U0*Lx/eta; write(5,*) "Magnetic Raynold's Number =", Rm
PM = mu/eta; write(5,*) "Prandtl Number =", PM

do i = 1, Np
xp(i) = rand(0)*Lx
yp(i) = rand(0)*Ly
zp(i) = rand(0)*Lz
uxp(i) = 0.0d0
uyp(i) = 0.0d0
uzp(i) = 0.0d0
end do

! Grid Generation.
do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      x(i)=0.0d0+real(i-1)*dx
      y(j)=0.0d0+real(j-1)*dy
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
      E(i,j,k) = P(i,j,k)/(spheat -1.0d0) + 0.50d0*rho(i,j,k) * (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)&
                 + 0.50d0*(Bx(i,j,k)**2 + By(i,j,k)**2 + Bz(i,j,k)**2)
      ! Initial Combined Variables.
      rho_ux(i,j,k) = rho(i,j,k) * ux(i,j,k)
      rho_uy(i,j,k) = rho(i,j,k) * uy(i,j,k)
      rho_uz(i,j,k) = rho(i,j,k) * uz(i,j,k)
      ! Keep Backup of the Arrays for FFTW. 
      !rho_dum(i,j,k) = rho(i,j,k)
      !ux_dum(i,j,k) = ux(i,j,k)
      !uy_dum(i,j,k) = uy(i,j,k)
      !uz_dum(i,j,k) = uz(i,j,k)
      !rho_ux_dum(i,j,k) = rho_ux(i,j,k)
      !rho_uy_dum(i,j,k) = rho_uy(i,j,k)
      !rho_uz_dum(i,j,k) = rho_uz(i,j,k)
      !P_dum(i,j,k) = P(i,j,k)
      !E_dum(i,j,k) = E(i,j,k) 
      !Bx_dum(i,j,k) = Bx(i,j,k)
      !By_dum(i,j,k) = By(i,j,k)
      !Bz_dum(i,j,k) = Bz(i,j,k)
      ! Write Initial Density and Velocity Distribution in File.
      !write(15,*) x(i),y(j),z(k),rho(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k),P(i,j,k),E(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    enddo
  enddo
enddo

 close(15)
write(*,*) "After initial grid generation"
!===================== INITIAL TIME DATA ===============================================
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_dum, rho_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_ux_dum, rho_ux_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_uy_dum, rho_uy_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_uz_dum, rho_uz_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, P_dum, P_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, E_dum, E_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Bx_dum, Bx_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, By_dum, By_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Bz_dum, Bz_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)

t1 = omp_get_wtime()

  call acc_set_device_num(0, 0)
  stream = acc_get_cuda_stream(acc_async_sync)
  call createCUFFTPlan3D(fftPlanD2ZMain, Nz, Ny, Nx, CUFFT_D2Z, stream) !indices are swapped as fortran is column major and c is row major
  call createCUFFTPlan3D(fftPlanZ2DMain, Nz, Ny, Nx, CUFFT_Z2D, stream) !indices are swapped as fortran is column major and c is row major
  
! !$acc data copy(rho, ux, uy, uz, rho_ux, rho_uy, rho_uz, P, E, Bx, By, Bz) & 
! !$acc & copyout(omega_x, omega_y, omega_z) &
!!$acc & create(Fx, Fy, Fz) &
!!$acc & create(d_ux_dx, d_uy_dy, d_uz_dz, d_Bx_dy, d_Bx_dz, d_By_dx, d_By_dz, d_Bz_dx, d_Bz_dy) &
!!$acc & create(Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3, Energy_x,Energy_y,Energy_z,E_Visc, Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2) &
!!$acc & create(x, y, z) &
!!$acc & copy(rho_k_new, rho_ux_k_new, rho_uy_k_new, rho_uz_k_new, E_k_new, Bx_k_new, By_k_new, Bz_k_new) &
!!$acc & copy(d_rho_k_dt_old, d_rho_ux_k_dt_old, d_rho_uy_k_dt_old, d_rho_uz_k_dt_old, d_E_k_dt_old, d_Bx_k_dt_old, d_By_k_dt_old, d_Bz_k_dt_old) &
!!$acc & copy(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new) &
!!$acc & copy(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new) copyout(P, div_B,B2,ux,uy,uz,omega2,Ek,Bk) 

 ! Inside this region the device data pointer will be used
 !$acc host_data use_device(rho,  ux, uy, uz, rho_ux, rho_uy, rho_uz, P, E, Bx, By, Bz)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho), C_LOC(rho), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(ux), C_LOC(ux), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(uy), C_LOC(uy), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(uz), C_LOC(uz), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho_ux), C_LOC(rho_ux), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho_uy), C_LOC(rho_uy), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho_uz), C_LOC(rho_uz), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(P), C_LOC(P), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(E), C_LOC(E), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Bx), C_LOC(Bx), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(By), C_LOC(By), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Bz), C_LOC(Bz), CUFFT_D2Z)

 !$acc end host_data
! Evaluate Initial Vorticity Spectra.
!$acc parallel firstprivate(Lx, Ly, Lz) present(ux, uy, uz, omega_x, omega_y, omega_z)
!$acc loop collapse(3)
do k = 1, Nz
  do j = 1, Ny 
    do i = 1, Nx/2+1
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      if (j <= Ny/2) then
        ky = 2.0d0*pi*dfloat(j-1)/Ly
      else
        ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      endif
      if (k <= Nz/2) then
        kz = 2.0d0*pi*dfloat(k-1)/Lz
      else
        kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      endif

      omega_x(2*i-1,j,k) = (-1.0d0)*ky*uz(2*i,j,k) - (-1.0d0)*kz*uy(2*i,j,k)
      omega_x(2*i,j,k) = (1.0d0)*ky*uz(2*i-1,j,k) - (1.0d0)*kz*uy(2*i-1,j,k)
      omega_y(2*i-1,j,k) = (-1.0d0)*kx*uz(2*i,j,k) - (-1.0d0)*kz*ux(2*i,j,k)
      omega_y(2*i,j,k) = (1.0d0)*kx*uz(2*i-1,j,k) - (1.0d0)*kz*ux(2*i-1,j,k)
      omega_z(2*i-1,j,k) = (-1.0d0)*kx*uy(2*i,j,k) - (-1.0d0)*ky*ux(2*i,j,k)
      omega_z(2*i,j,k) = (1.0d0)*kx*uy(2*i-1,j,k) - (1.0d0)*ky*ux(2*i-1,j,k)
      
!      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
!      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
!      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo
!$acc end parallel
 !$acc host_data use_device(omega_x, omega_y, omega_z)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_x_k_dum, omega_x, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(omega_x), C_LOC(omega_x), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_y_k_dum, omega_y, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(omega_y), C_LOC(omega_y), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_z_k_dum, omega_z, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(omega_z), C_LOC(omega_z), CUFFT_Z2D)

 ! call dfftw_destroy_plan_ (plan_forward)
 ! call dfftw_destroy_plan_ (plan_backward)

 !$acc end host_data
! Evaluate Initial Energy Spectra.
!do i = 1, Nx/2+1
  !do j = 1, Ny
    !do k = 1, Nz
      !write(30,*) i-1,j-1,k-1,abs(nk(i,j,k)),abs(ux_k(i,j,k)),abs(uy_k(i,j,k)),abs(uz_k(i,j,k)) 
    !enddo
  !enddo
!enddo

!Energy0 = 0.0d0; y_Energy0 = 0.0d0

!$acc parallel present(omega_x, omega_y, omega_z) 
!$acc loop collapse(3) 
do k = 1, Nz
  do j = 1, Ny
    do i = 1, Nx
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
    enddo
  enddo
enddo
!$acc end parallel
!!$acc end data
!NAGA - note: omega_x, omega_y and omega_z are not copied back as they dont seem to be used
! Write the Initial Energy in File.
!write(40,*) 0.0d0, Energy0,y_Energy0

!======================= MAIN PROGRAM =====================================================
write(*,*) "Before loop"
! Time Loop Starts.
do time = time_min,time_max,dt

write(*,*) "Iteration"
t = nint(time/dt) - dint(time_min/dt)

! Online Check of the Progress.
!if (mod(t,int(1.0/dt)) == 0) then
  !print*, t/int(1.0/dt)
!endif

!======================================================================================

!write(*,*) "Before derive"

!!$acc data copy(rho, rho_ux, rho_uy, rho_uz, Bx, By, Bz) copy(E) create(Fx, Fy, Fz) &
!!$acc & create(d_ux_dx, d_uy_dy, d_uz_dz, d_Bx_dy, d_Bx_dz, d_By_dx, d_By_dz, d_Bz_dx, d_Bz_dy) &
!!$acc & create(curl_x_B,curl_y_B,curl_z_B,Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3, Energy_x,Energy_y,Energy_z,E_Visc, Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2) &
!!$acc & create(i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k, i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k) &
!!$acc & create(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k,i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k) &
!!$acc & create(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k) &
!!$acc & create(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k, kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k) &
!!$acc & create(x, y, z) &
!!$acc & copy(rho_k_new, rho_ux_k_new, rho_uy_k_new, rho_uz_k_new, E_k_new, Bx_k_new, By_k_new, Bz_k_new) &
!!$acc & copy(d_rho_k_dt_old, d_rho_ux_k_dt_old, d_rho_uy_k_dt_old, d_rho_uz_k_dt_old, d_E_k_dt_old, d_Bx_k_dt_old, d_By_k_dt_old, d_Bz_k_dt_old) &
!!$acc & copy(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new) &
!!$acc & copy(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new) copyout(P, div_B,B2,ux,uy,uz,omega2,Ek,Bk) 

!$acc host_data use_device(rho, rho_ux, rho_uy, rho_uz, E, Bx, By, Bz)
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_k_dum, rho, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(rho), C_LOC(rho), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_ux_k_dum, rho_ux, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(rho_ux), C_LOC(rho_ux), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uy_k_dum, rho_uy, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(rho_uy), C_LOC(rho_uy), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uz_k_dum, rho_uz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(rho_uz), C_LOC(rho_uz), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, E_k_dum, E, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(E), C_LOC(E), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bx_k_dum, Bx, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(Bx), C_LOC(Bx), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, By_k_dum, By, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(By), C_LOC(By), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bz_k_dum, Bz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(Bz), C_LOC(Bz), CUFFT_Z2D)

!$acc end host_data

write(*,*) "Derive 1"
!if (mod(dfloat(t),0.0010d0/dt)==0.0d0) then 
 A = 1.0d0
 B = 1.0d0
 C = 1.0d0
!else
! A = 0.0d0
! B = 0.0d0
! C = 0.0d0
!endif

!!$OMP PARALLEL SHARED(dx,dy,dz,A,B,C,kf,x,y,z,rho,ux,uy,uz),&
!!$OMP & SHARED(rho_ux,rho_uy,rho_uz,Fx,Fy,Fz,E,Bx,By,Bz,B2) PRIVATE(i,j,k)
!!$OMP DO
!$acc parallel firstprivate(dx, dy, dz) present(rho, rho_ux, rho_uy, rho_uz, E, Bx, By, Bz, ux, uy, uz, Fx, Fy, Fz, B2, x, y, z) 
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx
    x(i)=0.0d0+real(i-1)*dx
    y(j)=0.0d0+real(j-1)*dy
    z(k)=0.0d0+real(k-1)*dz 
    ! FFTW Normalisation
    rho(i,j,k) = rho(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
    rho_ux(i,j,k) = rho_ux(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    rho_uy(i,j,k) = rho_uy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz)) 
    rho_uz(i,j,k) = rho_uz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
    E(i,j,k) = E(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
    Bx(i,j,k) = Bx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    By(i,j,k) = By(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    Bz(i,j,k) = Bz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
    ! Evaluate Velocity in Real Space. 
    ux(i,j,k) = rho_ux(i,j,k)/rho(i,j,k)
    uy(i,j,k) = rho_uy(i,j,k)/rho(i,j,k)
    uz(i,j,k) = rho_uz(i,j,k)/rho(i,j,k)

    ! Evaluate Forcing.
    Fx(i,j,k) = A*dsin(kf*z(k)) + C*dcos(kf*y(j))
    Fy(i,j,k) = B*dsin(kf*x(i)) + A*dcos(kf*z(k))
    Fz(i,j,k) = C*dsin(kf*y(j)) + B*dcos(kf*x(i))
  
    ! Evaluate Square of Magnetic Field.
    B2(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)
  
    ! Keep Backup of the Arrays for FFTW.
    !ux_dum(i,j,k) = ux(i,j,k)
    !uy_dum(i,j,k) = uy(i,j,k)
    !uz_dum(i,j,k) = uz(i,j,k)
    enddo
  enddo
enddo   
!$acc end parallel
!!$OMP END DO
!!$OMP END PARALLEL
write(*,*) "Derive 2"

!$acc host_data use_device(ux, uy,uz, Fx, Fy, Fz, Bx,By,Bz)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(ux), C_LOC(ux), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(uy), C_LOC(uy), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(uz), C_LOC(uz), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fx, Fx_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Fx), C_LOC(Fx), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fy, Fy_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Fy), C_LOC(Fy), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fz, Fz_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Fz), C_LOC(Fz), CUFFT_D2Z)
!Recomputing transform on Bx,By,Bz
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Bx), C_LOC(Bx), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(By), C_LOC(By), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Bz), C_LOC(Bz), CUFFT_D2Z)
!$acc end host_data

write(*,*) "Derive 3"
! Evaluate derivatives of Velocity and Magnetic Field.

!!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,Bx_k,By_k),&
!!$OMP & SHARED(i_kx_ux_k,i_ky_uy_k,i_kz_uz_k),&
!!$OMP & SHARED(i_ky_Bx_k,i_kz_Bx_k,i_kx_By_k,i_kz_By_k,i_kx_Bz_k,i_ky_Bz_k) PRIVATE(i,j,k,kx,ky,kz)
!!$OMP DO
!$acc parallel firstprivate(Lx, Ly, Lz) present(d_ux_dx, d_uy_dy, d_uz_dz, d_Bx_dy, d_Bx_dz, d_By_dx, d_By_dz, d_Bz_dx, d_Bz_dy) &
!$acc & present(ux, uy, uz, Bx, By, Bz) 
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx/2+1
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      if (j <= Ny/2) then
        ky = 2.0d0*pi*dfloat(j-1)/Ly
      else
        ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      endif
      if (k <= Nz/2) then
        kz = 2.0d0*pi*dfloat(k-1)/Lz
      else
        kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      endif

      d_ux_dx(2*i-1,j,k) = (-1.0d0)*kx*ux(2*i,j,k)
      d_ux_dx(2*i,j,k) = (1.0d0)*kx*ux(2*i-1,j,k)
      d_uy_dy(2*i-1,j,k) = (-1.0d0)*ky*uy(2*i,j,k)
      d_uy_dy(2*i,j,k) = (1.0d0)*ky*uy(2*i-1,j,k)
      d_uz_dz(2*i-1,j,k) = (-1.0d0)*kz*uz(2*i,j,k)
      d_uz_dz(2*i,j,k) = (1.0d0)*kz*uz(2*i-1,j,k)
      
      d_Bx_dy(2*i-1,j,k) = (-1.0d0)*ky*Bx(2*i,j,k)
      d_Bx_dy(2*i,j,k) = (1.0d0)*ky*Bx(2*i-1,j,k)
      d_Bx_dz(2*i-1,j,k) = (-1.0d0)*kz*Bx(2*i,j,k)
      d_Bx_dz(2*i,j,k) = (1.0d0)*kz*Bx(2*i-1,j,k)
      d_By_dx(2*i-1,j,k) = (-1.0d0)*kx*By(2*i,j,k)
      d_By_dx(2*i,j,k) = (1.0d0)*kx*By(2*i-1,j,k)
      d_By_dz(2*i-1,j,k) = (-1.0d0)*kz*By(2*i,j,k)
      d_By_dz(2*i,j,k) = (1.0d0)*kz*By(2*i-1,j,k)
      d_Bz_dx(2*i-1,j,k) = (-1.0d0)*kx*Bz(2*i,j,k)
      d_Bz_dx(2*i,j,k) = (1.0d0)*kx*Bz(2*i-1,j,k)
      d_Bz_dy(2*i-1,j,k) = (-1.0d0)*ky*Bz(2*i,j,k)
      d_Bz_dy(2*i,j,k) = (1.0d0)*ky*Bz(2*i-1,j,k)

    enddo 
  enddo
enddo
!$acc end parallel
!!$OMP END DO
!!$OMP END PARALLEL
write(*,*) "Derive 4"

!$acc host_data use_device(d_ux_dx, d_uy_dy, d_uz_dz, d_Bx_dy, d_Bx_dz, d_By_dx, d_By_dz, d_Bz_dx, d_Bz_dy, Bx,By,Bz)
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_ux_k, d_ux_dx, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(d_ux_dx), C_LOC(d_ux_dx), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_uy_k, d_uy_dy, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(d_uy_dy), C_LOC(d_uy_dy), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_uz_k, d_uz_dz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(d_uz_dz), C_LOC(d_uz_dz), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_Bx_k, d_Bx_dy, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(d_Bx_dy), C_LOC(d_Bx_dy), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_Bx_k, d_Bx_dz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(d_Bx_dz), C_LOC(d_Bx_dz), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_By_k, d_By_dx, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(d_By_dx), C_LOC(d_By_dx), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_By_k, d_By_dz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(d_By_dz), C_LOC(d_By_dz), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_Bz_k, d_Bz_dx, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(d_Bz_dx), C_LOC(d_Bz_dx), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_Bz_k, d_Bz_dy, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(d_Bz_dy), C_LOC(d_Bz_dy), CUFFT_Z2D)

!Recompute inverse transform of Bx,By,Bz
   call executeCUFFT3D(fftPlanZ2DMain, C_LOC(Bx), C_LOC(Bx), CUFFT_Z2D)
   call executeCUFFT3D(fftPlanZ2DMain, C_LOC(By), C_LOC(By), CUFFT_Z2D)
   call executeCUFFT3D(fftPlanZ2DMain, C_LOC(Bz), C_LOC(Bz), CUFFT_Z2D)

!$acc end host_data
write(*,*) "Derive 5"
!!$OMP PARALLEL SHARED(spheat,eta,ux,uy,uz,rho_ux,rho_uy,rho_uz,Bx,By,Bz),&
!!$OMP & SHARED(d_ux_dx,d_uy_dy,d_uz_dz,P,E),&
!!$OMP & SHARED(d_Bx_dy,d_Bx_dz,d_By_dx,d_By_dz,d_Bz_dx,d_Bz_dy,curl_x_B,curl_y_B,curl_z_B,B2),&
!!$OMP & SHARED(Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3),&
!!$OMP & SHARED(Energy_x,Energy_y,Energy_z,E_Visc),&
!!$OMP & SHARED(Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2) PRIVATE(i,j,k)
!!$OMP DO
! 
!$acc parallel firstprivate(spheat, eta) present(ux,uy,uz,rho_ux,rho_uy,rho_uz,Bx,By,Bz, d_ux_dx,d_uy_dy,d_uz_dz,E, B2, d_Bx_dy,d_Bx_dz,d_By_dx,d_By_dz,d_Bz_dx,d_Bz_dy) &
!$acc & private(curl_x_B, curl_y_B, curl_z_B) present(P, Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3, Energy_x,Energy_y,Energy_z,E_Visc, Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2) 
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx
!Recompute inverse transform of ux,uy,uz
    ux(i,j,k) = rho_ux(i,j,k)/rho(i,j,k)
    uy(i,j,k) = rho_uy(i,j,k)/rho(i,j,k)
    uz(i,j,k) = rho_uz(i,j,k)/rho(i,j,k)
    Bx(i,j,k) = Bx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    By(i,j,k) = By(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    Bz(i,j,k) = Bz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    ! FFTW Normalisation.
    d_ux_dx(i,j,k) = d_ux_dx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_uy_dy(i,j,k) = d_uy_dy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_uz_dz(i,j,k) = d_uz_dz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
   
    d_Bx_dy(i,j,k) = d_Bx_dy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_Bx_dz(i,j,k) = d_Bx_dz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_By_dx(i,j,k) = d_By_dx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_By_dz(i,j,k) = d_By_dz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_Bz_dx(i,j,k) = d_Bz_dx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_Bz_dy(i,j,k) = d_Bz_dy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
    ! Evaluate Curl of Magnetic Field.
    curl_x_B = d_Bz_dy(i,j,k) - d_By_dz(i,j,k)
    curl_y_B = d_Bz_dx(i,j,k) - d_Bx_dz(i,j,k)
    curl_z_B = d_By_dx(i,j,k) - d_Bx_dy(i,j,k)
    
    ! Evaluate ressure
    P(i,j,k) = CS*CS*rho(i,j,k) !( spheat - 1.0d0 ) * ( E(i,j,k) &
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
                      - eta * ( By(i,j,k) * curl_z_B + Bz(i,j,k) * curl_y_B) 
                      
    Energy_y(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * uy(i,j,k) &
                      - ux(i,j,k)*Bx(i,j,k)*By(i,j,k) - uy(i,j,k)*By(i,j,k)*By(i,j,k) - uz(i,j,k)*By(i,j,k)*Bz(i,j,k) &
                      + eta * ( Bx(i,j,k) * curl_z_B - Bz(i,j,k) * curl_x_B)
                      
    Energy_z(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * uz(i,j,k) &
                      - ux(i,j,k)*Bz(i,j,k)*Bx(i,j,k) - uy(i,j,k)*Bz(i,j,k)*By(i,j,k) - uz(i,j,k)*Bz(i,j,k)*Bz(i,j,k) & 
                      + eta * ( Bx(i,j,k) * curl_y_B + By(i,j,k) * curl_x_B)
                      
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
!$acc end parallel
write(*,*) "Derive 6"

!!$OMP END DO
!!$OMP END PARALLEL


!$acc host_data use_device(Mom_x_1, Mom_x_2, Mom_x_3, Mom_y_1, Mom_y_2, Mom_y_3, Mom_z_1, Mom_z_2, Mom_z_3, Energy_x, Energy_y, Energy_z, E_Visc, Mag_x_1, Mag_x_2, Mag_y_1, Mag_y_2) &
!$acc & use_device(Mag_z_1, Mag_z_2)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_1, Mom_x_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mom_x_1), C_LOC(Mom_x_1), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_2, Mom_x_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mom_x_2), C_LOC(Mom_x_2), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_3, Mom_x_3_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mom_x_3), C_LOC(Mom_x_3), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_1, Mom_y_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mom_y_1), C_LOC(Mom_y_1), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_2, Mom_y_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mom_y_2), C_LOC(Mom_y_2), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_3, Mom_y_3_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mom_y_3), C_LOC(Mom_y_3), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_1, Mom_z_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mom_z_1), C_LOC(Mom_z_1), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_2, Mom_z_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mom_z_2), C_LOC(Mom_z_2), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_3, Mom_z_3_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mom_z_3), C_LOC(Mom_z_3), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_x, Energy_x_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Energy_x), C_LOC(Energy_x), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_y, Energy_y_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Energy_y), C_LOC(Energy_y), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_z, Energy_z_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Energy_z), C_LOC(Energy_z), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, E_Visc, E_Visc_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(E_Visc), C_LOC(E_Visc), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_x_1, Mag_x_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mag_x_1), C_LOC(Mag_x_1), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_x_2, Mag_x_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mag_x_2), C_LOC(Mag_x_2), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_y_1, Mag_y_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mag_y_1), C_LOC(Mag_y_1), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_y_2, Mag_y_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mag_y_2), C_LOC(Mag_y_2), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_z_1, Mag_z_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mag_z_1), C_LOC(Mag_z_1), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_z_2, Mag_z_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Mag_z_2), C_LOC(Mag_z_2), CUFFT_D2Z)

  !call dfftw_destroy_plan_ (plan_backward)
  !call dfftw_destroy_plan_ (plan_forward)
!$acc end host_data


!Recompute tranforms of rho_ux,rho_uy,rho_uz,Bx,By,Bz,E
!!$acc parallel present(rho,rho_ux,rho_uy,rho_uz,Bx,By,Bz)
!!$acc loop collapse(3)
!do i = 1,Nx
!  do j = 1,Ny
!    do k = 1,Nz
!      rho(i,j,k) = rho(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      rho_ux(i,j,k) = rho_ux(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      rho_uy(i,j,k) = rho_uy(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      rho_uz(i,j,k) = rho_uz(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      Bx(i,j,k) = Bx(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      By(i,j,k) = By(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      Bz(i,j,k) = Bz(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      E(i,j,k) = E(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!    enddo
!   enddo
!enddo
!!$acc end parallel

!$acc host_data use_device(ux,uy,uz,rho,rho_ux,rho_uy,rho_uz,Bx,By,Bz,E)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho), C_LOC(rho), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho_ux), C_LOC(rho_ux), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho_uy), C_LOC(rho_uy), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho_uz), C_LOC(rho_uz), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Bx), C_LOC(Bx), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(By), C_LOC(By), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Bz), C_LOC(Bz), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(E), C_LOC(E), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(ux), C_LOC(ux), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(uy), C_LOC(uy), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(uz), C_LOC(uz), CUFFT_D2Z)

!$acc end host_data
write(*,*) "Derive 7"


!!$OMP PARALLEL SHARED(mu,eta,i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k),&
!!$OMP & SHARED(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k),&
!!$OMP & SHARED(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k),&
!!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,E_Visc_k)&
!!$OMP & SHARED(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k),&
!!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k,nu,Fx_k,Fy_k,Fz_k),&
!!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k),&
!!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new),&
!!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new) PRIVATE(i,j,k,kx,ky,kz)
!!$OMP DO 

!$acc parallel firstprivate(mu, eta,Lx,Ly,Lz) private(i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k) &
!$acc & private(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k) present(rho_ux,rho_uy,rho_uz,ux,uy,uz) &
!$acc & private(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k,i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k) &
!$acc & private(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k) &
!$acc & private(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k) &
!$acc & private(kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k) &
!$acc & present(Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3) &
!$acc & present(Mom_z_1,Mom_z_2,Mom_z_3,Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2) &
!$acc & present(E_Visc,Fx,Fy,Fz,Energy_x,Energy_y,Energy_z,Bx,By,Bz) &
!$acc & present(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new) &
!$acc & present(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new) 
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx/2+1
    



    kx = 2.0d0*pi*dfloat(i-1)/Lx
    if (j <= Ny/2) then
      ky = 2.0d0*pi*dfloat(j-1)/Ly
    else
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    endif
    if (k <= Nz/2) then
      kz = 2.0d0*pi*dfloat(k-1)/Lz
    else
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
    endif

    if (kx >= Nx/3 .and. ky >= Ny/3 .and. kz >= Nz/3) then!dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k = 0.0d0
        i_ky_rho_uy_k = 0.0d0
        i_kz_rho_uz_k = 0.0d0
        i_kx_Mom_x_1_k = 0.0d0
        i_ky_Mom_x_2_k = 0.0d0
        i_kz_Mom_x_3_k = 0.0d0
        i_kx_Mom_y_1_k = 0.0d0
        i_ky_Mom_y_2_k  = 0.0d0
        i_kz_Mom_y_3_k  = 0.0d0
        i_kx_Mom_z_1_k = 0.0d0
        i_ky_Mom_z_2_k = 0.0d0
        i_kz_Mom_z_3_k  = 0.0d0
        kx2_ux_k = 0.0d0
        ky2_ux_k = 0.0d0
        kz2_ux_k = 0.0d0
        kx2_uy_k = 0.0d0
        ky2_uy_k = 0.0d0
        kz2_uy_k = 0.0d0
        kx2_uz_k = 0.0d0
        ky2_uz_k = 0.0d0
        kz2_uz_k = 0.0d0
        i_kx_Energy_x_k = 0.0d0
        i_ky_Energy_y_k = 0.0d0
        i_kz_Energy_z_k = 0.0d0
        E_Visc(2*i-1,j,k) = 0.0d0
        E_Visc(2*i,j,k) = 0.0d0
        i_ky_Mag_x_1_k = 0.0d0
        i_kz_Mag_x_2_k = 0.0d0
        i_kx_Mag_y_1_k = 0.0d0
        i_kz_Mag_y_2_k = 0.0d0
        i_kx_Mag_z_1_k = 0.0d0
        i_ky_Mag_z_2_k = 0.0d0
        kx2_Bx_k = 0.0d0
        ky2_Bx_k = 0.0d0
        kz2_Bx_k = 0.0d0
        kx2_By_k = 0.0d0
        ky2_By_k = 0.0d0
        kz2_By_k = 0.0d0
        kx2_Bz_k = 0.0d0
        ky2_Bz_k = 0.0d0
        kz2_Bz_k = 0.0d0
    else
        i_kx_rho_ux_k = (0.0d0,1.0d0)*kx*cmplx(rho_ux(2*i-1,j,k),rho_ux(2*i,j,k),8) 
        i_ky_rho_uy_k = (0.0d0,1.0d0)*ky*cmplx(rho_uy(2*i-1,j,k),rho_uy(2*i,j,k),8) 
        i_kz_rho_uz_k = (0.0d0,1.0d0)*kz*cmplx(rho_uz(2*i-1,j,k),rho_uz(2*i,j,k),8) 
      i_kx_Mom_x_1_k = (0.0d0,1.0d0)*kx*cmplx(Mom_x_1(2*i-1,j,k),Mom_x_1(2*i,j,k),8)    
      i_ky_Mom_x_2_k = (0.0d0,1.0d0)*ky*cmplx(Mom_x_2(2*i-1,j,k),Mom_x_2(2*i,j,k),8)
      i_kz_Mom_x_3_k = (0.0d0,1.0d0)*kz*cmplx(Mom_x_3(2*i-1,j,k),Mom_x_3(2*i,j,k),8)   
      i_kx_Mom_y_1_k = (0.0d0,1.0d0)*kx*cmplx(Mom_y_1(2*i-1,j,k),Mom_y_1(2*i,j,k),8)
      i_ky_Mom_y_2_k = (0.0d0,1.0d0)*ky*cmplx(Mom_y_2(2*i-1,j,k),Mom_y_2(2*i,j,k),8)
      i_kz_Mom_y_3_k = (0.0d0,1.0d0)*kz*cmplx(Mom_y_3(2*i-1,j,k),Mom_y_3(2*i,j,k),8)
      i_kx_Mom_z_1_k = (0.0d0,1.0d0)*kx*cmplx(Mom_z_1(2*i-1,j,k), Mom_z_1(2*i,j,k),8) 
      i_ky_Mom_z_2_k = (0.0d0,1.0d0)*ky*cmplx(Mom_z_2(2*i-1,j,k),Mom_z_2(2*i,j,k),8) 
      i_kz_Mom_z_3_k = (0.0d0,1.0d0)*kz*cmplx(Mom_z_3(2*i-1,j,k),Mom_z_3(2*i,j,k),8)
      kx2_ux_k = kx*kx*cmplx(ux(2*i-1,j,k),ux(2*i,j,k),8)
      ky2_ux_k = ky*ky*cmplx(ux(2*i-1,j,k),ux(2*i,j,k),8)
      kz2_ux_k = kz*kz*cmplx(ux(2*i-1,j,k),ux(2*i,j,k),8)
      kx2_uy_k = kx*kx*cmplx(uy(2*i-1,j,k),uy(2*i,j,k),8)
      ky2_uy_k = ky*ky*cmplx(uy(2*i-1,j,k),uy(2*i,j,k),8)
      kz2_uy_k = kz*kz*cmplx(uy(2*i-1,j,k),uy(2*i,j,k),8)
      kx2_uz_k = kx*kx*cmplx(uz(2*i-1,j,k),uz(2*i,j,k),8)
      ky2_uz_k = ky*ky*cmplx(uz(2*i-1,j,k),uz(2*i,j,k),8)
      kz2_uz_k = kz*kz*cmplx(uz(2*i-1,j,k),uz(2*i,j,k),8)
      i_kx_Energy_x_k = (0.0d0,1.0d0)*kx*cmplx(Energy_x(2*i-1,j,k),Energy_x(2*i,j,k),8)
      i_ky_Energy_y_k = (0.0d0,1.0d0)*ky*cmplx(Energy_y(2*i-1,j,k),Energy_y(2*i,j,k),8)
      i_kz_Energy_z_k = (0.0d0,1.0d0)*kz*cmplx(Energy_z(2*i-1,j,k),Energy_z(2*i,j,k),8)
      i_ky_Mag_x_1_k = (0.0d0,1.0d0)*ky*cmplx(Mag_x_1(2*i-1,j,k),Mag_x_1(2*i,j,k),8)
      i_kz_Mag_x_2_k = (0.0d0,1.0d0)*kz*cmplx(Mag_x_2(2*i-1,j,k),Mag_x_2(2*i,j,k),8)
      i_kx_Mag_y_1_k = (0.0d0,1.0d0)*kx*cmplx(Mag_y_1(2*i-1,j,k),Mag_y_1(2*i,j,k),8)
      i_kz_Mag_y_2_k = (0.0d0,1.0d0)*kz*cmplx(Mag_y_2(2*i-1,j,k),Mag_y_2(2*i,j,k),8)
      i_kx_Mag_z_1_k = (0.0d0,1.0d0)*kx*cmplx(Mag_z_1(2*i-1,j,k),Mag_z_1(2*i,j,k),8)
      i_ky_Mag_z_2_k = (0.0d0,1.0d0)*ky*cmplx(Mag_z_2(2*i-1,j,k),Mag_z_2(2*i,j,k),8)
    
      kx2_Bx_k = kx*kx*cmplx(Bx(2*i-1,j,k),Bx(2*i,j,k),8)
      ky2_Bx_k = ky*ky*cmplx(Bx(2*i-1,j,k),Bx(2*i,j,k),8)
      kz2_Bx_k = kz*kz*cmplx(Bx(2*i-1,j,k),Bx(2*i,j,k),8)
      kx2_By_k = kx*kx*cmplx(By(2*i-1,j,k),By(2*i,j,k),8)
      ky2_By_k = ky*ky*cmplx(By(2*i-1,j,k),By(2*i,j,k),8)
      kz2_By_k = kz*kz*cmplx(By(2*i-1,j,k),By(2*i,j,k),8)
      kx2_Bz_k = kx*kx*cmplx(Bz(2*i-1,j,k),Bz(2*i,j,k),8)
      ky2_Bz_k = ky*ky*cmplx(Bz(2*i-1,j,k),Bz(2*i,j,k),8)
      kz2_Bz_k = kz*kz*cmplx(Bz(2*i-1,j,k),Bz(2*i,j,k),8)
    endif




      ! Density Equation.
      d_rho_k_dt_new(i,j,k) = - ( i_kx_rho_ux_k + i_ky_rho_uy_k + i_kz_rho_uz_k)
    
      ! Momentum Equation.
      d_rho_ux_k_dt_new(i,j,k) = - ( i_kx_Mom_x_1_k + i_ky_Mom_x_2_k + i_kz_Mom_x_3_k ) &
                                 - ( kx2_ux_k + ky2_ux_k + kz2_ux_k ) / mu !+ cmplx(Fx(2*i-1,j,k),Fx(2*i,j,k),8) 
    
      d_rho_uy_k_dt_new(i,j,k) = - ( i_kx_Mom_y_1_k + i_ky_Mom_y_2_k + i_kz_Mom_y_3_k ) &
                                 - ( kx2_uy_k + ky2_uy_k + kz2_uy_k ) / mu !+ cmplx(Fy(2*i-1,j,k),Fy(2*i,j,k),8)
    
      d_rho_uz_k_dt_new(i,j,k) = - ( i_kx_Mom_z_1_k + i_ky_Mom_z_2_k + i_kz_Mom_z_3_k ) &
                                 - ( kx2_uz_k + ky2_uz_k + kz2_uz_k ) / mu !+ cmplx(Fz(2*i-1,j,k),Fz(2*i,j,k),8) 
    
      ! Energy Equation.
      d_E_k_dt_new(i,j,k) = - ( i_kx_Energy_x_k + i_ky_Energy_y_k + i_kz_Energy_z_k ) &
                            + cmplx(E_Visc(2*i-1,j,k),E_Visc(2*i,j,k),8) / mu
    
      ! Magnetic Field Equation.
      d_Bx_k_dt_new(i,j,k) = + ( i_ky_Mag_x_1_k + i_kz_Mag_x_2_k ) &
                             - ( kx2_Bx_k + ky2_Bx_k + kz2_Bx_k ) / eta
                           
      d_By_k_dt_new(i,j,k) = - ( i_kx_Mag_y_1_k - i_kz_Mag_y_2_k ) &
                             - ( kx2_By_k + ky2_By_k + kz2_By_k ) / eta
                           
      d_Bz_k_dt_new(i,j,k) = - ( i_kx_Mag_z_1_k + i_ky_Mag_z_2_k ) &                       
                             - ( kx2_Bz_k + ky2_Bz_k + kz2_Bz_k ) / eta
    enddo                       
  enddo
enddo  
!$acc end parallel
!!$acc end data
!!$OMP END DO
!!$OMP END PARALLEL
write(*,*) "After derive"


! Non-Linear Term Evaluator.
!  call derive (Nx,Ny,Nz,Nh,pi,time,rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz, &
!               d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
!               d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
!               d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)

! Time Solvers.
! Adams-Bashforth
!  call ab (Nx,Ny,Nz,Nh,pi,time,rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz, &
!           rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new, &
!           d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
!           d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
!           d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)
!======================================================================================
!!$acc data copyin(rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz, d_rho_ux_k_dt_new, d_rho_uy_k_dt_new, d_rho_uz_k_dt_new, d_rho_k_dt_old, d_rho_ux_k_dt_old) &
!!$acc & copyin(d_rho_uy_k_dt_old, d_rho_uz_k_dt_old, d_E_k_dt_new, d_E_k_dt_old, d_Bx_k_dt_new, d_By_k_dt_new, d_Bz_k_dt_new, d_Bx_k_dt_old, d_By_k_dt_old, d_Bz_k_dt_old) &
!!$acc & copyout(rho_k_new, rho_ux_k_new, rho_uy_k_new, rho_uz_k_new, E_k_new, Bx_k_new, By_k_new, Bz_k_new)

!$acc parallel present(rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz, d_rho_ux_k_dt_new, d_rho_uy_k_dt_new, d_rho_uz_k_dt_new, d_rho_ux_k_dt_old) &
!$acc & present(d_rho_uy_k_dt_old, d_rho_uz_k_dt_old, d_E_k_dt_new, d_E_k_dt_old, d_Bx_k_dt_new, d_By_k_dt_new, d_Bz_k_dt_new, d_Bx_k_dt_old, d_By_k_dt_old, d_Bz_k_dt_old) &
!$acc & present(rho_k_new, rho_ux_k_new, rho_uy_k_new, rho_uz_k_new, E_k_new, Bx_k_new, By_k_new, Bz_k_new)

!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nh
      ! Density Equation Evolution.
      rho_k_new(i,j,k) = cmplx(rho(2*i-1,j,k),rho(2*i,j,k),8) + ( (3.0d0/2.0d0)*d_rho_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_k_dt_old(i,j,k) )*dt
      
      ! Momentum Equation Evolution.
      rho_ux_k_new(i,j,k) = cmplx(rho_ux(2*i-1,j,k),rho_ux(2*i,j,k),8) + ( (3.0d0/2.0d0)*d_rho_ux_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_ux_k_dt_old(i,j,k) )*dt
      rho_uy_k_new(i,j,k) = cmplx(rho_uy(2*i-1,j,k),rho_uy(2*i,j,k),8) + ( (3.0d0/2.0d0)*d_rho_uy_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_uy_k_dt_old(i,j,k) )*dt
      rho_uz_k_new(i,j,k) = cmplx(rho_uz(2*i-1,j,k),rho_uz(2*i,j,k),8) + ( (3.0d0/2.0d0)*d_rho_uz_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_uz_k_dt_old(i,j,k) )*dt
      
      ! Energy Equation Evolution.
      E_k_new(i,j,k) = cmplx(E(2*i-1,j,k),E(2*i,j,k),8) !+ ( (3.0d0/2.0d0)*d_E_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_E_k_dt_old(i,j,k) )*dt 
      
      ! Energy Equation Evolution.
      Bx_k_new(i,j,k) = cmplx(Bx(2*i-1,j,k),Bx(2*i,j,k),8) + ( (3.0d0/2.0d0)*d_Bx_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_Bx_k_dt_old(i,j,k) )*dt
      By_k_new(i,j,k) = cmplx(By(2*i-1,j,k),By(2*i,j,k),8) + ( (3.0d0/2.0d0)*d_By_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_By_k_dt_old(i,j,k) )*dt
      Bz_k_new(i,j,k) = cmplx(Bz(2*i-1,j,k),Bz(2*i,j,k),8) + ( (3.0d0/2.0d0)*d_Bz_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_Bz_k_dt_old(i,j,k) )*dt
    enddo
  enddo
enddo
!$acc end parallel
!!$acc end data
write(*,*) "After ab"
!!$OMP END DO
!!$OMP END PARALLEL

!!$acc data copyin(d_rho_k_dt_new, d_rho_ux_k_dt_new, d_rho_uy_k_dt_new, d_rho_uz_k_dt_new, d_E_k_dt_new, d_Bx_k_dt_new, d_By_k_dt_new, d_Bz_k_dt_new, rho_k_new, rho_ux_k_new, rho_uy_k_new, rho_uz_k_new, E_k_new, Bx_k_new, By_k_new, Bz_k_new) &
!!$acc & copyout(d_rho_k_dt_old, d_rho_ux_k_dt_old, d_rho_uy_k_dt_old, d_rho_uz_k_dt_old, d_E_k_dt_old, d_Bx_k_dt_old, d_By_k_dt_old, d_Bz_k_dt_old) &
!!$acc & copyout(rho, rho_ux, rho_uy, rho_uz, E, div_B, B2, j2, A2, u2, jx, jy, jz, Ax, Ay, Az) &
!!$acc & create(omega_x, omega_y, omega_z) &
!!$acc & copyout(P, ux, uy, uz, Bx, By, Bz, omega2, div_B)

!!$OMP PARALLEL SHARED(d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old),&
!!$OMP & SHARED(d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old),&
!!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new),& 
!!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new),& 
!!$OMP & SHARED(rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k),&
!!$OMP & SHARED(rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new),&
!!$OMP & PRIVATE(i,j,k)
!!$OMP DO 

!$acc parallel present(d_rho_k_dt_new, d_rho_ux_k_dt_new, d_rho_uy_k_dt_new, d_rho_uz_k_dt_new, d_E_k_dt_new, d_Bx_k_dt_new, d_By_k_dt_new, d_Bz_k_dt_new) &
!$acc & present(rho_k_new, rho_ux_k_new, rho_uy_k_new, rho_uz_k_new, E_k_new, Bx_k_new, By_k_new, Bz_k_new) &
!$acc & present(d_rho_k_dt_old, d_rho_ux_k_dt_old, d_rho_uy_k_dt_old, d_rho_uz_k_dt_old, d_E_k_dt_old, d_Bx_k_dt_old, d_By_k_dt_old, d_Bz_k_dt_old, rho, rho_ux, rho_uy, rho_uz, E, Bx, By, Bz)

!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny
     do i = 1,Nx/2+1
      ! Set the Variables in Proper Format for Next Time Iteration.
      d_rho_k_dt_old(i,j,k) = d_rho_k_dt_new(i,j,k)
    
      d_rho_ux_k_dt_old(i,j,k) = d_rho_ux_k_dt_new(i,j,k)
      d_rho_uy_k_dt_old(i,j,k) = d_rho_uy_k_dt_new(i,j,k)
      d_rho_uz_k_dt_old(i,j,k) = d_rho_uz_k_dt_new(i,j,k)
    
      d_E_k_dt_old(i,j,k) = d_E_k_dt_new(i,j,k)
    
      d_Bx_k_dt_old(i,j,k) = d_Bx_k_dt_new(i,j,k)
      d_By_k_dt_old(i,j,k) = d_By_k_dt_new(i,j,k)
      d_Bz_k_dt_old(i,j,k) = d_Bz_k_dt_new(i,j,k)
    
      rho(2*i-1,j,k) = real(rho_k_new(i,j,k))
      rho(2*i,j,k) = aimag(rho_k_new(i,j,k))
    
      rho_ux(2*i-1,j,k) = real(rho_ux_k_new(i,j,k))
      rho_ux(2*i,j,k) = aimag(rho_ux_k_new(i,j,k))
      rho_uy(2*i-1,j,k) = real(rho_uy_k_new(i,j,k))
      rho_uy(2*i,j,k) = aimag(rho_uy_k_new(i,j,k))
      rho_uz(2*i-1,j,k) = real(rho_uz_k_new(i,j,k))
      rho_uz(2*i,j,k) = aimag(rho_uz_k_new(i,j,k))
  
      E(2*i-1,j,k) = real(E_k_new(i,j,k))
      E(2*i,j,k) = aimag(E_k_new(i,j,k))
    
      Bx(2*i-1,j,k) = real(Bx_k_new(i,j,k))
      Bx(2*i,j,k) = aimag(Bx_k_new(i,j,k))
      By(2*i-1,j,k) = real(By_k_new(i,j,k))
      By(2*i,j,k) = aimag(By_k_new(i,j,k))
      Bz(2*i-1,j,k) = real(Bz_k_new(i,j,k))
      Bz(2*i,j,k) = aimag(Bz_k_new(i,j,k))
    
      ! Keep Backup of the Arrays for FFTW.
!      rho_k_dum(i,j,k) = rho_k(i,j,k)   
    
!      rho_ux_k_dum(i,j,k) = rho_ux_k(i,j,k)
!      rho_uy_k_dum(i,j,k) = rho_uy_k(i,j,k)
!      rho_uz_k_dum(i,j,k) = rho_uz_k(i,j,k)

!      E_k_dum(i,j,k) = E_k(i,j,k)
    
!      Bx_k_dum(i,j,k) = Bx_k(i,j,k)
!      By_k_dum(i,j,k) = By_k(i,j,k)
!      Bz_k_dum(i,j,k) = Bz_k(i,j,k)
    enddo
  enddo
enddo
!$acc end parallel
write(*,*) "Final loop1"

!!$OMP END DO
!!$OMP END PARALLEL 

! Evaluate Divergence of B in Spectral Space.

!!$OMP PARALLEL SHARED(Lx,Ly,Lz,Bx_k,By_k,Bz_k,div_B_k) PRIVATE(i,j,k,kx,ky,kz)
!!$OMP DO
!$acc parallel firstprivate(Lx, Ly, Lz) present(jx,jy,jz,Ax,Ay,Az,Bx,By,Bz,div_B) 
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx/2+1
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      if (j <= Ny/2) then
        ky = 2.0d0*pi*dfloat(j-1)/Ly
      else
        ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      endif
      if (k <= Nz/2) then
        kz = 2.0d0*pi*dfloat(k-1)/Lz
      else
        kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      endif
      jx(2*i-1,j,k) = (-1.0d0)*ky*Bz(2*i,j,k) + (-1.0d0)*kz*By(2*i,j,k)
      jx(2*i,j,k) = (1.0d0)*ky*Bz(2*i-1,j,k) + (1.0d0)*kz*By(2*i-1,j,k)
      jy(2*i-1,j,k) = (-1.0d0)*kz*Bx(2*i,j,k) + (-1.0d0)*kx*Bz(2*i,j,k)
      jy(2*i,j,k) = (1.0d0)*kz*Bx(2*i-1,j,k) + (1.0d0)*kx*Bz(2*i-1,j,k)
      jz(2*i-1,j,k) = (-1.0d0)*kx*By(2*i,j,k) + (-1.0d0)*ky*Bx(2*i,j,k)
      jz(2*i,j,k) = (1.0d0)*kx*By(2*i-1,j,k) + (1.0d0)*ky*Bx(2*i-1,j,k)
        if (i == 1 .and. j == 1 .and. k == 1) then
          Ax(2*i-1,j,k) = jx(2*i-1,j,k)
          Ax(2*i,j,k) = jx(2*i,j,k)
          Ay(2*i-1,j,k) = jy(2*i-1,j,k)
          Ay(2*i,j,k) = jy(2*i,j,k)
          Az(2*i-1,j,k) = jz(2*i-1,j,k)
          Az(2*i,j,k) = jz(2*i,j,k)
        else  
          Ax(2*i-1,j,k) = jx(2*i-1,j,k)/(kx*kx+ky*ky+kz*kz)
          Ax(2*i,j,k) = jx(2*i,j,k)/(kx*kx+ky*ky+kz*kz)
          Ay(2*i-1,j,k) = jy(2*i-1,j,k)/(kx*kx+ky*ky+kz*kz)
          Ay(2*i,j,k) = jy(2*i,j,k)/(kx*kx+ky*ky+kz*kz)
          Az(2*i-1,j,k) = jz(2*i-1,j,k)/(kx*kx+ky*ky+kz*kz)
          Az(2*i,j,k) = jz(2*i,j,k)/(kx*kx+ky*ky+kz*kz)
        endif
      div_B(2*i-1,j,k) = (-1.0d0)*kx*Bx(2*i,j,k) + (-1.0d0)*ky*By(2*i,j,k) + (-1.0d0)*kz*Bz(2*i,j,k)
      div_B(2*i,j,k) = (1.0d0)*kx*Bx(2*i-1,j,k) + (1.0d0)*ky*By(2*i-1,j,k) + (1.0d0)*kz*Bz(2*i-1,j,k)
    enddo
  enddo
enddo
!$acc end parallel
write(*,*) "Final loop2"

!!$OMP END DO
!!$OMP END PARALLEL 

 !$acc host_data use_device(rho, rho_ux, rho_uy, rho_uz, E, Bx, By, Bz, div_B, jx, jy, jz, Ax, Ay, Az)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_k_dum, rho, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(rho), C_LOC(rho), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_ux_k_dum, rho_ux, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(rho_ux), C_LOC(rho_ux), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uy_k_dum, rho_uy, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(rho_uy), C_LOC(rho_uy), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uz_k_dum, rho_uz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(rho_uz), C_LOC(rho_uz), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, E_k_dum, E, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(E), C_LOC(E), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bx_k_dum, Bx, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(Bx), C_LOC(Bx), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, By_k_dum, By, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(By), C_LOC(By), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bz_k_dum, Bz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(Bz), C_LOC(Bz), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, div_B_k, div_B, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(div_B), C_LOC(div_B), CUFFT_Z2D)

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, jx_k, jx, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(jx), C_LOC(jx), CUFFT_Z2D)

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, jy_k, jy, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(jy), C_LOC(jy), CUFFT_Z2D)

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, jz_k, jz, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(jz), C_LOC(jz), CUFFT_Z2D)

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Ax_k, Ax, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(Ax), C_LOC(Ax), CUFFT_Z2D)

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Ay_k, Ay, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(Ay), C_LOC(Ay), CUFFT_Z2D)

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Az_k, Az, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(Az), C_LOC(Az), CUFFT_Z2D)

 !$acc end host_data
write(*,*) "Final 3"

!!$OMP PARALLEL SHARED(spheat,rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz,div_B), &
!!$OMP & SHARED(ux,uy,uz,P) PRIVATE(i,j,k)
!!$OMP DO

!$acc parallel firstprivate(spheat) present(rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz,div_B, ux,uy,uz,jx,jy,jz,Ax,Ay,Az,u2,j2,A2,P,B2)
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx
      ! FFTW Normalisation.
      jx(i,j,k) = jx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      jy(i,j,k) = jy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      jz(i,j,k) = jz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      
      Ax(i,j,k) = Ax(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Ay(i,j,k) = Ay(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Az(i,j,k) = Az(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      
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
  
      ! Evaluate Square of Velocity and Magnetic Field.
      u2(i,j,k) = ux(i,j,k)*ux(i,j,k) + uy(i,j,k)*uy(i,j,k) + uz(i,j,k)*uz(i,j,k)
      B2(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)
      j2(i,j,k) = jx(i,j,k)*jx(i,j,k) + jy(i,j,k)*jy(i,j,k) + jz(i,j,k)*jz(i,j,k)
      A2(i,j,k) = Ax(i,j,k)*Ax(i,j,k) + Ay(i,j,k)*Ay(i,j,k) + Az(i,j,k)*Az(i,j,k)
  
      ! Keep Backup of the Arrays for FFTW.
!      ux_dum(i,j,k) = ux(i,j,k)
!      uy_dum(i,j,k) = uy(i,j,k)
!      uz_dum(i,j,k) = uz(i,j,k)
  
      ! Evaluate Pressure
      P(i,j,k) = ( spheat - 1.0d0 ) * ( E(i,j,k) - 0.50d0 * &
                 ( rho_ux(i,j,k)*ux(i,j,k)+rho_uy(i,j,k)*uy(i,j,k)+rho_uz(i,j,k)*uz(i,j,k) - B2(i,j,k) ) )
    enddo
  enddo
enddo   
!$acc end parallel
write(*,*) "Final 4"








!==================================================================================================================



! Velocity and Magnetic Field Interpolation from Grid to Tracer Particle...
!!$acc parallel loop vector_length(256) & 
!!$acc& present(ux,uy,uz,Bx,By,Bz,x,y,z,uxp,uyp,uzp,Bpx,Bpy,Bpz,x,y,z)
!$acc parallel loop present(ux,uy,uz,Bx,By,Bz,xp,yp,zp,uxp,uyp,uzp,Bpx,Bpy,Bpz,x,y,z)

do i = 1,Np

  m = int(xp(i)/dx) + 1

  n = int(yp(i)/dy) + 1
  
  q = int(zp(i)/dz) + 1


if (m==Nx .and. n/=Ny .and. q/=Nz) then

  uxp(i) = ux(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           ux(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           ux(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           ux(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           ux(m,n,q+1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           ux(1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m,n+1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           ux(1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpx(i) = Bx(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bx(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bx(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bx(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bx(m,n,q+1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m,n+1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bx(1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  uyp(i) = uy(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           uy(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           uy(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uy(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uy(m,n,q+1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uy(1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m,n+1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uy(1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  
  Bpy(i) = By(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           By(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           By(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           By(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           By(m,n,q+1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           By(1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           By(m,n+1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           By(1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  uzp(i) = uz(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           uz(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           uz(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uz(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uz(m,n,q+1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uz(1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m,n+1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uz(1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpz(i) = Bz(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bz(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bz(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bz(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bz(m,n,q+1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m,n+1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bz(1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

end if

if (m/=Nx .and. n==Ny .and. q/=Nz) then

  uxp(i) = ux(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           ux(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           ux(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           ux(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           ux(m,n,q+1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m+1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m,1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           ux(m+1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpx(i) = Bx(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bx(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bx(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bx(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bx(m,n,q+1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m+1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m,1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bx(m+1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  uyp(i) = uy(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           uy(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           uy(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uy(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uy(m,n,q+1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m+1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m,1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uy(m+1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  
  Bpy(i) = By(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           By(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           By(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           By(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           By(m,n,q+1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           By(m+1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           By(m,1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           By(m+1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  uzp(i) = uz(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           uz(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           uz(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uz(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uz(m,n,q+1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m+1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m,1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uz(m+1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpz(i) = Bz(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bz(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bz(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bz(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bz(m,n,q+1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m+1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m,1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bz(m+1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

endif

if (m/=Nx .and. n/=Ny .and. q==Nz) then

  uxp(i) = ux(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           ux(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           ux(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           ux(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           ux(m,n,1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m+1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m,n+1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           ux(m+1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpx(i) = Bx(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           Bx(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           Bx(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bx(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bx(m,n,1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m+1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m,n+1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bx(m+1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  uyp(i) = uy(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           uy(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           uy(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uy(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uy(m,n,1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m+1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m,n+1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uy(m+1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpy(i) = By(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           By(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           By(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           By(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           By(m,n,1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           By(m+1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           By(m,n+1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           By(m+1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h
  
  uzp(i) = uz(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           uz(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           uz(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uz(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uz(m,n,1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m+1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m,n+1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uz(m+1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpz(i) = Bz(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           Bz(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           Bz(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bz(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bz(m,n,1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m+1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m,n+1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bz(m+1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

endif

if (m==Nx .and. n==Ny .and. q/=Nz ) then

  uxp(i) = ux(m,n,q) * (Lx-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           ux(1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           ux(m,1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           ux(1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           ux(m,n,q+1) * (Lx-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           ux(1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m,1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           ux(1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpx(i) = Bx(m,n,q) * (Lx-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bx(1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bx(m,1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bx(1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bx(m,n,q+1) * (Lx-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m,1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bx(1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  uyp(i) = uy(m,n,q) * (Lx-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           uy(1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           uy(m,1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uy(1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uy(m,n,q+1) * (Lx-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uy(1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m,1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uy(1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpy(i) = By(m,n,q) * (Lx-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           By(1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           By(m,1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           By(1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           By(m,n,q+1) * (Lx-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           By(1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           By(m,1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           By(1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h
  
  uzp(i) = uz(m,n,q) * (Lx-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           uz(1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           uz(m,1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uz(1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uz(m,n,q+1) * (Lx-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uz(1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m,1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uz(1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpz(i) = Bz(m,n,q) * (Lx-xp(i)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bz(1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bz(m,1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bz(1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bz(m,n,q+1) * (Lx-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(1,n,q+1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m,1,q+1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bz(1,1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

endif

if (m==Nx .and. n/=Ny .and. q==Nz ) then

  uxp(i) = ux(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           ux(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           ux(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           ux(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           ux(m,n,1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           ux(1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m,n+1,1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           ux(1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpx(i) = Bx(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           Bx(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           Bx(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bx(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bx(m,n,1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m,n+1,1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bx(1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  uyp(i) = uy(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           uy(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           uy(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uy(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uy(m,n,1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uy(1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m,n+1,1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uy(1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpy(i) = By(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           By(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           By(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           By(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           By(m,n,1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           By(1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           By(m,n+1,1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           By(1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h
  
  uzp(i) = uz(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           uz(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           uz(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uz(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uz(m,n,1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uz(1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m,n+1,1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uz(1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpz(i) = Bz(m,n,q) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           Bz(1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (Lz-zp(i)) * h + &
           Bz(m,n+1,q) * (Lx-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bz(1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bz(m,n,1) * (Lx-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(1,n,1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m,n+1,1) * (Lx-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bz(1,n+1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

endif

if (m/=Nx .and. n==Ny .and. q==Nz ) then

  uxp(i) = ux(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           ux(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           ux(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           ux(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           ux(m,n,1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m+1,n,1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m,1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           ux(m+1,1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpx(i) = Bx(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           Bx(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           Bx(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bx(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bx(m,n,1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m+1,n,1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m,1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bx(m+1,1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  uyp(i) = uy(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           uy(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           uy(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uy(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uy(m,n,1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m+1,n,1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m,1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uy(m+1,1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpy(i) = By(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           By(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           By(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           By(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           By(m,n,1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           By(m+1,n,1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           By(m,1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           By(m+1,1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h
  
  uzp(i) = uz(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           uz(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           uz(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uz(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           uz(m,n,1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m+1,n,1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m,1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uz(m+1,1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpz(i) = Bz(m,n,q) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           Bz(m+1,n,q) * (xp(i)-x(m)) * (Ly-yp(i)) * (Lz-zp(i)) * h + &
           Bz(m,1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bz(m+1,1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (Lz-zp(i)) * h + &
           Bz(m,n,1) * (x(m+1)-xp(i)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m+1,n,1) * (xp(i)-x(m)) * (Ly-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m,1,1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bz(m+1,1,1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

endif

if (m/=Nx .and. n/=Ny .and. q/=Nz) then

  uxp(i) = ux(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           ux(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           ux(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           ux(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           ux(m,n,q+1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m+1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           ux(m,n+1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           ux(m+1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpx(i) = Bx(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bx(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bx(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bx(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bx(m,n,q+1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m+1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bx(m,n+1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bx(m+1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  uyp(i) = uy(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           uy(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           uy(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uy(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uy(m,n,q+1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m+1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uy(m,n+1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uy(m+1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpy(i) = By(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           By(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           By(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           By(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           By(m,n,q+1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           By(m+1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           By(m,n+1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           By(m+1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h
  
  uzp(i) = uz(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           uz(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           uz(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uz(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           uz(m,n,q+1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m+1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           uz(m,n+1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           uz(m+1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h

  Bpz(i) = Bz(m,n,q) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bz(m+1,n,q) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (z(q+1)-zp(i)) * h + &
           Bz(m,n+1,q) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bz(m+1,n+1,q) * (xp(i)-x(m)) * (yp(i)-y(n)) * (z(q+1)-zp(i)) * h + &
           Bz(m,n,q+1) * (x(m+1)-xp(i)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m+1,n,q+1) * (xp(i)-x(m)) * (y(n+1)-yp(i)) * (zp(i)-z(q)) * h + &
           Bz(m,n+1,q+1) * (x(m+1)-xp(i)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h + &
           Bz(m+1,n+1,q+1) * (xp(i)-x(m)) * (yp(i)-y(n)) * (zp(i)-z(q)) * h
end if

end do ! i
!$acc end parallel

! Boris Algorithm...
!!$acc parallel loop vector_length(256) present(vxp,vyp,vzp,uxp,uyp,uzp,x,y,z,Bpx,Bpy,Bpz) &
!$acc parallel loop present(vxp,vyp,vzp,uxp,uyp,uzp,x,y,z,Bpx,Bpy,Bpz) &
!$acc& private(vxd,vyd,vzd)
 
 do i=1,Np
  tx = 0.5d0*dt*Bpx(i)
  ty = 0.5d0*dt*Bpy(i)
  tz = 0.5d0*dt*Bpz(i)
  
  sx = 2.0d0*tx/(1.0d0+tx**2+ty**2+tz**2)
  sy = 2.0d0*ty/(1.0d0+tx**2+ty**2+tz**2)
  sz = 2.0d0*tz/(1.0d0+tx**2+ty**2+tz**2)

  vxd(i) = uxp(i) + (uyp(i)*tz-uzp(i)*ty)
  vyd(i) = uyp(i) - (uxp(i)*tz-uzp(i)*tx)
  vzd(i) = uzp(i) + (uxp(i)*ty-uyp(i)*tx)
  
  vxp(i) = uxp(i) + (vyd(i)*sz-vzd(i)*sy)
  vyp(i) = uyp(i) - (vxd(i)*sz-vzd(i)*sx)
  vzp(i) = uzp(i) + (vxd(i)*sy-vyd(i)*sx)
  
  xp(i) = xp(i) + vxp(i)*dt
  yp(i) = yp(i) + vyp(i)*dt
  zp(i) = zp(i) + vzp(i)*dt
  
  ! Periodic Boundary Condition Implemented...
  xp(i) = xp(i) - (int(xp(i)/Lx))*Lx
    if (xp(i) .lt. 0.0d0) then
    xp(i) = xp(i) + Lx
    endif
  
  yp(i) = yp(i) - (int(yp(i)/Ly))*Ly
    if (yp(i) .lt. 0.0d0) then
    yp(i) = yp(i) + Ly
    endif

  zp(i) = zp(i) - (int(zp(i)/Lz))*Lz
    if (zp(i) .lt. 0.0d0) then
    zp(i) = zp(i) + Lz
    endif

  enddo

!$acc end parallel

  if (mod(float(t),1000.0) == 0.0) then
  do i = 1,Np
  write(t+100,*) xp(i),yp(i),zp(i),uxp(i),uyp(i),uzp(i)
  end do
  endif
  
!==================================================================================================================















write(*,*) "Final 4a"

  
!!$OMP END DO
!!$OMP END PARALLEL 

!$acc host_data use_device(ux, uy, uz)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(ux), C_LOC(ux), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(uy), C_LOC(uy), CUFFT_D2Z)

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(uz), C_LOC(uz), CUFFT_D2Z)
! print*, "After fft2 in main" 

!$acc end host_data
write(*,*) "Final 5"

! Evaluate Vorticity in Spectral Space.

!!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,omega_x_k,omega_y_k,omega_z_k),&
!!$OMP & PRIVATE(i,j,k,kx,ky,kz)
!!$OMP DO 
!$acc parallel firstprivate(Lx, Ly, Lz) present(ux,uy,uz,omega_x,omega_y,omega_z)
!$acc loop collapse(3)
do k = 1, Nz
  do j = 1, Ny
    do i = 1, Nx/2+1
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      if (j <= Ny/2) then
        ky = 2.0d0*pi*dfloat(j-1)/Ly
      else
        ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      endif
      if (k <= Nz/2) then
        kz = 2.0d0*pi*dfloat(k-1)/Lz
      else
        kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      endif
     
      omega_x(2*i-1,j,k) = (-1.0d0)*ky*uz(2*i,j,k) - (-1.0d0)*kz*uy(2*i,j,k)
      omega_x(2*i,j,k) = (1.0d0)*ky*uz(2*i-1,j,k) - (1.0d0)*kz*uy(2*i-1,j,k)
      omega_y(2*i-1,j,k) = (-1.0d0)*kx*uz(2*i,j,k) - (-1.0d0)*kz*ux(2*i,j,k)
      omega_y(2*i,j,k) = (1.0d0)*kx*uz(2*i-1,j,k) - (1.0d0)*kz*ux(2*i-1,j,k)
      omega_z(2*i-1,j,k) = (-1.0d0)*kx*uy(2*i,j,k) - (-1.0d0)*ky*ux(2*i,j,k)
      omega_z(2*i,j,k) = (1.0d0)*kx*uy(2*i-1,j,k) - (1.0d0)*ky*ux(2*i-1,j,k)

     ! omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
     ! omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
     ! omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo
!$acc end parallel
write(*,*) "Final 6"

!!$OMP END DO
!!$OMP END PARALLEL

!$acc host_data use_device(omega_x, omega_y, omega_z)
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_x_k_dum, omega_x, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(omega_x), C_LOC(omega_x), CUFFT_Z2D)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_y_k_dum, omega_y, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(omega_y), C_LOC(omega_y), CUFFT_Z2D)
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_z_k_dum, omega_z, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(omega_z), C_LOC(omega_z), CUFFT_Z2D)

 !$acc end host_data
write(*,*) "Final 7"

  !call dfftw_destroy_plan_ (plan_forward)
  !call dfftw_destroy_plan_ (plan_backward)

! FFTW Normalisation and omega^2 Evaluation.

!!$OMP PARALLEL SHARED(omega_x,omega_y,omega_z,omega2) PRIVATE(i,j,k)
!!$OMP DO 
!$acc parallel present(omega_x,omega_y,omega_z,omega2)
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx

      omega_x(i,j,k) = omega_x(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega_y(i,j,k) = omega_y(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega_z(i,j,k) = omega_z(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega2(i,j,k) = omega_x(i,j,k)**2 + omega_y(i,j,k)**2 + omega_z(i,j,k)**2
    enddo
  enddo
enddo  
!$acc end parallel
!!$acc end data
write(*,*) "Final 8"

!!$OMP END DO
!!$OMP END PARALLEL

! Evaluate Energy Spectra at a Specified Time.
if (t >= int((time_max-time_max/10.0d0)/dt)) then

!!! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file. !!!
! !!$OMP PARALLEL SHARED(ux_k,uy_k,uz_k) PRIVATE(i,j,k) 
! !!$OMP DO REDUCTION (+:Ek) 

!$acc parallel loop collapse(3) present(ux,uy,uz,Ek) 

 do k = 1,Nz
    do j = 1,Ny
      do i = 1,Nx/2+1
        Ek(i,j,k) = Ek(i,j,k) + sqrt(sqrt(ux(2*i-1,j,k)**2+ux(2*i,j,k)**2)**2 + sqrt(uy(2*i-1,j,k)**2+uy(2*i,j,k)**2)**2 + sqrt(uz(2*i-1,j,k)**2+uz(2*i,j,k)**2)**2)&
                    /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    !    Bk(i,j,k) = Bk(i,j,k) + sqrt(abs(Bx(i,j,k))**2 + abs(By(i,j,k))**2 + abs(Bz(i,j,k))**2)&
    !                /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      enddo
   enddo
  enddo  
!$acc end parallel
! !!$OMP END PARALLEL

endif 

!Additional recompute of inverse transform of ux,uy,uz
!!$acc data copy(ux,uy,uz)
!$acc host_data use_device(ux,uy,uz) 
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(ux), C_LOC(ux), CUFFT_Z2D)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(uy), C_LOC(uy), CUFFT_Z2D)
  call executeCUFFT3D(fftPlanZ2DMain, C_LOC(uz), C_LOC(uz), CUFFT_Z2D)
!$acc end host_data
!$acc parallel present(ux,uy,uz)
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx
      ux(i,j,k) = ux(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      uy(i,j,k) = uy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      uz(i,j,k) = uz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    enddo
   enddo
enddo
!$acc end parallel
!!$acc end data

Pressure = 0.0d0; Energy = 0.0d0; T_Energy = 0.0d0; B_Field = 0.0d0
 C0 = 0.0d0; C1 = 0.0d0; C2 = 0.0d0; C3 = 0.0d0; div_B_Tot = 0.0d0
 Hf = 0.0d0; Hm = 0.0d0; HT = 0.0d0; HuB = 0.0d0; HAw = 0.0d0; Rayleigh = 0.0d0
   
! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file.
! !$OMP PARALLEL SHARED(mu_0,t,dt,x,y,z,omega_x,omega_y,omega_z,omega2,rho,ux,uy,uz) PRIVATE(i,j,k)
! !$OMP DO REDUCTION (+:Pressure,Energy,y_Energy,B_Field,C0,C1,C2,C3,div_B_Tot) 

!$acc parallel loop collapse(3) reduction(+:Pressure,Energy,T_Energy,B_Field,C0,C1,C2,C3,div_B_Tot,Hf,Hm,HT,HuB,HAw,Rayleigh) present(P,ux,uy,uz,rho,Bx,By,Bz,u2,B2,Ax,Ay,Az,jx,jy,jz,omega2,div_B)
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx
      ! Evaluate Pressure
      Pressure = Pressure + P(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Energy
      Energy = Energy + u2(i,j,k)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Growth Rate.
      !T_Energy = T_Energy + E(i,j,k)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Magnetic Field.
      B_Field = B_Field + B2(i,j,k)/(2.0d0*mu_0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
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
    enddo
  enddo
enddo  
      
!$acc end parallel
! !$OMP END PARALLEL

if (mod(float(t),1000.0) == 0.0) then
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx
      ! Write Grid Data in State Files.  
        write(t+110,*) ux(i,j,k),uy(i,j,k),uz(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    enddo
  enddo
enddo  
  close(t+110)
endif

if (mod(float(t),1000.0) == 0.0) then
do k = 1,Nz
  do j = 1,Ny
    do i = 1,Nx
        if (i == 1) then
        write(t+120,*) y(j),z(k),u2(i,j,k),B2(i,j,k),omega2(i,j,k),j2(i,j,k),A2(i,j,k)
        elseif (j == 1) then
        write(t+130,*) x(i),z(k),u2(i,j,k),B2(i,j,k),omega2(i,j,k),j2(i,j,k),A2(i,j,k)
        elseif (k == 1) then
        write(t+140,*) x(i),y(j),u2(i,j,k),B2(i,j,k),omega2(i,j,k),j2(i,j,k),A2(i,j,k)
        endif
    enddo
  enddo
enddo  
  close(t+120)
  close(t+130)
  close(t+140)
endif

if (mod(float(t),100.0) == 0.0) then
  write(40,*) time,Energy,B_Field!,C1,C2,C3,Hf,Hm,HT,HuB,HAw,Rayleigh,div_B_Tot
  call flush(40) 
endif

if (mod(float(t),1000.0) == 0.0) then
  close(t+110)
endif

!Recompute tranforms of rho,rho_ux,rho_uy,rho_uz,Bx,By,Bz,E
!!$acc data copy(rho,rho_ux,rho_uy,rho_uz,Bx,By,Bz,E)
!!$acc parallel present(rho,rho_ux,rho_uy,rho_uz,Bx,By,Bz,E)
!!$acc loop collapse(3)
!do i = 1,Nx
!  do j = 1,Ny
!    do k = 1,Nz
!      rho(i,j,k) = rho(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      rho_ux(i,j,k) = rho_ux(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      rho_uy(i,j,k) = rho_uy(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      rho_uz(i,j,k) = rho_uz(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      Bx(i,j,k) = Bx(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      By(i,j,k) = By(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      Bz(i,j,k) = Bz(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      E(i,j,k) = E(i,j,k) * (dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!    enddo
!   enddo
!enddo
!!$acc end parallel

!$acc host_data use_device(rho,rho_ux,rho_uy,rho_uz,Bx,By,Bz,E)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho), C_LOC(rho), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho_ux), C_LOC(rho_ux), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho_uy), C_LOC(rho_uy), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(rho_uz), C_LOC(rho_uz), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Bx), C_LOC(Bx), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(By), C_LOC(By), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(Bz), C_LOC(Bz), CUFFT_D2Z)
  call executeCUFFT3D(fftPlanD2ZMain, C_LOC(E), C_LOC(E), CUFFT_D2Z)

!$acc end host_data
!!$acc end data
!$acc parallel loop collapse(3) present(Bx,By,Bz,Bk)
 do k = 1,Nz
    do j = 1,Ny
      do i = 1,Nx/2+1
        Bk(i,j,k) = Bk(i,j,k) + sqrt(sqrt(Bx(2*i-1,j,k)**2+Bx(2*i,j,k)**2)**2 + sqrt(By(2*i-1,j,k)**2+By(2*i,j,k)**2)**2 + sqrt(Bz(2*i-1,j,k)**2+Bz(2*i,j,k)**2)**2)&
                    /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      enddo
   enddo
  enddo  

!$acc end parallel 
! close(t+100)

enddo ! time

!!$acc end data !======================================================================
t2 = omp_get_wtime()

!write(5,*) "Time taken for the run =",(t2 - t1)/(60.0d0*60.0d0),"Hours"
write(5,*) "Time taken for the run =",(t2 - t1),"Seconds"

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
    !write(35,*) i,j,k,sqrt(float((i-1)**2)+float((j-1)**2)+float((k-1)**2)),&
    !          abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*float(i-1)/Lx
    ky = 2.0d0*pi*float(j-1)/Ly
    kz = 2.0d0*pi*float((k-1)-Nz)/Lz
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
      Ek(i,j,k) = 0.0d0
      Bk(i,j,k) = 0.0d0
      endif
    !write(35,*) i,j,k-Nz,sqrt(float((i-1)**2)+float((j-1)**2)+float(((k-1)-Nz)**2)),&
    !          abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
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
    !write(35,*) i,j-Ny,k,sqrt(float((i-1)**2)+float(((j-1)-Ny)**2)+float((k-1)**2)),&
    !          abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*float(i-1)/Lx
    ky = 2.0d0*pi*float((j-1)-Ny)/Ly
    kz = 2.0d0*pi*float((k-1)-Nz)/Lz
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
      Ek(i,j,k) = 0.0d0
      Bk(i,j,k) = 0.0d0
      endif
    !write(35,*) i,j-Ny,k-Nz,sqrt(float((i-1)**2)+float(((j-1)-Ny)**2)+float(((k-1)-Nz)**2)),&
    !          abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
  enddo  
enddo 
 
 close(5)
 close(35)
 close(40)
 call destroyCUFFTPlan3D(fftPlanD2ZMain)
 call destroyCUFFTPlan3D(fftPlanZ2DMain)

!Perform deallocations
deallocate(x,y,z)
deallocate(rho,ux,uy,uz)
deallocate(rho_ux,rho_uy,rho_uz)
deallocate(omega_x,omega_y,omega_z,omega2,P,E,Bx,By,Bz,B2,div_B)
!deallocate(rho_k,ux_k,uy_k,uz_k,rho_ux_k,rho_uy_k,rho_uz_k,omega_x_k,omega_y_k,omega_z_k)
!deallocate(P_k,E_k,Ek,Bk,Bx_k,By_k,Bz_k,div_B_k)
deallocate(Ek,Bk)
deallocate(rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new)
deallocate(E_k_new,Bx_k_new,By_k_new,Bz_k_new)
deallocate(d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old)
deallocate(d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old)
deallocate(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new)
deallocate(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new)
!Variables used in derive
deallocate(Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3)

deallocate(d_ux_dx,d_uy_dy,d_uz_dz,Fx,Fy,Fz,Energy_x,Energy_y,Energy_z,E_Visc)

deallocate(Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2,d_Bx_dy,d_By_dx,d_Bx_dz,d_By_dz,d_Bz_dx,d_Bz_dy)

!deallocate(curl_x_B,curl_y_B,curl_z_B)

!End variables used in derive


end program PSMHD3  

