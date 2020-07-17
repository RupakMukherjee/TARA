!Module for invoking accfft

module accfft
    interface
    subroutine setGPUDevice(gpuid) bind (C, name = 'setGPUDevice')
         use iso_c_binding
         implicit none
         integer (c_int), value :: gpuid
     end subroutine setGPUDevice
     end interface

    interface
    subroutine accfftCreateComm(nprocs,comm) bind (C, name = 'accfftCreateComm')
          use iso_c_binding
          implicit none
          integer (c_int), value :: nprocs
          type (c_ptr) :: comm
    end subroutine accfftCreateComm
    end interface

    interface
    subroutine accfftDestroyComm(comm) bind (C, name = 'accfftDestroyComm')
          use iso_c_binding
          implicit none
          type (c_ptr) :: comm
    end subroutine accfftDestroyComm
    end interface

    interface
    subroutine accfftLocalSizeDFTR2CGPU(n,isize,istart,osize,ostart,comm, &
                                 allocsize) bind (C, name = 'accfftLocalSizeDFTR2CGPU')
          use iso_c_binding
          implicit none
          type (c_ptr), value ::n,isize,istart,osize,ostart,comm,allocsize
    end subroutine accfftLocalSizeDFTR2CGPU
    end interface

    interface
    subroutine accfftCreatePlan3DR2CGPU(n,idata,odata,comm,plan) &
     bind (C, name = 'accfftCreatePlan3DR2CGPU')
          use iso_c_binding
          implicit none
          type (c_ptr), value :: n,idata,odata,comm
          type (c_ptr) :: plan
    end subroutine accfftCreatePlan3DR2CGPU
    end interface

    interface
    subroutine accfftDestroyPlan3DR2CGPU(plan) bind (C, name = 'accfftDestroyPlan3DR2CGPU')
            use iso_c_binding
            implicit none
            type (c_ptr), value ::plan
    end subroutine accfftDestroyPlan3DR2CGPU
    end interface

    interface
    subroutine accfftExecuteR2CGPU(plan,idata,odata,accffttimer) &
      bind (C, name = 'accfftExecuteR2CGPU')
           use iso_c_binding
           implicit none
           type (c_ptr), value :: plan, idata, odata,accffttimer
    end subroutine accfftExecuteR2CGPU
    end interface

    interface
    subroutine accfftExecuteC2RGPU(plan,idata,odata,accffttimer) &
      bind (C, name = 'accfftExecuteC2RGPU')
           use iso_c_binding
           implicit none
           type (c_ptr), value :: plan, idata, odata, accffttimer
    end subroutine accfftExecuteC2RGPU
    end interface

end module accfft

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
use mpi 
use iso_c_binding
use accfft 
use openacc
implicit none
double precision rand !added by Naga to avoid the following compilation error:PGF90-S-0038-Symbol, rand, has not been explicitly declared (PSMHD3.f95)
! Define Grid Size.
integer ( kind = 4 ), parameter :: Np = 0 !10000
integer ( kind = 4 ), parameter :: Nx = 64
integer ( kind = 4 ), parameter :: Ny = 64
integer ( kind = 4 ), parameter :: Nz = 64
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

!multiGPU
integer ( kind = 4 ) :: gridSize(3),isize(3),osize(3),istart(3),ostart(3)
integer ( kind = 4 ) :: iglobal,jglobal,kglobal,lindex
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

integer ( kind = 4 ) i,j,k,t,m,n,q
real ( kind = 8 ) Lx,Ly,Lz,dx,dy,dz,kx,ky,kz,time,time_min,time_max,dt,G1,G2,G3,h,sx,sy,sz,tx,ty,tz
real ( kind = 8 ) spheat,rho0,U0,mu,kf,mode,sigma,uy0,MS,CS,Re,Rm,PM,P0,eta,mu_0,MA,VA,B0,ba,A,B,C
real ( kind = 8 ) Pressure,Energy,Energy0,T_Energy,y_Energy0,B_Field,C0,C1,C2,C3,div_B_tot,Hf,Hm,HT,HuB,HAw,Rayleigh
real ( kind = 8) Energyfinal, B_Fieldfinal

!real ( kind = 8 ) x(Nx),y(Ny),z(Nz),rho(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
real ( kind = 8 ), dimension (:), allocatable :: x,y,z
real (kind = 8), dimension (:,:,:), allocatable :: rho,ux,uy,uz,uxfull,uyfull,uzfull

!real ( kind = 8 ) rho_ux(Nx,Ny,Nz),rho_uy(Nx,Ny,Nz),rho_uz(Nx,Ny,Nz)
real ( kind = 8 ), dimension (:,:,:), allocatable :: rho_ux,rho_uy,rho_uz

!real ( kind = 8 ) omega_x(Nx,Ny,Nz),omega_y(Nx,Ny,Nz),omega_z(Nx,Ny,Nz),omega2(Nx,Ny,Nz)
real ( kind = 8 ), dimension (:,:,:), allocatable :: omega_x,omega_y,omega_z,omega2

!real ( kind = 8 ) P(Nx,Ny,Nz),E(Nx,Ny,Nz),Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz),div_B(Nx,Ny,Nz)
real ( kind = 8 ), dimension (:,:,:), allocatable :: P,E,Bx,By,Bz,B2,div_B,jx,jy,jz,Ax,Ay,Az,j2,A2,u2
real ( kind = 8 ), dimension (:,:,:), allocatable :: Bxfull,Byfull,Bzfull
complex ( kind = 8 ), dimension (:,:,:), allocatable :: jx_k,jy_k,jz_k,Ax_k,Ay_k,Az_k

real ( kind = 8 ), dimension (:), allocatable ::  xp,yp,zp,uxp,uyp,uzp,Bpx,Bpy,Bpz,vxp,vyp,vzp,vxd,vyd,vzd

!Backup variables for FFTW have been commented
!real ( kind = 8 ) rho_dum(Nx,Ny,Nz),ux_dum(Nx,Ny,Nz),uy_dum(Nx,Ny,Nz),uz_dum(Nx,Ny,Nz)
!real ( kind = 8 ) rho_ux_dum(Nx,Ny,Nz),rho_uy_dum(Nx,Ny,Nz),rho_uz_dum(Nx,Ny,Nz)
!real ( kind = 8 ) P_dum(Nx,Ny,Nz),E_dum(Nx,Ny,Nz),Bx_dum(Nx,Ny,Nz),By_dum(Nx,Ny,Nz),Bz_dum(Nx,Ny,Nz)

!complex ( kind = 8 ) rho_k(Nh,Ny,Nz),ux_k(Nh,Ny,Nz),uy_k(Nh,Ny,Nz),uz_k(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: rho_k,ux_k,uy_k,uz_k

!complex ( kind = 8 ) rho_ux_k(Nh,Ny,Nz),rho_uy_k(Nh,Ny,Nz),rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: rho_ux_k,rho_uy_k,rho_uz_k

!complex ( kind = 8 ) omega_x_k(Nh,Ny,Nz),omega_y_k(Nh,Ny,Nz),omega_z_k(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: omega_x_k,omega_y_k,omega_z_k

!complex ( kind = 8 ) P_k(Nh,Ny,Nz),E_k(Nh,Ny,Nz),Ek(Nh,Ny,Nz),Bk(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: P_k,E_k,Ek,Bk
complex ( kind = 8 ), dimension (:,:,:), allocatable :: Ekfull,Bkfull

!complex ( kind = 8 ) Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz),div_B_k(Nh,Ny,Nz)
complex ( kind = 8 ), dimension (:,:,:), allocatable :: Bx_k,By_k,Bz_k,div_B_k

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
complex ( kind = 8 ), dimension(:,:,:), allocatable :: Mom_x_1_k,Mom_x_2_k,Mom_x_3_k,Mom_y_1_k,Mom_y_2_k,Mom_y_3_k,Mom_z_1_k,Mom_z_2_k,Mom_z_3_k

real ( kind = 8 ), dimension(:,:,:), allocatable :: d_ux_dx,d_uy_dy,d_uz_dz,Fx,Fy,Fz
complex ( kind = 8 ), dimension(:,:,:), allocatable :: d_ux_dx_k,d_uy_dy_k,d_uz_dz_k,Fx_k,Fy_k,Fz_k
real ( kind = 8 ), dimension(:,:,:), allocatable :: Energy_x,Energy_y,Energy_z,E_Visc
complex ( kind = 8 ), dimension(:,:,:), allocatable :: Energy_x_k,Energy_y_k,Energy_z_k,E_Visc_k
real ( kind = 8 ), dimension(:,:,:), allocatable :: Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2
complex ( kind = 8 ), dimension(:,:,:), allocatable :: Mag_x_1_k,Mag_x_2_k,Mag_y_1_k,Mag_y_2_k,Mag_z_1_k,Mag_z_2_k
real ( kind = 8 ), dimension(:,:,:), allocatable :: d_Bx_dy,d_By_dx,d_Bx_dz
complex ( kind = 8 ), dimension(:,:,:), allocatable :: d_Bx_dy_k,d_By_dx_k,d_Bx_dz_k
real ( kind = 8 ), dimension(:,:,:), allocatable :: d_By_dz,d_Bz_dx,d_Bz_dy
complex ( kind = 8 ), dimension(:,:,:), allocatable :: d_By_dz_k,d_Bz_dx_k,d_Bz_dy_k
real ( kind = 8 ) :: curl_x_B,curl_y_B,curl_z_B
complex ( kind = 8 ) :: i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k
complex ( kind = 8 ) :: i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k
complex ( kind = 8 ) :: i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k
complex ( kind = 8 ) :: i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k
complex ( kind = 8 ) :: kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k
complex ( kind = 8 ) :: i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k
complex ( kind = 8 ) :: kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k

!End variables used in derive


!mpi and accfft
type (c_ptr) :: fftPlan
integer mpierr, nprocs, myrank,allocmax,ngpus,totgpus
character*4 :: numgpus
type (c_ptr) :: comm

real ( kind = 8 ) t1,t2,accffttimer(5)

common/comm/time_max,Lx,Ly,Lz,dx,dy,dz,spheat,mu,ms,CS,mu_0,eta,dt,kf

integer,parameter :: seed = 99999999

call srand(seed)

!Init MPI

call MPI_INIT(mpierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,mpierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
write(*,*) 'Nprocs: ',nprocs,' Rank: ',myrank

!Map MPI process to GPU
!call sleep(1000)
totgpus = acc_get_num_devices(acc_device_nvidia)
call getenv("NGPUS",numgpus)
write(*,*) 'Total GPUS :', totgpus, 'Used GPUs: ', numgpus
read(numgpus,*) ngpus
!ngpus = 2
ngpus = min(ngpus,totgpus)
write(*,*) 'GPUs used: ', ngpus
call setGPUDevice(mod(myrank,ngpus))
call acc_set_device_num(mod(myrank,ngpus),acc_device_nvidia)
write(*,*) 'Rank:',myrank,' setting gpu:', mod(myrank,ngpus)

!Init accfft
call accfftCreateComm(nprocs,comm)

!Do all the allocations
allocate(x(Nx))
allocate(y(Ny))
allocate(z(Nz))

!Get accfft local sizes
 gridSize(1) = Nz
 gridSize(2) = Ny
 gridSize(3) = Nx
 call accfftLocalSizeDFTR2CGPU(C_LOC(gridSize), C_LOC(isize), C_LOC(istart), &
             C_LOC(osize), C_LOC(ostart), comm, C_LOC(allocmax))

write(*,*) 'isize(1):',isize(1)
write(*,*) 'isize(2):',isize(2)
write(*,*) 'isize(3):',isize(3)
write(*,*) 'osize(1):',osize(1)
write(*,*) 'osize(2):',osize(2)
write(*,*) 'osize(3):',osize(3)

!Each process only allocated its portion of the data. Only rank 0 has the
!complete output arrays 
!Assume that Np is completely divisible by nprocs
allocate(xp(Np))
allocate(yp(Np))
allocate(zp(Np))
allocate(uxp(Np/nprocs))
allocate(uyp(Np/nprocs))
allocate(uzp(Np/nprocs))
allocate(Bpx(Np/nprocs))
allocate(Bpy(Np/nprocs))
allocate(Bpz(Np/nprocs))
allocate(vxp(Np/nprocs))
allocate(vyp(Np/nprocs))
allocate(vzp(Np/nprocs))
allocate(vxd(Np/nprocs))
allocate(vyd(Np/nprocs))
allocate(vzd(Np/nprocs))

!Flipping the sizes for fortran-c compatibility (column vs row)

allocate(rho(isize(3), isize(2), isize(1)))
allocate(ux(isize(3), isize(2), isize(1)))
allocate(uy(isize(3), isize(2), isize(1)))
allocate(uz(isize(3), isize(2), isize(1)))

allocate(rho_ux(isize(3), isize(2), isize(1)))
allocate(rho_uy(isize(3), isize(2), isize(1)))
allocate(rho_uz(isize(3), isize(2), isize(1)))

allocate(omega_x(isize(3), isize(2), isize(1)))
allocate(omega_y(isize(3), isize(2), isize(1)))
allocate(omega_z(isize(3), isize(2), isize(1)))
allocate(omega2(isize(3), isize(2), isize(1)))

allocate(P(isize(3), isize(2), isize(1)))
allocate(E(isize(3), isize(2), isize(1)))
allocate(Bx(isize(3), isize(2), isize(1)))
allocate(By(isize(3), isize(2), isize(1)))
allocate(Bz(isize(3), isize(2), isize(1)))
allocate(B2(isize(3), isize(2), isize(1)))
allocate(div_B(isize(3), isize(2), isize(1)))
allocate(div_B_k(osize(3), osize(2), osize(1)))
allocate(jx(isize(3), isize(2), isize(1)))
allocate(jy(isize(3), isize(2), isize(1)))
allocate(jz(isize(3), isize(2), isize(1)))
allocate(jx_k(osize(3), osize(2), osize(1)))
allocate(jy_k(osize(3), osize(2), osize(1)))
allocate(jz_k(osize(3), osize(2), osize(1)))
allocate(j2(isize(3), isize(2), isize(1)))
allocate(Ax(isize(3), isize(2), isize(1)))
allocate(Ay(isize(3), isize(2), isize(1)))
allocate(Az(isize(3), isize(2), isize(1)))
allocate(Ax_k(osize(3), osize(2), osize(1)))
allocate(Ay_k(osize(3), osize(2), osize(1)))
allocate(Az_k(osize(3), osize(2), osize(1)))
allocate(A2(isize(3), isize(2), isize(1)))
allocate(u2(isize(3), isize(2), isize(1)))


allocate(rho_k(osize(3),osize(2),osize(1)))
allocate(ux_k(osize(3),osize(2),osize(1)))
allocate(uy_k(osize(3),osize(2),osize(1)))
allocate(uz_k(osize(3),osize(2),osize(1)))

allocate(rho_ux_k(osize(3),osize(2),osize(1)))
allocate(rho_uy_k(osize(3),osize(2),osize(1)))
allocate(rho_uz_k(osize(3),osize(2),osize(1)))

allocate(omega_x_k(osize(3),osize(2),osize(1)))
allocate(omega_y_k(osize(3),osize(2),osize(1)))
allocate(omega_z_k(osize(3),osize(2),osize(1)))

allocate(P_k(osize(3),osize(2),osize(1)))
allocate(E_k(osize(3),osize(2),osize(1)))
allocate(Ek(osize(3), osize(2), osize(1)))
allocate(Bk(osize(3), osize(2), osize(1)))
Ek = 0.0
Bk = 0.0
if (myrank.eq.0) then
allocate(Ekfull(Nh,Ny,Nz))
allocate(Bkfull(Nh,Ny,Nz))
endif
allocate(uxfull(Nx,Ny,Nz))
allocate(uyfull(Nx,Ny,Nz))
allocate(uzfull(Nx,Ny,Nz))
allocate(Bxfull(Nx,Ny,Nz))
allocate(Byfull(Nx,Ny,Nz))
allocate(Bzfull(Nx,Ny,Nz))

allocate(Bx_k(osize(3),osize(2),osize(1)))
allocate(By_k(osize(3),osize(2),osize(1)))
allocate(Bz_k(osize(3),osize(2),osize(1)))

allocate(rho_k_new(osize(3), osize(2), osize(1)))
allocate(rho_ux_k_new(osize(3), osize(2), osize(1)))
allocate(rho_uy_k_new(osize(3), osize(2), osize(1)))
allocate(rho_uz_k_new(osize(3), osize(2), osize(1)))

allocate(E_k_new(osize(3), osize(2), osize(1)))
allocate(Bx_k_new(osize(3), osize(2), osize(1)))
allocate(By_k_new(osize(3), osize(2), osize(1)))
allocate(Bz_k_new(osize(3), osize(2), osize(1)))

allocate(d_rho_k_dt_old(osize(3), osize(2), osize(1)))
d_rho_k_dt_old = 0.0
allocate(d_rho_ux_k_dt_old(osize(3), osize(2), osize(1)))
d_rho_ux_k_dt_old = 0.0
allocate(d_rho_uy_k_dt_old(osize(3), osize(2), osize(1)))
d_rho_uy_k_dt_old = 0.0
allocate(d_rho_uz_k_dt_old(osize(3), osize(2), osize(1)))
d_rho_uz_k_dt_old = 0.0

allocate(d_E_k_dt_old(osize(3), osize(2), osize(1)))
d_E_k_dt_old = 0.0
allocate(d_Bx_k_dt_old(osize(3), osize(2), osize(1)))
d_Bx_k_dt_old = 0.0
allocate(d_By_k_dt_old(osize(3), osize(2), osize(1)))
d_By_k_dt_old = 0.0
allocate(d_Bz_k_dt_old(osize(3), osize(2), osize(1)))
d_Bz_k_dt_old = 0.0

allocate(d_rho_k_dt_new(osize(3), osize(2), osize(1)))
allocate(d_rho_ux_k_dt_new(osize(3), osize(2), osize(1)))
allocate(d_rho_uy_k_dt_new(osize(3), osize(2), osize(1)))
allocate(d_rho_uz_k_dt_new(osize(3), osize(2), osize(1)))

allocate(d_E_k_dt_new(osize(3), osize(2), osize(1)))
allocate(d_Bx_k_dt_new(osize(3), osize(2), osize(1)))
allocate(d_By_k_dt_new(osize(3), osize(2), osize(1)))
allocate(d_Bz_k_dt_new(osize(3), osize(2), osize(1)))

!Variables used in derive
allocate(Mom_x_1(isize(3),isize(2),isize(1)))
allocate(Mom_x_2(isize(3),isize(2),isize(1)))
allocate(Mom_x_3(isize(3),isize(2),isize(1)))
allocate(Mom_y_1(isize(3),isize(2),isize(1)))
allocate(Mom_y_2(isize(3),isize(2),isize(1)))
allocate(Mom_y_3(isize(3),isize(2),isize(1)))
allocate(Mom_z_1(isize(3),isize(2),isize(1)))
allocate(Mom_z_2(isize(3),isize(2),isize(1)))
allocate(Mom_z_3(isize(3),isize(2),isize(1)))


allocate(Mom_x_1_k(osize(3),osize(2),osize(1)))
allocate(Mom_x_2_k(osize(3),osize(2),osize(1)))
allocate(Mom_x_3_k(osize(3),osize(2),osize(1)))
allocate(Mom_y_1_k(osize(3),osize(2),osize(1)))
allocate(Mom_y_2_k(osize(3),osize(2),osize(1)))
allocate(Mom_y_3_k(osize(3),osize(2),osize(1)))
allocate(Mom_z_1_k(osize(3),osize(2),osize(1)))
allocate(Mom_z_2_k(osize(3),osize(2),osize(1)))
allocate(Mom_z_3_k(osize(3),osize(2),osize(1)))

allocate(d_ux_dx(isize(3),isize(2),isize(1)))
allocate(d_uy_dy(isize(3),isize(2),isize(1)))
allocate(d_uz_dz(isize(3),isize(2),isize(1)))
allocate(d_ux_dx_k(osize(3),osize(2),osize(1)))
allocate(d_uy_dy_k(osize(3),osize(2),osize(1)))
allocate(d_uz_dz_k(osize(3),osize(2),osize(1)))


allocate(Fx(isize(3),isize(2),isize(1)))
allocate(Fy(isize(3),isize(2),isize(1)))
allocate(Fz(isize(3),isize(2),isize(1)))
allocate(Fx_k(osize(3),osize(2),osize(1)))
allocate(Fy_k(osize(3),osize(2),osize(1)))
allocate(Fz_k(osize(3),osize(2),osize(1)))

allocate(Energy_x(isize(3),isize(2),isize(1)))
allocate(Energy_y(isize(3),isize(2),isize(1)))
allocate(Energy_z(isize(3),isize(2),isize(1)))
allocate(E_Visc(isize(3),isize(2),isize(1)))

allocate(Energy_x_k(osize(3),osize(2),osize(1)))
allocate(Energy_y_k(osize(3),osize(2),osize(1)))
allocate(Energy_z_k(osize(3),osize(2),osize(1)))
allocate(E_Visc_k(osize(3),osize(2),osize(1)))

allocate(Mag_x_1(isize(3),isize(2),isize(1)))
allocate(Mag_x_2(isize(3),isize(2),isize(1)))
allocate(Mag_y_1(isize(3),isize(2),isize(1)))
allocate(Mag_y_2(isize(3),isize(2),isize(1)))
allocate(Mag_z_1(isize(3),isize(2),isize(1)))
allocate(Mag_z_2(isize(3),isize(2),isize(1)))

allocate(Mag_x_1_k(osize(3),osize(2),osize(1)))
allocate(Mag_x_2_k(osize(3),osize(2),osize(1)))
allocate(Mag_y_1_k(osize(3),osize(2),osize(1)))
allocate(Mag_y_2_k(osize(3),osize(2),osize(1)))
allocate(Mag_z_1_k(osize(3),osize(2),osize(1)))
allocate(Mag_z_2_k(osize(3),osize(2),osize(1)))

allocate(d_Bx_dy(isize(3),isize(2),isize(1)))
allocate(d_By_dx(isize(3),isize(2),isize(1)))
allocate(d_Bx_dz(isize(3),isize(2),isize(1)))
allocate(d_By_dz(isize(3),isize(2),isize(1)))
allocate(d_Bz_dx(isize(3),isize(2),isize(1)))
allocate(d_Bz_dy(isize(3),isize(2),isize(1)))

allocate(d_Bx_dy_k(osize(3),osize(2),osize(1)))
allocate(d_By_dx_k(osize(3),osize(2),osize(1)))
allocate(d_Bx_dz_k(osize(3),osize(2),osize(1)))
allocate(d_By_dz_k(osize(3),osize(2),osize(1)))
allocate(d_Bz_dx_k(osize(3),osize(2),osize(1)))
allocate(d_Bz_dy_k(osize(3),osize(2),osize(1)))

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
!Only rank 0 does IO
if (myrank.eq.0) then
open(unit=5,file='System_information.dat',status='unknown')
!open(unit=100,file='omegaval-inplace.dat',status='unknown')
!open(unit=15,file='Initial_Grid_Data.dat',status='unknown')
!open(unit=20,file='INPUT.dat',status='old')  ! This input file is the file generated from Vorticity code. 
!open(unit=25,file='Initial_Grid_Data_Reproduced.dat',status='unknown')
!open(unit=30,file='Initial_Energy_Spectra.dat',status='unknown')
open(unit=35,file='Energy_Spectra.dat',status='unknown')
open(unit=40,file='Energy.dat',status='unknown')

endif
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
!time_max = 0.10d0
!time_max = 0.00020d0 !0.00002d0
!time_max = 0.020d0 !0.00002d0
time_max = 200.0d0 !0.00002d0
!time_max = 1.0d0 !0.00002d0
dt = 0.00010d0    

! Ratio of Specific Heats/ Adiabetic Index/ Heat Capacity Ratio/ Poisson Constant.
spheat = 1.0d0 !5.0d0/3.0d0

! Co-efficient of Viscosity.
mu = 0.20d0

! Co-efficient of Resistivity.
eta = 0.005d0

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
if (myrank.eq.0) then
write(5,*) "Sound Speed=", CS, "Initial Velocity=", U0
write(5,*) "Alfven Speed=", VA, "Initial Magnetic Field=", B0
endif
! Forcing Length Scale.
kf = 1.0d0

 A = 1.0d0
 B = 1.0d0
 C = 1.0d0

! Raynold's Number.
Re = 450.0d0!ms*rho0*U0*Lx/mu; if (myrank.eq.0) write(5,*) "Raynold's Number =", Re
Rm = 450.0d0 !U0*Lx/eta; if (myrank.eq.0) write(5,*) "Magnetic Raynold's Number =", Rm
PM = mu/eta; if (myrank.eq.0) write(5,*) "Prandtl Number =", PM

call srand(seed) !Naga: need to call this again here
do i = 1, Np
xp(i) = rand(0)*Lx
yp(i) = rand(0)*Ly
zp(i) = rand(0)*Lz
enddo
do i = 1, Np/nprocs
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
    enddo
  enddo
enddo

do k = 1, Nz/nprocs
  do j = 1, Ny
    do i = 1, Nx
      ! Initial Density Distribution.
      rho(i,j,k) = rho0
      ! Initial Velocity Distribution.
      ux(i,j,k) = 0.0d0 !U0 * ( A*dsin(kf*z(k+((Nz/nprocs)*myrank))) + C*dcos(kf*y(j)) )
      uy(i,j,k) = 0.0d0 !U0 * ( B*dsin(kf*x(i)) + A*dcos(kf*z(k+((Nz/nprocs)*myrank))) )
      uz(i,j,k) = U0 !* ( C*dsin(kf*y(j)) + B*dcos(kf*x(i)) )
      ! Initial Pressure Distribution.
      P(i,j,k) = P0
      ! Initial Magnetic Field Profile.
      Bx(i,j,k) = 0.0d0 + 0.10d0 * ( A*dsin(kf*z(k+((Nz/nprocs)*myrank))) + C*dcos(kf*y(j)) )
      By(i,j,k) = 0.0d0 + 0.10d0 * ( B*dsin(kf*x(i)) + A*dcos(kf*z(k+((Nz/nprocs)*myrank))) )
      Bz(i,j,k) = B0    + 0.10d0 * ( C*dsin(kf*y(j)) + B*dcos(kf*x(i)) )
      ! Initial Energy Distribution.
      E(i,j,k) = P(i,j,k) + 0.50d0*rho(i,j,k) * (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)&
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
!if (myrank.eq.0) then
! close(15)
!endif
write(*,*) "Rank: ",myrank, ":After initial grid generation"

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
accffttimer = 0.0
call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
if (myrank.eq.0)  then 
t1 =  MPI_WTIME() 
endif

!!$acc data copy(rho, ux, uy, uz, rho_ux, rho_uy, rho_uz, P, E, Bx, By, Bz,xp,yp,zp,uxp,uyp,uzp,vxp,vyp,vzp) & 
!!$acc & create(rho_k,ux_k,uy_k,uz_k,rho_ux_k,rho_uy_k,rho_uz_k,P_k,E_k,Bx_k,By_k,Bz_k) &
!!$acc & create(uxfull,uyfull,uzfull,Bxfull,Byfull,Bzfull) &
!!$acc & copyout(Bpx,Bpy,Bpz) &
!!$acc & copyout(omega_x, omega_y, omega_z) &
!!$acc & create(omega_x_k,omega_y_k,omega_z_k) &
!!$acc & create(Fx, Fy, Fz,Fx_k,Fy_k,Fz_k,d_ux_dx_k,d_uy_dy_k,d_uz_dz_k) &
!!$acc & create(d_ux_dx, d_uy_dy, d_uz_dz, d_Bx_dy, d_Bx_dz, d_By_dx, d_By_dz, d_Bz_dx, d_Bz_dy) &
!!$acc & create(d_Bx_dy_k,d_By_dx_k,d_Bx_dz_k,d_By_dz_k,d_Bz_dx_k,d_Bz_dy_k) &
!!$acc & create(Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3, Energy_x,Energy_y,Energy_z,E_Visc, Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2) &
!!$acc & create(Mom_x_1_k,Mom_x_2_k,Mom_x_3_k,Mom_y_1_k,Mom_y_2_k,Mom_y_3_k,Mom_z_1_k,Mom_z_2_k,Mom_z_3_k, Energy_x_k,Energy_y_k,Energy_z_k,E_Visc_k) &
!!$acc & create(Mag_x_1_k,Mag_x_2_k,Mag_y_1_k,Mag_y_2_k,Mag_z_1_k,Mag_z_2_k) &
!!$acc & copyin(x, y, z) &
!!$acc & copy(rho_k_new, rho_ux_k_new, rho_uy_k_new, rho_uz_k_new, E_k_new, Bx_k_new, By_k_new, Bz_k_new) &
!!$acc & copy(d_rho_k_dt_old, d_rho_ux_k_dt_old, d_rho_uy_k_dt_old, d_rho_uz_k_dt_old, d_E_k_dt_old, d_Bx_k_dt_old, d_By_k_dt_old, d_Bz_k_dt_old) &
!!$acc & copy(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new) &
!!$acc & create(jx_k,jy_k,jz_k,Ax_k,Ay_k,Az_k,div_B_k) &
!!$acc & copy(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new) copyout(div_B,B2,omega2,Ek,Bk, Ax,Ay,Az,jx,jy,jz,A2,j2,u2) 

 ! Inside this region the device data pointer will be used
 !$acc host_data use_device(rho,  ux, uy, uz, rho_ux, rho_uy, rho_uz, P, E, Bx, By, Bz) &
!$acc & use_device(rho_k,ux_k,uy_k,uz_k,rho_ux_k,rho_uy_k,rho_uz_k,P_k,E_k,Bx_k,By_k,Bz_k)

   call accfftCreatePlan3DR2CGPU(C_LOC(gridSize), C_LOC(rho), C_LOC(rho_k), &
                        comm, fftPlan)

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho), C_LOC(rho_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(ux), C_LOC(ux_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(uy), C_LOC(uy_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(uz), C_LOC(uz_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho_ux), C_LOC(rho_ux_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho_uy), C_LOC(rho_uy_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho_uz), C_LOC(rho_uz_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(P), C_LOC(P_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(E), C_LOC(E_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Bx), C_LOC(Bx_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(By), C_LOC(By_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Bz), C_LOC(Bz_k), &
                       C_LOC(accffttimer))

 !$acc end host_data

! Evaluate Initial Vorticity Spectra.
!$acc parallel firstprivate(Lx, Ly, Lz,nprocs,myrank) present(ux_k, uy_k, uz_k, omega_x_k, omega_y_k, omega_z_k)
!$acc loop collapse(3)
do k = 1, Nz
  do j = 1, Ny/nprocs 
    do i = 1, Nx/2+1
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      jglobal = j + myrank*(Ny/nprocs)
      if (jglobal <= Ny/2) then
        ky = 2.0d0*pi*dfloat(jglobal-1)/Ly
      else
        ky = 2.0d0*pi*dfloat((jglobal-1)-Ny)/Ly
      endif
      if (k <= Nz/2) then
        kz = 2.0d0*pi*dfloat(k-1)/Lz
      else
        kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      endif

      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)

!      omega_x(2*i-1,j,k) = (-1.0d0)*ky*uz(2*i,j,k) - (-1.0d0)*kz*uy(2*i,j,k)
!      omega_x(2*i,j,k) = (1.0d0)*ky*uz(2*i-1,j,k) - (1.0d0)*kz*uy(2*i-1,j,k)
!      omega_y(2*i-1,j,k) = (-1.0d0)*kx*uz(2*i,j,k) - (-1.0d0)*kz*ux(2*i,j,k)
!      omega_y(2*i,j,k) = (1.0d0)*kx*uz(2*i-1,j,k) - (1.0d0)*kz*ux(2*i-1,j,k)
!      omega_z(2*i-1,j,k) = (-1.0d0)*kx*uy(2*i,j,k) - (-1.0d0)*ky*ux(2*i,j,k)
!      omega_z(2*i,j,k) = (1.0d0)*kx*uy(2*i-1,j,k) - (1.0d0)*ky*ux(2*i-1,j,k)
      
!      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
!      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
!      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo
!$acc end parallel
 !$acc host_data use_device(omega_x, omega_y, omega_z,omega_x_k,omega_y_k,omega_z_k)

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_x_k_dum, omega_x, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(omega_x_k), C_LOC(omega_x), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_y_k_dum, omega_y, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(omega_y_k), C_LOC(omega_y), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_z_k_dum, omega_z, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(omega_z_k), C_LOC(omega_z), &
                       C_LOC(accffttimer))

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
do k = 1, Nz/nprocs
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
call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
write(*,*) "rank: ",myrank, ":Before loop"
! Time Loop Starts.
do time = time_min,time_max,dt

write(*,*) "rank: ",myrank, ":Iteration"
t = nint(time/dt) - dint(time_min/dt)

! Online Check of the Progress.
!if (mod(t,int(1.0/dt)) == 0) then
  !print*, t/int(1.0/dt)
!endif

!======================================================================================

write(*,*) "Before derive"

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

!$acc host_data use_device(rho, rho_ux, rho_uy, rho_uz, E, Bx, By, Bz) &
!$acc & use_device(rho_k, rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k)
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_k_dum, rho, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(rho_k), C_LOC(rho), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_ux_k_dum, rho_ux, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(rho_ux_k), C_LOC(rho_ux), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uy_k_dum, rho_uy, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(rho_uy_k), C_LOC(rho_uy), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uz_k_dum, rho_uz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(rho_uz_k), C_LOC(rho_uz), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, E_k_dum, E, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(E_k), C_LOC(E), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bx_k_dum, Bx, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(Bx_k), C_LOC(Bx), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, By_k_dum, By, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(By_k), C_LOC(By), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bz_k_dum, Bz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(Bz_k), C_LOC(Bz), &
                       C_LOC(accffttimer))

!$acc end host_data

write(*,*) "Rank: ",myrank, ":Derive 1"
!if (mod(dfloat(t),0.0010d0/dt)==0.0d0) then 
 A = 0.10d0
 B = 0.10d0
 C = 0.10d0
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
do k = 1,Nz/nprocs
  do j = 1,Ny
    do i = 1,Nx
     kglobal = k+(Nz/nprocs)*myrank
!    x(i)=0.0d0+real(i-1)*dx
!    y(j)=0.0d0+real(j-1)*dy
!    z(kglobal)=0.0d0+real(kglobal-1)*dz 
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
    Fx(i,j,k) = rho(i,j,k) * (A*dsin(kf*z(kglobal)) + C*dcos(kf*y(j)))
    Fy(i,j,k) = rho(i,j,k) * (B*dsin(kf*x(i)) + A*dcos(kf*z(kglobal)))
    Fz(i,j,k) = rho(i,j,k) * (C*dsin(kf*y(j)) + B*dcos(kf*x(i)))
  
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
write(*,*) "Rank: ",myrank, ":Derive 2"


!$acc host_data use_device(ux, uy,uz, Fx, Fy, Fz, ux_k,uy_k,uz_k,Fx_k,Fy_k,Fz_k,Bx,By,Bz,Bx_k,By_k,Bz_k,rho_ux_k,rho_uy_k,rho_uz_k,rho_ux,rho_uy,rho_uz,rho_k,rho,E_k,E)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(ux), C_LOC(ux_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(uy), C_LOC(uy_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(uz), C_LOC(uz_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fx, Fx_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Fx), C_LOC(Fx_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fy, Fy_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Fy), C_LOC(Fy_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fz, Fz_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Fz), C_LOC(Fz_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Bx), C_LOC(Bx_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(By), C_LOC(By_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Bz), C_LOC(Bz_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho_ux), C_LOC(rho_ux_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho_uy), C_LOC(rho_uy_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho_uz), C_LOC(rho_uz_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho), C_LOC(rho_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(E), C_LOC(E_k), &
                       C_LOC(accffttimer))
!$acc end host_data

write(*,*) "Rank: ",myrank,":Derive 3"
! Evaluate derivatives of Velocity and Magnetic Field.

!!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,Bx_k,By_k),&
!!$OMP & SHARED(i_kx_ux_k,i_ky_uy_k,i_kz_uz_k),&
!!$OMP & SHARED(i_ky_Bx_k,i_kz_Bx_k,i_kx_By_k,i_kz_By_k,i_kx_Bz_k,i_ky_Bz_k) PRIVATE(i,j,k,kx,ky,kz)
!!$OMP DO
!$acc parallel firstprivate(Lx, Ly, Lz) present(d_ux_dx_k, d_uy_dy_k, d_uz_dz_k, d_Bx_dy_k, d_Bx_dz_k, d_By_dx_k, d_By_dz_k, d_Bz_dx_k, d_Bz_dy_k) &
!$acc & present(ux_k, uy_k, uz_k, Bx_k, By_k, Bz_k) 
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny/nprocs
    do i = 1,Nx/2+1
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      jglobal = j + myrank*(Ny/nprocs)
      if (jglobal <= Ny/2) then
        ky = 2.0d0*pi*dfloat(jglobal-1)/Ly
      else
        ky = 2.0d0*pi*dfloat((jglobal-1)-Ny)/Ly
      endif
      if (k <= Nz/2) then
        kz = 2.0d0*pi*dfloat(k-1)/Lz
      else
        kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      endif

      d_ux_dx_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
      d_uy_dy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k) 
      d_uz_dz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k) 
      
      d_Bx_dy_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      d_Bx_dz_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
      d_By_dx_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
      d_By_dz_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
      d_Bz_dx_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      d_Bz_dy_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)

    enddo 
  enddo
enddo
!$acc end parallel
!!$OMP END DO
!!$OMP END PARALLEL
write(*,*) "Rank: ",myrank,":Derive 4"

!$acc host_data use_device(d_ux_dx, d_uy_dy, d_uz_dz, d_Bx_dy, d_Bx_dz, d_By_dx, d_By_dz, d_Bz_dx, d_Bz_dy, Bx,By,Bz) &
!$acc & use_device(d_ux_dx_k, d_uy_dy_k, d_uz_dz_k, d_Bx_dy_k, d_Bx_dz_k, d_By_dx_k,d_By_dz_k, d_Bz_dx_k, d_Bz_dy_k, Bx_k,By_k,Bz_k)
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_ux_k, d_ux_dx, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(d_ux_dx_k), C_LOC(d_ux_dx), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_uy_k, d_uy_dy, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(d_uy_dy_k), C_LOC(d_uy_dy), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_uz_k, d_uz_dz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(d_uz_dz_k), C_LOC(d_uz_dz), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_Bx_k, d_Bx_dy, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(d_Bx_dy_k), C_LOC(d_Bx_dy), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_Bx_k, d_Bx_dz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(d_Bx_dz_k), C_LOC(d_Bx_dz), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_By_k, d_By_dx, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(d_By_dx_k), C_LOC(d_By_dx), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_By_k, d_By_dz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(d_By_dz_k), C_LOC(d_By_dz), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_Bz_k, d_Bz_dx, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(d_Bz_dx_k), C_LOC(d_Bz_dx), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_Bz_k, d_Bz_dy, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(d_Bz_dy_k), C_LOC(d_Bz_dy), &
                       C_LOC(accffttimer))

!$acc end host_data
write(*,*) "rank:",myrank, ":Derive 5"

!Naga - with accfft, it seems like even with outof place fft, input is changed ,
!in C2R transforms. So, recomputing the transforms as needed

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
do k = 1,Nz/nprocs
  do j = 1,Ny
    do i = 1,Nx
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
write(*,*) "Rank: ",myrank,":Derive 6"


!!$OMP END DO
!!$OMP END PARALLEL


!$acc host_data use_device(Mom_x_1, Mom_x_2, Mom_x_3, Mom_y_1, Mom_y_2, Mom_y_3,Mom_z_1, Mom_z_2, Mom_z_3, Energy_x, Energy_y, Energy_z, E_Visc, Mag_x_1, Mag_x_2, Mag_y_1, Mag_y_2) &
!$acc & use_device(Mag_z_1, Mag_z_2, Mom_x_1_k, Mom_x_2_k, Mom_x_3_k, Mom_y_1_k,Mom_y_2_k, Mom_y_3_k, Mom_z_1_k, Mom_z_2_k, Mom_z_3_k, Energy_x_k, Energy_y_k, Energy_z_k, E_Visc_k, Mag_x_1_k, Mag_x_2_k, Mag_y_1_k, Mag_y_2_k, Mag_z_1_k,Mag_z_2_k)


  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_1, Mom_x_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mom_x_1), C_LOC(Mom_x_1_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_2, Mom_x_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mom_x_2), C_LOC(Mom_x_2_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_3, Mom_x_3_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mom_x_3), C_LOC(Mom_x_3_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_1, Mom_y_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mom_y_1), C_LOC(Mom_y_1_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_2, Mom_y_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mom_y_2), C_LOC(Mom_y_2_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_3, Mom_y_3_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mom_y_3), C_LOC(Mom_y_3_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_1, Mom_z_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mom_z_1), C_LOC(Mom_z_1_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_2, Mom_z_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mom_z_2), C_LOC(Mom_z_2_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_3, Mom_z_3_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mom_z_3), C_LOC(Mom_z_3_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_x, Energy_x_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Energy_x), C_LOC(Energy_x_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_y, Energy_y_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Energy_y), C_LOC(Energy_y_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_z, Energy_z_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Energy_z), C_LOC(Energy_z_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, E_Visc, E_Visc_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(E_Visc), C_LOC(E_Visc_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_x_1, Mag_x_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mag_x_1), C_LOC(Mag_x_1_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_x_2, Mag_x_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mag_x_2), C_LOC(Mag_x_2_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_y_1, Mag_y_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mag_y_1), C_LOC(Mag_y_1_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_y_2, Mag_y_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mag_y_2), C_LOC(Mag_y_2_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_z_1, Mag_z_1_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mag_z_1), C_LOC(Mag_z_1_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_z_2, Mag_z_2_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Mag_z_2), C_LOC(Mag_z_2_k), &
                       C_LOC(accffttimer))

  !call dfftw_destroy_plan_ (plan_backward)
  !call dfftw_destroy_plan_ (plan_forward)
!$acc end host_data


!!$OMP & SHARED(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k),&
!!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,E_Visc_k)&
!!$OMP & SHARED(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k),&
!!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k,nu,Fx_k,Fy_k,Fz_k),&
!!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k),&
!!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new),&
!!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new) PRIVATE(i,j,k,kx,ky,kz)
!!$OMP DO 
write(*,*) "Rank: ",myrank, ":Before derive"

!$acc parallel firstprivate(mu, eta,Lx,Ly,Lz) private(i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k) &
!$acc & private(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k) present(rho_ux_k,rho_uy_k,rho_uz_k,ux_k,uy_k,uz_k) &
!$acc & private(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k,i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k) &
!$acc & private(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k) &
!$acc & private(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k) &
!$acc & private(kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k) &
!$acc & present(Mom_x_1_k,Mom_x_2_k,Mom_x_3_k,Mom_y_1_k,Mom_y_2_k,Mom_y_3_k) &
!$acc & present(Mom_z_1_k,Mom_z_2_k,Mom_z_3_k,Mag_x_1_k,Mag_x_2_k,Mag_y_1_k,Mag_y_2_k,Mag_z_1_k,Mag_z_2_k) &
!$acc & present(E_Visc_k,Fx_k,Fy_k,Fz_k,Energy_x_k,Energy_y_k,Energy_z_k,Bx_k,By_k,Bz_k) &
!$acc & present(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new) &
!$acc & present(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new) 
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny/nprocs
    do i = 1,Nx/2+1
    
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    jglobal = j + myrank*(Ny/nprocs)
    if (jglobal <= Ny/2) then
      ky = 2.0d0*pi*dfloat(jglobal-1)/Ly
    else
      ky = 2.0d0*pi*dfloat((jglobal-1)-Ny)/Ly
    endif
    if (k <= Nz/2) then
      kz = 2.0d0*pi*dfloat(k-1)/Lz
    else
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
    endif

    !if (kx >= Nx/3 .and. ky >= Ny/3 .and. kz >= Nz/3) then!dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
     if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
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
        E_Visc_k(i,j,k) = 0.0d0
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
        i_kx_rho_ux_k = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k) 
        i_ky_rho_uy_k = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k)
        i_kz_rho_uz_k = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k)
      i_kx_Mom_x_1_k = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)
      i_ky_Mom_x_2_k = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
      i_kz_Mom_x_3_k = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)
      i_kx_Mom_y_1_k = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
      i_ky_Mom_y_2_k = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
      i_kz_Mom_y_3_k = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k)
      i_kx_Mom_z_1_k = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)
      i_ky_Mom_z_2_k = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)
      i_kz_Mom_z_3_k = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)
      kx2_ux_k = kx*kx*ux_k(i,j,k)
      ky2_ux_k = ky*ky*ux_k(i,j,k)
      kz2_ux_k = kz*kz*ux_k(i,j,k)
      kx2_uy_k = kx*kx*uy_k(i,j,k)
      ky2_uy_k = ky*ky*uy_k(i,j,k)
      kz2_uy_k = kz*kz*uy_k(i,j,k)
      kx2_uz_k = kx*kx*uz_k(i,j,k)
      ky2_uz_k = ky*ky*uz_k(i,j,k)
      kz2_uz_k = kz*kz*uz_k(i,j,k)
      i_kx_Energy_x_k = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
      i_ky_Energy_y_k = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
      i_kz_Energy_z_k = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)
      i_ky_Mag_x_1_k = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
      i_kz_Mag_x_2_k = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
      i_kx_Mag_y_1_k = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
      i_kz_Mag_y_2_k = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
      i_kx_Mag_z_1_k = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
      i_ky_Mag_z_2_k = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)
    
      kx2_Bx_k = kx*kx*Bx_k(i,j,k)
      ky2_Bx_k = ky*ky*Bx_k(i,j,k)
      kz2_Bx_k = kz*kz*Bx_k(i,j,k)
      kx2_By_k = kx*kx*By_k(i,j,k)
      ky2_By_k = ky*ky*By_k(i,j,k)
      kz2_By_k = kz*kz*By_k(i,j,k)
      kx2_Bz_k = kx*kx*Bz_k(i,j,k)
      ky2_Bz_k = ky*ky*Bz_k(i,j,k)
      kz2_Bz_k = kz*kz*Bz_k(i,j,k)
    endif




      ! Density Equation.
      d_rho_k_dt_new(i,j,k) = - ( i_kx_rho_ux_k + i_ky_rho_uy_k + i_kz_rho_uz_k)
    
      ! Momentum Equation.
      d_rho_ux_k_dt_new(i,j,k) = - ( i_kx_Mom_x_1_k + i_ky_Mom_x_2_k + i_kz_Mom_x_3_k ) &
                                 - ( kx2_ux_k + ky2_ux_k + kz2_ux_k ) / 450.0d0 !+ Fx_k(i,j,k)
    
      d_rho_uy_k_dt_new(i,j,k) = - ( i_kx_Mom_y_1_k + i_ky_Mom_y_2_k + i_kz_Mom_y_3_k ) &
                                 - ( kx2_uy_k + ky2_uy_k + kz2_uy_k ) / 450.0d0 !+ Fy_k(i,j,k)
    
      d_rho_uz_k_dt_new(i,j,k) = - ( i_kx_Mom_z_1_k + i_ky_Mom_z_2_k + i_kz_Mom_z_3_k ) &
                                 - ( kx2_uz_k + ky2_uz_k + kz2_uz_k ) / 450.0d0 !+ Fz_k(i,j,k) 
    
      ! Energy Equation.
      d_E_k_dt_new(i,j,k) = - ( i_kx_Energy_x_k + i_ky_Energy_y_k + i_kz_Energy_z_k ) &
                            + mu * E_Visc_k(i,j,k) !/ mu
    
      ! Magnetic Field Equation.
      d_Bx_k_dt_new(i,j,k) = + ( i_ky_Mag_x_1_k + i_kz_Mag_x_2_k ) &
                             - ( kx2_Bx_k + ky2_Bx_k + kz2_Bx_k ) / 450.0d0
                           
      d_By_k_dt_new(i,j,k) = - ( i_kx_Mag_y_1_k - i_kz_Mag_y_2_k ) &
                             - ( kx2_By_k + ky2_By_k + kz2_By_k ) / 450.0d0
                           
      d_Bz_k_dt_new(i,j,k) = - ( i_kx_Mag_z_1_k + i_ky_Mag_z_2_k ) &                       
                             - ( kx2_Bz_k + ky2_Bz_k + kz2_Bz_k ) / 450.0d0
    enddo                       
  enddo
enddo  
!$acc end parallel
!!$acc end data
!!$OMP END DO
!!$OMP END PARALLEL
write(*,*) "Rank: ",myrank, ":After derive"

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

!$acc parallel present(rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, d_rho_ux_k_dt_new, d_rho_uy_k_dt_new, d_rho_uz_k_dt_new, d_rho_ux_k_dt_old) &
!$acc & present(d_rho_uy_k_dt_old, d_rho_uz_k_dt_old, d_E_k_dt_new, d_E_k_dt_old, d_Bx_k_dt_new, d_By_k_dt_new, d_Bz_k_dt_new, d_Bx_k_dt_old, d_By_k_dt_old, d_Bz_k_dt_old) &
!$acc & present(rho_k_new, rho_ux_k_new, rho_uy_k_new, rho_uz_k_new, E_k_new, Bx_k_new, By_k_new, Bz_k_new)

!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny/nprocs
    do i = 1,Nh
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
  enddo
enddo
!$acc end parallel
!!$acc end data
write(*,*) "Rank: ",myrank, ":After ab"
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
!$acc & present(d_rho_k_dt_old, d_rho_ux_k_dt_old, d_rho_uy_k_dt_old, d_rho_uz_k_dt_old, d_E_k_dt_old, d_Bx_k_dt_old, d_By_k_dt_old, d_Bz_k_dt_old, rho_k, rho_ux_k, rho_uy_k, rho_uz_k, E_k, Bx_k, By_k, Bz_k)

!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny/nprocs
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
    
      rho_k(i,j,k) = rho_k_new(i,j,k)
   
      rho_ux_k(i,j,k) = rho_ux_k_new(i,j,k)
      rho_uy_k(i,j,k) = rho_uy_k_new(i,j,k)
      rho_uz_k(i,j,k) = rho_uz_k_new(i,j,k)
 
      E_k(i,j,k) = E_k_new(i,j,k)

      Bx_k(i,j,k) = Bx_k_new(i,j,k)
      By_k(i,j,k) = By_k_new(i,j,k)
      Bz_k(i,j,k) = Bz_k_new(i,j,k)
 
    
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
write(*,*) "Rank: ",myrank, ":Final loop1"

!!$OMP END DO
!!$OMP END PARALLEL 

! Evaluate Divergence of B in Spectral Space.

!!$OMP PARALLEL SHARED(Lx,Ly,Lz,Bx_k,By_k,Bz_k,div_B_k) PRIVATE(i,j,k,kx,ky,kz)
!!$OMP DO
!$acc parallel firstprivate(Lx, Ly, Lz) present(jx_k,jy_k,jz_k,Ax_k,Ay_k,Az_k,Bx_k,By_k,Bz_k,div_B_k) 
!$acc loop collapse(3)
do k = 1,Nz
  do j = 1,Ny/nprocs
    do i = 1,Nx/2+1
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      jglobal = j + myrank*(Ny/nprocs)

      if (jglobal <= Ny/2) then
        ky = 2.0d0*pi*dfloat(jglobal-1)/Ly
      else
        ky = 2.0d0*pi*dfloat((jglobal-1)-Ny)/Ly
      endif
      if (k <= Nz/2) then
        kz = 2.0d0*pi*dfloat(k-1)/Lz
      else
        kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      endif
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
!      jx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k) - (0.0d0,1.0d0)*kz*By_k(i,j,k)
!      jy_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k) - (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
!      jz_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k) - (0.0d0,1.0d0)*ky*Bx_k(i,j,k)

!        if (i == 1 .and. jglobal == 1 .and. k == 1) then
!          Ax_k(i,j,k) = jx_k(i,j,k)
!          Ay_k(i,j,k) = jy_k(i,j,k)
!          Az_k(i,j,k) = jz_k(i,j,k)
!
!        else  
!          Ax_k(i,j,k) = jx_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
!          Ay_k(i,j,k) = jy_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
!          Az_k(i,j,k) = jz_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
!
!        endif
    enddo
  enddo
enddo
!$acc end parallel
write(*,*) "Rank: ",myrank,":Final loop2"
!!$OMP END DO
!!$OMP END PARALLEL 
!$acc host_data use_device(rho, rho_ux, rho_uy, rho_uz, E, Bx, By, Bz, div_B, jx, jy, jz, Ax, Ay, Az) &
!$acc & use_device(rho_k, rho_ux_k, rho_uy_k, rho_uz_k, E_k, Bx_k, By_k, Bz_k, div_B_k, jx_k, jy_k, jz_k, Ax_k, Ay_k, Az_k) 

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_k_dum, rho, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(rho_k), C_LOC(rho), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_ux_k_dum, rho_ux, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(rho_ux_k), C_LOC(rho_ux), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uy_k_dum, rho_uy, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(rho_uy_k), C_LOC(rho_uy), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uz_k_dum, rho_uz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(rho_uz_k), C_LOC(rho_uz), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, E_k_dum, E, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(E_k), C_LOC(E), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bx_k_dum, Bx, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(Bx_k), C_LOC(Bx), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, By_k_dum, By, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(By_k), C_LOC(By), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bz_k_dum, Bz, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(Bz_k), C_LOC(Bz), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, div_B_k, div_B, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(div_B_k), C_LOC(div_B), &
                       C_LOC(accffttimer))

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, jx_k, jx, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(jx_k), C_LOC(jx), &
                       C_LOC(accffttimer))

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, jy_k, jy, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(jy_k), C_LOC(jy), &
                       C_LOC(accffttimer))

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, jz_k, jz, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(jz_k), C_LOC(jz), &
                       C_LOC(accffttimer))

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Ax_k, Ax, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(Ax_k), C_LOC(Ax), &
                       C_LOC(accffttimer))

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Ay_k, Ay, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(Ay_k), C_LOC(Ay), &
                       C_LOC(accffttimer))

!  call dfftw_init_threads(iret)
!  call dfftw_plan_with_nthreads(thread_num)
!  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Az_k, Az, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_backward)
!  call dfftw_destroy_plan_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(Az_k), C_LOC(Az), &
                       C_LOC(accffttimer))

 !$acc end host_data
write(*,*) "Rank: ",myrank, ":Final 3"

!!$OMP PARALLEL SHARED(spheat,rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz,div_B), &
!!$OMP & SHARED(ux,uy,uz,P) PRIVATE(i,j,k)
!!$OMP DO

!$acc parallel firstprivate(spheat) present(rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz,div_B, ux,uy,uz,jx,jy,jz,Ax,Ay,Az,u2,j2,A2,P,B2)
!$acc loop collapse(3)
do k = 1,Nz/nprocs
  do j = 1,Ny
    do i = 1,Nx
      ! FFTW Normalisation.
!      jx(i,j,k) = jx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      jy(i,j,k) = jy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      jz(i,j,k) = jz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      
!      Ax(i,j,k) = Ax(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      Ay(i,j,k) = Ay(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      Az(i,j,k) = Az(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      
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
!      j2(i,j,k) = jx(i,j,k)*jx(i,j,k) + jy(i,j,k)*jy(i,j,k) + jz(i,j,k)*jz(i,j,k)
!      A2(i,j,k) = Ax(i,j,k)*Ax(i,j,k) + Ay(i,j,k)*Ay(i,j,k) + Az(i,j,k)*Az(i,j,k)
  
      ! Keep Backup of the Arrays for FFTW.
!      ux_dum(i,j,k) = ux(i,j,k)
!      uy_dum(i,j,k) = uy(i,j,k)
!      uz_dum(i,j,k) = uz(i,j,k)
  
      ! Evaluate Pressure
      P(i,j,k) = CS*CS*rho(i,j,k)!( spheat - 1.0d0 ) * ( E(i,j,k) - 0.50d0 * &
                 !( rho_ux(i,j,k)*ux(i,j,k)+rho_uy(i,j,k)*uy(i,j,k)+rho_uz(i,j,k)*uz(i,j,k) - B2(i,j,k) ) )
    enddo
  enddo
enddo   
!$acc end parallel
write(*,*) "Rank: ",myrank, ":Final 4"

!==================================================================================================================

!Every process needs to have the complete ux,uy,uz,Bx,By,Bz before the following
!loop

!!$acc update self(ux,uy,uz,Bx,By,Bz)
!$acc host_data use_device(ux,uy,uz,uxfull,uyfull,uzfull,Bx,By,Bz,Bxfull,Byfull,Bzfull)
call MPI_ALLGATHER(ux, Nx*Ny*(Nz/nprocs), MPI_DOUBLE, uxfull, &
             Nx*Ny*(Nz/nprocs), MPI_DOUBLE, MPI_COMM_WORLD, mpierr)
call MPI_ALLGATHER(uy, Nx*Ny*(Nz/nprocs), MPI_DOUBLE, uyfull, &
             Nx*Ny*(Nz/nprocs), MPI_DOUBLE, MPI_COMM_WORLD, mpierr)
call MPI_ALLGATHER(uz, Nx*Ny*(Nz/nprocs), MPI_DOUBLE, uzfull, &
             Nx*Ny*(Nz/nprocs), MPI_DOUBLE, MPI_COMM_WORLD, mpierr)

call MPI_ALLGATHER(Bx, Nx*Ny*(Nz/nprocs), MPI_DOUBLE, Bxfull, &
             Nx*Ny*(Nz/nprocs), MPI_DOUBLE, MPI_COMM_WORLD, mpierr)
call MPI_ALLGATHER(By, Nx*Ny*(Nz/nprocs), MPI_DOUBLE, Byfull, &
             Nx*Ny*(Nz/nprocs), MPI_DOUBLE, MPI_COMM_WORLD, mpierr)
call MPI_ALLGATHER(Bz, Nx*Ny*(Nz/nprocs), MPI_DOUBLE, Bzfull, &
             Nx*Ny*(Nz/nprocs), MPI_DOUBLE, MPI_COMM_WORLD, mpierr)
!$acc end host_data
! Velocity and Magnetic Field Interpolation from Grid to Tracer Particle...
!!!$acc parallel loop vector_length(256) & 
!!!$acc& present(ux,uy,uz,Bx,By,Bz,x,y,z,uxp,uyp,uzp,Bpx,Bpy,Bpz,x,y,z)
!!$acc parallel loop present(uxfull,uyfull,uzfull,Bxfull,Byfull,Bzfull,xp,yp,zp,uxp,uyp,uzp,Bpx,Bpy,Bpz,x,y,z)

!do i = 1,Np/nprocs

!  iglobal = i + (Np/nprocs)*myrank

!  m = int(xp(iglobal)/dx) + 1

!  n = int(yp(iglobal)/dy) + 1
  
!  q = int(zp(iglobal)/dz) + 1


!if (m==Nx .and. n/=Ny .and. q/=Nz) then

!  uxp(i) = uxfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m,n,q+1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m,n+1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpx(i) = Bxfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m,n,q+1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m,n+1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  uyp(i) = uyfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m,n,q+1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m,n+1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

  
!  Bpy(i) = Byfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m,n,q+1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m,n+1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  uzp(i) = uzfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m,n,q+1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m,n+1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpz(i) = Bzfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m,n,q+1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m,n+1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!end if

!if (m/=Nx .and. n==Ny .and. q/=Nz) then

!  uxp(i) = uxfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m,1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m+1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpx(i) = Bxfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m,1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m+1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  uyp(i) = uyfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m,1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m+1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

  
!  Bpy(i) = Byfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m,1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m+1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  uzp(i) = uzfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m,1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m+1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpz(i) = Bzfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m,1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m+1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!endif

!if (m/=Nx .and. n/=Ny .and. q==Nz) then

!  uxp(i) = uxfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uxfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uxfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uxfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uxfull(m,n,1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m+1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m,n+1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m+1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpx(i) = Bxfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(m,n,1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m+1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m,n+1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m+1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  uyp(i) = uyfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uyfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uyfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uyfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uyfull(m,n,1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m+1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m,n+1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m+1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpy(i) = Byfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Byfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Byfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Byfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Byfull(m,n,1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m+1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m,n+1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m+1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h
  
!  uzp(i) = uzfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uzfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uzfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uzfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uzfull(m,n,1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m+1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m,n+1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m+1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpz(i) = Bzfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(m,n,1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m+1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m,n+1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m+1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!endif

!if (m==Nx .and. n==Ny .and. q/=Nz ) then

!  uxp(i) = uxfull(m,n,q) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m,1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m,n,q+1) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m,1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpx(i) = Bxfull(m,n,q) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m,1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m,n,q+1) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m,1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  uyp(i) = uyfull(m,n,q) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m,1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m,n,q+1) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m,1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpy(i) = Byfull(m,n,q) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m,1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m,n,q+1) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m,1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h
  
!  uzp(i) = uzfull(m,n,q) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m,1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m,n,q+1) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m,1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpz(i) = Bzfull(m,n,q) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m,1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m,n,q+1) * (Lx-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(1,n,q+1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m,1,q+1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(1,1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!endif

!if (m==Nx .and. n/=Ny .and. q==Nz ) then

!  uxp(i) = uxfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uxfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uxfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uxfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uxfull(m,n,1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m,n+1,1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpx(i) = Bxfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(m,n,1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m,n+1,1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  uyp(i) = uyfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uyfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uyfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uyfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uyfull(m,n,1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m,n+1,1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpy(i) = Byfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Byfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Byfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Byfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Byfull(m,n,1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m,n+1,1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h
  
!  uzp(i) = uzfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uzfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uzfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uzfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uzfull(m,n,1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m,n+1,1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpz(i) = Bzfull(m,n,q) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(m,n+1,q) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(m,n,1) * (Lx-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(1,n,1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m,n+1,1) * (Lx-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(1,n+1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!endif

!if (m/=Nx .and. n==Ny .and. q==Nz ) then

!  uxp(i) = uxfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uxfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uxfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uxfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uxfull(m,n,1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m+1,n,1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m,1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m+1,1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpx(i) = Bxfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bxfull(m,n,1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m+1,n,1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m,1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m+1,1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  uyp(i) = uyfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uyfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uyfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uyfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uyfull(m,n,1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m+1,n,1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m,1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m+1,1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpy(i) = Byfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Byfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Byfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Byfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Byfull(m,n,1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m+1,n,1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m,1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m+1,1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h
  
!  uzp(i) = uzfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uzfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           uzfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uzfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           uzfull(m,n,1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m+1,n,1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m,1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m+1,1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpz(i) = Bzfull(m,n,q) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(m+1,n,q) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(m,1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(m+1,1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (Lz-zp(iglobal)) * h + &
!           Bzfull(m,n,1) * (x(m+1)-xp(iglobal)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m+1,n,1) * (xp(iglobal)-x(m)) * (Ly-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m,1,1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m+1,1,1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!endif

!if (m/=Nx .and. n/=Ny .and. q/=Nz) then

!  uxp(i) = uxfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uxfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m,n+1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uxfull(m+1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpx(i) = Bxfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bxfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m,n+1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bxfull(m+1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  uyp(i) = uyfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uyfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m,n+1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uyfull(m+1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpy(i) = Byfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Byfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m,n+1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Byfull(m+1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h
  
!  uzp(i) = uzfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           uzfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m,n+1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           uzfull(m+1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h

!  Bpz(i) = Bzfull(m,n,q) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m+1,n,q) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m,n+1,q) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m+1,n+1,q) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (z(q+1)-zp(iglobal)) * h + &
!           Bzfull(m,n,q+1) * (x(m+1)-xp(iglobal)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m+1,n,q+1) * (xp(iglobal)-x(m)) * (y(n+1)-yp(iglobal)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m,n+1,q+1) * (x(m+1)-xp(iglobal)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h + &
!           Bzfull(m+1,n+1,q+1) * (xp(iglobal)-x(m)) * (yp(iglobal)-y(n)) * (zp(iglobal)-z(q)) * h
!end if

!end do ! i
!!$acc end parallel

!! Boris Algorithm...
!!!$acc parallel loop vector_length(256) present(vxp,vyp,vzp,uxp,uyp,uzp,x,y,z,Bpx,Bpy,Bpz) &
!!$acc parallel loop present(vxp,vyp,vzp,uxp,uyp,uzp,Bpx,Bpy,Bpz,xp,yp,zp) &
!!$acc& private(vxd,vyd,vzd)
 
! do i=1,Np/nprocs
!  tx = 0.5d0*dt*Bpx(i)
!  ty = 0.5d0*dt*Bpy(i)
!  tz = 0.5d0*dt*Bpz(i)

!  iglobal = i + myrank * (Np/nprocs)
  
!  sx = 2.0d0*tx/(1.0d0+tx**2+ty**2+tz**2)
!  sy = 2.0d0*ty/(1.0d0+tx**2+ty**2+tz**2)
!  sz = 2.0d0*tz/(1.0d0+tx**2+ty**2+tz**2)

!  vxd(i) = uxp(i) + (uyp(i)*tz-uzp(i)*ty)
!  vyd(i) = uyp(i) - (uxp(i)*tz-uzp(i)*tx)
!  vzd(i) = uzp(i) + (uxp(i)*ty-uyp(i)*tx)
  
!  vxp(i) =  + (vyd(i)*sz-vzd(i)*sy)
!  vyp(i) = uyp(i) - (vxd(i)*sz-vzd(i)*sx)
!  vzp(i) = uzp(i) + (vxd(i)*sy-vyd(i)*sx)
  
!  xp(iglobal) = xp(iglobal) + vxp(i)*dt
!  yp(iglobal) = yp(iglobal) + vyp(i)*dt
!  zp(iglobal) = zp(iglobal) + vzp(i)*dt
  
  ! Periodic Boundary Condition Implemented...
!  xp(iglobal) = xp(iglobal) - (int(xp(iglobal)/Lx))*Lx
!    if (xp(iglobal) .lt. 0.0d0) then
!    xp(iglobal) = xp(iglobal) + Lx
!    endif
  
!  yp(iglobal) = yp(iglobal) - (int(yp(iglobal)/Ly))*Ly
!    if (yp(iglobal) .lt. 0.0d0) then
!    yp(iglobal) = yp(iglobal) + Ly
!    endif

!  zp(iglobal) = zp(iglobal) - (int(zp(iglobal)/Lz))*Lz
!    if (zp(iglobal) .lt. 0.0d0) then
!    zp(iglobal) = zp(iglobal) + Lz
!    endif

!  enddo

!!$acc end parallel
!!$acc update host(xp,yp,zp,uxp,uyp,uzp)
!  if (mod(float(t),1000.0) == 0.0) then
!  do i = 1,Np/nprocs
!  iglobal = i + myrank * (Np/nprocs)
!  if (mod(iglobal,100) == 0) write(t*10+100+myrank,*) xp(iglobal),yp(iglobal),zp(iglobal),uxp(i),uyp(i),uzp(i)
!  end do
!  endif
  
!==================================================================================================================


write(*,*) "Rank: ",myrank,":Final 4a"

  
!!$OMP END DO
!!$OMP END PARALLEL 

!$acc host_data use_device(ux, uy, uz, ux_k,uy_k,uz_k,Bx,By,Bz,Bx_k,By_k,Bz_k,P,P_k,E,E_k,rho,rho_k,rho_ux,rho_ux_k,rho_uy,rho_uy_k,rho_uz,rho_uz_k)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(ux), C_LOC(ux_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(uy), C_LOC(uy_k), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(uz), C_LOC(uz_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Bx), C_LOC(Bx_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(By), C_LOC(By_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(Bz), C_LOC(Bz_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho), C_LOC(rho_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho_ux), C_LOC(rho_ux_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho_uy), C_LOC(rho_uy_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(rho_uz), C_LOC(rho_uz_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(P), C_LOC(P_k), &
                       C_LOC(accffttimer))
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteR2CGPU(fftPlan, C_LOC(E), C_LOC(E_k), &
                       C_LOC(accffttimer))
! print*, "After fft2 in main" 

!$acc end host_data
write(*,*) "Rank: ",myrank, ":Final 5"

! Evaluate Vorticity in Spectral Space.

!!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,omega_x_k,omega_y_k,omega_z_k),&
!!$OMP & PRIVATE(i,j,k,kx,ky,kz)
!!$OMP DO 
!$acc parallel firstprivate(Lx, Ly, Lz) present(ux_k,uy_k,uz_k,omega_x_k,omega_y_k,omega_z_k)
!$acc loop collapse(3)
do k = 1, Nz
  do j = 1, Ny/nprocs
    do i = 1, Nx/2+1
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      jglobal = j + myrank*(Ny/nprocs)

      if (jglobal <= Ny/2) then
        ky = 2.0d0*pi*dfloat(jglobal-1)/Ly
      else
        ky = 2.0d0*pi*dfloat((jglobal-1)-Ny)/Ly
      endif
      if (k <= Nz/2) then
        kz = 2.0d0*pi*dfloat(k-1)/Lz
      else
        kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      endif
    
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
 
     ! omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
     ! omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
     ! omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo
!$acc end parallel
write(*,*) "Rank: ",myrank, ":Final 6"

!!$OMP END DO
!!$OMP END PARALLEL

!$acc host_data use_device(omega_x, omega_y, omega_z,omega_x_k,omega_y_k,omega_z_k)
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_x_k_dum, omega_x, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(omega_x_k), C_LOC(omega_x), &
                       C_LOC(accffttimer))

  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_y_k_dum, omega_y, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(omega_y_k), C_LOC(omega_y), &
                       C_LOC(accffttimer))
  
  !call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_z_k_dum, omega_z, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_backward)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call accfftExecuteC2RGPU(fftPlan, C_LOC(omega_z_k), C_LOC(omega_z), &
                       C_LOC(accffttimer))

 !$acc end host_data
write(*,*) "Rank: ",myrank,":Final 7"

  !call dfftw_destroy_plan_ (plan_forward)
  !call dfftw_destroy_plan_ (plan_backward)

! FFTW Normalisation and omega^2 Evaluation.

!!$OMP PARALLEL SHARED(omega_x,omega_y,omega_z,omega2) PRIVATE(i,j,k)
!!$OMP DO 
!$acc parallel present(omega_x,omega_y,omega_z,omega2)
!$acc loop collapse(3)
do k = 1,Nz/nprocs
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
write(*,*) "Rank: ",myrank,":Final 8"

!!$OMP END DO
!!$OMP END PARALLEL

! Evaluate Energy Spectra at a Specified Time.
!if (t >= int((time_max-time_max/10.0d0)/dt)) then

!!! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file. !!!
! !!$OMP PARALLEL SHARED(ux_k,uy_k,uz_k) PRIVATE(i,j,k) 
! !!$OMP DO REDUCTION (+:Ek) 

!!$acc parallel loop collapse(3) present(ux_k,uy_k,uz_k,Ek) 

! do k = 1,Nz
!    do j = 1,Ny/nprocs
!      do i = 1,Nx/2+1
!      Ek(i,j,k) = Ek(i,j,k) + sqrt(abs(ux_k(i,j,k))**2 + abs(uy_k(i,j,k))**2 + abs(uz_k(i,j,k))**2)&
!                    /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!        Bk(i,j,k) = Bk(i,j,k) + sqrt(abs(Bx_k(i,j,k))**2 + abs(By_k(i,j,k))**2 + abs(Bz_k(i,j,k))**2)&
!                    /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

!      enddo
!   enddo
!  enddo  
!!$acc end parallel
! !!$OMP END PARALLEL

!endif 
!!$acc host_data use_device(ux,uy,uz,ux_k,uy_k,uz_k,Bx,By,Bz,Bx_k,By_k,Bz_k)

!  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
!  call accfftExecuteC2RGPU(fftPlan, C_LOC(Bx_k), C_LOC(Bx), &
!                       C_LOC(accffttimer))

!  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
!  call accfftExecuteC2RGPU(fftPlan, C_LOC(By_k), C_LOC(By), &
!                       C_LOC(accffttimer))

!  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
!  call accfftExecuteC2RGPU(fftPlan, C_LOC(Bz_k), C_LOC(Bz), &
!                       C_LOC(accffttimer))
!  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
!  call accfftExecuteC2RGPU(fftPlan, C_LOC(ux_k), C_LOC(ux), &
!                       C_LOC(accffttimer))
!  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
!  call accfftExecuteC2RGPU(fftPlan, C_LOC(uy_k), C_LOC(uy), &
!                       C_LOC(accffttimer))
!  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
!  call accfftExecuteC2RGPU(fftPlan, C_LOC(uz_k), C_LOC(uz), &
!                       C_LOC(accffttimer))
!!$acc end host_data
Pressure = 0.0d0; Energy = 0.0d0; T_Energy = 0.0d0; B_Field = 0.0d0
 C0 = 0.0d0; C1 = 0.0d0; C2 = 0.0d0; C3 = 0.0d0; div_B_Tot = 0.0d0
 Hf = 0.0d0; Hm = 0.0d0; HT = 0.0d0; HuB = 0.0d0; HAw = 0.0d0; Rayleigh = 0.0d0
   
! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file.
! !$OMP PARALLEL SHARED(mu_0,t,dt,x,y,z,omega_x,omega_y,omega_z,omega2,rho,ux,uy,uz) PRIVATE(i,j,k)
! !$OMP DO REDUCTION (+:Pressure,Energy,y_Energy,B_Field,C0,C1,C2,C3,div_B_Tot) 

!$acc parallel loop collapse(3) reduction(+:Pressure,Energy,T_Energy,B_Field,C0,C1,C2,C3,div_B_Tot,Hf,Hm,HT,HuB,HAw,Rayleigh) present(P,ux,uy,uz,rho,Bx,By,Bz,u2,B2,Ax,Ay,Az,jx,jy,jz,omega2,div_B)
do k = 1,Nz/nprocs
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
      !Hf = Hf + ( ux(i,j,k)*omega_x(i,j,k) + uy(i,j,k)*omega_y(i,j,k) + uz(i,j,k)*omega_z(i,j,k) ) &
      !          /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      !Hm = Hm + ( Ax(i,j,k)*Bx(i,j,k) + Ay(i,j,k)*By(i,j,k) + Az(i,j,k)*Bz(i,j,k) )/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))  
      !HT = HT + ( (ux(i,j,k)+Ax(i,j,k))*(omega_x(i,j,k)+Bx(i,j,k)) + (uy(i,j,k)+Ay(i,j,k))*(omega_y(i,j,k)+By(i,j,k)) &
      !          + (uz(i,j,k)+Az(i,j,k))*(omega_z(i,j,k)+Bz(i,j,k)) )/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz)) 
      !HuB = HuB + ( ux(i,j,k)*Bx(i,j,k) + uy(i,j,k)*By(i,j,k) + uz(i,j,k)*Bz(i,j,k) )/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      !HAw = HAw + ( Ax(i,j,k)*omega_x(i,j,k) + Ay(i,j,k)*omega_y(i,j,k) + Az(i,j,k)*omega_z(i,j,k) ) &
                !/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))            
      ! Check for Div B = 0
      div_B_Tot = div_B_Tot + div_B(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Rayleigh Quotient.
      !Rayleigh = Rayleigh + (omega2(i,j,k)+0.50d0*j2(i,j,k))/( (u2(i,j,k)+0.50d0*B2(i,j,k)) * dfloat(Nx)*dfloat(Ny)*dfloat(Nz) )
    enddo
  enddo
enddo  
      
!$acc end parallel
!!$acc update self(ux,uy,uz,Bx,By,Bz)
if (mod(float(t),100000.0) == 0.0) then
do i = 1,Nx
 do j = 1,Ny
    do k = 1,Nz/nprocs 
 ! Write Grid Data in State Files.  
        !write(t*10+100+myrank,*) ux(i,j,k),uy(i,j,k),uz(i,j,k)!,Bx(i,j,k),By(i,j,k),Bz(i,j,k)
        write(t+100+myrank,*) ux(i,j,k),uy(i,j,k),uz(i,j,k)
        write(t+110+myrank,*) Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    enddo
 enddo
enddo
endif  
! !$OMP END PARALLEL
!!$acc update self(ux,uy,uz,Bx,By,Bz)
!if (mod(float(t),1000.0) == 0.0) then
!do k = 1,Nz/nprocs,32
!  do j = 1,Ny,32
!    do i = 1,Nx,32
      ! Write Grid Data in State Files.  
        !write(t*10+110+myrank,*) ux(i,j,k),uy(i,j,k),uz(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
!    enddo
!  enddo
!enddo  
!  close(t*10+110+myrank)
!endif
!!$acc update self(u2,B2,omega2,j2,A2)
!if (mod(float(t),1000.0) == 0.0) then
!do k = 1,Nz/nprocs,32
!  do j = 1,Ny,32
!    do i = 1,Nx,32
!        kglobal = k + myrank * (Nz/nprocs)
!        if (i == 1) then
!        write(t*10+120+myrank,*) y(j),z(kglobal),u2(i,j,k),B2(i,j,k),omega2(i,j,k),j2(i,j,k),A2(i,j,k)
!        elseif (j == 1) then
!        write(t*10+130+myrank,*) x(i),z(kglobal),u2(i,j,k),B2(i,j,k),omega2(i,j,k),j2(i,j,k),A2(i,j,k)
!        elseif (kglobal == 1) then
!        write(t*10+140+myrank,*) x(i),y(j),u2(i,j,k),B2(i,j,k),omega2(i,j,k),j2(i,j,k),A2(i,j,k)
!        endif
!    enddo
!  enddo
!enddo  
!  close(t*10+120+myrank)
!  close(t*10+130+myrank)
!  close(t*10+140+myrank)
!endif

!Do mpi reduce on root 0
Energyfinal = 0.0
B_Fieldfinal = 0.0

call MPI_REDUCE(Energy, Energyfinal, 1, MPI_DOUBLE, MPI_SUM, 0, &
      MPI_COMM_WORLD, mpierr) 
call MPI_REDUCE(B_Field, B_Fieldfinal, 1, MPI_DOUBLE, MPI_SUM, 0, &
      MPI_COMM_WORLD, mpierr) 
if (myrank.eq.0) then
!if (mod(float(t),100.0) == 0.0) then
  write(40,*) time,Energyfinal,B_Fieldfinal!,C1,C2,C3,Hf,Hm,HT,HuB,HAw,Rayleigh,div_B_Tot
  call flush(40) 
!endif
endif
if (mod(float(t),100000.0) == 0.0) then
  close(t+100+myrank)
  close(t+110+myrank)
endif

! close(t+100)

enddo ! time

!!$acc end data !======================================================================
!t2 = omp_get_wtime()
call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
if (myrank.eq.0) then 
t2 =  MPI_WTIME() 

!write(5,*) "Time taken for the run =",(t2 - t1)/(60.0d0*60.0d0),"Hours"
write(5,*) "Time taken for the run =",(t2 - t1),"Seconds"
endif
!Collect Ek full array and Bk full array at process 0

!call MPI_GATHER(Ek, Nh*(Ny/nprocs)*Nz, MPI_DOUBLE_COMPLEX, Ekfull, &
!             Nh*(Ny/nprocs)*(Nz), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, mpierr)
!call MPI_GATHER(Bk, Nh*(Ny/nprocs)*Nz, MPI_DOUBLE_COMPLEX, Bkfull, &
!             Nh*(Ny/nprocs)*(Nz), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, mpierr)
!if (myrank.eq.0) then 
!do i = 1,Nx/2+1
!  do j = 1,Ny/2
!    do k = 1,Nz/2
!    kx = 2.0d0*pi*dfloat(i-1)/Lx
!    ky = 2.0d0*pi*dfloat(j-1)/Ly
!    kz = 2.0d0*pi*float(k-1)/Lz
!      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
!      Ekfull(i,j,k) = 0.0d0
!      Bkfull(i,j,k) = 0.0d0
!      endif
    !write(35,*) i,j,k,sqrt(float((i-1)**2)+float((j-1)**2)+float((k-1)**2)),&
    !          abs(Ekfull(i,j,k))/(time_max/(dt*10.0d0)), abs(Bkfull(i,j,k))/(time_max/(dt*10.0d0))
!    enddo
!    do k = Nz/2+1,Nz
!    kx = 2.0d0*pi*float(i-1)/Lx
!    ky = 2.0d0*pi*float(j-1)/Ly
!    kz = 2.0d0*pi*float((k-1)-Nz)/Lz
!      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
!      Ekfull(i,j,k) = 0.0d0
!      Bkfull(i,j,k) = 0.0d0
!      endif
    !write(35,*) i,j,k-Nz,sqrt(float((i-1)**2)+float((j-1)**2)+float(((k-1)-Nz)**2)),&
    !          abs(Ekfull(i,j,k))/(time_max/(dt*10.0d0)), abs(Bkfull(i,j,k))/(time_max/(dt*10.0d0))
!    enddo
!  enddo  
!  do j = Ny/2+1,Ny
!    do k = 1,Nz/2
!    kx = 2.0d0*pi*dfloat(i-1)/Lx
!    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
!    kz = 2.0d0*pi*float(k-1)/Lz
!      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
!      Ekfull(i,j,k) = 0.0d0
!      Bkfull(i,j,k) = 0.0d0
!      endif
    !write(35,*) i,j-Ny,k,sqrt(float((i-1)**2)+float(((j-1)-Ny)**2)+float((k-1)**2)),&
    !          abs(Ekfull(i,j,k))/(time_max/(dt*10.0d0)), abs(Bkfull(i,j,k))/(time_max/(dt*10.0d0))
!    enddo
!    do k = Nz/2+1,Nz
!    kx = 2.0d0*pi*float(i-1)/Lx
!    ky = 2.0d0*pi*float((j-1)-Ny)/Ly
!    kz = 2.0d0*pi*float((k-1)-Nz)/Lz
!      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
!      Ekfull(i,j,k) = 0.0d0
!      Bkfull(i,j,k) = 0.0d0
!      endif
    !write(35,*) i,j-Ny,k-Nz,sqrt(float((i-1)**2)+float(((j-1)-Ny)**2)+float(((k-1)-Nz)**2)),&
    !          abs(Ekfull(i,j,k))/(time_max/(dt*10.0d0)), abs(Bkfull(i,j,k))/(time_max/(dt*10.0d0))
!    enddo
!  enddo  
!enddo 
!endif 
if (myrank.eq.0) then
 close(5)
 close(35)
 close(40)
endif
! call destroyCUFFTPlan3D(fftPlanD2ZMain)
! call destroyCUFFTPlan3D(fftPlanZ2DMain)
call accfftDestroyPlan3DR2CGPU(fftPlan)
call accfftDestroyComm(comm)

!Perform deallocations
if (myrank.eq.0 ) then
deallocate(Ekfull)
deallocate(Bkfull)
endif
deallocate(uxfull,uyfull,uzfull)
deallocate(Bxfull,Byfull,Bzfull)
deallocate(x,y,z)
deallocate(rho,ux,uy,uz)
deallocate(rho_k,ux_k,uy_k,uz_k)
deallocate(rho_ux,rho_uy,rho_uz)
deallocate(rho_ux_k,rho_uy_k,rho_uz_k)
deallocate(omega_x,omega_y,omega_z,omega2,P,E,Bx,By,Bz,B2,div_B)
deallocate(omega_x_k,omega_y_k,omega_z_k,P_k,E_k,Bx_k,By_k,Bz_k,div_B_k)
!deallocate(rho_k,ux_k,uy_k,uz_k,rho_ux_k,rho_uy_k,rho_uz_k,omega_x_k,omega_y_k,omega_z_k)
!deallocate(P_k,E_k,Ek,Bk,Bx_k,By_k,Bz_k,div_B_k)
deallocate(Ek)
deallocate(Bk)
deallocate(rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new)
deallocate(E_k_new,Bx_k_new,By_k_new,Bz_k_new)
deallocate(d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old)
deallocate(d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old)
deallocate(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new)
deallocate(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new)
!Variables used in derive
deallocate(Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3)
deallocate(Mom_x_1_k,Mom_x_2_k,Mom_x_3_k,Mom_y_1_k,Mom_y_2_k,Mom_y_3_k,Mom_z_1_k,Mom_z_2_k,Mom_z_3_k)

deallocate(d_ux_dx,d_uy_dy,d_uz_dz,Fx,Fy,Fz,Energy_x,Energy_y,Energy_z,E_Visc)
deallocate(d_ux_dx_k,d_uy_dy_k,d_uz_dz_k,Fx_k,Fy_k,Fz_k,Energy_x_k,Energy_y_k,Energy_z_k,E_Visc_k)

deallocate(Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2,d_Bx_dy,d_By_dx,d_Bx_dz,d_By_dz,d_Bz_dx,d_Bz_dy)
deallocate(Mag_x_1_k,Mag_x_2_k,Mag_y_1_k,Mag_y_2_k,Mag_z_1_k,Mag_z_2_k,d_Bx_dy_k,d_By_dx_k,d_Bx_dz_k,d_By_dz_k,d_Bz_dx_k,d_Bz_dy_k)

!deallocate(curl_x_B,curl_y_B,curl_z_B)

!End variables used in derive

call MPI_FINALIZE(mpierr)
end program PSMHD3  

