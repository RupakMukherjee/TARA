! Read Input Parameters.

!===================================================================================
!================== SUBROUTINE Input - Parameter ===================================
!===================================================================================

subroutine Input_Parameter(Nx,Ny,Nz,Nh,pi,thread_num,Lx,Ly,Lz,dx,dy,dz,time_min,time_max,dt,&
                           spheat,mu,eta,mu_0,ms,rho0,P0,M,MA,kf,CS,U0,VA,B0)

use omp_lib
implicit none

integer ( kind = 4 ) G,Nx,Ny,Nz,Nh
real ( kind = 8 ) pi,threads,Lx,Ly,Lz,dx,dy,dz,time_min,time_max,dt
real ( kind = 8 ) spheat,mu,eta,mu_0,ms,rho0,P0,M,MA,kf,CS,U0,VA,B0,Re,Rm,PM

integer ( kind = 4 ) thread_num,num_thread,proc_num ! OMP

!===================== FILENAMES ==============================================	

open(unit=1,file='INPUT.DATA',status='old')
open(unit=5,file='System_information.dat',status='unknown')

!==============================================================================	

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

return

end subroutine Input_Parameter

!====================================================================================

