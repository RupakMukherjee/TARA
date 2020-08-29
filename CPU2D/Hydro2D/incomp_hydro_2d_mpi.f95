Program rupak
!
! This is MPI Benchmarked Two Dimensional Incompressible Viscous Fluid code, 
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
use, intrinsic :: iso_c_binding
implicit none
include 'mpif.h'

! Define Grid Size.
integer (C_INTPTR_T), parameter :: Nx = 128
integer (C_INTPTR_T), parameter :: Ny = 128
integer (C_INTPTR_T), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi = 3.14159265358979323846d0

include "fftw3-mpi.f03"

integer (C_INTPTR_T) :: i,j,t
REAL (C_DOUBLE) ::  Lx,Ly,dx,dy,kx,ky,time,time_min,time_max,dt,Rex,Rey,nux,nuy,W0,m,d,t1,t2,tmp
real (C_DOUBLE) ::  Energy,y_Energy
real (C_DOUBLE) ::  x(Nx),y(Ny)
real (C_DOUBLE), dimension(Nx,Ny) :: psi,omega,omega_dum,ux,uy,ux_dum,uy_dum
real (C_DOUBLE), dimension(Nx,Ny) :: domega_dx,domega_dy
real (C_DOUBLE), dimension(Nx,Ny) :: ux_domega_dx,uy_domega_dy
complex (C_DOUBLE_COMPLEX), dimension(Nh,Ny) :: omegak,omegak_dum,omegak_new,dt_omegak_old,dt_omegak_new
complex (C_DOUBLE_COMPLEX), dimension(Nh,Ny) :: psik,psik_dum,ukx,uky,ukx_dum,uky_dum,Ek
complex (C_DOUBLE_COMPLEX), dimension(Nh,Ny) :: i_kx_omegak,i_ky_omegak
complex (C_DOUBLE_COMPLEX), dimension(Nh,Ny) :: NLkx,NLky,NLk

integer :: ierr, myid, nproc, seed, root
type(C_PTR) :: planf, planb, cdatar, cdatac ! FFTW Variables.
integer(C_INTPTR_T) :: alloc_local, local_Ny, local_j_offset
real(C_DOUBLE), pointer :: idata(:,:)
complex(C_DOUBLE_complex), pointer :: odata(:,:)

open(unit=5,file='System_information.dat',status='unknown')
open(unit=10,file='Initial_Condition_Omega.dat',status='unknown')
open(unit=15,file='Initial_Condition_Velocity.dat',status='unknown')
open(unit=20,file='Initial_Condition_Velocity_Reproduced.dat',status='unknown')
open(unit=40,file='Energy.dat',status='unknown')

!===================== USER INPUTS ============================================		

  call mpi_init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

root = 0

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
    write(15,*) x(i),y(j),omega(i,j),ux(i,j),uy(i,j)
  end do ! j
end do ! i

 close (15)

!===================== FFTW MPI PLAN CREATION ==========================================

  call fftw_mpi_init()
  
  ! get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_2d(Ny, Nh, MPI_COMM_WORLD, local_Ny, local_j_offset)
  cdatar = fftw_alloc_real(2 * alloc_local)
  cdatac = fftw_alloc_complex(alloc_local)
  
  call c_f_pointer(cdatar, idata, [2*Nh,local_Ny])
  call c_f_pointer(cdatac, odata, [Nh,local_Ny])

  ! create MPI plan for out-of-place DFT (note dimension reversal)
  planf = fftw_mpi_plan_dft_r2c_2d(Ny, Nx, idata, odata, MPI_COMM_WORLD, FFTW_MEASURE)
  planb = fftw_mpi_plan_dft_c2r_2d(Ny, Nx, odata, idata, MPI_COMM_WORLD, FFTW_MEASURE)

!===================== INITIAL TIME DATA ===============================================

! Move to Spectral Space.
  do i = 1, Nx
    do j = 1, local_Ny
      idata(i, j) = omega(i, j + local_j_offset)
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, local_Ny
      omegak(i,j) = odata(i,j)
    end do
  end do

!======================= MAIN PROGRAM =====================================================

t1 = MPI_Wtime() 

do time = time_min,time_max,dt

t = nint(time/dt) - int(time_min/dt)

! Evaluation of Initial Stream-Function and Velocity Profile.
do i = 1,Nx/2+1
  do j = 1,local_Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j+local_j_offset-1)/Ly
      if (j+local_j_offset-1 <= Ny/2) then
        ky = 2.0d0*pi*dfloat(j+local_j_offset-1)/Ly
          if (i == 1 .and. j == 1) then
            psik(i,j) = (0.0d0,0.0d0)
            ukx(i,j) = (0.0d0,0.0d0)
            uky(i,j) = (0.0d0,0.0d0)
          else
            psik(i,j) = omegak(i,j)/( kx*kx + ky*ky ) 
            ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j) 
            uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j) 
          endif
      else
        ky = 2.0d0*pi*dfloat(j+local_j_offset-1-Ny)/Ly 
        psik(i,j) = omegak(i,j)/( kx*kx + ky*ky ) 
        ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j) 
        uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j) 
      endif
  end do ! j
end do ! i

! Evaluation of Derivatives in Spectral Space.
do i = 1,Nx/2+1
  do j = 1,local_Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j+local_j_offset-1)/Ly
      if (j+local_j_offset-1 <= Ny/2) then
        ky = 2.0d0*pi*dfloat(j+local_j_offset-1)/Ly
        i_kx_omegak(i,j) = (0.0d0,1.0d0)*kx*omegak(i,j)
        i_ky_omegak(i,j) = (0.0d0,1.0d0)*ky*omegak(i,j)
      else
        ky = 2.0d0*pi*dfloat(j+local_j_offset-1-Ny)/Ly 
        i_kx_omegak(i,j) = (0.0d0,1.0d0)*kx*omegak(i,j)
        i_ky_omegak(i,j) = (0.0d0,1.0d0)*ky*omegak(i,j)      
      endif
  end do ! j
end do ! i

! Move to Real Space.
  do i = 1, Nh
    do j = 1, local_Ny
      odata(i, j) = ukx(i, j)
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, local_Ny
      ux(i, j) = idata(i, j)/dfloat(Nx*Ny)
    end do
  end do

  do i = 1, Nh
    do j = 1, local_Ny
      odata(i, j) = uky(i, j)
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, local_Ny
      uy(i, j) = idata(i, j)/dfloat(Nx*Ny)
    end do
  end do

  do i = 1, Nh
    do j = 1, local_Ny
      odata(i, j) = i_kx_omegak(i, j)
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, local_Ny
      domega_dx(i, j) = idata(i, j)/dfloat(Nx*Ny)
    end do
  end do

  do i = 1, Nh
    do j = 1, local_Ny
      odata(i, j) = i_ky_omegak(i, j)
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, local_Ny
      domega_dy(i, j) = idata(i, j)/dfloat(Nx*Ny)
    end do
  end do

! Evaluation of Non-Linear Terms.  
do i = 1,Nx
  do j = 1,local_Ny
    ux_domega_dx(i,j) = ux(i,j)*domega_dx(i,j)
    uy_domega_dy(i,j) = uy(i,j)*domega_dy(i,j)
  end do ! j
end do ! i

! Move to Spectral Space.
  do i = 1, Nx
    do j = 1, local_Ny
      idata(i, j) = ux_domega_dx(i, j + local_j_offset)
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, local_Ny
      NLkx(i,j) = odata(i,j)
    end do
  end do

  do i = 1, Nx
    do j = 1, local_Ny
      idata(i, j) = uy_domega_dy(i, j + local_j_offset)
    end do
  end do

  call fftw_mpi_execute_dft_r2c(planf, idata, odata)

  do i = 1, Nh
    do j = 1, local_Ny
      NLky(i,j) = odata(i,j)
    end do
  end do

do i = 1,Nx/2+1
  do j = 1,local_Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (j+local_j_offset-1 <= Ny/2) then
        ky = 2.0d0*pi*float(j+local_j_offset-1)/Ly
        ! De - Aliazing Technique - 2/3 Truncation...
          if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0) then    
            NLkx(i,j) = 0.0d0
            NLky(i,j) = 0.0d0
          endif 
      else 
        ky = 2.0d0*pi*float(j+local_j_offset-1-Ny)/Ly
        ! De - Aliazing Technique - 2/3 Truncation...
          if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0) then
            NLkx(i,j) = 0.0d0
            NLky(i,j) = 0.0d0
          endif 
      endif    
  end do ! j
end do ! i 

! Evaluation of Nonlinear Term.  
do i = 1,Nx/2+1
  do j = 1,local_Ny  
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (j+local_j_offset-1 <= Ny/2) then
        ky = 2.0d0*pi*dfloat(j+local_j_offset-1)/Ly
        NLk(i,j) = ( NLkx(i,j) + NLky(i,j) )
        Nlk(i,j) = NLk(i,j) + (nux*kx*kx + nuy*ky*ky)*omegak(i,j)
      else 
        ky = 2.0d0*pi*dfloat(j+local_j_offset-1-Ny)/Ly
        NLk(i,j) = ( NLkx(i,j) + NLky(i,j) )
        Nlk(i,j) = NLk(i,j) + (nux*kx*kx + nuy*ky*ky)*omegak(i,j)
      endif  
  end do ! j
end do ! i

! Preparing for Time Evolution.
do i = 1,Nx/2+1
  do j = 1,local_Ny  
    dt_omegak_new(i,j) = - NLk(i,j)
  end do ! j
end do ! i 

! Adams-Bashforth Method for Time Evolution.
do i = 1,Nx/2+1
  do j = 1,local_Ny
    omegak_new(i,j) = omegak(i,j) + ( (3.0d0/2.0d0)*dt_omegak_new(i,j) - (1.0d0/2.0d0)*dt_omegak_old(i,j) )*dt
  end do ! j
end do ! i

! Reset the Values.
do i = 1,Nx/2+1
  do j = 1,local_Ny
    dt_omegak_old(i,j) = dt_omegak_new(i,j)
    omegak(i,j) = omegak_new(i,j)
  end do ! j
end do ! i

! Keep backup for FFTW.
do i = 1,Nx/2+1
  do j = 1,local_Ny
    omegak_dum(i,j) = omegak(i,j) 
  end do ! j
end do ! i

! Move to Real Space.
  do i = 1, Nh
    do j = 1, local_Ny
      odata(i, j) = omegak_dum(i, j)
    end do
  end do

  call fftw_mpi_execute_dft_c2r(planb, odata, idata)

  do i = 1, Nx
    do j = 1, local_Ny
      omega(i, j) = idata(i, j)/dfloat(Nx*Ny)
    end do
  end do

! FFTW Normalisation and Data Printing.   
Energy = 0.0d0

do i = 1,Nx
  do j = 1,local_Ny
    Energy = Energy + (omega(i,j)**2)/(2.0d0*dfloat(Nx*Ny))
  end do ! j
end do ! i

  CALL MPI_REDUCE(Energy, tmp, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  Energy = tmp

if (myid == root) then  
  if (mod(float(t),100.0) == 0.0) then
  write(40,*) time,Energy
  call flush(40) 
  endif
endif
  
end do ! time

t2 = MPI_Wtime()

write(5,*) "Elapsed time is", t2 - t1

  ! deallocate and destroy plans  
  call fftw_destroy_plan(planf)
  call fftw_destroy_plan(planb)
  call fftw_mpi_cleanup()
  call fftw_free(cdatar)
  call fftw_free(cdatac)

!====================================================================================

end program rupak  


