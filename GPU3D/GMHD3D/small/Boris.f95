Program Tracer
!
! Rupak, IPR, August 29, 2018
!
use openacc
implicit none

! Define Grid Size.
integer ( kind = 4 ), parameter :: Np = 10000000
integer ( kind = 4 ), parameter :: Nx = 128
integer ( kind = 4 ), parameter :: Ny = 256
integer ( kind = 4 ), parameter :: Nz = 512
integer ( kind = 4 ) i,j,k,t
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0
real ( kind = 8 ) Lx,Ly,Lz,dx,dy,dz,time,time_min,time_max,dt,h,U0,B0,A,B,C,kf
real ( kind = 8 ) x(Nx),y(Ny),z(Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
real ( kind = 8 ) Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz)
real ( kind = 8 ) xp(Np),yp(Np),zp(Np),uxp(Np),uyp(Np),uzp(Np)
real ( kind = 8 ) Bpx(Np), Bpy(Np), Bpz(Np),vxp(Np),vyp(Np),vzp(Np)
common/comm/Lx,Ly,Lz,dx,dy,dz,dt,h
integer,parameter :: seed = 99999999
double precision rand ! Added by Naga to avoid compilation error
real (kind = 8) starttime,stoptime
call srand(seed)
call acc_init(acc_device_nvidia)
!===================== FILENAMES ==============================================	

open(unit=15,file='Initial_Grid_Data.dat',status='unknown')
open(unit=10,file='data_multicore.dat',status='unknown')


!===================== USER INPUTS ============================================		

! System Size.
Lx = 2.0d0*pi; Ly = 4.0d0*pi; Lz = 8.0d0*pi 

! Grid Resolution.
dx = Lx/dfloat(Nx); dy = Ly/dfloat(Ny); dz = Lz/dfloat(Nz)

! Runtime Details and Time Step Width.
time_min = 0.0d0
time_max = 1.0d0
!time_max = 0.009d0
dt = 0.10d0
h = 1.0d0/(dx*dy*dz)

A = 1.0d0; B = 1.0d0; C = 1.0d0; kf = 4.0d0

! Maximum Velocity.
U0 = 0.10d0

! Initial Magnetic Field.
B0 = 0.10d0

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
  x(i) = dfloat(i-1)*dx
  do j = 1, Ny
    y(j) = dfloat(j-1)*dy
    do k = 1, Nz
      z(k) = dfloat(k-1)*dz
      ! Initial Velocity Distribution.
      ux(i,j,k) = U0*(A*dsin(kf*z(k)) + C*dcos(kf*y(j)))
      uy(i,j,k) = U0*(B*dsin(kf*x(i)) + A*dcos(kf*z(k)))
      uz(i,j,k) = U0*(C*dsin(kf*y(j)) + B*dcos(kf*x(i)))
      ! Initial Magnetic Field Profile.
      Bx(i,j,k) = B0*(A*dsin(kf*z(k)) + C*dcos(kf*y(j)))
      By(i,j,k) = B0*(B*dsin(kf*x(i)) + A*dcos(kf*z(k)))
      Bz(i,j,k) = B0*(C*dsin(kf*y(j)) + B*dcos(kf*x(i)))
      !write(15,*) x(i),y(j),z(k),ux(i,j,k),uy(i,j,k),uz(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    enddo
  end do
end do

 close(15)

!======================= MAIN PROGRAM =====================================================
! Time Loop Starts.
call cpu_time(starttime)
!$acc data copyin(ux,uy,uz,Bx,By,Bz,x,y,z)  copy(xp,yp,zp) &
!$acc& create(Bpx,Bpy,Bpz,uxp,uyp,uzp) copyout(vxp,vyp,vzp)

do time = time_min,time_max,dt

t = nint(time/dt) - dint(time_min/dt)

! Online Check of the Progress.
!if (mod(t,int(1.0/dt)) == 0) then
  !print*, t/int(1.0/dt)
!endif

  call tracer_particle (Np,Nx,Ny,Nz,pi,x,y,z,xp,yp,zp,ux,uy,uz,uxp,uyp,uzp,Bx,By,Bz,Bpx,Bpy,Bpz,vxp,vyp,vzp)

 ! if (time==time_max-dt) then
  
 
 ! !close((i-1)+100) 
!  endif
  
enddo ! time
!$acc end data
call cpu_time(stoptime)
write(*,*) 'Boris time in sec: ',stoptime-starttime
do i=1,NP
write(10,*) xp(i),yp(i),zp(i)
end do

contains

!====================================================================================

! This subroutine is implemented by Vinod Saini and Rupak Mukherjee, IPR on 7th April, 2018.

subroutine tracer_particle (Np,Nx,Ny,Nz,pi,xg,yg,zg,x,y,z,ux,uy,uz,uxp,uyp,uzp,Bx,By,Bz,Bpx,Bpy,Bpz,vxp,vyp,vzp)
implicit none
integer ( kind = 4 ) Np,Nx,Ny,Nz
integer ( kind = 4 ) i,j,k,m,n,p
real ( kind = 8 ) pi
double precision:: Lx,Ly,Lz,dx,dy,dz,dt,h,sx,sy,sz,tx,ty,tz
double precision:: x,y,z,uxp,uyp,uzp,Bx,By,Bz,Bpx,Bpy,Bpz
double precision:: xg,yg,zg,vxp,vyp,vzp,vxd,vyd,vzd,ux,uy,uz ! g stand for grid quantities
dimension:: x(Np),y(Np),z(Np),uxp(Np),uyp(Np),uzp(Np),Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),Bpx(Np),Bpy(Np),Bpz(Np)
dimension :: xg(Nx),yg(Ny),zg(Nz),vxp(Np),vyp(Np),vzp(Np)
dimension :: ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
common/comm/Lx,Ly,Lz,dx,dy,dz,dt,h


! Velocity and Magnetic Field Interpolation from Grid to Tracer Particle...
!$acc parallel loop vector_length(256) & 
!$acc& present(ux,uy,uz,Bx,By,Bz,x,y,z,uxp,uyp,uzp,Bpx,Bpy,Bpz,xg,yg,zg)
do i = 1,Np

  m = int(x(i)/dx) + 1

  n = int(y(i)/dy) + 1
  
  p = int(z(i)/dz) + 1

if (m==Nx .and. n/=Ny .and. p/=Nz) then

  uxp(i) = ux(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           ux(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           ux(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           ux(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           ux(m,n,p+1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           ux(1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           ux(m,n+1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           ux(1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpx(i) = Bx(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           Bx(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           Bx(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bx(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bx(m,n,p+1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bx(1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m,n+1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bx(1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  uyp(i) = uy(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           uy(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           uy(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uy(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uy(m,n,p+1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uy(1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uy(m,n+1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uy(1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  
  Bpy(i) = By(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           By(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           By(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           By(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           By(m,n,p+1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           By(1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           By(m,n+1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           By(1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  uzp(i) = uz(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           uz(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           uz(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uz(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uz(m,n,p+1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uz(1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uz(m,n+1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uz(1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpz(i) = Bz(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           Bz(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           Bz(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bz(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bz(m,n,p+1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bz(1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m,n+1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bz(1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

end if

if (m/=Nx .and. n==Ny .and. p/=Nz) then

  uxp(i) = ux(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           ux(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           ux(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           ux(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           ux(m,n,p+1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           ux(m+1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           ux(m,1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           ux(m+1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpx(i) = Bx(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           Bx(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           Bx(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bx(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bx(m,n,p+1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m+1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m,1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bx(m+1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  uyp(i) = uy(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           uy(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           uy(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uy(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uy(m,n,p+1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uy(m+1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uy(m,1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uy(m+1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  
  Bpy(i) = By(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           By(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           By(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           By(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           By(m,n,p+1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           By(m+1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           By(m,1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           By(m+1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  uzp(i) = uz(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           uz(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           uz(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uz(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uz(m,n,p+1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uz(m+1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uz(m,1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uz(m+1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpz(i) = Bz(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           Bz(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           Bz(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bz(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bz(m,n,p+1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m+1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m,1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bz(m+1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

endif

if (m/=Nx .and. n/=Ny .and. p==Nz) then

  uxp(i) = ux(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           ux(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           ux(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           ux(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           ux(m,n,1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           ux(m+1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           ux(m,n+1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           ux(m+1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpx(i) = Bx(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           Bx(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           Bx(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bx(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bx(m,n,1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m+1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m,n+1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bx(m+1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  uyp(i) = uy(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           uy(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           uy(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uy(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uy(m,n,1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uy(m+1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uy(m,n+1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uy(m+1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpy(i) = By(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           By(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           By(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           By(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           By(m,n,1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           By(m+1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           By(m,n+1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           By(m+1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h
  
  uzp(i) = uz(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           uz(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           uz(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uz(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uz(m,n,1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uz(m+1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uz(m,n+1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uz(m+1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpz(i) = Bz(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           Bz(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           Bz(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bz(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bz(m,n,1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m+1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m,n+1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bz(m+1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

endif

if (m==Nx .and. n==Ny .and. p/=Nz ) then

  uxp(i) = ux(m,n,p) * (Lx-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           ux(1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           ux(m,1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           ux(1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           ux(m,n,p+1) * (Lx-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           ux(1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           ux(m,1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           ux(1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpx(i) = Bx(m,n,p) * (Lx-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           Bx(1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           Bx(m,1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bx(1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bx(m,n,p+1) * (Lx-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bx(1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m,1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bx(1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  uyp(i) = uy(m,n,p) * (Lx-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           uy(1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           uy(m,1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uy(1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uy(m,n,p+1) * (Lx-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uy(1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uy(m,1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uy(1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpy(i) = By(m,n,p) * (Lx-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           By(1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           By(m,1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           By(1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           By(m,n,p+1) * (Lx-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           By(1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           By(m,1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           By(1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h
  
  uzp(i) = uz(m,n,p) * (Lx-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           uz(1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           uz(m,1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uz(1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uz(m,n,p+1) * (Lx-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uz(1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uz(m,1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uz(1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpz(i) = Bz(m,n,p) * (Lx-x(i)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           Bz(1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (zg(p+1)-z(i)) * h + &
           Bz(m,1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bz(1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bz(m,n,p+1) * (Lx-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bz(1,n,p+1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m,1,p+1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bz(1,1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

endif

if (m==Nx .and. n/=Ny .and. p==Nz ) then

  uxp(i) = ux(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           ux(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           ux(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           ux(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           ux(m,n,1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           ux(1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           ux(m,n+1,1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           ux(1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpx(i) = Bx(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           Bx(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           Bx(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bx(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bx(m,n,1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bx(1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m,n+1,1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bx(1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  uyp(i) = uy(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           uy(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           uy(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uy(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uy(m,n,1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uy(1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uy(m,n+1,1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uy(1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpy(i) = By(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           By(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           By(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           By(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           By(m,n,1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           By(1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           By(m,n+1,1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           By(1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h
  
  uzp(i) = uz(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           uz(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           uz(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uz(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uz(m,n,1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uz(1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uz(m,n+1,1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uz(1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpz(i) = Bz(m,n,p) * (Lx-x(i)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           Bz(1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (Lz-z(i)) * h + &
           Bz(m,n+1,p) * (Lx-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bz(1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bz(m,n,1) * (Lx-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bz(1,n,1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m,n+1,1) * (Lx-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bz(1,n+1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

endif

if (m/=Nx .and. n==Ny .and. p==Nz ) then

  uxp(i) = ux(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           ux(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           ux(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           ux(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           ux(m,n,1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           ux(m+1,n,1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           ux(m,1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           ux(m+1,1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpx(i) = Bx(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           Bx(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           Bx(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bx(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bx(m,n,1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m+1,n,1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m,1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bx(m+1,1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  uyp(i) = uy(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           uy(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           uy(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uy(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uy(m,n,1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uy(m+1,n,1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uy(m,1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uy(m+1,1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpy(i) = By(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           By(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           By(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           By(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           By(m,n,1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           By(m+1,n,1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           By(m,1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           By(m+1,1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h
  
  uzp(i) = uz(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           uz(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           uz(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uz(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           uz(m,n,1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uz(m+1,n,1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           uz(m,1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uz(m+1,1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpz(i) = Bz(m,n,p) * (xg(m+1)-x(i)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           Bz(m+1,n,p) * (x(i)-xg(m)) * (Ly-y(i)) * (Lz-z(i)) * h + &
           Bz(m,1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bz(m+1,1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (Lz-z(i)) * h + &
           Bz(m,n,1) * (xg(m+1)-x(i)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m+1,n,1) * (x(i)-xg(m)) * (Ly-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m,1,1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bz(m+1,1,1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

endif

if (m/=Nx .and. n/=Ny .and. p/=Nz) then

  uxp(i) = ux(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           ux(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           ux(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           ux(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           ux(m,n,p+1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           ux(m+1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           ux(m,n+1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           ux(m+1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpx(i) = Bx(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           Bx(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           Bx(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bx(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bx(m,n,p+1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m+1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bx(m,n+1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bx(m+1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  uyp(i) = uy(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           uy(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           uy(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uy(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uy(m,n,p+1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uy(m+1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uy(m,n+1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uy(m+1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpy(i) = By(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           By(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           By(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           By(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           By(m,n,p+1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           By(m+1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           By(m,n+1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           By(m+1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h
  
  uzp(i) = uz(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           uz(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           uz(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uz(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           uz(m,n,p+1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uz(m+1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           uz(m,n+1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           uz(m+1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h

  Bpz(i) = Bz(m,n,p) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           Bz(m+1,n,p) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (zg(p+1)-z(i)) * h + &
           Bz(m,n+1,p) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bz(m+1,n+1,p) * (x(i)-xg(m)) * (y(i)-yg(n)) * (zg(p+1)-z(i)) * h + &
           Bz(m,n,p+1) * (xg(m+1)-x(i)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m+1,n,p+1) * (x(i)-xg(m)) * (yg(n+1)-y(i)) * (z(i)-zg(p)) * h + &
           Bz(m,n+1,p+1) * (xg(m+1)-x(i)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h + &
           Bz(m+1,n+1,p+1) * (x(i)-xg(m)) * (y(i)-yg(n)) * (z(i)-zg(p)) * h
end if

end do ! i
!$acc end parallel

! Boris Algorithm...
!$acc parallel loop vector_length(256) present(vxp,vyp,vzp,uxp,uyp,uzp,x,y,z,Bpx,Bpy,Bpz) &
!$acc& private(vxd,vyd,vzd)
 
 do i=1,Np
  tx = 0.5d0*dt*Bpx(i)

  ty = 0.5d0*dt*Bpy(i)
  tz = 0.5d0*dt*Bpz(i)
  
  sx = 2.0d0*tx/(1.0d0+tx**2+ty**2+tz**2)
  sy = 2.0d0*ty/(1.0d0+tx**2+ty**2+tz**2)
  sz = 2.0d0*tz/(1.0d0+tx**2+ty**2+tz**2)

  vxd = uxp(i) + (uyp(i)*tz-uzp(i)*ty)
  vyd = uyp(i) - (uxp(i)*tz-uzp(i)*tx)
  vzd = uzp(i) + (uxp(i)*ty-uyp(i)*tx)
  
  vxp(i) = uxp(i) + (vyd*sz-vzd*sy)
  vyp(i) = uyp(i) - (vxd*sz-vzd*sx)
  vzp(i) = uzp(i) + (vxd*sy-vyd*sx)
  
  x(i) = x(i) + vxp(i)*dt
  y(i) = y(i) + vyp(i)*dt
  z(i) = z(i) + vzp(i)*dt
  
  ! Periodic Boundary Condition Implemented...
  x(i) = x(i) - (int(x(i)/Lx))*Lx
    if (x(i) .lt. 0.0d0) then
    x(i) = x(i) + Lx
    endif
  
  y(i) = y(i) - (int(y(i)/Ly))*Ly
    if (y(i) .lt. 0.0d0) then
    y(i) = y(i) + Ly
    endif

  z(i) = z(i) - (int(z(i)/Lz))*Lz
    if (z(i) .lt. 0.0d0) then
    z(i) = z(i) + Lz
    endif

  enddo
!$acc end parallel
end subroutine tracer_particle
  
!====================================================================================

end program Tracer

