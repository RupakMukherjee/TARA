!gcc dcst_3d.c -c && gfortran fec.f95 dcst_3d.o -lfftw3 && ./a.out
implicit none
! Define Grid Size.
integer ( kind = 4 ), parameter :: Nx = 64
integer ( kind = 4 ), parameter :: Ny = 64
integer ( kind = 4 ), parameter :: Nz = 64
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

integer (kind=4) i,j,k,t
real (kind=8) Lx,Ly,Lz,dx,dy,dz,time,time_min,time_max,dt,A,B,C,kf,Re,Rm
real (kind=8) MS,MA,U0,CS,VA,B0,Energy,B_Field,C1,C2,C3,Hf,Hm,HuB,HAw,R1,R2,Rayleigh,dBT
real (kind=8) x(Nx),y(Ny),z(Nz)
real (kind=8) ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz),P(Nx,Ny,Nz)
real (kind=8) wx(Nx,Ny,Nz),wy(Nx,Ny,Nz),wz(Nx,Ny,Nz),jx(Nx,Ny,Nz),jy(Nx,Ny,Nz),jz(Nx,Ny,Nz)
real (kind=8) Ax(Nx,Ny,Nz),Ay(Nx,Ny,Nz),Az(Nx,Ny,Nz)
real (kind=8) u2(Nx,Ny,Nz),w2(Nx,Ny,Nz),j2(Nx,Ny,Nz),dB(Nx,Ny,Nz)
real (kind=8) dux_dx(Nx,Ny,Nz),dux_dy(Nx,Ny,Nz),dux_dz(Nx,Ny,Nz)
real (kind=8) duy_dx(Nx,Ny,Nz),duy_dy(Nx,Ny,Nz),duy_dz(Nx,Ny,Nz)
real (kind=8) duz_dx(Nx,Ny,Nz),duz_dy(Nx,Ny,Nz),duz_dz(Nx,Ny,Nz)
real (kind=8) d2_ux_dx2(Nx,Ny,Nz),d2_ux_dy2(Nx,Ny,Nz),d2_ux_dz2(Nx,Ny,Nz)
real (kind=8) d2_uy_dx2(Nx,Ny,Nz),d2_uy_dy2(Nx,Ny,Nz),d2_uy_dz2(Nx,Ny,Nz)
real (kind=8) d2_uz_dx2(Nx,Ny,Nz),d2_uz_dy2(Nx,Ny,Nz),d2_uz_dz2(Nx,Ny,Nz)
real (kind=8) Mom_x_1(Nx,Ny,Nz),Mom_x_2(Nx,Ny,Nz),Mom_x_3(Nx,Ny,Nz)
real (kind=8) Mom_y_1(Nx,Ny,Nz),Mom_y_2(Nx,Ny,Nz),Mom_y_3(Nx,Ny,Nz)
real (kind=8) Mom_z_1(Nx,Ny,Nz),Mom_z_2(Nx,Ny,Nz),Mom_z_3(Nx,Ny,Nz)
real (kind=8) dMom_x_1_dx(Nx,Ny,Nz),dMom_x_2_dy(Nx,Ny,Nz),dMom_x_3_dz(Nx,Ny,Nz)
real (kind=8) dMom_y_1_dx(Nx,Ny,Nz),dMom_y_2_dy(Nx,Ny,Nz),dMom_y_3_dz(Nx,Ny,Nz)
real (kind=8) dMom_z_1_dx(Nx,Ny,Nz),dMom_z_2_dy(Nx,Ny,Nz),dMom_z_3_dz(Nx,Ny,Nz)
real (kind=8) Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz)
real (kind=8) dBx_dx(Nx,Ny,Nz),dBx_dy(Nx,Ny,Nz),dBx_dz(Nx,Ny,Nz)
real (kind=8) dBy_dx(Nx,Ny,Nz),dBy_dy(Nx,Ny,Nz),dBy_dz(Nx,Ny,Nz)
real (kind=8) dBz_dx(Nx,Ny,Nz),dBz_dy(Nx,Ny,Nz),dBz_dz(Nx,Ny,Nz)
real (kind=8) d2_Bx_dx2(Nx,Ny,Nz),d2_Bx_dy2(Nx,Ny,Nz),d2_Bx_dz2(Nx,Ny,Nz)
real (kind=8) d2_By_dx2(Nx,Ny,Nz),d2_By_dy2(Nx,Ny,Nz),d2_By_dz2(Nx,Ny,Nz)
real (kind=8) d2_Bz_dx2(Nx,Ny,Nz),d2_Bz_dy2(Nx,Ny,Nz),d2_Bz_dz2(Nx,Ny,Nz)
real (kind=8) Mag_x_1(Nx,Ny,Nz),Mag_x_2(Nx,Ny,Nz)
real (kind=8) Mag_y_1(Nx,Ny,Nz),Mag_y_2(Nx,Ny,Nz)
real (kind=8) Mag_z_1(Nx,Ny,Nz),Mag_z_2(Nx,Ny,Nz)
real (kind=8) dMag_x_1_dy(Nx,Ny,Nz),dMag_x_2_dz(Nx,Ny,Nz)
real (kind=8) dMag_y_1_dx(Nx,Ny,Nz),dMag_y_2_dz(Nx,Ny,Nz)
real (kind=8) dMag_z_1_dx(Nx,Ny,Nz),dMag_z_2_dy(Nx,Ny,Nz)
real (kind=8) dux_dt_new(Nx,Ny,Nz),duy_dt_new(Nx,Ny,Nz),duz_dt_new(Nx,Ny,Nz)
real (kind=8) dux_dt_old(Nx,Ny,Nz),duy_dt_old(Nx,Ny,Nz),duz_dt_old(Nx,Ny,Nz)
real (kind=8) dBx_dt_new(Nx,Ny,Nz),dBy_dt_new(Nx,Ny,Nz),dBz_dt_new(Nx,Ny,Nz)
real (kind=8) dBx_dt_old(Nx,Ny,Nz),dBy_dt_old(Nx,Ny,Nz),dBz_dt_old(Nx,Ny,Nz)
integer,parameter :: seed = 99999999
call srand(seed)

open(unit=5,file='IC.dat',status='unknown')
open(unit=10,file='Energy.dat',status='unknown')

! System Size.
Lx = 2.0d0*pi; Ly = 2.0d0*pi; Lz = 2.0d0*pi 

! Grid Resolution.
dx = Lx/dfloat(Nx); dy = Ly/dfloat(Ny); dz = Lz/dfloat(Nz)

! Runtime Details and Time Step Width.
time_min = 0.0d0
time_max = 100.d0
dt = 0.010d0

! Raynold's Numbers.
Re = 1000.0d0
Rm = 1000.0d0

! Alfven Mach Number.
MA = 100.0d0

! Maximum Velocity.
U0 = 0.10d0

! Alfven Speed.
VA = U0/MA

! Initial Magnetic Field.
B0 = VA

! Forcing Length Scale.
kf = 1.0d0

 A = 0.10d0
 B = 0.10d0
 C = 0.10d0
    
do i = 1, Nx
x(i) = dfloat(i-1)*dx
  do j = 1, Ny
  y(j) = dfloat(j-1)*dy
    do k = 1, Nz
    z(k) = dfloat(k-1)*dz
    if (x(i) .ge. Lx/3.0d0 .and. x(i) .le. 2.0d0*Lx/3.0d0 .and. &
        y(j) .ge. Ly/3.0d0 .and. y(j) .le. 2.0d0*Ly/3.0d0 .and. &
        z(k) .ge. Lz/3.0d0 .and. z(k) .le. 2.0d0*Lz/3.0d0) then
      ! Initial Velocity Distribution.
      ux(i,j,k) = U0 * ( A*dsin(kf*z(k)) + C*dcos(kf*y(j)) ) + U0 * 0.01 * (rand()-0.5)
      uy(i,j,k) = U0 * ( B*dsin(kf*x(i)) + A*dcos(kf*z(k)) ) + U0 * 0.01 * (rand()-0.5)
      uz(i,j,k) = U0 * ( C*dsin(kf*y(j)) + B*dcos(kf*x(i)) ) + U0 * 0.01 * (rand()-0.5)
      ! Initial Magnetic Field Distribution.
      Bx(i,j,k) = B0 * ( A*dsin(kf*z(k)) + C*dcos(kf*y(j)) ) + B0 * 0.01 * (rand()-0.5)
      By(i,j,k) = B0 * ( B*dsin(kf*x(i)) + A*dcos(kf*z(k)) ) + B0 * 0.01 * (rand()-0.5)
      Bz(i,j,k) = B0 * ( C*dsin(kf*y(j)) + B*dcos(kf*x(i)) ) + B0 * 0.01 * (rand()-0.5)
    else
      ux(i,j,k) = 0.0d0
      uy(i,j,k) = 0.0d0
      uz(i,j,k) = 0.0d0
      Bx(i,j,k) = 0.0d0
      By(i,j,k) = 0.0d0
      Bz(i,j,k) = 0.0d0
    endif
    ! Initial Pressure Distribution.
    P(i,j,k) = 0.0d0
    write(5,*) x(i),y(j),z(k),ux(i,j,k),uy(i,j,k),uz(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    ! Setting old variables to zero.
    dux_dt_old(i,j,k) = 0.0d0
    duy_dt_old(i,j,k) = 0.0d0
    duz_dt_old(i,j,k) = 0.0d0
    dBx_dt_old(i,j,k) = 0.0d0
    dBy_dt_old(i,j,k) = 0.0d0
    dBz_dt_old(i,j,k) = 0.0d0
    enddo
  enddo
enddo

  close(5)
  
!======================= MAIN PROGRAM =====================================================

! Time Loop Starts.
do time = time_min,time_max,dt

t = nint(time/dt) - dint(time_min/dt)

! Online Check of the Progress.
!if (mod(t,int(1.0/dt)) == 0) then
  !print*, t/int(1.0/dt)
!endif

  call dcsftdx(Nx, Ny, Nz, ux, dux_dx)
  call dcsftdy(Nx, Ny, Nz, uy, duy_dy)
  call dcsftdz(Nx, Ny, Nz, uz, duz_dz)
  
  call dcsftdxx(Nx, Ny, Nz, ux, d2_ux_dx2)
  call dcsftdyy(Nx, Ny, Nz, ux, d2_ux_dy2)
  call dcsftdzz(Nx, Ny, Nz, ux, d2_ux_dz2)

  call dcsftdxx(Nx, Ny, Nz, uy, d2_uy_dx2)
  call dcsftdyy(Nx, Ny, Nz, uy, d2_uy_dy2)
  call dcsftdzz(Nx, Ny, Nz, uy, d2_uy_dz2)

  call dcsftdxx(Nx, Ny, Nz, uz, d2_uz_dx2)
  call dcsftdyy(Nx, Ny, Nz, uz, d2_uz_dy2)
  call dcsftdzz(Nx, Ny, Nz, uz, d2_uz_dz2)

  call dcsftdxx(Nx, Ny, Nz, Bx, d2_Bx_dx2)
  call dcsftdyy(Nx, Ny, Nz, Bx, d2_Bx_dy2)
  call dcsftdzz(Nx, Ny, Nz, Bx, d2_Bx_dz2)

  call dcsftdxx(Nx, Ny, Nz, By, d2_By_dx2)
  call dcsftdyy(Nx, Ny, Nz, By, d2_By_dy2)
  call dcsftdzz(Nx, Ny, Nz, By, d2_By_dz2)

  call dcsftdxx(Nx, Ny, Nz, Bz, d2_Bz_dx2)
  call dcsftdyy(Nx, Ny, Nz, Bz, d2_Bz_dy2)
  call dcsftdzz(Nx, Ny, Nz, Bz, d2_Bz_dz2)

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
    P(i,j,k) = 0.0d0
    ! Evaluate Square of Magnetic Field.
    B2(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)
    ! Evaluate LHS of Momentum Equation.
    Mom_x_1(i,j,k) = ux(i,j,k)*ux(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 - Bx(i,j,k)*Bx(i,j,k)
    Mom_x_2(i,j,k) = ux(i,j,k)*uy(i,j,k) - Bx(i,j,k)*By(i,j,k)
    Mom_x_3(i,j,k) = ux(i,j,k)*uz(i,j,k) - Bx(i,j,k)*Bz(i,j,k)
  
    Mom_y_1(i,j,k) = ux(i,j,k)*uy(i,j,k) - Bx(i,j,k)*By(i,j,k) 
    Mom_y_2(i,j,k) = uy(i,j,k)*uy(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 - By(i,j,k)*By(i,j,k)
    Mom_y_3(i,j,k) = uy(i,j,k)*uz(i,j,k) - By(i,j,k)*Bz(i,j,k)
  
    Mom_z_1(i,j,k) = uz(i,j,k)*ux(i,j,k) - Bz(i,j,k)*Bx(i,j,k) 
    Mom_z_2(i,j,k) = uz(i,j,k)*uy(i,j,k) - Bz(i,j,k)*By(i,j,k) 
    Mom_z_3(i,j,k) = uz(i,j,k)*uz(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 - Bz(i,j,k)*Bz(i,j,k)
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

  call dcsftdx(Nx, Ny, Nz, Mom_x_1, dMom_x_1_dx)
  call dcsftdy(Nx, Ny, Nz, Mom_x_2, dMom_x_2_dy)
  call dcsftdz(Nx, Ny, Nz, Mom_x_3, dMom_x_3_dz)
  
  call dcsftdx(Nx, Ny, Nz, Mom_y_1, dMom_y_1_dx)
  call dcsftdy(Nx, Ny, Nz, Mom_y_2, dMom_y_2_dy)
  call dcsftdz(Nx, Ny, Nz, Mom_y_3, dMom_y_3_dz)

  call dcsftdx(Nx, Ny, Nz, Mom_z_1, dMom_z_1_dx)
  call dcsftdy(Nx, Ny, Nz, Mom_z_2, dMom_z_2_dy)
  call dcsftdz(Nx, Ny, Nz, Mom_z_3, dMom_z_3_dz)

  call dcsftdy(Nx, Ny, Nz, Mag_x_1, dMag_x_1_dy)
  call dcsftdz(Nx, Ny, Nz, Mag_x_2, dMag_x_2_dz)

  call dcsftdx(Nx, Ny, Nz, Mag_y_1, dMag_y_1_dx)
  call dcsftdz(Nx, Ny, Nz, Mag_y_2, dMag_y_2_dz)

  call dcsftdx(Nx, Ny, Nz, Mag_z_1, dMag_z_1_dx)
  call dcsftdy(Nx, Ny, Nz, Mag_z_2, dMag_z_2_dy)
  
do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
    dux_dt_new(i,j,k) = - dMom_x_1_dx(i,j,k) - dMom_x_2_dy(i,j,k) - dMom_x_3_dz(i,j,k) &
                            + ( d2_ux_dx2(i,j,k) + d2_ux_dy2(i,j,k) + d2_ux_dz2(i,j,k) ) / Re
    duy_dt_new(i,j,k) = - dMom_y_1_dx(i,j,k) - dMom_y_2_dy(i,j,k) - dMom_y_3_dz(i,j,k) &
                            + ( d2_uy_dx2(i,j,k) + d2_uy_dy2(i,j,k) + d2_uy_dz2(i,j,k) ) / Re
    duz_dt_new(i,j,k) = - dMom_z_1_dx(i,j,k) - dMom_z_2_dy(i,j,k) - dMom_z_3_dz(i,j,k) &
                            + ( d2_uz_dx2(i,j,k) + d2_uz_dy2(i,j,k) + d2_uz_dz2(i,j,k) ) / Re 

    dBx_dt_new(i,j,k) = + dMag_x_1_dy(i,j,k) + dMag_x_2_dz(i,j,k) &
                        + ( d2_Bx_dx2(i,j,k) + d2_Bx_dy2(i,j,k) + d2_Bx_dz2(i,j,k) ) / Rm
    dBy_dt_new(i,j,k) = - dMag_y_1_dx(i,j,k) + dMag_y_2_dz(i,j,k) &
                        + ( d2_By_dx2(i,j,k) + d2_By_dy2(i,j,k) + d2_By_dz2(i,j,k) ) / Rm
    dBz_dt_new(i,j,k) = - dMag_z_1_dx(i,j,k) - dMag_z_2_dy(i,j,k) &
                        + ( d2_Bz_dx2(i,j,k) + d2_Bz_dy2(i,j,k) + d2_Bz_dz2(i,j,k) ) / Rm
    enddo
  enddo
enddo   

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      ! Momentum Equation Evolution.
      ux(i,j,k) = ux(i,j,k) + ( (3.0d0/2.0d0)*dux_dt_new(i,j,k) - (1.0d0/2.0d0)*dux_dt_old(i,j,k) )*dt
      uy(i,j,k) = uy(i,j,k) + ( (3.0d0/2.0d0)*duy_dt_new(i,j,k) - (1.0d0/2.0d0)*duy_dt_old(i,j,k) )*dt
      uz(i,j,k) = uz(i,j,k) + ( (3.0d0/2.0d0)*duz_dt_new(i,j,k) - (1.0d0/2.0d0)*duz_dt_old(i,j,k) )*dt
      ! Magnetic Field Equation Evolution.
      Bx(i,j,k) = Bx(i,j,k) + ( (3.0d0/2.0d0)*dBx_dt_new(i,j,k) - (1.0d0/2.0d0)*dBx_dt_old(i,j,k) )*dt
      By(i,j,k) = By(i,j,k) + ( (3.0d0/2.0d0)*dBy_dt_new(i,j,k) - (1.0d0/2.0d0)*dBy_dt_old(i,j,k) )*dt
      Bz(i,j,k) = Bz(i,j,k) + ( (3.0d0/2.0d0)*dBz_dt_new(i,j,k) - (1.0d0/2.0d0)*dBz_dt_old(i,j,k) )*dt
    enddo
  end do
end do
      
do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      dux_dt_old(i,j,k) = dux_dt_new(i,j,k)
      duy_dt_old(i,j,k) = duy_dt_new(i,j,k)
      duz_dt_old(i,j,k) = duz_dt_new(i,j,k)
      dBx_dt_old(i,j,k) = dBx_dt_new(i,j,k)
      dBy_dt_old(i,j,k) = dBy_dt_new(i,j,k)
      dBz_dt_old(i,j,k) = dBz_dt_new(i,j,k)
    enddo
  end do
end do

  call dcsftdy(Nx, Ny, Nz, ux, dux_dy)
  call dcsftdz(Nx, Ny, Nz, ux, dux_dz)
  call dcsftdx(Nx, Ny, Nz, uy, duy_dx)
  call dcsftdz(Nx, Ny, Nz, uy, duy_dz)
  call dcsftdx(Nx, Ny, Nz, uz, duz_dx)
  call dcsftdy(Nx, Ny, Nz, uz, duz_dy)
  
  call dcsftdx(Nx, Ny, Nz, Bx, dBx_dx)
  call dcsftdy(Nx, Ny, Nz, Bx, dBx_dy)
  call dcsftdz(Nx, Ny, Nz, Bx, dBx_dz)
  call dcsftdx(Nx, Ny, Nz, By, dBy_dx)
  call dcsftdy(Nx, Ny, Nz, By, dBy_dy)  
  call dcsftdz(Nx, Ny, Nz, By, dBy_dz)
  call dcsftdx(Nx, Ny, Nz, Bz, dBz_dx)
  call dcsftdy(Nx, Ny, Nz, Bz, dBz_dy)
  call dcsftdz(Nx, Ny, Nz, Bz, dBz_dz)
    
do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
    u2(i,j,k) = ux(i,j,k)*ux(i,j,k) + uy(i,j,k)*uy(i,j,k) + uz(i,j,k)*uz(i,j,k)
    B2(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)
    wx(i,j,k) = duz_dy(i,j,k) - duy_dz(i,j,k)
    wy(i,j,k) = dux_dz(i,j,k) - duz_dx(i,j,k)
    wz(i,j,k) = duy_dx(i,j,k) - dux_dy(i,j,k)
    w2(i,j,k) = wx(i,j,k)*wx(i,j,k) + wy(i,j,k)*wy(i,j,k) + wz(i,j,k)*wz(i,j,k)
    jx(i,j,k) = dBz_dy(i,j,k) - dBy_dz(i,j,k)
    jy(i,j,k) = dBx_dz(i,j,k) - dBz_dx(i,j,k)
    jz(i,j,k) = dBy_dx(i,j,k) - dBx_dy(i,j,k)
    j2(i,j,k) = jx(i,j,k)*jx(i,j,k) + jy(i,j,k)*jy(i,j,k) + jz(i,j,k)*jz(i,j,k)
    dB(i,j,k) = dBx_dx(i,j,k) + dBy_dy(i,j,k) + dBz_dz(i,j,k)
    enddo
  end do
end do
  
  call dcsftdpoi(Nx, Ny, Nz, jx, Ax)
  call dcsftdpoi(Nx, Ny, Nz, jy, Ay)
  call dcsftdpoi(Nx, Ny, Nz, jz, Az)
  
if (mod(t,int(1.0/dt)) == 0) then
  do i = 1,Nx
    do j = 1,Ny
      do k = 1,Nz
      write(t+100,*) x(i),y(j),z(k),ux(i,j,k),uy(i,j,k),uz(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
      enddo
    end do
  end do
  close(t+100)
endif

Energy = 0.0d0; B_Field = 0.0d0; C1 = 0.0d0; C2 = 0.0d0; C3 = 0.0d0; Hf = 0.0d0; Hm = 0.0d0; HuB = 0.0d0; HAw = 0.0d0;
R1 = 0.0d0; R2 = 0.0d0; Rayleigh = 0.0d0; dBT = 0.0d0

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
    Energy = Energy + u2(i,j,k) / dfloat(Nx*Ny*Nz)
    B_Field = B_Field + B2(i,j,k) / dfloat(Nx*Ny*Nz)
    C1 = C1 + dsqrt( w2(i,j,k) ) / dfloat(Nx*Ny*Nz)
    C2 = C2 + w2(i,j,k) / dfloat(Nx*Ny*Nz)
    C3 = C3 + dsqrt( w2(i,j,k) )**3.0d0 / dfloat(Nx*Ny*Nz)
    Hf = Hf + ( ux(i,j,k)*wx(i,j,k) + uy(i,j,k)*wy(i,j,k) + uz(i,j,k)*wz(i,j,k) ) / dfloat(Nx*Ny*Nz)
    Hm = Hm + ( Ax(i,j,k)*Bx(i,j,k) + Ay(i,j,k)*By(i,j,k) + Az(i,j,k)*Bz(i,j,k) ) / dfloat(Nx*Ny*Nz)
    HuB = HuB + ( ux(i,j,k)*Bx(i,j,k) + uy(i,j,k)*By(i,j,k) + uz(i,j,k)*Bz(i,j,k) ) / dfloat(Nx*Ny*Nz)
    HAw = HAw + ( Ax(i,j,k)*wx(i,j,k) + Ay(i,j,k)*wy(i,j,k) + Az(i,j,k)*wz(i,j,k) ) / dfloat(Nx*Ny*Nz)
    R1 = R1 + w2(i,j,k) + 0.50d0 * j2(i,j,k) 
    R2 = R2 + u2(i,j,k) + 0.50d0 * B2(i,j,k)
    dBT = dBT + dB(i,j,k) / dfloat(Nx*Ny*Nz)
    enddo
  end do
end do

Rayleigh = R1 / R2

write (10,*) time, Energy, B_Field, C1, C2, C3, Hf, Hm, HuB, HAw, Rayleigh, dBT
  call flush(10)
  
enddo ! time

end

!===========================================================================================================

! Wrapper Subroutine to call fftw c functions for 3d arrays
    subroutine dcsftdx(n1, n2, n3, a, a1)
    implicit none
    integer n1, n2, n3, i1, i2, i3
    real*8 a(n1, n2, n3), a1(n1, n2, n3), b(n1*n2*n3)
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                b(((i1-1)*n2+(i2-1))*n3+i3)=a(i1,i2,i3)
            enddo
        enddo
    enddo
    
    call dcst3dx(n1, n2, n3, b)
    
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                a1(i1,i2,i3)=b(((i1-1)*n2+(i2-1))*n3+i3)
            enddo
        enddo
    enddo
    end subroutine
!======================================================
    subroutine dcsftdy(n1, n2, n3, a, a1)
    implicit none
    integer n1, n2, n3, i1, i2, i3
    real*8 a(n1, n2, n3), a1(n1, n2, n3), b(n1*n2*n3)
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                b(((i1-1)*n2+(i2-1))*n3+i3)=a(i1,i2,i3)
            enddo
        enddo
    enddo
    
    call dcst3dy(n1, n2, n3, b)
    
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                a1(i1,i2,i3)=b(((i1-1)*n2+(i2-1))*n3+i3)
            enddo
        enddo
    enddo
    end subroutine
!======================================================
    subroutine dcsftdz(n1, n2, n3, a, a1)
    implicit none
    integer n1, n2, n3, i1, i2, i3
    real*8 a(n1, n2, n3), a1(n1, n2, n3), b(n1*n2*n3)
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                b(((i1-1)*n2+(i2-1))*n3+i3)=a(i1,i2,i3)
            enddo
        enddo
    enddo
    
    call dcst3dz(n1, n2, n3, b)
    
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                a1(i1,i2,i3)=b(((i1-1)*n2+(i2-1))*n3+i3)
            enddo
        enddo
    enddo
    end subroutine
!======================================================
    subroutine dcsftdxx(n1, n2, n3, a, a1)
    implicit none
    integer n1, n2, n3, i1, i2, i3
    real*8 a(n1, n2, n3), a1(n1, n2, n3), b(n1*n2*n3)
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                b(((i1-1)*n2+(i2-1))*n3+i3)=a(i1,i2,i3)
            enddo
        enddo
    enddo
    
    call dcst3dxx(n1, n2, n3, b)
    
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                a1(i1,i2,i3)=b(((i1-1)*n2+(i2-1))*n3+i3)
            enddo
        enddo
    enddo
    end subroutine
!======================================================
    subroutine dcsftdyy(n1, n2, n3, a, a1)
    implicit none
    integer n1, n2, n3, i1, i2, i3
    real*8 a(n1, n2, n3), a1(n1, n2, n3), b(n1*n2*n3)
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                b(((i1-1)*n2+(i2-1))*n3+i3)=a(i1,i2,i3)
            enddo
        enddo
    enddo
    
    call dcst3dyy(n1, n2, n3, b)
    
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                a1(i1,i2,i3)=b(((i1-1)*n2+(i2-1))*n3+i3)
            enddo
        enddo
    enddo
    end subroutine
!======================================================
    subroutine dcsftdzz(n1, n2, n3, a, a1)
    implicit none
    integer n1, n2, n3, i1, i2, i3
    real*8 a(n1, n2, n3), a1(n1, n2, n3), b(n1*n2*n3)
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                b(((i1-1)*n2+(i2-1))*n3+i3)=a(i1,i2,i3)
            enddo
        enddo
    enddo
    
    call dcst3dzz(n1, n2, n3, b)
    
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                a1(i1,i2,i3)=b(((i1-1)*n2+(i2-1))*n3+i3)
            enddo
        enddo
    enddo
    end subroutine
!======================================================
    subroutine dcsftdpoi(n1, n2, n3, a, a1)
    implicit none
    integer n1, n2, n3, i1, i2, i3
    real*8 a(n1, n2, n3), a1(n1, n2, n3), b(n1*n2*n3)
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                b(((i1-1)*n2+(i2-1))*n3+i3)=a(i1,i2,i3)
            enddo
        enddo
    enddo
    
    call dcst3dpoi(n1, n2, n3, b)
    
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                a1(i1,i2,i3)=b(((i1-1)*n2+(i2-1))*n3+i3)
            enddo
        enddo
    enddo
    end subroutine
!====================================================== 
   subroutine dcsftdxy(n1, n2, n3, a, a1)
    implicit none
    integer n1, n2, n3, i1, i2, i3
    real*8 a(n1, n2, n3), a1(n1, n2, n3), b(n1*n2*n3)
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                b(((i1-1)*n2+(i2-1))*n3+i3)=a(i1,i2,i3)
            enddo
        enddo
    enddo
    
    call dcst3dxy(n1, n2, n3, b)
    
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                a1(i1,i2,i3)=b(((i1-1)*n2+(i2-1))*n3+i3)
            enddo
        enddo
    enddo
    end subroutine
!======================================================
    subroutine dcsftdyz(n1, n2, n3, a, a1)
    implicit none
    integer n1, n2, n3, i1, i2, i3
    real*8 a(n1, n2, n3), a1(n1, n2, n3), b(n1*n2*n3)
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                b(((i1-1)*n2+(i2-1))*n3+i3)=a(i1,i2,i3)
            enddo
        enddo
    enddo
    
    call dcst3dyz(n1, n2, n3, b)
    
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                a1(i1,i2,i3)=b(((i1-1)*n2+(i2-1))*n3+i3)
            enddo
        enddo
    enddo
    end subroutine
!======================================================
    subroutine dcsftdxz(n1, n2, n3, a, a1)
    implicit none
    integer n1, n2, n3, i1, i2, i3
    real*8 a(n1, n2, n3), a1(n1, n2, n3), b(n1*n2*n3)
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                b(((i1-1)*n2+(i2-1))*n3+i3)=a(i1,i2,i3)
            enddo
        enddo
    enddo
    
    call dcst3dxz(n1, n2, n3, b)
    
    do i1=1,n1
        do i2=1,n2
            do i3=1,n3
                a1(i1,i2,i3)=b(((i1-1)*n2+(i2-1))*n3+i3)
            enddo
        enddo
    enddo
    end subroutine
!======================================================

