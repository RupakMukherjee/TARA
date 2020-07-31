Program rupak
implicit none

integer ( kind = 4 ), parameter :: Nx = 512
integer ( kind = 4 ), parameter :: Ny = 512
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

include "fftw3.f"

integer ( kind = 4 ) i,j,t
real ( kind = 8 ) Lx,Ly,dx,dy,kx,ky,time,time_min,time_max,dt,nux,nuy,C0,C1,C2,C3,W0,m,d,G,Gb
real ( kind = 8 ) Energy,Energy0,y_Energy,y_Energy0
real ( kind = 8 ) x(Nx),y(Ny),psi(Nx,Ny),omega(Nx,Ny),omega_dum(Nx,Ny),ux(Nx,Ny),uy(Nx,Ny),ux_dum(Nx,Ny),uy_dum(Nx,Ny)
complex ( kind = 8 ) omegak(Nh,Ny),omegak_dum(Nh,Ny),omegak_new(Nh,Ny),dt_omegak_old(Nh,Ny),dt_omegak_new(Nh,Ny)
complex ( kind = 8 ) psik(Nh,Ny),psik_dum(Nh,Ny),ukx(Nh,Ny),uky(Nh,Ny),ukx_dum(Nh,Ny),uky_dum(Nh,Ny),Ek(Nh,Ny)
integer ( kind = 8 ) plan_forward,plan_backward
real ( kind = 8 ) t1,t2
open(unit=5,file='System_information.dat',status='unknown')
open(unit=10,file='Initial_Condition_Omega.dat',status='unknown')
open(unit=15,file='Initial_Condition_Velocity.dat',status='unknown')
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
time_max = 50.0d0 
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
    write(15,*) x(i),y(j),omega(i,j),ux(i,j),uy(i,j)
  end do
end do

!===================== INITIAL TIME DATA ===============================================

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, omega_dum, omegak, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

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

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, psik_dum, psi, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_destroy_plan_ (plan_forward)
  call dfftw_destroy_plan_ (plan_backward)

Energy0 = 0.0d0
y_Energy0 = 0.0d0
 C0 = 0.0d0
 C1 = 0.0d0
 C2 = 0.0d0
 C3 = 0.0d0
 
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
  write(10,*) x(i),y(j),omega(i,j),psi(i,j)
  write(20,*) x(i),y(j),omega(i,j),ux(i,j),uy(i,j)
  end do
end do

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

do i = 1,Nx/2+1
  do j = 1,Ny
    dt_omegak_old(i,j) = dt_omegak_new(i,j)
    omegak(i,j) = omegak_new(i,j)
  enddo
enddo

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

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_destroy_plan_ (plan_backward)

!if (t == 10000) then 
!do i = 1,Nx/2+1
!  do j = 1,Ny
!  Ek(i,j) = Ek(i,j) + sqrt(abs(ukx(i,j))**2 + abs(uky(i,j))**2)
!  write(35,*) i,j,sqrt(float((i-1)**2)+float((j-1)**2)),abs(Ek(i,j))
!  enddo
!enddo  
!endif 
 
Energy = 0.0d0 
y_Energy = 0.0d0
 C0 = 0.0d0
 C1 = 0.0d0
 C2 = 0.0d0
 C3 = 0.0d0
   
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

write(40,*) time,Energy,y_Energy/(2.0d0*dfloat(Nx)*dfloat(Ny))!/y_Energy0!, C0,C1,C2,C3 
  
enddo ! time
call cpu_time ( t2 )
write(5,*) "Time taken for the run = ",t2 - t1
 close(5)
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

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    i_kx_omegak(i,j) = (0.0d0,1.0d0)*kx*omegak_dum(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0)*ky*omegak_dum(i,j)
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    i_kx_omegak(i,j) = (0.0d0,1.0d0)*kx*omegak_dum(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0)*ky*omegak_dum(i,j)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_omegak, domega_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_omegak, domega_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
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

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_domega_dx, NLkx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_domega_dy, NLky, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  
  call dfftw_destroy_plan_ (plan_forward)
  call dfftw_destroy_plan_ (plan_backward)
  
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
  
do i = 1,Nx/2+1
  do j = 1,Ny/2  
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    NLk(i,j) = ( NLkx(i,j) + NLky(i,j) )
    Nlk(i,j) = NLk(i,j) + (nux*kx*kx + nuy*ky*ky)*omegak(i,j)
  enddo
  do j = Ny/2+1,Ny  
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    NLk(i,j) = ( NLkx(i,j) + NLky(i,j) )
    Nlk(i,j) = NLk(i,j) + (nux*kx*kx + nuy*ky*ky)*omegak(i,j)
  enddo
enddo   

do i = 1,Nx/2+1
  do j = 1,Ny  
    dt_omegak_new(i,j) = -NLk(i,j)
  enddo
enddo   

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

do i = 1,Nh
  do j = 1,Ny
    omegak_new(i,j) = omegak(i,j) + ( (3.0d0/2.0d0)*dt_omegak_new(i,j) - (1.0d0/2.0d0)*dt_omegak_old(i,j) )*dt
  end do
end do

return

end subroutine ab

!====================================================================================

end program rupak  

