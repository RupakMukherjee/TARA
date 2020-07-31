! To run the program : gfortran -I/usr/include <filename.f95> -lfftw3; ./a.out
! It uses Adam Bashforth Solver
program Burgers_Turbulence_AB
implicit none

include "fftw3.f"
	
double precision :: x,ux,dux,ddux
dimension x(1024),ux(1024),dux(1024),ddux(1024)
double precision :: abs_uk,a,b
dimension abs_uk(513),a(513),b(513)
double complex :: uk,duk,dduk,uk_backup
dimension uk(513),duk(513),dduk(513),uk_backup(513)
integer :: i
double precision :: dx
double precision, parameter :: pi=3.14159265358979323846d0
integer*8 :: plan
real*8 :: L

real*8 f,Random,time,time_min,time_max,dt,nu
real*8, dimension (20000,1024) :: nonlinear,force
integer*8 N,t
external f

integer,parameter :: seed = 99999999
call srand(seed)

		
!open(unit=10,file='fourier_data.dat',status='unknown')
!open(unit=12,file='nonlinear.dat',status='unknown')	
	
!L = 5.0*pi+0.5
L = 2.0*pi
dx = L/1024.0
   
N = 1024
time_min = 0.00d0
time_max = 2.00d0
dt = 0.00010d0

!nu = 0.0380d0
nu = 0.0d0

do i=1,N
x(i)=0.0+real(i-1)*dx
ux(i)= f(x(i))   
enddo

do time = time_min,time_max,dt

t = nint(time/dt) - int(time_min/dt)
!print*, t, time	
call dfftw_plan_dft_r2c_1d(plan,N,ux,uk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,ux,uk)
call dfftw_destroy_plan(plan)
	
do i=1,N/2+1
abs_uk(i)=abs(uk(i))
!if (i .ge. N/3+1) then     ! De - Aliazing Technique...
!uk(i) = (0.0d0,0.0d0)
!endif
!if (t == 10000) then 
!write(10,*) t, i-1.0, 2.0*abs_uk(i)/float(N)
!endif
enddo
	
do i=1,N/2+1
uk_backup(i) = uk(i)
duk(i) = uk(i) * (i-1) * (0.d0,1.d0)
dduk(i) = - uk(i) * (i-1) * (i-1)
enddo
if (mod(N,2)==0) uk(N/2)=cmplx(0.d0)
	
call dfftw_plan_dft_c2r_1d(plan,N,duk,dux,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,duk,dux)
call dfftw_destroy_plan(plan)
	
call dfftw_plan_dft_c2r_1d(plan,N,dduk,ddux,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,dduk,ddux)
call dfftw_destroy_plan(plan)

do i=1,N
dux(i)=2.0*pi*dux(i)/(L*float(N))
ddux(i)=2.0*pi*ddux(i)/(L*float(N))
if (t .ge. 1000 .and. mod(t,500) == 0) then
!if (t .ge. 11600 .and. mod(t,10) == 0 .and. t .le. 11700) then
write (t,*) t, x(i),ux(i)!,dux(i),ddux(i)
endif
enddo


do i = 1,N,1
!nonlinear(t,i) = 1.0*dux(i)                                  ! Acoustic Wave...
!nonlinear(t,i) = ux(i)*dux(i)                                ! Inviscid Burger Equation...
nonlinear(t,i) = ux(i)*dux(i) - nu*ddux(i)                    ! Burger Equation...
!nonlinear(t,i) = (nonlinear(t,i) + nonlinear(t-1,i))/2.0     ! Time Delay... 
!nonlinear(t,i) =  - nu*ddux(i)                               ! Diffusion Euation
!write(12,*) x(i),nonlinear(t,i)                 
enddo


do i = 1,N,1
Random = 0.0d0
!Random = rand()
force(t,i) = -nonlinear(t,i) + Random
ux(i) = ux(i) + ((3.0/2.0)*force(t,i) - (1.0/2.0)*force(t-1,i))*(dt)  
enddo
enddo

 close(10)
 close(12)

end program Burgers_Turbulence_AB



!==============================================================

function f(x)
implicit none
real*8 f,pi,x,L
pi = atan(1.0d0)*4.0d0
!L = 5.0*pi+0.5
L = 2.0*pi
f = sin(1.0*x)
!f= (x-L/2.0)**2
!f = sin(1.0*x) + sin(2.0*x+0.9) + sin(3.0*x)
!f = exp(-(x-L/2.0)**2)
!if (x .ge. 2.0 .and. x .le. 3.0) then
!f = 1.0
!else
!f = 0.0
!endif
return
end

!================================================================

