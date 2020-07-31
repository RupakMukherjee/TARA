! To run the program : gfortran -I/usr/include <filename.f95> -lfftw3; ./a.out
! It uses RK4 Solver
program Burgers_Turbulence_RK4
implicit none

include "fftw3.f"
	
double precision :: x,ux,NLx
dimension x(1024),ux(1024),NLx(1024)
double complex :: uk,uk_dum,duk,NLk,dtuk
dimension uk(513),uk_dum(513),duk(513),NLk(513),dtuk(513)
integer*8 :: i,j,k
double precision :: dx
double precision, parameter :: pi=3.14159265358979323846d0
integer*8 :: plan
real*8 :: L

real*8 time,time_min,time_max,dt
integer*8 N,t,Pint

!===================== USER INPUTS ============================================		

N = 1024
L = 2.0*pi
dx = L/float(N)

time_min = 0.0d0
time_max = 2.0d0 ! 2.0d0
dt = 0.00010d0    ! 0.00010d0

do i=1,N
x(i)=0.0+real(i-1)*dx
ux(i)= sin(1.0*x(i))
write(12,*) 0, x(i),ux(i)
enddo

!==================== MAIN PROGRAM ============================================

do time = time_min,time_max,dt

t = nint(time/dt) - int(time_min/dt)

call dfftw_plan_dft_r2c_1d(plan,N,ux,uk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,ux,uk)
call dfftw_destroy_plan(plan)

!Calculating the time evolution in k-space...
	
call derive(N,time,uk,NLk,dtuk)
call rk4(N,time,dt,uk,NLk,dtuk)

!Calculating the final velocity...

call dfftw_plan_dft_c2r_1d(plan,N,uk,ux,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,uk,ux)
call dfftw_destroy_plan(plan)

do i=1,N
ux(i)=2.0*pi*ux(i)/(L*float(N))
if (t .ge. 10000 .and. mod(t,100) == 0 .and. t .le. 11000) then
write (t,*) t, x(i),ux(i)
endif
enddo

enddo ! time

contains
!===================================================================================

subroutine derive(N,time,uk,NLk,dtuk)
real*8 pi,time
double precision :: ux,NLx
dimension ux(1024),NLx(1024)
double complex :: uk,uk_dum,dtuk,NLk
dimension uk(513),uk_dum(513),dtuk(513),NLk(513)
integer*8 i,N

pi = 4.0*atan(1.0)

do i = 1,N/2+1
uk_dum(i) = uk(i)/float(N)
enddo

call dfftw_plan_dft_c2r_1d(plan,N,uk_dum,ux,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan,uk_dum,ux)
call dfftw_destroy_plan(plan)

do i = 1,N
NLx(i) = ux(i)*ux(i)   
enddo

call dfftw_plan_dft_r2c_1d(plan,N,NLx,NLk,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan,NLx,NLk)
call dfftw_destroy_plan(plan)

do i = 1,N/2+1
NLk(i) = NLk(i) * (i-1) * (0.0d0,1.0d0) / 2.0
enddo

do i = 1,N/2
dtuk(i) = -NLk(i)
enddo

return
end 

!====================================================================================

subroutine rk4(N,time,dt,uk,NLk,dtuk)
real*8 time,dt
double complex :: uk,dtuk,NLk,k1,k2,k3,k4,dum
dimension uk(513),dtuk(513),NLk(513),k1(513),k2(513),k3(513),k4(513),dum(513)
integer*8 i,N
N = 1024

do i = 1,N/2
k1(i) = dtuk(i)
dum(i) = uk(i) + k1(i)*dt/2.0d0
enddo

call derive(N,time+dt/2.0d0,dum,NLk,dtuk)

do i = 1,N/2
k2(i) = dtuk(i)
dum(i) = uk(i) + k2(i)*dt/2.0d0
enddo

call derive(N,time+dt/2.0d0,dum,NLk,dtuk)

do i = 1,N/2
k3(i) = dtuk(i)
dum(i) = uk(i) + k3(i)*dt
enddo

call derive(N,time+dt,dum,NLk,dtuk)

do i = 1,N/2
k4(i) = dtuk(i)
uk(i) = uk(i) + dt/6.0d0*(k1(i) + 2.0d0*k2(i) + 2.0d0*k3(i) + k4(i))
enddo

return
end subroutine rk4

!=====================================================================================

end program Burgers_Turbulence_RK4
