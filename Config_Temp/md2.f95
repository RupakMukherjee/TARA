program molecular_dynamics
implicit none
real*8 pi,a,M,epsilon,wpd,L,Lx,Ly,Vxmax,Vymax,U1,U2,mu1,mu2,Sigma1,Sigma2,svx,svy,DGK,DES,VACFnum,VACFden,VACF
real*8 nbar,r,Gamma,Temp,K,fx,fy,KE,KEH,PE,TE,dP,BerendsenT,BerendsenP,Gaussian,Nose_Hoover,Config_Therm,MicroCan,Magnetic,B0
real*8 t,tmax,dt,xdiff,ydiff,scale,tauT,tauP,sumvx,sumvy,numx,denx,numy,deny,alphax,alphay,xi,xih,Q
real*8 virialx,virialy,virial,GF,force,gradforce,ConfigT,P0,Pressure,mu
integer*8 i,j,N,p,time
real*8, dimension (100000):: x,y,vx,vy,ax,ay,ux,uy,v,ac,x0,y0,vx0,vy0,x1,y1,vx1,vy1,GFT
integer,parameter :: seed = 99999999
pi = 4.0d0*atan(1.0d0)
call srand(seed)

!====================== User Inputs =============================

!Coupling Parameter...
Gamma = 170.0d0
Temp = 1.0d0/Gamma

!Screening Parameter...
K = 1.0d0

!Total Number of Particles...
N = 1225

!Areal*8 Density of Particles...
!nbar = float(N)/(2.0d0*Lx*2.0d0*Ly)
nbar = 1.0d0/pi

!Interparticle Distance...
a = sqrt(1.0d0/(pi*nbar))

!Normalized System Size...
L = sqrt(float(N)/(4.0d0*nbar))
Lx = L/a
Ly = L/a

!Initial Normalized Maximum Velocity of the Particles
Vxmax = 1.0d0*sqrt(2.0d0)/a
Vymax = 1.0d0*sqrt(2.0d0)/a
mu1 = 0.0d0
mu2 = 0.0d0
Sigma1 = sqrt(Temp)
Sigma2 = sqrt(Temp)

!Normalized Final Time and Time Steps...
tmax = 100.0d0!/sqrt(2.0d0)
dt = 0.010d0!/sqrt(2.0d0)
BerendsenT = tmax/2.0d0
!BerendsenP = tmax/20.0d0
!Nose_Hoover = tmax/20.0d0
!Config_Therm = tmax/20.0d0
!Gaussian = tmax/5.0d0
MicroCan = tmax
Magnetic = tmax

dP = (BerendsenP - BerendsenT)/4.0d0

xi = 0.0d0
time = 0

!====================== Output Filenames ===========================

open(unit=1,file='Initial_Configuration.dat',status='unknown')
!open(unit=2,file='Average_Velocity.dat',status='unknown')
open(unit=10,file='Energy.dat',status='unknown')
open(unit=20,file='Diff_VACF.dat',status='unknown')
open(unit=30,file='Pressure.dat',status='unknown')

!====================== Definition of initial state =================

call Initial_configuration (N,Lx,Ly,Sigma1,Sigma2,mu1,mu2,x,y,vx,vy,ux,uy,PE,KE)

!Initial Position, Velocity...
do i = 1,N,1
write(1,*) 0.0d0,x(i),y(i),vx(i),vy(i)
enddo

!Total Energy...
TE = KE + PE
write(10,*) 0.0d0,TE/float(N),KE/float(N),PE/float(N),0.0d0

!============== MD Time Evolution using velocity verlet method started ================

do t = dt,tmax,dt
B0 = 1.0d0
if (t .le. BerendsenT) then 
call Berendsen_Thermostat (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)

!elseif (t .gt. BerendsenT .and. t .le. BerendsenP) then
!call Berendsen_Thermostat (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)
!call Berendsen_Barostat (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure,BerendsenT,dP)

!elseif (t .gt. BerendsenP .and. t .le. Nose_Hoover) then
!call Nose_Hoover_Thermostat (Temp,K,N,t,dt,xi,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)

!elseif (t .gt. Nose_Hoover .and. t .le. Config_Therm) then
!call Configurational_Thermostat (Temp,K,N,t,dt,xi,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)

!elseif (t .gt. Config_Therm .and. t .le. Gaussian) then
!call Gaussian_Thermostat (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)

!elseif (t .gt. Gaussian .and. t .le. MicroCan) then
elseif (t .gt. BerendsenT) then 
call Micro_Canonical (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)
!call Diffusion_VACF (N,t,dt,x,y,vx,vy,Gaussian,time,DGK,DES,VACF) 
!write(20,*) time,DGK/(2.0d0*float(N)),DES/(4.0d0*float(N)*time),VACF

!elseif (t .gt. MicroCan .and. t .le. Magnetic) then
!call Magnetic_Field (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure,B0)

endif

!Total Energy...
TE = KE + PE

write(10,*) t,TE/float(N),KE/float(N),PE/float(N),1.0d0/ConfigT

write(30,*) t,Pressure,Lx,Ly

do i = 1,N,1
p = int(t/dt)
if (p  .ge. 5000 .and. mod(float(p),10.0d0) == 0.0d0) then
write(p,*) t,x(i),y(i),vx(i),vy(i)
endif
enddo

enddo

contains

!=======================================================================================







!       I N I T I A L   C O N F I G U R A T I O N 







!=======================================================================================

subroutine Initial_Configuration (N,Lx,Ly,Sigma1,Sigma2,mu1,mu2,x,y,vx,vy,ux,uy,PE,KE)
implicit none
real*8, intent (in) :: Lx,Ly,Sigma1,Sigma2,mu1,mu2
integer*8, intent (in) :: N
real*8, intent (out) :: PE,KE
real*8, dimension (10000), intent (out) :: x,y,vx,vy,ux,uy
real*8 U1,U2,svx,svy,xdiff,ydiff,fx,fy
real*8, dimension (10000) :: ax,ay
integer,parameter :: seed = 99999999
pi = 4.0d0*atan(1.0d0)
call srand(seed)

svx = 0.0d0
svy = 0.0d0
KE = 0.0d0
PE = 0.0d0

!Definition of the initial random positions and gaussian velocities in -Lx to Lx and -Ly to Ly rectangular box... 
do i = 1, N, 1
x(i) = (rand())*2.0d0*Lx - Lx
y(i) = (rand())*2.0d0*Ly - Ly
U1 = rand()
U2 = rand()
vx(i) = sqrt(-2.0d0*log(U1))*cos(2.0d0*pi*U2)*Sigma1 + mu1          ! Gaussian Velocity Distribution in x...
!vx(i) = (rand())*Vxmax - Vxmax/2.0d0
svx = svx + vx(i)                                               ! Center of Mass has Zero x Velocity...
vy(i) = sqrt(-2.0d0*log(U1))*sin(2.0d0*pi*U2)*Sigma2 + mu2          ! Gaussian Velocity Distribution in y...
!vy(i) = (rand())*Vymax - Vymax/2.0d0
svy = svy + vy(i)                                               ! Center of Mass has Zero y Velocity...
enddo 

!Definitions of corrected initial velocities...
do i = 1,N,1
vx(i) = vx(i) - svx/float(N)
vy(i) = vy(i) - svy/float(N)
KE = KE + (vx(i)**2 + vy(i)**2)/2.0d0                             ! Calculation of Kinetic Energy...
enddo

!write(2,*) 0.0d00,svx,svy

end subroutine Initial_Configuration

!===================================================================================







!       B E R E N D S E N   T H E R M O S T A T







!===================================================================================


subroutine Berendsen_Thermostat (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)
implicit none
real*8, intent (in) :: Temp,K,t,dt,Lx,Ly
real*8, intent (out) :: PE,KE,Pressure,ConfigT
integer*8, intent (in) :: N
real*8, dimension(10000), intent (inout) :: x,y,vx,vy,ux,uy
real*8 xdiff,ydiff,r,fx,fy,sumvx,sumvy,virialx,virialy,virial,tauT,scale
real*8, dimension (10000) :: ax,ay
integer*8 i,j

sumvx = 0.0d0
sumvy = 0.0d0
PE = 0.0d0
KE = 0.0d0
force = 0.0d0
gradforce = 0.0d0
virial = 0.0d0

!Calculating the acceleration at time dt

do i = 1,N,1
ax(i) = 0.0d0
ay(i) = 0.0d0
GFT(i) = 0.0d0
do j = 1,N,1
if (i .ne. j) then                                          ! .ne. is used, PE should be halved...
xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0d0*Lx))*2.0d0*Lx     ! Minimum Image Convension Introduced...
ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0d0*Ly))*2.0d0*Ly     ! Minimum Image Convension Introduced...
r = sqrt((xdiff)**2 + (ydiff)**2)
fx = (xdiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ax(i) = ax(i) + fx
virialx = xdiff*fx
fy = (ydiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ay(i) = ay(i) + fy
virialy = ydiff*fy
virial = virial + (virialx + virialy)/2.0d0                   ! Calculation of virial...
GF = (1.0d0 + K*r + K*K*r*r)*exp(-K*r)/r**3
GFT(i) = GFT(i) + GF
PE = PE + (exp(-k*r)/(2.0d0*r))                               ! Calculation of Potential Energy...
endif
enddo
force = force + ax(i)**2 + ay(i)**2
gradforce = gradforce + GFT(i)
enddo

ConfigT = force/gradforce                                   ! Calculation of Configurational Temperature...

!Calculating the velocity at time dt 
	  
do i = 1,N,1
sumvx = sumvx + vx(i)                                       ! Check for average x-velocity...
sumvy = sumvy + vy(i)                                       ! Check for average y-velocity...
vx(i) = ux(i) + ax(i)*dt/2.0d0
vy(i) = uy(i) + ay(i)*dt/2.0d0
KE = KE + (vx(i)**2 + vy(i)**2)/2.0d0                         ! Calculation of Kinetic Energy...
enddo

!write(2,*) t,sumvx,sumvy

Pressure = (KE + 0.5*virial)/(4.0d0*Lx*Ly)

!Berendsen Thermostat 

tauT = 10.0d0*dt
scale = sqrt(1.0d0 + (dt/tauT)*((Temp/(KE/float(N)))-1.0d0))

do i = 1,N,1
ux(i) = scale*ux(i)
uy(i) = scale*uy(i)
enddo

do i = 1,N,1
ux(i) = ux(i) + dt*ax(i)  
x(i) = x(i) + dt*ux(i)
x(i) = x(i) - (int(x(i)/Lx))*2.0d0*Lx                           ! Periodic Boundary Condition in x...

uy(i) = uy(i) + dt*ay(i)
y(i) = y(i) + dt*uy(i)
y(i) = y(i) - (int(y(i)/Ly))*2.0d0*Ly                           ! Periodic Boundary Condition in y...

enddo

end subroutine Berendsen_Thermostat

!=================================================================================







!         B E R E N D S E N   B A R O S T A T 







!=================================================================================

subroutine Berendsen_Barostat (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure,BerendsenT,dP)
implicit none
real*8, intent (in) :: Temp,K,t,dt,BerendsenT,dP
real*8, intent (out) :: PE,KE,Pressure,ConfigT
real*8, intent (inout) :: Lx,Ly
integer*8, intent (in) :: N
real*8, dimension(10000), intent (inout) :: x,y,vx,vy,ux,uy
real*8 pi,P0,xdiff,ydiff,r,fx,fy,sumvx,sumvy
real*8 virialx,virialy,virial,tauP,mu,GF,force,gradforce
real*8, dimension (10000) :: ax,ay
integer*8 i,j
pi = 4.0d0*atan(1.0d0)

sumvx = 0.0d0
sumvy = 0.0d0
PE = 0.0d0
KE = 0.0d0
force = 0.0d0
gradforce = 0.0d0
virial = 0.0d0

!Pressure Required...
if (t .lt. BerendsenT + dP) then
P0 = 1.0d0/pi
elseif (t .lt. BerendsenT + 2.0d0*dP) then
P0 = 2.0d0/pi
elseif (t .lt. BerendsenT + 3.0d0*dP) then
P0 = 1.0d0/(2.0d0*pi)
else
P0 = 1.0d0/pi
endif

!Calculating the acceleration at time dt

do i = 1,N,1
ax(i) = 0.0d0
ay(i) = 0.0d0
GFT(i) = 0.0d0
do j = 1,N,1
if (i .ne. j) then                                          ! .ne. is used, PE should be halved...
xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0d0*Lx))*2.0d0*Lx     ! Minimum Image Convension Introduced...
ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0d0*Ly))*2.0d0*Ly     ! Minimum Image Convension Introduced...
r = sqrt((xdiff)**2 + (ydiff)**2)
fx = (xdiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ax(i) = ax(i) + fx
virialx = xdiff*fx
fy = (ydiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ay(i) = ay(i) + fy
virialy = ydiff*fy
virial = virial + (virialx + virialy)/2.0d0                   ! Calculation of virial...
GF = (1.0d0 + K*r + K*K*r*r)*exp(-K*r)/r**3
GFT(i) = GFT(i) + GF
PE = PE + (exp(-k*r)/(2.0d0*r))                               ! Calculation of Potential Energy...
endif
enddo
force = force + ax(i)**2 + ay(i)**2
gradforce = gradforce + GFT(i)
enddo

ConfigT = force/gradforce                                   ! Calculation of Configurational Temperature...

!Calculating the velocity at time dt
	  
do i = 1,N,1
sumvx = sumvx + vx(i)                                       ! Check for average x-velocity...
sumvy = sumvy + vy(i)                                       ! Check for average y-velocity...
vx(i) = ux(i) + ax(i)*dt/2.0d0
vy(i) = uy(i) + ay(i)*dt/2.0d0
v(i) = sqrt(vx(i)**2 + vy(i)**2)
KE = KE + (vx(i)**2 + vy(i)**2)/2.0d0                         ! Calculation of Kinetic Energy...
enddo

!write(2,*) t,sumvx,sumvy

!Berendsen Barostat

tauP = 10.0d0*dt
Pressure = (KE + 0.5*virial)/(4.0d0*Lx*Ly)
mu = (1.0d0 - (dt/tauP)*(P0 - Pressure))**(1.0d0/3.0d0)

do i = 1,N,1
x(i) = mu*x(i)
y(i) = mu*y(i)
enddo
Lx = mu*Lx 
Ly = mu*Ly

do i = 1,N,1
ux(i) = ux(i) + dt*ax(i)  
x(i) = x(i) + dt*ux(i)
x(i) = x(i) - (int(x(i)/Lx))*2.0d0*Lx                           ! Periodic Boundary Condition in x...

uy(i) = uy(i) + dt*ay(i)
y(i) = y(i) + dt*uy(i)
y(i) = y(i) - (int(y(i)/Ly))*2.0d0*Ly                           ! Periodic Boundary Condition in y...

enddo

end subroutine Berendsen_Barostat

!======================================================================







!       N O S E - H O O V E R   T H E R M O S T A T 







!======================================================================

subroutine Nose_Hoover_Thermostat (Temp,K,N,t,dt,xi,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)
implicit none
real*8, intent (in) :: Temp,K,t,dt,Lx,Ly
integer*8, intent (in) :: N
real*8, intent (out) :: PE,KE,Pressure,ConfigT
real*8, intent (inout) :: xi
real*8, dimension (10000), intent (inout) :: x,y,vx,vy,ux,uy
real*8 pi,P0,xdiff,ydiff,r,fx,fy,sumvx,sumvy
real*8 virialx,virialy,virial,KEH,Q,GF,force,gradforce
real*8, dimension (10000) :: ax,ay,GFT
integer*8 i,j

sumvx = 0.0d0
sumvy = 0.0d0
PE = 0.0d0
KEH = 0.0d0
force = 0.0d0
gradforce = 0.0d0
virial = 0.0d0

Q = 6.0d0*(float(N)+0.50)*Temp

do i = 1,N,1

x(i) = x(i) + dt*vx(i) + (ax(i) - xi*vx(i))*dt*dt/2.0d0
x(i) = x(i) - (int(x(i)/Lx))*2.0d0*Lx                            ! Periodic Boundary Condition in x...
ux(i) = vx(i) + (ax(i) - xi*vx(i))*dt/2.0d0

y(i) = y(i) + dt*vy(i) + (ay(i) - xi*vy(i))*dt*dt/2.0d0
y(i) = y(i) - (int(y(i)/Ly))*2.0d0*Ly                            ! Periodic Boundary Condition in y...
uy(i) = vy(i) + (ay(i) - xi*vy(i))*dt/2.0d0

KEH = KEH + (ux(i)**2 + uy(i)**2)/2.0d0 

enddo

!Calculating the acceleration at time dt

do i = 1,N,1
ax(i) = 0.0d0
ay(i) = 0.0d0
GFT(i) = 0.0d0
do j = 1,N,1
if (i .ne. j) then                                              ! .ne. is used, PE should be halved...
xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0d0*Lx))*2.0d0*Lx     ! Minimum Image Convension Introduced...
ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0d0*Ly))*2.0d0*Ly     ! Minimum Image Convension Introduced...
r = sqrt((xdiff)**2 + (ydiff)**2)
fx = (xdiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ax(i) = ax(i) + fx
virialx = xdiff*fx
fy = (ydiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ay(i) = ay(i) + fy
virialy = ydiff*fy
virial = virial + (virialx + virialy)/2.0d0                   ! Calculation of virial...
GF = (1.0d0 + K*r + K*K*r*r)*exp(-K*r)/r**3
GFT(i) = GFT(i) + GF
PE = PE + (exp(-k*r)/(2.0d0*r))                               ! Calculation of Potential Energy...
endif
enddo
force = force + ax(i)**2 + ay(i)**2
gradforce = gradforce + GFT(i)
enddo

ConfigT = force/gradforce                                   ! Calculation of Configurational Temperature...

xih = xi + (KE - (float(N)+0.50)*Temp)*dt/(2.0d0*Q)
xi = xih + (KEH - (float(N)+0.50)*Temp)*dt/(2.0d0*Q)

!Calculating the velocity at time dt 

KE = 0.0d0
	  
do i = 1,N,1
sumvx = sumvx + vx(i)                                       ! Check for average x-velocity...
sumvy = sumvy + vy(i)                                       ! Check for average y-velocity...
vx(i) = (ux(i) + ax(i)*dt/2.0d0)/(1.0d0+xi*dt/2.0d0)
vy(i) = (uy(i) + ay(i)*dt/2.0d0)/(1.0d0+xi*dt/2.0d0)
KE = KE + (vx(i)**2 + vy(i)**2)/2.0d0                         ! Calculation of Kinetic Energy...
enddo

Pressure = (KE + 0.5*virial)/(4.0d0*Lx*Ly)

end subroutine Nose_Hoover_Thermostat

!======================================================================







!       C O N F I G U R A T I O N A L   T H E R M O S T A T 







!======================================================================

subroutine Configurational_Thermostat (Temp,K,N,t,dt,xi,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)
implicit none
real*8, intent (in) :: Temp,K,t,dt,Lx,Ly
integer*8, intent (in) :: N
real*8, intent (out) :: PE,KE,Pressure,ConfigT
real*8, intent (inout) :: xi
real*8, dimension (10000), intent (inout) :: x,y,vx,vy,ux,uy
real*8 pi,P0,xdiff,ydiff,r,fx,fy,sumvx,sumvy
real*8 virialx,virialy,virial,KEH,Q,GF,force,gradforce
real*8, dimension (10000) :: ax,ay,GFT
integer*8 i,j

sumvx = 0.0d0
sumvy = 0.0d0
PE = 0.0d0
KEH = 0.0d0
force = 0.0d0
gradforce = 0.0d0
virial = 0.0d0

Q = 6.0d0*(float(N)+0.50)*Temp

do i = 1,N,1

x(i) = x(i) + dt*vx(i) + (ax(i) - xi*vx(i))*dt*dt/2.0d0
x(i) = x(i) - (int(x(i)/Lx))*2.0d0*Lx                            ! Periodic Boundary Condition in x...
ux(i) = vx(i) + (ax(i) - xi*vx(i))*dt/2.0d0

y(i) = y(i) + dt*vy(i) + (ay(i) - xi*vy(i))*dt*dt/2.0d0
y(i) = y(i) - (int(y(i)/Ly))*2.0d0*Ly                            ! Periodic Boundary Condition in y...
uy(i) = vy(i) + (ay(i) - xi*vy(i))*dt/2.0d0

KEH = KEH + (ux(i)**2 + uy(i)**2)/2.0d0 

enddo

!Calculating the acceleration at time dt

do i = 1,N,1
ax(i) = 0.0d0
ay(i) = 0.0d0
GFT(i) = 0.0d0
do j = 1,N,1
if (i .ne. j) then                                          ! .ne. is used, PE should be halved...
xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0d0*Lx))*2.0d0*Lx     ! Minimum Image Convension Introduced...
ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0d0*Ly))*2.0d0*Ly     ! Minimum Image Convension Introduced...
r = sqrt((xdiff)**2 + (ydiff)**2)
fx = (xdiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ax(i) = ax(i) + fx
virialx = xdiff*fx
fy = (ydiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ay(i) = ay(i) + fy
virialy = ydiff*fy
virial = virial + (virialx + virialy)/2.0d0                   ! Calculation of virial...
GF = (1.0d0 + K*r + K*K*r*r)*exp(-K*r)/r**3
GFT(i) = GFT(i) + GF
PE = PE + (exp(-k*r)/(2.0d0*r))                               ! Calculation of Potential Energy...
endif
enddo
force = force + ax(i)**2 + ay(i)**2
gradforce = gradforce + GFT(i)
enddo

ConfigT = force/gradforce                                   ! Calculation of Configurational Temperature...

xih = xi + (ConfigT*float(N) - (float(N)+0.50)*Temp)*dt/(2.0d0*Q)
xi = xih + (KEH - (float(N)+0.50)*Temp)*dt/(2.0d0*Q)

!Calculating the velocity at time dt 

KE = 0.0d0
	  
do i = 1,N,1
sumvx = sumvx + vx(i)                                       ! Check for average x-velocity...
sumvy = sumvy + vy(i)                                       ! Check for average y-velocity...
vx(i) = (ux(i) + ax(i)*dt/2.0d0)/(1.0d0+xi*dt/2.0d0)
vy(i) = (uy(i) + ay(i)*dt/2.0d0)/(1.0d0+xi*dt/2.0d0)
KE = KE + (vx(i)**2 + vy(i)**2)/2.0d0                         ! Calculation of Kinetic Energy...
enddo

Pressure = (KE + 0.5*virial)/(4.0d0*Lx*Ly)

end subroutine Configurational_Thermostat

!======================================================================







!       G A U S S I A N   T H E R M O S T A T 







!======================================================================

subroutine Gaussian_Thermostat (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)
implicit none
real*8, intent (in) :: Temp,K,t,dt,Lx,Ly
integer*8, intent (in) :: N
real*8, intent (out) :: PE,KE,Pressure,ConfigT
real*8, dimension (10000), intent (inout) :: x,y,vx,vy,ux,uy
real*8 pi,P0,xdiff,ydiff,r,fx,fy,sumvx,sumvy,GF,force,gradforce
real*8 numx,numy,denx,deny,alphax,alphay,virialx,virialy,virial
real*8, dimension (10000) :: ax,ay,GFT
integer*8 i,j
pi = 4.0d0*atan(1.0d0)

sumvx = 0.0d0
sumvy = 0.0d0
PE = 0.0d0
KE = 0.0d0
numx = 0.0d0
denx = 0.0d0
numy = 0.0d0
deny = 0.0d0
force = 0.0d0
gradforce = 0.0d0
virial = 0.0d0

!Calculating the acceleration at time dt 

do i = 1,N,1
ax(i) = 0.0d0
ay(i) = 0.0d0
GFT(i) = 0.0d0
do j = 1,N,1
if (i .ne. j) then                                          ! .ne. is used, PE should be halved...
xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0d0*Lx))*2.0d0*Lx     ! Minimum Image Convension Introduced...
ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0d0*Ly))*2.0d0*Ly     ! Minimum Image Convension Introduced...
r = sqrt((xdiff)**2 + (ydiff)**2)
fx = (xdiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ax(i) = ax(i) + fx
virialx = xdiff*fx
fy = (ydiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ay(i) = ay(i) + fy
virialy = ydiff*fy
virial = virial + (virialx + virialy)/2.0d0                   ! Calculation of virial...
GF = (1.0d0 + K*r + K*K*r*r)*exp(-K*r)/r**3
GFT(i) = GFT(i) + GF
PE = PE + (exp(-k*r)/(2.0d0*r))                               ! Calculation of Potential Energy...
endif
enddo
force = force + ax(i)**2 + ay(i)**2
gradforce = gradforce + GFT(i)
enddo

ConfigT = force/gradforce                                   ! Calculation of Configurational Temperature...

!Calculating the velocity at time dt 
	  
do i = 1,N,1
sumvx = sumvx + vx(i)                                       ! Check for average x-velocity...
sumvy = sumvy + vy(i)                                       ! Check for average y-velocity...
vx(i) = ux(i) + ax(i)*dt/2.0d0
vy(i) = uy(i) + ay(i)*dt/2.0d0
numx = numx + vx(i)*ax(i)
denx = denx + vx(i)*vx(i)
numy = numy + vy(i)*ay(i)
deny = deny + vy(i)*vy(i)
KE = KE + (vx(i)**2 + vy(i)**2)/2.0d0                         ! Calculation of Kinetic Energy...
enddo

!write(2,*) t,sumvx,sumvy

Pressure = (KE + 0.5*virial)/(4.0d0*Lx*Ly)

alphax = -numx/denx                                           ! Canonical Run...
alphay = -numy/deny 

do i = 1,N,1

x(i) = x(i) + dt*vx(i) + ax(i)*dt*dt/2.0d0
x(i) = x(i) - (int(x(i)/Lx))*2.0d0*Lx                            ! Periodic Boundary Condition in x...
ux(i) = ux(i)*(1.0d0+alphax*dt) + (1.0d0+alphax*dt/2.0d0)*ax(i)*dt

y(i) = y(i) + dt*vy(i) + ay(i)*dt*dt/2.0d0
y(i) = y(i) - (int(y(i)/Ly))*2.0d0*Ly                            ! Periodic Boundary Condition in y...
uy(i) = uy(i)*(1.0d0+alphay*dt) + (1.0d0+alphay*dt/2.0d0)*ay(i)*dt

enddo

end subroutine Gaussian_Thermostat

!=====================================================================







!       M I C R O - C A N O N I C A L   R U N  







!======================================================================

subroutine Micro_Canonical (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure)
implicit none
real*8, intent (in) :: Temp,K,t,dt,Lx,Ly
integer*8, intent (in) :: N
real*8, intent (out) :: PE,KE,Pressure,ConfigT
real*8, dimension (10000), intent (inout) :: x,y,vx,vy,ux,uy
real*8 pi,P0,xdiff,ydiff,r,fx,fy,sumvx,sumvy
real*8 virialx,virialy,virial,GF,force,gradforce
real*8, dimension (10000) :: ax,ay,GFT
integer*8 i,j

sumvx = 0.0d0
sumvy = 0.0d0
PE = 0.0d0
KE = 0.0d0
force = 0.0d0
gradforce = 0.0d0
virial = 0.0d0

do i = 1,N,1

x(i) = x(i) + dt*vx(i) + ax(i)*dt*dt/2.0d0
x(i) = x(i) - (int(x(i)/Lx))*2.0d0*Lx                            ! Periodic Boundary Condition in x...
ux(i) = vx(i) + ax(i)*dt/2.0d0

y(i) = y(i) + dt*vy(i) + ay(i)*dt*dt/2.0d0
y(i) = y(i) - (int(y(i)/Ly))*2.0d0*Ly                            ! Periodic Boundary Condition in x...
uy(i) = vy(i) + ay(i)*dt/2.0d0

enddo

!Calculating the acceleration at time dt 

do i = 1,N,1
ax(i) = 0.0d0
ay(i) = 0.0d0
GFT(i) = 0.0d0
do j = 1,N,1
if (i .ne. j) then                                          ! .ne. is used, PE should be halved...
xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0d0*Lx))*2.0d0*Lx     ! Minimum Image Convension Introduced...
ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0d0*Ly))*2.0d0*Ly     ! Minimum Image Convension Introduced...
r = sqrt((xdiff)**2 + (ydiff)**2)
fx = (xdiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ax(i) = ax(i) + fx
virialx = xdiff*fx
fy = (ydiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ay(i) = ay(i) + fy
virialy = ydiff*fy
virial = virial + (virialx + virialy)/2.0d0                   ! Calculation of virial...
GF = (1.0d0 + K*r + K*K*r*r)*exp(-K*r)/r**3
GFT(i) = GFT(i) + GF
PE = PE + (exp(-k*r)/(2.0d0*r))                               ! Calculation of Potential Energy...
endif
enddo
force = force + ax(i)**2 + ay(i)**2
gradforce = gradforce + GFT(i)
enddo

ConfigT = force/gradforce                                   ! Calculation of Configurational Temperature...

do i = 1,N,1
sumvx = sumvx + vx(i)                                       ! Check for average x-velocity...
sumvy = sumvy + vy(i)                                       ! Check for average y-velocity...
vx(i) = ux(i) + ax(i)*dt/2.0d0
vy(i) = uy(i) + ay(i)*dt/2.0d0
KE = KE + (vx(i)**2 + vy(i)**2)/2.0d0                         ! Calculation of Kinetic Energy...
enddo

Pressure = (KE + 0.5*virial)/(4.0d0*Lx*Ly)

end subroutine Micro_Canonical

!======================================================================







!       M A G N E T I C   F I E L D   A P P L I E D  







!======================================================================

subroutine Magnetic_Field (Temp,K,N,t,dt,Lx,Ly,x,y,vx,vy,ux,uy,PE,KE,ConfigT,Pressure,B0)
implicit none
real*8, intent (in) :: Temp,K,t,dt,Lx,Ly,B0
integer*8, intent (in) :: N
real*8, intent (out) :: PE,KE,Pressure,ConfigT
real*8, dimension (10000), intent (inout) :: x,y,vx,vy,ux,uy
real*8 pi,P0,xdiff,ydiff,r,fx,fy,sumvx,sumvy
real*8 virialx,virialy,virial,GF,force,gradforce
real*8, dimension (10000) :: ax,ay,GFT
integer*8 i,j

sumvx = 0.0d0
sumvy = 0.0d0
PE = 0.0d0
KE = 0.0d0
force = 0.0d0
gradforce = 0.0d0
virial = 0.0d0

do i = 1,N,1

x(i) = x(i) + dt*vx(i) + (ax(i) + vy(i)*B0)*dt*dt/2.0d0
x(i) = x(i) - (int(x(i)/Lx))*2.0d0*Lx                            ! Periodic Boundary Condition in x...
ux(i) = vx(i) + (ax(i) + vy(i)*B0)*dt/2.0d0

y(i) = y(i) + dt*vy(i) + (ay(i) - vx(i)*B0)*dt*dt/2.0d0
y(i) = y(i) - (int(y(i)/Ly))*2.0d0*Ly                            ! Periodic Boundary Condition in x...
uy(i) = vy(i) + (ay(i) - vx(i)*B0)*dt/2.0d0

enddo 

!Calculating the acceleration at time dt 

do i = 1,N,1
ax(i) = 0.0d0
ay(i) = 0.0d0
GFT(i) = 0.0d0
do j = 1,N,1
if (i .ne. j) then                                          ! .ne. is used, PE should be halved...
xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0d0*Lx))*2.0d0*Lx     ! Minimum Image Convension Introduced...
ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0d0*Ly))*2.0d0*Ly     ! Minimum Image Convension Introduced...
r = sqrt((xdiff)**2 + (ydiff)**2)
fx = (xdiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ax(i) = ax(i) + fx
virialx = xdiff*fx
fy = (ydiff)*(1.0d0+k*r)*exp(-k*r)/r**3
ay(i) = ay(i) + fy
virialy = ydiff*fy
virial = virial + (virialx + virialy)/2.0d0                   ! Calculation of virial...
GF = (1.0d0 + K*r + K*K*r*r)*exp(-K*r)/r**3
GFT(i) = GFT(i) + GF
PE = PE + (exp(-k*r)/(2.0d0*r))                               ! Calculation of Potential Energy...
endif
enddo
force = force + ax(i)**2 + ay(i)**2
gradforce = gradforce + GFT(i)
enddo

ConfigT = force/gradforce                                   ! Calculation of Configurational Temperature...

do i = 1,N,1
sumvx = sumvx + vx(i)                                       ! Check for average x-velocity...
sumvy = sumvy + vy(i)                                       ! Check for average y-velocity...
vx(i) = ux(i) + ax(i)*dt/2.0d0
vy(i) = uy(i) + ay(i)*dt/2.0d0
KE = KE + (vx(i)**2 + vy(i)**2)/2.0d0                         ! Calculation of Kinetic Energy...
enddo

Pressure = (KE + 0.5*virial)/(4.0d0*Lx*Ly)

end subroutine Magnetic_Field

!======================================================================







!       G K - E S   D I F F    C O E F F   &   V A C F 







!======================================================================

subroutine Diffusion_VACF (N,t,dt,x,y,vx,vy,Gaussian,time,DGK,DES,VACF)
implicit none
real*8, intent (in) :: Gaussian,t,dt
real*8, dimension(10000), intent (in) :: x,y,vx,vy
real*8, intent (out) :: DGK,DES,VACF
integer*8, intent (in) :: N
integer*8, intent (inout) :: time
real*8 VACFnum,VACFden
real*8, dimension (10000) :: x0,x1,y0,y1,vx0,vx1,vy0,vy1

DGK = 0.0d0
DES = 0.0d0
VACFnum = 0.0d0
VACFden = 0.0d0

time = time + 1

do i = 1,N,1
if (t .ge. Gaussian .and. t .le. Gaussian+dt) then
x0(i) = x(i)
x1(i) = x(i)
y0(i) = y(i)
y1(i) = y(i)
vx0(i) = vx(i)
vx1(i) = vx(i)
vy0(i) = vy(i)
vy1(i) = vy(i)
DGK = DGK + vx(i)*vx1(i) + vy(i)*vy1(i)
DES = DES + sqrt( (x(i)-x0(i))**2.0d0 + (y(i)-y0(i))**2.0d0 )
VACFnum = VACFnum + vx0(i)*vx(i) + vy0(i)*vy(i)
VACFden = VACFden + vx0(i)*vx0(i) + vy0(i)*vy0(i)
VACF = VACFnum/VACFden
!elseif (time == 1000) then
!x1(i) = x(i)
!y1(t,i) = y(i)
!vx1(i) = vx(i)
!vy1(i) = vy(i)
!time = 0
else 
DGK = DGK + vx(i)*vx1(i) + vy(i)*vy1(i)
DES = DES + sqrt( (x(i)-x0(i))**2.0d0 + (y(i)-y0(i))**2.0d0 )
VACFnum = VACFnum + vx0(i)*vx(i) + vy0(i)*vy(i)
VACFden = VACFden + vx0(i)*vx0(i) + vy0(i)*vy0(i)
endif
enddo

VACF = VACFnum/VACFden

end subroutine Diffusion_VACF


end
