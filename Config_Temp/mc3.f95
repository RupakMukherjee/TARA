program phase_transition3d
implicit none
integer i,j,N,N1,N2,MCS,MCSmax,count
real pi,nbar,Lx,Ly,Lz,L,E,E1,dx,dy,dz,xdiff,ydiff,zdiff,r,rCut,R1,R2,ra,P,K,Gamma,Temp,Delta,Energy,Sum
real, dimension(100000) :: q,x,y,z,xmov,ymov,zmov,xdum,ydum,zdum,W,Wdum
integer,parameter :: seed = 99999999
pi = 4.0*atan(1.0)
 call srand(seed)

!====================== User Inputs =============================

!System Size in dimensionless quantities...
Lx = 7.105
Ly = 7.105
Lz = 7.105
L = Lx

!Coupling Parameter...
Gamma = 100.0
Temp = 1.0/Gamma

!Screening Parameter...
K = 1.0

N = 5000
nbar = float(N)/(2.0*Lx*2.0*Ly*2.0*Lz)
!nbar = 0.206 !1/pi
!L = (float(N)/(6.0*nbar))**(1.0/3.0)
!Lx = L
!Ly = L
!Lz = L

!Radius of the two types of particles...
R1 = 0.050
R2 = 0.050

!Number of particles having hard cores...
N1 = int(float(N)/2.0)
N2 = N-N1

!Charge of the particles...
do i = 1,N,1
q(i) = 1.0
enddo

!Declaration of the Total number of the Monte Carlo Steps...
MCSmax=150000

!Step width of a single particle for the operation of MC...
Delta = 0.068

!Cutoff of r...
rCut = 20.0 !L/2.0


Energy = 0.0

Sum = 0.0

!====================== Output Filenames ===========================

open(unit=10,file='Energy.dat',status='unknown')
open(unit=20,file='Acceptance_Ratio.dat',status='unknown')
open(unit=30,file='cross_check.dat',status='unknown')
open(unit=40,file='Initial_Configuration.dat',status='unknown')

!====================== Definition of initial state =================

   ! Definition of the initial random positions in -Lx to Lx and -Ly to Ly rectangular box...   
   do i = 1, N, 1

   x(i) = (rand())*2.0*Lx - Lx
   y(i) = (rand())*2.0*Ly - Ly
   z(i) = (rand())*2.0*Lz - Lz

    !read(40,*) MCS,x(i),y(i),z(i)

   enddo ! i

   ! Calculation of the Initial Potential Energy of the system...
   do i = 1, N, 1
   W(i) = 0.0
   do j = 1, N, 1

   ! Minimum Image Convension...
   xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0*Lx))*2.0*Lx     
!   if (xdiff >   Lx) xdiff = xdiff - 2.0*Lx
!   if (xdiff <= -Lx) xdiff = xdiff + 2.0*Lx

   ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0*Ly))*2.0*Ly
!   if (ydiff >   Ly) ydiff = ydiff - 2.0*Ly
!   if (ydiff <= -Ly) ydiff = ydiff + 2.0*Ly

   zdiff = (z(i)-z(j)) - nint((z(i)-z(j))/(2.0*Lz))*2.0*Lz
!   if (zdiff >   Lz) zdiff = zdiff - 2.0*Lz
!   if (zdiff <= -Lz) zdiff = zdiff + 2.0*Lz

   r = sqrt(xdiff**2 + ydiff**2 + zdiff**2) 


   if (i .ne. j .and. r*r .lt. rCut*rCut) then                 ! .ne. is used, Energy should be halved...
   E1 = (q(i)*q(j)/r)*exp(-K*r)
!   if (i .le. N1 .and. j .le. N1)then                                  
!   E1 = (q(i)*q(j)/(r+R1))*exp(-K*(r+R1))
!   elseif (i .gt. N1 .and. j .gt. N1)then
!   E1 = (q(i)*q(j)/(r+R2))*exp(-K*(r+R2))
!   elseif (i .gt. N1 .and. j .le. N1)then
!   E1 = (q(i)*q(j)/(r+(R1+R2)/2.0))*exp(-K*(r+(R1+R2)/2.0))
!   elseif (i .le. N1 .and. j .gt. N1)then
!   E1 = (q(i)*q(j)/(r+(R1+R2)/2.0))*exp(-K*(r+(R1+R2)/2.0))
!   else
!   E1 = 0.0
!   endif
   else
   E1 = 0.0
   endif
   
   W(i) = W(i) + E1

   enddo  ! j

   Energy = Energy + W(i)

   enddo ! i
   
   write(10,*) 0, Energy/(2.0*float(N))


!================ Monte Carlo started ===========================

do MCS = 1, MCSmax, 1

Energy = 0.0
count = 0

   do i = 1,N,1

   W(i) = 0.0

	!Recalculating the Energy for comparison with a Random Number...
	do j = 1,N,1
  
       ! Minimum Image Convension...
        xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0*Lx))*2.0*Lx
!        if (xdiff >   Lx) xdiff = xdiff - 2.0*Lx
!        if (xdiff <= -Lx) xdiff = xdiff + 2.0*Lx

        ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0*Ly))*2.0*Ly
!        if (ydiff >   Ly) ydiff = ydiff - 2.0*Ly
!        if (ydiff <= -Ly) ydiff = ydiff + 2.0*Ly

        zdiff = (z(i)-z(j)) - nint((z(i)-z(j))/(2.0*Lz))*2.0*Lz
!        if (zdiff >   Lz) zdiff = zdiff - 2.0*Lz
!        if (zdiff <= -Lz) zdiff = zdiff + 2.0*Lz

        r = sqrt(xdiff**2 + ydiff**2 + zdiff**2) 


        if (i .ne. j .and. r*r .lt. rCut*rCut)then               ! .ne. is used, Energy should be halved...
        E1 = (q(i)*q(j)/r)*exp(-K*r)
!        if (i .le. N1 .and. j .le. N1)then                               
!        E1 = (q(i)*q(j)/(r+R1))*exp(-K*(r+R1))
!        elseif (i .gt. N1 .and. j .gt. N1)then
!        E1 = (q(i)*q(j)/(r+R2))*exp(-K*(r+R2))
!        elseif (i .gt. N1 .and. j .le. N1)then
!        E1 = (q(i)*q(j)/(r+(R1+R2)/2.0))*exp(-K*(r+(R1+R2)/2.0))
!        elseif (i .le. N1 .and. j .gt. N1)then
!        E1 = (q(i)*q(j)/(r+(R1+R2)/2.0))*exp(-K*(r+(R1+R2)/2.0))
!        else
!        E1 = 0.0
!        endif
        else
	E1 = 0.0
        endif
   
        W(i) = W(i) + E1

	enddo  ! j

  
     ! Step Length of one displacement for a single particle...
     dx = rand()*Lx*Delta - Lx*Delta/2.0                      ! dx can be positive or negative...
     dy = rand()*Ly*Delta - Ly*Delta/2.0                      ! dy can be positive or negative...
     dz = rand()*Lz*Delta - Lz*Delta/2.0                      ! dz can be positive or negative...
     !print*, dx,dy

     xmov(i) = x(i) + dx
     xdum(i) = xmov(i) - (int(xmov(i)/Lx))*2.0*Lx                 ! Periodic Boundary Condition in x...

     ymov(i) = y(i) + dy          
     ydum(i) = ymov(i) - (int(ymov(i)/Ly))*2.0*Ly                 ! Periodic Boundary Condition in y...

     zmov(i) = z(i) + dz          
     zdum(i) = zmov(i) - (int(zmov(i)/Lz))*2.0*Lz                 ! Periodic Boundary Condition in z...
     
   Wdum(i) = 0.0

   ! Calculation of the Potential Energy for a probable configuration...
   do j = 1,N,1
 
     ! Minimum Image Convension...
     xdiff = (xdum(i)-x(j)) - nint((xdum(i)-x(j))/(2.0*Lx))*2.0*Lx
!     if (xdiff >   Lx) xdiff = xdiff - 2.0*Lx
!     if (xdiff <= -Lx) xdiff = xdiff + 2.0*Lx

     ydiff = (ydum(i)-y(j)) - nint((ydum(i)-y(j))/(2.0*Ly))*2.0*Ly
!     if (ydiff >   Ly) ydiff = ydiff - 2.0*Ly
!     if (ydiff <= -Ly) ydiff = ydiff + 2.0*Ly

     zdiff = (zdum(i)-z(j)) - nint((zdum(i)-z(j))/(2.0*Lz))*2.0*Lz
!     if (zdiff >   Lz) zdiff = zdiff - 2.0*Lz
!     if (zdiff <= -Lz) zdiff = zdiff + 2.0*Lz

     r = sqrt(xdiff**2 + ydiff**2 + zdiff**2) 


     if (i .ne. j .and. r*r .lt. rCut*rCut)then                    ! .ne. is used, Energy should be halved...
     E = (q(i)*q(j)/r)*exp(-K*r)
!     if (i .le. N1 .and. j .le. N1)then                                 
!     E = (q(i)*q(j)/(r+R1))*exp(-K*(r+R1))
!     elseif (i .gt. N1 .and. j .gt. N1)then
!     E = (q(i)*q(j)/(r+R2))*exp(-K*(r+R2))
!     elseif (i .gt. N1 .and. j .le. N1)then
!     E = (q(i)*q(j)/(r+(R1+R2)/2.0))*exp(-K*(r+(R1+R2)/2.0))
!     elseif (i .le. N1 .and. j .gt. N1)then
!     E = (q(i)*q(j)/(r+(R1+R2)/2.0))*exp(-K*(r+(R1+R2)/2.0))
!     else
!     E = 0.0
!     endif
     else
     E = 0.0
     endif
   
     Wdum(i) = Wdum(i) + E

   enddo  ! j 

     ! Calculation for the acceptibility of the probable configuration...
     P = exp(-(Wdum(i)-W(i))/Temp)
10   ra = rand()
     !print*, Wdum(i),W(i),ra,abs(Wdum(i)-W(i)),Temp

     if (ra .le. P)then                   ! Accept
     x(i) = xdum(i)
     y(i) = ydum(i)
     z(i) = zdum(i)
     count = count+1

     else                                 ! Reject
     x(i) = x(i)
     y(i) = y(i)
     z(i) = z(i)
     endif

   enddo ! i

   do i = 1,N,1

	if (MCS .ge. 80000 .and. mod(float(MCS),100.0) == 0.0) then
        write(MCS,*) MCS,x(i),y(i),z(i)
	endif

   W(i) = 0.0

	! Calculation of the Potential Energy after a single displacement (May be accepted or rejected)...
	do j = 1,N,1
 
       ! Minimum Image Convension...
        xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0*Lx))*2.0*Lx
!        if (xdiff >   Lx) xdiff = xdiff - 2.0*Lx
!        if (xdiff <= -Lx) xdiff = xdiff + 2.0*Lx

        ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0*Ly))*2.0*Ly
!        if (ydiff >   Ly) ydiff = ydiff - 2.0*Ly
!        if (ydiff <= -Ly) ydiff = ydiff + 2.0*Ly

        zdiff = (z(i)-z(j)) - nint((z(i)-z(j))/(2.0*Lz))*2.0*Lz
!        if (zdiff >   Lz) zdiff = zdiff - 2.0*Lz
!        if (zdiff <= -Lz) zdiff = zdiff + 2.0*Lz

        r = sqrt(xdiff**2 + ydiff**2 + zdiff**2) 


        if (i .ne. j .and. r*r .lt. rCut*rCut)then           ! .ne. is used, Energy should be halved...
        E1 = (q(i)*q(j)/r)*exp(-K*r)
!        if (i .le. N1 .and. j .le. N1)then                                  
!        E1 = (q(i)*q(j)/(r+R1))*exp(-K*(r+R1))
!        elseif (i .gt. N1 .and. j .gt. N1)then
!        E1 = (q(i)*q(j)/(r+R2))*exp(-K*(r+R2))
!        elseif (i .gt. N1 .and. j .le. N1)then
!        E1 = (q(i)*q(j)/(r+(R1+R2)/2.0))*exp(-K*(r+(R1+R2)/2.0))
!        elseif (i .le. N1 .and. j .gt. N1)then
!        E1 = (q(i)*q(j)/(r+(R1+R2)/2.0))*exp(-K*(r+(R1+R2)/2.0))
!        else
!        E1 = 0.0
!        endif
        else
	E1 = 0.0
        endif
   
        W(i) = W(i) + E1

	enddo  ! j

   Energy = Energy + W(i)

   enddo ! i

write(10,*) MCS, Energy/(2.0*float(N))

Sum = Sum + count/float(N)
write(20,*) MCS, count/float(N), Sum

enddo ! MCS

! Sum/tmax is the acceptance ratio...
write(30,*) "# of Particles","Density of Particles","Average Acceptance Ratio"
write(30,*) N,float(N)/(Lx*2.0*Ly*2.0),Sum/MCSmax
 
 close(10)
 close(20)
 close(30)
 close(40)

end program phase_transition3d

