program Configurational_Temperature
implicit none
integer i,j,N,MCS,MCSmin,MCSmax,MCSdiff,G
real pi,nbar,Lx,Ly,Lz,L,xdiff,ydiff,zdiff,r,rCut,K,Gamma,Temp,Fx,Fy,Fz,GF,force,gradforce,grandtemp
real, dimension(686) :: q,x,y,z,FxT,FyT,FzT,GFT
pi = 4.0*atan(1.0)
!pi = 3.141592653589793238460d0


!====================== User Inputs =============================

!Coupling Parameter...
Gamma = 300.0 ! MAKE YOUR GAMMA LOW (RUGH & BUTLER)...
Temp = 1.0/Gamma

!Screening Parameter...
K = 1.0

N = 686

Lx = 7.105
Ly = 7.105
Lz = 7.105
L = Lx

nbar = float(N)/(2.0*Lx*2.0*Ly*2.0*Lz)

!Charge of the particles...
do i = 1,N,1
q(i) = 1.0
enddo

!Declaration of the Total number of the Monte Carlo Steps...
MCSmin = 80000
MCSmax = 100000
MCSdiff = 100

!Cutoff of r...
rCut = 20.0 !L/2.0


!====================== Output Filenames ===========================

open(unit=10,file='Temp.dat',status='unknown')

!====================== Energy & Config T calculation ===============

do MCS = MCSmin, MCSmax, MCSdiff


force = 0.0
gradforce = 0.0


   do i = 1, N, 1
   read(MCS,*) G,x(i),y(i),z(i)
   enddo ! i


   do i = 1,N,1

   FxT(i) = 0.0
   FyT(i) = 0.0
   FzT(i) = 0.0
   GFT(i) = 0.0

	! Calculation of the Potential Energy & Config T...
	
	do j = 1,N,1
 
       ! Minimum Image Convension...
        xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0*Lx))*2.0*Lx
        ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0*Ly))*2.0*Ly
        zdiff = (z(i)-z(j)) - nint((z(i)-z(j))/(2.0*Lz))*2.0*Lz
        
        r = sqrt(xdiff**2 + ydiff**2 + zdiff**2) 


        if (i .ne. j .and. r*r .lt. rCut*rCut)then           ! .ne. is used, Energy should be halved...
	    Fx = (q(i)*q(j))*(xdiff)*(1.0 + K*r)*exp(-K*r)/r**3
	    Fy = (q(i)*q(j))*(ydiff)*(1.0 + K*r)*exp(-K*r)/r**3
	    Fz = (q(i)*q(j))*(zdiff)*(1.0 + K*r)*exp(-K*r)/r**3
	    GF = (q(i)*q(j))*(1.0 + K*r + K*K*r*r)*exp(-K*r)/r**3
        else
	    Fx = 0.0
	    Fy = 0.0
	    Fz = 0.0
	    GF = 0.0
        endif
   
	    FxT(i) = FxT(i) + Fx
	    FyT(i) = FyT(i) + Fy
	    FzT(i) = FzT(i) + Fz
	    GFT(i) = GFT(i) + GF
	enddo  ! j

   force = force + (FxT(i)*FxT(i) + FyT(i)*FyT(i) + FzT(i)*FzT(i))
   gradforce = gradforce + GFT(i)

   enddo ! i

write(10,*) MCS, 5.0*force/(3.0*gradforce), 3.0*gradforce/(5.0*force), Gamma

grandtemp = grandtemp + force/gradforce

enddo ! MCS

write(10,*) grandtemp/((float(MCSmax) - float(MCSmin))/float(MCSdiff) + 1.0), Temp
 
 close(10)

end program Configurational_Temperature

