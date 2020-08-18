! Evaluate the Multiplications in Real Space.

!===================================================================================
!================== SUBROUTINE Vel - Force - Mag ===================================
!===================================================================================

subroutine Vel_Force_Mag(Nx,Ny,Nz,Nh,pi,dx,dy,dz,A,B,C,kf,rho,ux,uy,uz,ux_dum,uy_dum,uz_dum,&
                         rho_ux,rho_uy,rho_uz,Fx,Fy,Fz,E,Bx,By,Bz,B2)

use omp_lib
implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,dx,dy,dz,A,B,C,kf

real ( kind = 8 ) x(Nx),y(Ny),z(Nz),rho(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz),E(Nx,Ny,Nz)
real ( kind = 8 ) ux_dum(Nx,Ny,Nz),uy_dum(Nx,Ny,Nz),uz_dum(Nx,Ny,Nz)
real ( kind = 8 ) Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz)
real ( kind = 8 ) rho_ux(Nx,Ny,Nz),rho_uy(Nx,Ny,Nz),rho_uz(Nx,Ny,Nz)
real ( kind = 8 ) Fx(Nx,Ny,Nz),Fy(Nx,Ny,Nz),Fz(Nx,Ny,Nz)

!if (mod(dfloat(t),0.0010d0/dt)==0.0d0) then 
 A = 0.01d0
 B = 0.01d0
 C = 0.01d0
!else
! A = 0.0d0
! B = 0.0d0
! C = 0.0d0
!endif

!$OMP PARALLEL SHARED(dx,dy,dz,A,B,C,kf,x,y,z,rho,ux,uy,uz,ux_dum,uy_dum,uz_dum),&
!$OMP & SHARED(rho_ux,rho_uy,rho_uz,Fx,Fy,Fz,E,Bx,By,Bz,B2) PRIVATE(i,j,k)
!$OMP DO
 
do i = 1,Nx
x(i)=0.0d0+real(i-1)*dx
  do j = 1,Ny
  y(j)=0.0d0+real(j-1)*dy
    do k = 1,Nz
    z(k)=0.0d0+real(k-1)*dz 
    ! FFTW Normalisation
    rho(i,j,k) = rho(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
    rho_ux(i,j,k) = rho_ux(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    rho_uy(i,j,k) = rho_uy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz)) 
    rho_uz(i,j,k) = rho_uz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
    E(i,j,k) = E(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
    Bx(i,j,k) = Bx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    By(i,j,k) = By(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    Bz(i,j,k) = Bz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
    ! Evaluate Velocity in Real Space. 
    ux(i,j,k) = rho_ux(i,j,k)/rho(i,j,k)
    uy(i,j,k) = rho_uy(i,j,k)/rho(i,j,k)
    uz(i,j,k) = rho_uz(i,j,k)/rho(i,j,k)

    ! Evaluate Forcing.
    Fx(i,j,k) = A*dsin(kf*z(k)) + C*dcos(kf*y(j))
    Fy(i,j,k) = B*dsin(kf*x(i)) + A*dcos(kf*z(k))
    Fz(i,j,k) = C*dsin(kf*y(j)) + B*dcos(kf*x(i))
  
    ! Evaluate Square of Magnetic Field.
    B2(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)
  
    ! Keep Backup of the Arrays for FFTW.
    ux_dum(i,j,k) = ux(i,j,k)
    uy_dum(i,j,k) = uy(i,j,k)
    uz_dum(i,j,k) = uz(i,j,k)
    enddo
  enddo
enddo   

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine Vel_Force_Mag

!====================================================================================

