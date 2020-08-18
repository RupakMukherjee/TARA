! Evaluate Divergence of Magnetic Field.

!===================================================================================
!================== SUBROUTINE DIVERGENCE_B ========================================
!===================================================================================

subroutine U_B2_P (Nx,Ny,Nz,Nh,spheat,rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz,B2,&
                   ux,uy,uz,ux_dum,uy_dum,uz_dum,P)

use omp_lib
implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) spheat,rho(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
real ( kind = 8 ) rho_ux(Nx,Ny,Nz),rho_uy(Nx,Ny,Nz),rho_uz(Nx,Ny,Nz)
real ( kind = 8 ) ux_dum(Nx,Ny,Nz),uy_dum(Nx,Ny,Nz),uz_dum(Nx,Ny,Nz)
real ( kind = 8 ) P(Nx,Ny,Nz),E(Nx,Ny,Nz),Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz)

!$OMP PARALLEL SHARED(spheat,rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz,B2), &
!$OMP & SHARED(ux,uy,uz,ux_dum,uy_dum,uz_dum,P) PRIVATE(i,j,k)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      ! FFTW Normalisation.
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
  
      ! Evaluate Square of Magnetic Field.
      B2(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)
  
      ! Keep Backup of the Arrays for FFTW.
      ux_dum(i,j,k) = ux(i,j,k)
      uy_dum(i,j,k) = uy(i,j,k)
      uz_dum(i,j,k) = uz(i,j,k)
  
      ! Evaluate Pressure
      P(i,j,k) = ( spheat - 1.0d0 ) * ( E(i,j,k) - 0.50d0 * &
                 ( rho_ux(i,j,k)*ux(i,j,k)+rho_uy(i,j,k)*uy(i,j,k)+rho_uz(i,j,k)*uz(i,j,k) ) )! - B2(i,j) ) )
    enddo
  enddo
enddo 
  
!$OMP END DO
!$OMP END PARALLEL 
      
return

end subroutine U_B2_P

!====================================================================================

