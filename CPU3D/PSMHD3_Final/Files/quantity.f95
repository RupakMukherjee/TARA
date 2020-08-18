! Evaluate Divergence of Magnetic Field.

!===================================================================================
!================== SUBROUTINE DIVERGENCE_B ========================================
!===================================================================================

subroutine E_Casimir_Div_B (Nx,Ny,Nz,Nh,pi,mu_0,omega_x,omega_y,omega_z,omega2,rho,ux,uy,uz,P,E,Bx,By,Bz,div_B,&
                            Pressure,Energy,y_Energy,B_Field,C0,C1,C2,C3,div_B_Tot)

use omp_lib
implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,mu_0
real ( kind = 8 ) rho(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
real ( kind = 8 ) P(Nx,Ny,Nz),E(Nx,Ny,Nz),Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz)
real ( kind = 8 ) omega_x(Nx,Ny,Nz),omega_y(Nx,Ny,Nz),omega_z(Nx,Ny,Nz),omega2(Nx,Ny,Nz),div_B(Nx,Ny,Nz)
real ( kind = 8 ) Pressure,Energy,y_Energy,B_Field,C0,C1,C2,C3,div_B_tot

Pressure = 0.0d0; Energy = 0.0d0; y_Energy = 0.0d0; B_Field = 0.0d0
 C0 = 0.0d0; C1 = 0.0d0; C2 = 0.0d0; C3 = 0.0d0; div_B_Tot = 0.0d0
   
!$OMP PARALLEL SHARED(mu_0,omega_x,omega_y,omega_z,omega2,rho,ux,uy,uz,P,Bx,By,Bz,div_B) PRIVATE(i,j,k)
!$OMP DO REDUCTION (+:Pressure,Energy,y_Energy,B_Field,C0,C1,C2,C3,div_B_Tot) 

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      ! Evaluate Pressure
      Pressure = Pressure + P(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Energy
      !Energy = Energy + E(i,j,k)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Energy = Energy + (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Growth Rate.
      y_Energy = y_Energy + rho(i,j,k)*(uy(i,j,k)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Magnetic Field.
      B_Field = B_Field + (Bx(i,j,k)**2+By(i,j,k)**2+Bz(i,j,k)**2)/(2.0d0*mu_0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Casimirs.
      C0 = C0 + dsqrt(omega2(i,j,k))**0.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C1 = C1 + dsqrt(omega2(i,j,k))**1.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C2 = C2 + dsqrt(omega2(i,j,k))**2.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C3 = C3 + dsqrt(omega2(i,j,k))**3.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Check for Div B = 0
      div_B_Tot = div_B_Tot + div_B(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    enddo
  enddo
enddo  
      
!$OMP END PARALLEL
      
return

end subroutine E_Casimir_Div_B

!====================================================================================


