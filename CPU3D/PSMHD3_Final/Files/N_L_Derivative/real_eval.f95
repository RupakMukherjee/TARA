! Evaluate the Multiplications in Real Space.

!===================================================================================
!================== SUBROUTINE Real - Evaluation ===================================
!===================================================================================

subroutine Real_Evaluation(Nx,Ny,Nz,Nh,pi,spheat,eta,ux,uy,uz,rho_ux,rho_uy,rho_uz,Bx,By,Bz,&
                           d_ux_dx,d_uy_dy,d_uz_dz,P,E,&
                           d_Bx_dy,d_Bx_dz,d_By_dx,d_By_dz,d_Bz_dx,d_Bz_dy,curl_x_B,curl_y_B,curl_z_B,B2,&
                           Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3,&
                           Energy_x,Energy_y,Energy_z,E_Visc,&
                           Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2)

use omp_lib
implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,eta,spheat
real ( kind = 8 ) ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz),E(Nx,Ny,Nz),P(Nx,Ny,Nz)
real ( kind = 8 ) rho_ux(Nx,Ny,Nz),rho_uy(Nx,Ny,Nz),rho_uz(Nx,Ny,Nz)
real ( kind = 8 ) Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz)

real ( kind = 8 ) d_ux_dx(Nx,Ny,Nz),d_uy_dy(Nx,Ny,Nz),d_uz_dz(Nx,Ny,Nz)
real ( kind = 8 ) d_Bx_dy(Nx,Ny,Nz),d_By_dx(Nx,Ny,Nz),d_Bx_dz(Nx,Ny,Nz)
real ( kind = 8 ) d_By_dz(Nx,Ny,Nz),d_Bz_dx(Nx,Ny,Nz),d_Bz_dy(Nx,Ny,Nz)
real ( kind = 8 ) curl_x_B(Nx,Ny,Nz),curl_y_B(Nx,Ny,Nz),curl_z_B(Nx,Ny,Nz)

real ( kind = 8 ) Mom_x_1(Nx,Ny,Nz),Mom_x_2(Nx,Ny,Nz),Mom_x_3(Nx,Ny,Nz)
real ( kind = 8 ) Mom_y_1(Nx,Ny,Nz),Mom_y_2(Nx,Ny,Nz),Mom_y_3(Nx,Ny,Nz)
real ( kind = 8 ) Mom_z_1(Nx,Ny,Nz),Mom_z_2(Nx,Ny,Nz),Mom_z_3(Nx,Ny,Nz)
real ( kind = 8 ) Energy_x(Nx,Ny,Nz),Energy_y(Nx,Ny,Nz),Energy_z(Nx,Ny,Nz),E_Visc(Nx,Ny,Nz)
real ( kind = 8 ) Mag_x_1(Nx,Ny,Nz),Mag_x_2(Nx,Ny,Nz),Mag_y_1(Nx,Ny,Nz)
real ( kind = 8 ) Mag_y_2(Nx,Ny,Nz),Mag_z_1(Nx,Ny,Nz),Mag_z_2(Nx,Ny,Nz)

!$OMP PARALLEL SHARED(spheat,eta,ux,uy,uz,rho_ux,rho_uy,rho_uz,Bx,By,Bz),&
!$OMP & SHARED(d_ux_dx,d_uy_dy,d_uz_dz,P,E),&
!$OMP & SHARED(d_Bx_dy,d_Bx_dz,d_By_dx,d_By_dz,d_Bz_dx,d_Bz_dy,curl_x_B,curl_y_B,curl_z_B,B2),&
!$OMP & SHARED(Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3),&
!$OMP & SHARED(Energy_x,Energy_y,Energy_z,E_Visc),&
!$OMP & SHARED(Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2) PRIVATE(i,j,k)
!$OMP DO
 
do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
    ! FFTW Normalisation.
    d_ux_dx(i,j,k) = d_ux_dx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_uy_dy(i,j,k) = d_uy_dy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_uz_dz(i,j,k) = d_uz_dz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
   
    d_Bx_dy(i,j,k) = d_Bx_dy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_Bx_dz(i,j,k) = d_Bx_dz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_By_dx(i,j,k) = d_By_dx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_By_dz(i,j,k) = d_By_dz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_Bz_dx(i,j,k) = d_Bz_dx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_Bz_dy(i,j,k) = d_Bz_dy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
  
    ! Evaluate Curl of Magnetic Field.
    curl_x_B(i,j,k) = d_Bz_dy(i,j,k) - d_By_dz(i,j,k)
    curl_y_B(i,j,k) = d_Bz_dx(i,j,k) - d_Bx_dz(i,j,k)
    curl_z_B(i,j,k) = d_By_dx(i,j,k) - d_Bx_dy(i,j,k)
    
    ! Evaluate Pressure
    P(i,j,k) = ( spheat - 1.0d0 ) * ( E(i,j,k) &
               - 0.50d0 * ( rho_ux(i,j,k)*ux(i,j,k)+rho_uy(i,j,k)*uy(i,j,k)+rho_uz(i,j,k)*uz(i,j,k) )) !&
               !- B2(i,j,k) ) )
  
    ! Evaluate LHS of Momentum Equation.
    Mom_x_1(i,j,k) = rho_ux(i,j,k)*ux(i,j,k) !+ P(i,j,k) + B2(i,j,k)/2.0d0 - Bx(i,j,k)*Bx(i,j,k)
    Mom_x_2(i,j,k) = rho_ux(i,j,k)*uy(i,j,k) !- Bx(i,j,k)*By(i,j,k)
    Mom_x_3(i,j,k) = rho_ux(i,j,k)*uz(i,j,k) !- Bx(i,j,k)*Bz(i,j,k)
  
    Mom_y_1(i,j,k) = rho_ux(i,j,k)*uy(i,j,k) !- Bx(i,j,k)*By(i,j,k) 
    Mom_y_2(i,j,k) = rho_uy(i,j,k)*uy(i,j,k) !+ P(i,j,k) + B2(i,j,k)/2.0d0 - By(i,j,k)*By(i,j,k)
    Mom_y_3(i,j,k) = rho_uy(i,j,k)*uz(i,j,k) !- By(i,j,k)*Bz(i,j,k)
  
    Mom_z_1(i,j,k) = rho_uz(i,j,k)*ux(i,j,k) !- Bz(i,j,k)*Bx(i,j,k) 
    Mom_z_2(i,j,k) = rho_uz(i,j,k)*uy(i,j,k) !- Bz(i,j,k)*By(i,j,k) 
    Mom_z_3(i,j,k) = rho_uz(i,j,k)*uz(i,j,k) !+ P(i,j,k) + B2(i,j,k)/2.0d0 - Bz(i,j,k)*Bz(i,j,k)
  
    ! Evaluate LHS of Energy Equation.
    Energy_x(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * ux(i,j,k) &
                      - ux(i,j,k)*Bx(i,j,k)*Bx(i,j,k) - uy(i,j,k)*Bx(i,j,k)*By(i,j,k) - uz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) &
                      - eta * ( By(i,j,k) * curl_z_B(i,j,k) + Bz(i,j,k) * curl_y_B(i,j,k)) 
                      
    Energy_y(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * uy(i,j,k) &
                      - ux(i,j,k)*Bx(i,j,k)*By(i,j,k) - uy(i,j,k)*By(i,j,k)*By(i,j,k) - uz(i,j,k)*By(i,j,k)*Bz(i,j,k) &
                      + eta * ( Bx(i,j,k) * curl_z_B(i,j,k) - Bz(i,j,k) * curl_x_B(i,j,k))
                      
    Energy_z(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * uz(i,j,k) &
                      - ux(i,j,k)*Bz(i,j,k)*Bx(i,j,k) - uy(i,j,k)*Bz(i,j,k)*By(i,j,k) - uz(i,j,k)*Bz(i,j,k)*Bz(i,j,k) & 
                      + eta * ( Bx(i,j,k) * curl_y_B(i,j,k) + By(i,j,k) * curl_x_B(i,j,k))
                      
    ! Evaluate RHS of Energy Equation.
    E_Visc(i,j,k) = ( d_ux_dx(i,j,k) + d_uy_dy(i,j,k) + d_uz_dz(i,j,k) )**2
  
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

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine Real_Evaluation

!====================================================================================

