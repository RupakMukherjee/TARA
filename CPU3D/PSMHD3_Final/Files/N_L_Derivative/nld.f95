! Evaluate the Non Linear Derivatives in Spectral Space.

!===================================================================================
!================== SUBROUTINE Non Linear Derivative ===============================
!===================================================================================

subroutine Non_Linear_Derivative(Nx,Ny,Nz,Nh,pi,mu,eta,i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k,&
                                 i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,&
                                 i_kz_Mom_y_3_k,i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k,&
                                 i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,E_Visc_k,&
                                 i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,&
                                 i_kx_Mag_z_1_k,i_ky_Mag_z_2_k,&
                                 kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,&
                                 ky2_uz_k,kz2_uz_k,nu,Fx_k,Fy_k,Fz_k,&
                                 kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k,&
                                 d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new,&
                                 d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new)

use omp_lib
implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,mu,eta

complex ( kind = 8 ) nu(Nh,Ny,Nz),Fx_k(Nh,Ny,Nz),Fy_k(Nh,Ny,Nz),Fz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_rho_ux_k(Nh,Ny,Nz),i_ky_rho_uy_k(Nh,Ny,Nz),i_kz_rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_x_1_k(Nh,Ny,Nz),i_ky_Mom_x_2_k(Nh,Ny,Nz),i_kz_Mom_x_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_y_1_k(Nh,Ny,Nz),i_ky_Mom_y_2_k(Nh,Ny,Nz),i_kz_Mom_y_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_z_1_k(Nh,Ny,Nz),i_ky_Mom_z_2_k(Nh,Ny,Nz),i_kz_Mom_z_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Energy_x_k(Nh,Ny,Nz),i_ky_Energy_y_k(Nh,Ny,Nz),i_kz_Energy_z_k(Nh,Ny,Nz),E_Visc_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_ky_Mag_x_1_k(Nh,Ny,Nz),i_kz_Mag_x_2_k(Nh,Ny,Nz),i_kx_Mag_y_1_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kz_Mag_y_2_k(Nh,Ny,Nz),i_kx_Mag_z_1_k(Nh,Ny,Nz),i_ky_Mag_z_2_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_ky_Bx_k(Nh,Ny,Nz),i_kx_By_k(Nh,Ny,Nz),i_kx_Bz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kz_Bx_k(Nh,Ny,Nz),i_ky_Bz_k(Nh,Ny,Nz),i_kz_By_k(Nh,Ny,Nz)

complex ( kind = 8 ) kx2_ux_k(Nh,Ny,Nz),ky2_ux_k(Nh,Ny,Nz),kz2_ux_k(Nh,Ny,Nz),kx2_uy_k(Nh,Ny,Nz),ky2_uy_k(Nh,Ny,Nz)
complex ( kind = 8 ) kz2_uy_k(Nh,Ny,Nz),kx2_uz_k(Nh,Ny,Nz),ky2_uz_k(Nh,Ny,Nz),kz2_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) kx2_Bx_k(Nh,Ny,Nz),ky2_Bx_k(Nh,Ny,Nz),kz2_Bx_k(Nh,Ny,Nz),kx2_By_k(Nh,Ny,Nz),ky2_By_k(Nh,Ny,Nz)
complex ( kind = 8 ) kz2_By_k(Nh,Ny,Nz),kx2_Bz_k(Nh,Ny,Nz),ky2_Bz_k(Nh,Ny,Nz),kz2_Bz_k(Nh,Ny,Nz)

complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny,Nz),d_rho_ux_k_dt_old(Nh,Ny,Nz),d_rho_uy_k_dt_old(Nh,Ny,Nz),d_rho_uz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_old(Nh,Ny,Nz),d_Bx_k_dt_old(Nh,Ny,Nz),d_By_k_dt_old(Nh,Ny,Nz),d_Bz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny,Nz),d_rho_ux_k_dt_new(Nh,Ny,Nz),d_rho_uy_k_dt_new(Nh,Ny,Nz),d_rho_uz_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_new(Nh,Ny,Nz),d_Bx_k_dt_new(Nh,Ny,Nz),d_By_k_dt_new(Nh,Ny,Nz),d_Bz_k_dt_new(Nh,Ny,Nz)


!$OMP PARALLEL SHARED(mu,eta,i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k),&
!$OMP & SHARED(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k),&
!$OMP & SHARED(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k),&
!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,E_Visc_k)&
!$OMP & SHARED(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k),&
!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k,nu,Fx_k,Fy_k,Fz_k),&
!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k),&
!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new),&
!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new) PRIVATE(i,j,k)
!$OMP DO 

do i = 1,Nx/2+1
  do j = 1,Ny
    do k = 1,Nz
      ! Density Equation.
      d_rho_k_dt_new(i,j,k) = - ( i_kx_rho_ux_k(i,j,k) + i_ky_rho_uy_k(i,j,k) + i_kz_rho_uz_k(i,j,k) )
    
      ! Momentum Equation.
      d_rho_ux_k_dt_new(i,j,k) = - ( i_kx_Mom_x_1_k(i,j,k) + i_ky_Mom_x_2_k(i,j,k) + i_kz_Mom_x_3_k(i,j,k) ) &
                                 - mu * ( kx2_ux_k(i,j,k) + ky2_ux_k(i,j,k) + kz2_ux_k(i,j,k) ) + Fx_k(i,j,k) 
    
      d_rho_uy_k_dt_new(i,j,k) = - ( i_kx_Mom_y_1_k(i,j,k) + i_ky_Mom_y_2_k(i,j,k) + i_kz_Mom_y_3_k(i,j,k) ) &
                                 - mu * ( kx2_uy_k(i,j,k) + ky2_uy_k(i,j,k) + kz2_uy_k(i,j,k) ) + Fy_k(i,j,k)
    
      d_rho_uz_k_dt_new(i,j,k) = - ( i_kx_Mom_z_1_k(i,j,k) + i_ky_Mom_z_2_k(i,j,k) + i_kz_Mom_z_3_k(i,j,k) ) &
                                 - mu * ( kx2_uz_k(i,j,k) + ky2_uz_k(i,j,k) + kz2_uz_k(i,j,k) ) + Fz_k(i,j,k) 
    
      ! Energy Equation.
      d_E_k_dt_new(i,j,k) = - ( i_kx_Energy_x_k(i,j,k) + i_ky_Energy_y_k(i,j,k) + i_kz_Energy_z_k(i,j,k) ) &
                            + mu * E_Visc_k(i,j,k)
    
      ! Magnetic Field Equation.
      d_Bx_k_dt_new(i,j,k) = + ( i_ky_Mag_x_1_k(i,j,k) + i_kz_Mag_x_2_k(i,j,k) ) &
                             - eta * ( kx2_Bx_k(i,j,k) + ky2_Bx_k(i,j,k) + kz2_Bx_k(i,j,k) )
                           
      d_By_k_dt_new(i,j,k) = - ( i_kx_Mag_y_1_k(i,j,k) - i_kz_Mag_y_2_k(i,j,k) ) &
                             - eta * ( kx2_By_k(i,j,k) + ky2_By_k(i,j,k) + kz2_By_k(i,j,k) ) 
                           
      d_Bz_k_dt_new(i,j,k) = - ( i_kx_Mag_z_1_k(i,j,k) + i_ky_Mag_z_2_k(i,j,k) ) &                       
                             - eta * ( kx2_Bz_k(i,j,k) + ky2_Bz_k(i,j,k) + kz2_Bz_k(i,j,k) ) 
    enddo                       
  enddo
enddo  

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine Non_Linear_Derivative

!====================================================================================

