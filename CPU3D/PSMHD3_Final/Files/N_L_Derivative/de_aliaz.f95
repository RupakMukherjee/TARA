! De - Aliazing Technique With 2/3 Rule for All the Non-Linear Terms.  

!===================================================================================
!================== SUBROUTINE DE - ALIAZING =======================================
!===================================================================================

subroutine DE_ALIAZING (Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k,&
                        i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k,&
                        i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k,&
                        i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,E_Visc_k,&
                        i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k,&
                        kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k,nu,&
                        kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k)

use omp_lib
implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,Lx,Ly,Lz,kx,ky,kz

complex ( kind = 8 ) i_kx_rho_ux_k(Nh,Ny,Nz),i_ky_rho_uy_k(Nh,Ny,Nz),i_kz_rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_x_1_k(Nh,Ny,Nz),i_ky_Mom_x_2_k(Nh,Ny,Nz),i_kz_Mom_x_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_y_1_k(Nh,Ny,Nz),i_ky_Mom_y_2_k(Nh,Ny,Nz),i_kz_Mom_y_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_z_1_k(Nh,Ny,Nz),i_ky_Mom_z_2_k(Nh,Ny,Nz),i_kz_Mom_z_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Energy_x_k(Nh,Ny,Nz),i_ky_Energy_y_k(Nh,Ny,Nz),i_kz_Energy_z_k(Nh,Ny,Nz),E_Visc_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_ky_Mag_x_1_k(Nh,Ny,Nz),i_kz_Mag_x_2_k(Nh,Ny,Nz),i_kx_Mag_y_1_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kz_Mag_y_2_k(Nh,Ny,Nz),i_kx_Mag_z_1_k(Nh,Ny,Nz),i_ky_Mag_z_2_k(Nh,Ny,Nz)
complex ( kind = 8 ) kx2_ux_k(Nh,Ny,Nz),ky2_ux_k(Nh,Ny,Nz),kz2_ux_k(Nh,Ny,Nz),kx2_uy_k(Nh,Ny,Nz),ky2_uy_k(Nh,Ny,Nz)
complex ( kind = 8 ) kz2_uy_k(Nh,Ny,Nz),kx2_uz_k(Nh,Ny,Nz),ky2_uz_k(Nh,Ny,Nz),kz2_uz_k(Nh,Ny,Nz),nu(Nh,Ny,Nz)
complex ( kind = 8 ) kx2_Bx_k(Nh,Ny,Nz),ky2_Bx_k(Nh,Ny,Nz),kz2_Bx_k(Nh,Ny,Nz),kx2_By_k(Nh,Ny,Nz),ky2_By_k(Nh,Ny,Nz)
complex ( kind = 8 ) kz2_By_k(Nh,Ny,Nz),kx2_Bz_k(Nh,Ny,Nz),ky2_Bz_k(Nh,Ny,Nz),kz2_Bz_k(Nh,Ny,Nz)

!$OMP PARALLEL SHARED(Lx,Ly,Lz,i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k),&
!$OMP & SHARED(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k),&
!$OMP & SHARED(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k),&
!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,E_Visc_k)&
!$OMP & SHARED(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k),&
!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k,nu),&
!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO 

do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
        if (kx >= Nx/3 .and. ky >= Ny/3 .and. kz >= Nz/3) then!dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k(i,j,k) = 0.0d0
        i_ky_rho_uy_k(i,j,k) = 0.0d0
        i_kz_rho_uz_k(i,j,k) = 0.0d0
        i_kx_Mom_x_1_k(i,j,k) = 0.0d0
        i_ky_Mom_x_2_k(i,j,k) = 0.0d0
        i_kz_Mom_x_3_k(i,j,k) = 0.0d0
        i_kx_Mom_y_1_k(i,j,k) = 0.0d0
        i_ky_Mom_y_2_k(i,j,k) = 0.0d0
        i_kz_Mom_y_3_k(i,j,k) = 0.0d0
        i_kx_Mom_z_1_k(i,j,k) = 0.0d0
        i_ky_Mom_z_2_k(i,j,k) = 0.0d0
        i_kz_Mom_z_3_k(i,j,k) = 0.0d0
        kx2_ux_k(i,j,k) = 0.0d0
        ky2_ux_k(i,j,k) = 0.0d0
        kz2_ux_k(i,j,k) = 0.0d0
        kx2_uy_k(i,j,k) = 0.0d0
        ky2_uy_k(i,j,k) = 0.0d0
        kz2_uy_k(i,j,k) = 0.0d0
        kx2_uz_k(i,j,k) = 0.0d0
        ky2_uz_k(i,j,k) = 0.0d0
        kz2_uz_k(i,j,k) = 0.0d0
        i_kx_Energy_x_k(i,j,k) = 0.0d0
        i_ky_Energy_y_k(i,j,k) = 0.0d0
        i_kz_Energy_z_k(i,j,k) = 0.0d0
        E_Visc_k(i,j,k) = 0.0d0
        i_ky_Mag_x_1_k(i,j,k) = 0.0d0
        i_kz_Mag_x_2_k(i,j,k) = 0.0d0
        i_kx_Mag_y_1_k(i,j,k) = 0.0d0
        i_kz_Mag_y_2_k(i,j,k) = 0.0d0
        i_kx_Mag_z_1_k(i,j,k) = 0.0d0
        i_ky_Mag_z_2_k(i,j,k) = 0.0d0
        kx2_Bx_k(i,j,k) = 0.0d0
        ky2_Bx_k(i,j,k) = 0.0d0
        kz2_Bx_k(i,j,k) = 0.0d0
        kx2_By_k(i,j,k) = 0.0d0
        ky2_By_k(i,j,k) = 0.0d0
        kz2_By_k(i,j,k) = 0.0d0
        kx2_Bz_k(i,j,k) = 0.0d0
        ky2_Bz_k(i,j,k) = 0.0d0
        kz2_Bz_k(i,j,k) = 0.0d0
        nu(i,j,k) = 0.0d0
        endif
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float((k-1)-Nz)/Lz
        if (kx >= Nx/3 .and. ky >= Ny/3 .and. kz >= Nz/3) then!dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k(i,j,k) = 0.0d0
        i_ky_rho_uy_k(i,j,k) = 0.0d0
        i_kz_rho_uz_k(i,j,k) = 0.0d0
        i_kx_Mom_x_1_k(i,j,k) = 0.0d0
        i_ky_Mom_x_2_k(i,j,k) = 0.0d0
        i_kz_Mom_x_3_k(i,j,k) = 0.0d0
        i_kx_Mom_y_1_k(i,j,k) = 0.0d0
        i_ky_Mom_y_2_k(i,j,k) = 0.0d0
        i_kz_Mom_y_3_k(i,j,k) = 0.0d0
        i_kx_Mom_z_1_k(i,j,k) = 0.0d0
        i_ky_Mom_z_2_k(i,j,k) = 0.0d0
        i_kz_Mom_z_3_k(i,j,k) = 0.0d0
        kx2_ux_k(i,j,k) = 0.0d0
        ky2_ux_k(i,j,k) = 0.0d0
        kz2_ux_k(i,j,k) = 0.0d0
        kx2_uy_k(i,j,k) = 0.0d0
        ky2_uy_k(i,j,k) = 0.0d0
        kz2_uy_k(i,j,k) = 0.0d0
        kx2_uz_k(i,j,k) = 0.0d0
        ky2_uz_k(i,j,k) = 0.0d0
        kz2_uz_k(i,j,k) = 0.0d0
        i_kx_Energy_x_k(i,j,k) = 0.0d0
        i_ky_Energy_y_k(i,j,k) = 0.0d0
        i_kz_Energy_z_k(i,j,k) = 0.0d0
        E_Visc_k(i,j,k) = 0.0d0
        i_ky_Mag_x_1_k(i,j,k) = 0.0d0
        i_kz_Mag_x_2_k(i,j,k) = 0.0d0
        i_kx_Mag_y_1_k(i,j,k) = 0.0d0
        i_kz_Mag_y_2_k(i,j,k) = 0.0d0
        i_kx_Mag_z_1_k(i,j,k) = 0.0d0
        i_ky_Mag_z_2_k(i,j,k) = 0.0d0
        kx2_Bx_k(i,j,k) = 0.0d0
        ky2_Bx_k(i,j,k) = 0.0d0
        kz2_Bx_k(i,j,k) = 0.0d0
        kx2_By_k(i,j,k) = 0.0d0
        ky2_By_k(i,j,k) = 0.0d0
        kz2_By_k(i,j,k) = 0.0d0
        kx2_Bz_k(i,j,k) = 0.0d0
        ky2_Bz_k(i,j,k) = 0.0d0
        kz2_Bz_k(i,j,k) = 0.0d0
        nu(i,j,k) = 0.0d0
        endif
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k-1)/Lz
        if (kx >= Nx/3 .and. ky >= Ny/3 .and. kz >= Nz/3) then!dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k(i,j,k) = 0.0d0
        i_ky_rho_uy_k(i,j,k) = 0.0d0
        i_kz_rho_uz_k(i,j,k) = 0.0d0
        i_kx_Mom_x_1_k(i,j,k) = 0.0d0
        i_ky_Mom_x_2_k(i,j,k) = 0.0d0
        i_kz_Mom_x_3_k(i,j,k) = 0.0d0
        i_kx_Mom_y_1_k(i,j,k) = 0.0d0
        i_ky_Mom_y_2_k(i,j,k) = 0.0d0
        i_kz_Mom_y_3_k(i,j,k) = 0.0d0
        i_kx_Mom_z_1_k(i,j,k) = 0.0d0
        i_ky_Mom_z_2_k(i,j,k) = 0.0d0

        i_kz_Mom_z_3_k(i,j,k) = 0.0d0
        kx2_ux_k(i,j,k) = 0.0d0
        ky2_ux_k(i,j,k) = 0.0d0
        kz2_ux_k(i,j,k) = 0.0d0
        kx2_uy_k(i,j,k) = 0.0d0
        ky2_uy_k(i,j,k) = 0.0d0
        kz2_uy_k(i,j,k) = 0.0d0
        kx2_uz_k(i,j,k) = 0.0d0
        ky2_uz_k(i,j,k) = 0.0d0
        kz2_uz_k(i,j,k) = 0.0d0
        i_kx_Energy_x_k(i,j,k) = 0.0d0
        i_ky_Energy_y_k(i,j,k) = 0.0d0
        i_kz_Energy_z_k(i,j,k) = 0.0d0
        E_Visc_k(i,j,k) = 0.0d0
        i_ky_Mag_x_1_k(i,j,k) = 0.0d0
        i_kz_Mag_x_2_k(i,j,k) = 0.0d0
        i_kx_Mag_y_1_k(i,j,k) = 0.0d0
        i_kz_Mag_y_2_k(i,j,k) = 0.0d0
        i_kx_Mag_z_1_k(i,j,k) = 0.0d0
        i_ky_Mag_z_2_k(i,j,k) = 0.0d0
        kx2_Bx_k(i,j,k) = 0.0d0
        ky2_Bx_k(i,j,k) = 0.0d0
        kz2_Bx_k(i,j,k) = 0.0d0
        kx2_By_k(i,j,k) = 0.0d0
        ky2_By_k(i,j,k) = 0.0d0
        kz2_By_k(i,j,k) = 0.0d0
        kx2_Bz_k(i,j,k) = 0.0d0
        ky2_Bz_k(i,j,k) = 0.0d0
        kz2_Bz_k(i,j,k) = 0.0d0
        nu(i,j,k) = 0.0d0
        endif
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float((k-1)-Nz)/Lz
        if (kx >= Nx/3 .and. ky >= Ny/3 .and. kz >= Nz/3) then!dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k(i,j,k) = 0.0d0
        i_ky_rho_uy_k(i,j,k) = 0.0d0
        i_kz_rho_uz_k(i,j,k) = 0.0d0
        i_kx_Mom_x_1_k(i,j,k) = 0.0d0
        i_ky_Mom_x_2_k(i,j,k) = 0.0d0
        i_kz_Mom_x_3_k(i,j,k) = 0.0d0
        i_kx_Mom_y_1_k(i,j,k) = 0.0d0
        i_ky_Mom_y_2_k(i,j,k) = 0.0d0
        i_kz_Mom_y_3_k(i,j,k) = 0.0d0
        i_kx_Mom_z_1_k(i,j,k) = 0.0d0
        i_ky_Mom_z_2_k(i,j,k) = 0.0d0
        i_kz_Mom_z_3_k(i,j,k) = 0.0d0
        kx2_ux_k(i,j,k) = 0.0d0
        ky2_ux_k(i,j,k) = 0.0d0
        kz2_ux_k(i,j,k) = 0.0d0
        kx2_uy_k(i,j,k) = 0.0d0
        ky2_uy_k(i,j,k) = 0.0d0
        kz2_uy_k(i,j,k) = 0.0d0
        kx2_uz_k(i,j,k) = 0.0d0
        ky2_uz_k(i,j,k) = 0.0d0
        kz2_uz_k(i,j,k) = 0.0d0
        i_kx_Energy_x_k(i,j,k) = 0.0d0
        i_ky_Energy_y_k(i,j,k) = 0.0d0
        i_kz_Energy_z_k(i,j,k) = 0.0d0
        E_Visc_k(i,j,k) = 0.0d0
        i_ky_Mag_x_1_k(i,j,k) = 0.0d0
        i_kz_Mag_x_2_k(i,j,k) = 0.0d0
        i_kx_Mag_y_1_k(i,j,k) = 0.0d0
        i_kz_Mag_y_2_k(i,j,k) = 0.0d0
        i_kx_Mag_z_1_k(i,j,k) = 0.0d0
        i_ky_Mag_z_2_k(i,j,k) = 0.0d0
        kx2_Bx_k(i,j,k) = 0.0d0
        ky2_Bx_k(i,j,k) = 0.0d0
        kz2_Bx_k(i,j,k) = 0.0d0
        kx2_By_k(i,j,k) = 0.0d0
        ky2_By_k(i,j,k) = 0.0d0
        kz2_By_k(i,j,k) = 0.0d0
        kx2_Bz_k(i,j,k) = 0.0d0
        ky2_Bz_k(i,j,k) = 0.0d0
        kz2_Bz_k(i,j,k) = 0.0d0
        nu(i,j,k) = 0.0d0
        endif
    enddo
  enddo
enddo  

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine DE_ALIAZING

!====================================================================================

