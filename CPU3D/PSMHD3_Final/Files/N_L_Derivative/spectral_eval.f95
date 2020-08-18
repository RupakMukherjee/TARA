! Evaluate the Derivatives in Spectral Space.

!===================================================================================
!================== SUBROUTINE Spectral - Evaluation ===============================
!===================================================================================

subroutine Spectral_Evaluation(Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,ux_k,uy_k,uz_k,rho_ux_k,rho_uy_k,rho_uz_k,Bx_k,By_k,Bz_k,&
                               Mom_x_1_k,Mom_x_2_k,Mom_x_3_k,Mom_y_1_k,Mom_y_2_k,Mom_y_3_k,Mom_z_1_k,Mom_z_2_k,Mom_z_3_k,&
                               Energy_x_k,Energy_y_k,Energy_z_k,Mag_x_1_k,Mag_x_2_k,Mag_y_1_k,Mag_y_2_k,Mag_z_1_k,Mag_z_2_k,&
                               i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k,&
                               i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k,&
                               i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k,&
                               i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,&
                               i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k,&
                               kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k,&
                               kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k)

use omp_lib
implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,Lx,Ly,Lz,kx,ky,kz

complex ( kind = 8 ) ux_k(Nh,Ny,Nz),uy_k(Nh,Ny,Nz),uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_ux_k(Nh,Ny,Nz),rho_uy_k(Nh,Ny,Nz),rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mom_x_1_k(Nh,Ny,Nz),Mom_x_2_k(Nh,Ny,Nz),Mom_x_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mom_y_1_k(Nh,Ny,Nz),Mom_y_2_k(Nh,Ny,Nz),Mom_y_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mom_z_1_k(Nh,Ny,Nz),Mom_z_2_k(Nh,Ny,Nz),Mom_z_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) Energy_x_k(Nh,Ny,Nz),Energy_y_k(Nh,Ny,Nz),Energy_z_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mag_x_1_k(Nh,Ny,Nz),Mag_x_2_k(Nh,Ny,Nz),Mag_y_1_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mag_y_2_k(Nh,Ny,Nz),Mag_z_1_k(Nh,Ny,Nz),Mag_z_2_k(Nh,Ny,Nz)

complex ( kind = 8 ) i_kx_rho_ux_k(Nh,Ny,Nz),i_ky_rho_uy_k(Nh,Ny,Nz),i_kz_rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_x_1_k(Nh,Ny,Nz),i_ky_Mom_x_2_k(Nh,Ny,Nz),i_kz_Mom_x_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_y_1_k(Nh,Ny,Nz),i_ky_Mom_y_2_k(Nh,Ny,Nz),i_kz_Mom_y_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_z_1_k(Nh,Ny,Nz),i_ky_Mom_z_2_k(Nh,Ny,Nz),i_kz_Mom_z_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Energy_x_k(Nh,Ny,Nz),i_ky_Energy_y_k(Nh,Ny,Nz),i_kz_Energy_z_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_ky_Mag_x_1_k(Nh,Ny,Nz),i_kz_Mag_x_2_k(Nh,Ny,Nz),i_kx_Mag_y_1_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kz_Mag_y_2_k(Nh,Ny,Nz),i_kx_Mag_z_1_k(Nh,Ny,Nz),i_ky_Mag_z_2_k(Nh,Ny,Nz)
complex ( kind = 8 ) kx2_ux_k(Nh,Ny,Nz),ky2_ux_k(Nh,Ny,Nz),kz2_ux_k(Nh,Ny,Nz),kx2_uy_k(Nh,Ny,Nz),ky2_uy_k(Nh,Ny,Nz)
complex ( kind = 8 ) kz2_uy_k(Nh,Ny,Nz),kx2_uz_k(Nh,Ny,Nz),ky2_uz_k(Nh,Ny,Nz),kz2_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) kx2_Bx_k(Nh,Ny,Nz),ky2_Bx_k(Nh,Ny,Nz),kz2_Bx_k(Nh,Ny,Nz),kx2_By_k(Nh,Ny,Nz),ky2_By_k(Nh,Ny,Nz)
complex ( kind = 8 ) kz2_By_k(Nh,Ny,Nz),kx2_Bz_k(Nh,Ny,Nz),ky2_Bz_k(Nh,Ny,Nz),kz2_Bz_k(Nh,Ny,Nz)

!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,rho_ux_k,rho_uy_k,rho_uz_k,Bx_k,By_k,Bz_k),&
!$OMP & SHARED(Mom_x_1_k,Mom_x_2_k,Mom_x_3_k,Mom_y_1_k,Mom_y_2_k,Mom_y_3_k),&
!$OMP & SHARED(Mom_z_1_k,Mom_z_2_k,Mom_z_3_k),&
!$OMP & SHARED(Energy_x_k,Energy_y_k,Energy_z_k,Mag_x_1_k,Mag_x_2_k,Mag_y_1_k,Mag_y_2_k,Mag_z_1_k,Mag_z_2_k),&
!$OMP & SHARED(i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k),&
!$OMP & SHARED(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k),&
!$OMP & SHARED(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k),&
!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k),&
!$OMP & SHARED(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k),&
!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k),&
!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO 

do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*float(k-1)/Lz
  
      i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
      i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k) 
      i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k) 
        
      i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)    
      i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
      i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)   
      i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
      i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
      i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k) 
      i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)    
      i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)    
      i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)
     
      kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
      ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
      kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
      kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
      ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k) 
      kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
      kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
      ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
      kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

      i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
      i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
      i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)
    
      i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
      i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
      i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
      i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
      i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
      i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)
    
      kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
      ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
      kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
      kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
      ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k) 
      kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
      kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
      ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
      kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*float(i-1)/Lx
    ky = 2.0d0*pi*float(j-1)/Ly
    kz = 2.0d0*pi*float((k-1)-Nz)/Lz
    
      i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
      i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k) 
      i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k) 
        
      i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)    
      i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
      i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)   
      i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
      i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
      i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k) 
      i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)    
      i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)    
      i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)
     
      kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
      ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
      kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
      kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
      ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k) 
      kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
      kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
      ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
      kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

      i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
      i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
      i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)
    
      i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
      i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
      i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
      i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
      i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
      i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)
    
      kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
      ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
      kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
      kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
      ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k) 
      kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
      kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
      ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
      kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
    enddo
  enddo
  
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
    kx = 2.0d0*pi*float(i-1)/Lx
    ky = 2.0d0*pi*float((j-1)-Ny)/Ly
    kz = 2.0d0*pi*float(k-1)/Lz

      i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
      i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k) 
      i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k) 
        
      i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)    
      i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
      i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)   
      i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
      i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
      i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k) 
      i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)    
      i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)    
      i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)
     
      kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
      ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
      kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
      kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
      ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k) 
      kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
      kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
      ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
      kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

      i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
      i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
      i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)
    
      i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
      i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
      i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
      i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
      i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
      i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)
    
      kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
      ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
      kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
      kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
      ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k) 
      kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
      kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
      ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
      kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*float(i-1)/Lx
    ky = 2.0d0*pi*float((j-1)-Ny)/Ly
    kz = 2.0d0*pi*float((k-1)-Nz)/Lz

      i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
      i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k) 
      i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k) 
        
      i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)    
      i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
      i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)   
      i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
      i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
      i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k) 
      i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)    
      i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)    
      i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)
     
      kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
      ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
      kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
      kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
      ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k) 
      kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
      kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
      ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
      kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

      i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
      i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
      i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)
    
      i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
      i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
      i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
      i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
      i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
      i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)
    
      kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
      ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
      kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
      kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
      ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k) 
      kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
      kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
      ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
      kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine Spectral_Evaluation

!====================================================================================

