! Evaluate derivatives of Velocity and Magnetic Field.

!===================================================================================
!================== SUBROUTINE Velocity - Magnetic Field Derivative ================
!===================================================================================

subroutine Vel_Mag_Derive(Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,ux_k,uy_k,uz_k,Bx_k,By_k,Bz_k,&
                          i_kx_ux_k,i_ky_uy_k,i_kz_uz_k,&
                          i_ky_Bx_k,i_kz_Bx_k,i_kx_By_k,i_kz_By_k,i_kx_Bz_k,i_ky_Bz_k)
                          
use omp_lib
implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,Lx,Ly,Lz,kx,ky,kz

complex ( kind = 8 ) ux_k(Nh,Ny,Nz),uy_k(Nh,Ny,Nz),uz_k(Nh,Ny,Nz)                          
complex ( kind = 8 ) Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_ux_k(Nh,Ny,Nz),i_ky_uy_k(Nh,Ny,Nz),i_kz_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_ky_Bx_k(Nh,Ny,Nz),i_kx_By_k(Nh,Ny,Nz),i_kx_Bz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kz_Bx_k(Nh,Ny,Nz),i_ky_Bz_k(Nh,Ny,Nz),i_kz_By_k(Nh,Ny,Nz)

!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,Bx_k,By_k,Bz_k),&
!$OMP & SHARED(i_kx_ux_k,i_ky_uy_k,i_kz_uz_k),&
!$OMP & SHARED(i_ky_Bx_k,i_kz_Bx_k,i_kx_By_k,i_kz_By_k,i_kx_Bz_k,i_ky_Bz_k) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO
 
do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*float(k-1)/Lz
      i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
      i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
      i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)
      
      i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
      i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
      i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
      i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
    enddo 
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float((k-1)-Nz)/Lz
      i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
      i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
      i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)
      
      i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
      i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
      i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
      i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k-1)/Lz
      i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
      i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
      i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)
      
      i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
      i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
      i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
      i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float((k-1)-Nz)/Lz
      i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
      i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
      i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)
      
      i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
      i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
      i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
      i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
    enddo
  enddo   
enddo

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine Vel_Mag_Derive

!====================================================================================

