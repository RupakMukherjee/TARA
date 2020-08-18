! Evaluate Vorticity in Spectral Space.

!===================================================================================
!================== SUBROUTINE OMEGA ===============================================
!===================================================================================

subroutine OMEGA (Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,ux_k,uy_k,uz_k,omega_x_k,omega_y_k,omega_z_k,&
                  omega_x,omega_y,omega_z,omega2)
                  
use omp_lib
implicit none

include "fftw3.f"

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,Lx,Ly,Lz,kx,ky,kz

real ( kind = 8 ) omega_x(Nx,Ny,Nz),omega_y(Nx,Ny,Nz),omega_z(Nx,Ny,Nz),omega2(Nx,Ny,Nz)

complex ( kind = 8 ) ux_k(Nh,Ny,Nz),uy_k(Nh,Ny,Nz),uz_k(Nh,Ny,Nz)   
complex ( kind = 8 ) omega_x_k(Nh,Ny,Nz),omega_y_k(Nh,Ny,Nz),omega_z_k(Nh,Ny,Nz) 
complex ( kind = 8 ) omega_x_k_dum(Nh,Ny,Nz),omega_y_k_dum(Nh,Ny,Nz),omega_z_k_dum(Nh,Ny,Nz)  
         
integer ( kind = 8 ) plan_backward ! FFTW

!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,omega_x_k,omega_y_k,omega_z_k),&
!$OMP & SHARED(omega_x_k_dum,omega_y_k_dum,omega_z_k_dum) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO 

do i = 1, Nx/2+1
  do j = 1, Ny/2
    do k = 1, Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
    do k = Nz/2+1, Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo 
  enddo
  do j = Ny/2+1, Ny
    do k = 1, Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
    do k = Nz/2+1, Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)
      
      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo 

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_x_k_dum, omega_x, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_y_k_dum, omega_y, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_z_k_dum, omega_z, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_destroy_plan_ (plan_backward)

! FFTW Normalisation and omega^2 Evaluation.

!$OMP PARALLEL SHARED(omega_x,omega_y,omega_z,omega2) PRIVATE(i,j,k)
!$OMP DO 

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      omega_x(i,j,k) = omega_x(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega_y(i,j,k) = omega_y(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega_z(i,j,k) = omega_z(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega2(i,j,k) = omega_x(i,j,k)**2 + omega_y(i,j,k)**2 + omega_z(i,j,k)**2
    end do
  end do
enddo

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine OMEGA

!====================================================================================

