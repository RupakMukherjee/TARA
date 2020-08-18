! Evaluate Divergence of Magnetic Field.

!===================================================================================
!================== SUBROUTINE DIVERGENCE_B ========================================
!===================================================================================

subroutine DIVERGENCE_B (Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,Bx_k,By_k,Bz_k,div_B)

use omp_lib
implicit none

include "fftw3.f"

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,Lx,Ly,Lz,kx,ky,kz,div_B(Nx,Ny,Nz)

complex ( kind = 8 ) Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz),div_B_k(Nh,Ny,Nz)

integer ( kind = 8 ) plan_backward ! FFTW

!$OMP PARALLEL SHARED(Lx,Ly,Lz,Bx_k,By_k,Bz_k,div_B_k) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO
 
do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k) 
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL 

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, div_B_k, div_B, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

!$OMP PARALLEL SHARED(div_B) PRIVATE(i,j,k)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      ! FFTW Normalisation.
      div_B(i,j,k) = div_B(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    enddo
  enddo
enddo   

!$OMP END DO
!$OMP END PARALLEL 
      
return

end subroutine DIVERGENCE_B

!====================================================================================

