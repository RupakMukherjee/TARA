!===================================================================================
!================== SUBROUTINE NONLINEAR DERIVATIVE ================================
!===================================================================================

subroutine DERIVE(Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
                  d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
                  d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
                  d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)
use omp_lib
implicit none

include "fftw3.f"

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,time,dt,Lx,Ly,Lz,dx,dy,dz,kx,ky,kz,mu,ms,CS,mu_0,eta,spheat,kf,A,B,C,time_max

real ( kind = 8 ) x(Nx),y(Ny),z(Nz),rho(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz),E(Nx,Ny,Nz),P(Nx,Ny,Nz)
real ( kind = 8 ) ux_dum(Nx,Ny,Nz),uy_dum(Nx,Ny,Nz),uz_dum(Nx,Ny,Nz)
real ( kind = 8 ) Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz)
real ( kind = 8 ) rho_ux(Nx,Ny,Nz),rho_uy(Nx,Ny,Nz),rho_uz(Nx,Ny,Nz)
real ( kind = 8 ) Mom_x_1(Nx,Ny,Nz),Mom_x_2(Nx,Ny,Nz),Mom_x_3(Nx,Ny,Nz),Mom_y_1(Nx,Ny,Nz),Mom_y_2(Nx,Ny,Nz),Mom_y_3(Nx,Ny,Nz)
real ( kind = 8 ) Mom_z_1(Nx,Ny,Nz),Mom_z_2(Nx,Ny,Nz),Mom_z_3(Nx,Ny,Nz)

real ( kind = 8 ) d_ux_dx(Nx,Ny,Nz),d_uy_dy(Nx,Ny,Nz),d_uz_dz(Nx,Ny,Nz),Fx(Nx,Ny,Nz),Fy(Nx,Ny,Nz),Fz(Nx,Ny,Nz)
real ( kind = 8 ) Energy_x(Nx,Ny,Nz),Energy_y(Nx,Ny,Nz),Energy_z(Nx,Ny,Nz),E_Visc(Nx,Ny,Nz)
real ( kind = 8 ) Mag_x_1(Nx,Ny,Nz),Mag_x_2(Nx,Ny,Nz),Mag_y_1(Nx,Ny,Nz),Mag_y_2(Nx,Ny,Nz),Mag_z_1(Nx,Ny,Nz),Mag_z_2(Nx,Ny,Nz)
real ( kind = 8 ) d_Bx_dy(Nx,Ny,Nz),d_By_dx(Nx,Ny,Nz),d_Bx_dz(Nx,Ny,Nz)
real ( kind = 8 ) d_By_dz(Nx,Ny,Nz),d_Bz_dx(Nx,Ny,Nz),d_Bz_dy(Nx,Ny,Nz)
real ( kind = 8 ) curl_x_B(Nx,Ny,Nz),curl_y_B(Nx,Ny,Nz),curl_z_B(Nx,Ny,Nz)

complex ( kind = 8 ) rho_k(Nh,Ny,Nz),ux_k(Nh,Ny,Nz),uy_k(Nh,Ny,Nz),uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) nu(Nh,Ny,Nz),Fx_k(Nh,Ny,Nz),Fy_k(Nh,Ny,Nz),Fz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_ux_k(Nh,Ny,Nz),rho_uy_k(Nh,Ny,Nz),rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_ux_k_dum(Nh,Ny,Nz),rho_uy_k_dum(Nh,Ny,Nz),rho_uz_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) P_k(Nh,Ny,Nz),E_k(Nh,Ny,Nz),Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_k_dum(Nh,Ny,Nz),ux_k_dum(Nh,Ny,Nz),uy_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) E_k_dum(Nh,Ny,Nz),Bx_k_dum(Nh,Ny,Nz),By_k_dum(Nh,Ny,Nz),Bz_k_dum(Nh,Ny,Nz)

complex ( kind = 8 ) i_kx_rho_ux_k(Nh,Ny,Nz),i_ky_rho_uy_k(Nh,Ny,Nz),i_kz_rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mom_x_1_k(Nh,Ny,Nz),Mom_x_2_k(Nh,Ny,Nz),Mom_x_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mom_y_1_k(Nh,Ny,Nz),Mom_y_2_k(Nh,Ny,Nz),Mom_y_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mom_z_1_k(Nh,Ny,Nz),Mom_z_2_k(Nh,Ny,Nz),Mom_z_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_x_1_k(Nh,Ny,Nz),i_ky_Mom_x_2_k(Nh,Ny,Nz),i_kz_Mom_x_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_y_1_k(Nh,Ny,Nz),i_ky_Mom_y_2_k(Nh,Ny,Nz),i_kz_Mom_y_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_z_1_k(Nh,Ny,Nz),i_ky_Mom_z_2_k(Nh,Ny,Nz),i_kz_Mom_z_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) kx2_ux_k(Nh,Ny,Nz),ky2_ux_k(Nh,Ny,Nz),kz2_ux_k(Nh,Ny,Nz),kx2_uy_k(Nh,Ny,Nz),ky2_uy_k(Nh,Ny,Nz)
complex ( kind = 8 ) kz2_uy_k(Nh,Ny,Nz),kx2_uz_k(Nh,Ny,Nz),ky2_uz_k(Nh,Ny,Nz),kz2_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_ux_k(Nh,Ny,Nz),i_ky_uy_k(Nh,Ny,Nz),i_kz_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) Energy_x_k(Nh,Ny,Nz),Energy_y_k(Nh,Ny,Nz),Energy_z_k(Nh,Ny,Nz),E_Visc_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mag_x_1_k(Nh,Ny,Nz),Mag_x_2_k(Nh,Ny,Nz),Mag_y_1_k(Nh,Ny,Nz),Mag_y_2_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mag_z_1_k(Nh,Ny,Nz),Mag_z_2_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Energy_x_k(Nh,Ny,Nz),i_ky_Energy_y_k(Nh,Ny,Nz),i_kz_Energy_z_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_ky_Mag_x_1_k(Nh,Ny,Nz),i_kz_Mag_x_2_k(Nh,Ny,Nz),i_kx_Mag_y_1_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kz_Mag_y_2_k(Nh,Ny,Nz),i_kx_Mag_z_1_k(Nh,Ny,Nz),i_ky_Mag_z_2_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_ky_Bx_k(Nh,Ny,Nz),i_kx_By_k(Nh,Ny,Nz),i_kx_Bz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kz_Bx_k(Nh,Ny,Nz),i_ky_Bz_k(Nh,Ny,Nz),i_kz_By_k(Nh,Ny,Nz)
complex ( kind = 8 ) kx2_Bx_k(Nh,Ny,Nz),ky2_Bx_k(Nh,Ny,Nz),kz2_Bx_k(Nh,Ny,Nz),kx2_By_k(Nh,Ny,Nz),ky2_By_k(Nh,Ny,Nz)
complex ( kind = 8 ) kz2_By_k(Nh,Ny,Nz),kx2_Bz_k(Nh,Ny,Nz),ky2_Bz_k(Nh,Ny,Nz),kz2_Bz_k(Nh,Ny,Nz)

complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny,Nz),d_rho_ux_k_dt_old(Nh,Ny,Nz),d_rho_uy_k_dt_old(Nh,Ny,Nz),d_rho_uz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_old(Nh,Ny,Nz),d_Bx_k_dt_old(Nh,Ny,Nz),d_By_k_dt_old(Nh,Ny,Nz),d_Bz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny,Nz),d_rho_ux_k_dt_new(Nh,Ny,Nz),d_rho_uy_k_dt_new(Nh,Ny,Nz),d_rho_uz_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_new(Nh,Ny,Nz),d_Bx_k_dt_new(Nh,Ny,Nz),d_By_k_dt_new(Nh,Ny,Nz),d_Bz_k_dt_new(Nh,Ny,Nz)

integer ( kind = 8 ) plan_forward,plan_backward ! FFTW

common/comm/time_max,Lx,Ly,Lz,dx,dy,dz,spheat,mu,ms,CS,mu_0,eta,dt,kf

! Keep Backup of the Arrays for FFTW.

!$OMP PARALLEL SHARED(rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k),&
!$OMP & SHARED(rho_k_dum,rho_ux_k_dum,rho_uy_k_dum,rho_uz_k_dum),&
!$OMP & SHARED(E_k_dum,Bx_k_dum,By_k_dum,Bz_k_dum) PRIVATE(i,j,k)
!$OMP DO
 
do i = 1,Nx/2+1
  do j = 1,Ny
    do k = 1,Nz
    rho_k_dum(i,j,k) = rho_k(i,j,k)  
    rho_ux_k_dum(i,j,k) = rho_ux_k(i,j,k)
    rho_uy_k_dum(i,j,k) = rho_uy_k(i,j,k)
    rho_uz_k_dum(i,j,k) = rho_uz_k(i,j,k)
    E_k_dum(i,j,k) = E_k(i,j,k)
    Bx_k_dum(i,j,k) = Bx_k(i,j,k)
    By_k_dum(i,j,k) = By_k(i,j,k)
    Bz_k_dum(i,j,k) = Bz_k(i,j,k)
    enddo
  enddo
enddo  

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_k_dum, rho, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_ux_k_dum, rho_ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uy_k_dum, rho_uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uz_k_dum, rho_uz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, E_k_dum, E, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bx_k_dum, Bx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, By_k_dum, By, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bz_k_dum, Bz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

! Evaluate Velocity Force and Magnetic Field.

  call Vel_Force_Mag(Nx,Ny,Nz,Nh,pi,dx,dy,dz,A,B,C,kf,rho,ux,uy,uz,ux_dum,uy_dum,uz_dum,&
                     rho_ux,rho_uy,rho_uz,Fx,Fy,Fz,E,Bx,By,Bz,B2)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fx, Fx_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fy, Fy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fz, Fz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

! Evaluate derivatives of Velocity and Magnetic Field.

  call Vel_Mag_Derive(Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,ux_k,uy_k,uz_k,Bx_k,By_k,Bz_k,&
                      i_kx_ux_k,i_ky_uy_k,i_kz_uz_k,&
                      i_ky_Bx_k,i_kz_Bx_k,i_kx_By_k,i_kz_By_k,i_kx_Bz_k,i_ky_Bz_k) 

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_ux_k, d_ux_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_uy_k, d_uy_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_uz_k, d_uz_dz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_Bx_k, d_Bx_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_Bx_k, d_Bx_dz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_By_k, d_By_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_By_k, d_By_dz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_Bz_k, d_Bz_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_Bz_k, d_Bz_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

! Evaluate the Multiplications in Real Space.

  call Real_Evaluation(Nx,Ny,Nz,Nh,pi,spheat,eta,ux,uy,uz,rho_ux,rho_uy,rho_uz,Bx,By,Bz,&
                       d_ux_dx,d_uy_dy,d_uz_dz,P,E,&
                       d_Bx_dy,d_Bx_dz,d_By_dx,d_By_dz,d_Bz_dx,d_Bz_dy,curl_x_B,curl_y_B,curl_z_B,B2,&
                       Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3,&
                       Energy_x,Energy_y,Energy_z,E_Visc,&
                       Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_1, Mom_x_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_2, Mom_x_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_3, Mom_x_3_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_1, Mom_y_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_2, Mom_y_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_3, Mom_y_3_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_1, Mom_z_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_2, Mom_z_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_3, Mom_z_3_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_x, Energy_x_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_y, Energy_y_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_z, Energy_z_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, E_Visc, E_Visc_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_x_1, Mag_x_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_x_2, Mag_x_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_y_1, Mag_y_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_y_2, Mag_y_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_z_1, Mag_z_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_z_2, Mag_z_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_destroy_plan_ (plan_backward)
  call dfftw_destroy_plan_ (plan_forward)

! Evaluate the Derivatives in Spectral Space.

  call Spectral_Evaluation(Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,ux_k,uy_k,uz_k,rho_ux_k,rho_uy_k,rho_uz_k,Bx_k,By_k,Bz_k,&
                           Mom_x_1_k,Mom_x_2_k,Mom_x_3_k,Mom_y_1_k,Mom_y_2_k,Mom_y_3_k,Mom_z_1_k,Mom_z_2_k,Mom_z_3_k,&
                           Energy_x_k,Energy_y_k,Energy_z_k,Mag_x_1_k,Mag_x_2_k,Mag_y_1_k,Mag_y_2_k,Mag_z_1_k,Mag_z_2_k,&
                           i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k,&
                           i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k,&
                           i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k,&
                           i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,&
                           i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k,&
                           kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k,&
                           kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k)


! De - Aliazing Technique With 2/3 Rule for All the Non-Linear Terms. 
 
  call DE_ALIAZING(Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k,&
                   i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k,&
                   i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k,&
                   i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,E_Visc_k,&
                   i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k,&
                   kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k,nu,&
                   kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k)

! Evaluate the Non Linear Derivatives in Spectral Space.

  call Non_Linear_Derivative(Nx,Ny,Nz,Nh,pi,mu,eta,i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k,&
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
    
return

end subroutine DERIVE

!====================================================================================

