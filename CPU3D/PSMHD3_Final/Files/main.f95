Program PSMHD3
 
use omp_lib
implicit none

! Define Grid Size.
integer ( kind = 4 ), parameter :: Nx = 32
integer ( kind = 4 ), parameter :: Ny = 32
integer ( kind = 4 ), parameter :: Nz = 32
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

include "fftw3.f"

integer ( kind = 4 ) G,i,j,k,t
real ( kind = 8 ) Lx,Ly,Lz,dx,dy,dz,kx,ky,kz,time,time_min,time_max,dt,G1,G2,G3
real ( kind = 8 ) spheat,rho0,ms,U0,mu,kf,mode,a,sigma,uy0,M,CS,Re,Rm,PM,P0,eta,mu_0,MA,VA,B0,ba
real ( kind = 8 ) Pressure,Energy,Energy0,y_Energy,y_Energy0,B_Field,C0,C1,C2,C3,div_B_tot

real ( kind = 8 ) x(Nx),y(Ny),z(Nz),rho(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
real ( kind = 8 ) rho_ux(Nx,Ny,Nz),rho_uy(Nx,Ny,Nz),rho_uz(Nx,Ny,Nz)
real ( kind = 8 ) omega_x(Nx,Ny,Nz),omega_y(Nx,Ny,Nz),omega_z(Nx,Ny,Nz),omega2(Nx,Ny,Nz)
real ( kind = 8 ) P(Nx,Ny,Nz),E(Nx,Ny,Nz),Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz),div_B(Nx,Ny,Nz)
real ( kind = 8 ) rho_dum(Nx,Ny,Nz),ux_dum(Nx,Ny,Nz),uy_dum(Nx,Ny,Nz),uz_dum(Nx,Ny,Nz)
real ( kind = 8 ) rho_ux_dum(Nx,Ny,Nz),rho_uy_dum(Nx,Ny,Nz),rho_uz_dum(Nx,Ny,Nz)
real ( kind = 8 ) P_dum(Nx,Ny,Nz),E_dum(Nx,Ny,Nz),Bx_dum(Nx,Ny,Nz),By_dum(Nx,Ny,Nz),Bz_dum(Nx,Ny,Nz)

complex ( kind = 8 ) rho_k(Nh,Ny,Nz),ux_k(Nh,Ny,Nz),uy_k(Nh,Ny,Nz),uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_ux_k(Nh,Ny,Nz),rho_uy_k(Nh,Ny,Nz),rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) omega_x_k(Nh,Ny,Nz),omega_y_k(Nh,Ny,Nz),omega_z_k(Nh,Ny,Nz)
complex ( kind = 8 ) P_k(Nh,Ny,Nz),E_k(Nh,Ny,Nz),Ek(Nh,Ny,Nz),Bk(Nh,Ny,Nz)
complex ( kind = 8 ) Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz),div_B_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_k_dum(Nh,Ny,Nz),ux_k_dum(Nh,Ny,Nz),uy_k_dum(Nh,Ny,Nz),uz_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) omega_x_k_dum(Nh,Ny,Nz),omega_y_k_dum(Nh,Ny,Nz),omega_z_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) rho_ux_k_dum(Nh,Ny,Nz),rho_uy_k_dum(Nh,Ny,Nz),rho_uz_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) E_k_dum(Nh,Ny,Nz),Bx_k_dum(Nh,Ny,Nz),By_k_dum(Nh,Ny,Nz),Bz_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) rho_k_new(Nh,Ny,Nz),rho_ux_k_new(Nh,Ny,Nz),rho_uy_k_new(Nh,Ny,Nz),rho_uz_k_new(Nh,Ny,Nz)
complex ( kind = 8 ) E_k_new(Nh,Ny,Nz),Bx_k_new(Nh,Ny,Nz),By_k_new(Nh,Ny,Nz),Bz_k_new(Nh,Ny,Nz)

complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny,Nz),d_rho_ux_k_dt_old(Nh,Ny,Nz),d_rho_uy_k_dt_old(Nh,Ny,Nz),d_rho_uz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_old(Nh,Ny,Nz),d_Bx_k_dt_old(Nh,Ny,Nz),d_By_k_dt_old(Nh,Ny,Nz),d_Bz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny,Nz),d_rho_ux_k_dt_new(Nh,Ny,Nz),d_rho_uy_k_dt_new(Nh,Ny,Nz),d_rho_uz_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_new(Nh,Ny,Nz),d_Bx_k_dt_new(Nh,Ny,Nz),d_By_k_dt_new(Nh,Ny,Nz),d_Bz_k_dt_new(Nh,Ny,Nz)

integer ( kind = 8 ) plan_forward,plan_backward ! FFTW

integer ( kind = 4 ) thread_num,num_thread,proc_num ! OMP
real ( kind = 8 ) threads,t1,t2

common/comm/time_max,Lx,Ly,Lz,dx,dy,dz,spheat,mu,ms,CS,mu_0,eta,dt,kf

integer,parameter :: seed = 99999999
call srand(seed)

!===================== FILENAMES ==============================================	

open(unit=15,file='Initial_Grid_Data.dat',status='unknown')
open(unit=35,file='Energy_Spectra.dat',status='unknown')
open(unit=40,file='Energy.dat',status='unknown')

!==============================================================================	

! Read Input Parameters.		

  call Input_Parameter (Nx,Ny,Nz,Nh,pi,thread_num,Lx,Ly,Lz,dx,dy,dz,time_min,time_max,dt,&
                        spheat,mu,eta,mu_0,ms,rho0,P0,M,MA,kf,CS,U0,VA,B0) 

  call Initial_Condition (Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,dx,dy,dz,&
                          spheat,rho0,U0,P0,P,B0,rho,x,y,z,ux,uy,uz,Bx,By,Bz)

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      ! Initial Combined Variables.
      rho_ux(i,j,k) = rho(i,j,k) * ux(i,j,k)
      rho_uy(i,j,k) = rho(i,j,k) * uy(i,j,k)
      rho_uz(i,j,k) = rho(i,j,k) * uz(i,j,k)
      ! Keep Backup of the Arrays for FFTW. 
      rho_dum(i,j,k) = rho(i,j,k)
      ux_dum(i,j,k) = ux(i,j,k)
      uy_dum(i,j,k) = uy(i,j,k)
      uz_dum(i,j,k) = uz(i,j,k)
      rho_ux_dum(i,j,k) = rho_ux(i,j,k)
      rho_uy_dum(i,j,k) = rho_uy(i,j,k)
      rho_uz_dum(i,j,k) = rho_uz(i,j,k)
      P_dum(i,j,k) = P(i,j,k)
      E_dum(i,j,k) = E(i,j,k) 
      Bx_dum(i,j,k) = Bx(i,j,k)
      By_dum(i,j,k) = By(i,j,k)
      Bz_dum(i,j,k) = Bz(i,j,k)
      ! Write Initial Density and Velocity Distribution in File.
      write(15,*) x(i),y(j),z(k),rho(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k),P(i,j,k),E(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    enddo
  end do
end do

 close(15)

!===================== INITIAL TIME DATA ===============================================

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_dum, rho_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_ux_dum, rho_ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_uy_dum, rho_uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_uz_dum, rho_uz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, P_dum, P_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, E_dum, E_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Bx_dum, Bx_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, By_dum, By_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Bz_dum, Bz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_destroy_plan_ (plan_forward)

! Evaluate Vorticity in Spectral Space.

  call OMEGA (Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,ux_k,uy_k,uz_k,omega_x_k,omega_y_k,omega_z_k,&
              omega_x,omega_y,omega_z,omega2) 

!======================= MAIN PROGRAM =====================================================
t1 = omp_get_wtime()

! Time Loop Starts.
do time = time_min,time_max,dt

t = nint(time/dt) - dint(time_min/dt)

! Online Check of the Progress.
!if (mod(t,int(1.0/dt)) == 0) then
  !print*, t/int(1.0/dt)
!endif

!======================================================================================

! Non-Linear Term Evaluator.
  call DERIVE (Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
               d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
               d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
               d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)
  
! Time Solvers.
! Adams-Bashforth
  call ADAMS_BASHFORTH (Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
           rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new, &
           d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
           d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
           d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)
  
!======================================================================================

! Reset the variables for next time iteration.

  call RESET (Nx,Ny,Nz,Nh,d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old,&
              d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old,&
              d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new,&
              d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new,&
              rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k,&
              rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new,&
              rho_k_dum,rho_ux_k_dum,rho_uy_k_dum,rho_uz_k_dum,E_k_dum,Bx_k_dum,By_k_dum,Bz_k_dum)

! Evaluate Divergence of B in Spectral Space.

  call DIVERGENCE_B (Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,Bx_k,By_k,Bz_k,div_B)

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
  
  call U_B2_P (Nx,Ny,Nz,Nh,spheat,rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz,B2,&
               ux,uy,uz,ux_dum,uy_dum,uz_dum,P)
  
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_destroy_plan_ (plan_forward)
  call dfftw_destroy_plan_ (plan_backward)

! Evaluate Vorticity in Spectral Space.

  call OMEGA (Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,ux_k,uy_k,uz_k,omega_x_k,omega_y_k,omega_z_k,&
              omega_x,omega_y,omega_z,omega2) 

! Evaluate Energy Spectra at a Specified Time.
if (t >= int((time_max-time_max/10.0d0)/dt)) then

  !$OMP PARALLEL SHARED(ux_k,uy_k,uz_k) PRIVATE(i,j,k) 
  !$OMP DO REDUCTION (+:Ek,Bk) 

  do i = 1,Nx/2+1
    do j = 1,Ny
      do k = 1,Nz
        Ek(i,j,k) = Ek(i,j,k) + sqrt(abs(ux_k(i,j,k))**2 + abs(uy_k(i,j,k))**2 + abs(uz_k(i,j,k))**2)&
                    /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
        Bk(i,j,k) = Bk(i,j,k) + sqrt(abs(Bx_k(i,j,k))**2 + abs(By_k(i,j,k))**2 + abs(Bz_k(i,j,k))**2)&
                    /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      enddo
    enddo
  enddo  

  !$OMP END PARALLEL

endif 
 
  call E_Casimir_Div_B (Nx,Ny,Nz,Nh,pi,mu_0,omega_x,omega_y,omega_z,omega2,rho,ux,uy,uz,P,E,Bx,By,Bz,div_B,&
                        Pressure,Energy,y_Energy,B_Field,C0,C1,C2,C3,div_B_Tot)
  
!do i = 1,Nx
!  do j = 1,Ny
!    do k = 1,Nz
!      ! Write Grid Data in State Files.  
!        if (mod(float(t),100.0) == 0.0) then
!        write(t+100,*) x(i),y(j),z(k),omega_x(i,j,k),omega_y(i,j,k),omega_z(i,j,k), &
!                        rho(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k), &
!                        P(i,j,k),E(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
!        endif
!    enddo
!  enddo
!enddo  

write(40,*) time,Pressure,Energy,y_Energy/2.0d0,B_Field,C0,C1,C2,C3,div_B_Tot
  
! close(t+100)
 
enddo ! time

t2 = omp_get_wtime()

write(5,*) "Time taken for the run =",(t2 - t1)/(60.0d0*60.0d0),"Hours"

  call E_B_Spectra (Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,time_max,dt,Ek,Bk)  
 
  close(5)
  close(40)

!====================================================================================

end program PSMHD3  

