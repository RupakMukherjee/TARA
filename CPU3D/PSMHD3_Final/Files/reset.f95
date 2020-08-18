! Reset the variables for next time iteration.

!===================================================================================
!================== SUBROUTINE RESET ===============================================
!===================================================================================

subroutine RESET (Nx,Ny,Nz,Nh,d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old,&
                  d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old,&
                  d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new,&
                  d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new,&
                  rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k,&
                  rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new,&
                  rho_k_dum,rho_ux_k_dum,rho_uy_k_dum,rho_uz_k_dum,E_k_dum,Bx_k_dum,By_k_dum,Bz_k_dum)

implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k

complex ( kind = 8 ) rho_k(Nh,Ny,Nz),rho_ux_k(Nh,Ny,Nz),rho_uy_k(Nh,Ny,Nz),rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) E_k(Nh,Ny,Nz),Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_k_new(Nh,Ny,Nz),rho_ux_k_new(Nh,Ny,Nz),rho_uy_k_new(Nh,Ny,Nz),rho_uz_k_new(Nh,Ny,Nz)
complex ( kind = 8 ) E_k_new(Nh,Ny,Nz),Bx_k_new(Nh,Ny,Nz),By_k_new(Nh,Ny,Nz),Bz_k_new(Nh,Ny,Nz)
complex ( kind = 8 ) rho_k_dum(Nh,Ny,Nz),rho_ux_k_dum(Nh,Ny,Nz),rho_uy_k_dum(Nh,Ny,Nz),rho_uz_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) E_k_dum(Nh,Ny,Nz),Bx_k_dum(Nh,Ny,Nz),By_k_dum(Nh,Ny,Nz),Bz_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny,Nz),d_rho_ux_k_dt_old(Nh,Ny,Nz),d_rho_uy_k_dt_old(Nh,Ny,Nz),d_rho_uz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_old(Nh,Ny,Nz),d_Bx_k_dt_old(Nh,Ny,Nz),d_By_k_dt_old(Nh,Ny,Nz),d_Bz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny,Nz),d_rho_ux_k_dt_new(Nh,Ny,Nz),d_rho_uy_k_dt_new(Nh,Ny,Nz),d_rho_uz_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_new(Nh,Ny,Nz),d_Bx_k_dt_new(Nh,Ny,Nz),d_By_k_dt_new(Nh,Ny,Nz),d_Bz_k_dt_new(Nh,Ny,Nz)

!$OMP PARALLEL SHARED(d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old),&
!$OMP & SHARED(d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old),&
!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new),& 
!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new),& 
!$OMP & SHARED(rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k),&
!$OMP & SHARED(rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new),&
!$OMP & SHARED(rho_k_dum,rho_ux_k_dum,rho_uy_k_dum,rho_uz_k_dum,E_k_dum,Bx_k_dum,By_k_dum,Bz_k_dum) PRIVATE(i,j,k)
!$OMP DO 

do i = 1,Nx/2+1
  do j = 1,Ny
    do k = 1,Nz
      ! Set the Variables in Proper Format for Next Time Iteration.
      d_rho_k_dt_old(i,j,k) = d_rho_k_dt_new(i,j,k)
    
      d_rho_ux_k_dt_old(i,j,k) = d_rho_ux_k_dt_new(i,j,k)
      d_rho_uy_k_dt_old(i,j,k) = d_rho_uy_k_dt_new(i,j,k)
      d_rho_uz_k_dt_old(i,j,k) = d_rho_uz_k_dt_new(i,j,k)
    
      d_E_k_dt_old(i,j,k) = d_E_k_dt_new(i,j,k)
    
      d_Bx_k_dt_old(i,j,k) = d_Bx_k_dt_new(i,j,k)
      d_By_k_dt_old(i,j,k) = d_By_k_dt_new(i,j,k)
      d_Bz_k_dt_old(i,j,k) = d_Bz_k_dt_new(i,j,k)
    
      rho_k(i,j,k) = rho_k_new(i,j,k)
    
      rho_ux_k(i,j,k) = rho_ux_k_new(i,j,k)
      rho_uy_k(i,j,k) = rho_uy_k_new(i,j,k)
      rho_uz_k(i,j,k) = rho_uz_k_new(i,j,k)
  
      E_k(i,j,k) = E_k_new(i,j,k)
    
      Bx_k(i,j,k) = Bx_k_new(i,j,k)
      By_k(i,j,k) = By_k_new(i,j,k)
      Bz_k(i,j,k) = Bz_k_new(i,j,k)
    
      ! Keep Backup of the Arrays for FFTW.
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
      
return

end subroutine RESET

!====================================================================================

