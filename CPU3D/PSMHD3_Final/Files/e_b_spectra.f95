! Evaluate Energy and Magnetic Energy Spectra.

!===================================================================================
!================== SUBROUTINE E - B - Spectra =====================================
!===================================================================================

subroutine E_B_Spectra (Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,time_max,dt,Ek,Bk)

implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,Lx,Ly,Lz,kx,ky,kz,time_max,dt

complex ( kind = 8 ) Ek(Nh,Ny,Nz),Bk(Nh,Ny,Nz)

do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    kz = 2.0d0*pi*dfloat(k-1)/Lz
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
      Ek(i,j,k) = 0.0d0
      Bk(i,j,k) = 0.0d0
      endif
    write(35,*) i,j,k,sqrt(float((i-1)**2)+float((j-1)**2)+float((k-1)**2)),&
              abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
      Ek(i,j,k) = 0.0d0
      Bk(i,j,k) = 0.0d0
      endif
    write(35,*) i,j,k-Nz,sqrt(float((i-1)**2)+float((j-1)**2)+float(((k-1)-Nz)**2)),&
              abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
  enddo  
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    kz = 2.0d0*pi*dfloat(k-1)/Lz
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
      Ek(i,j,k) = 0.0d0
      Bk(i,j,k) = 0.0d0
      endif
    write(35,*) i,j-Ny,k,sqrt(float((i-1)**2)+float(((j-1)-Ny)**2)+float((k-1)**2)),&
              abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 1) then!kx >= Nx/3 .and. ky >= Ny/3) then
      Ek(i,j,k) = 0.0d0
      Bk(i,j,k) = 0.0d0
      endif
    write(35,*) i,j-Ny,k-Nz,sqrt(float((i-1)**2)+float(((j-1)-Ny)**2)+float(((k-1)-Nz)**2)),&
              abs(Ek(i,j,k))/(time_max/(dt*10.0d0)), abs(Bk(i,j,k))/(time_max/(dt*10.0d0))
    enddo
  enddo  
enddo 

 close(35)

return

end subroutine E_B_Spectra 

!====================================================================================

