program spectra3
! gfortran -I/soft/fftw-3.3.3/include -L/home/rupak/bin/lib spectra3.f95 -lfftw3 -lm
implicit none

include "fftw3.f"

integer ( kind = 4 ), parameter :: Nx = 64
integer ( kind = 4 ), parameter :: Ny = 64
integer ( kind = 4 ), parameter :: Nz = 64
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1
    
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

real( kind = 8 ) Lx,Ly,Lz,kx,ky,kz,Es,Bs,r,r_min,r_max,dr,G
integer ( kind = 4 ) i,j,k,t
real ( kind = 8 ) x(Nx),y(Ny),z(Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz),Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz)
real ( kind = 8 ) E(Nx,Ny,Nz),B(Nx,Ny,Nz)
complex ( kind = 8 ) Ek(Nh,Ny,Nz),Bk(Nh,Ny,Nz)
integer ( kind = 8 ) plan_forward
    
!open(10,file="Spectra.dat")

Lx = 2.0d0*pi; Ly = 2.0d0*pi; Lz = 2.0d0*pi

r_min=1.0d0
r_max=int(sqrt(dfloat(Nx**2 + Ny**2 + Nz**2)/3.0d0))
dr = 1.0d0
    
do t = 100, 6100, 1000    
do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
    read(t,*) G,G,G,ux(i,j,k),uy(i,j,k),uz(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    enddo
  enddo
enddo    

  close(t)

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
    E(i,j,k) = ux(i,j,k)*ux(i,j,k) + uy(i,j,k)*uy(i,j,k) + uz(i,j,k)*uz(i,j,k)
    B(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)
    enddo
  enddo
enddo 
   
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, E, Ek, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, B, Bk, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

do i = 1,Nh
  do j = 1,Ny
    do k = 1,Nz
    Ek(i,j,k) = Ek(i,j,k)/dfloat(Nx*Ny*Nz)
    Bk(i,j,k) = Bk(i,j,k)/dfloat(Nx*Ny*Nz)
    enddo
  enddo
enddo 

do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    kz = 2.0d0*pi*dfloat(k-1)/Lz
      if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
      !if (kx .ge. 2*Nx/3 .or. ky .ge. 2*Ny/3 .or. kz .ge. 2*Nz/3) then
      !Ek(i,j,k) = 0.0d0
      !Bk(i,j,k) = 0.0d0
      endif
    !write(35,*) i,j,k,dsqrt(dfloat((i-1)**2)+dfloat((j-1)**2)+dfloat((k-1)**2)),&
    !            abs(Ek(i,j,k)), abs(Bk(i,j,k))
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
      !if (kx .ge. 2*Nx/3 .or. ky .ge. 2*Ny/3 .or. kz .ge. 2*Nz/3) then
      !Ek(i,j,k) = 0.0d0
      !Bk(i,j,k) = 0.0d0
      endif
    !write(35,*) i,j,k-Nz,dsqrt(dfloat((i-1)**2)+dfloat((j-1)**2)+dfloat(((k-1)-Nz)**2)),&
    !            abs(Ek(i,j,k)), abs(Bk(i,j,k))
    enddo
  enddo  
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    kz = 2.0d0*pi*dfloat(k-1)/Lz
      if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
      !if (kx .ge. 2*Nx/3 .or. ky .ge. 2*Ny/3 .or. kz .ge. 2*Nz/3) then
      !Ek(i,j,k) = 0.0d0
      !Bk(i,j,k) = 0.0d0
      endif
    !write(35,*) i,j-Ny,k,dsqrt(dfloat((i-1)**2)+dfloat(((j-1)-Ny)**2)+dfloat((k-1)**2)),&
    !            abs(Ek(i,j,k)), abs(Bk(i,j,k))
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
      !if (kx .ge. 2*Nx/3 .or. ky .ge. 2*Ny/3 .or. kz .ge. 2*Nz/3) then
      !Ek(i,j,k) = 0.0d0
      !Bk(i,j,k) = 0.0d0
      endif
    !write(35,*) i,j-Ny,k-Nz,dsqrt(dfloat((i-1)**2)+dfloat(((j-1)-Ny)**2)+dfloat(((k-1)-Nz)**2)),&
    !            abs(Ek(i,j,k)), abs(Bk(i,j,k))
    enddo
  enddo  
enddo 

do r = r_min,r_max,dr
Es = 0.0d0
Bs = 0.0d0
  do i=1,Nh
    do j=1,Ny/2
      do k=1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
        if(dsqrt(kx**2+ky**2+kz**2)<=r .and. dsqrt(kx**2+ky**2+kz**2)>(r-dr)) then
        Es = Es + abs (Ek(i,j,k))
        Bs = Bs + abs (Bk(i,j,k))
        endif
      enddo
      do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
        if(dsqrt(kx**2+ky**2+kz**2)<=r .and. dsqrt(kx**2+ky**2+kz**2)>(r-dr)) then
        Es = Es + abs (Ek(i,j,k))
        Bs = Bs + abs (Bk(i,j,k))
        endif
      enddo 
    enddo
    do j = Ny/2+1,Ny
      do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
        if(dsqrt(kx**2+ky**2+kz**2)<=r .and. dsqrt(kx**2+ky**2+kz**2)>(r-dr)) then
        Es = Es + abs (Ek(i,j,k))
        Bs = Bs + abs (Bk(i,j,k))
        endif
      enddo 
      do k = Nz/2+1,Nz   
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
        if(dsqrt(kx**2+ky**2+kz**2)<=r .and. dsqrt(kx**2+ky**2+kz**2)>(r-dr)) then
        Es = Es + abs (Ek(i,j,k))
        Bs = Bs + abs (Bk(i,j,k))
        endif
      enddo
    enddo
  enddo
  write(t+50,*) r-1,Es,Bs

  if (r == r_min+dr) then
  write(10,*) t,Es,Bs
  elseif (r == r_min+2.0d0*dr) then
  write(20,*) t,Es,Bs
  elseif (r == r_min+3.0d0*dr) then
  write(30,*) t,Es,Bs
  elseif (r == r_min+4.0d0*dr) then
  write(40,*) t,Es,Bs
  elseif (r == r_min+5.0d0*dr) then
  write(50,*) t,Es,Bs
  endif
  
  call flush (10)
  call flush (20)
  call flush (30)
  call flush (40)
  call flush (50)
  
enddo
  close(t+50)
enddo ! t
end program spectra3
