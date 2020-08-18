! Provide Initial Condition.

!===================================================================================
!================== SUBROUTINE Initial Condition ===================================
!===================================================================================

subroutine Initial_Condition (Nx,Ny,Nz,Nh,pi,Lx,Ly,Lz,dx,dy,dz,&
                              spheat,rho0,U0,P0,P,B0,rho,x,y,z,ux,uy,uz,Bx,By,Bz)

implicit none

integer ( kind = 4 ) Nx,Ny,Nz,Nh,i,j,k
real ( kind = 8 ) pi,Lx,Ly,Lz,dx,dy,dz,spheat,rho0,U0,P0,B0

real ( kind = 8 ) x(Nx),y(Ny),z(Nz),rho(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz)
real ( kind = 8 ) P(Nx,Ny,Nz),E(Nx,Ny,Nz),Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz)

integer,parameter :: seed = 99999999
call srand(seed)

! Grid Generation.
do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0d0+real(j-1)*dy
    do k = 1, Nz
      z(k)=0.0d0+real(k-1)*dz
      ! Initial Density Distribution.
      rho(i,j,k) = rho0
      ! Initial Velocity Distribution.
      ux(i,j,k) = U0 * rand() 
      uy(i,j,k) = U0 * rand() 
      uz(i,j,k) = U0 * rand() 
      ! Initial Pressure Distribution.
      P(i,j,k) = P0
      ! Initial Magnetic Field Profile.
      Bx(i,j,k) = B0 * rand() 
      By(i,j,k) = B0 * rand() 
      Bz(i,j,k) = B0 * rand() 
      ! Initial Energy Distribution.
      E(i,j,k) = P(i,j,k)/(spheat -1.0d0) + 0.50d0*rho(i,j,k) * (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)!&
      ! + 0.50d0*(Bx(i,j)**2 + By(i,j)**2 + Bz(i,j)**2)
    enddo
  end do
end do

return

end subroutine Initial_Condition

!====================================================================================

