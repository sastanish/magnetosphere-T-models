module NATsyganenko_magnetosphere
contains

  subroutine run_TS05(parmod,ps,x,y,z,Bx,By,Bz,nx,ny,nz)
    integer :: nx, ny, nz
    integer :: i, j, k
    real :: xx,yy,zz,bbx,bby,bbz !dummy variables
    real :: x(nx), y(ny), z(nz)
    real :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
    real :: parmod(10), ps

    use T04_s

    bbx = 0.0
    bby = 0.0
    bbz = 0.0
    do i 1,nx
      xx = x(i)
      do j 1,ny
        yy = y(j)
        do k 1,nz
          zz = z(k)
          call T04_s(0,parmod,ps,xx,yy,zz,bbx,bby,bbz)
          Bx(i,j,k) = bbx
          By(i,j,k) = bby
          Bz(i,j,k) = bbz
        end do
      end do
    end do
    
  end subroutine run_TS05

  subroutine run_dipole()
  end subroutine run_dipole

  subroutine run_igrf_dipole()
  end subroutine run_igrf_dipole

end module
