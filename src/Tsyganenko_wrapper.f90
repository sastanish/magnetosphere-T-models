module models
contains

  subroutine run_TA16(parmod,ps,x,y,z,Bx,By,Bz,nx,ny,nz)

    use TA16, only : RBF_MODEL_2016

    integer :: nx, ny, nz
    integer :: i, j, k
    real(8) :: xx,yy,zz,bbx,bby,bbz !dummy variables
    real(8) :: x(nx), y(ny), z(nz)
    real(8) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
    real(8) :: parmod(10), ps

    bbx = 0.0
    bby = 0.0
    bbz = 0.0
    do i = 1,nx
      xx = x(i)
      do j = 1,ny
        yy = y(j)
        do k = 1,nz
          zz = z(k)
          call RBF_MODEL_2016(0,parmod,ps,xx,yy,zz,bbx,bby,bbz)
          Bx(i,j,k) = bbx
          By(i,j,k) = bby
          Bz(i,j,k) = bbz
        end do
      end do
    end do
    
  end subroutine run_TA16

  subroutine run_TS05(parmod,ps,x,y,z,Bx,By,Bz,nx,ny,nz)

    use TS04c, only : T04_s

    integer :: nx, ny, nz
    integer :: i, j, k
    real(8) :: xx,yy,zz,bbx,bby,bbz !dummy variables
    real(8) :: x(nx), y(ny), z(nz)
    real(8) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
    real(8) :: parmod(10), ps

    bbx = 0.0
    bby = 0.0
    bbz = 0.0
    do i = 1,nx
      xx = x(i)
      do j = 1,ny
        yy = y(j)
        do k = 1,nz
          zz = z(k)
          call T04_s(0,parmod,ps,xx,yy,zz,bbx,bby,bbz)
          Bx(i,j,k) = bbx
          By(i,j,k) = bby
          Bz(i,j,k) = bbz
        end do
      end do
    end do
    
  end subroutine run_TS05

  subroutine run_igrf_dipole(iyear, iday, ihour, imin, isec, vx_gse, vy_gse, vz_gse, &
                             x_gsw, y_gsw, z_gsw, hx_gsw, hy_gsw, hz_gsw, nx, ny, nz)

    use geopack, only : RECALC_08, IGRF_GSW_08

    integer :: nx, ny, nz
    integer :: i, j, k
    integer :: iyear, iday, ihour, imin, isec
    real(8) :: vx_gse, vy_gse, vz_gse
    real(8) :: xx,yy,zz,hhx,hhy,hhz !dummy variables
    real(8) :: x_gsw(nx), y_gsw(ny), z_gsw(nz)
    real(8) :: hx_gsw(nx,ny,nz), hy_gsw(nx,ny,nz), hz_gsw(nx,ny,nz)
    
    hhx = 0.0
    hhy = 0.0
    hhz = 0.0
    call RECALC_08(iyear, iday, ihour, imin, isec, vx_gse, vy_gse, vz_gse)
    do i = 1,nx
      xx = x_gsw(i)
      do j = 1,ny
        yy = y_gsw(j)
        do k = 1,nz
          zz = z_gsw(k)
          call IGRF_GSW_08(xx,yy,zz,hhx,hhy,hhz)
          hx_gsw(i,j,k) = hhx
          hy_gsw(i,j,k) = hhy
          hz_gsw(i,j,k) = hhz
        end do
      end do
    end do
    
  end subroutine run_igrf_dipole

  subroutine run_dipole(iyear, iday, ihour, imin, isec, vx_gse, vy_gse, vz_gse, &
                             x_gsw, y_gsw, z_gsw, hx_gsw, hy_gsw, hz_gsw, nx, ny, nz)

    use geopack, only : DIP_08, RECALC_08

    integer :: nx, ny, nz
    integer :: i, j, k
    integer :: iyear, iday, ihour, imin, isec
    real(8) :: vx_gse, vy_gse, vz_gse
    real(8) :: xx,yy,zz,hhx,hhy,hhz !dummy variables
    real(8) :: x_gsw(nx), y_gsw(ny), z_gsw(nz)
    real(8) :: hx_gsw(nx,ny,nz), hy_gsw(nx,ny,nz), hz_gsw(nx,ny,nz)

    hhx = 0.0
    hhy = 0.0
    hhz = 0.0
    call RECALC_08(iyear, iday, ihour, imin, isec, vx_gse, vy_gse, vz_gse)
    do i = 1,nx
      xx = x_gsw(i)
      do j = 1,ny
        yy = y_gsw(j)
        do k = 1,nz
          zz = z_gsw(k)
          call DIP_08(xx,yy,zz,hhx,hhy,hhz)
          hx_gsw(i,j,k) = hhx
          hy_gsw(i,j,k) = hhy
          hz_gsw(i,j,k) = hhz
        end do
      end do
    end do
    
  end subroutine run_dipole


end module
