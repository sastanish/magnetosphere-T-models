module compute

  implicit none

contains

  subroutine run_TA16(parmod,ps,x,y,z,Bx,By,Bz,nx,ny,nz)

    use TA16, only : RBF_MODEL_2016
    use omp_lib

    integer :: nx, ny, nz
    integer :: i, j, k
    real(8) :: xx,yy,zz,bbx,bby,bbz !dummy variables
    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
    real(8), intent(in) :: parmod(10), ps

    !$OMP PARALLEL SHARED(Bx,By,Bz,nx,ny,nz,parmod,ps) PRIVATE(bbx,bby,bbz,xx,yy,zz,i,j,k)
    !$ print *, OMP_GET_MAX_THREADS()
    !$OMP DO COLLAPSE(3)
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          xx = x(i)
          yy = y(j)
          zz = z(k)
          call RBF_MODEL_2016(0,parmod,ps,xx,yy,zz,bbx,bby,bbz)
          Bx(i,j,k) = bbx
          By(i,j,k) = bby
          Bz(i,j,k) = bbz
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
  end subroutine run_TA16

  subroutine run_TS05(parmod,ps,x,y,z,Bx,By,Bz,nx,ny,nz)

    use TS05, only : T04_s
    use omp_lib

    integer :: nx, ny, nz
    integer :: i, j, k
    real(8) :: xx,yy,zz,bbx,bby,bbz !dummy variables
    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
    real(8), intent(in) :: parmod(10), ps

    !$OMP PARALLEL SHARED(Bx,By,Bz,nx,ny,nz,parmod,ps) PRIVATE(bbx,bby,bbz,xx,yy,zz,i,j,k)
    !$ print *, OMP_GET_MAX_THREADS()
    !$OMP DO COLLAPSE(3)
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          xx = x(i)
          yy = y(j)
          zz = z(k)
          call T04_s(0,parmod,ps,xx,yy,zz,bbx,bby,bbz)
          Bx(i,j,k) = bbx
          By(i,j,k) = bby
          Bz(i,j,k) = bbz
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
  end subroutine run_TS05

  subroutine run_igrf_dipole(iyear, iday, ihour, imin, isec, vx_gse, vy_gse, vz_gse, &
                             x_gsw, y_gsw, z_gsw, hx_gsw, hy_gsw, hz_gsw, nx, ny, nz)

    use geopack, only : RECALC_08, IGRF_GSW_08

    integer :: nx, ny, nz
    integer :: i, j, k
    real(8) :: xx,yy,zz,hhx,hhy,hhz !dummy variables
    integer, intent(in) :: iyear, iday, ihour, imin, isec
    real(8), intent(in) :: vx_gse, vy_gse, vz_gse
    real(8), intent(in) :: x_gsw(nx), y_gsw(ny), z_gsw(nz)
    real(8), intent(out) :: hx_gsw(nx,ny,nz), hy_gsw(nx,ny,nz), hz_gsw(nx,ny,nz)
    
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
    real(8) :: xx,yy,zz,hhx,hhy,hhz !dummy variables
    integer, intent(in) :: iyear, iday, ihour, imin, isec
    real(8), intent(in) :: vx_gse, vy_gse, vz_gse
    real(8), intent(in) :: x_gsw(nx), y_gsw(ny), z_gsw(nz)
    real(8), intent(out) :: hx_gsw(nx,ny,nz), hy_gsw(nx,ny,nz), hz_gsw(nx,ny,nz)

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

  subroutine metrics(x,y,z,bx,by,bz,Mout,nx,ny,nz)

    use reconnectionMetrics, only: compute_all_metrics

    integer :: nx, ny, nz

    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), dimension(nx,ny,nz), intent(in) :: bx, by, bz
    real(8), intent(out) :: Mout(23,nx,ny,nz)

    call compute_all_metrics(x,y,z,bx,by,bz,Mout)

  end subroutine metrics

  subroutine total_rate(x,y,z,bx,by,bz,Mout,nx,ny,nz)

    use reconnectionMetrics, only: compute_total_rate

    integer :: nx, ny, nz

    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), dimension(nx,ny,nz), intent(in) :: bx, by, bz
    real(8), intent(out) :: Mout(nx,ny,nz)

    call compute_total_rate(x,y,z,bx,by,bz,Mout)

  end subroutine total_rate

end module compute
