module reconnection_metrics
  contains
    subroutine compute_current(x,y,z,bx,by,bz,nx,ny,nz,jx,jy,jz,mag_j)

      implicit none

      integer, intent(in) :: nx,ny,nz
      !f2py intent(hide) :: nx, ny, nz
      real(8), intent(in) :: x(nx), y(ny), z(nz)
      !f2py intent(in) :: x(nx), y(ny), z(nz)
      real(8), intent(in) :: bx(nx,ny,nz), by(nx,ny,nz), bz(nx,ny,nz)
      !f2py intent(in) :: bx(nx,ny,nz), by(nx,ny,nz), bz(nx,ny,nz)
      real(8), intent(out) :: jx(nx,ny,nz), jy(nx,ny,nz), jz(nx,ny,nz), mag_j(nx,ny,nz)
      !f2py intent(out) :: jx(nx,ny,nz), jy(nx,ny,nz), jz(nx,ny,nz), mag_j(nx,ny,nz)

      integer :: ix, iy, iz
      real(8) :: xx, yy, zz, hx, hy, hz

      ! We assume that the grid spacing is uniform
      hx = abs(x(2) - x(1))
      hy = abs(y(2) - y(1))
      hz = abs(z(2) - z(1))

      do ix = 2,nx-1
      xx = x(ix)
        do iy = 2,ny-1
        yy = y(iy)
          do iz = 2,nz-1
          zz = z(iz)
          jx(ix,iy,iz) = o2_central_diff(hy,bz(ix,iy+1,iz),bz(ix,iy,iz),bz(ix,iy-1,iz))&
                       - o2_central_diff(hz,by(ix,iy,iz+1),by(ix,iy,iz),by(ix,iy,iz-1))
          jy(ix,iy,iz) = o2_central_diff(hz,bx(ix,iy,iz+1),bx(ix,iy,iz),bx(ix,iy,iz-1))&
                       - o2_central_diff(hx,bz(ix+1,iy,iz),bz(ix,iy,iz),bz(ix-1,iy,iz))
          jz(ix,iy,iz) = o2_central_diff(hx,by(ix+1,iy,iz),by(ix,iy,iz),by(ix-1,iy,iz))&
                       - o2_central_diff(hy,bx(ix,iy+1,iz),bx(ix,iy,iz),bx(ix,iy-1,iz))
          mag_j(ix,iy,iz) = sqrt(jx(ix,iy,iz)**2 + jy(ix,iy,iz)**2 + jz(ix,iy,iz)**2 )
          end do
        end do
      end do
    end subroutine compute_current

    function o2_central_diff(h,fp,f,fm)
      implicit none
      real(8), intent(in) :: h, fp, f, fm
      real(8) :: o2_central_diff

      o2_central_diff = (fp - 2.0*f + fm) / ( h**2 )
      
    end function o2_central_diff
      
end module reconnection_metrics
