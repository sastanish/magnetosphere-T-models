module reconnection_metrics
  implicit none

  integer :: nx,ny,nz
  real(8), dimension(:,:,:), allocatable :: bx,by,bz
  real(8), dimension(:), allocatable :: x,y,z

contains

  subroutine read_netcdf_field_file(filename)
    ! Note, this subroutine assumes that the layout of the netcdf file
    ! as generated from the Tsynenko model.

    use netcdf

    character(*), intent(in) :: filename
    character(len=50) :: dummy
    character(len=50) :: xname, yname, zname
    character(len=50) :: bxname, byname, bzname
    integer :: nc_id

    ! read the file and fill reference id, nc_id
    call check(nf90_open(filename, nf90_nowrite, nc_id))

    ! get indices of dimensions

    ! Read size of dimensions
    call check(nf90_inquire_dimension(nc_id,1,xname,nx))
    call check(nf90_inquire_dimension(nc_id,2,yname,ny))
    call check(nf90_inquire_dimension(nc_id,3,zname,nz))
    
    print *, "(nx, ny, nz): ", nx, ny, nz

    ! Allocate data in ram
    allocate(x(nx))
    allocate(y(ny))
    allocate(z(nz))

    ! netcdf files are backwards indexed
    allocate(bx(nz,ny,nx))
    allocate(by(nz,ny,nx))
    allocate(bz(nz,ny,nx))

    ! Read all data (dimension and variables) into memory
    call check(nf90_get_var(nc_id,4,x))
    call check(nf90_get_var(nc_id,5,y))
    call check(nf90_get_var(nc_id,6,z))
    call check(nf90_get_var(nc_id,1,bx))
    call check(nf90_get_var(nc_id,2,by))
    call check(nf90_get_var(nc_id,3,bz))

    print *, shape(bx)

    print *, 'close infile'
    call check(nf90_close(nc_id))

  end subroutine read_netcdf_field_file

  subroutine write_netcdf_metric_file(filename,dat_name,datx,daty,datz,magdat)
    ! Assumes you wish to write out some calculated vector quantity and it's magnitude
    ! on the same mesh as the original magnetic field file.  

    use netcdf

    real(8), dimension(nz,ny,nx), intent(in) :: datx,daty,datz,magdat
    character(*), intent(in) :: filename
    character(*), intent(in) :: dat_name
    integer :: nc_id

    integer :: id_nc,xdim_id,ydim_id,zdim_id
    integer :: x_id,y_id,z_id,dx_id,dy_id,dz_id,dm_id
    integer, dimension(3) :: dim_ids

    call check(nf90_create(filename, NF90_CLOBBER, id_nc))

    ! setup dimensions
    call check(nf90_def_dim(id_nc, "x", nx, xdim_id))
    call check(nf90_def_dim(id_nc, "y", ny, ydim_id))
    call check(nf90_def_dim(id_nc, "z", nz, zdim_id))

    ! create variables for data
    call check(nf90_def_var(id_nc, "x", NF90_DOUBLE, xdim_id, x_id))
    call check(nf90_def_var(id_nc, "y", NF90_DOUBLE, ydim_id, y_id))
    call check(nf90_def_var(id_nc, "z", NF90_DOUBLE, zdim_id, z_id))

    dim_ids = (/ zdim_id,ydim_id,xdim_id /)

    call check(nf90_def_var(id_nc, dat_name//"x", NF90_DOUBLE, dim_ids, dx_id))
    call check(nf90_def_var(id_nc, dat_name//"y", NF90_DOUBLE, dim_ids, dy_id))
    call check(nf90_def_var(id_nc, dat_name//"z", NF90_DOUBLE, dim_ids, dz_id))
    call check(nf90_def_var(id_nc, "mag_"//dat_name, NF90_DOUBLE, dim_ids, dm_id))

    call check(nf90_enddef(id_nc))

    ! Write data

    print *, 'write data'
    call check(nf90_put_var(id_nc, x_id, x))
    call check(nf90_put_var(id_nc, y_id, y))
    call check(nf90_put_var(id_nc, z_id, z))
    call check(nf90_put_var(id_nc, dx_id, datx))
    call check(nf90_put_var(id_nc, dy_id, daty))
    call check(nf90_put_var(id_nc, dz_id, datz))
    call check(nf90_put_var(id_nc, dm_id, magdat))

    print *, 'close'
    call check(nf90_close(nc_id))

  end subroutine write_netcdf_metric_file

  subroutine compute_current(jx, jy, jz, mag_j)

    real(8), intent(inout), dimension(nz,ny,nx) :: jx, jy, jz, mag_j

    integer :: ix, iy, iz
    integer :: ixp, ixm, iyp, iym, izp, izm
    real(8) :: hx, hy, hz

    ! We assume that the grid spacing is uniform
    print *, "grid spacing"
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

    print *, "start loop"
    ! compute current via second order central diff
    do ix = 1,nx
      if (ix==1)  then; ixm=ix; else; ixm=ix-1; endif
      if (ix==nx) then; ixp=ix; else; ixp=ix+1; endif
      do iy = 1,ny
        if (iy==1)  then; iym=iy; else; iym=iy-1; endif
        if (iy==ny) then; iyp=iy; else; iyp=iy+1; endif
        do iz = 1,nz
          if (iz==1)  then; izm=iz; else; izm=iz-1; endif
          if (iz==nz) then; izp=iz; else; izp=iz+1; endif

          print *, "(ix, iy, iz): ", ix, iy, iz
          ! compute
          print *, "jx"
          jx(iz,iy,ix) = o2_central_diff(hy,bz(iz ,iyp,ix ),bz(iz,iy,ix),bz(iz ,iym,ix ))&
                       - o2_central_diff(hz,by(izp,iy ,ix ),by(iz,iy,ix),by(izm,iy ,ix ))
          print *, "jy"
          jy(iz,iy,ix) = o2_central_diff(hz,bx(izp,iy ,ix ),bx(iz,iy,ix),bx(izm,iy ,ix ))&
                       - o2_central_diff(hx,bz(iz ,iy ,ixp),bz(iz,iy,ix),bz(iz ,iy ,ixm))
          print *, "jz"
          jz(iz,iy,ix) = o2_central_diff(hx,by(iz ,iy ,ixp),by(iz,iy,ix),by(iz ,iy ,ixm))&
                       - o2_central_diff(hy,bx(iz ,iyp,ix ),bx(iz,iy,ix),bx(iz ,iym,ix ))
          print *, "magj"
          mag_j(iz,iy,ix) = sqrt(jx(iz ,iy ,ix)**2 + jy(iz ,iy ,ix)**2 + jz(iz ,iy ,ix)**2 )
        end do
      end do
    end do
    print *, "end loop"

  end subroutine compute_current

  subroutine calculate_metrics(filename)

    character(*), intent(in) :: filename
    real(8), dimension(:,:,:), allocatable :: jx, jy, jz, mag_j

    integer :: name_len

    call read_netcdf_field_file(filename)

    print *, 'allocate jx'
    allocate(jx(nz,ny,nx))
    allocate(jy(nz,ny,nx))
    allocate(jz(nz,ny,nx))
    allocate(mag_j(nz,ny,nx))

    print *, 'compute current'
    call compute_current(jx,jy,jz,mag_j)

    name_len = len(filename)
    print *, 'start writing'
    call write_netcdf_metric_file(trim(filename(1:name_len-3))//"_j.nc","j",jx,jy,jz,mag_j)

  end subroutine calculate_metrics

  subroutine check(istatus)
    use netcdf
    implicit none
    integer, intent (in) :: istatus
    if (istatus /= nf90_noerr) then
    write(*,*) trim(adjustl(nf90_strerror(istatus)))
    end if
  end subroutine check

  function o2_central_diff(h,fp,f,fm)
    real(8), intent(in) :: h, fp, f, fm
    real(8) :: o2_central_diff

    o2_central_diff = (fp - 2.0*f + fm) / ( h**2 )
    
  end function o2_central_diff
    
end module reconnection_metrics
