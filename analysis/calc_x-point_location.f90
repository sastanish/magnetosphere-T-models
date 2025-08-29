program main

  implicit none

  !Data vars
  real(8), dimension(:,:,:), allocatable :: Bx, By, Bz, rate
  real(8), dimension(:,:), allocatable :: pressure
  real(8), dimension(:), allocatable :: x,y,z,xr,yr,zr,out_rate,out_pos
  integer :: year, day, hour, minn

  !Indices
  integer :: f,i,j,xind,yind,zind,max_x,max_z

  !File handling
  integer :: fileind, start_ind, end_ind, outfile, omnifile
  character(4) :: str_ind, start_str, end_str
  character(50) :: dir

  call GET_COMMAND_ARGUMENT(1,dir)
  call GET_COMMAND_ARGUMENT(2,start_str)
  call GET_COMMAND_ARGUMENT(3,end_str)
  read(start_str,'(I5)') start_ind
  read(end_str,'(I5)') end_ind

  allocate(out_rate(start_ind:end_ind))
  allocate(out_pos(start_ind:end_ind))

  do f = start_ind,end_ind

    write( str_ind, '(I4)' ) f
    print *, 'reading: '//trim(dir)//'*_'//trim(adjustl(str_ind))//'.nc'

    call load_rate_from_netcdf(trim(dir)//'rate_'//trim(adjustl(str_ind))//'.nc',xr,yr,zr,rate)
    deallocate(xr)
    deallocate(yr)
    deallocate(zr)
    call load_field_from_netcdf(trim(dir)//'output_'//trim(adjustl(str_ind))//'.nc',x,y,z,bx,by,bz)

    ! find index of y=0 or close to.
    yind = 8
    !do i=1,size(y)
    !  if (abs(y(i)) .le. 1e-3) yind=i
    !end do

    allocate(pressure(size(x),size(z)))
    do j=1,size(z)
    do i=1,size(x)
      pressure(i,j) = sqrt(bx(i,yind,j)**2 + by(i,yind,j)**2 + bz(i,yind,j)**2) 
    end do
    end do

    call closest_critical_point(x,z,pressure,xind,zind)

    out_rate(f) = rate(xind,yind,zind)
    out_pos(f)  = sqrt(x(xind)**2 + y(yind)**2 + z(zind)**2)

    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(bx)
    deallocate(by)
    deallocate(bz)
    deallocate(rate)
    deallocate(pressure)
  end do

  open(newunit=outfile, file=trim(dir)//"x-point_location.lst", status="new")
  open(newunit=omnifile, file=trim(dir)//"input_data.lst", status="old")
  write(outfile,"(A)") "X-point locations."
  write(outfile,"(A)") "FORMAT:"
  write(outfile,"(A)") "YEAR    DAY    HOUR    MIN    RATE    DISTANCE"

  ! skip omni file to correct line
  do i = 1,start_ind
    read(omnifile,*)
  end do
  do i = start_ind,end_ind
    read(omnifile,*) year, day, hour, minn
    write(outfile,'(I4, 4x, I3, 4x, I4, 4x, I3, 4x, f8.2, 4x, f6.3)') year, day, hour, minn, out_rate(i), out_pos(i)
  end do
  close(outfile)
  close(omnifile)

  deallocate(out_rate)
  deallocate(out_pos)

contains

subroutine closest_critical_point(x,y,array,out_xind,out_yind)

  real(8), intent(in) :: array(:,:), x(:), y(:)
  integer, intent(out) :: out_xind, out_yind

  integer :: i,j
  real(8) :: dxp,dxm,dyp,dym,pos,out_pos

  out_pos = 50
  out_xind = 1
  out_yind = 1

  do j = 2,size(array,2)-1
  do i = 2,size(array,1)-1
    dxp = array(i+1,j)-array(i,j)
    dxm = array(i,j)-array(i-1,j)
    dyp = array(i,j+1)-array(i,j)
    dym = array(i,j)-array(i,j-1)
    if ((dxp*dxm <= 0) .and. (dyp*dym <= 0)) then ! ie. not / / or \ \ => same signed derivatives
      pos = sqrt(x(i)**2 + y(j)**2)
      if ((pos < out_pos) .and. (pos > 2)) then
        out_pos = pos
        out_xind = i
        out_yind = j
      end if
    end if
  end do
  end do

end subroutine closest_critical_point

subroutine save_field_to_netcdf(x,y,z,bx,by,bz,filename)

  use netcdf

  real(8), intent(in), dimension(:) :: x, y, z
  real(8), intent(in), dimension(:,:,:) :: bx, by, bz
  character(*), intent(in) :: filename

  integer :: i
  integer :: nx, ny, nz

  integer :: id_nc,xdim_id,ydim_id,zdim_id
  integer :: x_id,y_id,z_id
  integer :: id_bx, id_by, id_bz

  integer, dimension(3) :: dim_ids

  call check(nf90_create(filename, NF90_NETCDF4, id_nc))

  ! setup dimensions
  call check(nf90_def_dim(id_nc, "x", size(x), xdim_id))
  call check(nf90_def_dim(id_nc, "y", size(y), ydim_id))
  call check(nf90_def_dim(id_nc, "z", size(z), zdim_id))

  ! create variables for data
  call check(nf90_def_var(id_nc, "x", NF90_DOUBLE, xdim_id, x_id))
  call check(nf90_def_var(id_nc, "y", NF90_DOUBLE, ydim_id, y_id))
  call check(nf90_def_var(id_nc, "z", NF90_DOUBLE, zdim_id, z_id))

  dim_ids = (/ xdim_id,ydim_id,zdim_id /)

  call check(nf90_def_var(id_nc, "bx", NF90_DOUBLE, dim_ids, id_bx , deflate_level=7))
  call check(nf90_def_var(id_nc, "by", NF90_DOUBLE, dim_ids, id_by , deflate_level=7))
  call check(nf90_def_var(id_nc, "bz", NF90_DOUBLE, dim_ids, id_bz , deflate_level=7))

  call check(nf90_enddef(id_nc))

  call check(nf90_put_var(id_nc, x_id, x))
  call check(nf90_put_var(id_nc, y_id, y))
  call check(nf90_put_var(id_nc, z_id, z))

  call check(nf90_put_var(id_nc, id_bx, bx))
  call check(nf90_put_var(id_nc, id_by, by))
  call check(nf90_put_var(id_nc, id_bz, bz))
  
  call check(nf90_close(id_nc))

end subroutine save_field_to_netcdf

subroutine save_rate_to_netcdf(x,y,z,rate,filename)

  use netcdf

  real(8), intent(in), dimension(:) :: x, y, z
  real(8), intent(in), dimension(:,:,:) :: rate
  character(*), intent(in) :: filename

  integer :: i
  integer :: nx, ny, nz

  integer :: id_nc,xdim_id,ydim_id,zdim_id
  integer :: x_id,y_id,z_id
  integer :: id_rate

  integer, dimension(3) :: dim_ids

  call check(nf90_create(filename, NF90_NETCDF4, id_nc))

  ! setup dimensions
  call check(nf90_def_dim(id_nc, "x", size(x), xdim_id))
  call check(nf90_def_dim(id_nc, "y", size(y), ydim_id))
  call check(nf90_def_dim(id_nc, "z", size(z), zdim_id))

  ! create variables for data
  call check(nf90_def_var(id_nc, "x", NF90_DOUBLE, xdim_id, x_id))
  call check(nf90_def_var(id_nc, "y", NF90_DOUBLE, ydim_id, y_id))
  call check(nf90_def_var(id_nc, "z", NF90_DOUBLE, zdim_id, z_id))

  dim_ids = (/ xdim_id,ydim_id,zdim_id /)

  call check(nf90_def_var(id_nc, "rate", NF90_DOUBLE, dim_ids, id_rate , deflate_level=7))

  call check(nf90_enddef(id_nc))

  call check(nf90_put_var(id_nc, x_id, x))
  call check(nf90_put_var(id_nc, y_id, y))
  call check(nf90_put_var(id_nc, z_id, z))

  call check(nf90_put_var(id_nc, id_rate, rate))
  
  call check(nf90_close(id_nc))

end subroutine save_rate_to_netcdf

subroutine load_field_from_netcdf(filename,x,y,z,Bx,By,Bz)

  use netcdf

  character(*), intent(in) :: filename
  real(8), allocatable, dimension(:,:,:), intent(inout) :: bx,by,bz
  real(8), allocatable, dimension(:), intent(inout) :: x,y,z

  character(len=50) :: dummy
  character(len=50) :: xname, yname, zname
  character(len=50) :: bxname, byname, bzname
  integer :: nc_id
  integer :: nx, ny, nz

  ! read the file and fill reference id, nc_id
  call check(nf90_open(trim(filename), nf90_nowrite, nc_id))

  ! Read size of dimensions
  call check(nf90_inquire_dimension(nc_id,1,xname,nx))
  call check(nf90_inquire_dimension(nc_id,2,yname,ny))
  call check(nf90_inquire_dimension(nc_id,3,zname,nz))
  
  ! Allocate data in ram
  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))

  allocate(Bx(nx,ny,nz))
  allocate(By(nx,ny,nz))
  allocate(Bz(nx,ny,nz))

  ! Read all data (dimension and variables) into memory
  call check(nf90_get_var(nc_id,1,x))
  call check(nf90_get_var(nc_id,2,y))
  call check(nf90_get_var(nc_id,3,z))
  call check(nf90_get_var(nc_id,4,Bx))
  call check(nf90_get_var(nc_id,5,By))
  call check(nf90_get_var(nc_id,6,Bz))

  call check(nf90_close(nc_id))

end subroutine load_field_from_netcdf

subroutine load_rate_from_netcdf(filename,x,y,z,rate)

  use netcdf

  character(*), intent(in) :: filename
  real(8), allocatable, dimension(:,:,:), intent(inout) :: rate
  real(8), allocatable, dimension(:), intent(inout) :: x,y,z

  character(len=50) :: dummy
  character(len=50) :: xname, yname, zname
  character(len=50) :: bxname, byname, bzname
  integer :: nc_id
  integer :: nx, ny, nz

  ! read the file and fill reference id, nc_id
  call check(nf90_open(trim(filename), nf90_nowrite, nc_id))

  ! Read size of dimensions
  call check(nf90_inquire_dimension(nc_id,1,xname,nx))
  call check(nf90_inquire_dimension(nc_id,2,yname,ny))
  call check(nf90_inquire_dimension(nc_id,3,zname,nz))
  
  ! Allocate data in ram
  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))

  allocate(rate(nx,ny,nz))

  ! Read all data (dimension and variables) into memory
  call check(nf90_get_var(nc_id,1,x))
  call check(nf90_get_var(nc_id,2,y))
  call check(nf90_get_var(nc_id,3,z))
  call check(nf90_get_var(nc_id,4,rate))

  call check(nf90_close(nc_id))

end subroutine load_rate_from_netcdf

subroutine check(istatus)
  use netcdf
  implicit none
  integer, intent (in) :: istatus
  if (istatus /= nf90_noerr) then
  write(*,*) trim(adjustl(nf90_strerror(istatus)))
  end if
end subroutine check

end program main

