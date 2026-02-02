module inputOutput

  implicit none
  private
  public :: save_field_to_netcdf, save_rate_to_netcdf, load_field_from_netcdf, load_rate_from_netcdf

contains

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

end module inputOutput
