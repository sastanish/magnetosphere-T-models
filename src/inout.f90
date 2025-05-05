module inputOutput

  implicit none

contains

  subroutine save_to_netcdf(x,y,z,bx,by,bz,reconMetrics,outputFlags,filename)

    use netcdf

    real(8), intent(in), dimension(:,:,:) :: x, y, z
    real(8), intent(in), dimension(:,:,:) :: bx, by, bz
    real(8), intent(in), dimension(:,:,:,:) :: reconMetrics
    logical, intent(in), dimension(:) :: outputFlags

    integer :: i

    integer :: id_nc,xdim_id,ydim_id,zdim_id
    integer :: x_id,y_id,z_id
    integer :: id_bx, id_by, id_bz, id_b_mag 
    integer :: id_t1_x, id_t1_y, id_t1_z, id_t1_mag 
    integer :: id_t2_x, id_t2_y, id_t2_z, id_t2_mag 
    integer :: id_t3_x, id_t3_y, id_t3_z, id_t3_mag 

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

    do i = 1,size(outputFlags)
      if (outputFlags(i) == .True.) then
        if (i == 1) then
          call check(nf90_def_var(id_nc, "bx", NF90_DOUBLE, dim_ids, id_bx , deflate_level=7))
          call check(nf90_def_var(id_nc, "by", NF90_DOUBLE, dim_ids, id_by , deflate_level=7))
          call check(nf90_def_var(id_nc, "bz", NF90_DOUBLE, dim_ids, id_bz , deflate_level=7))
        end if
        if (i == 2) call check(nf90_def_var(id_nc, "b_mag" , NF90_DOUBLE, dim_ids, id_b_mag , deflate_level=7))
        if (i == 3) then
          call check(nf90_def_var(id_nc, "c2_t1_x", NF90_DOUBLE, dim_ids, id_t1_x , deflate_level=7))
          call check(nf90_def_var(id_nc, "c2_t1_y", NF90_DOUBLE, dim_ids, id_t1_y , deflate_level=7))
          call check(nf90_def_var(id_nc, "c2_t1_z", NF90_DOUBLE, dim_ids, id_t1_z , deflate_level=7))
        end if
        if (i == 4) call check(nf90_def_var(id_nc, "c2_t1_mag" , NF90_DOUBLE, dim_ids, id_t1_mag , deflate_level=7))
        if (i == 5) then
          call check(nf90_def_var(id_nc, "c2_t2_x", NF90_DOUBLE, dim_ids, id_t2_x , deflate_level=7))
          call check(nf90_def_var(id_nc, "c2_t2_y", NF90_DOUBLE, dim_ids, id_t2_y , deflate_level=7))
          call check(nf90_def_var(id_nc, "c2_t2_z", NF90_DOUBLE, dim_ids, id_t2_z , deflate_level=7))
        end if
        if (i == 6) call check(nf90_def_var(id_nc, "c2_t2_mag" , NF90_DOUBLE, dim_ids, id_t2_mag , deflate_level=7))
        if (i == 7) then
          call check(nf90_def_var(id_nc, "c2_t3_x", NF90_DOUBLE, dim_ids, id_t3_x , deflate_level=7))
          call check(nf90_def_var(id_nc, "c2_t3_y", NF90_DOUBLE, dim_ids, id_t3_y , deflate_level=7))
          call check(nf90_def_var(id_nc, "c2_t3_z", NF90_DOUBLE, dim_ids, id_t3_z , deflate_level=7))
        end if
        if (i == 8) call check(nf90_def_var(id_nc, "c2_t3_mag" , NF90_DOUBLE, dim_ids, id_t3_mag , deflate_level=7))
      end if 

    call check(nf90_enddef(id_nc))

    call check(nf90_put_var(id_nc, x_id, x))
    call check(nf90_put_var(id_nc, y_id, y))
    call check(nf90_put_var(id_nc, z_id, z))

    ! Write data
    do i = 1,size(outputFlags)
      if (outputFlags(i) == .True.) then
        if (i == 1) then
          call check(nf90_put_var(id_nc, id_bx, bx))
          call check(nf90_put_var(id_nc, id_by, by))
          call check(nf90_put_var(id_nc, id_bz, bz))
        end if
        if (i == 2) call check(nf90_put_var(id_nc, id_b_mag, sqrt(bx**2+by**2+bz**2)))
        if (i == 3) then
          call check(nf90_put_var(id_nc, id_t1_x, reconMetrics(:,:,:)))
          call check(nf90_put_var(id_nc, id_t1_y, by))
          call check(nf90_put_var(id_nc, id_t1_z, bz))
        end if
        if (i == 4) call check(nf90_def_var(id_nc, "c2_t1_mag" , NF90_DOUBLE, dim_ids, id_t1_mag , deflate_level=7))
        if (i == 5) then
          call check(nf90_def_var(id_nc, "c2_t2_x", NF90_DOUBLE, dim_ids, id_t2_x , deflate_level=7))
          call check(nf90_def_var(id_nc, "c2_t2_y", NF90_DOUBLE, dim_ids, id_t2_y , deflate_level=7))
          call check(nf90_def_var(id_nc, "c2_t2_z", NF90_DOUBLE, dim_ids, id_t2_z , deflate_level=7))
        end if
        if (i == 6) call check(nf90_def_var(id_nc, "c2_t2_mag" , NF90_DOUBLE, dim_ids, id_t2_mag , deflate_level=7))
        if (i == 7) then
          call check(nf90_def_var(id_nc, "c2_t3_x", NF90_DOUBLE, dim_ids, id_t3_x , deflate_level=7))
          call check(nf90_def_var(id_nc, "c2_t3_y", NF90_DOUBLE, dim_ids, id_t3_y , deflate_level=7))
          call check(nf90_def_var(id_nc, "c2_t3_z", NF90_DOUBLE, dim_ids, id_t3_z , deflate_level=7))
        end if
        if (i == 8) call check(nf90_def_var(id_nc, "c2_t3_mag" , NF90_DOUBLE, dim_ids, id_t3_mag , deflate_level=7))
      end if 


    call 
    call check(nf90_put_var(id_nc, dy_id, daty))
    call check(nf90_put_var(id_nc, dz_id, datz))

    call check(nf90_close(id_nc))

  end subroutine save_to_netcdf
  
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

    call check(nf90_close(nc_id))

  end subroutine read_netcdf_field_file

  subroutine write_netcdf_vector(filename,dat_name,datx,daty,datz)
    ! Assumes you wish to write out some calculated vector quantity and it's magnitude
    ! on the same mesh as the original magnetic field file.  

    use netcdf

    real(8), dimension(nz,ny,nx), intent(in) :: datx,daty,datz
    character(*), intent(in) :: filename
    character(*), intent(in) :: dat_name
    integer :: nc_id

    integer :: id_nc,xdim_id,ydim_id,zdim_id
    integer :: x_id,y_id,z_id,dx_id,dy_id,dz_id
    integer, dimension(3) :: dim_ids

    call check(nf90_create(filename, NF90_NETCDF4, id_nc))

    ! setup dimensions
    call check(nf90_def_dim(id_nc, "x", nx, xdim_id))
    call check(nf90_def_dim(id_nc, "y", ny, ydim_id))
    call check(nf90_def_dim(id_nc, "z", nz, zdim_id))

    ! create variables for data
    call check(nf90_def_var(id_nc, "x", NF90_DOUBLE, xdim_id, x_id))
    call check(nf90_def_var(id_nc, "y", NF90_DOUBLE, ydim_id, y_id))
    call check(nf90_def_var(id_nc, "z", NF90_DOUBLE, zdim_id, z_id))

    dim_ids = (/ zdim_id,ydim_id,xdim_id /)

    call check(nf90_def_var(id_nc, dat_name//"x", NF90_DOUBLE, dim_ids, dx_id, deflate_level=7))
    call check(nf90_def_var(id_nc, dat_name//"y", NF90_DOUBLE, dim_ids, dy_id, deflate_level=7))
    call check(nf90_def_var(id_nc, dat_name//"z", NF90_DOUBLE, dim_ids, dz_id, deflate_level=7))

    call check(nf90_enddef(id_nc))

    ! Write data

    call check(nf90_put_var(id_nc, x_id, x))
    call check(nf90_put_var(id_nc, y_id, y))
    call check(nf90_put_var(id_nc, z_id, z))
    call check(nf90_put_var(id_nc, dx_id, datx))
    call check(nf90_put_var(id_nc, dy_id, daty))
    call check(nf90_put_var(id_nc, dz_id, datz))

    call check(nf90_close(id_nc))

  end subroutine write_netcdf_vector

  subroutine write_netcdf_scalar(filename,data_name,data)
    ! Writes a scalar quantity on the same mesh as the 
    ! original magnetic field file.  

    use netcdf

    real(8), dimension(nz,ny,nx), intent(in) :: data
    character(*), intent(in) :: filename, data_name
    integer :: nc_id

    integer :: id_nc,xdim_id,ydim_id,zdim_id
    integer :: x_id,y_id,z_id,data_id
    integer, dimension(3) :: dim_ids

    call check(nf90_create(filename, NF90_NETCDF4, id_nc))

    ! setup dimensions
    call check(nf90_def_dim(id_nc, "x", nx, xdim_id))
    call check(nf90_def_dim(id_nc, "y", ny, ydim_id))
    call check(nf90_def_dim(id_nc, "z", nz, zdim_id))

    ! create variables for data
    call check(nf90_def_var(id_nc, "x", NF90_DOUBLE, xdim_id, x_id))
    call check(nf90_def_var(id_nc, "y", NF90_DOUBLE, ydim_id, y_id))
    call check(nf90_def_var(id_nc, "z", NF90_DOUBLE, zdim_id, z_id))

    dim_ids = (/ zdim_id,ydim_id,xdim_id /)

    call check(nf90_def_var(id_nc, data_name, NF90_DOUBLE, dim_ids, data_id,deflate_level=7))

    call check(nf90_enddef(id_nc))

    ! Write data

    call check(nf90_put_var(id_nc, x_id, x))
    call check(nf90_put_var(id_nc, y_id, y))
    call check(nf90_put_var(id_nc, z_id, z))
    call check(nf90_put_var(id_nc, data_id, data))

    call check(nf90_close(id_nc))

  end subroutine write_netcdf_scalar

  subroutine check(istatus)
    use netcdf
    implicit none
    integer, intent (in) :: istatus
    if (istatus /= nf90_noerr) then
    write(*,*) trim(adjustl(nf90_strerror(istatus)))
    end if
  end subroutine check

end module inputOutput
