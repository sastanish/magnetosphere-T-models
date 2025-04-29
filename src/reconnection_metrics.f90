module reconnection_metrics
  implicit none

  integer :: nx,ny,nz
  real(8), dimension(:,:,:), allocatable :: bx,by,bz
  real(8), dimension(:), allocatable :: x,y,z

contains

  subroutine init_no_file(xin,yin,zin,bxin,byin,bzin)
    real(8), dimension(:,:,:), intent(in) :: bxin,byin,bzin
    real(8), dimension(:), intent(in) :: xin,yin,zin

    nx = size(xin)
    ny = size(yin)
    nz = size(zin)

    ! Allocate data in ram
    allocate(x(nx))
    allocate(y(ny))
    allocate(z(nz))
    allocate(bx(nz,ny,nx))
    allocate(by(nz,ny,nx))
    allocate(bz(nz,ny,nx))

    x = xin
    y = yin
    z = zin
    bx = bxin
    by = byin
    bz = bzin

  end subroutine init_no_file

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

  subroutine compute_current(filename,jx,jy,jz,write_data)

    character(*), intent(in) :: filename
    logical, optional, intent(inout) :: write_data

    real(8), dimension(:,:,:), intent(out) :: jx, jy, jz

    integer :: ix, iy, iz
    real(8) :: hx, hy, hz
    integer :: r(3,3) !9-point stencil

    ! Check for writing data flag
    if (.not.present(write_data)) then
      write_data = .false.
    end if

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

    ! compute internal current via second order central diff
    do ix = 1,nx
      if (ix==1) then; r(:,3) = [1, 2, 3];
      else if (ix==nx) then; r(:,3) = [nx-2, nx-1, nx];
      else; r(:,3) = [ix-1, ix, ix+1];
      end if

      do iy = 1,ny
        if (iy==1) then; r(:,2) = [1, 2, 3];
        else if (iy==ny) then; r(:,2) = [ny-2, ny-1, ny];
        else; r(:,2) = [iy-1, iy, iy+1];
        end if

        do iz = 1,nz
          if (iz==1) then; r(:,1) = [1, 2, 3];
          else if (iz==nz) then; r(:,1) = [nz-2, nz-1, nz];
          else; r(:,1) = [iz-1, iz, iz+1];
          end if

          ! compute J=curl(B)
          jx(iz,iy,ix) = (bz(r(2,1),r(3,2),r(2,3))     &
                         -bz(r(2,1),r(1,2),r(2,3)))/(2*hy) &
                       - (by(r(3,1),r(2,2),r(2,3))     &
                         -by(r(1,1),r(2,2),r(2,3)))/(2*hz)

          jy(iz,iy,ix) = (bx(r(3,1),r(2,2),r(2,3))     &
                         -bx(r(1,1),r(2,2),r(2,3)))/(2*hz) &
                       - (bz(r(2,1),r(2,2),r(3,3))     &
                         -bz(r(2,1),r(2,2),r(1,3)))/(2*hx) 

          jz(iz,iy,ix) = (by(r(2,1),r(2,2),r(3,3))     &
                         -by(r(2,1),r(2,2),r(1,3)))/(2*hx) &
                       - (bx(r(2,1),r(3,2),r(2,3))     &
                         -bx(r(2,1),r(1,2),r(2,3)))/(2*hy)
        end do
      end do
    end do

    if (write_data) then
      ix = len(filename) !reusing variable here
      call write_netcdf_vector(trim(filename(1:ix-3))//"_j.nc","j",jx,jy,jz)
    end if

  end subroutine compute_current

  subroutine compute_lorrenz(filename, jx,jy,jz, write_data)

    ! Takes the field current as input. Note, this subroutine modifies
    ! the current arrays as output

    character(*), intent(in) :: filename
    real(8), intent(inout) :: jx(:,:,:), jy(:,:,:), jz(:,:,:)
    logical, optional, intent(inout) :: write_data

    integer :: ix, iy, iz
    real(8) :: dx, dy, dz !Dummy variables for writing jx,y,z

    ! Check for writing data flag
    if (.not.present(write_data)) then
      write_data = .false.
    end if

    ! F = J cross B (writes back out to jx,y,z)
    do ix = 1,nx
      do iy = 1,ny
        do iz = 1,nz

          dx = jx(iz,iy,ix)
          dy = jy(iz,iy,ix)
          dz = jz(iz,iy,ix)

          jx(iz,iy,ix) = dy * bz(iz,iy,ix) - dz * by(iz,iy,ix)
          jy(iz,iy,ix) = dz * bx(iz,iy,ix) - dx * bz(iz,iy,ix)
          jz(iz,iy,ix) = dx * by(iz,iy,ix) - dy * bx(iz,iy,ix)

        end do
      end do
    end do

    if (write_data) then
      ix = len(filename) !reusing variable here
      call write_netcdf_vector(trim(filename(1:ix-3))//"_F.nc","F",jx,jy,jz)
    end if

  end subroutine compute_lorrenz

  subroutine compute_B_fp(filename, Fx, Fy, Fz, write_data)

    ! Takes the lorrenz force as input. Note, this subroutine modifies
    ! the Lorrenz arrays as output

    character(*), intent(in) :: filename
    real(8), intent(inout) :: Fx(:,:,:), Fy(:,:,:), Fz(:,:,:)
    logical, optional, intent(inout) :: write_data

    integer :: ix, iy, iz
    real(8) :: hx, hy, hz, mag_F
    real(8) :: dx, dy, dz !Dummy variables for writing jx,y,z

    ! Check for writing data flag
    if (.not.present(write_data)) then
      write_data = .false.
    end if

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

    ! B_fperp = B cross F/|F| (writes back out to jx,y,z)
    do ix = 1,nx
      do iy = 1,ny
        do iz = 1,nz

          dx = Fx(iz,iy,ix)
          dy = Fy(iz,iy,ix)
          dz = Fz(iz,iy,ix)

          mag_F = sqrt(dx**2 + dy**2 + dz**2)

          Fx(iz,iy,ix) = (dz * by(iz,iy,ix) - dy * bz(iz,iy,ix))/mag_F
          Fy(iz,iy,ix) = (dx * bz(iz,iy,ix) - dz * bx(iz,iy,ix))/mag_F
          Fz(iz,iy,ix) = (dy * bx(iz,iy,ix) - dx * by(iz,iy,ix))/mag_F

        end do
      end do
    end do

    if (write_data) then
      ix = len(filename) !reusing variable here
      call write_netcdf_vector(trim(filename(1:ix-3))//"_B_fp.nc","B_fp",Fx,Fy,Fz)
    end if

  end subroutine compute_B_fp

  subroutine compute_c2_t1_coeff(filename,c2_t1,write_data)

    character(*), intent(in) :: filename
    real(8), intent(out) :: c2_t1(nz,ny,nx)
    logical, optional, intent(inout) :: write_data

    real(8), dimension(nz,ny,nx) :: jx, jy, jz
    real(8), dimension(nz,ny,nx) :: Fx, Fy, Fz
    real(8), dimension(nz,ny,nx) :: B_fpx, B_fpy, B_fpz

    integer :: ix, iy, iz, r(3,3)
    real(8) :: hx, hy, hz
    real(8) :: mag_b(nz,ny,nx), mag_F(nz,ny,nx), lambda(nz,ny,nx)
    real(8) :: curl_B_fpx, curl_B_fpy, curl_B_fpz
    logical :: dummy

    ! Check for writing data flag
    if (.not.present(write_data)) then
      write_data = .false.
    end if

    dummy = .false.

    call compute_current(filename,jx,jy,jz,dummy)
    Fx = jx
    Fy = jy
    Fz = jz
    call compute_lorrenz(filename,Fx,Fy,Fz,dummy)
    B_fpx = Fx
    B_fpy = Fy
    B_fpz = Fz
    call compute_B_fp(filename,B_fpx,B_fpy,B_fpz,dummy)

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

    mag_b = sqrt(bx**2 + by**2 + bz**2)
    mag_F = sqrt(Fx**2 + Fy**2 + Fz**2)
    lambda = (jx*B_fpx + jy*B_fpy + jz*B_fpz)/mag_b**2

    do ix = 1,nx
      if (ix==1) then; r(:,3) = [1, 2, 3];
      else if (ix==nx) then; r(:,3) = [nx-2, nx-1, nx];
      else; r(:,3) = [ix-1, ix, ix+1];
      end if

      do iy = 1,ny
        if (iy==1) then; r(:,2) = [1, 2, 3];
        else if (iy==ny) then; r(:,2) = [ny-2, ny-1, ny];
        else; r(:,2) = [iy-1, iy, iy+1];
        end if

        do iz = 1,nz
          if (iz==1) then; r(:,1) = [1, 2, 3];
          else if (iz==nz) then; r(:,1) = [nz-2, nz-1, nz];
          else; r(:,1) = [iz-1, iz, iz+1];
          end if

          !curl(B_fperp)
          curl_B_fpx = (B_fpz(r(2,1),r(3,2),r(2,3))     &
                       -B_fpz(r(2,1),r(1,2),r(2,3)))/(2*hy) &
                     - (B_fpy(r(3,1),r(2,2),r(2,3))     &
                       -B_fpy(r(1,1),r(2,2),r(2,3)))/(2*hz)

          curl_B_fpy = (B_fpx(r(3,1),r(2,2),r(2,3))     &
                       -B_fpx(r(1,1),r(2,2),r(2,3)))/(2*hz) &
                     - (B_fpz(r(2,1),r(2,2),r(3,3))     &
                       -B_fpz(r(2,1),r(2,2),r(1,3)))/(2*hx) 

          curl_B_fpz = (B_fpy(r(2,1),r(2,2),r(3,3))     &
                       -B_fpy(r(2,1),r(2,2),r(1,3)))/(2*hx) &
                     - (B_fpx(r(2,1),r(3,2),r(2,3))     &
                       -B_fpx(r(2,1),r(1,2),r(2,3)))/(2*hy)

          ! (lambda * curl(B_fperp) dot F/|F| - Grad(lambda) dot B)/|B|
          c2_t1(iz,iy,ix) = (&
                            lambda(iz,iy,ix)*(curl_B_fpx*Fx(iz,iy,ix) + curl_B_fpy*Fy(iz,iy,ix) + curl_B_fpz*Fz(iz,iy,ix))/mag_F(iz,iy,ix)&
                          - ((lambda(r(2,1),r(2,2),r(3,3)) - lambda(r(2,1),r(2,2),r(1,3)))/(2*hx)*bx(iz,iy,ix) &
                            +(lambda(r(2,1),r(3,2),r(2,3)) - lambda(r(2,1),r(1,2),r(2,3)))/(2*hy)*by(iz,iy,ix) &
                            +(lambda(r(3,1),r(2,2),r(2,3)) - lambda(r(1,1),r(2,2),r(2,3)))/(2*hz)*bz(iz,iy,ix) &
                            )&
                            )/mag_b(iz,iy,ix)
        end do
      end do
    end do

    if (write_data) then
      ix = len(filename) !reusing variable here
      call write_netcdf_scalar(trim(filename(1:ix-3))//"_c2_t1.nc","c2_t1",c2_t1)
    end if

  end subroutine compute_c2_t1_coeff

  subroutine compute_c2_t2_coeff(filename,c2_t2,write_data)

    character(*), intent(in) :: filename
    real(8), intent(out) :: c2_t2(nz,ny,nx)
    logical, optional, intent(inout) :: write_data

    real(8), dimension(nz,ny,nx) :: jx, jy, jz
    real(8), dimension(nz,ny,nx) :: Fx, Fy, Fz
    real(8), dimension(nz,ny,nx) :: B_fpx, B_fpy, B_fpz

    integer :: ix, iy, iz, r(3,3)
    real(8) :: hx, hy, hz
    real(8) :: mag_b(nz,ny,nx), mag_F(nz,ny,nx), lambda(nz,ny,nx), alpha(nz,ny,nx)
    real(8) :: curl_B_fpx, curl_B_fpy, curl_B_fpz
    logical :: dummy

    ! Check for writing data flag
    if (.not.present(write_data)) then
      write_data = .false.
    end if

    dummy = .false.

    call compute_current(filename,jx,jy,jz,dummy)
    Fx = jx
    Fy = jy
    Fz = jz
    call compute_lorrenz(filename,Fx,Fy,Fz,dummy)
    B_fpx = Fx
    B_fpy = Fy
    B_fpz = Fz
    call compute_B_fp(filename,B_fpx,B_fpy,B_fpz,dummy)

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

    mag_b = sqrt(bx**2 + by**2 + bz**2)
    mag_F = sqrt(Fx**2 + Fy**2 + Fz**2)
    lambda = (jx*B_fpx + jy*B_fpy + jz*B_fpz)/mag_b**2
    alpha = (jx*bx + jy*by + jz*bz)/mag_b**2

    do ix = 1,nx
      if (ix==1) then; r(:,3) = [1, 2, 3];
      else if (ix==nx) then; r(:,3) = [nx-2, nx-1, nx];
      else; r(:,3) = [ix-1, ix, ix+1];
      end if

      do iy = 1,ny
        if (iy==1) then; r(:,2) = [1, 2, 3];
        else if (iy==ny) then; r(:,2) = [ny-2, ny-1, ny];
        else; r(:,2) = [iy-1, iy, iy+1];
        end if

        do iz = 1,nz
          if (iz==1) then; r(:,1) = [1, 2, 3];
          else if (iz==nz) then; r(:,1) = [nz-2, nz-1, nz];
          else; r(:,1) = [iz-1, iz, iz+1];
          end if

          !curl(B_fperp)
          curl_B_fpx = (B_fpz(r(2,1),r(3,2),r(2,3))     &
                       -B_fpz(r(2,1),r(1,2),r(2,3)))/(2*hy) &
                     - (B_fpy(r(3,1),r(2,2),r(2,3))     &
                       -B_fpy(r(1,1),r(2,2),r(2,3)))/(2*hz)

          curl_B_fpy = (B_fpx(r(3,1),r(2,2),r(2,3))     &
                       -B_fpx(r(1,1),r(2,2),r(2,3)))/(2*hz) &
                     - (B_fpz(r(2,1),r(2,2),r(3,3))     &
                       -B_fpz(r(2,1),r(2,2),r(1,3)))/(2*hx) 

          curl_B_fpz = (B_fpy(r(2,1),r(2,2),r(3,3))     &
                       -B_fpy(r(2,1),r(2,2),r(1,3)))/(2*hx) &
                     - (B_fpx(r(2,1),r(3,2),r(2,3))     &
                       -B_fpx(r(2,1),r(1,2),r(2,3)))/(2*hy)

          ! lambda (curl(B_fperp) dot B_fperp/|B|^2 + alpha)/|B|
          c2_t2(iz,iy,ix) = lambda(iz,iy,ix)*( &
                  ( curl_B_fpx*B_fpx(iz,iy,ix) &
                  + curl_B_fpy*B_fpy(iz,iy,ix) &
                  + curl_B_fpz*B_fpz(iz,iy,ix) )/mag_b(iz,iy,ix)**2 &
              + alpha(iz,iy,ix))/mag_b(iz,iy,ix)
        end do
      end do
    end do

    if (write_data) then
      ix = len(filename) !reusing variable here
      call write_netcdf_scalar(trim(filename(1:ix-3))//"_c2_t2.nc","c2_t2",c2_t2)
    end if

  end subroutine compute_c2_t2_coeff

  subroutine compute_c2_t3_coeff(filename,c2_t3,write_data)

    character(*), intent(in) :: filename
    real(8), intent(out) :: c2_t3(nz,ny,nx)
    logical, optional, intent(inout) :: write_data

    real(8), dimension(nz,ny,nx) :: jx, jy, jz
    real(8), dimension(nz,ny,nx) :: Fx, Fy, Fz
    real(8), dimension(nz,ny,nx) :: B_fpx, B_fpy, B_fpz

    integer :: ix, iy, iz, r(3,3)
    real(8) :: hx, hy, hz
    real(8) :: mag_b(nz,ny,nx), mag_F(nz,ny,nx), lambda(nz,ny,nx), alpha(nz,ny,nx)
    real(8) :: alpha_dx(nz,ny,nx), alpha_dy(nz,ny,nx), alpha_dz(nz,ny,nx)
    logical :: dummy

    ! Check for writing data flag
    if (.not.present(write_data)) then
      write_data = .false.
    end if

    dummy = .false.

    call compute_current(filename,jx,jy,jz,dummy)
    Fx = jx
    Fy = jy
    Fz = jz
    call compute_lorrenz(filename,Fx,Fy,Fz,dummy)
    B_fpx = Fx
    B_fpy = Fy
    B_fpz = Fz
    call compute_B_fp(filename,B_fpx,B_fpy,B_fpz,dummy)

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

    mag_b = sqrt(bx**2 + by**2 + bz**2)
    alpha = (jx*bx + jy*by + jz*bz)/mag_b**2

    do ix = 1,nx
      if (ix==1) then; r(:,3) = [1, 2, 3];
      else if (ix==nx) then; r(:,3) = [nx-2, nx-1, nx];
      else; r(:,3) = [ix-1, ix, ix+1];
      end if

      do iy = 1,ny
        if (iy==1) then; r(:,2) = [1, 2, 3];
        else if (iy==ny) then; r(:,2) = [ny-2, ny-1, ny];
        else; r(:,2) = [iy-1, iy, iy+1];
        end if

        do iz = 1,nz
          if (iz==1) then; r(:,1) = [1, 2, 3];
          else if (iz==nz) then; r(:,1) = [nz-2, nz-1, nz];
          else; r(:,1) = [iz-1, iz, iz+1];
          end if

          alpha_dx(iz,iy,ix) = (alpha(r(2,1),r(2,2),r(3,3))-alpha(r(2,1),r(2,2),r(1,3)))/(2*hx)
          alpha_dy(iz,iy,ix) = (alpha(r(2,1),r(3,2),r(2,3))-alpha(r(2,1),r(1,2),r(2,3)))/(2*hy)
          alpha_dz(iz,iy,ix) = (alpha(r(3,1),r(2,2),r(2,3))-alpha(r(1,1),r(2,2),r(2,3)))/(2*hz)

        end do
      end do
    end do

    c2_t3 = sqrt( (alpha_dy*bz-alpha_dz*by)**2 + (alpha_dz*bx-alpha_dx*bz)**2 + (alpha_dx*by-alpha_dy*bx)**2 )/mag_b

    if (write_data) then
      ix = len(filename) !reusing variable here
      call write_netcdf_scalar(trim(filename(1:ix-3))//"_c2_t3.nc","c2_t3",c2_t3)
    end if

  end subroutine compute_c2_t3_coeff

  subroutine cleanup()
    deallocate(bx)
    deallocate(by)
    deallocate(bz)
    deallocate(x)
    deallocate(y)
    deallocate(z)
  end subroutine cleanup

  subroutine check(istatus)
    use netcdf
    implicit none
    integer, intent (in) :: istatus
    if (istatus /= nf90_noerr) then
    write(*,*) trim(adjustl(nf90_strerror(istatus)))
    end if
  end subroutine check

end module reconnection_metrics
