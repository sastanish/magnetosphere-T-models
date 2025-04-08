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

    call check(nf90_put_var(id_nc, x_id, x))
    call check(nf90_put_var(id_nc, y_id, y))
    call check(nf90_put_var(id_nc, z_id, z))
    call check(nf90_put_var(id_nc, dx_id, datx))
    call check(nf90_put_var(id_nc, dy_id, daty))
    call check(nf90_put_var(id_nc, dz_id, datz))
    call check(nf90_put_var(id_nc, dm_id, magdat))

    call check(nf90_close(id_nc))

  end subroutine write_netcdf_metric_file

  subroutine compute_current(filename, jx,jy,jz,mag_j)

    real(8), intent(inout), dimension(nz,ny,nx) :: jx, jy, jz, mag_j

    integer :: ix, iy, iz
    integer :: ixp, ixm, iyp, iym, izp, izm
    real(8) :: hx, hy, hz

    character(*), intent(in) :: filename

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

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

          ! compute
          jx(iz,iy,ix) = o2_central_diff(hy,bz(iz ,iyp,ix ),bz(iz,iy,ix),bz(iz ,iym,ix ))&
                       - o2_central_diff(hz,by(izp,iy ,ix ),by(iz,iy,ix),by(izm,iy ,ix ))
          jy(iz,iy,ix) = o2_central_diff(hz,bx(izp,iy ,ix ),bx(iz,iy,ix),bx(izm,iy ,ix ))&
                       - o2_central_diff(hx,bz(iz ,iy ,ixp),bz(iz,iy,ix),bz(iz ,iy ,ixm))
          jz(iz,iy,ix) = o2_central_diff(hx,by(iz ,iy ,ixp),by(iz,iy,ix),by(iz ,iy ,ixm))&
                       - o2_central_diff(hy,bx(iz ,iyp,ix ),bx(iz,iy,ix),bx(iz ,iym,ix ))
          mag_j(iz,iy,ix) = sqrt(jx(iz ,iy ,ix)**2 + jy(iz ,iy ,ix)**2 + jz(iz ,iy ,ix)**2 )

        end do
      end do
    end do

    ix = len(filename) !reusing variable here
    call write_netcdf_metric_file(trim(filename(1:ix-3))//"_j.nc","j",jx,jy,jz,mag_j)

  end subroutine compute_current

  subroutine compute_first_term(filename, jx,jy,jz, mag_j)

    real(8), intent(in), dimension(nz,ny,nx) :: jx, jy, jz, mag_j
    character(*), intent(in) :: filename

    real(8), dimension(nz,ny,nx) :: cx, cy, cz, mag_c

    integer :: ix, iy, iz, ixp, ixm, iyp, iym, izp, izm
    real(8) :: hx, hy, hz
    real(8) :: mag_b
    real(8) :: Fx, Fy, Fz, mag_F
    real(8) :: B_fpx, B_fpy, B_fpz, B_fpx_xm, B_fpy_xm, B_fpz_xm, B_fpx_xp, B_fpy_xp, B_fpz_xp, B_fpx_ym, B_fpy_ym, B_fpz_ym, B_fpx_yp, B_fpy_yp, B_fpz_yp, B_fpx_zm, B_fpy_zm, B_fpz_zm, B_fpx_zp, B_fpy_zp, B_fpz_zp
    real(8) :: curl_B_fpx, curl_B_fpy, curl_B_fpz
    real(8) :: lambda, lambda_zp, lambda_zm, lambda_yp, lambda_ym, lambda_xp, lambda_xm
    real(8) :: omega_1, coeff

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

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

          ! compute
          Fx = jy(iz,iy,ix) * bz(iz,iy,ix) - jz(iz,iy,ix) * by(iz,iy,ix)
          Fy = jz(iz,iy,ix) * bx(iz,iy,ix) - jx(iz,iy,ix) * bz(iz,iy,ix)
          Fz = jx(iz,iy,ix) * by(iz,iy,ix) - jy(iz,iy,ix) * bx(iz,iy,ix)
          mag_F = sqrt(Fx**2 + Fy**2 + Fz**2)
          mag_b = sqrt(bx(iz ,iy ,ix)**2 + by(iz ,iy ,ix)**2 + bz(iz ,iy ,ix)**2 )

          B_fpx = (by(iz,iy,ix)*Fz - bz(iz,iy,ix)*Fy)/mag_F
          B_fpy = (bz(iz,iy,ix)*Fx - bx(iz,iy,ix)*Fz)/mag_F
          B_fpz = (bx(iz,iy,ix)*Fy - by(iz,iy,ix)*Fx)/mag_F

          B_fpx_zp = (by(izp,iy,ix)*Fz - bz(izp,iy,ix)*Fy)/mag_F
          B_fpx_zm = (by(izm,iy,ix)*Fz - bz(izm,iy,ix)*Fy)/mag_F
          B_fpx_yp = (by(iz,iyp,ix)*Fz - bz(iz,iyp,ix)*Fy)/mag_F
          B_fpx_ym = (by(iz,iym,ix)*Fz - bz(iz,iym,ix)*Fy)/mag_F
          B_fpy_zp = (bz(izp,iy,ix)*Fx - bx(izp,iy,ix)*Fz)/mag_F
          B_fpy_zm = (bz(izm,iy,ix)*Fx - bx(izm,iy,ix)*Fz)/mag_F
          B_fpy_xp = (bz(iz,iy,ixp)*Fx - bx(iz,iy,ixp)*Fz)/mag_F
          B_fpy_xm = (bz(iz,iy,ixm)*Fx - bx(iz,iy,ixm)*Fz)/mag_F
          B_fpz_yp = (bx(iz,iyp,ix)*Fy - by(iz,iyp,ix)*Fx)/mag_F
          B_fpz_ym = (bx(iz,iym,ix)*Fy - by(iz,iym,ix)*Fx)/mag_F
          B_fpz_xp = (bx(iz,iy,ixp)*Fy - by(iz,iy,ixp)*Fx)/mag_F
          B_fpz_xm = (bx(iz,iy,ixm)*Fy - by(iz,iy,ixm)*Fx)/mag_F

          curl_B_fpx = o2_central_diff(hy,B_fpz_yp,B_fpz,B_fpz_ym)&
                     - o2_central_diff(hz,B_fpy_zp,B_fpy,B_fpy_zm)
          curl_B_fpy = o2_central_diff(hz,B_fpx_zp,B_fpx,B_fpx_zm)&
                     - o2_central_diff(hx,B_fpz_xp,B_fpz,B_fpz_xm)
          curl_B_fpz = o2_central_diff(hx,B_fpy_xp,B_fpy,B_fpy_xm)&
                     - o2_central_diff(hy,B_fpx_yp,B_fpx,B_fpx_ym)

          lambda = (jx(iz,iy,ix)*B_fpx + jy(iz,iy,ix)*B_fpy + jz(iz,iy,ix)*B_fpz )/mag_b**2
          lambda_zp = (jx(izp,iy,ix)*B_fpx + jy(izp,iy,ix)*B_fpy + jz(izp,iy,ix)*B_fpz )/mag_b**2
          lambda_zm = (jx(izm,iy,ix)*B_fpx + jy(izm,iy,ix)*B_fpy + jz(izm,iy,ix)*B_fpz )/mag_b**2
          lambda_yp = (jx(iz,iyp,ix)*B_fpx + jy(iz,iyp,ix)*B_fpy + jz(iz,iyp,ix)*B_fpz )/mag_b**2
          lambda_ym = (jx(iz,iym,ix)*B_fpx + jy(iz,iym,ix)*B_fpy + jz(iz,iym,ix)*B_fpz )/mag_b**2
          lambda_xp = (jx(iz,iy,ixp)*B_fpx + jy(iz,iy,ixp)*B_fpy + jz(iz,iy,ixp)*B_fpz )/mag_b**2
          lambda_xm = (jx(iz,iy,ixm)*B_fpx + jy(iz,iy,ixm)*B_fpy + jz(iz,iy,ixm)*B_fpz )/mag_b**2

          omega_1 = (curl_B_fpx*Fx + curl_B_fpy*Fy + curl_B_fpz*Fz)/mag_F

          coeff = lambda*omega_1 - (&
                    o2_central_diff(hx,lambda_xp,lambda,lambda_xm)*bx(ix,iy,iz)&
                  + o2_central_diff(hy,lambda_yp,lambda,lambda_ym)*by(ix,iy,iz)&
                  + o2_central_diff(hz,lambda_zp,lambda,lambda_zm)*bz(ix,iy,iz))

          cx(iz,iy,ix) = coeff* Fx/(mag_F*mag_b)
          cy(iz,iy,ix) = coeff* Fy/(mag_F*mag_b)
          cz(iz,iy,ix) = coeff* Fz/(mag_F*mag_b)
          mag_c(iz,iy,ix) = coeff

        end do
      end do
    end do

    ix = len(filename)! reusing variable here
    call write_netcdf_metric_file(trim(filename(1:ix-3))//"_C2_t1.nc","t1",cx,cy,cz,mag_c)

  end subroutine compute_first_term

  subroutine compute_second_term(filename, jx,jy,jz, mag_j)

    real(8), intent(in), dimension(nz,ny,nx) :: jx, jy, jz, mag_j
    character(*), intent(in) :: filename

    real(8), dimension(nz,ny,nx) :: cx, cy, cz, mag_c

    integer :: ix, iy, iz, ixp, ixm, iyp, iym, izp, izm
    real(8) :: hx, hy, hz
    real(8) :: mag_b
    real(8) :: Fx, Fy, Fz, mag_F
    real(8) :: B_fpx, B_fpy, B_fpz, B_fpx_xm, B_fpy_xm, B_fpz_xm, B_fpx_xp, B_fpy_xp, B_fpz_xp, B_fpx_ym, B_fpy_ym, B_fpz_ym, B_fpx_yp, B_fpy_yp, B_fpz_yp, B_fpx_zm, B_fpy_zm, B_fpz_zm, B_fpx_zp, B_fpy_zp, B_fpz_zp
    real(8) :: curl_B_fpx, curl_B_fpy, curl_B_fpz
    real(8) :: omega_2, coeff, alpha, lambda

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

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

          ! compute
          Fx = jy(iz,iy,ix) * bz(iz,iy,ix) - jz(iz,iy,ix) * by(iz,iy,ix)
          Fy = jz(iz,iy,ix) * bx(iz,iy,ix) - jx(iz,iy,ix) * bz(iz,iy,ix)
          Fz = jx(iz,iy,ix) * by(iz,iy,ix) - jy(iz,iy,ix) * bx(iz,iy,ix)
          mag_F = sqrt(Fx**2 + Fy**2 + Fz**2)
          mag_b = sqrt(bx(iz ,iy ,ix)**2 + by(iz ,iy ,ix)**2 + bz(iz ,iy ,ix)**2 )

          B_fpx = (by(iz,iy,ix)*Fz - bz(iz,iy,ix)*Fy)/mag_F
          B_fpy = (bz(iz,iy,ix)*Fx - bx(iz,iy,ix)*Fz)/mag_F
          B_fpz = (bx(iz,iy,ix)*Fy - by(iz,iy,ix)*Fx)/mag_F

          B_fpx_zp = (by(izp,iy,ix)*Fz - bz(izp,iy,ix)*Fy)/mag_F
          B_fpx_zm = (by(izm,iy,ix)*Fz - bz(izm,iy,ix)*Fy)/mag_F
          B_fpx_yp = (by(iz,iyp,ix)*Fz - bz(iz,iyp,ix)*Fy)/mag_F
          B_fpx_ym = (by(iz,iym,ix)*Fz - bz(iz,iym,ix)*Fy)/mag_F
          B_fpy_zp = (bz(izp,iy,ix)*Fx - bx(izp,iy,ix)*Fz)/mag_F
          B_fpy_zm = (bz(izm,iy,ix)*Fx - bx(izm,iy,ix)*Fz)/mag_F
          B_fpy_xp = (bz(iz,iy,ixp)*Fx - bx(iz,iy,ixp)*Fz)/mag_F
          B_fpy_xm = (bz(iz,iy,ixm)*Fx - bx(iz,iy,ixm)*Fz)/mag_F
          B_fpz_yp = (bx(iz,iyp,ix)*Fy - by(iz,iyp,ix)*Fx)/mag_F
          B_fpz_ym = (bx(iz,iym,ix)*Fy - by(iz,iym,ix)*Fx)/mag_F
          B_fpz_xp = (bx(iz,iy,ixp)*Fy - by(iz,iy,ixp)*Fx)/mag_F
          B_fpz_xm = (bx(iz,iy,ixm)*Fy - by(iz,iy,ixm)*Fx)/mag_F

          curl_B_fpx = o2_central_diff(hy,B_fpz_yp,B_fpz,B_fpz_ym)&
                     - o2_central_diff(hz,B_fpy_zp,B_fpy,B_fpy_zm)
          curl_B_fpy = o2_central_diff(hz,B_fpx_zp,B_fpx,B_fpx_zm)&
                     - o2_central_diff(hx,B_fpz_xp,B_fpz,B_fpz_xm)
          curl_B_fpz = o2_central_diff(hx,B_fpy_xp,B_fpy,B_fpy_xm)&
                     - o2_central_diff(hy,B_fpx_yp,B_fpx,B_fpx_ym)

          lambda = (jx(iz,iy,ix)*B_fpx + jy(iz,iy,ix)*B_fpy + jz(iz,iy,ix)*B_fpz )/mag_b**2
          alpha  = (jx(iz,iy,ix) *bx(iz,iy,ix)  + jy(iz,iy,ix) *by(iz,iy,ix)  + jz(iz,iy,ix) * bz(iz,iy,ix) )/mag_b**2

          omega_2 = (curl_B_fpx*B_fpx + curl_B_fpy*B_fpy + curl_B_fpz*B_fpz)/mag_b**2

          coeff = lambda*(alpha + omega_2)

          cx(iz,iy,ix) = coeff* B_fpx/mag_b
          cy(iz,iy,ix) = coeff* B_fpy/mag_b
          cz(iz,iy,ix) = coeff* B_fpz/mag_b
          mag_c(iz,iy,ix) = coeff

        end do
      end do
    end do

    ix = len(filename)! reusing variable here
    call write_netcdf_metric_file(trim(filename(1:ix-3))//"_C2_t2.nc","t2",cx,cy,cz,mag_c)

  end subroutine compute_second_term

  subroutine compute_third_term(filename, jx,jy,jz, mag_j)

    real(8), intent(in), dimension(nz,ny,nx) :: jx, jy, jz, mag_j
    character(*), intent(in) :: filename

    real(8), dimension(nz,ny,nx) :: cx, cy, cz, mag_c

    integer :: ix, iy, iz, ixp, ixm, iyp, iym, izp, izm
    real(8) :: hx, hy, hz
    real(8) :: mag_b
    real(8) :: alpha, alpha_zp, alpha_zm, alpha_yp, alpha_ym, alpha_xp, alpha_xm

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

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

          mag_b = sqrt(bx(iz ,iy ,ix)**2 + by(iz ,iy ,ix)**2 + bz(iz ,iy ,ix)**2 )

          alpha =    (jx(iz,iy,ix) *bx(iz,iy,ix)  + jy(iz,iy,ix) *by(iz,iy,ix)  + jz(iz,iy,ix) * bz(iz,iy,ix) )/mag_b**2
          alpha_zp = (jx(izp,iy,ix)*bx(izp,iy,ix) + jy(izp,iy,ix)*by(izp,iy,ix) + jz(izp,iy,ix)* bz(izp,iy,ix))/mag_b**2
          alpha_zm = (jx(izm,iy,ix)*bx(izm,iy,ix) + jy(izm,iy,ix)*by(izm,iy,ix) + jz(izm,iy,ix)* bz(izm,iy,ix))/mag_b**2
          alpha_yp = (jx(iz,iyp,ix)*bx(iz,iyp,ix) + jy(iz,iyp,ix)*by(iz,iyp,ix) + jz(iz,iyp,ix)* bz(iz,iyp,ix))/mag_b**2
          alpha_ym = (jx(iz,iym,ix)*bx(iz,iym,ix) + jy(iz,iym,ix)*by(iz,iym,ix) + jz(iz,iym,ix)* bz(iz,iym,ix))/mag_b**2
          alpha_xp = (jx(iz,iy,ixp)*bx(iz,iy,ixp) + jy(iz,iy,ixp)*by(iz,iy,ixp) + jz(iz,iy,ixp)* bz(iz,iy,ixp))/mag_b**2
          alpha_xm = (jx(iz,iy,ixm)*bx(iz,iy,ixm) + jy(iz,iy,ixm)*by(iz,iy,ixm) + jz(iz,iy,ixm)* bz(iz,iy,ixm))/mag_b**2

          cx(iz,iy,ix) =( o2_central_diff(hy,alpha_yp,alpha,alpha_ym)*bz(iz,iy,ix)&
                        - o2_central_diff(hz,alpha_zp,alpha,alpha_zm)*by(iz,iy,ix))/mag_b**2
          cy(iz,iy,ix) =( o2_central_diff(hz,alpha_zp,alpha,alpha_zm)*bx(iz,iy,ix)&
                        - o2_central_diff(hx,alpha_xp,alpha,alpha_xm)*bz(iz,iy,ix))/mag_b**2
          cz(iz,iy,ix) =( o2_central_diff(hx,alpha_xp,alpha,alpha_xm)*by(iz,iy,ix)&
                        - o2_central_diff(hy,alpha_yp,alpha,alpha_ym)*bx(iz,iy,ix))/mag_b**2
          mag_c(iz,iy,iz) = sqrt(cx(iz ,iy ,ix)**2 + cy(iz ,iy ,ix)**2 + cz(iz ,iy ,ix)**2 )

        end do
      end do
    end do

    ix = len(filename)! reusing variable here
    call write_netcdf_metric_file(trim(filename(1:ix-3))//"_C2_t3.nc","t3",cx,cy,cz,mag_c)

  end subroutine compute_third_term

  subroutine calculate_metrics(filename)

    character(*), intent(in) :: filename
    real(8), dimension(:,:,:), allocatable :: jx, jy, jz, mag_j

    call read_netcdf_field_file(filename)

    allocate(jx(nz,ny,nx))
    allocate(jy(nz,ny,nx))
    allocate(jz(nz,ny,nx))
    allocate(mag_j(nz,ny,nx))

    call compute_current(filename,jx,jy,jz,mag_j)

    call compute_first_term(filename, jx, jy, jz, mag_j)
    call compute_second_term(filename, jx, jy, jz, mag_j)
    call compute_third_term(filename, jx, jy, jz, mag_j)

    deallocate(jx)
    deallocate(jy)
    deallocate(jz)
    deallocate(mag_j)

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
