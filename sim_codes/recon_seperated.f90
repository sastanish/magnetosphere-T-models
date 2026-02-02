program main

  !$ use omp_lib

  implicit none

  real(8), dimension(:,:,:), allocatable :: Bx,By,Bz,magB,Fx,Fy,Fz,magF
  real(8), dimension(:,:,:), allocatable :: Jx,Jy,Jz,Bfx,Bfy,Bfz
  real(8), dimension(:,:,:), allocatable :: lambda,alpha,rate,Lorr_rate,Wind_rate,Alig_rate

  real(8), dimension(:), allocatable :: x,y,z

  integer :: i,j,k,nx,ny,nz
  real(8) :: hx,hy,hz
  real(8) :: Gax,Gay,Gaz,omega1,omega2,BLG,cBfx,cBfy,cBfz,rate_x,rate_y,rate_z
  real(8) :: lrate_x,lrate_y,lrate_z,wrate_x,wrate_y,wrate_z,arate_x,arate_y,arate_z

  !File handling
  integer :: fileind, start_ind, end_ind
  character(4) :: str_ind, start_str, end_str

  call GET_COMMAND_ARGUMENT(1,start_str)
  call GET_COMMAND_ARGUMENT(2,end_str)
  read(start_str,'(I5)') start_ind
  read(end_str,'(I5)') end_ind

  do fileind = start_ind,end_ind

    write( str_ind, '(I4)' ) fileind
    print *, 'reading: '//'output_'//trim(adjustl(str_ind))//'.nc'
  
    call load_field_from_netcdf('output_'//trim(adjustl(str_ind))//'.nc',x,y,z,Bx,By,Bz)

    nx = size(x)
    ny = size(y)
    nz = size(z)

    ! ALLOCATIONS
    allocate(magB(nx,ny,nz))
    allocate(Fx(nx,ny,nz))
    allocate(Fy(nx,ny,nz))
    allocate(Fz(nx,ny,nz))
    allocate(magF(nx,ny,nz))
    allocate(Jx(nx,ny,nz))
    allocate(Jy(nx,ny,nz))
    allocate(Jz(nx,ny,nz))
    allocate(Bfx(nx,ny,nz))
    allocate(Bfy(nx,ny,nz))
    allocate(Bfz(nx,ny,nz))
    allocate(alpha(nx,ny,nz))
    allocate(lambda(nx,ny,nz))
    allocate(rate(nx,ny,nz))
    allocate(Lorr_rate(nx,ny,nz))
    allocate(Wind_rate(nx,ny,nz))
    allocate(Alig_rate(nx,ny,nz))

    !!!!!!!!!!!!!!!!!!
    ! Calculate rate !
    !!!!!!!!!!!!!!!!!!

    hx = abs(x(2)-x(1))
    hy = abs(y(2)-y(1))
    hz = abs(z(2)-z(1))
    !$OMP Parallel shared(Bx,By,Bz,magB,Fx,Fy,Fz,Jx,Jy,Jz,Bfx,Bfy,Bfz,lambda,alpha) private(i,j,k)
    !$OMP DO COLLAPSE(3)
    do k=2,nz-1
    do j=2,ny-1
    do i=2,nx-1
      ! |B|
      magB(i,j,k) = sqrt(Bx(i,j,k)**2 + By(i,j,k)**2 + Bz(i,j,k)**2)
      ! J = Curl(B)
      Jx(i,j,k) = 1/(2*hy) * (Bz(i,j+1,k) - Bz(i,j-1,k)) - 1/(2*hz) * (By(i,j,k+1) - By(i,j,k-1))
      Jy(i,j,k) = 1/(2*hz) * (Bx(i,j,k+1) - Bx(i,j,k-1)) - 1/(2*hx) * (Bz(i+1,j,k) - Bz(i-1,j,k))
      Jz(i,j,k) = 1/(2*hx) * (By(i+1,j,k) - By(i-1,j,k)) - 1/(2*hy) * (Bx(i,j+1,k) - Bx(i,j-1,k))
      ! F = J x B
      Fx(i,j,k) = Jy(i,j,k)*Bz(i,j,k) - Jz(i,j,k)*By(i,j,k) 
      Fy(i,j,k) = Jz(i,j,k)*Bx(i,j,k) - Jx(i,j,k)*Bz(i,j,k) 
      Fz(i,j,k) = Jx(i,j,k)*By(i,j,k) - Jy(i,j,k)*Bx(i,j,k) 
      magF(i,j,k) = sqrt(Fx(i,j,k)**2 + Fy(i,j,k)**2 + Fz(i,j,k)**2)
      ! B_f = B x F/|F|
      Bfx(i,j,k) = (By(i,j,k)*Fz(i,j,k) - Bz(i,j,k)*Fy(i,j,k) )/magF(i,j,k)
      Bfy(i,j,k) = (Bz(i,j,k)*Fx(i,j,k) - Bx(i,j,k)*Fz(i,j,k) )/magF(i,j,k)
      Bfz(i,j,k) = (Bx(i,j,k)*Fy(i,j,k) - By(i,j,k)*Fx(i,j,k) )/magF(i,j,k)
      ! alpha = j * b / |B|^2
      alpha(i,j,k) = (Jx(i,j,k)*Bx(i,j,k) + Jy(i,j,k)*By(i,j,k) + Jz(i,j,k)*Bz(i,j,k) )/magB(i,j,k)**2
      !lamb = j * B_f / |B|^2
      lambda(i,j,k) = (Jx(i,j,k)*Bfx(i,j,k) + Jy(i,j,k)*Bfy(i,j,k) + Jz(i,j,k)*Bfz(i,j,k) )/magB(i,j,k)**2
    end do
    end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP Parallel shared(rate,Lorr_rate,Wind_rate,Alig_rate) private(i,j,k,cBfx,cBfy,cBfz,omega1,omega2,lrate_x,lrate_y,lrate_z,wrate_x,wrate_y,wrate_z,arate_x,arate_y,arate_z,Gax,Gay,Gaz)
    !$OMP DO COLLAPSE(3)
    do k=3,nz-2
    do j=3,ny-2
    do i=3,nx-2
      ! Curl(B_f)
      cBfx = 1/(2*hy) * (Bfz(i,j+1,k) - Bfz(i,j-1,k)) - 1/(2*hz) * (Bfy(i,j,k+1) - Bfy(i,j,k-1))
      cBfy = 1/(2*hz) * (Bfx(i,j,k+1) - Bfx(i,j,k-1)) - 1/(2*hx) * (Bfz(i+1,j,k) - Bfz(i-1,j,k))
      cBfz = 1/(2*hx) * (Bfy(i+1,j,k) - Bfy(i-1,j,k)) - 1/(2*hy) * (Bfx(i,j+1,k) - Bfx(i,j-1,k))
      ! w1 = Curl(B_f) * F/|F|
       omega1 = ( cBfx * Fx(i,j,k) + cBfy * Fy(i,j,k) + cBfz * Fz(i,j,k) ) / magF(i,j,k)
      ! w2 = curl(B_f) * B_f / |B|^2
       omega2 = ( cBfx * Bfx(i,j,k) + cBfy * Bfy(i,j,k) + cBfz * Bfz(i,j,k) ) / magB(i,j,k)**2
      ! GLB = Grad(lambda) * B
      BLG = ( 1/(2*hx) * (lambda(i+1,j,k) - lambda(i-1,j,k)) * Bx(i,j,k) &
            + 1/(2*hy) * (lambda(i,j+1,k) - lambda(i,j-1,k)) * By(i,j,k) &
            + 1/(2*hz) * (lambda(i,j,k+1) - lambda(i,j,k-1)) * Bz(i,j,k) )
      ! Ga = Grad(alpha)
      Gax = 1/(2*hx) * (alpha(i+1,j,k) - alpha(i-1,j,k))
      Gay = 1/(2*hy) * (alpha(i,j+1,k) - alpha(i,j-1,k))
      Gaz = 1/(2*hz) * (alpha(i,j,k+1) - alpha(i,j,k-1))

      ! Lorr_rate = -1/|B|  (lamb w1 - G(lamb) * B) F/|F| 
      lrate_x = - 1/magB(i,j,k) * ( (lambda(i,j,k) * omega1 - BLG) * Fx(i,j,k) / magF(i,j,k) )
      lrate_y = - 1/magB(i,j,k) * ( (lambda(i,j,k) * omega1 - BLG) * Fy(i,j,k) / magF(i,j,k) )
      lrate_z = - 1/magB(i,j,k) * ( (lambda(i,j,k) * omega1 - BLG) * Fz(i,j,k) / magF(i,j,k) )

      Lorr_rate(i,j,k) = sqrt(lrate_x**2 + lrate_y**2 + lrate_z**2)

      ! Wind_rate = -1/|B| ( lamb (alph + w2)B_f )
      wrate_x = - 1/magB(i,j,k) * ( lambda(i,j,k) * (alpha(i,j,k) + omega2) * Bfx(i,j,k) )
      wrate_y = - 1/magB(i,j,k) * ( lambda(i,j,k) * (alpha(i,j,k) + omega2) * Bfy(i,j,k) )
      wrate_z = - 1/magB(i,j,k) * ( lambda(i,j,k) * (alpha(i,j,k) + omega2) * Bfz(i,j,k) )

      Wind_rate(i,j,k) = sqrt(wrate_x**2 + wrate_y**2 + wrate_z**2)

      ! Alig_rate = -1/|B| ( G(alph) x B )
      arate_x = - 1/magB(i,j,k) * ( Gay * Bz(i,j,k) - Gaz * By(i,j,k) )
      arate_y = - 1/magB(i,j,k) * ( Gaz * Bx(i,j,k) - Gax * Bz(i,j,k) )
      arate_z = - 1/magB(i,j,k) * ( Gax * By(i,j,k) - Gay * Bx(i,j,k) )

      Alig_rate(i,j,k) = sqrt(arate_x**2 + arate_y**2 + arate_z**2)

      rate(i,j,k) = sqrt( (arate_x+wrate_x+lrate_x)**2 &
                         +(arate_y+wrate_y+lrate_y)**2 &
                         +(arate_z+wrate_z+lrate_z)**2 )

    end do
    end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Write to file
    print *, 'Writing '//'rate_'//trim(adjustl(str_ind))//'.nc'
    call save_rates_to_netcdf(x,y,z,rate,Lorr_rate,Wind_rate,Alig_rate,'rates_'//trim(adjustl(str_ind))//'.nc')

    ! ALLOCATIONS
    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(Bx)
    deallocate(By)
    deallocate(Bz)
    deallocate(magB)
    deallocate(Fx)
    deallocate(Fy)
    deallocate(Fz)
    deallocate(magF)
    deallocate(Jx)
    deallocate(Jy)
    deallocate(Jz)
    deallocate(Bfx)
    deallocate(Bfy)
    deallocate(Bfz)
    deallocate(alpha)
    deallocate(lambda)
    deallocate(rate)
    deallocate(Lorr_rate)
    deallocate(Wind_rate)
    deallocate(Alig_rate)

  end do

contains

  subroutine save_rates_to_netcdf(x,y,z,rate,lrate,wrate,arate,filename)

    use netcdf

    real(8), intent(in), dimension(:) :: x, y, z
    real(8), intent(in), dimension(:,:,:) :: rate, lrate, wrate, arate
    character(*), intent(in) :: filename

    integer :: i
    integer :: nx, ny, nz

    integer :: id_nc,xdim_id,ydim_id,zdim_id
    integer :: x_id,y_id,z_id
    integer :: id_rate, id_lrate, id_wrate, id_arate

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

    call check(nf90_def_var(id_nc, "Total rate",   NF90_DOUBLE, dim_ids, id_rate,  deflate_level=7))
    call check(nf90_def_var(id_nc, "Lorrenz rate", NF90_DOUBLE, dim_ids, id_lrate, deflate_level=7))
    call check(nf90_def_var(id_nc, "Winding rate", NF90_DOUBLE, dim_ids, id_wrate, deflate_level=7))
    call check(nf90_def_var(id_nc, "Aligned rate", NF90_DOUBLE, dim_ids, id_arate, deflate_level=7))

    call check(nf90_enddef(id_nc))

    call check(nf90_put_var(id_nc, x_id, x))
    call check(nf90_put_var(id_nc, y_id, y))
    call check(nf90_put_var(id_nc, z_id, z))

    call check(nf90_put_var(id_nc, id_rate,   rate))
    call check(nf90_put_var(id_nc, id_lrate, lrate))
    call check(nf90_put_var(id_nc, id_wrate, wrate))
    call check(nf90_put_var(id_nc, id_arate, arate))
    
    call check(nf90_close(id_nc))

  end subroutine save_rates_to_netcdf
  
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

  subroutine check(istatus)
    use netcdf
    implicit none
    integer, intent (in) :: istatus
    if (istatus /= nf90_noerr) then
    write(*,*) trim(adjustl(nf90_strerror(istatus)))
    end if
  end subroutine check


end program main
