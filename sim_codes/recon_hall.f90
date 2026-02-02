program main

  !$ use omp_lib

  implicit none

  real(8), dimension(:,:,:), allocatable :: Bx,By,Bz,magB
  real(8), dimension(:,:,:), allocatable :: Jx,Jy,Jz,Bex,Bey,Bez,JxB_x,JxB_y,JxB_z
  real(8), dimension(:,:,:), allocatable :: rate

  real(8), dimension(:), allocatable :: x,y,z

  integer :: i,j,k,nx,ny,nz
  real(8) :: hx,hy,hz
  real(8) :: rate_x,rate_y,rate_z,cJxB_x,cJxB_y,cJxB_z,JxB_dot_b

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
    allocate(rate(nx,ny,nz))
    allocate(magB(nx,ny,nz))
    allocate(Jx(nx,ny,nz))
    allocate(Jy(nx,ny,nz))
    allocate(Jz(nx,ny,nz))
    allocate(Bex(nx,ny,nz))
    allocate(Bey(nx,ny,nz))
    allocate(Bez(nx,ny,nz))
    allocate(JxB_x(nx,ny,nz))
    allocate(JxB_y(nx,ny,nz))
    allocate(JxB_z(nx,ny,nz))

    rate = 0.0

    !!!!!!!!!!!!!!!!!!
    ! Calculate rate !
    !!!!!!!!!!!!!!!!!!

    hx = abs(x(2)-x(1))
    hy = abs(y(2)-y(1))
    hz = abs(z(2)-z(1))
    !$OMP Parallel shared(Bx,By,Bz,magB,Jx,Jy,Jz,Bex,Bey,Bez,JxB_x,JxB_y,JxB_z) private(i,j,k)
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

      ! B_e = B/|B|
      Bex(i,j,k) = By(i,j,k)/magB(i,j,k)
      Bey(i,j,k) = Bz(i,j,k)/magB(i,j,k)
      Bez(i,j,k) = Bx(i,j,k)/magB(i,j,k)

      ! JxB = J x B
      JxB_x(i,j,k) = Jy(i,j,k)*Bz(i,j,k) - Jz(i,j,k)*By(i,j,k)
      JxB_y(i,j,k) = Jz(i,j,k)*Bx(i,j,k) - Jx(i,j,k)*Bz(i,j,k)
      JxB_z(i,j,k) = Jx(i,j,k)*By(i,j,k) - Jy(i,j,k)*Bx(i,j,k)

    end do
    end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP Parallel shared(Bx,By,Bz,magB,Jx,Jy,Jz,Bex,Bey,Bez,JxB_x,JxB_y,JxB_z,rate) private(i,j,k,cJxB_x,cJxB_y,cJxB_z,rate_x,rate_y,rate_z,JxB_dot_b)
    !$OMP DO COLLAPSE(3)
    do k=3,nz-2
    do j=3,ny-2
    do i=3,nx-2

      ! Curl(JxB)
      cJxB_x = 1/(2*hy) * (JxB_z(i,j+1,k) - JxB_z(i,j-1,k)) - 1/(2*hz) * (JxB_y(i,j,k+1) - JxB_y(i,j,k-1))
      cJxB_y = 1/(2*hz) * (JxB_x(i,j,k+1) - JxB_x(i,j,k-1)) - 1/(2*hx) * (JxB_z(i+1,j,k) - JxB_z(i-1,j,k))
      cJxB_z = 1/(2*hx) * (JxB_y(i+1,j,k) - JxB_y(i-1,j,k)) - 1/(2*hy) * (JxB_x(i,j+1,k) - JxB_x(i,j-1,k))

      JxB_dot_b = Bex(i,j,k)*cJxB_x + Bey(i,j,k)*cJxB_y + Bez(i,j,k)*cJxB_z

      rate_x = - (cJxB_x - JxB_dot_b * Bex(i,j,k))/magB(i,j,k)
      rate_y = - (cJxB_y - JxB_dot_b * Bey(i,j,k))/magB(i,j,k)
      rate_z = - (cJxB_z - JxB_dot_b * Bez(i,j,k))/magB(i,j,k)

      rate(i,j,k) = sqrt( (rate_x)**2 &
                         +(rate_y)**2 &
                         +(rate_z)**2 )

    end do
    end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Write to file
    print *, 'Writing '//'rate_'//trim(adjustl(str_ind))//'.nc'
    call save_rate_to_netcdf(x,y,z,rate,'hall_rate_'//trim(adjustl(str_ind))//'.nc')

    ! ALLOCATIONS
    deallocate(rate)
    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(Bx)
    deallocate(By)
    deallocate(Bz)
    deallocate(magB)
    deallocate(Jx)
    deallocate(Jy)
    deallocate(Jz)
    deallocate(Bex)
    deallocate(Bey)
    deallocate(Bez)
    deallocate(JxB_x)
    deallocate(JxB_y)
    deallocate(JxB_z)

  end do

contains

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

    call check(nf90_def_var(id_nc, "Hall rate",   NF90_DOUBLE, dim_ids, id_rate,  deflate_level=7))

    call check(nf90_enddef(id_nc))

    call check(nf90_put_var(id_nc, x_id, x))
    call check(nf90_put_var(id_nc, y_id, y))
    call check(nf90_put_var(id_nc, z_id, z))

    call check(nf90_put_var(id_nc, id_rate,   rate))
    
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

  subroutine check(istatus)
    use netcdf
    implicit none
    integer, intent (in) :: istatus
    if (istatus /= nf90_noerr) then
    write(*,*) trim(adjustl(nf90_strerror(istatus)))
    end if
  end subroutine check


end program main
