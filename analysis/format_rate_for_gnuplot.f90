program main

  use netcdf

  implicit none

  integer i,k,ofile_id,yind
  real(8), allocatable, dimension(:,:,:) :: rate
  real(8), allocatable, dimension(:) :: x,y,z

  character(50) ifilename

  call GET_COMMAND_ARGUMENT(1,ifilename)
  call load_rate_from_netcdf(trim(ifilename),x,y,z,rate)

  do i = 1,size(y)
    if (abs(y(i)) <= 0.01) then
      yind = i
      exit
    end if
  end do

  open(newunit=ofile_id,file="rate_table.lst")
  write(ofile_id,"('x', 9x, 'z', 9x, 'rate')")
  do k = 3,size(z)-2
  do i = 3,size(x)-2
    if (sqrt(x(i)**2 + z(k)**2) <= 1) then
      rate(i,yind,k) = -1
    end if
    write(ofile_id,"(F8.2, 4x, F8.2, 4x, F12.2)") x(i), z(k), rate(i,yind,k)
  end do 
  write(ofile_id,*)
  end do 
  close(ofile_id)

  deallocate(x)
  deallocate(z)
  deallocate(rate)

contains

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

