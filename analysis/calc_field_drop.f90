program main

  use netcdf

  implicit none

  !Data vars
  real(8), dimension(:,:,:), allocatable :: Bx, By, Bz
  real(8), dimension(:), allocatable :: x,y,z,strength
  real(8) :: Area
  integer :: year, day, hour, minn

  !Indices
  integer :: f,i,j,k,eind

  !File handling
  integer :: fileind, start_ind, end_ind, outfile, omnifile
  character(4) :: str_ind, start_str, end_str
  character(50) :: dir

  call GET_COMMAND_ARGUMENT(1,dir)
  call GET_COMMAND_ARGUMENT(2,start_str)
  call GET_COMMAND_ARGUMENT(3,end_str)
  read(start_str,'(I5)') start_ind
  read(end_str,'(I5)') end_ind

  allocate(strength(start_ind:end_ind))

  do f = start_ind,end_ind

    write( str_ind, '(I4)' ) f
    print *, 'reading: '//trim(dir)//'*_'//trim(adjustl(str_ind))//'.nc'

    call load_from_netcdf(trim(dir)//'output_'//trim(adjustl(str_ind))//'.nc',x,y,z,bx,by,bz)
    Area = size(z)*size(y)*eind
    do i = 1,size(x)
      if ( x(i) > -2) then
        eind = i
        exit
      end if
    end do

    strength(f) = 0
    do k=1,size(z)
    do j=1,size(y)
    do i=1,eind
      strength(f) = strength(f) + sqrt(bx(i,j,k)**2 + by(i,j,k)**2 + bz(i,j,k)**2)/Area 
    end do
    end do
    end do

    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(bx)
    deallocate(by)
    deallocate(bz)
  end do

  open(newunit=outfile, file=trim(dir)//"field_strength.lst", status="new")
  open(newunit=omnifile, file=trim(dir)//"input_data.lst", status="old")
  write(outfile,"(A)") "|B| in magnetotail where x < -2."
  write(outfile,"(A)") "FORMAT:"
  write(outfile,"(A)") "YEAR    DAY    HOUR    MIN    |B|"

  ! skip omni file to correct line
  do i = 1,start_ind
    read(omnifile,*)
  end do
  print *, strength
  do i = start_ind,end_ind
    read(omnifile,*) year, day, hour, minn
    write(outfile,'(I4, 4x, I3, 4x, I4, 4x, I3, 4x, F8.4)') year, day, hour, minn, strength(i)
  end do
  close(outfile)
  close(omnifile)
  deallocate(strength)

contains

  subroutine load_from_netcdf(filename,x,y,z,Bx,By,Bz)

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

  end subroutine load_from_netcdf
 
  subroutine check(istatus)
    use netcdf
    implicit none
    integer, intent (in) :: istatus
    if (istatus /= nf90_noerr) then
    write(*,*) trim(adjustl(nf90_strerror(istatus)))
    end if
  end subroutine check


end program main
