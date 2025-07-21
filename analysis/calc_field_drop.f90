program main

  use inputOutput, only : load_field_from_netcdf

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

    call load_field_from_netcdf(trim(dir)//'output_'//trim(adjustl(str_ind))//'.nc',x,y,z,bx,by,bz)
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

end program main
