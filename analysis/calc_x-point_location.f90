program main

  use inputOutput, only : load_rate_from_netcdf, load_field_from_netcdf

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
    yind = 0
    do i=1,size(y)
      if (abs(y(i)) .le. 1e-3) yind=i
    end do

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

end program main
