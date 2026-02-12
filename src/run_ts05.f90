program main

  use control, only : setup_grid,  read_TS05_input_data, calculate_TS05, progress_bar
  use inputOutput, only : save_field_to_netcdf

  implicit none

  real(8), dimension(:,:,:), allocatable :: Bx, By, Bz
  real(8), dimension(:), allocatable :: x,y,z

  ! Input File params
  integer, dimension(:), allocatable :: year,day,hour,mint

  ! TS05 params
  real(8), dimension(:), allocatable :: ivx,ivy,ivz,tilt,pydn,symh,iby,ibz,w1,w2,w3,w4,w5,w6

  ! vars
  integer :: n
  character(4) :: str_ind

  ! Command line arguments
  integer :: start_ind, end_ind
  character(4) :: start_str, end_str

  ! Convert command line args to integers for control
  call GET_COMMAND_ARGUMENT(1,start_str)
  call GET_COMMAND_ARGUMENT(2,end_str)
  read(start_str,'(I5)') start_ind
  read(end_str,'(I5)') end_ind

  call setup_grid('input_parameters.txt',x,y,z,Bx,By,Bz)
  call read_TS05_input_data('input_data.lst',year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symh,iby,ibz,w1,w2,w3,w4,w5,w6)

  do n = start_ind,end_ind

    call calculate_TS05(year(n),day(n),hour(n),mint(n),ivx(n),ivy(n),ivz(n),tilt(n),pydn(n),symh(n),iby(n),ibz(n),w1(n),w2(n),w3(n),w4(n),w5(n),w6(n),x,y,z,Bx,By,Bz)

    write( str_ind, '(I4)' ) n
    call save_field_to_netcdf(x,y,z,Bx,By,Bz,'output_'//trim(adjustl(str_ind))//'.nc')
    call progress_bar(n,end_ind-start_ind,'output_'//trim(adjustl(str_ind))//'.nc')

  end do

  deallocate(Bx)
  deallocate(By)
  deallocate(Bz)
  deallocate(x)
  deallocate(y)
  deallocate(z)

end program main
