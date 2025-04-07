program main
  use reconection_metrics

  implicit none
  character(*) :: in_file

  in_file = "test.nc"
  call calculate_metrics(in_file)

end program main
