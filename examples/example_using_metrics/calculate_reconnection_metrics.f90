program main
  use reconnection_metrics

  implicit none
  character(len=50) :: in_file
  integer :: io, i
  
  ! loops through a list of files to compute over.  
  ! generate via `ls ./data > file_list.txt`
  open(unit=1,file="file_list.txt",status="OLD",action="READ")
  do
    read(1, '(A)', end=200) in_file
    call calculate_metrics(trim(in_file))
  end do
  200 continue

end program main
