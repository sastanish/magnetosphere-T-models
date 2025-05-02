program main

  use parser

  implicit none
  integer :: nx, ny, nz

  !$OMP PARALLEL DO
  do lineNumber = 1,nLines
    call read_data_line(lineNumber,dataFileName,inData)
    call run_models(inData,controlParams,gridDims)
  end do
  !$OMP END PARALLEL DO



  call final_cleanup()

end program main
