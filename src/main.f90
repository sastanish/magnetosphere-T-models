program main

  use magnetosphereModels, only : run_models
  use parser

  implicit none
  integer :: nx, ny, nz

  ! Get a command line argument for controlFileName
  call parse_input_file(controlFileName,controlParams,gridDims)
  ! controlParams holds all of our control info about the sim inside it.

  !$OMP PARALLEL DO
  do lineNumber = 1,nLines
    call read_data_line(lineNumber,dataFileName,inData)
    call run_models(inData,controlParams,gridDims)
  end do
  !$OMP END PARALLEL DO

  call final_cleanup()

end program main
