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

  contains
    subroutine compute_metrics(x,y,z,bx,by,bz,Mout)
      
      use reconnectionMetrics, only : compute_parameters, compute_forces,&
                                      compute_c2_t1, compute_c2_t2, compute_c2_t3

      real(8), intent(in), dimension(:) :: x,y,z
      real(8), intent(in), dimension(:,:,:) :: bx,by,bz

      real(8), intent(out), dimension(:,:,:,:), allocatable :: Mout
      ! First dimension is [jx, jy, jz, fx, fy, fz, bfpx, bfpy, bfpz,&
      !                     alpha, lambda, c2_t1, c2_t2, c2_t3]

      call compute_forces(x,y,z,bx,by,bz,Mout(1,:,:,:),&
                                         Mout(2,:,:,:),&
                                         Mout(3,:,:,:),&
                                         Mout(4,:,:,:),&
                                         Mout(5,:,:,:),&
                                         Mout(6,:,:,:) &
                                         )

      call compute_parameters(x,y,z,bx,by,bz,Mout(1,:,:,:),&
                                             Mout(2,:,:,:),&
                                             Mout(3,:,:,:),&
                                             Mout(4,:,:,:),&
                                             Mout(5,:,:,:),&
                                             Mout(6,:,:,:),&
                                             Mout(7,:,:,:),&
                                             Mout(8,:,:,:),&
                                             Mout(9,:,:,:),&
                                             Mout(10,:,:,:),&
                                             Mout(11,:,:,:)&
                                             )

      call compute_c2_t1(x,y,z,bx,by,bz,Mout(4,:,:,:),&
                                        Mout(5,:,:,:),&
                                        Mout(6,:,:,:),&
                                        Mout(7,:,:,:),&
                                        Mout(8,:,:,:),&
                                        Mout(9,:,:,:),&
                                        Mout(11,:,:,:),&
                                        Mout(12,:,:,:)&
                                        )

      call compute_c2_t2(x,y,z,bx,by,bz,Mout(4,:,:,:),&
                                        Mout(5,:,:,:),&
                                        Mout(6,:,:,:),&
                                        Mout(7,:,:,:),&
                                        Mout(8,:,:,:),&
                                        Mout(9,:,:,:),&
                                        Mout(10,:,:,:),&
                                        Mout(11,:,:,:),&
                                        Mout(13,:,:,:)&
                                        )

      call compute_c2_t3(x,y,z,bx,by,bz,Mout(10,:,:,:),Mout(14,:,:,:))
      
    end subroutine compute_metrics



end program main
