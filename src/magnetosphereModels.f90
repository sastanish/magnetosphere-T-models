module MagnetosphereModels

  implicit none

  real(8), allocatable, dimension(:,:,:) :: x, y, z
  real(8), allocatable, dimension(:,:,:) :: Bfpx, Bfpy, Bfpz, bx, by, bz
  real(8), allocatable, dimension(:,:,:) :: Fx, Fy, Fz, jx, jy, jz
  real(8), allocatable, dimension(:,:,:) :: c2_t1_x, c2_t1_y, c2_t1_z
  real(8), allocatable, dimension(:,:,:) :: c2_t2_x, c2_t2_y, c2_t2_z
  real(8), allocatable, dimension(:,:,:) :: c2_t3_x, c2_t3_y, c2_t3_z
  real(8), allocatable, dimension(:,:,:) :: alpha, lambda

  public :: run_models, cleanup

  contains

    ! First function that gets called from main program
    subroutine run_models(xin,yin,zin,parmod,ps,controlParams,outFileName)
      real(8), dimension(:), intent(in) :: xin, yin, zin
      real(8), dimension(:), intent(in) :: parmod
      real(8), intent(in) :: ps
      integer, dimension(:), intent(in) :: controlParams
      character(:), intent(in) :: outFileName

      allocate(x(SIZE(xin)))
      allocate(y(SIZE(yin)))
      allocate(z(SIZE(zin)))

      x = xin
      y = yin
      z = zin

      call init_data_vars(controlParams(4:))
      call compute_field(parmod,ps,controlParams)
      if (sum(controlParams(4:)) > 10) call compute_secondary_vars
      if (sum(controlParams(4:)) > 200) call compute_reconnection_metrics
      call write_data(controlParams,outFileName)

    end subroutine run_models

    subroutine init_data_vars(target_vars)

      integer, intent(in) :: nx,ny,nz
      integer, intent(in) :: target_vars(:)
      integer :: id

      do id = 1,size(target_vars)
        if (target_vars(id) == 1) allocate(bx(nx,ny,nz))
        if (target_vars(id) == 2) allocate(by(nx,ny,nz))
        if (target_vars(id) == 3) allocate(bz(nx,ny,nz))
        if (target_vars(id) == 11) allocate(jx(nx,ny,nz))
        if (target_vars(id) == 12) allocate(jy(nx,ny,nz))
        if (target_vars(id) == 13) allocate(jz(nx,ny,nz))
        if (target_vars(id) == 21) allocate(fx(nx,ny,nz))
        if (target_vars(id) == 22) allocate(fy(nx,ny,nz))
        if (target_vars(id) == 23) allocate(fz(nx,ny,nz))
        if (target_vars(id) == 201) allocate(c2_t1_x(nx,ny,nz))
        if (target_vars(id) == 202) allocate(c2_t1_y(nx,ny,nz))
        if (target_vars(id) == 203) allocate(c2_t1_z(nx,ny,nz))
        if (target_vars(id) == 301) allocate(c2_t2_x(nx,ny,nz))
        if (target_vars(id) == 302) allocate(c2_t2_y(nx,ny,nz))
        if (target_vars(id) == 303) allocate(c2_t2_z(nx,ny,nz))
        if (target_vars(id) == 401) allocate(c2_t3_x(nx,ny,nz))
        if (target_vars(id) == 402) allocate(c2_t3_y(nx,ny,nz))
        if (target_vars(id) == 403) allocate(c2_t3_z(nx,ny,nz))
      end do

    end subroutine init_data_vars

    subroutine cleanup

      ! free all allocated arrays

    end subroutine cleanup

end module magnetosphereModels
