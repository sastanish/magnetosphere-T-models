module MagnetosphereModels

  implicit none

  real(8), allocatable, dimension(:,:,:) :: x, y, z
  real(8), allocatable, dimension(:,:,:) :: Bfpx, Bfpy, Bfpz, bx, by, bz
  real(8), allocatable, dimension(:,:,:) :: Fx, Fy, Fz, jx, jy, jz
  real(8), allocatable, dimension(:,:,:) :: c2_t1_x, c2_t1_y, c2_t1_z
  real(8), allocatable, dimension(:,:,:) :: c2_t2_x, c2_t2_y, c2_t2_z
  real(8), allocatable, dimension(:,:,:) :: c2_t3_x, c2_t3_y, c2_t3_z

  public :: run_models

  contains

    ! First function that gets called from main program
    subroutine run_models(inData,controlParams,gridDims)
      real(8), dimension(:), intent(in) :: inData
      real(8), dimension(:), intent(in) :: gridDims
      integer, dimension(:), intent(in) :: controlParams

      call


    end subroutine run_models

    subroutine init_coordinates(nx,ny,nz,xp,xm,yp,ym,zp,zm)

      integer, intent(in) :: nx,ny,nz
      real(8), intent(in) :: xp,xm,yp,ym,zp,zm
      integer :: i

      allocate(x(nx))
      allocate(y(ny))
      allocate(z(nz))

      do i = 1,nx
        x(i) = xm + ((xp-xm)/nx) * (i-1)
      end do
      do i = 1,ny
        y(i) = ym + ((yp-ym)/ny) * (i-1)
      end do
      do i = 1,nz
        z(i) = zm + ((zp-zm)/nz) * (i-1)
      end do

    end subroutine init_coordinates

    subroutine init_data_vars(nx,ny,nz,target_vars)

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
        if (target_vars(id) == 101) allocate(c2_t1_x(nx,ny,nz))
        if (target_vars(id) == 102) allocate(c2_t1_y(nx,ny,nz))
        if (target_vars(id) == 103) allocate(c2_t1_z(nx,ny,nz))
        if (target_vars(id) == 201) allocate(c2_t2_x(nx,ny,nz))
        if (target_vars(id) == 202) allocate(c2_t2_y(nx,ny,nz))
        if (target_vars(id) == 203) allocate(c2_t2_z(nx,ny,nz))
        if (target_vars(id) == 301) allocate(c2_t3_x(nx,ny,nz))
        if (target_vars(id) == 302) allocate(c2_t3_y(nx,ny,nz))
        if (target_vars(id) == 303) allocate(c2_t3_z(nx,ny,nz))
      end do

    end subroutine init_data_vars

end module magnetosphereModels
