module TsyganenkoWrapper

  implicit none

contains
  subroutine field(x,y,z,dateInfo,velocity,parmod,ps,modelNumber,dipoleNumber,bx,by,bz,nx,ny,nz)

    use magnetosphereModels, only: compute_field

    integer :: nx, ny, nz

    real(8), intent(in) :: x(nx), y(ny), z(nz)
    integer, intent(in) :: dateinfo(5)
    real(8), intent(in) :: velocity(3)
    real(8), intent(in) :: parmod(10)
    real(8), intent(in) :: ps
    integer, intent(in) :: modelNumber, dipoleNumber

    real(8), dimension(nx,ny,nz), intent(out) :: bx, by, bz

    call compute_field(x,y,z,dateInfo,velocity,parmod,ps,modelNumber,dipoleNumber,bx,by,bz)

  end subroutine field

  subroutine metrics(x,y,z,bx,by,bz,Mout,nx,ny,nz)

    use reconnectionMetrics, only: compute_all_metrics

    integer :: nx, ny, nz

    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), dimension(nx,ny,nz), intent(in) :: bx, by, bz
    real(8), intent(out) :: Mout(14,nx,ny,nz)

    call compute_all_metrics(x,y,z,bx,by,bz,Mout)

  end subroutine metrics

end module TsyganenkoWrapper
