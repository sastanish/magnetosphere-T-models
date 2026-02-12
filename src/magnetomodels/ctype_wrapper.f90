module magnetomodels
  use iso_c_binding
  implicit none
  contains
subroutine TA16(date,data,x,y,z,Bx,By,Bz,nx,ny,nz) bind(c, name='ta16')

  use control, only : calculate_TA16, setup_TA16_model

  !! Input !!
  integer(c_int), intent(in) :: nx, ny, nz
  real(c_double), intent(in) :: x(nx), y(ny), z(nz)
  real(c_double), intent(in) :: data(8)
  integer(c_int), intent(in) :: date(4)

  !! Output !!
  real(c_double), intent(out) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)

  !! Vars for exploding arrays !!
  integer :: year,day,hour,mint
  real(8) :: ivx,ivy,ivz,tilt,pydn,symhc,nind,aby

  year = date(1)
  day = date(2)
  hour = date(3)
  mint = date(4)
  ivx = data(1)
  ivy = data(2)
  ivz = data(3)
  tilt = data(4)
  pydn = data(5)
  symhc = data(6)
  nind = data(7)
  aby = data(8)

  call setup_TA16_model()
  call calculate_TA16(year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symhc,nind,aby,x,y,z,Bx,By,Bz)
  
end subroutine TA16

end module magnetomodels
