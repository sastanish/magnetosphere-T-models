! Module compute defined in file python_wrapper.f90

subroutine f90wrap_compute__ta16(f90wrap_n0, f90wrap_n1, f90wrap_n2, f90wrap_n3, f90wrap_n4, f90wrap_n5, f90wrap_n6, &
    f90wrap_n7, f90wrap_n8, f90wrap_n9, f90wrap_n10, f90wrap_n11, date, data, x, y, z, bx, by, bz, nx, ny, nz)
    use compute
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    
    integer(c_int) :: f90wrap_n0
    !f2py intent(hide), depend(x) :: f90wrap_n0 = shape(x,0)
    integer(c_int) :: f90wrap_n1
    !f2py intent(hide), depend(y) :: f90wrap_n1 = shape(y,0)
    integer(c_int) :: f90wrap_n2
    !f2py intent(hide), depend(z) :: f90wrap_n2 = shape(z,0)
    integer(c_int) :: f90wrap_n3
    !f2py intent(hide), depend(bx) :: f90wrap_n3 = shape(bx,0)
    integer(c_int) :: f90wrap_n4
    !f2py intent(hide), depend(bx) :: f90wrap_n4 = shape(bx,1)
    integer(c_int) :: f90wrap_n5
    !f2py intent(hide), depend(bx) :: f90wrap_n5 = shape(bx,2)
    integer(c_int) :: f90wrap_n6
    !f2py intent(hide), depend(by) :: f90wrap_n6 = shape(by,0)
    integer(c_int) :: f90wrap_n7
    !f2py intent(hide), depend(by) :: f90wrap_n7 = shape(by,1)
    integer(c_int) :: f90wrap_n8
    !f2py intent(hide), depend(by) :: f90wrap_n8 = shape(by,2)
    integer(c_int) :: f90wrap_n9
    !f2py intent(hide), depend(bz) :: f90wrap_n9 = shape(bz,0)
    integer(c_int) :: f90wrap_n10
    !f2py intent(hide), depend(bz) :: f90wrap_n10 = shape(bz,1)
    integer(c_int) :: f90wrap_n11
    !f2py intent(hide), depend(bz) :: f90wrap_n11 = shape(bz,2)
    integer, intent(in), dimension(4) :: date
    real(8), intent(in), dimension(8) :: data
    real(8), intent(in), dimension(f90wrap_n0) :: x
    real(8), intent(in), dimension(f90wrap_n1) :: y
    real(8), intent(in), dimension(f90wrap_n2) :: z
    real(8), intent(inout), dimension(f90wrap_n3,f90wrap_n4,f90wrap_n5) :: bx
    real(8), intent(inout), dimension(f90wrap_n6,f90wrap_n7,f90wrap_n8) :: by
    real(8), intent(inout), dimension(f90wrap_n9,f90wrap_n10,f90wrap_n11) :: bz
    integer(c_int), intent(in) :: nx
    integer(c_int), intent(in) :: ny
    integer(c_int), intent(in) :: nz
    integer :: f90wrap_f90wrap_n0_default
    integer :: f90wrap_f90wrap_n1_default
    integer :: f90wrap_f90wrap_n2_default
    integer :: f90wrap_f90wrap_n3_default
    integer :: f90wrap_f90wrap_n4_default
    integer :: f90wrap_f90wrap_n5_default
    integer :: f90wrap_f90wrap_n6_default
    integer :: f90wrap_f90wrap_n7_default
    integer :: f90wrap_f90wrap_n8_default
    integer :: f90wrap_f90wrap_n9_default
    integer :: f90wrap_f90wrap_n10_default
    integer :: f90wrap_f90wrap_n11_default
    integer :: f90wrap_nx_default
    integer :: f90wrap_ny_default
    integer :: f90wrap_nz_default
    f90wrap_f90wrap_n0_default = f90wrap_n0
    f90wrap_f90wrap_n1_default = f90wrap_n1
    f90wrap_f90wrap_n2_default = f90wrap_n2
    f90wrap_f90wrap_n3_default = f90wrap_n3
    f90wrap_f90wrap_n4_default = f90wrap_n4
    f90wrap_f90wrap_n5_default = f90wrap_n5
    f90wrap_f90wrap_n6_default = f90wrap_n6
    f90wrap_f90wrap_n7_default = f90wrap_n7
    f90wrap_f90wrap_n8_default = f90wrap_n8
    f90wrap_f90wrap_n9_default = f90wrap_n9
    f90wrap_f90wrap_n10_default = f90wrap_n10
    f90wrap_f90wrap_n11_default = f90wrap_n11
    f90wrap_nx_default = nx
    f90wrap_ny_default = ny
    f90wrap_nz_default = nz
    call TA16(date=date, data=data, x=x, y=y, z=z, Bx=bx, By=by, Bz=bz, nx=f90wrap_nx_default, ny=f90wrap_ny_default, &
        nz=f90wrap_nz_default)
    f90wrap_n0 = f90wrap_f90wrap_n0_default
    f90wrap_n1 = f90wrap_f90wrap_n1_default
    f90wrap_n2 = f90wrap_f90wrap_n2_default
    f90wrap_n3 = f90wrap_f90wrap_n3_default
    f90wrap_n4 = f90wrap_f90wrap_n4_default
    f90wrap_n5 = f90wrap_f90wrap_n5_default
    f90wrap_n6 = f90wrap_f90wrap_n6_default
    f90wrap_n7 = f90wrap_f90wrap_n7_default
    f90wrap_n8 = f90wrap_f90wrap_n8_default
    f90wrap_n9 = f90wrap_f90wrap_n9_default
    f90wrap_n10 = f90wrap_f90wrap_n10_default
    f90wrap_n11 = f90wrap_f90wrap_n11_default
end subroutine f90wrap_compute__ta16

! End of module compute defined in file python_wrapper.f90

