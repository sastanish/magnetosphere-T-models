program main

  use iso_fortran_env, only: real64

  implicit none

  integer, parameter :: dp=real64
  real(dp) :: out_error

  call test_rbf_centers(out_error)
  print *, "rbf centers subroutine, abs error;"
  print *, "     ", out_error

contains

  subroutine test_rbf_centers(errors)
    ! Tests new rbf_centers subroutine against the original 
    ! TA16 calculation. New variables have an n appended
    ! in front. 
    !
    ! Returns number of mismatched values

    use Original_TA16
    use TA16, only : calculate_rbf_centers

    real(dp), dimension(1296) :: nx,ny,nz,nst,nrho,nzsp,nzcp,nrhbr
    real(dp) :: hold_parmod(10), hbx,hby,hbz
    real(dp), intent(out) :: errors
    integer :: i

    hold_parmod = 1.0_dp
    call calculate_rbf_centers(nx,ny,nz,nst,nrho,nzsp,nzcp,nrhbr)
    call RBF_MODEL_2016(0,hold_parmod, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, hbx, hby, hbz)

    do i = 1,size(nx)
      print *, nx(i), " ", XX(i)
    end do

    errors = 0.0
    errors = errors + sum(abs(nx-XX))
    !errors = errors + sum(abs(ny-YY))
    !errors = errors + sum(abs(nz-ZZ))
    !errors = errors + sum(abs(nst-ST))
    !errors = errors + sum(abs(nrho-RHO))

  end subroutine test_rbf_centers

end program main
