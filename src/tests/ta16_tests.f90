program main

  use iso_fortran_env, only: real64

  implicit none

  integer, parameter :: dp=real64

  call test_rbf_centers

contains

  subroutine test_rbf_centers
    ! Tests new rbf_centers subroutine against the original 
    ! TA16 calculation. New variables have an n appended
    ! in front. 
    
    use Original_TA16
    use TA16, only : calculate_rbf_centers

    real(dp), dimension(1296) :: nx,ny,nz,nst,nrho,nzsp,nzcp,nrhbr
    real(dp) :: hold_parmod(10), hbx,hby,hbz
    integer :: i

    hold_parmod = 1.0_dp
    call calculate_rbf_centers(nx,ny,nz,nst,nrho,nzsp,nzcp,nrhbr)
    call RBF_MODEL_2016(0,hold_parmod, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, hbx, hby, hbz)

    print *, "rbf centers subroutine errors;"
    print *, "  x -> ", sum(abs(nx-XX))
    print *, "  y -> ", sum(abs(ny-YY))
    print *, "  z -> ", sum(abs(nz-ZZ))
    print *, "  st -> ", sum(abs(nst-ST))
    print *, "  rho -> ", sum(abs(nrho-RHO))
    print *, "  zsp -> ", sum(abs(nzsp-ZSP))
    print *, "  zcp -> ", sum(abs(nzcp-ZCP))
    print *, "  rhbr -> ", sum(abs(nrhbr-RHBR))

  end subroutine test_rbf_centers

end program main
