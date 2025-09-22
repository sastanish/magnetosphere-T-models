program main

  use iso_fortran_env, only: real64

  implicit none

  integer, parameter :: dp=real64

  call test_rbf_centers
  call test_ta16_calc

contains

  subroutine test_ta16_calc

    use Original_TA16
    use TA16

    real(dp) :: hold_parmod(10), hbx,hby,hbz, nbx, nby, nbz
    integer :: i

    hold_parmod(1) = 1.82_dp
    hold_parmod(2) = -16.5_dp
    hold_parmod(3) = 0.0026_dp
    hold_parmod(4) = 0.62_dp

    call read_ta16_parameters('TA16_RBF.par')
    call calculate_ta16_field(hold_parmod(1), hold_parmod(2), hold_parmod(3), hold_parmod(4), 0.0904_dp, -1.0_dp, 1.0_dp, 0.0_dp, nbx, nby, nbz)
    call RBF_MODEL_2016(0,hold_parmod, 0.0904_dp, -1.0_dp, 1.0_dp, 0.0_dp, hbx, hby, hbz)

    print *, "ta16 field calculation errors"
    print *, "  bx -> ", abs(nbx-hbx)
    print *, "  by -> ", abs(nby-hby)
    print *, "  bz -> ", abs(nbz-hbz)
    print *, "ta16 field calculation Mags"
    print *, "old bx -> ", abs(hbx), "new bx -> ", abs(nbx)
    print *, "old by -> ", abs(hby), "new by -> ", abs(nby)
    print *, "old bz -> ", abs(hbz), "new bz -> ", abs(nbz)

  end subroutine test_ta16_calc

  subroutine test_rbf_centers
    ! Tests new rbf_centers subroutine against the original 
    ! TA16 calculation. New variables have an n appended
    ! in front. 
    
    use Original_TA16
    use TA16

    real(dp) :: hold_parmod(10), hbx,hby,hbz
    integer :: i

    hold_parmod = 1.0_dp
    call calculate_rbf_centers
    call RBF_MODEL_2016(0,hold_parmod, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, hbx, hby, hbz)

    print *, "rbf centers subroutine errors;"
    print *, "  x -> ", sum(abs(xrb-XX))
    print *, "  y -> ", sum(abs(yrb-YY))
    print *, "  z -> ", sum(abs(zrb-ZZ))
    print *, "  st -> ", sum(abs(strb-ST))
    print *, "  rho -> ", sum(abs(rhorb-RHO))
    print *, "  zsp -> ", sum(abs(zsprb-ZSP))
    print *, "  zcp -> ", sum(abs(zcprb-ZCP))
    print *, "  rhbr -> ", sum(abs(rhbrrb-RHBR))

  end subroutine test_rbf_centers

end program main
