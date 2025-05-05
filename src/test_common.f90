module testing
  implicit none
contains

  subroutine common_from_geo(iy,id,ih,im,is,vx,vy,vz)

    use geopack
    use test_old

    integer, intent(in) :: iy, id, ih, im, is
    real(8), intent(in) :: vx, vy, vz
    real(8) :: A(105)

    call RECALC_08(iy,id,ih,im,is,vx,vy,vz)
    call get_common(A)

    print *, A
  end subroutine common_from_geo

end module testing

program main

  use testing, only : common_from_geo

  implicit none
  real(8) :: vx, vy, vz
  vx = -415.1
  vy = -28.4
  vz = -27.2

  call common_from_geo(2020,50,12,00,00,vx,vy,vz)

end program main



  
