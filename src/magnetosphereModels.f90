module magnetosphereModels

  implicit none

contains

  subroutine compute_field(x,y,z,dateInfo,velocity,parmod,ps,modelNumber,dipoleNumber,bx,by,bz)

    use TA16, only : RBF_MODEL_2016
    use TS05, only : T04_s
    use geopack, only : RECALC_08, DIP_08, IGRF_GSW_08

    real(8), dimension(:), intent(in) :: x,y,z

    integer, dimension(:), intent(in) :: dateInfo ! [year, dayNumber, hour, min, sec]
    real(8), dimension(:), intent(in) :: velocity ! GSW Solar wind components [vx, vy, vz]
    real(8), dimension(:), intent(in) :: parmod
    real(8), intent(in) :: ps ! Geo-dipole tilt angle (in radians)
    integer, intent(in) :: modelNumber, dipoleNumber
    ! control parameters, modelNumber = 5 (TS05), 16 (TA16)
    !                     dipoleNumber = 1 (DIP_08), 2(IGRF_GSW_08)
    
    real(8), intent(out), dimension(:,:,:) :: bx,by,bz

    ! dummy vars
    real(8) :: xx,yy,zz,ebx,eby,ebz,ibx,iby,ibz
    integer :: i, j, k, nx, ny, nz

    integer :: iyear, iday, ihour, imin, isec
    real(8) :: vx, vy, vz

    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)

    iyear = dateInfo(1)
    iday = dateInfo(2)
    ihour = dateInfo(3)
    isec = dateInfo(4)

    vx = velocity(1)
    vy = velocity(2)
    vz = velocity(3)
  
    call RECALC_08(iyear, iday, ihour, imin, isec, vx, vy, vz)
    do k = 1,nz
      zz = z(k)
      do j = 1,ny
        yy = y(j)
        do i = 1,nx
          xx = x(i)

          if ( dipoleNumber == 1 ) call DIP_08(xx,yy,zz,ibx,iby,ibz)
          if ( dipoleNumber == 2 ) call IGRF_GSW_08(xx,yy,zz,ibx,iby,ibz)
          if ( modelNumber == 5  ) call T04_s(0,parmod,ps,xx,yy,zz,ebx,eby,ebz)
          if ( modelNumber == 16 ) call RBF_MODEL_2016(0,parmod,ps,xx,yy,zz,ebx,eby,ebz)

          bx(i,j,k) = ibx+ebx
          by(i,j,k) = iby+eby
          bz(i,j,k) = ibz+ebz

        end do
      end do
    end do

  end subroutine compute_field

end module magnetosphereModels
