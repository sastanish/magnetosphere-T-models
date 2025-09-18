module TA16

implicit none

integer, parameter :: dp=selected_real_kind(1.d0)

real(dp) :: A(23328)

contains

subroutine calculate_rbf_centers

  integer :: klat, nlat

  real(dp) :: rlow_grid, rhigh_grid
  real(dp) :: dlat_deg, dlon_deg, xcolatd
  real(dp) :: Radius

  real(dp), parameter :: rh=8.0,alpha=3.0,pd_tr=0.5
  real(dp), parameter :: psi=0.0,pm=0.0
  real(dp), parameter :: pi=4.d0*atan(1.d0)

  real(dp), parameter, dimension(22) :: An = (12.544D0,-0.194D0,0.305D0,0.0573D0,2.178D0,0.0571D0,-0.999D0,16.473D0,0.00152D0,0.382D0,0.0431D0,-0.00763D0,-0.210D0,0.0405D0,-4.430D0,-0.636D0,-2.600D0,0.832D0,-5.328D0,1.103D0,-0.907D0,1.450D0)

  integer :: i,j,k,l

  klat = 8
  nlat = klat + 1

  rlow_grid = 3.3
  rhigh_grid = 16.0

  dlat_deg = 90.0/(nlat-0.5)

  l = 0

  Radius = rlow_grid

  An(10:14) = 0.0 ! Simplified version of Lin et al model

  
  do j = 1,100

    do i = 1,nlat

      nlon = 4*(i-1)

      if (i /= 1) then
        dlon_deg = 360.0/nlon
      else
        nlon = 1
        dlon_deg = 0.0
      end if

      do k = 1,nlon

        ! Coordinate transform
        xlond = (k-1)*dlon_deg

        xcolat = (2.0*pi/360.0) * dlat_deg * (i-1)
        xlon = (2.0*pi/360.0) * dlon_deg * (k-1)

        xp = R*sin(xcolat)*cos(xlon)
        yp = R*sin(xcolat)*sin(xlon)
        zp = R*cos(xcolat)

        ! Calculate rbf nodes


      end do

    end do
    l = l+1 !THIS IS ON THE BREAK CONDITION on R
  end do

end subroutine calculate_rbf_centers

subroutine calculate_ta16_field
end subroutine calculate_ta16_field

subroutine read_ta16_parameters(filename)

  character(*), intent(in) :: filename
  integer :: ta_file

  open(newunit=ta_file,file=filename,action="read")
  do i=1,size(A)
    read(ta_file,"(G15.6)") A(i)
  end do
  close(ta_file)

end subroutine read_ta16_parameters

end module TA16
