  module TA16

  use iso_fortran_env, only: real64

  implicit none

  integer, parameter :: dp=real64

  !real(dp) :: A(23328)

  contains

  subroutine calculate_rbf_centers(x,y,z,st,rho,zsp,zcp,rhbr)

    real(dp), intent(out), dimension(1296) :: x,y,z,st,rho,zsp,zcp,rhbr

    integer :: klat, nlat, nlon

    real(dp) :: rlow_grid, rhigh_grid
    real(dp) :: dlat_deg, dlon_deg, xcolat
    real(dp) :: Radius
    real(dp) :: theta_n, theta_s, ctn, cts, stn, sts, sp, cp
    real(dp) :: zp, yp, xp, xlond, xlon, T, stt, rm, R1, R0
    real(dp) :: psin, psis, last_radius, f1, ds, dn, ctt, cn, b2, b0

    real(dp), parameter :: rh=8.0_dp,alpha=3.0_dp,pd_tr=0.5_dp
    real(dp), parameter :: psi=0.0_dp,pm=0.0_dp,Bzimf_tr=0.0_dp
  real(dp), parameter :: pi=4.0_dp*atan(1.0_dp)

  real(dp), dimension(22) :: An(0:21) = [12.544_dp,-0.194_dp,0.305_dp,0.0573_dp,2.178_dp,0.0571_dp,-0.999_dp,16.473_dp,0.00152_dp,0.382_dp,0.0431_dp,-0.00763_dp,-0.210_dp,0.0405_dp,-4.430_dp,-0.636_dp,-2.600_dp,0.832_dp,-5.328_dp,1.103_dp,-0.907_dp,1.45_dp]

  integer :: i,j,k,l

  klat = 8
  nlat = klat + 1

  rlow_grid = 3.3_dp
  rhigh_grid = 16.0_dp

  dlat_deg = 90.0_dp/(nlat-0.5_dp)

  l = 0

  Radius = rlow_grid

  An(10) = 0.0_dp ! Simplified version of Lin et al model
  An(11) = 0.0_dp
  An(13) = 0.0_dp
  An(14) = 0.0_dp
  
  do while (Radius <= rhigh_grid)

lat_loop : do i = 1,nlat

      nlon = 4*(i-1)

      if (i /= 1) then
        dlon_deg = 360.0_dp/nlon
      else
        nlon = 1
        dlon_deg = 0.0_dp
      end if

lon_loop : do k = 1,nlon

        ! Coordinate transform
        xlond = (k-1)*dlon_deg

        xcolat = (2.0_dp*pi/360.0_dp) * dlat_deg * (i-1.0_dp)
        xlon = (2.0_dp*pi/360.0_dp) * dlon_deg * (k-1.0_dp)

        xp = Radius*sin(xcolat)*cos(xlon)
        yp = Radius*sin(xcolat)*sin(xlon)
        zp = Radius*cos(xcolat)

        ! Calculate rbf nodes
        theta_n = An(19)+An(20)*psi
        theta_s = An(19)-An(20)*psi

        ctn = cos(theta_n)
        cts = cos(theta_s)
        stn = sin(theta_n)
        sts = sin(theta_s)

        R1 = sqrt(xp**2 + yp**2 + zp**2)

        ctt = xp/R1
        stt = sqrt(yp**2 + zp**2)/R1
        T = atan2(stt,ctt)

        if (sqrt(yp**2 + zp**2) > 1.d-5) then
          sp = zp/R1
          cp = yp/R1
        else if ( xp > 0.0_dp) then
          sp = 0.0_dp
          cp = 1.0_dp
        else
          Radius = 1000.0_dp
          return
        end if

        psin = acos(ctt*ctn+stt*stn*sp)
        psis = acos(ctt*cts-stt*sts*sp)

        dn = An(16) + (An(17)+An(18)*psi)*psi
        ds = An(16) - (An(17)-An(18)*psi)*psi

        cn = An(14)*(pd_tr+pm)**An(15)

        b0 = An(6) + An(7)*(exp(An(8)*Bzimf_tr)-1.0)/(exp(An(9)*Bzimf_tr)+1.0)
        b2 = An(11)+An(12)*psi

        f1 = ( sqrt(0.5_dp * (1.0_dp+ctt)) + An(5)*2.0_dp*stt*ctt*(1.0_dp-exp(-T)))**(b0+An(10)*cp+b2*sp+An(13)*sp**2)

        R0 = An(0)*(pd_tr+pm)**An(1) * ( 1.0_dp + An(2)*(exp(An(3)*Bzimf_tr)-1.0_dp)/(exp(An(4)*Bzimf_tr)+1.0_dp) )
        Rm = R0*f1 + cn*exp(dn*psin**An(21)) + cn*exp(ds*psis**An(21))

        if (Radius .GT. Rm) cycle

        l = l+1
        x(l) = xp
        y(l) = yp
        z(l) = zp
        st(l) = sin(xcolat)
        rho(l) = Radius*st(l)
        zsp(l) = z(l)*sin(xlon)
        zcp(l) = z(l)*cos(xlon)
        rhbr(l) = rh/Radius * (1.0_dp - (1.0_dp + (Radius/rh)**alpha)**(1.0_dp/alpha))

      end do lon_loop
    end do lat_loop

    Radius = Radius*(nlat - 0.5_dp + pi/4.0_dp)/(nlat - 0.5_dp - pi/4.0_dp)

  end do

end subroutine calculate_rbf_centers

subroutine calculate_ta16_field
end subroutine calculate_ta16_field

!subroutine read_ta16_parameters(filename)

! character(*), intent(in) :: filename
! integer :: ta_file, i
!
! open(newunit=ta_file,file=filename,action="read")
! do i=1,size(A)
!   read(ta_file,"(G15.6)") A(i)
! end do
! close(ta_file)

!end subroutine read_ta16_parameters

end module TA16
