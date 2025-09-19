module TA16

use iso_fortran_env, only: real64

implicit none

integer, parameter :: dp=real64

real(dp) :: A(23328)

contains

subroutine calculate_rbf_centers(x,y,z,st,rho,zsp,zcp,rhbr)

  real(dp), intent(out), dimension(1296) :: x,y,z,st,rho,zsp,zcp,rhbr

  integer :: klat, nlat, nlon

  real(dp) :: rlow_grid, rhigh_grid
  real(dp) :: dlat_deg, dlon_deg, xcolat
  real(dp) :: Radius
  real(dp) :: theta_n, theta_s, ctn, cts, stn, sts, sp, cp
  real(dp) :: zp, yp, xp, xlond, xlon, T, stt, rm, R1, R0
  real(dp) :: psin, psis, last_radius, f1, es, en, ds, dn, ctt, cn, b2, b0

  real(dp), parameter :: rh=8.0,alpha=3.0,pd_tr=0.5
  real(dp), parameter :: psi=0.0,pm=0.0,Bzimf_tr=0.0
  real(dp), parameter :: pi=4.0_dp*atan(1.0_dp)

  real(dp), dimension(22) :: An(0:21) = [12.544,-0.194,0.305,0.0573,2.178,0.0571,-0.999,16.473,0.00152,0.382,0.0431,-0.00763,-0.210,0.0405,-4.430,-0.636,-2.600,0.832,-5.328,1.103,-0.907,1.45]

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
        else if ( xp > 0.0) then
          sp = 0.0
          cp = 1.0
        else
          Radius = 1000.0
          return
        end if

        psin = acos(ctt*ctn+stt*stn*sp)
        psis = acos(ctt*cts-stt*sts*sp)

        dn = An(16) + (An(17)+An(18)*psi)*psi
        ds = An(16) - (An(17)-An(18)*psi)*psi

        cn = An(14)*(pd_tr+pm)**An(15)

        b0 = An(6) + An(7)*(exp(An(8)*Bzimf_tr)-1.0)/(exp(An(9)*Bzimf_tr)+1.0)
        b2 = An(11)+An(12)*psi

        f1 = ( sqrt(0.50 * (1.0+ctt)) + An(5)*2.0*stt*ctt*(1.0-exp(-T)))**(b0+An(10)*cp+b2*sp+An(13)*sp**2)

        R0 = An(0)*(pd_tr+pm)**An(1) * ( 1.0 + An(2)*(exp(An(3)*Bzimf_tr)-1.0)/(exp(An(4)*Bzimf_tr)+1.0) )
        Rm = R0*f1 + cn*exp(dn*psin**en) + cn*exp(ds*psis**es)

        if (Radius > Rm) exit

        l = l+1

        x(l) = xp
        y(l) = yp
        z(l) = zp
        st(l) = sin(xcolat)
        rho(l) = Radius*st(l)
        zsp(l) = z(l)*sin(xlon)
        zcp(l) = z(l)*cos(xlon)
        rhbr(l) = rh/Radius * (1.0 - (1.0 + (Radius/rh)**alpha)**(1.0/alpha))

      end do

      last_radius = Radius
      Radius = Radius*(nlat - 0.5 + pi/4.0)/(nlat - 0.5 - pi/4.0)

      if (Radius > rhigh_grid) return

    end do
  end do

end subroutine calculate_rbf_centers

subroutine calculate_ta16_field
end subroutine calculate_ta16_field

subroutine read_ta16_parameters(filename)

  character(*), intent(in) :: filename
  integer :: ta_file, i

  open(newunit=ta_file,file=filename,action="read")
  do i=1,size(A)
    read(ta_file,"(G15.6)") A(i)
  end do
  close(ta_file)

end subroutine read_ta16_parameters

end module TA16
