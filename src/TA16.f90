module TA16

implicit none

integer, parameter :: dp=selected_real_kind(1.d0)

real(dp) :: A(23328)

contains

subroutine calculate_rbf_centers

  real(dp), intent(out), dimension(1296) :: x,y,z,st,rho,zsp,zcp,rhbr

  integer :: klat, nlat

  real(dp) :: rlow_grid, rhigh_grid
  real(dp) :: dlat_deg, dlon_deg, xcolatd
  real(dp) :: Radius
  real(dp) :: theta_n, theta_s, ctn, cts, stn, sts, sp, cp

  real(dp), parameter :: rh=8.0,alpha=3.0,pd_tr=0.5
  real(dp), parameter :: psi=0.0,pm=0.0,Bzimf_tr=0.0
  real(dp), parameter :: pi=4.d0*atan(1.d0)

  real(dp), parameter, dimension(22) :: An(0:21) = (12.544D0,-0.194D0,0.305D0,0.0573D0,2.178D0,0.0571D0,-0.999D0,16.473D0,0.00152D0,0.382D0,0.0431D0,-0.00763D0,-0.210D0,0.0405D0,-4.430D0,-0.636D0,-2.600D0,0.832D0,-5.328D0,1.103D0,-0.907D0,1.450D0)

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
        else if ( xp > 0.d0) then
          sp = 0.d0
          cp = 1.d0
        else
          Radius = 1000.d0
          return
        end if

        psin = acos(ctt*ctn+stt*stn*sp)
        psis = acos(ctt*cts-stt*sts*sp)

        dn = An(16) + (An(17)+An(18)*psi)*psi
        ds = An(16) - (An(17)-An(18)*psi)*psi

        cn = An(14)*(pd_tr+pm)**An(15)

        b0 = An(6) + An(7)*(exp(An(8)*Bzimf_tr)-1.d0)/(exp(An(9)*Bzimf_tr)+1.d0)
        b2 = An(11)+An(12)*psi

        f1 = ( sqrt(0.5d0 * (1.d0+ctt)) + An(5)*2.d0*stt*ctt*(1.d0-exp(-T)))**(b0+An(10)*cp+b2*sp+An(13)*sp**2)

        R0 = An(0)*(pd_tr+pm)**An(1) * ( 1.d0 + An(2)*(exp(An(3)*Bzimf_tr)-1.d0)/(exp(An(4)*Bzimf_tr)+1.d0) )
        Rm = R0*f1 + cn*exp(dn*psin**en) + cs*exp(ds*psis**es)

        if (Radius > Rm) exit

        l = l+1

        x(l) = xp
        y(l) = yp
        z(l) = zp
        st(l) = sin(xcolat)
        rho(l) = Radius*st(l)
        zsp(l) = zp(l)*sin(xlon)
        zcp(l) = zp(l)*cos(xlon)
        rhbr(l) = rh/Radius * (1.d0 - (1.d0 + (Radius/rh)**alpha)**(1.d0/alpha))

      end do

      last_radius = Radius
      Radius = Radius*(nlat - 0.5d0 + pi/4.d0)/(nlat - 0.5d0 - pi/4.d0)

      if (Radius > rhigh_grid) return

    end do
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
