  module TA16

  use iso_fortran_env, only: real64

  implicit none

  integer, parameter :: dp=real64

  real(dp) :: A(23328) ! These are the TA16 fitting parameters
  real(dp), dimension(1296) :: xrb,yrb,zrb,strb,rhorb,zsprb,zcprb,rhbrrb !rbf coordinates

  contains

  subroutine calculate_rbf_centers

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

        if (Radius <= Rm) then
          l = l+1
          xrb(l) = xp
          yrb(l) = yp
          zrb(l) = zp
          strb(l) = sin(xcolat)
          rhorb(l) = Radius*strb(l)
          zsprb(l) = zrb(l)*sin(xlon)
          zcprb(l) = zrb(l)*cos(xlon)
          rhbrrb(l) = rh/Radius * (1.0_dp - (1.0_dp + (Radius/rh)**alpha)**(1.0_dp/alpha))
        end if

      end do lon_loop
    end do lat_loop

    Radius = Radius*(nlat - 0.5_dp + pi/4.0_dp)/(nlat - 0.5_dp - pi/4.0_dp)

  end do

end subroutine calculate_rbf_centers

subroutine calculate_ta16_field(pydn, symHc, Nind, Byimf, tilt, x, y, z, bx, by, bz)
  
  real(dp), intent(out) :: bx, by, bz
  real(dp), intent(in) :: pydn, symHc, Nind, Byimf, tilt, x, y, z
  real(dp) :: symH

  real(dp) :: xsm, ysm, zsm, fpd, cps, sps, tps
  real(dp) :: acp, act, ap, asp, ast, at, cdtm1, cm, cp, cpx, cpy, cpz, ctx, cty, ctz
  real(dp) :: dcmx, dcmx2, dcmxy, dcmxz, dcmy, dcmy2, dcmyz, dcmz, dcmz2, dcpx, dcpx2, dcpxy, dcpxz, dcpy, dcpy2, dcpyz, dcpz, dcpz2
  real(dp) :: delta_zr, dtheta, dxm, dxp, dym, dyp, dzm, dzp
  real(dp) :: pxcm, pxcp, pycm, pycp, pzcm, pzcp, sdt, spx, spy, spz, stx, sty, stz
  real(dp) :: txcm, txcp, tycm, tycp, tzcm, tzcp, xm, xp, ym, yp, zm, zp

  integer :: i
  real(dp), parameter :: d=4.0_dp

  xsm = x*cos(tilt) - z*sin(tilt)
  ysm = y
  zsm = z*cos(tilt) + x*sin(tilt)

  fpd = sqrt(pydn/2.0_dp) - 1.0_dp
  symH = symHc/50.0_dp

  cps = cos(tilt)
  sps = sin(tilt)
  tps = sps/cps

  bx = 0.0_dp
  by = 0.0_dp
  bz = 0.0_dp

  do i = 1,1296

    xp = xrb(i)
    yp = yrb(i)
    zp = zrb(i)
    xm = xp
    ym = yp
    zm =-zp

    delta_zr=rhbrrb(i)*tps
    dtheta =-asin(delta_zr)*strb(i)
    sdt=sin(dtheta)
    cdtm1=cos(dtheta)-1.0_dp
    dxp=xp*cdtm1+sdt*zcprb(i)
    dyp=yp*cdtm1+sdt*zsprb(i)
    dzp=zp*cdtm1-rhorb(i)*sdt
    dxm=xm*cdtm1-sdt*zcprb(i)
    dym=ym*cdtm1-sdt*zsprb(i)
    dzm=zm*cdtm1-rhorb(i)*sdt

    cp=sqrt((xsm-xp-dxp)**2+(ysm-yp-dyp)**2+(zsm-zp-dzp)**2+d**2)    ! rbf ch_i+
    cm=sqrt((xsm-xm-dxm)**2+(ysm-ym-dym)**2+(zsm-zm-dzm)**2+d**2)    ! rbf ch_i-
    dcpx=(xsm-xp-dxp)/cp
    dcmx=(xsm-xm-dxm)/cm
    dcpy=(ysm-yp-dyp)/cp
    dcmy=(ysm-ym-dym)/cm
    dcpz=(zsm-zp-dzp)/cp
    dcmz=(zsm-zm-dzm)/cm

    dcpx2=1.0_dp/cp-dcpx**2/cp
    dcmx2=1.0_dp/cm-dcmx**2/cm
    dcpy2=1.0_dp/cp-dcpy**2/cp
    dcmy2=1.0_dp/cm-dcmy**2/cm
    dcpz2=1.0_dp/cp-dcpz**2/cp
    dcmz2=1.0_dp/cm-dcmz**2/cm
    dcpxy=-dcpx*dcpy/cp
    dcmxy=-dcmx*dcmy/cm
    dcpxz=-dcpx*dcpz/cp
    dcmxz=-dcmx*dcmz/cm
    dcpyz=-dcpy*dcpz/cp
    dcmyz=-dcmy*dcmz/cm

    txcp=zsm*dcpy-ysm*dcpz
    tycp=xsm*dcpz-zsm*dcpx
    tzcp=ysm*dcpx-xsm*dcpy
    txcm=zsm*dcmy-ysm*dcmz
    tycm=xsm*dcmz-zsm*dcmx
    tzcm=ysm*dcmx-xsm*dcmy
 
    pxcp=2.0_dp*dcpx-xsm*(dcpy2+dcpz2)+ysm*dcpxy+zsm*dcpxz
    pycp=2.0_dp*dcpy-ysm*(dcpx2+dcpz2)+zsm*dcpyz+xsm*dcpxy
    pzcp=2.0_dp*dcpz-zsm*(dcpx2+dcpy2)+xsm*dcpxz+ysm*dcpyz
    pxcm=2.0_dp*dcmx-xsm*(dcmy2+dcmz2)+ysm*dcmxy+zsm*dcmxz
    pycm=2.0_dp*dcmy-ysm*(dcmx2+dcmz2)+zsm*dcmyz+xsm*dcmxy
    pzcm=2.0_dp*dcmz-zsm*(dcmx2+dcmy2)+xsm*dcmxz+ysm*dcmyz

    ctx = cps*(txcp+txcm)
    cty = cps*(tycp+tycm)
    ctz = cps*(tzcp+tzcm)

    stx = sps*(txcp-txcm)
    sty = sps*(tycp-tycm)
    stz = sps*(tzcp-tzcm)

    cpx = cps*(pxcp-pxcm)
    cpy = cps*(pycp-pycm)
    cpz = cps*(pzcp-pzcm)

    spx = sps*(pxcp+pxcm)
    spy = sps*(pycp+pycm)
    spz = sps*(pzcp+pzcm)

    ! Total field calculation
    act=A(i)+A(i+5184)*fpd+A(i+10368)*symH+A(i+15552)*Nind
    ast=A(i+1296)+A(i+6480)*fpd+A(i+11664)*symH+A(i+16848)*Nind
    at =A(i+20736)*Byimf
    acp=A(i+2592)+A(i+7776)*fpd+A(i+12960)*symH+A(i+18144)*Nind
    asp=A(i+3888)+A(i+9072)*fpd+A(i+14256)*symH+A(i+19440)*Nind
    ap =A(i+22032)*Byimf

    bx=bx+ctx*act+stx*ast+(txcp-txcm)*at+cpx*acp+spx*asp+(pxcp+pxcm)*ap
    by=by+cty*act+sty*ast+(tycp-tycm)*at+cpy*acp+spy*asp+(pycp+pycm)*ap
    bz=bz+ctz*act+stz*ast+(tzcp-tzcm)*at+cpz*acp+spz*asp+(pzcp+pzcm)*ap

  end do

  ! Convert back from SM to GSM
  bx=bx*cps+bz*sps
  by=by
  bz=bz*cps-bx*sps

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
