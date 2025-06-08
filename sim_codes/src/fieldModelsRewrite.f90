module compute

  implicit none
contains

subroutine rbf_model(pdyn,symh,xind,byimf,ps,&
                     xin,yin,zin,bx,by,bz,nx,ny,nz)

  implicit none
  integer :: nx,ny,nz
  real, intent(in) :: pdyn, symh, xind, byimf, ps
  real, intent(in) :: xin(nx),yin(ny),zin(nz)
  real, intent(out), dimension(nx,ny,nz) :: bx,by,bz

  real, dimension(1296) :: xx,yy,zz,st,rho,zsp,zcp,rhbr

  integer :: i,ix,iy,iz
  real :: xsm,ysm,zsm,fpd,x,y,z,cps,sps,tps,bxsm,bysm,bzsm
  real :: xp,yp,zp,xm,ym,zm
  real :: delta_zr,dtheta,sdt,cdtm1,dxp,dyp,dzp,dxm,dym,dzm
  real :: cp,cm,dcpx,dcmx,dcpy,dcmy,dcpz,dcmz
  real :: dcpx2,dcmx2,dcpy2,dcmy2,dcpz2,dcmz2,dcpxy,dcmxy,dcpxz,dcmxz,dcpyz,dcmyz
  real :: txcp,tycp,tzcp,txcm,tycm,tzcm,pxcp,pycp,pzcp,pxcm,pycm,pzcm
  real :: ctx,cty,ctz,stx,sty,stz,cpx,cpy,cpz,spx,spy,spz
  real :: act,ast,at,acp,asp,ap
  real, dimension(:), allocatable :: A

  integer, parameter :: Nlat = 8
  real, parameter :: Rlow = 3.3, Rhigh = 16.0, alpha = 3.0, rh = 8.0
  real, parameter :: pi = 3.14159265359, D = 4.0, pd_tr = 0.5

  allocate(A(23328))
  print *, "reading params"
  OPEN (UNIT=777,FILE='TA16_RBF.par')
  READ (777,100) A
  100 FORMAT (G15.6)
  CLOSE(777)
  print *, "finished, starting rbf_centers"
  call rbf_centers(xx,yy,zz,st,rho,zsp,zcp,rhbr)
  print *, "done rbf_centers, starting loop over x,y,z"

  do iz = 1,nz
  z = zin(iz)
  do iy = 1,ny
  y = yin(iy)
  do ix = 1,nx
  x = xin(ix)

!--------------------  START CALCULATING THE MODEL B-FIELD  ---------------------------------------
    xsm=x*cos(ps)-z*sin(ps)         !  rbf expansions are in sm coordinates
    ysm=y                           !  ->  convert x,y,z from gsw to sm 
    zsm=z*cos(ps)+x*sin(ps)

    fpd=sqrt(pdyn/2.0)-1.0

    cps=cos(ps)
    sps=sin(ps)
    tps=sps/cps

    bxsm=0.0
    bysm=0.0
    bzsm=0.0

    do i=1,1296

      xp = xx(i)
      yp = yy(i)
      zp = zz(i)
      xm = xp
      ym = yp
      zm =-zp

      delta_zr=rhbr(i)*tps
      dtheta =-asin(delta_zr)*st(i)
      sdt=sin(dtheta)
      cdtm1=cos(dtheta)-1.0
      dxp=xp*cdtm1+sdt*zcp(i)
      dyp=yp*cdtm1+sdt*zsp(i)
      dzp=zp*cdtm1-rho(i)*sdt
      dxm=xm*cdtm1-sdt*zcp(i)
      dym=ym*cdtm1-sdt*zsp(i)
      dzm=zm*cdtm1-rho(i)*sdt

      cp=sqrt((xsm-xp-dxp)**2+(ysm-yp-dyp)**2+(zsm-zp-dzp)**2+d**2)    ! rbf ch_i+
      cm=sqrt((xsm-xm-dxm)**2+(ysm-ym-dym)**2+(zsm-zm-dzm)**2+d**2)    ! rbf ch_i-
      dcpx=(xsm-xp-dxp)/cp
      dcmx=(xsm-xm-dxm)/cm
      dcpy=(ysm-yp-dyp)/cp
      dcmy=(ysm-ym-dym)/cm
      dcpz=(zsm-zp-dzp)/cp
      dcmz=(zsm-zm-dzm)/cm

      dcpx2=1.0/cp-dcpx**2/cp
      dcmx2=1.0/cm-dcmx**2/cm
      dcpy2=1.0/cp-dcpy**2/cp
      dcmy2=1.0/cm-dcmy**2/cm
      dcpz2=1.0/cp-dcpz**2/cp
      dcmz2=1.0/cm-dcmz**2/cm
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

      pxcp=2.0*dcpx-xsm*(dcpy2+dcpz2)+ysm*dcpxy+zsm*dcpxz    !  correction of 11/23/2021: z -> zsm
      pycp=2.0*dcpy-ysm*(dcpx2+dcpz2)+zsm*dcpyz+xsm*dcpxy
      pzcp=2.0*dcpz-zsm*(dcpx2+dcpy2)+xsm*dcpxz+ysm*dcpyz
      pxcm=2.0*dcmx-xsm*(dcmy2+dcmz2)+ysm*dcmxy+zsm*dcmxz
      pycm=2.0*dcmy-ysm*(dcmx2+dcmz2)+zsm*dcmyz+xsm*dcmxy
      pzcm=2.0*dcmz-zsm*(dcmx2+dcmy2)+xsm*dcmxz+ysm*dcmyz

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

      !-----------------   total field:    -----------------------------------

      act=a(i)+a(i+5184)*fpd+a(i+10368)*symh/50.0+a(i+15552)*xind
      ast=a(i+1296)+a(i+6480)*fpd+a(i+11664)*symh/50.0+a(i+16848)*xind
      at =a(i+20736)*byimf
      acp=a(i+2592)+a(i+7776)*fpd+a(i+12960)*symh/50.0+a(i+18144)*xind
      asp=a(i+3888)+a(i+9072)*fpd+a(i+14256)*symh/50.0+a(i+19440)*xind
      ap =a(i+22032)*byimf

      bxsm=bxsm+ctx*act+stx*ast+(txcp-txcm)*at+cpx*acp+spx*asp+(pxcp+pxcm)*ap
      bysm=bysm+cty*act+sty*ast+(tycp-tycm)*at+cpy*acp+spy*asp+(pycp+pycm)*ap
      bzsm=bzsm+ctz*act+stz*ast+(tzcp-tzcm)*at+cpz*acp+spz*asp+(pzcp+pzcm)*ap

    end do
    !----------------------------------------------------------------------------------------
    !   NOW CONVERT THE OBTAINED MAGNETIC FIELD VECTOR BACK FROM SM TO GSM (GSW) SYSTEM:
    !-----------------------------------------------------------------------------------------
    bx(ix,iy,iz)=bxsm*cps+bzsm*sps
    by(ix,iy,iz)=bysm
    bz(ix,iy,iz)=bzsm*cps-bxsm*sps

  end do
  end do
  end do
  !$OMP END PARALLEL DO 
end subroutine rbf_model

subroutine rbf_centers(xx,yy,zz,st,rho,zsp,zcp,rhbr)
  
  implicit none
  real, intent(out), dimension(:) :: xx,yy,zz,st,rho,zsp,zcp,rhbr

  integer, parameter :: Nlat = 8
  real, parameter :: Rlow = 3.3, Rhigh = 16.0, alpha = 3.0, rh = 8.0
  real, parameter :: pi = 3.14159265359, D = 4.0, pd_tr = 0.5
  real, dimension(22) :: A = [12.544,-0.194,0.305,0.0573,2.178,&
                              0.0571,-0.999,16.473,0.00152,0.382,&
                              0.0,0.0,0.0,0.0,-4.430,&
                              -0.636,-2.600,0.832,-5.328,1.103,&
                              -0.907,1.450]

  real :: dlat_deg,r,xcolatd, xlond, xcolat, xlon, xxxx, yyyy, zzzz
  real :: r1, rho1, ctt, stt, t, sp, cp, rm, psin, psis, f1, r0, rlast, dlon_deg
  integer :: nlon, i, k, l

  dlat_deg = 90.0/(Nlat-0.5)
  l=0
  R=Rlow
  print *, "A is"
  print *, A

  do while ( R.LE.Rhigh )

    do i=1,Nlat

      xcolatd = dlat_deg*(i-1)
      nlon = 4*(i-1)

      if (i.ne.1) then; dlon_deg = 360.0/nlon
      else; nlon=1; dlon_deg=0.0
      end if

      do k=1,nlon
        xlond = (k-1)*dlon_deg
        xcolat = xcolatd * 0.01745329252
        xlon = xlond * 0.01745329252
        xxxx = R*sin(xcolat)*cos(xlon)
        yyyy = R*sin(xcolat)*sin(xlon)
        zzzz = R*cos(xcolat)

        r1=sqrt(xxxx**2+(yyyy**2+zzzz**2))
        rho1=sqrt(yyyy**2+zzzz**2)

        ctt=xxxx/r1
        stt=sqrt(yyyy**2+zzzz**2)/r1
        t=atan2(stt,ctt)

        if (rho1.GT.1e-5) then
          sp = zzzz/rho1
          cp = yyyy/rho1
        else if (xxxx.GT.0.0) then
          sp = 0.0
          cp = 1.0
        else
          rm = 1000.0
          return
        end if

        psin=acos(ctt*cos(A(20))+stt*sin(A(20))*sp)
        psis=acos(ctt*cos(A(20))-stt*sin(A(20))*sp)

        f1=(sqrt(0.5*(1.0+ctt))+a(6)*2.0*stt*ctt*(1.0-exp(-t)))**(a(7)+a(11)*cp+a(12)*sp+a(14)*sp**2)
        r0=a(1)*(pd_tr)**a(2)*(1.0/2.0)
        rm=r0*f1+a(15)*(pd_tr)**a(16)*exp(a(17)*psin**a(22))+a(15)*(pd_tr)**a(16)*exp(a(17)*psis**a(22))    ! POSITION OF THE MODEL MAGNETOPAUSE

        if (r.gt.rm) then
          l=l+1                          !  counter of the rbf centers
          xx(l)=xxxx
          yy(l)=yyyy
          zz(l)=zzzz
          st(l)=sin(xcolat)
          rho(l)=r*st(l)
          zsp(l)=zz(l)*sin(xlon)
          zcp(l)=zz(l)*cos(xlon)
          rhbr(l)=rh/r*(1.0-(1.0+(r/rh)**alpha)**(1.0/alpha))
        end if
      end do
    end do
  rlast=r
  r=r*(nlat-0.5+pi/4.0)/(nlat-0.5-pi/4.0)  ! increment r by a fixed factor
  end do

end subroutine rbf_centers

end module compute
