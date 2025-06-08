program main

  use TS04c, only : T04_s
  use geopack, only : IGRF_GSW_08, RECALC_08
  use reconnection_metrics, only : init_no_file, write_netcdf_vector, compute_c2_t1_coeff

  implicit none

  !Input parameters
  real(8), parameter :: Lxm = -20
  real(8), parameter :: Lxp = 2
  real(8), parameter :: Lym = -4
  real(8), parameter :: Lyp = 4
  real(8), parameter :: Lzm = -4
  real(8), parameter :: Lzp = 4

  !Target File
  character(*), parameter :: omni_file="selected_times.txt"

  !Output Grid
  integer, parameter :: nx=400, ny=100, nz=100
  real(8), dimension(nz,ny,nx) :: bx,by,bz
  real(8), dimension(nx) :: x
  real(8), dimension(ny) :: y
  real(8), dimension(nz) :: z

  !Tsyganenko parameters
  integer :: iyear,iday,ihour,imin,isec
  character(4) :: syear
  character(4) :: sday
  character(3) :: shour, smin
  real(8) :: bx_in, by_in,bz_in,vx_in,vy_in,vz_in,density,temp,symh
  integer :: imfflag,iswflag
  real(8) :: tilt,pdyn,w1,w2,w3,w4,w5,w6
  real(8) :: parmod(10)

  !Dummy variables
  real(8) :: dx,dy,dz, ddx, ddy, ddz
  real(8) :: xx,yy,zz
  integer :: ix,iy,iz

  character(50) :: outfile

  ! Setup dimension grids
  do ix=1,nx
    x(ix) = Lxm + ((Lxp-Lxm)/nx)*ix
  end do
  do iy=1,ny
    y(iy) = Lym + ((Lyp-Lym)/ny)*iy
  end do
  do iz=1,nz
    z(iz) = Lzm + ((Lzp-Lzm)/nz)*iz
  end do

  bx=0.0
  by=0.0
  bz=0.0

  call init_no_file(x,y,z,bx,by,bz)

  ! Note the bad format method. Using as the is the intended input for the
  ! old f77 code we are using in TS04c
  505 FORMAT (2I4,2I3,3F8.2,3F8.1,F7.2,F9.0,F7.1,2(3X,I2),F8.4,7F7.2)
  open(unit=1,file=omni_file,status="OLD",action="READ")

  do
    read(1, 505, end=200) iyear,iday,ihour,imin,bx_in, &
                          by_in,bz_in,vx_in,vy_in,vz_in,density, &
                          temp,symh,imfflag,iswflag, &
                          tilt,pdyn,w1,w2,w3,w4,w5,w6

    parmod = [pdyn, symh, by_in, bz_in, w1, w2, w3, w4, w5, w6]

    isec = 0
    !adjust v
    vy_in = vy_in + 29.78

    call RECALC_08(iyear, iday, ihour, imin, isec, vx_in, vy_in, vz_in)

    do ix = 1,nx
    xx = x(i)
    do iy = 1,ny
    yy = y(i)
    do iz = 1,nz
    zz = z(i)
      call IGRF_GSW_08(xx,yy,zz,dx,dy,dz)
      call T04_s(0,parmod,tilt,xx,yy,zz,ddx,ddy,ddz)
      bx(iz,iy,ix) = dx + ddx
      by(iz,iy,ix) = dy + ddy
      bz(iz,iy,ix) = dz + ddz
    end do
    end do
    end do

    write(syear,"(I4)") iyear
    write(sday,"(I3)") iday
    write(shour,"(I2)") ihour
    write(smin,"(I2)") imin

    outfile = trim(syear)//"_"//adjustl(sday)//"_"//adjustl(shour)//"_"//adjustl(smin)//".nc"
    call write_netcdf_vector(trim(outfile),"b",bx,by,bz)

  end do

  200 continue

end program main
