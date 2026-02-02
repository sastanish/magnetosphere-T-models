module compute

  implicit none

contains

  subroutine TA16(date,data,x,y,z,Bx,By,Bz,nx,ny,nz)

    use control, only : calculate_TA16, setup_TA16_model

    !! Input !!
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), intent(in) :: data(8)
    integer, intent(in) :: date(4)

    !! Output !!
    real(8), intent(out) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)

    !! Vars for exploding arrays !!
    integer :: year,day,hour,mint
    real(8) :: ivx,ivy,ivz,tilt,pydn,symhc,nind,aby

    year = date(1)
    day = date(2)
    hour = date(3)
    mint = date(4)
    ivx = data(1)
    ivy = data(2)
    ivz = data(3)
    tilt = data(4)
    pydn = data(5)
    symhc = data(6)
    nind = data(7)
    aby = data(8)

    call setup_TA16_model()
    call calculate_TA16(year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symhc,nind,aby,x,y,z,Bx,By,Bz)
    
  end subroutine TA16

  subroutine TS05(date,data,x,y,z,Bx,By,Bz,nx,ny,nz)

    use control, only : calculate_TA16, setup_TA16_model

    !! Input !!
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), intent(in) :: data(14)
    integer, intent(in) :: date(4)

    !! Output !!
    real(8), intent(out) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)

    !! Vars for exploding arrays !!
    integer :: year,day,hour,mint
    real(8) :: ivx,ivy,ivz,tilt,pydn,symh,iby,ibz,w1,w2,w3,w4,w5,w6

    year = date(1)
    day = date(2)
    hour = date(3)
    mint = date(4)
    ivx = data(1)
    ivy = data(2)
    ivz = data(3)
    tilt = data(4)
    pydn = data(5)
    symh = data(6)
    iby = data(7)
    ibz = data(8)
    w1 = data(9)
    w2 = data(10)
    w3 = data(11)
    w4 = data(12)
    w5 = data(13)
    w6 = data(14)

    call calculate_TS05(year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symh,iby,ibz,w1,w2,w3,w4,w5,w6,x,y,z,Bx,By,Bz)
    
  end subroutine TS05

! subroutine reconnection(x,y,z,bx,by,bz,nx,ny,nz,recon)
!
!   integer :: nx, ny, nz
!
!   real(8), intent(in) :: x(nx), y(ny), z(nz)
!   real(8), dimension(nx,ny,nz), intent(in) :: bx, by, bz
!   real(8), intent(out) :: Mout(23,nx,ny,nz)
!
!   call compute_all_metrics(x,y,z,bx,by,bz,Mout)
!
! end subroutine reconnection

end module compute
