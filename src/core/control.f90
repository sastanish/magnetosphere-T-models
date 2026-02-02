module control
  implicit none

contains

  subroutine calculate_TA16(year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symhc,nind,aby,x,y,z,Bx,By,Bz)

    use TA16, only : RBF_MODEL_2016
    use geopack, only : RECALC_08, IGRF_GSW_08
    !$ use omp_lib

    implicit none

    !! Inputs !!
    ! Grid
    real(8), intent(in), dimension(:) :: x,y,z
    ! Input File params
    integer, intent(in) :: year,day,hour,mint
    real(8), intent(in) :: ivx,ivy,ivz,tilt,pydn,symhc,nind,aby

    !! Outputs !!
    real(8), intent(out), dimension(:,:,:) :: Bx, By, Bz

    !! Vars !!
    integer :: n, i, j, k
    real(8) :: xx,yy,zz,bbx,bby,bbz,hhx,hhy,hhz
    real(8) :: parmod(10)

    ivy = ivy + 29.78 !velocity correction

    call RECALC_08(year, day, hour, mint, 0, ivx, ivy, ivz)

    parmod(1) = pydn
    parmod(2) = symhc
    parmod(3) = nind
    parmod(4) = aby

    !$OMP PARALLEL PRIVATE(parmod,xx,yy,zz,hhx,hhy,hhz,bbx,bby,bbz,i,j,k) SHARED(x,y,z,Bx,By,Bz)
    !$OMP DO COLLAPSE(3)
    do k = 1,size(z)
      do j = 1,size(y)
        do i = 1,size(x)

          xx = x(i)
          yy = y(j)
          zz = z(k)
          hhx = 0
          hhy = 0
          hhz = 0
          bbx = 0
          bby = 0
          bbz = 0

          ! External field
          call RBF_MODEL_2016(0,parmod,tilt,xx,yy,zz,hhx,hhy,hhz)

          ! Internal field
          call IGRF_GSW_08(xx,yy,zz,bbx,bby,bbz)

          Bx(i,j,k) = hhx + bbx
          By(i,j,k) = hhy + bby
          Bz(i,j,k) = hhz + bbz

        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine calculate_TA16

  subroutine calculate_TS05(year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symh,iby,ibz,w1,w2,w3,w4,w5,w6,x,y,z,Bx,By,Bz)

  use TS05, only : T04_s
  use geopack, only : RECALC_08, IGRF_GSW_08
  !$ use omp_lib

  implicit none

  !! Inputs !!
  ! Grid
  real(8), intent(in), dimension(:) :: x,y,z
  ! Input File params
  integer, intent(in) :: year,day,hour,mint
  real(8), intent(in) :: ivx,ivy,ivz,tilt,pydn,symh,iby,ibz,w1,w2,w3,w4,w5,w6

  !! Outputs !!
  real(8), intent(out), dimension(:,:,:) :: Bx, By, Bz

  !! Vars !!
  integer :: n, i, j, k
  real(8) :: xx,yy,zz,bbx,bby,bbz,hhx,hhy,hhz
  real(8) :: parmod(10)

  ivy = ivy + 29.78 !velocity correction

  call RECALC_08(year, day, hour, mint, 0, ivx, ivy, ivz)

  !$OMP PARALLEL PRIVATE(parmod,xx,yy,zz,hhx,hhy,hhz,bbx,bby,bbz,i,j,k) SHARED(x,y,z,Bx,By,Bz)
  parmod(1) = pydn
  parmod(2) = symh
  parmod(3) = iby
  parmod(4) = ibz
  parmod(5) = w1
  parmod(6) = w2
  parmod(7) = w3
  parmod(8) = w4
  parmod(9) = w5
  parmod(10) = w6

  !$OMP DO COLLAPSE(3)
  do k = 1,size(z)
    do j = 1,size(y)
      do i = 1,size(x)

        xx = x(i)
        yy = y(j)
        zz = z(k)
        hhx = 0
        hhy = 0
        hhz = 0
        bbx = 0
        bby = 0
        bbz = 0

        ! External field
        call T04_s(0,parmod,tilt,xx,yy,zz,hhx,hhy,hhz)

        ! Internal field
        call IGRF_GSW_08(xx,yy,zz,bbx,bby,bbz)

        Bx(i,j,k) = hhx + bbx
        By(i,j,k) = hhy + bby
        Bz(i,j,k) = hhz + bbz

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  end subroutine calculate_TS05

  subroutine setup_TA16_model()

    use TA16, only : CALCULATE_RBF_CENTERS,READ_TA16_PARS

    call read_TA16_pars
    call calculate_rbf_centers

  end subroutine setup_TA16_model

  subroutine setup_grid(InFileName,x,y,z,Bx,By,Bz)

    implicit none

    integer :: nx, ny, nz, i
    integer :: file

    real(8) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(8), dimension(:), allocatable, intent(inout) :: x,y,z
    real(8), dimension(:,:,:), allocatable, intent(inout) :: Bx, By, Bz

    character(*), intent(in) :: InFileName

    ! Get input parameters for grid
    open(newunit=file, file=InFileName, status='old', action='read')
    read(file, *) !skip the header line
    read(file, *) xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz
    close(file)

    allocate(x(nx))
    allocate(y(ny))
    allocate(z(nz))
    allocate(Bx(nx,ny,nz))
    allocate(By(nx,ny,nz))
    allocate(Bz(nx,ny,nz))

    ! Setup grid
    do i = 1,nx
      x(i) = xmin + (i-1)*(xmax-xmin)/nx
    end do
    do i = 1,ny
      y(i) = ymin + (i-1)*(ymax-ymin)/ny
    end do
    do i = 1,nz
      z(i) = zmin + (i-1)*(zmax-zmin)/nz
    end do

  end subroutine setup_grid

  subroutine read_TA16_input_data(filename,year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symhc,nind,aby)

    implicit none

    integer :: Nlines, file, i
    character(*), intent(in) :: filename

    integer, dimension(:), allocatable, intent(inout) :: year,day,hour,mint
    real(8), dimension(:), allocatable, intent(inout) :: ivx,ivy,ivz,tilt,pydn,symhc,nind,aby

    ! un-used input params
    integer :: aeind, symh, imf
    real(8) :: ibx, iby, ibz, den, temp, sw, rp, abx, abz

    call get_num_lines(filename,nlines)

    allocate(year(nlines))
    allocate(day(nlines))
    allocate(hour(nlines))
    allocate(mint(nlines))
    allocate(ivx(nlines))
    allocate(ivy(nlines))
    allocate(ivz(nlines))
    allocate(tilt(nlines))
    allocate(pydn(nlines))
    allocate(symhc(nlines))
    allocate(nind(nlines))
    allocate(aby(nlines))

    open(newunit=file, file=filename, status='old', action='read')
    do i = 1,nlines
      read(file,*) year(i),day(i),hour(i),mint(i),ibx,iby,ibz,ivx(i),ivy(i),ivz(i),den,temp,pydn(i),aeind,symh,imf,sw,tilt(i),rp,abx,aby(i),abz,nind(i),symhc(i)
    end do
    close(file)
  end subroutine read_TA16_input_data

  subroutine read_TS05_input_data(filename,year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symh,iby,ibz,w1,w2,w3,w4,w5,w6)

    implicit none

    integer :: Nlines, file, i
    character(*), intent(in) :: filename

    integer, dimension(:), allocatable, intent(inout) :: year,day,hour,mint
    real(8), dimension(:), allocatable, intent(inout) :: ivx,ivy,ivz,tilt,pydn,w1,w2,w3,w4,w5,w6,iby,ibz,symh

    ! un-used input params
    integer :: imf, sw
    real(8) :: ibx, den, temp

    call get_num_lines(filename,nlines)

    allocate(year(nlines))
    allocate(day(nlines))
    allocate(hour(nlines))
    allocate(mint(nlines))
    allocate(ivx(nlines))
    allocate(ivy(nlines))
    allocate(ivz(nlines))
    allocate(pydn(nlines))
    allocate(symh(nlines))
    allocate(tilt(nlines))
    allocate(iby(nlines))
    allocate(ibz(nlines))
    allocate(w1(nlines))
    allocate(w2(nlines))
    allocate(w3(nlines))
    allocate(w4(nlines))
    allocate(w5(nlines))
    allocate(w6(nlines))

    open(newunit=file, file=filename, status='old', action='read')
    do i = 1,nlines
      read(file,*) year(i),day(i),hour(i),mint(i),ibx,iby(i),ibz(i),ivx(i),ivy(i),ivz(i),den,temp,symh(i),imf,sw,tilt(i),pydn(i),w1(i),w2(i),w3(i),w4(i),w5(i),w6(i)
    end do
    close(file)
  end subroutine read_TS05_input_data

  subroutine get_num_lines(filename,nlines)
    implicit none

    character(*), intent(in) :: filename
    integer, intent(out) :: nlines
    integer :: file, i, stat

    nlines = 0
    open(newunit=file, file=filename, status='old', action='read')
    do 
      read(file, *, iostat=stat)
      if (stat < 0) exit !checking for end of file
      nlines = nlines + 1
    end do
    close(file)
  end subroutine get_num_lines

end module control
