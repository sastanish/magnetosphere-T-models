program main

  use TA16, only : RBF_MODEL_2016,CALCULATE_RBF_CENTERS,READ_TA16_PARS
  use geopack, only : RECALC_08, IGRF_GSW_08

  !$ use omp_lib

  implicit none

  real(8), dimension(:), allocatable :: x_points,y,z

  ! Input File params
  integer, dimension(:), allocatable :: year,day,hour,mint
  real(8), dimension(:), allocatable :: ivx,ivy,ivz,tilt,pydn,symhc,nind,aby
  real(8) :: pressure, flux

  ! vars
  integer :: n, i, j, k, f_file, p_file

  character(4) :: str_ind

  real(8) :: xx,yy,zz,bbx,bby,bbz,hhx,hhy,hhz !dummy variables
  real(8) :: parmod(10)

  !command line args
  integer :: start_ind, end_ind
  character(4) :: start_str, end_str

  call GET_COMMAND_ARGUMENT(1,start_str)
  call GET_COMMAND_ARGUMENT(2,end_str)
  read(start_str,'(I5)') start_ind
  read(end_str,'(I5)') end_ind

  call setup_grid(x_points,y,z)
  call read_input_data('input_data.lst',year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symhc,nind,aby)
  call read_TA16_pars
  call calculate_rbf_centers

  ! Setup output files
  open(newunit=p_file, file="pressure.lst", status="new")
  open(newunit=f_file, file="flux.lst", status="new")
  write(p_file,"(A)") "#Avg Pressure across y-z plane in magnetotail"
  write(f_file,"(A)") "#Avg |Bx| across y-z plane in magnetotail"
  write(p_file,"(A)") "#FORMAT:"
  write(f_file,"(A)") "#FORMAT:"
  write(p_file,"(A, F6.2, 4x)",advance="no") "YEAR    DAY    HOUR    MIN    x:", x_points(1)
  write(f_file,"(A, F6.2, 4x)",advance="no") "YEAR    DAY    HOUR    MIN    x:", x_points(1)
  do i = 2,size(x_points)
    write(p_file,"(F6.2, 6x)",advance="no") x_points(i)
    write(f_file,"(F6.2, 6x)",advance="no") x_points(i)
  end do
  write(p_file,*)
  write(f_file,*)

  ivy = ivy + 29.78 !velocity correction

  do n = start_ind,end_ind

    write(f_file,"(I4, 6x, I3, 6x, I4, 6x, I3, 6x)",advance="no") year(n), day(n), hour(n), mint(n)
    write(p_file,"(I4, 6x, I3, 6x, I4, 6x, I3, 6x)",advance="no") year(n), day(n), hour(n), mint(n)
    print *, n

    call RECALC_08(year(n), day(n), hour(n), mint(n), 0, ivx(n), ivy(n), ivz(n))

    parmod(1) = pydn(n)
    parmod(2) = symhc(n)
    parmod(3) = nind(n)
    parmod(4) = aby(n)

    do i=1,size(x_points)
      xx = x_points(i)
      pressure = 0
      flux = 0
      !$omp parallel private(parmod,yy,zz,hhx,hhy,hhz,bbx,bby,bbz,j,k) shared(pressure,flux,xx)
      !$omp do collapse(2)
      do j=1,size(y)
        do k=1,size(z)
          yy = y(j)
          zz = z(k)
          hhx = 0
          hhy = 0
          hhz = 0
          bbx = 0
          bby = 0
          bbz = 0

          ! External field
          call RBF_MODEL_2016(0,parmod,tilt(n),xx,yy,zz,hhx,hhy,hhz)

          ! Internal field
          call IGRF_GSW_08(xx,yy,zz,bbx,bby,bbz)

          !$omp critical
          flux = flux + abs(hhx + bbx)
          pressure = pressure + (bbx*hhx)**2 + (bby+hhy)**2 + (bbz+hhz)**2
          !$omp end critical

        end do
      end do
      !$omp end do
      !$omp end parallel

    write(p_file,"(F14.0, 2x)",advance="no") pressure/(size(y)*size(z))
    write(f_file,"(F12.0, 2x)",advance="no") flux/(size(y)*size(z))

    end do

    write(p_file,*)
    write(f_file,*)

  end do

  close(p_file)
  close(f_file)

contains

  subroutine setup_grid(x,y,z)

    implicit none

    integer :: nx, ny, nz, i
    integer :: file

    real(8) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(8), dimension(:), allocatable, intent(inout) :: x,y,z

    ! Get input parameters for grid
    open(newunit=file, file='input_parameters.txt', status='old', action='read')
    read(file, *) !skip the header line
    read(file, *) xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz
    close(file)

    allocate(x(nx))
    allocate(y(ny))
    allocate(z(nz))

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

  subroutine read_input_data(filename,year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symhc,nind,aby)

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
  end subroutine read_input_data

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

end program main
