program main

  use TA16, only : RBF_MODEL_2016,CALCULATE_RBF_CENTERS,READ_TA16_PARS
  use TS05, only : T04_s
  use geopack, only : RECALC_08, IGRF_GSW_08
  use inputOutput, only : save_field_to_netcdf
  !$ use omp_lib

  implicit none

  real(8), dimension(:,:,:), allocatable :: Bx, By, Bz
  real(8), dimension(:), allocatable :: x,y,z

  ! Input File params
  integer, dimension(:), allocatable :: year,day,hour,mint
  real(8), dimension(:), allocatable :: ivx,ivy,ivz,tilt,pydn,symhc,nind,aby

  ! vars
  integer :: n, i, j, k

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

  call setup_grid(x,y,z,Bx,By,Bz)
  call read_input_data('input_data.lst',year,day,hour,mint,ivx,ivy,ivz,tilt,pydn,symhc,nind,aby)
  call read_TA16_pars
  call calculate_rbf_centers

  ivy = ivy + 29.78 !velocity correction

  do n = start_ind,end_ind

    call RECALC_08(year(n), day(n), hour(n), mint(n), 0, ivx(n), ivy(n), ivz(n))

    !$OMP PARALLEL PRIVATE(parmod,xx,yy,zz,hhx,hhy,hhz,bbx,bby,bbz,i,j,k) SHARED(n,x,y,z,Bx,By,Bz)
    parmod(1) = pydn(n)
    parmod(2) = symhc(n)
    parmod(3) = nind(n)
    parmod(4) = aby(n)

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
          call RBF_MODEL_2016(0,parmod,tilt(n),xx,yy,zz,hhx,hhy,hhz)

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

    ! Write to file
    write( str_ind, '(I4)' ) n
    print *, 'my name is '//'output_'//trim(adjustl(str_ind))//'.nc'
    call save_field_to_netcdf(x,y,z,Bx,By,Bz,'output_'//trim(adjustl(str_ind))//'.nc')

  end do

  deallocate(Bx)
  deallocate(By)
  deallocate(Bz)
  deallocate(x)
  deallocate(y)
  deallocate(z)

contains

  subroutine setup_grid(x,y,z,Bx,By,Bz)

    implicit none

    integer :: nx, ny, nz, i
    integer :: file

    real(8) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(8), dimension(:), allocatable, intent(inout) :: x,y,z
    real(8), dimension(:,:,:), allocatable, intent(inout) :: Bx, By, Bz

    ! Get input parameters for grid
    open(newunit=file, file='input_parameters.txt', status='old', action='read')
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

