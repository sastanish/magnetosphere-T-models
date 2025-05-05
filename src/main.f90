program main

  use parser, only : parse_input_file
  use magnetosphereModels, only : calculate_magnetic_field
  use reconnectionMetrics
  use inputOutput

  implicit none
  character(len=*) :: controlFile

  integer :: nx, ny, nz
  real(8), allocatable, dimension(:,:,:) :: x, y, z
  real(8), allocatable, dimension(:,:,:) :: bx, by, bz
  real(8), allocatable, dimension(:,:,:,:) :: reconMetrics

  ! parser outputs
  integer :: model, dipole
  logical :: outputFlags
  character(len=*) :: outputFilePath
  character(len=*) :: dataFilePath
  integer :: cells(3), bounds(6)

  ! data outputs
  integer :: dateInfo(:,:)
  real(8) :: velocity(:,:), parmod(:,:), ps(:)

  ! misc
  integer :: l

  controlFile = './control.nml'

  call parse_input_file(controlFile,model,dipole,cells,bounds,&
                        dataFilePath,outputFilePath,outputFlags)
  nx = cells(1)
  ny = cells(2)
  nz = cells(3)

  ! make grid
  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))

  call make_grid(bounds,x,y,z)
  call setup_data_file(dataFilePath,nLines,dataFileUnit)

  ! allocate field
  allocate(bx(nx,ny,nz))
  allocate(by(nx,ny,nz))
  allocate(bz(nx,ny,nz))

  ! allocate recon metrics
  allocate(reconMetrics(23,nx,ny,nz))

  if (model == 1) call read_TS04_data(dataFilePath,ps,parmod,velocity,dateinfo)
  if (model == 2) call read_TA16_data(dataFilePath,ps,parmod,velocity,dateinfo)

!  !$OMP PARALLEL DO

  do l = 1,size(ps)
    call calculate_magnetic_field(x,y,z,model,dipole,ps(l),parmod(l,:),velocity(l,:),dateinfo(l,:),bx,by,bz)
    call calculate_reconnection_metrics(x,y,z,bx,by,bz,reconMetrics)
    FILENAME = 
    call save_to_netcdf(x,y,z,bx,by,bz,reconMetrics,outputFlags,FILENAME)
  end do

!  !$OMP END PARALLEL DO

  ! allocate field
  deallocate(bx)
  deallocate(by)
  deallocate(bz)

  ! allocate recon metrics
  deallocate(reconMetrics)

contains

  subroutine read_TA16_data(dataFilePath,ps,parmod,velocity,dateinfo)

    character(:), intent(in) :: dataFilePath
    integer, intent(out), allocatable :: dateInfo(:,:)
    real(8), intent(out), allocatable :: velocity(:,:), parmod(:,:), ps(:)
    integer l,nLines, dataFileUnit
    real(8) :: BX, BY, BZ, VX, VY, VZ, DEN, TEMP, SYMH, TILT, Pdyn, &
               Bindex, Nindex, symHC
    integer :: IYEAR, IDAY, IHOUR, IMIN, IMFFLAG, ISWFLAG 

    call setup_data_file(dataFilePath,nLines,dataFileUnit)

    allocate(dateInfo(l,4))
    allocate(velocity(l,3))
    allocate(parmod(l,10))
    allocate(ps(l))

    do l = 1,nLines
      read(dataFileUnit,*) IYEAR, IDAY, IHOUR, IMIN, &
        BX, BY, BZ, VX, VY, VZ, &
        DEN, TEMP, SYMH, IMFFLAG, ISWFLAG, TILT, Pdyn, &
        Nindex, Bindex, symHC
      dateInfo(l,1) = IYEAR
      dateInfo(l,2) = IDAY
      dateInfo(l,3) = IHOUR
      dateInfo(l,4) = IMIN
      velocity(l,1) = VX
      velocity(l,2) = VY
      velocity(l,3) = VZ
      ps(l) = TILT
      parmod(l,1) = PDYN
      parmod(l,2) = symHC
      parmod(l,3) = Nindex
      parmod(l,4) = BY
      parmod(l,5) = 0
      parmod(l,6) = 0
      parmod(l,7) = 0
      parmod(l,8) = 0
      parmod(l,9) = 0
      parmod(l,10) = 0
    end do
    close(dataFileUnit)

  end subroutine read_TA16_data

  subroutine read_TS05_data(dataFilePath,ps,parmod,velocity,dateinfo)

    character(:), intent(in) :: dataFilePath
    integer, intent(out), allocatable :: dateInfo(:,:)
    real(8), intent(out), allocatable :: velocity(:,:), parmod(:,:), ps(:)
    integer l,nLines, dataFileUnit
    real(8) :: BX, BY, BZ, VX, VY, VZ, DEN, TEMP, SYMH, TILT, Pdyn, &
               W1, W2, W3, W4, W5, W6
    integer :: IYEAR, IDAY, IHOUR, IMIN, IMFFLAG, ISWFLAG 

    call setup_data_file(dataFilePath,nLines,dataFileUnit)

    allocate(dateInfo(l,4))
    allocate(velocity(l,3))
    allocate(parmod(l,10))
    allocate(ps(l))

    do l = 1,nLines
      read(dataFileUnit,*) IYEAR, IDAY, IHOUR, IMIN, &
        BX, BY, BZ, VX, VY, VZ, &
        DEN, TEMP, SYMH, IMFFLAG, ISWFLAG, TILT, Pdyn, &
        W1, W2, W3, W4, W5, W6
      dateInfo(l,1) = IYEAR
      dateInfo(l,2) = IDAY
      dateInfo(l,3) = IHOUR
      dateInfo(l,4) = IMIN
      velocity(l,1) = VX
      velocity(l,2) = VY
      velocity(l,3) = VZ
      ps(l) = TILT
      parmod(l,1) = PDYN
      parmod(l,2) = TEMP
      parmod(l,3) = BY
      parmod(l,4) = BZ
      parmod(l,5) = W1
      parmod(l,6) = W2
      parmod(l,7) = W3
      parmod(l,8) = W4
      parmod(l,9) = W5
      parmod(l,10) = W6
    end do
    close(dataFileUnit)

  end subroutine read_TS05_data


  subroutine setup_data_file(dataFilePath,nLines,dataFileUnit)
    character(:), intent(in) :: dataFilePath
    integer, intent(out) :: nLines, dataFileUnit
    integer iostat

    nlines = 1
    inquire (file=dataFilePath, iostat=iostat)
    if (iostat /= 0) then
      write (stderr, '(3a)') 'Error: file "', trim(dataFilePath), '" not found!'
    end if
    open(action="read",file=dataFilePath,iostat=iostat,newunit=dataFileUnit)
    do
      read(dataFileUnit,*,end=10)
      nLines = nLines + 1
    end do
    10 close(dataFileUnit)
    open(action="read",file=dataFilePath,iostat=iostat,newunit=dataFileUnit)

  end subroutine setup_data_file

  subroutine make_grid(bounds,x,y,z)

    real(8), intent(in), dimension(:)  :: bounds
    real(8), intent(out), dimension(:) :: x,y,z
    real(8) :: h
    integer :: i, nx,ny,nz
    
    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)

    h = (bounds(2)-bounds(1))/nx
    do i = 1,nx
      x(i) = bounds(1) + h*(i-1)
    end do

    h = (bounds(2)-bounds(1))/ny
    do i = 1,ny
      y(i) = bounds(1) + h*(i-1)
    end do

    h = (bounds(2)-bounds(1))/nz
    do i = 1,nz
      z(i) = bounds(1) + h*(i-1)
    end do

  end subroutine make_grid




end program main
