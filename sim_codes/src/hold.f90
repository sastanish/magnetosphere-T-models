program main

  use compute, only : run_TA16, run_igrf_dipole
  use inputOutput, only : save_field_to_netcdf

  implicit none
  integer :: nx, ny, nz
  integer :: i, j, k
  integer :: ip, id, stat

  integer :: fileind
  character(4) :: str_ind

  real(8) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(8), dimension(:), allocatable :: x,y,z
  real(8), dimension(:,:,:), allocatable :: Bx, By, Bz
  real(8), dimension(:,:,:), allocatable :: Hx, Hy, Hz
  real(8) :: parmod(10)

  ! Input File params
  integer :: year, day, hour, mint, aeind, symh, imf
  real(8) :: ibx, iby, ibz, ivx, ivy, ivz, den, temp, p, sw, tilt, rp
  real(8) :: nind, symhc, abx, aby, abz

  ! Get input parameters for grid
  open(newunit=ip, file='input_parameters.txt', status='old', action='read')
  read(ip, *) !skip the header line
  read(ip, *) xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz
  close(ip)

  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))
  allocate(Bx(nx,ny,nz))
  allocate(By(nx,ny,nz))
  allocate(Bz(nx,ny,nz))
  allocate(Hx(nx,ny,nz))
  allocate(Hy(nx,ny,nz))
  allocate(Hz(nx,ny,nz))


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

  ! Now read the next line of data from the input data
  open(newunit=id, file='input_data.txt', status='old', action='read')
  fileind = 1
  do 
    read(id, *, iostat=stat) year, day, hour, mint, ibx, iby, ibz, ivx, ivy, ivz, &
                             den, temp, p, aeind, symh, imf, sw, tilt, rp, abx, aby, abz,&
                             nind, symhc
    if (stat < 0) exit !checking for end of file
    parmod(1) = p
    parmod(2) = symhc
    parmod(3) = nind
    parmod(4) = aby 
    ivy = ivx + 29.78

    ! External field
    call run_igrf_dipole(year, day, hour, mint, 0, ivx, ivy, ivz, &
                         x, y, z, Bx, By, Bz, nx, ny, nz)

    print *, 'do we get to here?'
    ! Internal field
    call run_TA16(parmod,tilt,x,y,z,Hx,Hy,Hz,nx,ny,nz)

    print *, 'How bout here?'
    Bx = Bx + Hx
    By = By + Hy
    Bz = Bz + Hz

    ! Write to file
    write( str_ind, '(I4)' ) fileind
    print *, 'my name is '//'output_'//trim(adjustl(str_ind))//'.nc'
    call save_field_to_netcdf(x,y,z,Bx,By,Bz,'output_'//trim(adjustl(str_ind))//'.nc')
    fileind = fileind + 1

  end do

  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(Bx)
  deallocate(By)
  deallocate(Bz)
  deallocate(Hx)
  deallocate(Hy)
  deallocate(Hz)


end program main
