program compute

  use TA16, only : RBF_MODEL_2016
  use geopack, only : RECALC_08, IGRF_GSW_08
  use inputOutput, only : save_field_to_netcdf
  !$ use omp_lib

  implicit none
  integer :: nx, ny, nz
  integer :: i, j, k
  integer :: ip, id, stat

  integer :: fileind
  character(4) :: str_ind

  real(8) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(8) :: xx,yy,zz,bbx,bby,bbz,hhx,hhy,hhz !dummy variables
  real(8) :: x(nx), y(ny), z(nz)
  real(8) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
  real(8) :: parmod(10), ps

  ! Input File params
  integer :: year, day, hour, mint, aeind, symh, imf
  real(8) :: ibx, iby, ibz, ivx, ivy, ivz, den, temp, p, sw, tilt, rp
  real(8) :: nind, symhc, abx, aby, abz

  ! Get input parameters for grid
  open(newunit=ip, file='input_parameters.txt', status='old', action='read')
  read(ip, *) !skip the header line
  read(ip, *) xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz
  close(ip)

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
                             den, temp, p, symh, imf, sw, tilt, rp, abx, aby, abz,&
                             nind, symhc
    if (stat .le. 0) exit !checking for end of file
    parmod(1) = p
    parmod(2) = symhc
    parmod(3) = nind
    parmod(4) = aby

    call RECALC_08(year, day, hour, mint, 0, ivx, ivy, ivz)

    !$OMP PARALLEL SHARED(Bx,By,Bz,nx,ny,nz,parmod,ps) PRIVATE(bbx,bby,bbz,xx,yy,zz,i,j,k)
    !$OMP DO COLLAPSE(3)
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz

          xx = x(i)
          yy = y(j)
          zz = z(k)

          ! External field
          call RBF_MODEL_2016(0,parmod,ps,xx,yy,zz,hhx,hhy,hhz)

          ! Internal field
          call IGRF_GSW_08(xx,yy,zz,bbx,bby,bbz)

          Bx(i,j,k) = bbx + hhx
          By(i,j,k) = bby + hhy
          Bz(i,j,k) = bbz + hhz

        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Write to file
    write( str_ind, '(I4)' ) fileind
    call save_field_to_netcdf(x,y,z,Bx,By,Bz,'output_'//str_ind//'.nc')

  end do

end program main
