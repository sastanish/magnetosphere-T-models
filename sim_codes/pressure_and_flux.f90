program main

use TA16, only : RBF_MODEL_2016,CALCULATE_RBF_CENTERS,READ_TA16_PARS
use geopack, only : RECALC_08, IGRF_GSW_08

implicit none

! Image Handling
integer :: img, nimgs, io_img, io_stat

! Input Data
integer :: flux_file, pres_file, omni_file
integer :: date(4)[*]
real(8) :: field_data(9)[*],index_data(3)[*],derived_parameters(8)[*],parmod(10)

! Grid
integer :: Area,Nx,Ny,Nz,Nypc
real(8) :: xmin, xmax, ymin, ymax, zmin, zmax, xx,yy,zz

! Computation 
real(8) :: flux[*], pres[*]
real(8) :: bx,by,bz,hx,hy,hz 

! Dummy Vars
integer :: i,j,k, line

! Parallel setup
io_img = 1
img = this_image()
nimgs = num_images()

if (img==io_img) then
  print *, "Initalizing enviorment on ", nimgs
end if

sync all
print *, "Proc: ", img, " online"

! Parameters TODO: write input script for better ingestion
Nx = 10
Ny = 300
Nz = 300
xmin = -12
xmax = -2
ymin = -6
ymax = 6
zmin = -6
zmax = 6
Area = Ny*Nz

Nypc = Ny/nimgs !Careful that Nimgs divides Ny else we will miss points. TODO: Make error msg

! Setup output files
open(newunit=flux_file, file="flux.lst")
open(newunit=pres_file, file="pressure.lst")

! Open input file on image 1 and write headers
if (img == io_img) then
  open(newunit=omni_file, file="input_data.lst")
  write(pres_file,"(A)") "#Avg Pressure across y-z plane in magnetotail"
  write(flux_file,"(A)") "#Avg |Bx| across y-z plane in magnetotail"
  write(pres_file,"(A)") "#FORMAT:"
  write(flux_file,"(A)") "#FORMAT:"
  write(pres_file,"(A, F6.2, 4x)",advance="no") "YEAR    DAY    HOUR    MIN    x:", xmin
  write(flux_file,"(A, F6.2, 4x)",advance="no") "YEAR    DAY    HOUR    MIN    x:", xmin
  do i = 1,nx+1
    write(pres_file,"(F6.2, 6x)",advance="no")  xmin + i
    write(flux_file,"(F6.2, 6x)",advance="no")  xmin + i
  end do
  write(pres_file,*)
  write(flux_file,*)
end if

! Setup grid
call read_TA16_pars
call calculate_rbf_centers

! Computational loop
io_stat = 0
line = 1
do while ( io_stat == 0 )
  ! Read in data
  if (img == io_img) then
    read(omni_file,*,iostat=io_stat) date,field_data,index_data,derived_parameters
    write(flux_file,"(I4, 6x, I3, 6x, I4, 6x, I3, 6x)",advance="no") date
    write(pres_file,"(I4, 6x, I3, 6x, I4, 6x, I3, 6x)",advance="no") date
    !! Date = [Year, Day, Hour, Min]
    !! Field_data = [Bx,By,Bz,Vx,Vy,Vz,den,temp,pydn]
    !! Index_data = [Ae index, SymH, IMF_flag]
    !! Derived_parameters = [sw,tilt,rp,<bx>,<by>,<bz>,Nindex,SymHc]
    print *, "Line number: ", line
  end if

  sync all ! Wait all for read

  ! Broadcast data to arrays
  call co_broadcast(io_stat,io_img)
  call co_broadcast(date,io_img)
  call co_broadcast(field_data,io_img)
  call co_broadcast(derived_parameters,io_img)

  ! Start calculation                                                  Velocity correction
  call RECALC_08(date(1), date(2), date(3), date(4), 0, field_data(4), field_data(5)+29.78, field_data(6))

  parmod(1) = field_data(9) !Pydn (pressure)
  parmod(2) = derived_parameters(8) !SymHc (Modified SymH)
  parmod(3) = derived_parameters(7) !N index
  parmod(4) = derived_parameters(5) !<By>

  do i = 0,Nx

    flux = 0
    pres = 0

    xx = xmin + i

    do j = (img-1)*Nypc+1,img*Nypc
    yy = ymin + j*(ymax - ymin)/Ny
    do k = 1,Nz

      zz = zmin + k*(zmax - zmin)/Nz

      ! Compute field
      call RBF_MODEL_2016(0,parmod,derived_parameters(2),xx,yy,zz,hx,hy,hz)
      call IGRF_GSW_08(xx,yy,zz,bx,by,bz)

      flux = flux + abs(hx+bx)
      pres = pres + (hx+bx)**2+(hy+by)**2+(hz+bz)**2

    end do
    end do

    sync all
    call co_sum(pres, io_img)
    call co_sum(flux, io_img)

    if (img == io_img) then
      write(flux_file,"(F12.0, 2x)",advance="no") flux/Area
      write(pres_file,"(F12.0, 2x)",advance="no") pres/Area
    end if

  end do

  if (img == io_img) then
    write(flux_file,*)
    write(pres_file,*)
  end if

  line = line + 1

end do

close(flux_file)
close(pres_file)
if (img==io_img) close(omni_file)

end program main
