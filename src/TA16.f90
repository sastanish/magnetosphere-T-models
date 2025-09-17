module TA16

implicit none

integer, parameter :: dp=selected_real_kind(1.d0)

real(dp) :: A(23328)

contains

subroutine calculate_rbf_centers

  integer :: klat, nlat
  real(dp) :: rlow_grid, rhigh_grid
  real(dp) :: dlat_deg
  real(dp) :: Radius

  real(dp), parameter :: rh=8.0,alpha=3.0

  klat = 8
  nlat = klat + 1

  rlow_grid = 3.3
  rhigh_grid = 16.0

  dlat_deg = 90.0/(nlat-0.5)

  l = 0
  do j = 1,100
    do i = 1,nlat
      do k = 1,nlon
      end do
    end do
    l = l+1 !THIS IS ON THE BREAK CONDITION on R
  end do
      


end subroutine calculate_rbf_centers

subroutine read_ta16_parameters(filename)

  character(*), intent(in) :: filename
  integer :: ta_file

  open(newunit=ta_file,file=filename,action="read")
  do i=1,size(A)
    read(ta_file,"(G15.6)") A(i)
  end do
  close(ta_file)

end subroutine read_ta16_parameters

subroutine calculate_ta16_field
end subroutine calculate_ta16_field

end module TA16
