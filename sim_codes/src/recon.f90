program main

  use inputOutput, only : save_rate_to_netcdf, load_field_from_netcdf

  implicit none

  real(8), dimension(:,:,:), allocatable :: Bx,By,Bz,magB,Fx,Fy,Fz,magF
  real(8), dimension(:,:,:), allocatable :: Jx,Jy,Jz,Bfx,Bfy,Bfz
  real(8), dimension(:,:,:), allocatable :: lambda,alpha,rate

  real(8), dimension(:), allocatable :: x,y,z

  integer :: i,j,k,nx,ny,nz
  real(8) :: hx,hy,hz
  real(8) :: Gax,Gay,Gaz,omega1,omega2,BLG,cBfx,cBfy,cBfz,rate_x,rate_y,rate_z

  !File handling
  integer :: fileind, start_ind, end_ind
  character(4) :: str_ind

  GET_COMMAND_ARGUMENT(1,start_ind)
  GET_COMMAND_ARGUMENT(2,end_ind)
  do fileind = start_ind,end_ind

    write( str_ind, '(I4)' ) fileind
    print *, 'reading: '//'output_'//trim(adjustl(str_ind))//'.nc'
  
    call load_field_from_netcdf('output_'//trim(adjustl(str_ind))//'.nc',x,y,z,Bx,By,Bz)

    nx = size(x)
    ny = size(y)
    nz = size(z)

    ! ALLOCATIONS
    allocate(magB(nx,ny,nz))
    allocate(Fx(nx,ny,nz))
    allocate(Fy(nx,ny,nz))
    allocate(Fz(nx,ny,nz))
    allocate(magF(nx,ny,nz))
    allocate(Jx(nx,ny,nz))
    allocate(Jy(nx,ny,nz))
    allocate(Jz(nx,ny,nz))
    allocate(Bfx(nx,ny,nz))
    allocate(Bfy(nx,ny,nz))
    allocate(Bfz(nx,ny,nz))
    allocate(alpha(nx,ny,nz))
    allocate(lambda(nx,ny,nz))
    allocate(rate(nx,ny,nz))

    !!!!!!!!!!!!!!!!!!
    ! Calculate rate !
    !!!!!!!!!!!!!!!!!!

    !$OMP Parallel public(Bx,By,Bz,magB,Fx,Fy,Fz,Jx,Jy,Jz,Bfx,Bfy,Bfz,lambda,alpha) private(i,j,k)
    !$OMP DO COLLAPSE(3)
    do k=2,nz-1
    do j=2,ny-1
    do i=2,nx-1
      ! |B|
      magB(i,j,k) = sqrt(Bx(i,j,k)**2 + By(i,j,k)**2 + Bz(i,j,k)**2)
      ! F = u x B
      Fx(i,j,k) = uy(i,j,k)*Bz(i,j,k) - uz(i,j,k)*By(i,j,k) 
      Fy(i,j,k) = uz(i,j,k)*Bx(i,j,k) - ux(i,j,k)*Bz(i,j,k) 
      Fz(i,j,k) = ux(i,j,k)*By(i,j,k) - uy(i,j,k)*Bx(i,j,k) 
      magF(i,j,k) = sqrt(Fx(i,j,k)**2 + Fy(i,j,k)**2 + Fz(i,j,k)**2)
      ! J = Curl(B)
      Jx(i,j,k) = 1/(2*hy) * (Bz(i,j+1,k) - Bz(i,j-1,k)) - 1/(2*hz) * (By(i,j,k+1) - By(i,j,k-1))
      Jy(i,j,k) = 1/(2*hz) * (Bx(i,j,k+1) - Bx(i,j,k-1)) - 1/(2*hx) * (Bz(i+1,j,k) - Bz(i-1,j,k))
      Jz(i,j,k) = 1/(2*hx) * (By(i+1,j,k) - By(i-1,j,k)) - 1/(2*hy) * (Bx(i,j+1,k) - Bx(i,j-1,k))
      magJ(i,j,k) = sqrt(Jx(i,j,k)**2 + Jy(i,j,k)**2 + Jz(i,j,k)**2)
      ! B_f = B x F/|F|
      Bfx(i,j,k) = (By(i,j,k)*Fz(i,j,k) - Fz(i,j,k)*By(i,j,k) )/magF(i,j,k)
      Bfy(i,j,k) = (Bz(i,j,k)*Fx(i,j,k) - Fx(i,j,k)*Bz(i,j,k) )/magF(i,j,k)
      Bfz(i,j,k) = (Bx(i,j,k)*Fy(i,j,k) - Fy(i,j,k)*Bx(i,j,k) )/magF(i,j,k)
      ! alpha = j * b / |B|^2
      alpha(i,j,k) = (Jx(i,j,k)*Bx(i,j,k) + Jy(i,j,k)*By(i,j,k) + Jz(i,j,k)*Bz(i,j,k) )/magB(i,j,k)**2
      !lamb = j * B_f / |B|^2
      lambda(i,j,k) = (Jx(i,j,k)*Bfx(i,j,k) + Jy(i,j,k)*Bfy(i,j,k) + Jz(i,j,k)*Bfz(i,j,k) )/magB(i,j,k)**2
    end do
    end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP Parallel public(rate) private(i,j,k,cBfx,cBfy,cBfz,omega1,omega2,rate_x,rate_y,rate_z,Gax,Gay,Gaz)
    !$OMP DO COLLAPSE(3)
    do k=3,nz-2
    do j=3,ny-2
    do i=3,nx-2
      ! Curl(B_f)
      cBfx = 1/(2*hy) * (Bfz(i,j+1,k) - Bfz(i,j-1,k)) - 1/(2*hz) * (Bfy(i,j,k+1) - Bfy(i,j,k-1))
      cBfy = 1/(2*hz) * (Bfx(i,j,k+1) - Bfx(i,j,k-1)) - 1/(2*hx) * (Bfz(i+1,j,k) - Bfz(i-1,j,k))
      cBfz = 1/(2*hx) * (Bfy(i+1,j,k) - Bfy(i-1,j,k)) - 1/(2*hy) * (Bfx(i,j+1,k) - Bfx(i,j-1,k))
      ! w1 = Curl(B_f) * F/|F|
       omega1 = ( cBfx * Fx + cBfy * Fy + cBfz * Fz ) / magF(i,j,k)
      ! w2 = curl(B_f) * B_f / |B|^2
       omega2 = ( cBfx * Bfx + cBfy * Bfy + cBfz * Bfz ) / magB(i,j,k)**2
      ! GLB = Grad(lambda) * B
      BLG = ( 1/(2*hx) * (lambda(i+1,j,k) - lambda(i-1,j,k)) * Bx(i,j,k) &
            + 1/(2*hy) * (lambda(i,j+1,k) - lambda(i,j-1,k)) * By(i,j,k) &
            + 1/(2*hz) * (lambda(i,j,k+1) - lambda(i,j,k-1)) * Bz(i,j,k) )
      ! Ga = Grad(alpha)
      Gax = 1/(2*hx) * (alpha(i+1,j,k) - alpha(i-1,j,k))
      Gay = 1/(2*hy) * (alpha(i,j+1,k) - alpha(i,j-1,k))
      Gaz = 1/(2*hz) * (alpha(i,j,k+1) - alpha(i,j,k-1))
      ! Rate_vec = -1/|B| ( (lamb w1 - G(lamb) * B) F/|F| + lamb (alph + w2)B_f + G(alph) x B)
      rate_x = - 1/magB(i,j,k) * ( (lambda(i,j,k) * omega1 - BLG) * Fx(i,j,k) / magF(i,j,k) &
                                 + lambda(i,j,k) * (alpha(i,j,k) + omega2) * Bfx(i,j,k) &
                                 + ( Gay * Bz(i,j,k) - Gaz * By(i,j,k) ) &
                                 )
      rate_y = - 1/magB(i,j,k) * ( (lambda(i,j,k) * omega1 - BLG) * Fy(i,j,k) / magF(i,j,k) &
                                 + lambda(i,j,k) * (alpha(i,j,k) + omega2) * Bfy(i,j,k) &
                                 + ( Gaz * Bx(i,j,k) - Gax * Bz(i,j,k) ) &
                                 )
      rate_z = - 1/magB(i,j,k) * ( (lambda(i,j,k) * omega1 - BLG) * Fz(i,j,k) / magF(i,j,k) &
                                 + lambda(i,j,k) * (alpha(i,j,k) + omega2) * Bfz(i,j,k) &
                                 + ( Gax * By(i,j,k) - Gay * Bx(i,j,k) ) &
                                 )
      ! Rate = sqrt( ratex^2 + ratey^2 + ratez^2)
      rate(i,j,k) = sqrt(rate_x**2 + rate_y**2 + rate_z**2)
    end do
    end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Write to file
    print *, 'Writing '//'rate_'//trim(adjustl(str_ind))//'.nc'
    call save_rate_to_netcdf(x,y,z,rate,'rate_'//trim(adjustl(str_ind))//'.nc')

  end do

end program main
