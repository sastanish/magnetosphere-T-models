module reconnectionMetrics

  implicit none

contains

  subroutine compute_all_metrics(x,y,z,bx,by,bz,Mout)
    
    real(8), intent(in), dimension(:) :: x,y,z
    real(8), intent(in), dimension(:,:,:) :: bx,by,bz

    real(8), intent(out), dimension(:,:,:,:) :: Mout
    ! First dimension is [jx, jy, jz, fx, fy, fz, bfpx, bfpy, bfpz, alpha, lambda,&
    !                     mag_c2_t1, c2_t1_x, c2_t1_y, c2_t1_z,
    !                     mag_c2_t2, c2_t2_x, c2_t2_y, c2_t2_z,
    !                     mag_c2_t3, c2_t3_x, c2_t3_y, c2_t3_z]

    call compute_forces(x,y,z,bx,by,bz,Mout(1,:,:,:),&
                                       Mout(2,:,:,:),&
                                       Mout(3,:,:,:),&
                                       Mout(4,:,:,:),&
                                       Mout(5,:,:,:),&
                                       Mout(6,:,:,:) &
                                       )

    call compute_parameters(x,y,z,bx,by,bz,Mout(1,:,:,:),&
                                           Mout(2,:,:,:),&
                                           Mout(3,:,:,:),&
                                           Mout(4,:,:,:),&
                                           Mout(5,:,:,:),&
                                           Mout(6,:,:,:),&
                                           Mout(7,:,:,:),&
                                           Mout(8,:,:,:),&
                                           Mout(9,:,:,:),&
                                           Mout(10,:,:,:),&
                                           Mout(11,:,:,:)&
                                           )

    call compute_c2_t1(x,y,z,bx,by,bz,Mout(4,:,:,:),&
                                      Mout(5,:,:,:),&
                                      Mout(6,:,:,:),&
                                      Mout(7,:,:,:),&
                                      Mout(8,:,:,:),&
                                      Mout(9,:,:,:),&
                                      Mout(11,:,:,:),&
                                      Mout(12,:,:,:),&
                                      Mout(13,:,:,:),&
                                      Mout(14,:,:,:),&
                                      Mout(15,:,:,:)&
                                      )

    call compute_c2_t2(x,y,z,bx,by,bz,Mout(4,:,:,:),&
                                      Mout(5,:,:,:),&
                                      Mout(6,:,:,:),&
                                      Mout(7,:,:,:),&
                                      Mout(8,:,:,:),&
                                      Mout(9,:,:,:),&
                                      Mout(10,:,:,:),&
                                      Mout(11,:,:,:),&
                                      Mout(16,:,:,:),&
                                      Mout(17,:,:,:),&
                                      Mout(18,:,:,:),&
                                      Mout(19,:,:,:)&
                                      )

    call compute_c2_t3(x,y,z,bx,by,bz,Mout(10,:,:,:),Mout(20,:,:,:),&
                                                     Mout(21,:,:,:),&
                                                     Mout(22,:,:,:),&
                                                     Mout(23,:,:,:)&
                                                     )
    

  end subroutine compute_all_metrics

  subroutine compute_forces(x,y,z,bx,by,bz,jx,jy,jz,fx,fy,fz)
    real(8), intent(in), dimension(:) :: x,y,z
    real(8), intent(in), dimension(:,:,:) :: bx,by,bz

    real(8), dimension(:,:,:), intent(out) :: jx,jy,jz
    real(8), dimension(:,:,:), intent(out) :: fx,fy,fz

    integer nx,ny,nz
    integer :: ix, iy, iz
    real(8) :: hx, hy, hz, jjx, jjy, jjz
    integer :: r(3,3) !9-point stencil

    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

    ! compute internal current via second order central diff
    ! and lorrentz force
    do iz = 1,nz
      if (iz==1) then; r(:,3) = [1, 2, 3];
      else if (iz==nz) then; r(:,3) = [nz-2, nz-1, nz];
      else; r(:,3) = [iz-1, iz, iz+1];
      end if

      do iy = 1,ny
        if (iy==1) then; r(:,2) = [1, 2, 3];
        else if (iy==ny) then; r(:,2) = [ny-2, ny-1, ny];
        else; r(:,2) = [iy-1, iy, iy+1];
        end if

        do ix = 1,nx
          if (ix==1) then; r(:,1) = [1, 2, 3];
          else if (ix==nx) then; r(:,1) = [nx-2, nx-1, nx];
          else; r(:,1) = [ix-1, ix, ix+1];
          end if

          ! compute J=curl(B)
          jjx = (bz(r(2,1),r(3,2),r(2,3))     &
                -bz(r(2,1),r(1,2),r(2,3)))/(2*hy) &
              - (by(r(2,1),r(2,2),r(3,3))     &
                -by(r(2,1),r(2,2),r(1,3)))/(2*hz)

          jjy = (bx(r(2,1),r(2,2),r(3,3))     &
                -bx(r(2,1),r(2,2),r(1,3)))/(2*hz) &
              - (bz(r(3,1),r(2,2),r(2,3))     &
                -bz(r(1,1),r(2,2),r(2,3)))/(2*hx) 

          jjz = (by(r(3,1),r(2,2),r(2,3))     &
                -by(r(1,1),r(2,2),r(2,3)))/(2*hx) &
              - (bx(r(2,1),r(3,2),r(2,3))     &
                -bx(r(2,1),r(1,2),r(2,3)))/(2*hy)

          ! compute F= J cross B
          fx(ix,iy,iz) = jjy * bz(ix,iy,iz) - jjz * by(ix,iy,iz)
          fy(ix,iy,iz) = jjz * bx(ix,iy,iz) - jjx * bz(ix,iy,iz)
          fz(ix,iy,iz) = jjx * by(ix,iy,iz) - jjy * bx(ix,iy,iz)

          ! fill in j
          jx(ix,iy,iz) = jjx
          jy(ix,iy,iz) = jjy
          jz(ix,iy,iz) = jjz

        end do
      end do
    end do

  end subroutine compute_forces

  subroutine compute_parameters(x,y,z,bx,by,bz,jx,jy,jz,fx,fy,fz,&
                                bfpx,bfpy,bfpz,alpha,lambda)

    real(8), intent(in), dimension(:) :: x,y,z
    real(8), intent(in), dimension(:,:,:) :: bx,by,bz
    real(8), intent(in), dimension(:,:,:) :: jx,jy,jz
    real(8), intent(in), dimension(:,:,:) :: fx,fy,fz

    real(8), intent(out), dimension(:,:,:) :: bfpx,bfpy,bfpz,alpha,lambda

    real(8), dimension(:,:,:), allocatable :: magF, magB
    integer nx,ny,nz

    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)

    allocate(magF(nx,ny,nz))
    allocate(magB(nx,ny,nz))

    magF = sqrt(fx**2 + fy**2 + fz**2)
    magB = sqrt(bx**2 + by**2 + bz**2)

    bfpx = (fz * by - fy * bz)/magF
    bfpy = (fx * bz - fz * bx)/magF
    bfpz = (fy * bx - fx * by)/magF

    lambda = (jx*bfpx + jy*bfpy + jz*bfpz)/magB**2
    alpha = (jx*bx + jy*by + jz*bz)/magB**2

    deallocate(magF)
    deallocate(magB)

  end subroutine compute_parameters

  subroutine compute_c2_t1(x,y,z,bx,by,bz,fx,fy,fz,bfpx,bfpy,bfpz,lambda,&
                           mag,c2_t1_x,c2_t1_y,c2_t1_z)

    real(8), intent(in), dimension(:) :: x,y,z
    real(8), intent(in), dimension(:,:,:) :: bx,by,bz
    real(8), intent(in), dimension(:,:,:) :: fx,fy,fz
    real(8), intent(in), dimension(:,:,:) :: bfpx,bfpy,bfpz
    real(8), intent(in), dimension(:,:,:) :: lambda

    real(8), intent(out), dimension(:,:,:) :: mag,c2_t1_x,c2_t1_y,c2_t1_z

    integer :: ix, iy, iz, r(3,3)
    real(8) :: hx, hy, hz
    real(8) :: curl_bfpx, curl_bfpy, curl_bfpz

    real(8), dimension(:,:,:), allocatable :: magF, magB
    integer nx,ny,nz

    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)

    allocate(magF(nx,ny,nz))
    allocate(magB(nx,ny,nz))

    magF = sqrt(fx**2 + fy**2 + fz**2)
    magB = sqrt(bx**2 + by**2 + bz**2)

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

    do iz = 1,nz
      if (iz==1) then; r(:,3) = [1, 2, 3];
      else if (iz==nz) then; r(:,3) = [nz-2, nz-1, nz];
      else; r(:,3) = [iz-1, iz, iz+1];
      end if

      do iy = 1,ny
        if (iy==1) then; r(:,2) = [1, 2, 3];
        else if (iy==ny) then; r(:,2) = [ny-2, ny-1, ny];
        else; r(:,2) = [iy-1, iy, iy+1];
        end if

        do ix = 1,nx
          if (ix==1) then; r(:,1) = [1, 2, 3];
          else if (ix==nx) then; r(:,1) = [nx-2, nx-1, nx];
          else; r(:,1) = [ix-1, ix, ix+1];
          end if

          !curl(B_fperp)
          curl_bfpx = (bfpz(r(2,1),r(3,2),r(2,3))     &
                      -bfpz(r(2,1),r(1,2),r(2,3)))/(2*hy) &
                    - (bfpy(r(2,1),r(2,2),r(3,3))     &
                      -bfpy(r(2,1),r(2,2),r(1,3)))/(2*hz)

          curl_bfpy = (bfpx(r(2,1),r(2,2),r(3,3))     &
                      -bfpx(r(2,1),r(2,2),r(1,3)))/(2*hz) &
                    - (bfpz(r(3,1),r(2,2),r(2,3))     &
                      -bfpz(r(1,1),r(2,2),r(2,3)))/(2*hx) 

          curl_bfpz = (bfpy(r(3,1),r(2,2),r(2,3))     &
                      -bfpy(r(1,1),r(2,2),r(2,3)))/(2*hx) &
                    - (bfpx(r(2,1),r(3,2),r(2,3))     &
                      -bfpx(r(2,1),r(1,2),r(2,3)))/(2*hy)

          ! (lambda * curl(B_fperp) dot F/|F| - Grad(lambda) dot B)/|B|
          mag(ix,iy,iz) = (&
          lambda(ix,iy,iz)*(curl_bfpx*fx(ix,iy,iz) + curl_bfpy*fy(ix,iy,iz) + curl_bfpz*fz(ix,iy,iz))/magF(ix,iy,iz)&
        - ((lambda(r(3,1),r(2,2),r(2,3)) - lambda(r(1,1),r(2,2),r(2,3)))/(2*hx)*bx(ix,iy,iz) &
          +(lambda(r(2,1),r(3,2),r(2,3)) - lambda(r(2,1),r(1,2),r(2,3)))/(2*hy)*by(ix,iy,iz) &
          +(lambda(r(2,1),r(2,2),r(3,3)) - lambda(r(2,1),r(2,2),r(1,3)))/(2*hz)*bz(ix,iy,iz) &
          )&
          )/magB(ix,iy,iz)
        end do
      end do
    end do

    c2_t1_x = -1.0*mag * fx/magF
    c2_t1_y = -1.0*mag * fy/magF
    c2_t1_z = -1.0*mag * fz/magF

    mag = abs(mag)

    deallocate(magF)
    deallocate(magB)

  end subroutine compute_c2_t1

  subroutine compute_c2_t2(x,y,z,bx,by,bz,fx,fy,fz,bfpx,bfpy,bfpz,lambda,alpha,&
                           mag,c2_t2_x,c2_t2_y,c2_t2_z)

    real(8), intent(in), dimension(:) :: x,y,z
    real(8), intent(in), dimension(:,:,:) :: bx,by,bz
    real(8), intent(in), dimension(:,:,:) :: fx,fy,fz
    real(8), intent(in), dimension(:,:,:) :: bfpx,bfpy,bfpz
    real(8), intent(in), dimension(:,:,:) :: alpha,lambda

    real(8), intent(out), dimension(:,:,:) :: mag,c2_t2_x,c2_t2_y,c2_t2_z

    integer :: ix, iy, iz, r(3,3)
    real(8) :: hx, hy, hz
    real(8) :: curl_bfpx, curl_bfpy, curl_bfpz

    real(8), dimension(:,:,:), allocatable :: magF, magB
    integer nx,ny,nz

    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)

    allocate(magF(nx,ny,nz))
    allocate(magB(nx,ny,nz))

    magF = sqrt(fx**2 + fy**2 + fz**2)
    magB = sqrt(bx**2 + by**2 + bz**2)

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

    do iz = 1,nz
      if (iz==1) then; r(:,3) = [1, 2, 3];
      else if (iz==nz) then; r(:,3) = [nz-2, nz-1, nz];
      else; r(:,3) = [iz-1, iz, iz+1];
      end if

      do iy = 1,ny
        if (iy==1) then; r(:,2) = [1, 2, 3];
        else if (iy==ny) then; r(:,2) = [ny-2, ny-1, ny];
        else; r(:,2) = [iy-1, iy, iy+1];
        end if

        do ix = 1,nx
          if (ix==1) then; r(:,1) = [1, 2, 3];
          else if (ix==nx) then; r(:,1) = [nx-2, nx-1, nx];
          else; r(:,1) = [ix-1, ix, ix+1];
          end if

          !curl(B_fperp)
          curl_bfpx = (bfpz(r(2,1),r(3,2),r(2,3))     &
                      -bfpz(r(2,1),r(1,2),r(2,3)))/(2*hy) &
                    - (bfpy(r(2,1),r(2,2),r(3,3))     &
                      -bfpy(r(2,1),r(2,2),r(1,3)))/(2*hz)

          curl_bfpy = (bfpx(r(2,1),r(2,2),r(3,3))     &
                      -bfpx(r(2,1),r(2,2),r(1,3)))/(2*hz) &
                    - (bfpz(r(3,1),r(2,2),r(2,3))     &
                      -bfpz(r(1,1),r(2,2),r(2,3)))/(2*hx) 

          curl_bfpz = (bfpy(r(3,1),r(2,2),r(2,3))     &
                      -bfpy(r(1,1),r(2,2),r(2,3)))/(2*hx) &
                    - (bfpx(r(2,1),r(3,2),r(2,3))     &
                      -bfpx(r(2,1),r(1,2),r(2,3)))/(2*hy)

          ! lambda (curl(B_fperp) dot B_fperp/|B|^2 + alpha)/|B|
          mag(ix,iy,iz) = lambda(ix,iy,iz)*( &
                  ( curl_bfpx*bfpx(ix,iy,iz) &
                  + curl_bfpy*bfpy(ix,iy,iz) &
                  + curl_bfpz*bfpz(ix,iy,iz) )/magB(ix,iy,iz)**2 &
              + alpha(ix,iy,iz))/magB(ix,iy,iz)

        end do
      end do
    end do

    c2_t2_x = -1.0*mag * bfpx
    c2_t2_y = -1.0*mag * bfpy
    c2_t2_z = -1.0*mag * bfpz

    mag = abs(mag)

    deallocate(magF)
    deallocate(magB)

  end subroutine compute_c2_t2

  subroutine compute_c2_t3(x,y,z,bx,by,bz,alpha,mag,c2_t3_x,c2_t3_y,c2_t3_z)

    real(8), intent(in), dimension(:) :: x,y,z
    real(8), intent(in), dimension(:,:,:) :: bx,by,bz
    real(8), intent(in), dimension(:,:,:) :: alpha

    real(8), intent(out), dimension(:,:,:) :: mag,c2_t3_x,c2_t3_y,c2_t3_z

    integer :: ix, iy, iz, r(3,3)
    real(8) :: hx, hy, hz
    real(8) :: dxa, dya, dza
    real(8), dimension(:,:,:), allocatable :: magB
    integer nx,ny,nz

    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)

    allocate(magB(nx,ny,nz))

    magB = sqrt(bx**2 + by**2 + bz**2)

    ! We assume that the grid spacing is uniform
    hx = abs(x(2) - x(1))
    hy = abs(y(2) - y(1))
    hz = abs(z(2) - z(1))

    do iz = 1,nz
      if (iz==1) then; r(:,3) = [1, 2, 3];
      else if (iz==nz) then; r(:,3) = [nz-2, nz-1, nz];
      else; r(:,3) = [iz-1, iz, iz+1];
      end if

      do iy = 1,ny
        if (iy==1) then; r(:,2) = [1, 2, 3];
        else if (iy==ny) then; r(:,2) = [ny-2, ny-1, ny];
        else; r(:,2) = [iy-1, iy, iy+1];
        end if

        do ix = 1,nx
          if (ix==1) then; r(:,1) = [1, 2, 3];
          else if (ix==nx) then; r(:,1) = [nx-2, nx-1, nx];
          else; r(:,1) = [ix-1, ix, ix+1];
          end if

          dxa = (alpha(r(3,1),r(2,2),r(2,3))-alpha(r(1,1),r(2,2),r(2,3)))/(2*hx)
          dya = (alpha(r(2,1),r(3,2),r(2,3))-alpha(r(2,1),r(1,2),r(2,3)))/(2*hy)
          dza = (alpha(r(2,1),r(2,2),r(3,3))-alpha(r(2,1),r(2,2),r(1,3)))/(2*hz)

          !c2_t3 = (Grad(alpha) cross B)/|B|
          c2_t3_x(ix,iy,iz) = -1.0*(dya*bz(ix,iy,iz)-dza*by(ix,iy,iz))/magB(ix,iy,iz)
          c2_t3_y(ix,iy,iz) = -1.0*(dza*bx(ix,iy,iz)-dxa*bz(ix,iy,iz))/magB(ix,iy,iz)
          c2_t3_z(ix,iy,iz) = -1.0*(dxa*by(ix,iy,iz)-dya*bx(ix,iy,iz))/magB(ix,iy,iz)

        end do
      end do
    end do

    mag = sqrt(c2_t3_x**2 + c2_t3_y**2 + c2_t3_z**2)

    deallocate(magB)

    end subroutine compute_c2_t3

end module reconnectionMetrics
