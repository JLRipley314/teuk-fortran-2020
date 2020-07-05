module mod_bkgrd
!-----------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: ip => int64, rp => real64

  use mod_sim_params, only: &
    nx, ny, cl, &
    bhs, bhm, &
    r_pts, y_pts, cs, sn

  implicit none
! all variables are public
!-----------------------------------------------------------------------------
  integer(ip) :: dims(2)

  complex(rp), allocatable :: &
    mu_0(:,:), ta_0(:,:), pi_0(:,:), rh_0(:,:), &
    thorn_prime_ta_0(:,:), &
    psi2_0(:,:)
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
  subroutine init_bkgrd()

    integer(ip) :: i, j

    allocate(mu_0(nx,ny))
    allocate(ta_0(nx,ny))
    allocate(pi_0(nx,ny))
    allocate(rh_0(nx,ny))

    allocate(thorn_prime_ta_0(nx,ny))

    allocate(psi2_0(nx,ny))

    y_loop: do j=1,ny
      x_loop: do i=1,nx

        mu_0(i,j) = cmplx( &
        - (cl**2) &
          / ((cl**4) + (bhs**2)*(cs(j)**2)*(r_pts(i)**2)) &
          , &
        - (bhs*cs(j)*r_pts(i)) & 
          / ((cl**4) + (bhs**2)*(cs(j)**2)*(r_pts(i)**2)) &
          , & 
          kind=rp)

        ta_0(i,j) = cmplx( &
        - (sqrt(2.0_rp)*(bhs**2)*cs(j)*(cl**2)*r_pts(i)*sn(j)) &
          / (((cl**4) + (bhs**2)*(cs(j)**2)*(r_pts(i)**2))**2) &
          , &
          (bhs*sn(j)*(4.0_rp*(cl**4) + (bhs**2)*(r_pts(i)**2)*(-1.0_rp - 3.0_rp*(cs(j)**2) + (sn(j)**2)))) &
          / (sqrt(2.0_rp)*((2*(cl**4) + (bhs**2)*(r_pts(i)**2)*(1.0_rp + (cs(j)**2) - (sn(j)**2)))**2)) &
          ,&
          kind=rp)

        pi_0(i,j) = cmplx( &
          0.0_rp &
          , &
        - (bhs*sn(j)) &
          / (sqrt(2.0_rp)*(cl**4) + (bhs**2)*(cs(j)**2)*(r_pts(i)**2)) &
          , &
          kind=rp)

        rh_0(i,j) = cmplx( &
          (-2.0_rp*((cl**6) - 2.0_rp*(cl**4)*bhm*r_pts(i) + (bhs**2)*(cl**2)*(r_pts(i)**2))) &
          / ((2.0_rp*(cl**4) + (bhs**2)*(r_pts(i)**2)*(1.0_rp + (cs(j)**2) - (sn(j)**2)))**2) &
          , & 
          (-2.0_rp*bhs*cs(j)*r_pts(i)*((cl**4) - 2.0_rp*(cl**2)*bhm*r_pts(i) + (bhs**2)*(r_pts(i)**2))) &
          / ((2.0_rp*(cl**4) + (bhs**2)*(r_pts(i)**2)*(1.0_rp + (cs(j)**2) - (sn(j)**2)))**2) &
          , &
          kind=rp)

        thorn_prime_ta_0(i,j) = 0

        psi2_0(i,j) = 0

      end do x_loop
    end do y_loop

  end subroutine init_bkgrd
!-----------------------------------------------------------------------------
end module mod_bkgrd
