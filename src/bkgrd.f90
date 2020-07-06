module mod_bkgrd_np
!-----------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: ip => int64, rp => real64

  use mod_sim_params, only: &
    nx, ny, cl, &
    bhs, bhm, &
    R, cy, sy

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
  subroutine init_bkgrd_np()

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
          / ((cl**4) + (bhs**2)*(cy(j)**2)*(R(i)**2)) &
          , &
        - (bhs*cy(j)*R(i)) & 
          / ((cl**4) + (bhs**2)*(cy(j)**2)*(R(i)**2)) &
          , & 
          kind=rp)

        ta_0(i,j) = cmplx( &
        - (sqrt(2.0_rp)*(bhs**2)*cy(j)*(cl**2)*R(i)*sy(j)) &
          / (((cl**4) + (bhs**2)*(cy(j)**2)*(R(i)**2))**2) &
          , &
          (bhs*sy(j)*(4.0_rp*(cl**4) + (bhs**2)*(R(i)**2)*(-1.0_rp - 3.0_rp*(cy(j)**2) + (sy(j)**2)))) &
          / (sqrt(2.0_rp)*((2*(cl**4) + (bhs**2)*(R(i)**2)*(1.0_rp + (cy(j)**2) - (sy(j)**2)))**2)) &
          ,&
          kind=rp)

        pi_0(i,j) = cmplx( &
          0.0_rp &
          , &
        - (bhs*sy(j)) &
          / (sqrt(2.0_rp)*(cl**4) + (bhs**2)*(cy(j)**2)*(R(i)**2)) &
          , &
          kind=rp)

        rh_0(i,j) = cmplx( &
          (-2.0_rp*((cl**6) - 2.0_rp*(cl**4)*bhm*R(i) + (bhs**2)*(cl**2)*(R(i)**2))) &
          / ((2.0_rp*(cl**4) + (bhs**2)*(R(i)**2)*(1.0_rp + (cy(j)**2) - (sy(j)**2)))**2) &
          , & 
          (-2.0_rp*bhs*cy(j)*R(i)*((cl**4) - 2.0_rp*(cl**2)*bhm*R(i) + (bhs**2)*(R(i)**2))) &
          / ((2.0_rp*(cl**4) + (bhs**2)*(R(i)**2)*(1.0_rp + (cy(j)**2) - (sy(j)**2)))**2) &
          , &
          kind=rp)

        thorn_prime_ta_0(i,j) = cmplx( &
          (4.0_rp*sqrt(2.0_rp)*(bhs**2)*cy(j)*R(i)*sy(j)*(-6.0_rp*(cl**4) &
            + (bhs**2)*(R(i)**2)*(1.0_rp + (cy(j)**2) - (sy(j)**2))) &
          ) &
          / ((2.0_rp*(cl**4) + (bhs**2)*(R(i)**2)*(1.0_rp + (cy(j)**2) - (sy(j)**2)))**3) &
          , &
          (2.0_rp*sqrt(2.0_rp)*bhs*(cl**2)*sy(j)*( &
              4.0_rp*(cl**4) &
            + 3.0_rp*(bhs**2)*(R(i)**2)*(-1.0_rp - 3.0_rp*(cy(j)**2) + (sy(j)**2)) &
            )) &
          / ((2.0_rp*(cl**4) + (bhs**2)*(R(i)**2)*(1.0_rp + (cy(j)**2) - (sy(j)**2)))**3) &
          , &
          kind=rp)

        psi2_0(i,j) = cmplx( &
          (4.0_rp*(cl**2)*bhm*( &
            - 2.0_rp*(cl**4) &
            + 3.0_rp*(bhs**2)*(R(i)**2)*(1.0_rp + (cy(j)**2) - (sy(j)**2)) &
            )) &
          / ((2.0_rp*(cl**4) + (bhs**2)*(R(i)**2)*(1.0_rp + (cy(j)**2) - (sy(j)**2)))**3) &
          , &
          (2.0_rp*bhs*cy(j)*bhm*R(i)*( &
            - 12.0_rp*(cl**4) &
            + (bhs**2)*(R(i)**2)*(3.0_rp + (cy(j)**2) - 3.0_rp*(sy(j)**2)) &
            )) &
          / ((2.0_rp*(cl**4) + (bhs**2)*(R(i)**2)*(1.0_rp + (cy(j)**2) - (sy(j)**2)))**3) &
          , &
          kind=rp)

      end do x_loop
    end do y_loop

  end subroutine init_bkgrd_np

  ! allocatable arrays are freed at end of main
  ! but clean is here so valgrind doesn't get confused
  subroutine clear_bkgrd_np

    deallocate(mu_0)
    deallocate(ta_0)
    deallocate(pi_0)
    deallocate(rh_0)

    deallocate(thorn_prime_ta_0)

    deallocate(psi2_0)
  end subroutine clear_bkgrd_np
!-----------------------------------------------------------------------------
end module mod_bkgrd_np
