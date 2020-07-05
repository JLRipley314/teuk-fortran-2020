module mod_bkgrd
!-----------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: ip => int64, rp => real64

  use mod_sim_params, only: &
    nx, ny, cl, &
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
    complex(rp) :: val

    allocate(mu_0(nx,ny))
    allocate(ta_0(nx,ny))
    allocate(pi_0(nx,ny))
    allocate(rh_0(nx,ny))

    allocate(thorn_prime_ta_0(nx,ny))

    allocate(psi2_0(nx,ny))

    y_loop: do j=1,ny
      x_loop: do i=1,nx

        val = (2*cl,cl)

        mu_0(i,j) = (cl,cl) !( &
!        - ((cl**2) / ((cl**4) + (bhs**2)*(cs(j)**2)*(R(i)**2))), &
!        - ((bhs*cs(j)*R(i)) / ((cl**4) + (bhs**2)*(cs(j)**2)*(R(i)**2))) & 
!        )
        ta_0(i,j) = 0
        pi_0(i,j) = 0
        rh_0(i,j) = 0

        thorn_prime_ta_0(i,j) = 0

        psi2_0(i,j) = 0

      end do x_loop
    end do y_loop

  end subroutine init_bkgrd
!-----------------------------------------------------------------------------
end module mod_bkgrd
