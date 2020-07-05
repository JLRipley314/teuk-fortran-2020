module mod_bkgrd
!-----------------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env, only: ip => int64, rp => real64

    implicit none
!-----------------------------------------------------------------------------
    integer(ip) :: dims(2)

    real(rp), allocatable :: &
        r_pts(:), y_pts(:), &
        cs(:), sn(:) ! cosine and sine and gauss points

    complex(rp), allocatable :: &
        mu_0(:,:), ta_0(:,:), pi_0(:,:), rh_0(:,:), &
        thorn_prime_ta_0(:,:), &
        psi2_0(:,:)
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
    subroutine init_bkgrd(dims)
        integer(ip), intent(in) :: dims(2) ! domain size in x and y

        integer(ip) :: i, j

        j_loop: do j=1,dims(2)

            y_pts(j) = 0
            cs(j)    = cos(y_pts(j))
            sn(j)    = sin(y_pts(j))

            i_loop: do i=1,dims(2)
                r_pts(i) = 0

                mu_0(i,j) = 0
                ta_0(i,j) = 0
                pi_0(i,j) = 0
                rh_0(i,j) = 0

                thorn_prime_ta_0(i,j) = 0

                psi2_0(i,j) = 0

            end do i_loop
        end do j_loop

    end subroutine init_bkgrd
!-----------------------------------------------------------------------------
end module mod_bkgrd
