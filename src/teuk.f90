module teuk
    use, intrinsic :: iso_fortran_env, only: rp=>real64, ip=>int64
    use mod_field 
    implicit none
!-----------------------------------------------------------------------------
    real(rp), allocatable :: &
        A_pp(:,:), A_pq(:,:), A_pf(:,:), &
        A_qp(:,:), A_qq(:,:), A_qf(:,:), &
        A_fp(:,:), A_fq(:,:), A_ff(:,:), &
        B_pp(:,:), B_pq(:,:), B_pf(:,:), &
        B_qp(:,:), B_qq(:,:), B_qf(:,:), &
        B_fp(:,:), B_fq(:,:), B_ff(:,:)
!-----------------------------------------------------------------------------
    contains
!-----------------------------------------------------------------------------
    subroutine init_teuk(nx,ny)
        integer(ip), intent(in) :: nx, ny

        integer(ip) :: i, j

        j_loop: do j=1,ny
            i_loop: do i=1,nx

                A_pp(i,j) = 1
                A_pq(i,j) = 1

            end do i_loop
        end do j_loop

    end subroutine init_teuk
!-----------------------------------------------------------------------------
end module teuk

