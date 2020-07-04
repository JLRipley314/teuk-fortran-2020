module teuk
    use, intrinsic :: iso_fortran_env, only: prec=>real64
    implicit none

    real(prec), allocatable :: &
        A_pp(:,:), A_pq(:,:), A_pf(:,:), &
        A_qp(:,:), A_qq(:,:), A_qf(:,:), &
        A_fp(:,:), A_fq(:,:), A_ff(:,:), &
        B_pp(:,:), B_pq(:,:), B_pf(:,:), &
        B_qp(:,:), B_qq(:,:), B_qf(:,:), &
        B_fp(:,:), B_fq(:,:), B_ff(:,:)

end module teuk

