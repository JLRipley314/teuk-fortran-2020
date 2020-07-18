module mod_ghp
!=============================================================================
   use mod_prec
   use mod_cheb,   only: Rvec=>R, compute_DR
   use mod_swal,   only: Yvec=>Y, cy, sy, swal_lower, swal_raise
   use mod_bkgrd,  only: ep_0
   use mod_params, only: &
      dt, nx, ny, lmax, &
      spin, min_m, max_m, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin
!=============================================================================
   implicit none
   private

   public :: ghp_edth
!=============================================================================
   pure subroutine ghp_edth(level, f, f_edth)
      integer(ip), intent(in)    :: level
      type(field), intent(inout) :: f
      complex(rp), intent(out)   :: f_edth

      call compute_DR(f

   end subroutine ghp_edth
!=============================================================================
end module mod_ghp
