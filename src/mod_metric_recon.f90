!
! Metric reconstruction evolution equations
!
!=============================================================================
module metric_recon
!=============================================================================
   use mod_prec
   use mod_cheb,     only: Rvec=>R, compute_DR
   use mod_field,    only: field, set_level, set_DT
   use mod_ghp,      only: ghp_edth, ghp_edth_prime, ghp_thorn, ghp_thorn_prime
   use mod_bkgrd_np, only: mu_0, ta_0, pi_0, rh_0, thorn_prime_ta_0, psi2_0
   use mod_params,   only: &
      dt, nx, ny, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin

!=============================================================================
   implicit none
   private

   public :: set_k 

   complex(rp), dimension(nx,ny) :: &
      psi4,  psi3,    psi2, &
      la,      pi,          &
      muhll, hlmb,   hmbmb, &
      edth_psi4, edth_psi3, &
      edth_pi,              &
      edth_hlmb, edth_hmbmb
!=============================================================================
   contains
!=============================================================================
   pure subroutine set_k(fname, kl)
      character(*), intent(in) :: fname
      complex(rp), dimension(nx,ny), intent(inout) :: kl

      integer(ip) :: i,j

      select_field: select case (fname)
         !-------------------------------------------------------------------
         case ("psi3")
            do j=1,ny
            do i=1,nx
               kl(i,j) = &
               -  4.0_rp*mu_0(i,j)*psi3(i,j) &
               -  ta_0(i,j)*psi4(i,j) &
               +  edth_psi4(i,j)
            end do 
            end do
         !--------------------------------------------------------------------
         case ("la")
            do j=1,ny
            do i=1,nx
               kl(i,j) = &
               -  (mu_0(i,j) + conjg(mu_0(i,j))) * la(i,j) &
               -  psi4(i,j)
            end do
            end do
         !--------------------------------------------------------------------
         case ("psi2")
            do j=1,ny
            do i=1,nx
               kl(i,j) = &
               -  (3.0_rp*mu_0(i,j)) * psi2(i,j) &
               -  2.0_rp*ta_0(i,j)*psi3(i,j) &
               +  edth_psi3(i,j)
            end do
            end do
         !--------------------------------------------------------------------
         case ("hmbmb")
            do j=1,ny
            do i=1,nx
               kl(i,j) = & 
                  (mu_0(i,j) - conjg(mu_0(i,j))) * hmbmb(i,j) &
               -  2.0_rp*la(i,j)
            end do
            end do
         !--------------------------------------------------------------------
         case ("pi")
            do j=1,ny
            do i=1,nx
               kl(i,j) = &
               -  (conjg(pi_0(i,j)) + ta_0(i,j)) * la(i,j) &
               +  0.5_rp*mu_0(i,j)*(conjg(pi_0(i,j))+ta_0(i,j))*hmbmb(i,j) &
               -  psi3(i,j)
            end do
            end do
         !--------------------------------------------------------------------
         case ("hlmb")
            do j=1,ny
            do i=1,nx
               kl(i,j) = &
               -  2.0_rp*pi(i,j) &
               -  ta_0(i,j)*hmbmb(i,j) &
               -  conjg(mu_0(i,j))*hlmb(i,j)
            end do
            end do
         !--------------------------------------------------------------------
         case ("muhll")
            do j=1,ny
            do i=1,nx
               kl(i,j) = &
               -  conjg(mu_0(i,j))*muhll(i,j) &
               -  mu_0(i,j)*edth_hlmb(i,j) &
               -  mu_0(i,j)*(conjg(pi_0(i,j))+2.0_rp*ta_0(i,j))*hlmb(i,j) &
               -  2.0_rp*edth_pi(i,j) &
               -  2.0_rp*conjg(pi_0(i,j))*pi(i,j) &
               -  2.0_rp*psi2(i,j) &
               -  pi_0(i,j)*conjg(edth_hmbmb(i,j)) &
               +  (pi_0(i,j)**2)*conjg(hmbmb(i,j)) &
               +  mu_0(i,j)*conjg(edth_hlmb(i,j))  &
               -  3.0_rp*mu_0(i,j)*pi_0(i,j)*conjg(hlmb(i,j)) &
               +  2.0_rp*conjg(mu_0(i,j))*pi_0(i,j)*conjg(hlmb(i,j)) &
               -  2.0_rp*mu_0(i,j)*conjg(ta_0(i,j))*conjg(hlmb(i,j)) &
               -  2.0_rp*pi_0(i,j)*conjg(pi(i,j))
            end do
            end do
         !--------------------------------------------------------------------
         case default
            kl = -1.0_rp

      end select select_field

   end subroutine set_k
!=============================================================================
end module metric_recon 
