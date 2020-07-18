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
      dt, nx, ny, min_m, max_m, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin

!=============================================================================
   implicit none
   private

   type(field), public :: &
      psi4, psi3, psi2,   &
      la, pi,             &
      muhll, hlmb, hmbmb  
!=============================================================================
   contains
!=============================================================================
   pure subroutine set_k(m_ang, fname, falloff, level, level_DR, kl)
      integer(ip),  intent(in) :: m_ang
      character(*), intent(in) :: fname
      integer(ip),  intent(in) :: falloff
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)    :: level
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: level_DR
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: kl

      real(rp) :: R

      integer(ip) :: i,j

      select_field: select case (fname)
         !-------------------------------------------------------------------
         case ("psi3")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -  R*4.0_rp*mu_0(i,j)*level(i,j,m_ang) &
               -  R*ta_0(i,j)*(psi4%level(i,j,m_ang)) &
               +    (psi4%edth(i,j,m_ang))
            end do 
            end do
         !--------------------------------------------------------------------
         case ("la")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -  R*(mu_0(i,j) + conjg(mu_0(i,j)))*level(i,j,m_ang) &
               -    (psi4%level(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("psi2")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -  R*(3.0_rp*mu_0(i,j))*level(i,j,m_ang) &
               -  R*2.0_rp*ta_0(i,j)*(psi3%level(i,j,m_ang)) &
               +    (psi3%edth(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("hmbmb")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = & 
                  R*(mu_0(i,j) - conjg(mu_0(i,j)))*level(i,j,m_ang) &
               -    2.0_rp*(la%level(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("pi")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -       R*(conjg(pi_0(i,j)) + ta_0(i,j))*(la%level(i,j,m_ang)) &
               +  (R**2)*0.5_rp*mu_0(i,j)*(conjg(pi_0(i,j))+ta_0(i,j))*(hmbmb%level(i,j,m_ang)) &
               -         (psi3%level(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("hlmb")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -  R*conjg(mu_0(i,j))*level(i,j,m_ang) &
               -    2.0_rp*(pi%level(i,j,m_ang)) &
               -  R*ta_0(i,j)*(hmbmb%level(i,j,m_ang)) 
            end do
            end do
         !--------------------------------------------------------------------
         case ("muhll")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -       R*conjg(mu_0(i,j))*level(i,j,m_ang) &
               -       R*mu_0(i,j)*(hlmb%edth(i,j,m_ang)) &
               -  (R**2)*mu_0(i,j)*(conjg(pi_0(i,j))+2.0_rp*ta_0(i,j))*(hlmb%level(i,j,m_ang)) &
               -         2.0_rp*(pi%edth(i,j,m_ang)) &
               -       R*2.0_rp*conjg(pi_0(i,j))*(pi%level(i,j,m_ang)) &
               -         2.0_rp*(psi2%level(i,j,m_ang)) &
               -  (R**2)*pi_0(i,j)*conjg(hmbmb%edth(i,j,-m_ang)) &
               +  (R**2)*(pi_0(i,j)**2)*conjg(hmbmb%level(i,j,-m_ang)) &
               +       R*mu_0(i,j)*conjg(hlmb%edth(i,j,-m_ang))  &
               -  (R**2)*3.0_rp*mu_0(i,j)*pi_0(i,j)*conjg(hlmb%level(i,j,-m_ang)) &
               +  (R**2)*2.0_rp*conjg(mu_0(i,j))*pi_0(i,j)*conjg(hlmb%level(i,j,-m_ang)) &
               -  (R**2)*2.0_rp*mu_0(i,j)*conjg(ta_0(i,j))*conjg(hlmb%level(i,j,-m_ang)) &
               -       R*2.0_rp*pi_0(i,j)*conjg(pi%level(i,j,-m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case default
            kl = -1.0_rp

      end select select_field
      !-----------------------------------------------------------------------
      call compute_DR(m_ang, level, level_DR)

      do j=1,ny
      do i=1,nx
         R = Rvec(i)

         kl(i,j,m_ang) = kl(i,j,m_ang) &
         -  ((R/cl)**2)*level_DR(i,j,m_ang) &
         -  (falloff*R/(cl**2))*level(i,j,m_ang)
      
         kl(i,j,m_ang) = kl(i,j,m_ang) / (2.0_rp*(1.0_rp+(2.0_rp*bhm*R/(cl**2))))
      end do
      end do

   end subroutine set_k
!=============================================================================
   pure subroutine take_step(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      !-----------------------------------------------------------------------
      select_field: select case (step)
         !--------------------------------------------------------------------
         case (1)
            call set_k(m_ang, f%name, f%falloff, f%n,   f%DR, f%k1)
         !--------------------------------------------------------------------
         case (2)
            call set_k(m_ang, f%name, f%falloff, f%l2,  f%DR, f%k2)
         !--------------------------------------------------------------------
         case (3)
            call set_k(m_ang, f%name, f%falloff, f%l3,  f%DR, f%k3)
         !--------------------------------------------------------------------
         case (4)
            call set_k(m_ang, f%name, f%falloff, f%l4,  f%DR, f%k4)
         !--------------------------------------------------------------------
         case (5)
            call set_k(m_ang, f%name, f%falloff, f%np1, f%DR, f%k5)
         !--------------------------------------------------------------------
         case default
            f%error = -1.0_ip
      end select select_field
      !-----------------------------------------------------------------------
   end subroutine take_step
!=============================================================================
   pure subroutine step_all_fields(step, m_ang)
      integer(ip), intent(in) :: step, m_ang
   end subroutine step_all_fields
!=============================================================================
end module metric_recon 
