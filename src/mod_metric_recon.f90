!
! Metric reconstruction evolution equations
!
!=============================================================================
module mod_metric_recon
!=============================================================================
   use mod_prec

   use mod_cheb,     only: R=>Rarr, compute_DR
   use mod_field,    only: field, set_level
   use mod_ghp,      only: set_edth, set_edth_prime, set_thorn, set_thorn_prime
   use mod_bkgrd_np, only: mu_0, ta_0, pi_0, rh_0, psi2_0
   use mod_params,   only: &
      dt, nx, ny, min_m, max_m, &
      cl=>compactification_length, &
      bhm=>black_hole_mass

   use mod_fields_list, only: &
      psi4_lin_f, &
      psi3, psi2,         &
      la, pi,             &
      muhll, hlmb, hmbmb, &
      res_bianchi3, res_bianchi2, res_hll 
!=============================================================================
   implicit none
   private

   public :: metric_recon_time_step, metric_recon_indep_res
!=============================================================================
   contains
!=============================================================================
   subroutine set_k(step, m_ang, f, kl)
      integer(ip),  intent(in) :: step, m_ang
      type(field), intent(inout) :: f
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: kl
      !----------------------------------------------------------------------
      call set_level(step, m_ang, f)
      !----------------------------------------------------------------------
      select_field: select case (f%fname)
         !-------------------------------------------------------------------
         case ("psi3")
            kl(:,:,m_ang) = &
            -  4.0_rp*R*mu_0*f%level(:,:,m_ang) &

            -  R*ta_0*(psi4_lin_f%level(:,:,m_ang)) &

            +  (psi4_lin_f%edth(:,:,m_ang))
         !--------------------------------------------------------------------
         case ("la")
            kl(:,:,m_ang) = &
            -  R*( &
                  mu_0 &
               +  conjg(mu_0) &
               )*f%level(:,:,m_ang) &

            - (psi4_lin_f%level(:,:,m_ang))
         !--------------------------------------------------------------------
         case ("psi2")
            kl(:,:,m_ang) = &
            -  3.0_rp*R*mu_0*f%level(:,:,m_ang) &

            -  2.0_rp*R*ta_0*(psi3%level(:,:,m_ang)) &

            +  (psi3%edth(:,:,m_ang))
         !--------------------------------------------------------------------
         case ("hmbmb")
            kl(:,:,m_ang) = & 
               R*( &
                  mu_0 &
               -  conjg(mu_0) &
               )*f%level(:,:,m_ang) &

            -  2.0_rp*(la%level(:,:,m_ang))
         !--------------------------------------------------------------------
         case ("pi")
            kl(:,:,m_ang) = &
            -  R*( &
                  conjg(pi_0) &
               +  ta_0 &
               )*(la%level(:,:,m_ang)) &

            +  (R**2)*0.5_rp*mu_0*( &
                  conjg(pi_0) &
               +  ta_0 &
               )*(hmbmb%level(:,:,m_ang)) &
            -  (psi3%level(:,:,m_ang))
         !--------------------------------------------------------------------
         case ("hlmb")
               kl(:,:,m_ang) = &
               -  R*conjg(mu_0)*f%level(:,:,m_ang) &

               -  2.0_rp*(pi%level(:,:,m_ang)) &

               -  R*ta_0*(hmbmb%level(:,:,m_ang)) 
         !--------------------------------------------------------------------
         case ("muhll")
            kl(:,:,m_ang) = &
            -  R*conjg(mu_0)*f%level(:,:,m_ang) &

            -  R*mu_0*(hlmb%edth(:,:,m_ang)) &

            -  (R**2)*mu_0*( &
                  conjg(pi_0) &
               +  2.0_rp*ta_0 &
               )*(hlmb%level(:,:,m_ang)) &

            -  2.0_rp*(pi%edth(:,:,m_ang)) &

            -  2.0_rp*R*conjg(pi_0)*(pi%level(:,:,m_ang)) &

            -  2.0_rp*(psi2%level(:,:,m_ang)) &
            !
            ! complex conjugate of fields
            !
            -  2.0_rp*R*pi_0*conjg(pi%level(:,:,-m_ang)) &

            -  R*pi_0*conjg(hmbmb%edth(:,:,-m_ang)) &

            +  (R**2)*(pi_0**2)*conjg(hmbmb%level(:,:,-m_ang)) &

            +  R*mu_0*conjg(hlmb%edth(:,:,-m_ang))  &

            +  (R**2)*( &
               -  3.0_rp*mu_0*pi_0 &
               +  2.0_rp*conjg(mu_0)*pi_0 &
               -  2.0_rp*mu_0*conjg(ta_0) &
               )*conjg(hlmb%level(:,:,-m_ang)) 
         !--------------------------------------------------------------------
         case default
            write (*,*) "ERROR: set_k, " // f%fname // ", not in list"
            stop
      end select select_field
      !-----------------------------------------------------------------------
      call compute_DR(step, m_ang, f)

      kl(:,:,m_ang) = ( &
            kl(:,:,m_ang) &

         -  ((R/cl)**2)*f%DR(:,:,m_ang) &

         -  (f%falloff*R/(cl**2))*f%level(:,:,m_ang) &
         )/( &
            2.0_rp+(4.0_rp*bhm*R/(cl**2)) &
         )

   end subroutine set_k
!=============================================================================
! RK4 evolution
!=============================================================================
   subroutine take_step(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f 
      !-----------------------------------------------------------------------
      select_step: select case (step)
         !--------------------------------------------------------------------
         case (1)
            !--------------------------------------------------------
            ! if first time then k1 has not been set from k5
            !--------------------------------------------------------
            if (f%first_time) then
               call set_k(step, m_ang, f, f%k1)

               f%first_time = .false.
            end if

            f%l2(:,:,m_ang)= f%n(:,:,m_ang)+0.5_rp*dt*f%k1(:,:,m_ang)
         !--------------------------------------------------------------------
         case (2)
            call set_k(step, m_ang, f, f%k2)

            f%l3(:,:,m_ang)= f%n(:,:,m_ang)+0.5_rp*dt*f%k2(:,:,m_ang)
         !--------------------------------------------------------------------
         case (3)
            call set_k(step, m_ang, f, f%k3)

            f%l4(:,:,m_ang)= f%n(:,:,m_ang)+dt*f%k3(:,:,m_ang)
         !--------------------------------------------------------------------
         case (4)
            call set_k(step, m_ang, f, f%k4)

            f%np1(:,:,m_ang)= f%n(:,:,m_ang) &
            +  (dt/6.0_rp)*(f%k1(:,:,m_ang)+2.0_rp*f%k2(:,:,m_ang)+2.0_rp*f%k3(:,:,m_ang)+f%k4(:,:,m_ang))
         !--------------------------------------------------------------------
         case (5)
            !-----------------------------------------------------------------
            ! need k5 for source term and independent residual
            !-----------------------------------------------------------------
            call set_k(step, m_ang, f, f%k5)
         !--------------------------------------------------------------------
         case default
            write(*,*) "Error: take_step, " // f%fname // ", step != 1,..,5"
            stop
      end select select_step
      !-----------------------------------------------------------------------
   end subroutine take_step
!=============================================================================
   subroutine set_indep_res(m_ang, fname)
      integer(ip),  intent(in) :: m_ang
      character(*), intent(in) :: fname
      !-----------------------------------------------------------------------
      select_field: select case (fname)
         !--------------------------------------------------------------------
         case ("res_bianchi3")

            call set_level(5_ip,m_ang,psi4_lin_f)
            call set_level(5_ip,m_ang,psi3)
            call set_level(5_ip,m_ang,la)

            call set_thorn(     5_ip,m_ang,psi4_lin_f)
            call set_edth_prime(5_ip,m_ang,psi3)
            !-----------------------------------------------   
            res_bianchi3%np1(:,:,m_ang) = &
               R*psi3%edth_prime(:,:,m_ang) &

            +  4.0_rp*(R**2)*pi_0*psi3%level(:,:,m_ang) &

            -  psi4_lin_f%thorn(:,:,m_ang) &

            +  rh_0*psi4_lin_f%level(:,:,m_ang) &

            -  3.0_rp*(R**2)*psi2_0*la%level(:,:,m_ang) 
         !--------------------------------------------------------------------
         case ("res_bianchi2")

            call set_level(5_ip,m_ang,pi)

            call set_level(5_ip,m_ang,psi2)
            call set_level(5_ip,m_ang,psi3)

            call set_level(5_ip,m_ang,hlmb)
            call set_level(5_ip,m_ang,hmbmb)

            call set_thorn(     5_ip,m_ang,psi3)
            call set_edth_prime(5_ip,m_ang,psi2)
            !-----------------------------------------------   
            res_bianchi2%np1(:,:,m_ang) = &
               psi2_0*( &
               -  3.0_rp*(R**3)*mu_0*hlmb%level(:,:,m_ang) &
               -  1.5_rp*(R**3)*ta_0*hmbmb%level(:,:,m_ang) &
               -       3.0_rp*(R**2)*pi%level(:,:,m_ang) &
               ) &

            -  R*psi2%edth_prime(:,:,m_ang) &

            -  3.0_rp*(R**2)*pi_0*psi2%level(:,:,m_ang) &

            +  psi3%thorn(:,:,m_ang) &

            -  2.0_rp*rh_0*psi3%level(:,:,m_ang)
         !--------------------------------------------------------------------
         ! hll should be real
         !--------------------------------------------------------------------
         case ("res_hll")

            call set_level(5_ip,-m_ang,muhll)
            call set_level(5_ip, m_ang,muhll)
            !-----------------------------------------------   
            res_hll%np1(:,:,m_ang) = &
                    (muhll%level(:,:, m_ang) / mu_0) &
            -  conjg(muhll%level(:,:,-m_ang) / mu_0) 
         !-------------------------------------------------------------------- 
         case default
            continue
      end select select_field
   end subroutine set_indep_res
!=============================================================================
   subroutine step_all_fields_preserve_m_ang(step, m_ang)
      integer(ip), intent(in) :: step, m_ang
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,psi4_lin_f)
      call set_edth( step,m_ang,psi4_lin_f)

      call take_step(step,m_ang,psi3)
      call take_step(step,m_ang,la)
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,psi3)
      call set_edth( step,m_ang,psi3)

      call take_step(step,m_ang,psi2)
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,la)

      call take_step(step,m_ang,hmbmb)
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,hmbmb)

      call take_step(step,m_ang,pi)
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,pi)

      call take_step(step,m_ang,hlmb)
   end subroutine step_all_fields_preserve_m_ang
!=============================================================================
! evolve +/- m_ang as muhll requires both
!=============================================================================
   subroutine step_all_fields_mix_m_ang(step, m_ang)
      integer(ip), intent(in) :: step, m_ang
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,   pi); call set_level(step,-m_ang,   pi)
      call set_level(step,m_ang, psi2); call set_level(step,-m_ang, psi2)
      call set_level(step,m_ang, hlmb); call set_level(step,-m_ang, hlmb)
      call set_level(step,m_ang,hmbmb); call set_level(step,-m_ang,hmbmb)

      call set_edth(step,m_ang,   pi); call set_edth(step,-m_ang,   pi)
      call set_edth(step,m_ang, hlmb); call set_edth(step,-m_ang, hlmb)
      call set_edth(step,m_ang,hmbmb); call set_edth(step,-m_ang,hmbmb)

      call take_step(step,m_ang,muhll); call take_step(step,-m_ang,muhll)
   end subroutine step_all_fields_mix_m_ang
!=============================================================================
! evolve +/- m angular number fields
!=============================================================================
   subroutine metric_recon_time_step(m_ang)
      integer(ip), intent(in) :: m_ang

      call step_all_fields_preserve_m_ang(1_ip,m_ang);
      call step_all_fields_preserve_m_ang(2_ip,m_ang);
      call step_all_fields_preserve_m_ang(3_ip,m_ang);
      call step_all_fields_preserve_m_ang(4_ip,m_ang);
      call step_all_fields_preserve_m_ang(5_ip,m_ang);

      if (m_ang>0) then
         call step_all_fields_preserve_m_ang(1_ip,-m_ang);
         call step_all_fields_preserve_m_ang(2_ip,-m_ang);
         call step_all_fields_preserve_m_ang(3_ip,-m_ang);
         call step_all_fields_preserve_m_ang(4_ip,-m_ang);
         call step_all_fields_preserve_m_ang(5_ip,-m_ang);
      end if

      call step_all_fields_mix_m_ang(1_ip,m_ang);
      call step_all_fields_mix_m_ang(2_ip,m_ang);
      call step_all_fields_mix_m_ang(3_ip,m_ang);
      call step_all_fields_mix_m_ang(4_ip,m_ang);
      call step_all_fields_mix_m_ang(5_ip,m_ang);

   end subroutine metric_recon_time_step
!=============================================================================
   subroutine metric_recon_indep_res(m_ang)
      integer(ip), intent(in) :: m_ang

      call set_indep_res(m_ang,"res_bianchi3")
      call set_indep_res(m_ang,"res_bianchi2")
      call set_indep_res(m_ang,"res_hll")

   end subroutine metric_recon_indep_res
!=============================================================================
end module mod_metric_recon 
