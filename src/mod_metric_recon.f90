!
! Metric reconstruction evolution equations
!
!=============================================================================
module mod_metric_recon
!=============================================================================
   use mod_prec
   use mod_cheb,     only: R, compute_DR
   use mod_field,    only: field, set_field, set_level
   use mod_ghp,      only: set_edth, set_edth_prime, set_thorn, set_thorn_prime
   use mod_bkgrd_np, only: mu_0, ta_0, pi_0, rh_0, psi2_0
   use mod_teuk,     only: psi4_f
   use mod_params,   only: &
      dt, nx, ny, min_m, max_m, &
      cl=>compactification_length, &
      bhm=>black_hole_mass

!=============================================================================
   implicit none
   private

   public :: metric_recon_time_step, metric_recon_indep_res

   type(field), public :: &
      psi3, psi2,         &
      la, pi,             &
      muhll, hlmb, hmbmb, &
      bianchi3_res, bianchi2_res, hll_res 
!=============================================================================
   contains
!=============================================================================
   pure subroutine set_k(error, m_ang, fname, falloff, level, level_DR, kl)
      character(*), intent(inout) :: error
      integer(ip),  intent(in) :: m_ang
      character(*), intent(in) :: fname
      integer(ip),  intent(in) :: falloff
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: level
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: level_DR
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: kl

      integer(ip) :: i,j

      select_field: select case (fname)
         !-------------------------------------------------------------------
         case ("psi3")
            do j=1,ny
            do i=1,nx
               kl(i,j,m_ang) = &
               -  4.0_rp*R(i)*mu_0(i,j)*level(i,j,m_ang) &
               -         R(i)*ta_0(i,j)*(psi4_f%level(i,j,m_ang)) &
               +                        (psi4_f%edth(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("la")
            do j=1,ny
            do i=1,nx
               kl(i,j,m_ang) = &
               -  R(i)*(mu_0(i,j)+conjg(mu_0(i,j)))*level(i,j,m_ang) &
               -                                    (psi4_f%level(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("psi2")
            do j=1,ny
            do i=1,nx
               kl(i,j,m_ang) = &
               -  3.0_rp*R(i)*mu_0(i,j)*level(i,j,m_ang) &
               -  2.0_rp*R(i)*ta_0(i,j)*(psi3%level(i,j,m_ang)) &
               +                        (psi3%edth(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("hmbmb")
            do j=1,ny
            do i=1,nx
               kl(i,j,m_ang) = & 
                  R(i)*(mu_0(i,j)-conjg(mu_0(i,j)))*level(i,j,m_ang) &
               -                             2.0_rp*(la%level(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("pi")
            do j=1,ny
            do i=1,nx
               kl(i,j,m_ang) = &
               -  R(i)     *(conjg(pi_0(i,j))+ta_0(i,j))*(la%level(i,j,m_ang)) &
               +  (R(i)**2)*(0.5_rp*mu_0(i,j)*(conjg(pi_0(i,j))+ta_0(i,j))) &
                                                        *(hmbmb%level(i,j,m_ang)) &
               -                                         (psi3%level(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("hlmb")
            do j=1,ny
            do i=1,nx
               kl(i,j,m_ang) = &
               -  R(i)*conjg(mu_0(i,j))*level(i,j,m_ang) &
               -                 2.0_rp*(pi%level(i,j,m_ang)) &
               -         R(i)*ta_0(i,j)*(hmbmb%level(i,j,m_ang)) 
            end do
            end do
         !--------------------------------------------------------------------
         case ("muhll")
            do j=1,ny
            do i=1,nx
               kl(i,j,m_ang) = &
               -       R(i)*conjg(mu_0(i,j))*level(i,j,m_ang) &
               -              R(i)*mu_0(i,j)*(hlmb%edth(i,j,m_ang)) &
               -  (R(i)**2)*mu_0(i,j)*(conjg(pi_0(i,j))+2.0_rp*ta_0(i,j)) &
                                            *(hlmb%level(i,j,m_ang)) &
               -                      2.0_rp*(psi2%level(i,j,m_ang)) &
               -            R(i)*(pi_0(i,j))*conjg(hmbmb%edth( i,j,-m_ang)) &
               +    (R(i)**2)*(pi_0(i,j)**2)*conjg(hmbmb%level(i,j,-m_ang)) &
               +              R(i)*mu_0(i,j)*conjg(hlmb%edth(i,j,-m_ang))  &
               -  (R(i)**2)*3.0_rp*mu_0(i,j)*pi_0(i,j) &
                                            *conjg(hlmb%level(i,j,-m_ang)) &
               +  (R(i)**2)*2.0_rp*conjg(mu_0(i,j))*pi_0(i,j) &
                                            *conjg(hlmb%level(i,j,-m_ang)) &
               -  (R(i)**2)*2.0_rp*mu_0(i,j)*conjg(ta_0(i,j)) &
                                            *conjg(hlmb%level(i,j,-m_ang)) &
               -                      2.0_rp*(pi%edth(i,j,m_ang)) &
               -       R(i)*2.0_rp*conjg(pi_0(i,j)) &
                                            *(pi%level(i,j,m_ang)) &
               -       R(i)*2.0_rp*pi_0(i,j)*conjg(pi%level(i,j,-m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case default
            write(error,*) "ERROR: set_k, " // fname // ", not in list"

      end select select_field
      !-----------------------------------------------------------------------
      call compute_DR(m_ang, level, level_DR)

      do i=1,nx
         kl(i,:,m_ang) = kl(i,:,m_ang) &
         -          ((R(i)/cl)**2)*level_DR(i,:,m_ang) &
         -  (falloff*R(i)/(cl**2))*level(i,:,m_ang)
      
         kl(i,:,m_ang) = kl(i,:,m_ang) / (2.0_rp+(4.0_rp*bhm*R(i)/(cl**2)))
      end do

   end subroutine set_k
!=============================================================================
! RK4 evolution
!=============================================================================
   pure subroutine take_step(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f
      
      integer(ip) :: i, j

      !-----------------------------------------------------------------------
      select_step: select case (step)
         !--------------------------------------------------------------------
         case (1)
            !--------------------------------------------------------
            ! if first time then k1 has not been set from k5
            !--------------------------------------------------------
            if (f%first_time) then
               call set_k(f%error, m_ang, f%name, f%falloff, f%n, f%DR, f%k1)

               f%first_time = .false.
            end if

            do j=1,ny
            do i=1,nx
               f%l2(i,j,m_ang)= f%n(i,j,m_ang)+0.5_rp*dt*f%k1(i,j,m_ang)
            end do
            end do
         !--------------------------------------------------------------------
         case (2)
            call set_k(f%error, m_ang, f%name, f%falloff, f%l2, f%DR, f%k2)

            do j=1,ny
            do i=1,nx
               f%l3(i,j,m_ang)= f%n(i,j,m_ang)+0.5_rp*dt*f%k2(i,j,m_ang)
            end do
            end do
         !--------------------------------------------------------------------
         case (3)
            call set_k(f%error, m_ang, f%name, f%falloff, f%l3, f%DR, f%k3)

            do j=1,ny
            do i=1,nx
               f%l4(i,j,m_ang)= f%n(i,j,m_ang)+dt*f%k3(i,j,m_ang)
            end do
            end do
         !--------------------------------------------------------------------
         case (4)
            call set_k(f%error, m_ang, f%name, f%falloff, f%l4, f%DR, f%k4)

            do j=1,ny
            do i=1,nx
               f%np1(i,j,m_ang)= f%n(i,j,m_ang) &
               +  (dt/6.0_rp)*(f%k1(i,j,m_ang)+2.0_rp*f%k2(i,j,m_ang)+2.0_rp*f%k3(i,j,m_ang)+f%k4(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case (5)
            !-----------------------------------------------------------------
            ! need k5 for source term and independent residual
            !-----------------------------------------------------------------
            call set_k(f%error, m_ang, f%name, f%falloff, f%np1, f%DR, f%k5)
         !--------------------------------------------------------------------
         case default
            write(f%error,*) "Error: take_step, " // f%name // ", step != 1,..,5"
      end select select_step
      !-----------------------------------------------------------------------
   end subroutine take_step
!=============================================================================
   subroutine set_indep_res(m_ang, fname)
      integer(ip),  intent(in) :: m_ang
      character(*), intent(in) :: fname

      integer(ip) :: i, j
      !-----------------------------------------------------------------------
      select_field: select case (fname)
         !--------------------------------------------------------------------
         case ("bianchi3_res")

            call set_level(5_ip,m_ang,psi4_f)
            call set_level(5_ip,m_ang,psi3)
            call set_level(5_ip,m_ang,la)

            call set_thorn(     5_ip,m_ang,psi4_f)
            call set_edth_prime(5_ip,m_ang,psi3)

            do j=1,ny
            do i=1,nx
               bianchi3_res%np1(i,j,m_ang) = &
                                                psi3%edth_prime(i,j,m_ang) !&
         !                                 R(i)*psi3%edth_prime(i,j,m_ang) !&
         !      +    4.0_rp*(R(i)**2)*pi_0(i,j)*psi3%level(i,j,m_ang) &
         !      -                               psi4_f%thorn(i,j,m_ang) &
         !      +                     rh_0(i,j)*psi4_f%level(i,j,m_ang) &
         !      -  3.0_rp*(R(i)**2)*psi2_0(i,j)*la%level(i,j,m_ang) 
            end do
            end do
         !--------------------------------------------------------------------
         case ("bianchi2_res")
            do j=1,ny
            do i=1,nx
               bianchi2_res%np1(i,j,m_ang) = 0.0_rp 
            end do
            end do
         !--------------------------------------------------------------------
         case ("hll_res")
            do j=1,ny
            do i=1,nx
               hll_res%np1(i,j,m_ang) = 0.0_rp 
            end do
            end do
         !-------------------------------------------------------------------- 
         case default
            continue
      end select select_field
   end subroutine set_indep_res
!=============================================================================
   subroutine step_all_fields(step, m_ang)
      integer(ip), intent(in) :: step, m_ang
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,psi4_f)
      call set_edth( step,m_ang,psi4_f)

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
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,psi2)
      call set_level(step,m_ang,hlmb)

      call set_edth( step,m_ang,pi)
      call set_edth( step,m_ang,hlmb)

      call set_level(step,-m_ang,hmbmb)
      call set_edth( step,-m_ang,hmbmb)
      call set_level(step,-m_ang, hlmb)
      call set_level(step,-m_ang,   pi)

      call take_step(step,m_ang,muhll)
      !-----------------------------------------------------------------------
      ! check if there were any errors in take_step routines
      !-----------------------------------------------------------------------
      if (psi4_f%error/="") then
         write (*,*) psi4_f%error
         stop
      else if (psi3%error/="") then
         write (*,*) psi3%error
         stop
      else if (psi2%error/="") then
         write (*,*) psi2%error
         stop
      else if (la%error/="") then
         write (*,*) la%error
         stop
      else if (pi%error/="") then
         write (*,*) pi%error
         stop
      else if (hmbmb%error/="") then
         write (*,*) hmbmb%error 
         stop
      else if (hlmb%error/="") then
         write (*,*) hlmb%error 
         stop
      else if (muhll%error/="") then
         write (*,*) muhll%error
         stop
      else
         continue
      end if

   end subroutine step_all_fields
!=============================================================================
   subroutine metric_recon_time_step(m_ang)
      integer(ip), intent(in) :: m_ang

      call step_all_fields(1_ip,m_ang)
      call step_all_fields(2_ip,m_ang)
      call step_all_fields(3_ip,m_ang)
      call step_all_fields(4_ip,m_ang)
      call step_all_fields(5_ip,m_ang)

   end subroutine metric_recon_time_step
!=============================================================================
   subroutine metric_recon_indep_res(m_ang)
      integer(ip), intent(in) :: m_ang
      !-----------------------------------------------------------------------
      call set_indep_res(m_ang,"bianchi3_res")
      !-----------------------------------------------------------------------
      call set_indep_res(m_ang,"bianchi2_res")

   end subroutine metric_recon_indep_res
!=============================================================================
end module mod_metric_recon 
