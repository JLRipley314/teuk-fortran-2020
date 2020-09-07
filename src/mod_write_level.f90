!
! module for writing fields to file
!
!=============================================================================
module mod_write_level
!=============================================================================
   use, intrinsic :: iso_fortran_env, only: &
      stdout=>output_unit, stdin=>input_unit, stderr=>error_unit

   use mod_prec

   use mod_params, only: nx, ny 

   use mod_io,   only: write_csv

   use mod_cheb, only: cheb_real_to_coef
   use mod_swal, only: swal_real_to_coef

   use mod_teuk, only: compute_res_q

   use mod_metric_recon, only: metric_recon_indep_res

   use mod_params, only: &
      lin_write_m, scd_write_m, &
      len_lin_write_m, len_scd_write_m, &
      pm1_ang, pm2_ang, &
      metric_recon, scd_order, &
      write_metric_recon_fields, &
      write_indep_res, write_scd_order_source, &
      write_coefs, constrained_evo

   use mod_fields_list, only: &
      psi4_lin_p, psi4_lin_q, psi4_lin_f, &
      res_lin_q, & 

      psi4_scd_p, psi4_scd_q, psi4_scd_f, &
      res_scd_q, & 

      psi3, psi2, la, pi, muhll, hlmb, hmbmb, &
      res_bianchi3, res_bianchi2, res_hll

   use mod_scd_order_source, only: source
!=============================================================================
   implicit none
   private

   public :: write_level, write_diagnostics
!=============================================================================
contains
!=============================================================================
   subroutine write_diagnostics(time)
      real(rp), intent(in) :: time
      real(rp), dimension(12) :: diag 

      complex(rp) :: mean
      integer(ip) :: m_ang 
      !-----------------------------------------------------------------------
      ! field values at future null infinity and horizon
      !-----------------------------------------------------------------------
      m_ang = maxval(lin_write_m)

      mean = sum(psi4_lin_f%np1(1,:,m_ang)) / ny
      diag(1) = real( mean,kind=rp)
      diag(2) = aimag(mean)

      mean = sum(psi4_lin_f%np1(nx,:,m_ang)) / ny
      diag(3) = real( mean,kind=rp)
      diag(4) = aimag(mean)
      !-----------------------------------------------------------------------
      m_ang = maxval(scd_write_m)

      mean = sum(psi4_scd_f%np1(1,:,m_ang)) / ny
      diag(5) = real( mean,kind=rp)
      diag(6) = aimag(mean)

      mean = sum(psi4_scd_f%np1(nx,:,m_ang)) / ny
      diag(7) = real( mean,kind=rp)
      diag(8) = aimag(mean)
      !-----------------------------------------------------------------------
      m_ang = 0 

      mean = sum(psi4_scd_f%np1(1,:,0)) / ny
      diag(9)  = real( mean,kind=rp)
      diag(10) = aimag(mean)

      mean = sum(psi4_scd_f%np1(nx,:,0)) / ny
      diag(11) = real( mean,kind=rp)
      diag(12) = aimag(mean)
      !-----------------------------------------------------------------------
      write(stdout,'(e14.6)',advance='no') time
      write(stdout,'(e14.6)',advance='no') diag(1)
      write(stdout,'(e14.6)',advance='no') diag(2)
      write(stdout,'(e14.6)',advance='no') diag(3)
      write(stdout,'(e14.6)',advance='no') diag(4)
      write(stdout,'(e14.6)',advance='no') diag(5)
      write(stdout,'(e14.6)',advance='no') diag(6)
      write(stdout,'(e14.6)',advance='no') diag(7)
      write(stdout,'(e14.6)',advance='no') diag(8)
      write(stdout,'(e14.6)',advance='no') diag(9)
      write(stdout,'(e14.6)',advance='no') diag(10)
      write(stdout,'(e14.6)',advance='no') diag(11)
      write(stdout,'(e14.6)',advance='no') diag(12)
      ! line break
      write(stdout,*)
   end subroutine write_diagnostics
!=============================================================================
   subroutine write_level(time)
      real(rp), intent(in) :: time

      integer(ip) :: i
      !-----------------------------------------------------------------------
      ! \Psi_4^{(1)} and linear metric reconstruction 
      !--------------------------------------------------------------------------
      do i=1,len_lin_write_m
         call write_csv(time,lin_write_m(i),psi4_lin_f)
      end do 

      if (write_metric_recon_fields) then
         do i=1,len_lin_write_m
            call write_csv(time,lin_write_m(i),psi3)
            call write_csv(time,lin_write_m(i),psi2)
            call write_csv(time,lin_write_m(i),hmbmb)
            call write_csv(time,lin_write_m(i),hlmb)
            call write_csv(time,lin_write_m(i),muhll)
         end do
      end if

      if (write_indep_res) then
         if (.not. constrained_evo) then
            do i=1,len_lin_write_m
               call compute_res_q( lin_write_m(i),psi4_lin_q,psi4_lin_f,res_lin_q)
               call write_csv(time,lin_write_m(i),res_lin_q)
            end do
         end if
         if (metric_recon) then
            do i=1,len_lin_write_m
               call metric_recon_indep_res(lin_write_m(i))

               call write_csv(time,lin_write_m(i),res_bianchi3)
               call write_csv(time,lin_write_m(i),res_bianchi2)
               call write_csv(time,lin_write_m(i),res_hll)
            end do
         end if
      end if
      !-----------------------------------------------------------------------
      ! \Psi_4^{(2)} and 2nd order source term 
      !-----------------------------------------------------------------------
      if (scd_order) then

         do i=1,len_scd_write_m
            call write_csv(time,scd_write_m(i),psi4_scd_f)
         end do

         if (write_indep_res) then
            if (.not. constrained_evo) then
               do i=1,len_scd_write_m
                  call compute_res_q( scd_write_m(i),psi4_scd_q,psi4_scd_f,res_scd_q)
                  call write_csv(time,scd_write_m(i),res_scd_q)
               end do
            end if
         end if

         if (write_scd_order_source) then 
            do i=1,len_scd_write_m
               call write_csv( &
                  source%fname, &
                  time, &
                  scd_write_m(i), &
                  source%np1(:,:,scd_write_m(i)) &
               )
            end do
         end if

      end if
      !-----------------------------------------------------------------------
      if (write_coefs) then
         !--------------------------------------------------------------------
         do i=1,len_lin_write_m
            call cheb_real_to_coef( &
               lin_write_m(i), &
               psi4_lin_f%np1, &
               psi4_lin_f%coefs_cheb &
            )
            call swal_real_to_coef( &
               psi4_lin_f%spin, &
               lin_write_m(i), &
               psi4_lin_f%coefs_cheb, &
               psi4_lin_f%coefs_both &
            )
            call write_csv( &
               "coefs_"//psi4_lin_f%fname, &
               time, &
               lin_write_m(i), &
               psi4_lin_f%coefs_both(:,:,lin_write_m(i)) &
            )
         end do 
         !--------------------------------------------------------------------
         if (scd_order) then
            do i=1,len_scd_write_m
               call cheb_real_to_coef( &
                  scd_write_m(i), &
                  psi4_scd_f%np1, &
                  psi4_scd_f%coefs_cheb &
               )
               call swal_real_to_coef( &
                  psi4_scd_f%spin, &
                  scd_write_m(i), &
                  psi4_scd_f%coefs_cheb, &
                  psi4_scd_f%coefs_both &
               )
               call write_csv( &
                  "coefs_"//psi4_scd_f%fname, &
                  time, &
                  scd_write_m(i), &
                  psi4_scd_f%coefs_both(:,:,scd_write_m(i)) &
               )
            end do 
         end if
         !--------------------------------------------------------------------
      end if
      !--------------------------------------------------------------------
   end subroutine write_level
!=============================================================================
end module mod_write_level
