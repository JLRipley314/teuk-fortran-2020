!
! Evolve linear Teukolsky field, reconstruct metric, and evolve
! second order Teukolsky field
!
!=============================================================================
program main
!=============================================================================
   use, intrinsic :: iso_fortran_env, only: &
      stdout=>output_unit, stdin=>input_unit, stderr=>error_unit

   use mod_prec
   use mod_params,       only: nt, dt, t_step_save, black_hole_mass, &
                           psi_spin, psi_boost, pm_ang, &
                           metric_recon, &
                           write_metric_recon_fields, &
                           write_indep_res!, write_source
   use mod_field,        only: set_field, shift_time_step
   use mod_cheb,         only: cheb_init, cheb_filter, cheb_test
   use mod_swal,         only: swal_init, swal_filter, swal_test_orthonormal
   use mod_io,           only: write_csv
   use mod_teuk,         only: teuk_init, teuk_time_step, compute_res_q
   use mod_initial_data, only: set_initial_data
   use mod_bkgrd_np,     only: bkgrd_np_init
   use mod_metric_recon, only: 

   use mod_fields_list, only: &
      psi4_lin_p, psi4_lin_q, psi4_lin_f, &
      res_lin_q, & 

      psi3, psi2, la, pi, muhll, hlmb, hmbmb, &
      res_bianchi3, res_bianchi2, res_hll

   implicit none
!=============================================================================
! so valgrind doesn't get confused about automatically deallocated memory
clean_memory: block
!=============================================================================
! declare and initialize variables, fields, etc.
!=============================================================================
   integer(ip) :: t_step
   real(rp)    :: time
!=============================================================================
   write (*,*) "Initializing fields"   
!-----------------------------------------------------------------------------
! first order metric field
!-----------------------------------------------------------------------------
   call set_field(name="lin_p",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_lin_p)
   call set_field(name="lin_q",spin=psi_spin,boost=psi_boost,falloff=2_ip,f=psi4_lin_q)
   call set_field(name="lin_f",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_lin_f)
!-----------------------------------------------------------------------------
! metric reconstructed fields
!-----------------------------------------------------------------------------
   call set_field(name="psi3",spin=-1_ip,boost=-1_ip,falloff=2_ip,f=psi3)
   call set_field(name="psi2",spin= 0_ip,boost= 0_ip,falloff=3_ip,f=psi2)

   call set_field(name="la",spin=-2_ip,boost=-1_ip,falloff=1_ip,f=la)
   call set_field(name="pi",spin=-1_ip,boost= 0_ip,falloff=2_ip,f=pi)

   call set_field(name="muhll",spin= 0_ip,boost=1_ip,falloff=3_ip,f=muhll)
   call set_field(name="hlmb" ,spin=-1_ip,boost=1_ip,falloff=2_ip,f= hlmb)
   call set_field(name="hmbmb",spin=-2_ip,boost=0_ip,falloff=1_ip,f=hmbmb)
!-----------------------------------------------------------------------------
! independent residual fields
!-----------------------------------------------------------------------------
   call set_field(name="res_lin_q",spin=-2_ip,boost=-2_ip,falloff=2_ip,f=res_lin_q)

   call set_field(name="res_bianchi3",spin=-2_ip,boost=-1_ip,falloff=2_ip,f=res_bianchi3)
   call set_field(name="res_bianchi2",spin=-1_ip,boost= 0_ip,falloff=2_ip,f=res_bianchi2)
   call set_field(name="res_hll",     spin= 0_ip,boost= 2_ip,falloff=2_ip,f=res_hll)
!-----------------------------------------------------------------------------
! initialize chebyshev diff matrices, swal matrices, etc.
!-----------------------------------------------------------------------------
   call cheb_init()
   call swal_init()
   call bkgrd_np_init()
   call teuk_init()
!=============================================================================
! initial data 
!=============================================================================
   write (stdout,*) "Setting up initial data"
!-----------------------------------------------------------------------------
   time = 0.0_rp

   call set_initial_data(-pm_ang, psi4_lin_p, psi4_lin_q, psi4_lin_f)
   call set_initial_data( pm_ang, psi4_lin_p, psi4_lin_q, psi4_lin_f)
   !--------------------------------------------------------------------------
   ! write to file
   !--------------------------------------------------------------------------
   call write_csv(time,-pm_ang,psi4_lin_f)
   call write_csv(time, pm_ang,psi4_lin_f)
   !--------------------------------------------------------------------------
   if (write_metric_recon_fields) then
      call write_csv(time,-pm_ang,psi3)
      call write_csv(time, pm_ang,psi3)

      call write_csv(time,-pm_ang,psi2)
      call write_csv(time, pm_ang,psi2)

      call write_csv(time,-pm_ang,hmbmb)
      call write_csv(time, pm_ang,hmbmb)

      call write_csv(time,-pm_ang,hlmb)
      call write_csv(time, pm_ang,hlmb)

      call write_csv(time,-pm_ang,muhll)
      call write_csv(time, pm_ang,muhll)
   end if
   !--------------------------------------------------------------------------
   if (write_indep_res) then
      call compute_res_q(psi4_lin_q,psi4_lin_f,res_lin_q)

      call write_csv(time,-pm_ang,res_lin_q)
      call write_csv(time, pm_ang,res_lin_q)

      if (metric_recon) then
         call write_csv(time,-pm_ang,res_bianchi3)
         call write_csv(time, pm_ang,res_bianchi3)

         call write_csv(time,-pm_ang,res_bianchi2)
         call write_csv(time, pm_ang,res_bianchi2)

         call write_csv(time,-pm_ang,res_hll)
         call write_csv(time, pm_ang,res_hll)
      end if
   end if
!=============================================================================
! integrate in time 
!=============================================================================
!-----------------------------------------------------------------------------
   write (stdout,*) "Beginning time evolution"
!-----------------------------------------------------------------------------
   time_evolve: do t_step=1,nt
      time = t_step*dt
      !-----------------------------------------------------------------------
      call teuk_time_step(-pm_ang, psi4_lin_p, psi4_lin_q, psi4_lin_f)
      call teuk_time_step( pm_ang, psi4_lin_p, psi4_lin_q, psi4_lin_f)
      !-----------------------------------------------------------------------
      if (metric_recon) then 
         call metric_recon_time_step(pm_ang)
      end if
      !-----------------------------------------------------------------------
      if (mod(t_step,t_step_save)==0) then
         !--------------------------------------------------------------------
         write (stdout,*) time / black_hole_mass
         flush (stdout)
         !--------------------------------------------------------------------
         call write_csv(time,-pm_ang,psi4_lin_f)
         call write_csv(time, pm_ang,psi4_lin_f)
         !--------------------------------------------------------------------
         if (write_metric_recon_fields) then
            call write_csv(time,-pm_ang,psi3)
            call write_csv(time, pm_ang,psi3)

            call write_csv(time,-pm_ang,psi2)
            call write_csv(time, pm_ang,psi2)

            call write_csv(time,-pm_ang,hmbmb)
            call write_csv(time, pm_ang,hmbmb)

            call write_csv(time,-pm_ang,hlmb)
            call write_csv(time, pm_ang,hlmb)

            call write_csv(time,-pm_ang,muhll)
            call write_csv(time, pm_ang,muhll)
         end if
         !--------------------------------------------------------------------
         if (write_indep_res) then
            call compute_res_q(psi4_lin_q,psi4_lin_f,res_lin_q)

            call write_csv(time,-pm_ang,res_lin_q)
            call write_csv(time, pm_ang,res_lin_q)
            !-----------------------------------------------------------------
            if (metric_recon) then
               call metric_recon_indep_res(-pm_ang)
               call metric_recon_indep_res( pm_ang)

               call write_csv(time,-pm_ang,res_bianchi3)
               call write_csv(time, pm_ang,res_bianchi3)

               call write_csv(time,-pm_ang,res_bianchi2)
               call write_csv(time, pm_ang,res_bianchi2)

               call write_csv(time,-pm_ang,res_hll)
               call write_csv(time, pm_ang,res_hll)
            end if
         end if
         !--------------------------------------------------------------------
      end if
      !-----------------------------------------------------------------------
      call cheb_filter(psi4_lin_p%np1,psi4_lin_p%coefs_cheb)
      call cheb_filter(psi4_lin_q%np1,psi4_lin_q%coefs_cheb)
      call cheb_filter(psi4_lin_f%np1,psi4_lin_f%coefs_cheb)

      call swal_filter(psi4_lin_p%spin,psi4_lin_p%np1,psi4_lin_p%coefs_swal)
      call swal_filter(psi4_lin_q%spin,psi4_lin_q%np1,psi4_lin_q%coefs_swal)
      call swal_filter(psi4_lin_f%spin,psi4_lin_f%np1,psi4_lin_f%coefs_swal)

      if (metric_recon) then
         call cheb_filter(psi3%np1,psi3%coefs_cheb)
         call cheb_filter(psi2%np1,psi2%coefs_cheb)

         call cheb_filter(la%np1,la%coefs_cheb)
         call cheb_filter(pi%np1,pi%coefs_cheb)

         call cheb_filter(hmbmb%np1,hmbmb%coefs_cheb)
         call cheb_filter( hlmb%np1, hlmb%coefs_cheb)
         call cheb_filter(muhll%np1,muhll%coefs_cheb)
      end if
      !-----------------------------------------------------------------------
      call shift_time_step(psi4_lin_p)
      call shift_time_step(psi4_lin_q)
      call shift_time_step(psi4_lin_f)

      if (metric_recon) then
         call shift_time_step(psi3)
         call shift_time_step(psi2)

         call shift_time_step(la)
         call shift_time_step(pi)

         call shift_time_step(hmbmb)
         call shift_time_step(hlmb)
         call shift_time_step(muhll)
      end if
   end do time_evolve
!=============================================================================
end block clean_memory
!=============================================================================
end program main
