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
   use mod_params, only: &
      nt, dt, t_step_save, black_hole_mass, &
      lin_m, scd_m, &
      psi_spin, psi_boost, &
      metric_recon, scd_order, &
      scd_order_start_time

   use mod_field,        only: set_field, shift_time_step
   use mod_cheb,         only: cheb_init, cheb_filter, cheb_test
   use mod_swal,         only: swal_init, swal_filter, swal_test_orthonormal
   use mod_teuk,         only: teuk_init, teuk_time_step
   use mod_initial_data, only: set_initial_data
   use mod_bkgrd_np,     only: bkgrd_np_init
   use mod_metric_recon, only: metric_recon_time_step 
   use mod_write_level,  only: write_level

   use mod_fields_list, only: &
      psi4_lin_p, psi4_lin_q, psi4_lin_f, &
      res_lin_q, & 

      psi4_scd_p, psi4_scd_q, psi4_scd_f, &
      res_scd_q, & 

      psi3, psi2, la, pi, muhll, hlmb, hmbmb, &
      res_bianchi3, res_bianchi2, res_hll

   use mod_scd_order_source, only: &
      scd_order_source, &
      source, &
      scd_order_source_init, &
      scd_order_source_compute, &
      scd_order_source_shift_time_step

   implicit none
!=============================================================================
! so valgrind doesn't get confused about automatically deallocated memory
clean_memory: block
!=============================================================================
! declare and initialize variables, fields, etc.
!=============================================================================
   integer(ip) :: i, t_step
   real(rp)    :: time
!=============================================================================
   write (*,*) "Initializing fields"   
!-----------------------------------------------------------------------------
! first order metric field
!-----------------------------------------------------------------------------
   call set_field(name="lin_p",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_lin_p)
   call set_field(name="lin_q",spin=psi_spin,boost=psi_boost,falloff=2_ip,f=psi4_lin_q)
   call set_field(name="lin_f",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_lin_f)

   if (scd_order) then
      call set_field(name="scd_p",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_scd_p)
      call set_field(name="scd_q",spin=psi_spin,boost=psi_boost,falloff=2_ip,f=psi4_scd_q)
      call set_field(name="scd_f",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_scd_f)
   end if
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

   if (scd_order) then
      call set_field(name="res_scd_q",spin=-2_ip,boost=-2_ip,falloff=2_ip,f=res_scd_q)
   end if

   call set_field(name="res_bianchi3",spin=-2_ip,boost=-1_ip,falloff=2_ip,f=res_bianchi3)
   call set_field(name="res_bianchi2",spin=-1_ip,boost= 0_ip,falloff=2_ip,f=res_bianchi2)
   call set_field(name="res_hll",     spin= 0_ip,boost= 2_ip,falloff=2_ip,f=res_hll)
!-----------------------------------------------------------------------------
! source term for \psi_4^{(2)}
!-----------------------------------------------------------------------------
   if (scd_order) then
      call scd_order_source_init(name="scd_order_source",sf=source)
   end if
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

   do i=1,size(lin_m)
      call set_initial_data(lin_m(i), psi4_lin_p, psi4_lin_q, psi4_lin_f)
   end do
   call write_level(time)
!=============================================================================
! integrate in time 
!=============================================================================
   write (stdout,*) "Beginning time evolution"
!-----------------------------------------------------------------------------
   time_evolve: do t_step=1,nt
      time = t_step*dt
      !-----------------------------------------------------------------------
      ! \Psi_4^{(1)} evolution 
      !-----------------------------------------------------------------------
      do i=1,size(lin_m)
         call teuk_time_step(lin_m(i),psi4_lin_p,psi4_lin_q,psi4_lin_f)
      end do
      !-----------------------------------------------------------------------
      ! metric recon evolves +/- m_ang
      !-----------------------------------------------------------------------
      if (metric_recon) then 
         do i=1,size(lin_m)
            if (lin_m(i)>=0) then
               call metric_recon_time_step(lin_m(i))
            end if
         end do
      end if
      !-----------------------------------------------------------------------
      ! \Psi_4^{(2)} evolution 
      !-----------------------------------------------------------------------
      if (scd_order) then
         do i=1,size(scd_m)
            call scd_order_source_compute(scd_m(i),source) 
         end do
         if (time>=scd_order_start_time) then
            do i=1,size(scd_m)
               call teuk_time_step(scd_m(i),source,psi4_scd_p,psi4_scd_q,psi4_scd_f)
            end do
         end if
      end if
      !-----------------------------------------------------------------------
      ! save to file 
      !-----------------------------------------------------------------------
      save_level: if (mod(t_step,t_step_save)==0) then
         !--------------------------------------------------------------------
         write (stdout,*) time / black_hole_mass
         flush (stdout)
         call write_level(time)
      end if save_level
      !-----------------------------------------------------------------------
      ! low pass filter (in spectral space)
      !-----------------------------------------------------------------------
      do i=1,size(lin_m)
         call cheb_filter(lin_m(i),psi4_lin_p)
         call cheb_filter(lin_m(i),psi4_lin_q)
         call cheb_filter(lin_m(i),psi4_lin_f)
      end do

      do i=1,size(scd_m)
         call swal_filter(scd_m(i),psi4_lin_p)
         call swal_filter(scd_m(i),psi4_lin_q)
         call swal_filter(scd_m(i),psi4_lin_f)
      end do

      if (metric_recon) then
         do i=1,size(lin_m)
            call cheb_filter(lin_m(i),psi3)
            call cheb_filter(lin_m(i),psi2)

            call cheb_filter(lin_m(i),la)
            call cheb_filter(lin_m(i),pi)

            call cheb_filter(lin_m(i),hmbmb)
            call cheb_filter(lin_m(i), hlmb)
            call cheb_filter(lin_m(i),muhll)
         end do
      end if
      !-----------------------------------------------------------------------
      ! shift time steps
      !-----------------------------------------------------------------------
      do i=1,size(lin_m)
         call shift_time_step(lin_m(i),psi4_lin_p)
         call shift_time_step(lin_m(i),psi4_lin_q)
         call shift_time_step(lin_m(i),psi4_lin_f)
      end do

      if (metric_recon) then
         do i=1,size(lin_m)
            call shift_time_step(lin_m(i),psi3)
            call shift_time_step(lin_m(i),psi2)

            call shift_time_step(lin_m(i),la)
            call shift_time_step(lin_m(i),pi)

            call shift_time_step(lin_m(i),hmbmb)
            call shift_time_step(lin_m(i),hlmb)
            call shift_time_step(lin_m(i),muhll) 
         end do 
      end if

      if (scd_order) then
         do i=1,size(scd_m)
            call shift_time_step(scd_m(i),psi4_scd_p)
            call shift_time_step(scd_m(i),psi4_scd_q)
            call shift_time_step(scd_m(i),psi4_scd_f)

            call scd_order_source_shift_time_step(scd_m(i),source)
         end do
      end if

   end do time_evolve
!=============================================================================
end block clean_memory
!=============================================================================
end program main
