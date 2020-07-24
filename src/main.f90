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
   use mod_params,       only: nt, dt, t_step_save, black_hole_mass, pm_ang, &
                           write_fields, write_indep_res
   use mod_field,        only: field, set_field, shift_time_step
   use mod_cheb,         only: cheb_init, cheb_test
   use mod_swal,         only: swal_init, swal_test_orthonormal
   use mod_io,           only: write_csv
   use mod_teuk,         only: psi4_f, psi4_p, psi4_q, q_res, &
                               teuk_init, teuk_time_step, compute_q_res
   use mod_initial_data, only: set_initial_data
   use mod_bkgrd_np,     only: bkgrd_np_init
   use mod_metric_recon, only: psi3, psi2, la, pi, muhll, hlmb, hmbmb, &
                               bianchi3_res, bianchi2_res, hll_res,    &
                               metric_recon_time_step, &
                               metric_recon_indep_res

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
   call set_field(name="p",spin=-2_ip,boost=-2_ip,falloff=1_ip,f=psi4_p)
   call set_field(name="q",spin=-2_ip,boost=-2_ip,falloff=2_ip,f=psi4_q)
   call set_field(name="f",spin=-2_ip,boost=-2_ip,falloff=1_ip,f=psi4_f)
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
   call set_field(name="res_q",spin=-2_ip,boost=-2_ip,falloff=2_rp,f=q_res)

   call set_field(name="bianchi3_res",spin=-2_ip,boost=-1_ip,falloff=2_rp,f=bianchi3_res)
   call set_field(name="bianchi2_res",spin=-1_ip,boost= 0_ip,falloff=2_rp,f=bianchi2_res)
   call set_field(name="hll_res",     spin= 0_ip,boost= 2_ip,falloff=2_rp,f=hll_res)
!-----------------------------------------------------------------------------
! initialize chebyshev diff matrices, swal matrices, etc.
!-----------------------------------------------------------------------------
   call cheb_init()
   call swal_init()
   call bkgrd_np_init()
   call teuk_init()
!=============================================================================
! integrate Teukolsky field in time
!=============================================================================
   write (stdout,*) "Setting up initial data"
!-----------------------------------------------------------------------------
   time = 0.0_rp

   call set_initial_data(pm_ang,psi4_p, psi4_q, psi4_f)
   !--------------------------------------------------------------------------
   ! write to file
   !--------------------------------------------------------------------------
   call write_csv(time,pm_ang,psi4_p)
   call write_csv(time,pm_ang,psi4_q)
   call write_csv(time,pm_ang,psi4_f)
   !--------------------------------------------------------------------------
   if (write_fields) then
      call write_csv(time,pm_ang,psi3)
      call write_csv(time,pm_ang,psi2)

      call write_csv(time,pm_ang,muhll)
      call write_csv(time,pm_ang, hlmb)
      call write_csv(time,pm_ang,hmbmb)
   end if
   !--------------------------------------------------------------------------
   if (write_indep_res) then
      call compute_q_res(psi4_q,psi4_f,q_res)

      call write_csv(time,pm_ang,q_res)

      call write_csv(time,pm_ang,bianchi3_res)
      call write_csv(time,pm_ang,bianchi2_res)
      call write_csv(time,pm_ang,hll_res)
   end if
!-----------------------------------------------------------------------------
   write (stdout,*) "Beginning time evolution"
!-----------------------------------------------------------------------------
   time_evolve: do t_step=1,nt
      time = t_step*dt

      call teuk_time_step(pm_ang, psi4_p, psi4_q, psi4_f)
   
      call metric_recon_time_step(pm_ang)
      !-----------------------------------------------------------------------
      if (mod(t_step,t_step_save)==0) then
         !--------------------------------------------------------------------
         write (stdout,*) time / black_hole_mass
         flush (stdout)
         !--------------------------------------------------------------------
         call write_csv(time,pm_ang,psi4_f)
         !--------------------------------------------------------------------
         if (write_fields) then
            call write_csv(time,pm_ang,psi3)
            call write_csv(time,pm_ang,psi2)

            call write_csv(time,pm_ang,muhll)
            call write_csv(time,pm_ang, hlmb)
            call write_csv(time,pm_ang,hmbmb)
         end if
         !--------------------------------------------------------------------
         if (write_indep_res) then
            call compute_q_res(psi4_q,psi4_f,q_res)

            call write_csv(time,pm_ang,q_res)
            !-----------------------------------------------------------------
            call metric_recon_indep_res(pm_ang)

            call write_csv(time,pm_ang,bianchi3_res)
            call write_csv(time,pm_ang,bianchi2_res)
            call write_csv(time,pm_ang,hll_res)
         end if
         !--------------------------------------------------------------------
      end if
      !-----------------------------------------------------------------------
      call shift_time_step(psi4_p)
      call shift_time_step(psi4_q)
      call shift_time_step(psi4_f)

      call shift_time_step(psi3)
      call shift_time_step(psi2)

      call shift_time_step(la)
      call shift_time_step(pi)

      call shift_time_step(hmbmb)
      call shift_time_step(hlmb)
      call shift_time_step(muhll)

   end do time_evolve
!=============================================================================
end block clean_memory
!=============================================================================
end program main
