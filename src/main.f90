!
! Evolve linear Teukolsky field, reconstruct metric, and evolve
! second order Teukolsky field
!
!=============================================================================
program main
!=============================================================================
   use mod_prec
   use mod_params,       only: nt, dt, t_step_save, black_hole_mass, pm_ang
   use mod_field,        only: field, set_field, shift_time_step
   use mod_cheb,         only: cheb_init, cheb_test
   use mod_swal,         only: swal_init
   use mod_io,           only: write_csv
   use mod_teuk,         only: psi4_f, psi4_p, psi4_q, &
                               teuk_init, teuk_time_step, compute_q_res
   use mod_initial_data, only: set_initial_data
   use mod_bkgrd_np,     only: bkgrd_np_init
   use mod_metric_recon, only: psi3, psi2, la, pi, muhll, hlmb, hmbmb, &
                               bianci3_res, bianci2_res, hll_res,      &
                               metric_recon_time_step

   implicit none
!=============================================================================
! so valgrind doesn't get confused about automatically deallocated memory
clean_memory: block
!=============================================================================
! declare variables, fields, etc.
!=============================================================================
   integer(ip) :: t_step
   real(rp) :: time
   type(field) :: q_res
!=============================================================================
! fields initialization format: name, spin, boost, falloff, field
!-----------------------------------------------------------------------------
   write (*,*) "Initializing fields"   
!-----------------------------------------------------------------------------
! first order metric field
!-----------------------------------------------------------------------------
   call set_field("p",-2_ip,2_ip,1_ip,psi4_p)
   call set_field("q",-2_ip,2_ip,2_ip,psi4_q)
   call set_field("f",-2_ip,2_ip,1_ip,psi4_f)
!-----------------------------------------------------------------------------
! metric reconstructed fields
!-----------------------------------------------------------------------------
   call set_field("psi3",-1_ip,1_ip,2_ip,psi3)
   call set_field("psi2", 0_ip,0_ip,3_ip,psi2)

   call set_field("la",-2_ip,-1_ip,1_ip,la)
   call set_field("pi",-1_ip, 0_ip,2_ip,pi)

   call set_field("muhll", 0_ip,2_ip,3_ip,muhll)
   call set_field("hlmb" ,-1_ip,1_ip,2_ip, hlmb)
   call set_field("hmbmb",-2_ip,0_ip,1_ip,hmbmb)
!-----------------------------------------------------------------------------
! independent residual fields
!-----------------------------------------------------------------------------
   call set_field("q_res",-2_ip,2_ip,2_rp,q_res)

   call set_field("bianci3_res",-2_ip,2_ip,2_rp,bianci3_res)
   call set_field("bianci2_res",-2_ip,2_ip,2_rp,bianci2_res)
   call set_field("ll_res",     -2_ip,2_ip,2_rp,hll_res)
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
   write (*,*) "Setting up initial data"
!-----------------------------------------------------------------------------
   time = 0.0_rp

   call set_initial_data(pm_ang,psi4_p, psi4_q, psi4_f)

   call compute_q_res(psi4_q,psi4_f,q_res)

   call write_csv(time,pm_ang,psi4_p)
   call write_csv(time,pm_ang,psi4_q)
   call write_csv(time,pm_ang,psi4_f)
   
   call write_csv(time,pm_ang,q_res)
!-----------------------------------------------------------------------------
!   write (*,*) "Performing tests"
!   call cheb_test()
!   stop
!-----------------------------------------------------------------------------
   write (*,*) "Beginning time evolution"
!-----------------------------------------------------------------------------
   time_evolve: do t_step=1,nt
      time = t_step*dt

      call teuk_time_step(pm_ang, psi4_p, psi4_q, psi4_f)
   
      call metric_recon_time_step(pm_ang)
      !-----------------------------------------------------------------------
      if (mod(t_step,t_step_save)==0) then
         write (*,*) time / black_hole_mass

         call compute_q_res(psi4_q,psi4_f,q_res)

         call write_csv(time,pm_ang,psi4_p)
         call write_csv(time,pm_ang,psi4_q)
         call write_csv(time,pm_ang,psi4_f)

         call write_csv(time,pm_ang,q_res)
      end if
      !-----------------------------------------------------------------------
      call shift_time_step(psi4_p)
      call shift_time_step(psi4_q)
      call shift_time_step(psi4_f)
   end do time_evolve
!=============================================================================
end block clean_memory
!=============================================================================
end program main
