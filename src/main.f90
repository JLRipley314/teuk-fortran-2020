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
   use mod_bkgrd_np,     only: bkgrd_np_init
   use mod_teuk,         only: teuk_init, teuk_time_step, compute_q_indep_res
   use mod_initial_data, only: set_initial_data
   use mod_io,           only: write_csv

   implicit none
!=============================================================================
! so valgrind doesn't get confused about automatically deallocated memory
clean_memory: block
!=============================================================================
! declare variables, fields, etc.
!=============================================================================
   integer(ip) :: t_step
   real(rp) :: time
   type(field) :: p, q, f, q_indep_res
!=============================================================================
   write (*,*) "Initializing fields"   
!-----------------------------------------------------------------------------
   call set_field("p",-2_ip,2_ip,1_ip,p)
   call set_field("q",-2_ip,2_ip,2_ip,q)
   call set_field("f",-2_ip,2_ip,1_ip,f)

   call set_field("q_indep_res",-2_ip,2_ip,2_rp,q_indep_res)
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

   call set_initial_data(pm_ang,p, q, f)

   call compute_q_indep_res(q,f,q_indep_res)

   call write_csv(time,pm_ang,p)
   call write_csv(time,pm_ang,q)
   call write_csv(time,pm_ang,f)
   
   call write_csv(time,pm_ang,q_indep_res)
!-----------------------------------------------------------------------------
!   write (*,*) "Performing tests"
!   call cheb_test()
!   stop
!-----------------------------------------------------------------------------
   write (*,*) "Beginning time evolution"
!-----------------------------------------------------------------------------
   time_evolve: do t_step=1,nt
      time = t_step*dt

      call teuk_time_step(pm_ang, p, q, f)
      !-----------------------------------------------------------------------
      if (mod(t_step,t_step_save)==0) then
         write (*,*) time / black_hole_mass

         call compute_q_indep_res(q,f,q_indep_res)

         call write_csv(time,pm_ang,p)
         call write_csv(time,pm_ang,q)
         call write_csv(time,pm_ang,f)

         call write_csv(time,pm_ang,q_indep_res)
      end if
      !-----------------------------------------------------------------------
      call shift_time_step(p)
      call shift_time_step(q)
      call shift_time_step(f)
   end do time_evolve
!=============================================================================
end block clean_memory
!=============================================================================
end program main
