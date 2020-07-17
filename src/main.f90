!
! Evolve linear Teukolsky field, reconstruct metric, and evolve
! second order Teukolsky field
!
!=============================================================================
program main
!=============================================================================
   use mod_prec
   use mod_params,       only: nt, dt, t_step_save, pm_ang
   use mod_field,        only: field, shift_time_step
   use mod_cheb,         only: cheb_init, cheb_test
   use mod_swal,         only: swal_init
   use mod_bkgrd_np,     only: bkgrd_np_init
   use mod_teuk,         only: teuk_init, teuk_time_step
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
   type(field) :: p, q, f
!=============================================================================
   write (*,*) "Initializing fields"   
!-----------------------------------------------------------------------------
   p = field("p")
   q = field("q")
   f = field("f")
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

   call set_initial_data(p, q, f)

   call write_csv(time,p)
   call write_csv(time,q)
   call write_csv(time,f)
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
         write (*,*) time
         call write_csv(time,p)
         call write_csv(time,q)
         call write_csv(time,f)
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
