!
! Evolve linear Teukolsky field, reconstruct metric, and evolve
! second order Teukolsky field
!
!=============================================================================
program main
!=============================================================================
   use mod_prec
   use mod_params,       only: nt, t_step_save, pm_ang
   use mod_field,        only: field
   use mod_cheb,         only: cheb_init
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
   integer(ip) :: t

   type(field) :: p, q, f
!=============================================================================
   p = field("p")
   q = field("q")
   f = field("f")
!=============================================================================
   call cheb_init()
   call swal_init()
   call bkgrd_np_init()
   call teuk_init()
!=============================================================================
! integrate Teukolsky field in time
!=============================================================================
   call set_initial_data(p, q, f)

   call write_csv(p)
   call write_csv(q)
   call write_csv(f)
   !--------------------------------------------------------------------------
   time_evolve: do t=1,nt

      call teuk_time_step(pm_ang, p, q, f)
      !-----------------------------------------------------------------------
      if (mod(t,t_step_save)==0) then
         write (*,*) t
         call write_csv(p)
         call write_csv(q)
         call write_csv(f)
      end if
      !-----------------------------------------------------------------------
   end do time_evolve
!=============================================================================
end block clean_memory
!=============================================================================
end program main
