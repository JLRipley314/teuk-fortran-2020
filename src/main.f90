!
! Evolve linear Teukolsky field, reconstruct metric, and evolve
! second order Teukolsky field
!
!=============================================================================
program main
!=============================================================================
   use mod_prec
   use mod_cheb
   use mod_swal
   use mod_field
   use mod_params
   use mod_teuk, only: teuk_init, teuk_time_step
   use mod_bkgrd_np

   implicit none
!=============================================================================
! so valgrind doesn't get confused about automatically deallocated memory
clean_memory: block
!=============================================================================
   integer(ip) :: t
!=============================================================================
   call cheb_init()
   call swal_init()
   call bkgrd_np_init()
   call teuk_init()
!=============================================================================
! integrate Teukolsky field in time
!=============================================================================
   time_evolve: do t=1,nt

   end do time_evolve
!=============================================================================
end block clean_memory
!=============================================================================
end program main
