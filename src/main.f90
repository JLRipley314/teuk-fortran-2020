program main
!-----------------------------------------------------------------------------
   use mod_prec
   use mod_cheb
   use mod_swal
   use mod_field
   use mod_params
   use mod_teuk, only: Teuk, teuk_constructor
   use mod_bkgrd_np

   implicit none
!-----------------------------------------------------------------------------
! so valgrind doesn't get confused about automatically deallocated memory
clean_memory: block
!-----------------------------------------------------------------------------
   integer(ip) :: t

   type(Teuk) :: psi4_pm
   psi4_pm = Teuk(pm_ang)
!-----------------------------------------------------------------------------
   call cheb_init()
   call swal_init()
   call bkgrd_np_init()
!-----------------------------------------------------------------------------
! integrate Teukolsky field in time
!-----------------------------------------------------------------------------
   time_evolve: do t=1,nt

   end do time_evolve
!-----------------------------------------------------------------------------
end block clean_memory
!-----------------------------------------------------------------------------
end program main
