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
!   integer(ip) :: j

   type(Teuk) :: psi4_m
   psi4_m = teuk_constructor()
!-----------------------------------------------------------------------------
   call cheb_init()
   call swal_init()
   call bkgrd_np_init()
!-----------------------------------------------------------------------------
   call swal_write()

!   y_loop: do j=1,ny
!      write(*,*) mu_0(:,j)
!   end do y_loop 
!-----------------------------------------------------------------------------
end block clean_memory
!-----------------------------------------------------------------------------
end program main
