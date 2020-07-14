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
   integer(ip) :: i, j

   real(rp), dimension(2,3), parameter :: test = reshape([ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ],[2,3])

   type(Teuk) :: psi4_m
   psi4_m = teuk_constructor()
!-----------------------------------------------------------------------------
   call cheb_init()
   call swal_init()
   call bkgrd_np_init()
!-----------------------------------------------------------------------------
   y_loop: do j=1,ny
   x_loop: do i=1,nx
      write(*,*) mu_0(i,j)
   end do x_loop
   end do y_loop 

   write (*,*) test
!-----------------------------------------------------------------------------
end block clean_memory
!-----------------------------------------------------------------------------
end program main
