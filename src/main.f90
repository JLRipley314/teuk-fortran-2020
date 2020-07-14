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

   real(rp), dimension(0:1,3), parameter :: test = reshape([&
      1.0_rp,&
      2.0_rp,&
      3.0_rp,&
      4.0_rp,&
      5.0_rp,&
      6.0_rp&
      ],[2,3])

   type(Teuk) :: psi4_m
   psi4_m = teuk_constructor()
!-----------------------------------------------------------------------------
!   call cheb_init()
!   call swal_init()
   call bkgrd_np_init()
!-----------------------------------------------------------------------------
   y_loop: do j=1,ny
   x_loop: do i=1,nx
      write(*,*) mu_0(i,j)
   end do x_loop
   end do y_loop 

   write (*,*) test(0,1)
   write (*,*) test(1,1)

   write (*,*) test(0,2)
   write (*,*) test(1,2)

   write (*,*) test(0,3)
   write (*,*) test(1,3)
!-----------------------------------------------------------------------------
end block clean_memory
!-----------------------------------------------------------------------------
end program main
