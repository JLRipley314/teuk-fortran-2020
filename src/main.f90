program main
!-----------------------------------------------------------------------------
  use mod_prec
  use mod_field
  use mod_sim_params
  use mod_teuk
  use mod_bkgrd_np

  implicit none
!-----------------------------------------------------------------------------
! so valgrind doesn't get confused about automatically deallocated memory
clean_memory: block
!-----------------------------------------------------------------------------
  integer(ip) :: i, j
  complex(rp) :: val
!-----------------------------------------------------------------------------
  call init_bkgrd_np()
!-----------------------------------------------------------------------------
  y_loop: do j=1,ny
    x_loop: do i=1,nx

      val = mu_0(i,j)

      print *, val

    end do x_loop
  end do y_loop 
!-----------------------------------------------------------------------------
  call clear_bkgrd_np()
!-----------------------------------------------------------------------------
end block clean_memory
!-----------------------------------------------------------------------------
end program main
