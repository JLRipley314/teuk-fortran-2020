module mod_cheb 
!-----------------------------------------------------------------------------
   use mod_prec
   use mod_io, only: set_arr
   use mod_sim_params, only: dir_tables, nx, ny, R 

   implicit none
!-----------------------------------------------------------------------------
   private
   public :: cheb_init, cheb_der

   ! Chebyshev matrix and
   ! Chebyshev differentiation matrix  
   real(rp), dimension(nx,nx) ::    cheb = 0
   real(rp), dimension(nx,nx) ::  D_cheb = 0
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
   subroutine cheb_init()
      call set_arr('cheb.txt',  nx,nx,  cheb)
      call set_arr('D_cheb.txt',nx,nx,D_cheb)
   end subroutine cheb_init
!-----------------------------------------------------------------------------
   pure subroutine cheb_der(vals,D_vals)
      complex(rp), dimension(nx,nx), intent(in)  :: vals
      complex(rp), dimension(nx,nx), intent(out) :: D_vals 
      integer(ip) :: i, j

      D_vals = 0
      do j=1,nx
      do i=1,nx
         D_vals(i,:) = D_vals(i,:) + (D_cheb(i,j) * vals(j,:))
      end do
      end do
   end subroutine cheb_der 
!-----------------------------------------------------------------------------
end module mod_cheb
