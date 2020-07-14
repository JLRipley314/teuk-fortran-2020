module mod_cheb 
!-----------------------------------------------------------------------------
   use mod_prec
   use mod_io, only: set_arr
   use mod_params, only: dir_tables, r_max, nx, ny 

   implicit none
!-----------------------------------------------------------------------------
   private
   public :: R, cheb_init, cheb_der

   ! radial points
   real(rp), protected, dimension(nx) :: R = 0

   ! Chebyshev differentiation matrix  
   real(rp), dimension(nx,nx) ::  D_cheb = 0
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
   subroutine cheb_init()
      implicit none
      integer(ip) :: i

      call set_arr('pts_cheb.txt', nx,     R)
      call set_arr('D_cheb.txt',nx,nx,D_cheb)

      D_cheb = (2.0_rp/r_max) * D_cheb

      do i=1,nx
         R(i) = (r_max/2.0_rp) * (R(i) + 1.0_rp)
      end do

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
