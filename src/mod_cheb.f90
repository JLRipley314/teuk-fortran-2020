module mod_cheb 
!-----------------------------------------------------------------------------
   use mod_prec
   use mod_io, only: set_arr
   use mod_params, only: tables_dir, R_max, nx, ny 

   implicit none
!-----------------------------------------------------------------------------
   private

   ! radial points
   real(rp), dimension(nx), protected, public :: R = 0

   ! subroutines
   public :: cheb_init, compute_DR

   ! Chebyshev differentiation matrix  
   real(rp), dimension(nx,nx) :: D_cheb = 0
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
   subroutine cheb_init()
      implicit none

      call set_arr('cheb_pts.txt', nx,     R)
      call set_arr('cheb_D.txt',nx,nx,D_cheb)

      D_cheb = (2.0_rp/R_max) * D_cheb

      R = (R_max/2.0_rp) * (R + 1.0_rp)

   end subroutine cheb_init
!-----------------------------------------------------------------------------
   pure subroutine compute_DR(vals,D_vals)
      complex(rp), dimension(nx,ny), intent(in)  :: vals
      complex(rp), dimension(nx,ny), intent(out) :: D_vals 
      integer(ip) :: i, j

      D_vals = 0
      do j=1,nx
      do i=1,nx
         D_vals(i,:) = D_vals(i,:) + (D_cheb(i,j) * vals(j,:))
      end do
      end do
   end subroutine compute_DR
!-----------------------------------------------------------------------------
end module mod_cheb
