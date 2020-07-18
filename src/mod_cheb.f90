module mod_cheb 
!-----------------------------------------------------------------------------
!   use, intrinsic :: iso_c_binding ! for fftw

   use mod_prec
   use mod_io, only: set_arr
   use mod_params, only: tables_dir, R_max, nx, ny, min_m, max_m 

   implicit none
!-----------------------------------------------------------------------------
!   include 'fftw3.f03'
!-----------------------------------------------------------------------------
   private

   ! radial points
   real(rp), dimension(nx), protected, public :: R = 0

   ! subroutines
   public :: cheb_init, compute_DR, cheb_test

   ! Chebyshev differentiation matrix  
   real(rp), dimension(nx,nx) :: D_cheb = 0
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
   subroutine cheb_init()
      call set_arr('cheb_pts.txt', nx,     R)
      call set_arr('cheb_D.txt',nx,nx,D_cheb)

      D_cheb = (2.0_rp/R_max) * D_cheb

      R = (R_max/2.0_rp) * (R + 1.0_rp)

   end subroutine cheb_init
!-----------------------------------------------------------------------------
   pure subroutine compute_DR(m_ang,vals,D_vals)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: D_vals 
      integer(ip) :: i, j

      D_vals = 0
      do j=1,nx
      do i=1,nx
         D_vals(i,:,m_ang) = D_vals(i,:,m_ang) + (D_cheb(i,j) * vals(j,:,m_ang))
      end do
      end do
   end subroutine compute_DR
!-----------------------------------------------------------------------------
   subroutine cheb_test()
      complex(rp), dimension(nx,ny,min_m:max_m) :: vals, DR_vals, computed_DR_vals
      integer(ip) :: i

      do i=1,nx
         vals(i,:,:) = sin(R(i))
         DR_vals(i,:,:) = cos(R(i))
      end do

      call compute_DR(0_ip, vals, computed_DR_vals)

      do i=1,nx
         write (*,*) computed_DR_vals(i,:,0_ip) - DR_vals(i,:,0_ip)
      end do

   end subroutine cheb_test
!-----------------------------------------------------------------------------
end module mod_cheb
