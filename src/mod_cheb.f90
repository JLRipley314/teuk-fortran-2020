module mod_cheb 
!=============================================================================
   use mod_fftw3

   use mod_prec
   use mod_io, only: set_arr
   use mod_params, only: tables_dir, R_max, nx, ny, min_m, max_m 

   implicit none
!=============================================================================
   private

   ! fftw3 objects
   type(c_ptr) :: plan

   ! radial points
   real(rp), dimension(nx), protected, public :: R = 0

   ! subroutines
   public :: cheb_init, compute_DR, radial_smooth, cheb_test

   ! Chebyshev differentiation matrix  
   real(rp), dimension(nx,nx) :: D_cheb  = 0
   real(rp), dimension(nx)    :: weights = 0

   ! For radial smoothing 
   real(rp), parameter :: width = R_max/real(nx,rp)

!=============================================================================
contains
!=============================================================================
   subroutine cheb_init()

      integer(ip), parameter :: N = 10
      complex(rp), dimension(N) :: test_in,test_out

      call set_arr('cheb_pts.txt', nx,     R)
      call set_arr('cheb_D.txt',nx,nx,D_cheb)
      call set_arr('weights_clenshaw_curtis.txt', nx, weights)

      D_cheb = (2.0_rp/R_max) * D_cheb

      R = (R_max/2.0_rp) * (R + 1.0_rp)

      weights = (2.0_rp/R_max) * weights

      call dfftw_plan_dft_1d(plan,N,test_in,test_out,FFTW_FORWARD,FFTW_ESTIMATE)

   end subroutine cheb_init
!=============================================================================
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
!=============================================================================
   real(rp) pure function smooth(i,j) result(res)
      integer(ip), intent(in) :: i, j

      res = exp(abs((R(i)-R(j))/width)**20)

   end function smooth
!=============================================================================
! Convolve functions with Gaussian in radial direction.
! Radial integration is done via Clenshaw-Curtis quadrature
!=============================================================================
   pure subroutine radial_smooth(m_ang, tmp, vals)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: tmp
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: vals

      integer(ip) :: i, j
      real(rp) :: norm

      tmp  = vals
      vals = 0.0_rp

      do i=1,nx
         norm = 0.0_rp
         do j=1,nx

            vals(i,:,m_ang) = vals(i,:,m_ang) &
            +  weights(abs(i-j)+1_ip)*smooth(i,j)*tmp(j,:,m_ang)

            norm = norm + weights(abs(i-j)+1_ip)*smooth(i,j)

         end do
         vals(i,:,m_ang) = vals(i,:,m_ang) / norm
      end do

   end subroutine radial_smooth
!=============================================================================
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
!=============================================================================
end module mod_cheb
