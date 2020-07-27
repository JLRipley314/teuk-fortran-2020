module mod_cheb 
!=============================================================================
   use mod_prec
   use mod_io, only: set_arr
   use mod_params, only: tables_dir, R_max, nx, ny, min_m, max_m 

   implicit none
!=============================================================================
   private

   ! fftw3 objects
!   type(c_ptr) :: plan_forward, plan_backward

   ! radial points
   real(rp), dimension(nx), protected, public :: R = 0

   ! subroutines
   public :: cheb_init, compute_DR, cheb_filter, cheb_test

   ! Chebyshev matrices  
   real(rp), dimension(nx,nx) :: D_cheb  = 0
   real(rp), dimension(nx,nx) :: to_cheb = 0
   real(rp), dimension(nx,nx) :: to_real = 0
   real(rp), dimension(nx)    :: weights = 0

   ! For radial smoothing 
   real(rp), parameter :: width = R_max/real(nx,rp)

!=============================================================================
contains
!=============================================================================
   subroutine cheb_init()

!      complex(rp), dimension(nx) :: test_in,test_out

      call set_arr('cheb_pts.txt', nx,     R)
      call set_arr('cheb_D.txt',nx,nx,D_cheb)
      call set_arr('weights_clenshaw_curtis.txt', nx, weights)

      call set_arr('real_to_cheb.txt',nx,nx,to_cheb)
      call set_arr('cheb_to_real.txt',nx,nx,to_real)

      D_cheb = (2.0_rp/R_max) * D_cheb

      R = (R_max/2.0_rp) * (R + 1.0_rp)

      weights = (2.0_rp/R_max) * weights

!      call dfftw_plan_dft_1d(plan_forward, nx,test_in,test_out,FFTW_FORWARD, FFTW_ESTIMATE)
!      call dfftw_plan_dft_1d(plan_backward,nx,test_in,test_out,FFTW_BACKWARD,FFTW_ESTIMATE)

   end subroutine cheb_init
!=============================================================================
   pure subroutine cheb_real_to_coef(m_ang,vals,coefs)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: coefs

      integer(ip) :: i,j 

      coefs = 0.0_rp

      do i=1,nx
      do j=1,nx
         coefs(i,:,m_ang) = coefs(i,:,m_ang) + to_cheb(i,j)*vals(j,:,m_ang) 
      end do
      end do
   end subroutine cheb_real_to_coef
!=============================================================================
   pure subroutine cheb_coef_to_real(m_ang,vals,coefs)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: coefs
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: vals 

      integer(ip) :: i,j 

      vals = 0.0_rp

      do i=1,nx
      do j=1,nx
         vals(i,:,m_ang) = vals(i,:,m_ang) + to_real(i,j)*coefs(j,:,m_ang) 
      end do
      end do
   end subroutine cheb_coef_to_real
!=============================================================================
   pure subroutine compute_DR(m_ang,vals,D_vals)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)    :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: D_vals 
      integer(ip) :: i, j

      D_vals(:,:,m_ang) = 0

      do j=1,nx
      do i=1,nx
         D_vals(i,:,m_ang) = D_vals(i,:,m_ang) + (D_cheb(i,j) * vals(j,:,m_ang))
      end do
      end do
   end subroutine compute_DR
!=============================================================================
! Low pass filter. A smooth filter appears to help prevent Gibbs-like ringing
!=============================================================================
   pure subroutine cheb_filter(vals,coefs)
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: coefs

      integer(ip) :: m_ang, i

      do m_ang=min_m,max_m

         call cheb_real_to_coef(m_ang,vals,coefs) 
         do i=1,nx
            coefs(i,:,m_ang) = &
               exp(-36.0_rp*(real(i-1,rp)/real(nx-1,rp))**25)*coefs(i,:,m_ang)
         end do
         call cheb_coef_to_real(m_ang,coefs,vals) 

      end do
   end subroutine cheb_filter
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
