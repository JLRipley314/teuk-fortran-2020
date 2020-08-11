module mod_cheb 
!=============================================================================
   use mod_prec
   use mod_field,  only: field, set_level
   use mod_io,     only: set_arr
   use mod_params, only: tables_dir, R_max, nx, ny, min_m, max_m 

   implicit none
!=============================================================================
   private

   ! radial points
   real(rp), dimension(nx),    protected, public :: Rvec = 0
   real(rp), dimension(nx,ny), protected, public :: Rarr = 0

   ! subroutines
   public :: cheb_init, compute_DR, cheb_filter, cheb_test

   interface compute_DR
      module procedure compute_DR_arr, compute_DR_field
   end interface

   ! Chebyshev matrices  
   real(rp), dimension(nx,nx) :: D_cheb  = 0
   real(rp), dimension(nx,nx) :: to_cheb = 0
   real(rp), dimension(nx,nx) :: to_real = 0

   ! For radial smoothing 
   real(rp), parameter :: width = R_max/real(nx,rp)

!=============================================================================
contains
!=============================================================================
   subroutine cheb_init()
      integer :: j

      call set_arr('cheb_pts.txt', nx,  Rvec)
      call set_arr('cheb_D.txt',nx,nx,D_cheb)

      call set_arr('real_to_cheb.txt',nx,nx,to_cheb)
      call set_arr('cheb_to_real.txt',nx,nx,to_real)

      D_cheb = (2.0_rp/R_max) * D_cheb

      Rvec = (R_max/2.0_rp) * (Rvec + 1.0_rp)

      do j=1,ny
         Rarr(:,j) = Rvec
      end do
   end subroutine cheb_init
!=============================================================================
   subroutine cheb_real_to_coef(m_ang,vals,coefs)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: coefs

      integer(ip) :: i,j 

      coefs(:,:,m_ang) = 0.0_rp

      do i=1,nx
      do j=1,nx
         coefs(i,:,m_ang) = coefs(i,:,m_ang) + to_cheb(i,j)*vals(j,:,m_ang) 
      end do
      end do
   end subroutine cheb_real_to_coef
!=============================================================================
   subroutine cheb_coef_to_real(m_ang,vals,coefs)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: coefs
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: vals 

      integer(ip) :: i,j 

      vals(:,:,m_ang) = 0.0_rp

      do i=1,nx
      do j=1,nx
         vals(i,:,m_ang) = vals(i,:,m_ang) + to_real(i,j)*coefs(j,:,m_ang) 
      end do
      end do
   end subroutine cheb_coef_to_real
!=============================================================================
   subroutine compute_DR_arr(m_ang,vals,D_vals)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: D_vals 
      integer(ip) :: i, j

      D_vals(:,:,m_ang) = 0

      do j=1,nx
      do i=1,nx
         D_vals(i,:,m_ang) = D_vals(i,:,m_ang) + (D_cheb(i,j) * vals(j,:,m_ang))
      end do
      end do
   end subroutine compute_DR_arr
!=============================================================================
   subroutine compute_DR_field(step,m_ang,f)
      integer(ip), intent(in) :: step, m_ang
      type(field), intent(inout) :: f

      call set_level(step,m_ang,f)

      call compute_DR_arr(m_ang, f%level, f%DR)

   end subroutine compute_DR_field
!=============================================================================
! Low pass filter. A smooth filter appears to help prevent Gibbs-like ringing
!=============================================================================
   subroutine cheb_filter(m_ang, f)
      integer(ip), intent(in) :: m_ang
      type(field), intent(inout) :: f

      integer(ip) :: i

      call cheb_real_to_coef(m_ang,f%np1,f%coefs_cheb) 

      do i=1,nx
         f%coefs_cheb(i,:,m_ang) = &
            exp(-36.0_rp*(real(i-1,rp)/real(nx-1,rp))**25)*f%coefs_cheb(i,:,m_ang)
      end do

      call cheb_coef_to_real(m_ang,f%coefs_cheb,f%np1) 
   end subroutine cheb_filter
!=============================================================================
   subroutine cheb_test()
      complex(rp), dimension(nx,ny,min_m:max_m) :: vals, DR_vals, computed_DR_vals
      integer(ip) :: i

      do i=1,nx
         vals(i,:,:) = sin(Rvec(i))
         DR_vals(i,:,:) = cos(Rvec(i))
      end do

      call compute_DR(0_ip, vals, computed_DR_vals)

      do i=1,nx
         write (*,*) computed_DR_vals(i,:,0_ip) - DR_vals(i,:,0_ip)
      end do

   end subroutine cheb_test
!=============================================================================
end module mod_cheb
