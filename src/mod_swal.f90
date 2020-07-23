!
! Spin-weighted associated Legendre function routines, along with
! Legendre polynomial roots and weights for Gaussian qudrature to
! go to/from coefficient and real space
!
module mod_swal 
!=============================================================================
   use mod_prec
   use mod_io, only: set_arr
   use mod_params, only: &
      nx, ny, nl, lmax, &
      min_m, max_m, min_s, max_s

   implicit none
!=============================================================================
   private

   ! gauss points y, cos(y), sin(y)
   real(rp), dimension(ny), protected, public :: Y, cy, sy

   ! subroutines 
   public :: swal_init, swal_laplacian, swal_lower, swal_raise, swal_filter, swal_test

   ! weights for Gaussian integration 
   real(rp), dimension(ny) :: weights

   ! Note: range of indices
   real(rp), dimension(1:ny, 0:lmax, min_m:max_m, min_s:max_s), protected, public :: swal = 0
!=============================================================================
contains
!=============================================================================
   integer(ip) pure function compute_lmin(spin, m_ang) result(lmin)
      integer(ip), intent(in)  :: spin, m_ang

      lmin = max(abs(spin),abs(m_ang))
   end function compute_lmin
!=============================================================================
   subroutine swal_init()
      integer(ip) :: m_ang = 0
      integer(ip) :: spin  = 0
      character(:), allocatable :: mstr, sstr

      ! weights for Gaussian quadrature
      call set_arr('roots_legendre.txt', ny, Y)

      cy = cos(Y)
      sy = sin(Y)

      ! weights for Gaussian quadrature
      call set_arr('weights_legendre.txt', ny, weights)

      ! spin-weighted spherical associated Legendre polynomials 
      do spin=min_s,max_s

         sstr = '     '
         write (sstr,'(i5)') spin
         sstr = trim(adjustl(sstr))

         do m_ang=min_m,max_m

            mstr = '     '
            write (mstr,'(i5)') m_ang
            mstr = trim(adjustl(mstr))

            call set_arr('s_'//sstr//'_m_' //mstr//'.txt', ny,nl,swal(:,:,m_ang,spin))
         end do
      end do
   end subroutine swal_init
!=============================================================================
! Gaussian quadrature.
! Integrate against complex conjugate cc[s_Y^m_l]=(-1)^{m+2} {-s}_Y^{-m}_l
!=============================================================================
   pure subroutine swal_real_to_coef(spin,m_ang,vals,coefs)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,    ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,0:lmax,min_m:max_m), intent(out) :: coefs

      integer(ip) :: lmin, j, k

      lmin  = compute_lmin(spin,m_ang)
      coefs = 0.0_rp

      do k=lmin,lmax
      do j=1,ny
         coefs(:,k,m_ang) =  &
            coefs(:,k,m_ang) &
         +  (vals(:,j,m_ang) * weights(j) * ((-1.0_rp)**(spin+m_ang))*swal(j,k,-m_ang,-spin))
      end do
      end do
   end subroutine swal_real_to_coef
!=============================================================================
! coefficient synthesis
!=============================================================================
   pure subroutine swal_coef_to_real(spin,m_ang,coefs,vals)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,0:lmax,min_m:max_m), intent(in)  :: coefs
      complex(rp), dimension(nx,    ny,min_m:max_m), intent(out) :: vals

      integer(ip) :: lmin, j, k

      lmin = compute_lmin(spin,m_ang)
      vals = 0.0_rp

      do j=1,ny
      do k=lmin,lmax
         vals(:,j,m_ang) = &
            vals(:,j,m_ang) &
         + (coefs(:,k,m_ang) * swal(j,k,m_ang,spin))
      end do
      end do
   end subroutine swal_coef_to_real
!=============================================================================
   pure subroutine swal_laplacian(spin,m_ang,vals,coefs,vals_lap)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,    ny,min_m:max_m), intent(in)    :: vals
      complex(rp), dimension(nx,0:lmax,min_m:max_m), intent(inout) :: coefs
      complex(rp), dimension(nx,    ny,min_m:max_m), intent(out)   :: vals_lap

      integer(ip) :: k

      call swal_real_to_coef(spin,m_ang,vals,coefs) 

      do k=0,lmax
         coefs(:,k,m_ang) = &
         - real(k-spin,rp)*real(k+spin+1.0_rp,rp)*coefs(:,k,m_ang)
      end do

      call swal_coef_to_real(spin,m_ang,coefs,vals_lap) 
   end subroutine swal_laplacian
!=============================================================================
   pure subroutine swal_lower(spin,m_ang,vals,coefs,vals_lowered)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,    ny,min_m:max_m), intent(in)    :: vals
      complex(rp), dimension(nx,0:lmax,min_m:max_m), intent(inout) :: coefs
      complex(rp), dimension(nx,    ny,min_m:max_m), intent(out)   :: vals_lowered 
      integer(ip) :: lmin, k

      call swal_real_to_coef(spin,m_ang,vals,coefs) 

      lmin = max(-spin,spin-1_ip)

      do k=lmin,lmax
         coefs(:,k,m_ang) = &
         -  sqrt(real(k+spin,rp)*real(k-spin+1.0_rp,rp))*coefs(:,k,m_ang)
      end do

      call swal_coef_to_real(spin-1,m_ang,coefs,vals_lowered) 
   end subroutine swal_lower
!=============================================================================
   pure subroutine swal_raise(spin,m_ang,vals,coefs,vals_raised)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,    ny,min_m:max_m), intent(in)    :: vals
      complex(rp), dimension(nx,0:lmax,min_m:max_m), intent(inout) :: coefs
      complex(rp), dimension(nx,    ny,min_m:max_m), intent(out)   :: vals_raised

      integer(ip) :: lmin, k

      call swal_real_to_coef(spin,m_ang,vals,coefs) 

      lmin = max(spin,-spin-1_ip)

      do k=lmin,lmax
         coefs(:,k,m_ang) = &
            sqrt(real(k-spin,rp)*real(k+spin+1.0_rp,rp))*coefs(:,k,m_ang)
      end do

      call swal_coef_to_real(spin+1,m_ang,coefs,vals_raised) 
   end subroutine swal_raise
!=============================================================================
! Low pass filter. A smooth filter helps prevent Gibbs-like ringing
!=============================================================================
   pure subroutine swal_filter(spin,m_ang,vals,coefs)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,    ny,min_m:max_m), intent(inout) :: vals
      complex(rp), dimension(nx,0:lmax,min_m:max_m), intent(inout) :: coefs

      integer(ip) :: k

      call swal_real_to_coef(spin,m_ang,vals,coefs) 

      do k=0,lmax
         coefs(:,k,m_ang) = &
            exp(-36.0_rp*(real(k,rp)/real(lmax,rp))**25)*coefs(:,k,m_ang)
      end do

      call swal_coef_to_real(spin,m_ang,coefs,vals) 
   end subroutine swal_filter
!=============================================================================
! test that the swal are orthogonal
!=============================================================================
   subroutine swal_test()

      integer(ip) :: s, m, l1, l2, j
      real(rp)    :: integral

      do s=-2,2
      do m=min_m,max_m
      do l1=0,lmax
      do l2=0,lmax
         integral = 0.0_rp
         do j=1,ny 
            integral = integral + weights(j)*swal(j,l1,m,s)*(((-1.0_rp)**(m+s))*swal(j,l2,-m,-s))
         end do
         write (*,*) s,m,l1,l2,integral
      end do 
      end do 
      end do 
      end do 

   end subroutine swal_test
!=============================================================================
end module mod_swal
