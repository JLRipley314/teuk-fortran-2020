module mod_swal 
!-----------------------------------------------------------------------------
   use mod_prec
   use mod_io, only: set_arr
   use mod_params, only: &
      nx, ny, nl, &
      min_m, max_m, min_s, max_s

   implicit none
!-----------------------------------------------------------------------------
   private

   ! maximum l angular value
   integer(ip), parameter :: lmax = nl-1

   ! gauss points y, cos(y), sin(y)
   real(rp), dimension(ny), protected, public :: Y, cy, sy

   ! subroutines 
   public :: swal_init, swal_laplacian, swal_lower, swal_raise, swal_write

   ! weights for Gaussian integration 
   real(rp), dimension(ny) :: weights

   ! Note: range of indices
   real(rp), dimension(1:ny, 0:lmax, min_m:max_m, min_s:max_s), protected, public :: swal = 0
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
   subroutine swal_init()
      integer(ip) :: m = 0
      character(:), allocatable :: mstr

      ! weights for Gaussian quadrature
      call set_arr('roots_legendre.txt', ny, Y)

      cy = cos(Y)
      sy = sin(Y)

      ! weights for Gaussian quadrature
      call set_arr('weights_legendre.txt', ny, weights)

      ! spin-weighted spherical associated Legendre polynomials 
      do m=min_m,max_m

         ! inelegant int to str conversion
         mstr = '     '
         write (mstr,'(i5)') m
         mstr = trim(adjustl(mstr))

         call set_arr('s_2_m_' //mstr//'.txt', ny, nl, swal(:,:,m,  2))
         call set_arr('s_1_m_' //mstr//'.txt', ny, nl, swal(:,:,m,  1))
         call set_arr('s_0_m_' //mstr//'.txt', ny, nl, swal(:,:,m,  0))
         call set_arr('s_-1_m_'//mstr//'.txt', ny, nl, swal(:,:,m, -1))
         call set_arr('s_-2_m_'//mstr//'.txt', ny, nl, swal(:,:,m, -2))
         call set_arr('s_-3_m_'//mstr//'.txt', ny, nl, swal(:,:,m, -3))
      end do
   end subroutine swal_init
!-----------------------------------------------------------------------------
   pure subroutine swal_real_to_coef(spin,m_ang,vals,coefs)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny),     intent(in)  :: vals
      complex(rp), dimension(nx,0:lmax), intent(out) :: coefs

      integer(ip) :: j, k

      ! Gaussian quadrature
      coefs = 0
      do k=0,lmax
      do j=1,ny
         coefs(:,k) =  &
            coefs(:,k) &
         +  (vals(:,j) * weights(j) * swal(j,k,m_ang,spin))
      end do
      end do
   end subroutine swal_real_to_coef
!-----------------------------------------------------------------------------
   pure subroutine swal_coef_to_real(spin,m_ang,coefs,vals)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,0:lmax), intent(in)  :: coefs
      complex(rp), dimension(nx,ny),     intent(out) :: vals

      integer(ip) :: j, k

      ! coefficient synthesis
      vals = 0
      do j=1,ny
      do k=0,lmax
         vals(:,j) = &
            vals(:,j) &
         + (coefs(:,k) * swal(j,k,m_ang,spin))
      end do
      end do
   end subroutine swal_coef_to_real
!-----------------------------------------------------------------------------
   pure subroutine swal_laplacian(spin,m_ang,vals, vals_lap)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny), intent(in)  :: vals
      complex(rp), dimension(nx,ny), intent(out) :: vals_lap

      real(rp)    :: pre
      integer(ip) :: k

      complex(rp), dimension(nx,0:lmax) :: coefs

      call swal_real_to_coef(spin,m_ang,vals,coefs) 

      do k=0,lmax
         pre = - (k - spin) * (k + spin + 1.0_rp)
         coefs(:,k) = pre*coefs(:,k)
      end do

      call swal_coef_to_real(spin,m_ang,coefs,vals_lap) 
   end subroutine swal_laplacian
!-----------------------------------------------------------------------------
   pure subroutine swal_lower(spin,m_ang,coefs,vals)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,0:lmax), intent(inout) :: coefs
      complex(rp), dimension(nx,ny),     intent(inout) :: vals

      real(rp)    :: pre
      integer(ip) :: k

      call swal_real_to_coef(spin,m_ang,vals,coefs) 

      do k=0,lmax
         pre = -sqrt( &
               (real(k,rp)+real(spin,rp)) &
            *  (real(k,rp)-real(spin,rp)+1.0_rp) &
            ) 
         coefs(:,k) = pre*coefs(:,k)
      end do

      call swal_coef_to_real(spin-1,m_ang,coefs,vals) 
   end subroutine swal_lower
!-----------------------------------------------------------------------------
   pure subroutine swal_raise(spin,m_ang,coefs,vals)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,0:lmax), intent(inout) :: coefs
      complex(rp), dimension(nx,ny),     intent(inout) :: vals

      real(rp)    :: pre
      integer(ip) :: k

      call swal_real_to_coef(spin,m_ang,vals,coefs) 

      do k=0,lmax
         pre = sqrt( &
               (real(k,rp)-real(spin,rp)) &
            *  (real(k,rp)+real(spin,rp)+1.0_rp) &
            ) 
         coefs(:,k) = pre*coefs(:,k)
      end do

      call swal_coef_to_real(spin+1,m_ang,coefs,vals) 
   end subroutine swal_raise
!-----------------------------------------------------------------------------
   subroutine swal_write()
      integer(ip) :: j

      do j=1,ny
         write (*,*) swal(j,:,0,0)
      end do
   end subroutine swal_write
!-----------------------------------------------------------------------------
end module mod_swal
