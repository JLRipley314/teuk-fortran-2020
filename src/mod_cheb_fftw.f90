module mod_cheb 
!=============================================================================
   use mod_prec
   use mod_field,  only: field, set_level
   use mod_io,     only: set_arr
   use mod_params, only: tables_dir, R_max, nx, ny, min_m, max_m 

!=============================================================================
   use, intrinsic :: iso_c_binding

   implicit none
!=============================================================================
   include 'fftw3.f03'
!=============================================================================
   private

   real(rp), parameter :: dx_over_dR = 2.0_rp / R_max

   integer(ip), parameter :: N = nx - 1

   ! fftw3 objects
   type(c_ptr) :: plan_dct

   ! radial points
   real(rp), dimension(nx),    protected, public :: Rvec
   real(rp), dimension(nx,ny), protected, public :: Rarr

   ! subroutines
   public :: cheb_init, compute_DR, cheb_filter, cheb_test

   interface compute_DR
      module procedure compute_DR_arr, compute_DR_field
   end interface
!=============================================================================
contains
!=============================================================================
   subroutine cheb_init()
      integer :: j
      complex(rp), dimension(nx) :: test_in,test_out

      call set_arr('cheb_pts.txt', nx, Rvec)

      Rvec = (R_max/2.0_rp) * (Rvec + 1.0_rp)

      do j=1,ny
         Rarr(:,j) = Rvec
      end do

      ! setup dct fftw plan
      call dfftw_plan_r2r_1d( &
         plan_dct, &
         nx, &
         test_in,test_out, &
         FFTW_REDFT00,FFTW_PATIENT)

   end subroutine cheb_init
!=============================================================================
   subroutine cheb_real_to_coef(m_ang,vals,coefs)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: coefs

      integer(ip) :: i,j 
      real(rp), dimension(nx) :: rep, imp

      do j=1,ny
         rep = real( vals(:,j,m_ang),kind=rp)
         imp = aimag(vals(:,j,m_ang))

         call dfftw_execute_r2r(plan_dct,rep,rep)
         call dfftw_execute_r2r(plan_dct,imp,imp)

         rep(1) = rep(1)/(2.0_rp*N)
         imp(1) = imp(1)/(2.0_rp*N)
         
         rep(nx) = rep(nx)/(2.0_rp*N)
         imp(nx) = imp(nx)/(2.0_rp*N)

         do i=2,N
            rep(i) = rep(i)/N
            imp(i) = imp(i)/N
         end do

         coefs(:,j,m_ang) = cmplx(rep,imp,kind=rp)
      end do

   end subroutine cheb_real_to_coef
!=============================================================================
   subroutine cheb_coef_to_real(m_ang,coefs,vals)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: coefs
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: vals

      integer(ip) :: i, j 
      real(rp), dimension(nx) :: rep, imp

      do j=1,ny
         rep = real( coefs(:,j,m_ang),kind=rp)
         imp = aimag(coefs(:,j,m_ang))

         do i=2,N
            rep(i) = rep(i) / 2.0_rp
            imp(i) = imp(i) / 2.0_rp
         end do

         call dfftw_execute_r2r(plan_dct,rep,rep)
         call dfftw_execute_r2r(plan_dct,imp,imp)

         vals(:,j,m_ang) = cmplx(rep,imp,kind=rp)
      end do

   end subroutine cheb_coef_to_real
!=============================================================================
   subroutine compute_DR_arr(m_ang,vals,coefs,DR)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: coefs
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: DR

      integer(ip) :: i
   
      call cheb_real_to_coef(m_ang,vals,coefs)

      coefs(nx,  :,m_ang) = 0
      coefs(nx-1,:,m_ang) = 0

      ! use D_vals as a temporary array
      DR(:,:,m_ang) = coefs(:,:,m_ang)

      ! note fortran indexing 1..nx (and not 0..nx-1)
      do i=nx-1, 2, -1
         coefs(i-1,:,m_ang) = 2.0_rp*(i-1)*DR(i,:,m_ang) + coefs(i+1,:,m_ang)
      end do

      coefs(1,:,m_ang) = coefs(1,:,m_ang) / 2.0_rp

      call cheb_coef_to_real(m_ang,coefs,DR)

      DR(:,:,m_ang) = dx_over_dR * DR(:,:,m_ang)

   end subroutine compute_DR_arr
!=============================================================================
   subroutine compute_DR_field(step,m_ang,f)
      integer(ip), intent(in) :: step, m_ang
      type(field), intent(inout) :: f

      call set_level(step,m_ang,f)
      call compute_DR_arr(m_ang,f%level,f%coefs_cheb,f%DR)

   end subroutine compute_DR_field
!=============================================================================
! Low pass filter. A smooth filter appears to help prevent Gibbs-like ringing
!=============================================================================
   subroutine cheb_filter(m_ang, f)
      integer(ip), intent(in)    :: m_ang
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
      complex(rp), dimension(nx,ny,min_m:max_m) :: vals, coefs, DR_vals, computed_DR_vals
      integer(ip) :: i

      integer(ip), parameter :: m_ang = 0_ip

      do i=1,nx
         vals(i,:,:)    = sin(Rvec(i))**2
         DR_vals(i,:,:) = 2*sin(Rvec(i))*cos(Rvec(i))
      end do

      call compute_DR(m_ang, vals, coefs, computed_DR_vals)

      do i=1,nx
         write (*,*) computed_DR_vals(i,:,m_ang) - DR_vals(i,:,m_ang)
      end do


   end subroutine cheb_test
!=============================================================================
end module mod_cheb
