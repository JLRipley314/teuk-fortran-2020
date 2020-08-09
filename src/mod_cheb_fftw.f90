module mod_cheb 
!=============================================================================
   use mod_prec
   use mod_field,  only: field, set_level, get_level
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

   ! fftw3 objects
   type(c_ptr) :: plan_dct

   ! radial points
   real(rp), dimension(nx),    protected, public :: Rvec = 0
   real(rp), dimension(nx,ny), protected, public :: Rarr = 0

   ! subroutines
   public :: cheb_init, compute_DR, cheb_filter, cheb_test

   ! For radial smoothing 
   real(rp), parameter :: width = R_max/real(nx,rp)
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
   subroutine cheb_real_to_coef(step,m_ang,f)
      integer(ip), intent(in) :: step,m_ang
      type(field), intent(inout) :: f

      integer(ip) :: i,j,N 
      real(rp), dimension(nx) :: rep, imp

      call set_level(step,m_ang,f)

      N = nx-1

      do j=1,ny
         rep = real( f%level(:,j,m_ang),kind=rp)
         imp = aimag(f%level(:,j,m_ang))

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

         f%coefs(:,j,m_ang) = cmplx(rep,imp,kind=rp)
      end do

   end subroutine cheb_real_to_coef
!=============================================================================
   subroutine cheb_coef_to_real(step,m_ang,f)
      integer(ip), intent(in) :: step,m_ang

      integer(ip) :: j,N 
      real(rp), dimension(nx) :: rep, imp

      N = nx-1

      do j=1,ny
         rep = real( f%coefs(:,j,m_ang),kind=rp)
         imp = aimag(f%coefs(:,j,m_ang))

         rep = rep / 2.0_rp
         imp = imp / 2.0_rp

         call dfftw_execute_r2r(plan_dct,rep,rep)
         call dfftw_execute_r2r(plan_dct,imp,imp)

         f%level(:,j,m_ang) = cmplx(rep,imp,kind=rp)

         call get_level(step,m_ang,f)
      end do
   end subroutine cheb_coef_to_real
!=============================================================================
   subroutine compute_DR(step,m_ang,f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      integer(ip) :: i
   
      call cheb_real_to_coef(step,m_ang,f)

      ! use D_vals as a temporary array
      f%DR(:,:,m_ang) = f%coefs_cheb(:,:,m_ang)

      f%coefs_cheb(nx,  :,m_ang) = 0
      f%coefs_cheb(nx-1,:,m_ang) = 0

      do i=nx-1, 2, -1
         f%coefs_cheb(i-1,:,m_ang) = 2.0_rp*i*f%DR(i,:,m_ang) + f%coefs_cheb(i+1,:,m_ang)
      end do

      f%coefs_cheb(1,:,m_ang) = f%coefs_cheb(1,:,m_ang) / 2.0_rp

      call cheb_coef_to_real(step,m_ang,f)

      f%DR(:,:,m_ang) = dx_over_dR * f%DR(:,:,m_ang)

   end subroutine compute_DR
!=============================================================================
! Low pass filter. A smooth filter appears to help prevent Gibbs-like ringing
!=============================================================================
   subroutine cheb_filter(f)
      type(field), intent(inout) :: f

      integer(ip), parameter :: step = 5_ip

      integer(ip) :: m_ang, i
      
      do m_ang=min_m,max_m

         call cheb_real_to_coef(step,m_ang,f) 

         do i=1,nx
            f%coefs(i,:,m_ang) = &
               exp(-36.0_rp*(real(i-1,rp)/real(nx-1,rp))**25)*f%coefs(i,:,m_ang)
         end do

         call cheb_coef_to_real(step,m_ang,f) 

      end do
   end subroutine cheb_filter
!=============================================================================
end module mod_cheb
