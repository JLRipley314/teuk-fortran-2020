module mod_cheb 
!=============================================================================
   use mod_prec
   use mod_fftw3
   use mod_io, only: set_arr
   use mod_params, only: tables_dir, R_max, nx, ny, min_m, max_m 

   implicit none
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
         nx,test_in,test_out, &
         FFTW_REDFT00,FFTW_PATIENT)

   end subroutine cheb_init
!=============================================================================
   subroutine cheb_real_to_coef(m_ang,vals,coefs)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: coefs

      integer(ip) :: i,j,N 
      real(rp), dimension(nx) :: rep, imp

      N = nx-1

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

      integer(ip) :: j,N 
      real(rp), dimension(nx) :: rep, imp

      N = nx-1

      do j=1,ny
         rep = real( coefs(:,j,m_ang),kind=rp)
         imp = aimag(coefs(:,j,m_ang))

         rep = rep / 2.0_rp
         imp = imp / 2.0_rp

         call dfftw_execute_r2r(plan_dct,rep,rep)
         call dfftw_execute_r2r(plan_dct,imp,imp)

         vals(:,j,m_ang) = cmplx(rep,imp,kind=rp)
      end do
   end subroutine cheb_coef_to_real
!=============================================================================
   subroutine compute_DR(m_ang,vals,coefs,D_vals)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)    :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: coefs
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: D_vals 

      integer(ip) :: i
   
      call cheb_real_to_coef(m_ang,vals,coefs)

      ! use D_vals as a temporary array
      D_vals(:,:,m_ang) = coefs(:,:,m_ang)

      coefs(nx,  :,m_ang) = 0
      coefs(nx-1,:,m_ang) = 0

      do i=nx-1, 2, -1
         coefs(i-1,:,m_ang) = 2.0_rp*i*D_vals(i,:,m_ang) + coefs(i+1,:,m_ang)
      end do

      coefs(1,:,m_ang) = coefs(1,:,m_ang) / 2.0_rp

      call cheb_coef_to_real(m_ang,coefs,D_vals)

      D_vals(:,:,m_ang) = dx_over_dR * D_vals(:,:,m_ang)

   end subroutine compute_DR
!=============================================================================
! Low pass filter. A smooth filter appears to help prevent Gibbs-like ringing
!=============================================================================
   subroutine cheb_filter(vals,coefs)
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
      complex(rp), dimension(nx,ny,min_m:max_m) :: vals, coefs, DR_vals, computed_DR_vals
      integer(ip) :: i

      do i=1,nx
         vals(i,:,:) = sin(Rvec(i))
         DR_vals(i,:,:) = cos(Rvec(i))
      end do

      call compute_DR(0_ip, vals, coefs, computed_DR_vals)

      do i=1,nx
         write (*,*) computed_DR_vals(i,:,0_ip) - DR_vals(i,:,0_ip)
      end do

   end subroutine cheb_test
!=============================================================================
end module mod_cheb
