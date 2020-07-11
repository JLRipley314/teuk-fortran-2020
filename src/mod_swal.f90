module mod_swal 
!-----------------------------------------------------------------------------
   use mod_prec

   use mod_sim_params, only: dir_tables, &
   nx, ny, lmax, &
   min_m, max_m, min_s, max_s

   implicit none
!-----------------------------------------------------------------------------
   private
   public :: swal_init, swal_lower, swal_raise

   ! weights for Gaussian integration 
   real(rp), dimension(ny) :: weights

   ! Note: range of indices
   real(rp), dimension(1:ny, 0:lmax, min_m:max_m, min_s:max_s) :: swal = 0
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
   subroutine swal_init()
      call swal_set('spin_p2.txt',swal(:,:,:, 2))
      call swal_set('spin_p1.txt',swal(:,:,:, 1))
      call swal_set('spin_0.txt', swal(:,:,:, 0))
      call swal_set('spin_n1.txt',swal(:,:,:,-1))
      call swal_set('spin_n2.txt',swal(:,:,:,-2))
      call swal_set('spin_n3.txt',swal(:,:,:,-3))
   end subroutine swal_init
!-----------------------------------------------------------------------------
   subroutine swal_set(fn, arr)
      character(*),                   intent(in)  :: fn
      real(rp), dimension(ny,0:lmax), intent(out) :: arr

      character(:), allocatable :: rn
      integer(ip) :: ierror
      integer(ip) :: uf = 3
      ! set the file name to read from
      rn = dir_tables // fn

      ! Note: here we ASSUME the input file is correctly formatted
      open(unit=uf,file=rn,status='old',action='read',iostat=ierror)
         if (ierror/=0) then
            write (*,*) "Error: ierror=", ierror
            write (*,*) "file = ", rn
            stop
         end if
         read (uf,*,iostat=ierror) arr
      close(uf)
   end subroutine swal_set
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

      ! synthesis
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
end module mod_swal
