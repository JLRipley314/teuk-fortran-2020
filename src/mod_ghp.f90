!
! Geroch-Held-Penrose (GHP) operators
!
module mod_ghp
!=============================================================================
   use mod_prec
   use mod_cheb,   only: Rvec=>R, compute_DR
   use mod_swal,   only: cyvec=>cy, syvec=>sy, swal_lower, swal_raise
   use mod_field,  only: field, set_level, set_DT

   use mod_bkgrd_np, only: ep_0
   use mod_params,   only: &
      dt, nx, ny, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin
!=============================================================================
   implicit none
   private
   
   complex(rp), parameter :: ci = (0.0_rp, 1.0_rp) 

   public :: ghp_edth, ghp_edth_prime, ghp_thorn, ghp_thorn_prime
!=============================================================================
   contains
!=============================================================================
   pure subroutine ghp_edth(step, f, f_edth)
      integer(ip), intent(in)    :: step
      type(field), intent(inout) :: f
      complex(rp), dimension(nx,ny), intent(out)   :: f_edth

      integer(ip) :: i, j
      real(rp) :: R, cy, sy, p

      p = (f%spin+f%boost)

      call set_level(step, f)
      call set_DT(   step, f)

      call compute_DR(f%level, f%DR)

      call swal_lower(f%spin, f%m_ang, f%level, f%coefs, f%lowered)

      y_loop: do j=1,ny
      x_loop: do i=1,nx
   
         R = Rvec(i)
         cy = cyvec(j)
         sy = syvec(j)

         f_edth(i,j) = &
            (r/sqrt(2.0_rp)) * (1.0_rp/((cl**2)-ci*bhs*r*cy)) * ( &
               - ci*bhs*sy*(f%DT(i,j)) &
               + (f%lowered(i,j)) &
            ) &
         -  ( &
               (ci*p/sqrt(2.0_rp)) * bhs * (R**2) * sy  &
            /  (((cl**2) - ci*bhs*R*cy)**2) &
            ) * (f%level(i,j))

      end do x_loop
      end do y_loop
   end subroutine ghp_edth
!=============================================================================
   pure subroutine ghp_edth_prime(step, f, f_edth_prime)
      integer(ip), intent(in)    :: step
      type(field), intent(inout) :: f
      complex(rp), dimension(nx,ny), intent(out)   :: f_edth_prime

      integer(ip) :: i, j
      real(rp) :: R, cy, sy, q

      q = (-f%spin+f%boost)

      call set_level(step, f)
      call set_DT(   step, f)

      call compute_DR(f%level, f%DR)

      call swal_raise(f%spin, f%m_ang, f%level, f%coefs, f%raised)

      y_loop: do j=1,ny
      x_loop: do i=1,nx

         R = Rvec(i)
         cy = cyvec(j)
         sy = syvec(j)

         f_edth_prime(i,j) = &
            (r/sqrt(2.0_rp)) * (1.0_rp/((cl**2)+ci*bhs*r*cy)) * ( &
                 ci*bhs*sy*(f%DT(i,j)) &
               + (f%raised(i,j)) &
            ) &
         -  ( &
               (ci*q/sqrt(2.0_rp)) * bhs * (R**2) * sy  &
            /  (((cl**2) + ci*bhs*R*cy)**2) &
            ) * (f%level(i,j))

      end do x_loop
      end do y_loop
   end subroutine ghp_edth_prime
!=============================================================================
   pure subroutine ghp_thorn(step, f, f_thorn)
      integer(ip), intent(in)    :: step
      type(field), intent(inout) :: f
      complex(rp), dimension(nx,ny), intent(out)   :: f_thorn

      integer(ip) :: i, j, m_ang
      real(rp)    :: R, cy, p, q
      complex(rp) :: ep, ep_cc

      m_ang = f%m_ang

      p = ( f%spin+f%boost)
      q = (-f%spin+f%boost)

      call set_level(step, f)
      call set_DT(   step, f)

      call compute_DR(f%level, f%DR)

      y_loop: do j=1,ny
      x_loop: do i=1,nx

         R = Rvec(i)
         cy = cyvec(j)

         ep = ep_0(i,j)
         ep_cc = conjg(ep)

         f_thorn(i,j) = &
            (r**2)/((cl**4)+((bhs*r*cy)**2)) * ( &
               2.0_rp*bhm*(2.0_rp*bhm-((bhs/cl)**2)*r) * f%DT(i,j) &
            -  0.5_rp*((cl**2)-(2.0_rp*bhm*r) + ((bhs/cl)**2)*r) * f%DR(i,j) &
            ) &
            +  (ci*bhs*m_ang - p*ep - q*ep_cc)*f%level(i,j)

      end do x_loop
      end do y_loop
   end subroutine ghp_thorn
!=============================================================================
   pure subroutine ghp_thorn_prime(step, f, f_thorn_prime)
      integer(ip), intent(in)    :: step
      type(field), intent(inout) :: f
      complex(rp), dimension(nx,ny), intent(out)   :: f_thorn_prime

      integer(ip) :: i, j
      real(rp)    :: R

      call set_level(step, f)
      call set_DT(   step, f)

      call compute_DR(f%level, f%DR)

      y_loop: do j=1,ny
      x_loop: do i=1,nx

         R = Rvec(i)

         f_thorn_prime(i,j) = &
            (2.0_rp + (4.0_rp*bhm*r/(cl**2))) * f%DT(i,j) &
         +  ((r/cl)**2)*f%DR(i,j)

      end do x_loop
      end do y_loop
   end subroutine ghp_thorn_prime
!=============================================================================
end module mod_ghp
