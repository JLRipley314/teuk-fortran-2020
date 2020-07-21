!
! Geroch-Held-Penrose (GHP) operators
!
module mod_ghp
!=============================================================================
   use mod_prec
   use mod_cheb,   only: R, compute_DR
   use mod_swal,   only: cy, sy, swal_lower, swal_raise
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

   public :: set_edth, set_edth_prime, set_thorn, set_thorn_prime
!=============================================================================
   contains
!=============================================================================
! edth rescaled by R
!=============================================================================
   pure subroutine set_edth(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      integer(ip) :: i, j
      real(rp) :: p

      p = (f%spin+f%boost)

      call set_level(step, m_ang, f)
      call set_DT(   step, m_ang, f)

      call swal_raise(f%spin, m_ang, f%level, f%coefs, f%raised)

      do j=1,ny
      do i=1,nx
         f%edth(i,j,m_ang) = &
            (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2)-ci*bhs*R(i)*cy(j))) * ( &
               - ci*bhs*sy(j)*(f%DT(i,j,m_ang)) &
               + (f%raised(i,j,m_ang)) &
            ) &
         -  ( &
               (ci*p*bhs*R(i)*sy(j)/sqrt(2.0_rp)) &
            /  (((cl**2) - ci*bhs*R(i)*cy(j))**2) &
            )*(f%level(i,j,m_ang))
      end do 
      end do 
   end subroutine set_edth
!=============================================================================
! edth_prime rescaled by R
!=============================================================================
   pure subroutine set_edth_prime(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      integer(ip) :: i, j
      real(rp) :: q

      q = (-f%spin+f%boost)

      call set_level(step, m_ang, f)
      call set_DT(   step, m_ang, f)

      call swal_lower(f%spin, m_ang, f%level, f%coefs, f%lowered)

      do j=1,ny
      do i=1,nx
         f%edth_prime(i,j,m_ang) = &
            (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2)+ci*bhs*R(i)*cy(j))) * ( &
                 ci*bhs*sy(j)*(f%DT(i,j,m_ang)) &
               + (f%lowered(i,j,m_ang)) &
            ) &
         +  ( &
               (ci*q*bhs*R(i)*sy(j)/sqrt(2.0_rp)) &
            /  (((cl**2) + ci*bhs*R(i)*cy(j))**2) &
            )*(f%level(i,j,m_ang))

      end do
      end do
   end subroutine set_edth_prime
!=============================================================================
! thorn rescaled by R
!=============================================================================
   pure subroutine set_thorn(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      integer(ip) :: i, j
      real(rp)    :: p, q

      p = ( f%spin+f%boost)
      q = (-f%spin+f%boost)

      call set_level(step, m_ang, f)
      call set_DT(   step, m_ang, f)

      call compute_DR(m_ang, f%level, f%DR)

      do j=1,ny
      do i=1,nx
         f%thorn(i,j,m_ang) = &
            (1.0_rp/((cl**4)+((bhs*R(i)*cy(j))**2)))*( &
               R(i)*2.0_rp*bhm*(2.0_rp*bhm-((bhs/cl)**2)*R(i))*(f%DT(i,j,m_ang)) &
            -  0.5_rp*((cl**2)-(2.0_rp*bhm*R(i)) + ((bhs*R(i)/cl)**2))*( &
                  R(i)*f%DR(i,j,m_ang) &
               +  (f%falloff)*(f%level(i,j,m_ang)) &
               ) &
            +  R(i)*(ci*bhs*m_ang)*(f%level(i,j,m_ang)) &
            ) &
          - R(i)*(p*ep_0(i,j) + q*conjg(ep_0(i,j)))*(f%level(i,j,m_ang))
      end do
      end do
   end subroutine set_thorn
!=============================================================================
! no rescaling in R for thorn prime
!=============================================================================
   pure subroutine set_thorn_prime(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      integer(ip) :: i, j

      call set_level(step, m_ang, f)
      call set_DT(   step, m_ang, f)

      call compute_DR(m_ang, f%level, f%DR)

      do j=1,ny
      do i=1,nx
         f%thorn_prime(i,j,m_ang) = &
            (2.0_rp + (4.0_rp*bhm*R(i)/(cl**2)))*(f%DT(i,j,m_ang)) &
         +  ((1.0_rp/cl)**2)*( &
               (R(i)**2)*(f%DR(i,j,m_ang)) &
            +  R(i)*(f%falloff)*(f%level(i,j,m_ang)) &
         )
      end do
      end do
   end subroutine set_thorn_prime
!=============================================================================
end module mod_ghp
