!
! Geroch-Held-Penrose (GHP) operators
!
! edth:        rescaled by R
! edth_prime:  rescaled by R
! thorn:       rescaled by R
! thorn_prime: not rescaled (is already nonzero at future null infinity)
!
module mod_ghp
!=============================================================================
   use mod_prec
   use mod_cheb,   only: R=>Rvec, compute_DR
   use mod_swal,   only: cy=>cyvec, sy=>syvec, swal_lower, swal_raise
   use mod_field,  only: field, set_level, set_DT

   use mod_bkgrd_np, only: ep_0
   use mod_params,   only: &
      dt, nx, ny, min_m, max_m, max_l, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin
!=============================================================================
   implicit none
   private
   
   public :: set_edth, set_edth_prime, set_thorn, set_thorn_prime
  
   complex(rp), parameter :: ZI = (0.0_rp, 1.0_rp) 
!=============================================================================
   interface set_edth
      module procedure set_edth_field, set_edth_arr
   end interface

   interface set_edth_prime
      module procedure set_edth_prime_field, set_edth_prime_arr
   end interface

   interface set_thorn
      module procedure set_thorn_field, set_thorn_arr
   end interface

   interface set_thorn_prime
      module procedure set_thorn_prime_field, set_thorn_prime_arr
   end interface
!=============================================================================
   contains
!=============================================================================
! edth rescaled by R
!=============================================================================
   pure subroutine set_edth_field(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      integer(ip) :: i, j
      real(rp) :: p

      p = (f%spin+f%boost)

      call set_level(step, m_ang, f)
      call set_DT(   step, m_ang, f)

      call swal_raise(f%spin, m_ang, f%level, f%coefs_swal, f%raised)

      do j=1,ny
      do i=1,nx
         f%edth(i,j,m_ang) = &
            (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2) - ZI*bhs*R(i)*cy(j))) * ( &
               -  ZI*bhs*sy(j)*(f%DT(i,j,m_ang)) &
               +  (f%raised(i,j,m_ang)) &
            ) &
         +  ( &
               (ZI*p*bhs*R(i)*sy(j)/sqrt(2.0_rp)) &
            /  ((ZI*(cl**2) + bhs*R(i)*cy(j))**2) &
            )*(f%level(i,j,m_ang))
      end do 
      end do 
   end subroutine set_edth_field
!=============================================================================
   pure subroutine set_edth_arr(m_ang, spin, boost, level, DT, raised, edth_arr)
      integer(ip), intent(in) :: m_ang, spin, boost
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: level, DT
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: raised
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: edth_arr

      integer(ip) :: i, j
      real(rp) :: p

      p = (spin+boost)

      do j=1,ny
      do i=1,nx
         edth_arr(i,j,m_ang) = &
            (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2) - ZI*bhs*R(i)*cy(j))) * ( &
               -  ZI*bhs*sy(j)*(DT(i,j,m_ang)) &
               +  (raised(i,j,m_ang)) &
            ) &
         +  ( &
               (ZI*p*bhs*R(i)*sy(j)/sqrt(2.0_rp)) &
            /  ((ZI*(cl**2) + bhs*R(i)*cy(j))**2) &
            )*(level(i,j,m_ang))
      end do 
      end do 
   end subroutine set_edth_arr
!=============================================================================
! edth_prime rescaled by R
!=============================================================================
   pure subroutine set_edth_prime_field(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      integer(ip) :: i, j
      real(rp) :: q

      q = (-f%spin+f%boost)

      call set_level(step, m_ang, f)
      call set_DT(   step, m_ang, f)

      call swal_lower(f%spin, m_ang, f%level, f%coefs_swal, f%lowered)

      do j=1,ny
      do i=1,nx
         f%edth_prime(i,j,m_ang) = &
            (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2) + ZI*bhs*R(i)*cy(j))) * ( &
                  ZI*bhs*sy(j)*(f%DT(i,j,m_ang)) &
               +  (f%lowered(i,j,m_ang)) &
            ) &
         +  ( &
               (ZI*q*bhs*R(i)*sy(j)/sqrt(2.0_rp)) &
            /  ((cl**2 + ZI*bhs*R(i)*cy(j))**2) &
            )*(f%level(i,j,m_ang))

      end do
      end do
   end subroutine set_edth_prime_field
!=============================================================================
   pure subroutine set_edth_prime_arr(m_ang, spin, boost, level, DT, lowered, edth_prime_arr)
      integer(ip), intent(in) :: m_ang, spin, boost
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: level, DT
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: lowered 
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: edth_prime_arr

      integer(ip) :: i, j
      real(rp) :: q

      q = (-spin+boost)

      do j=1,ny
      do i=1,nx
         edth_prime_arr(i,j,m_ang) = &
            (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2) + ZI*bhs*R(i)*cy(j))) * ( &
                  ZI*bhs*sy(j)*(DT(i,j,m_ang)) &
               +  (lowered(i,j,m_ang)) &
            ) &
         +  ( &
               (ZI*q*bhs*R(i)*sy(j)/sqrt(2.0_rp)) &
            /  ((cl**2 + ZI*bhs*R(i)*cy(j))**2) &
            )*(level(i,j,m_ang))

      end do
      end do
   end subroutine set_edth_prime_arr
!=============================================================================
! thorn rescaled by R
!=============================================================================
   pure subroutine set_thorn_field(step, m_ang, f)
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
            +  R(i)*(ZI*m_ang*bhs)*(f%level(i,j,m_ang)) &
            ) &
          - R(i)*(p*ep_0(i,j) + q*conjg(ep_0(i,j)))*(f%level(i,j,m_ang))
      end do
      end do
   end subroutine set_thorn_field
!=============================================================================
   pure subroutine set_thorn_arr(m_ang, spin, boost, falloff, level, DT, DR, thorn_arr)
      integer(ip), intent(in)    :: m_ang, spin, boost, falloff
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in) :: level, DT, DR
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: thorn_arr

      integer(ip) :: i, j
      real(rp)    :: p, q

      p = ( spin+boost)
      q = (-spin+boost)

      do j=1,ny
      do i=1,nx
         thorn_arr(i,j,m_ang) = &
            (1.0_rp/((cl**4)+((bhs*R(i)*cy(j))**2)))*( &
               R(i)*2.0_rp*bhm*(2.0_rp*bhm-((bhs/cl)**2)*R(i))*(DT(i,j,m_ang)) &
            -  0.5_rp*((cl**2)-(2.0_rp*bhm*R(i)) + ((bhs*R(i)/cl)**2))*( &
                  R(i)*DR(i,j,m_ang) &
               +  (falloff)*(level(i,j,m_ang)) &
               ) &
            +  R(i)*(ZI*m_ang*bhs)*(level(i,j,m_ang)) &
            ) &
          - R(i)*(p*ep_0(i,j) + q*conjg(ep_0(i,j)))*(level(i,j,m_ang))
      end do
      end do
   end subroutine set_thorn_arr
!=============================================================================
! no rescaling in R for thorn prime
!=============================================================================
   pure subroutine set_thorn_prime_field(step, m_ang, f)
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
   end subroutine set_thorn_prime_field
!=============================================================================
   pure subroutine set_thorn_prime_arr(m_ang, falloff, level, DT, DR, thorn_prime_arr)
      integer(ip), intent(in) :: m_ang, falloff
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in) :: level, DT, DR
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: thorn_prime_arr

      integer(ip) :: i, j

      do j=1,ny
      do i=1,nx
         thorn_prime_arr(i,j,m_ang) = &
            (2.0_rp + (4.0_rp*bhm*R(i)/(cl**2)))*(DT(i,j,m_ang)) &
         +  ((1.0_rp/cl)**2)*( &
               (R(i)**2)*(DR(i,j,m_ang)) &
            +  R(i)*(falloff)*(level(i,j,m_ang)) &
         )
      end do
      end do
   end subroutine set_thorn_prime_arr
!=============================================================================
end module mod_ghp
