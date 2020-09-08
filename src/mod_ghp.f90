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
   use mod_cheb,   only: R=>Rarr, compute_DR
   use mod_swal,   only: cy=>cyarr, sy=>syarr, swal_lower, swal_raise
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
   subroutine set_edth_arr(m_ang, spin, boost, level, DT, raised, edth_arr)
      integer(ip), intent(in) :: m_ang, spin, boost
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: level, DT
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: raised
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: edth_arr

      real(rp) :: p

      p = (spin+boost)

      edth_arr(:,:,m_ang) = &
         (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2) - ZI*bhs*R*cy)) * ( &
            -  ZI*bhs*sy*(DT(:,:,m_ang)) &
            +  (raised(:,:,m_ang)) &
         ) &
      +  ( &
            (ZI*p*bhs*R*sy/sqrt(2.0_rp)) &
         /  ((ZI*(cl**2) + bhs*R*cy)**2) &
         )*(level(:,:,m_ang))

   end subroutine set_edth_arr
!=============================================================================
   subroutine set_edth_field(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      call set_level(step, m_ang, f)
      call set_DT(   step, m_ang, f)

      call swal_raise(f%spin, m_ang, f%level, f%raised)

      call set_edth_arr(m_ang, f%spin, f%boost, f%level, f%DT, f%raised, f%edth)

   end subroutine set_edth_field
!=============================================================================
! edth_prime rescaled by R
!=============================================================================
   subroutine set_edth_prime_arr(m_ang, spin, boost, level, DT, lowered, edth_prime_arr)
      integer(ip), intent(in) :: m_ang, spin, boost
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: level, DT
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: lowered 
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: edth_prime_arr

      real(rp) :: q

      q = (-spin+boost)

      edth_prime_arr(:,:,m_ang) = &
         (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2) + ZI*bhs*R*cy)) * ( &
               ZI*bhs*sy*(DT(:,:,m_ang)) &
            +  (lowered(:,:,m_ang)) &
         ) &
      +  ( &
            (ZI*q*bhs*R*sy/sqrt(2.0_rp)) &
         /  ((cl**2 + ZI*bhs*R*cy)**2) &
         )*(level(:,:,m_ang))

   end subroutine set_edth_prime_arr
!=============================================================================
   subroutine set_edth_prime_field(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      call set_level(step, m_ang, f)
      call set_DT(   step, m_ang, f)

      call swal_lower(f%spin, m_ang, f%level, f%lowered)

      call set_edth_prime_arr(m_ang, f%spin, f%boost, f%level, f%DT, f%lowered, f%edth_prime)

   end subroutine set_edth_prime_field
!=============================================================================
! thorn rescaled by R
!=============================================================================
   subroutine set_thorn_arr(m_ang, spin, boost, falloff, level, DT, DR, thorn_arr)
      integer(ip), intent(in)    :: m_ang, spin, boost, falloff
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in) :: level, DT, DR
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: thorn_arr

      real(rp) :: p, q

      p = ( spin+boost)
      q = (-spin+boost)

      thorn_arr(:,:,m_ang) = &
         (1.0_rp/((cl**4)+((bhs*R*cy)**2)))*( &
            R*2.0_rp*bhm*(2.0_rp*bhm-((bhs/cl)**2)*R)*(DT(:,:,m_ang)) &
         -  0.5_rp*((cl**2)-(2.0_rp*bhm*R) + ((bhs*R/cl)**2))*( &
               R*DR(:,:,m_ang) &
            +  (falloff)*(level(:,:,m_ang)) &
            ) &
         +  R*(ZI*m_ang*bhs)*(level(:,:,m_ang)) &
         ) &
       - R*(p*ep_0 + q*conjg(ep_0))*(level(:,:,m_ang))

   end subroutine set_thorn_arr
!=============================================================================
   subroutine set_thorn_field(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      call set_level( step, m_ang, f)
      call set_DT(    step, m_ang, f)
      call compute_DR(step, m_ang, f)

      call set_thorn_arr(m_ang, f%spin, f%boost, f%falloff, f%level, f%DT, f%DR, f%thorn)

   end subroutine set_thorn_field
!=============================================================================
! no rescaling in R for thorn prime
!=============================================================================
   subroutine set_thorn_prime_arr(m_ang, falloff, level, DT, DR, thorn_prime_arr)
      integer(ip), intent(in) :: m_ang, falloff
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in) :: level, DT, DR
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: thorn_prime_arr


      thorn_prime_arr(:,:,m_ang) = &
         (2.0_rp + (4.0_rp*bhm*R/(cl**2)))*(DT(:,:,m_ang)) &
      +  ((1.0_rp/cl)**2)*( &
            (R**2)*(DR(:,:,m_ang)) &
         +  R*(falloff)*(level(:,:,m_ang)) &
      )

   end subroutine set_thorn_prime_arr
!=============================================================================
   subroutine set_thorn_prime_field(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      call set_level( step, m_ang, f)
      call set_DT(    step, m_ang, f)
      call compute_DR(step, m_ang, f)

      call set_thorn_prime_arr(m_ang, f%falloff, f%level, f%DT, f%DR, f%thorn_prime)

   end subroutine set_thorn_prime_field
!=============================================================================
end module mod_ghp
