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
   
   public :: ghp_init, set_edth, set_edth_prime, set_thorn, set_thorn_prime
  
   complex(rp), parameter :: ZI = (0.0_rp, 1.0_rp) 

   complex(rp), allocatable :: &
      pre_edth_DT(:,:), &
      pre_edth_raised(:,:), &
      pre_edth_level(:,:), &
      pre_edth_prime_DT(:,:), &
      pre_edth_prime_lowered(:,:), &
      pre_edth_prime_level(:,:), &
      pre_thorn_DT(:,:), &
      pre_thorn_DR(:,:), &
      pre_thorn_level(:,:), &
      pre_thorn_prime_DT(:,:), &
      pre_thorn_prime_DR(:,:) 
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
! precompute some of the prefactors for ghp operators
!=============================================================================
   subroutine ghp_init()
      allocate(pre_edth_DT(nx,ny))
      allocate(pre_edth_raised(nx,ny))
      allocate(pre_edth_level(nx,ny))
      allocate(pre_edth_prime_DT(nx,ny))
      allocate(pre_edth_prime_lowered(nx,ny))
      allocate(pre_edth_prime_level(nx,ny))
      allocate(pre_thorn_DT(nx,ny))
      allocate(pre_thorn_DR(nx,ny))
      allocate(pre_thorn_level(nx,ny))
      allocate(pre_thorn_prime_DT(nx,ny))
      allocate(pre_thorn_prime_DR(nx,ny)) 
   !--------------------------------------------------------------------------
   ! divided by R
   !--------------------------------------------------------------------------
      pre_edth_DT = &
         (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2) - ZI*bhs*R*cy))*(-ZI*bhs*sy)

      pre_edth_raised = &
         (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2) - ZI*bhs*R*cy))

      pre_edth_level = &
         ( &
            (ZI*bhs*R*sy/sqrt(2.0_rp)) &
         /  ((ZI*(cl**2) + bhs*R*cy)**2) &
         )
   !--------------------------------------------------------------------------
   ! divided by R
   !--------------------------------------------------------------------------
      pre_edth_prime_DT = &
         (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2) + ZI*bhs*R*cy)) * ZI*bhs*sy

      pre_edth_prime_lowered = &
         (1.0_rp/sqrt(2.0_rp))*(1.0_rp/((cl**2) + ZI*bhs*R*cy))

      pre_edth_prime_level = &
         ( &
            (ZI*bhs*R*sy/sqrt(2.0_rp)) &
         /  ((cl**2 + ZI*bhs*R*cy)**2) &
         )
   !--------------------------------------------------------------------------
   ! divided by R
   !--------------------------------------------------------------------------
      pre_thorn_DT = &
         (1.0_rp/((cl**4)+((bhs*R*cy)**2))) &
         *R*2.0_rp*bhm*(2.0_rp*bhm-((bhs/cl)**2)*R)

      !--------------
      ! divided by R
      !--------------
      pre_thorn_DR = &
         (1.0_rp/((cl**4)+((bhs*R*cy)**2)))*( &
         -  0.5_rp*((cl**2)-(2.0_rp*bhm*R) + ((bhs*R/cl)**2)) &
         )

      pre_thorn_level = &
         (1.0_rp/((cl**4)+((bhs*R*cy)**2)))*R*(ZI*bhs)
   !--------------------------------------------------------------------------
   ! NOT divided by R
   !--------------------------------------------------------------------------
      pre_thorn_prime_DT = &
         (2.0_rp + (4.0_rp*bhm*R/(cl**2)))

      !--------------
      ! divided by R
      !--------------
      pre_thorn_prime_DR = &
         ((1.0_rp/cl)**2)*R
   !--------------------------------------------------------------------------
   end subroutine ghp_init
!=============================================================================
   subroutine set_edth_arr(m_ang, spin, boost, level, DT, raised, edth_arr)
      integer(ip), intent(in) :: m_ang, spin, boost
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: level, DT
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: raised
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: edth_arr

      real(rp) :: p

      p = (spin+boost)
      
      edth_arr(:,:,m_ang) = &
         pre_edth_DT     *DT(:,:,m_ang) &
      +  pre_edth_raised *raised(:,:,m_ang) &
      +  p*pre_edth_level*level(:,:,m_ang)

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
   subroutine set_edth_prime_arr(m_ang, spin, boost, level, DT, lowered, edth_prime_arr)
      integer(ip), intent(in) :: m_ang, spin, boost
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: level, DT
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: lowered 
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: edth_prime_arr

      real(rp) :: q

      q = (-spin+boost)

      edth_prime_arr(:,:,m_ang) = &
         pre_edth_prime_DT     *DT(:,:,m_ang) &
      +  pre_edth_prime_lowered*lowered(:,:,m_ang) &
      +  q*pre_edth_prime_level*level(:,:,m_ang)

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
   subroutine set_thorn_arr(m_ang, spin, boost, falloff, level, DT, DR, thorn_arr)
      integer(ip), intent(in)    :: m_ang, spin, boost, falloff
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in) :: level, DT, DR
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: thorn_arr

      real(rp) :: p, q

      p = ( spin+boost)
      q = (-spin+boost)

      thorn_arr(:,:,m_ang) = &
         pre_thorn_DT*DT(:,:,m_ang) &
      +  pre_thorn_DR*(R*DR(:,:,m_ang) + falloff*level(:,:,m_ang)) &
      +  pre_thorn_level*m_ang     *level(:,:,m_ang) &
      -  R*(p*ep_0 + q*conjg(ep_0))*level(:,:,m_ang)

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
         pre_thorn_prime_DT*DT(:,:,m_ang) &
      +  pre_thorn_prime_DR*(R*DR(:,:,m_ang) + falloff*level(:,:,m_ang))

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
