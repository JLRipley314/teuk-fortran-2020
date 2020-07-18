!
! Provides the field type, and some basic utility functions 
!
module mod_field
!=============================================================================
   use mod_prec
   use mod_params, only: nx, ny, lmax
   implicit none
!=============================================================================
   private
   public :: field, set_field, set_level, set_DT, shift_time_step
!=============================================================================
   type :: field

   character(:), allocatable :: name
   ! always have two levels: n, np1 and intermediate levels for RK4 time
   ! evolution. k1, etc. are the time derivatives (note 1==n, 5==np1).
   ! These functions are complex as we work in NP formalism

   integer(ip) :: m_ang, spin, boost, falloff

   complex(rp) :: &
   n( nx,ny), l2(nx,ny), l3(nx,ny), l4(nx,ny), np1(nx,ny), &
   k1(nx,ny), k2(nx,ny), k3(nx,ny), k4(nx,ny), k5( nx,ny), &

   level(nx,ny), &
   DT(nx,ny), DR(nx,ny), &
   raised(nx,ny), lowered(nx,ny), lap(nx,ny), coefs(nx,0:lmax)

   end type field

   interface set_level
      module procedure set_level_interior, set_level_exterior 
   end interface set_level
!=============================================================================
contains
!=============================================================================
   pure subroutine set_field(name, m_ang, spin, boost, falloff, f)
      character(*), intent(in)  :: name ! field name
      integer(ip),  intent(in)  :: m_ang, spin, boost, falloff
      type(field),  intent(out) :: f

      f % name = name 

      f % m_ang = m_ang
      f % spin  = spin
      f % boost = boost
      f % falloff = falloff

      f % n   = 0.0_rp
      f % l2  = 0.0_rp
      f % l3  = 0.0_rp
      f % l4  = 0.0_rp
      f % np1 = 0.0_rp

      f % k1 = 0.0_rp
      f % k2 = 0.0_rp
      f % k3 = 0.0_rp
      f % k4 = 0.0_rp
      f % k5 = 0.0_rp

      f % level = 0.0_rp
      f % DT = 0.0_rp
      f % DR = 0.0_rp

      f % raised  = 0.0_rp
      f % lowered = 0.0_rp
      f % lap     = 0.0_rp

      f % coefs = 0.0_rp

   end subroutine set_field
!=============================================================================
   pure subroutine set_level_interior(step, f)
      integer(ip), intent(in)    :: step
      type(field), intent(inout) :: f

      select case (step)
         case (1)
            f % level = f % n 
         case (2)
            f % level = f % l2 
         case (3)
            f % level = f % l3 
         case (4)
            f % level = f % l4 
         case (5)
            f % level = f % np1 
         case default
            f % level = -1.0_rp
      end select
   end subroutine set_level_interior
!=============================================================================
   pure subroutine set_level_exterior(step, f, ex)
      integer(ip), intent(in)  :: step
      type(field), intent(in)  :: f
      complex(rp), dimension(nx,ny), intent(out) :: ex

      select case (step)
         case (1)
            ex = f % n 
         case (2)
            ex = f % l2 
         case (3)
            ex = f % l3 
         case (4)
            ex = f % l4 
         case (5)
            ex = f % np1 
         case default
            ex = -1.0_rp
      end select
   end subroutine set_level_exterior
!=============================================================================
   pure subroutine set_DT(step, f)
      integer(ip), intent(in)    :: step
      type(field), intent(inout) :: f

      select case (step)
         case (1)
            f % DT = f % k1 
         case (2)
            f % DT = f % k2 
         case (3)
            f % DT = f % k3 
         case (4)
            f % DT = f % k4 
         case (5)
            f % DT = f % k5
         case default
            f % DT = -1.0_rp
      end select
   end subroutine set_DT
!=============================================================================
   pure subroutine shift_time_step(f)
      type(field), intent(inout) :: f 

      f % n  = f % np1
      f % k1 = f % k5
   end subroutine shift_time_step
!=============================================================================
end module mod_field
