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
   public :: field, shift_time_step, copy_field
!=============================================================================
   type :: field

   character(:), allocatable :: name
   ! always have two levels: n, np1 and intermediate levels for RK4 time
   ! evolution. k1, etc. are the time derivatives (note 1==n, 5==np1).
   ! These functions are complex as we work in NP formalism
   complex(rp) :: &
   n( nx,ny), l2(nx,ny), l3(nx,ny), l4(nx,ny), np1(nx,ny), &
   k1(nx,ny), k2(nx,ny), k3(nx,ny), k4(nx,ny), k5( nx,ny), &
   DR(nx,ny),lap(nx,ny), coefs(nx,0:lmax)

   end type field
!=============================================================================
   interface field
      module procedure :: field_constructor
   end interface field
!=============================================================================
contains
!=============================================================================
   type(field) function field_constructor(name) result(self)
      character(*), intent(in) :: name ! field name
      self % name = name 

      self % n   = 0.0_rp
      self % l2  = 0.0_rp 
      self % l3  = 0.0_rp
      self % l4  = 0.0_rp
      self % np1 = 0.0_rp
      self % k1  = 0.0_rp
      self % k2  = 0.0_rp
      self % k3  = 0.0_rp
      self % k4  = 0.0_rp
      self % k5  = 0.0_rp

      self % DR  = 0.0_rp
      self % lap = 0.0_rp
      self % coefs = 0.0_rp
   end function field_constructor
!=============================================================================
   pure subroutine copy_field(target, source)
      ! Initializes field instance target using components
      ! from field instance source. Used to initialize a 
      ! field from another field without invoking the 
      ! assignment operator.
      type(field), intent(in out) :: target
      type(field), intent(in) :: source
      target % name = source % name

      target % n   = source % n
      target % l2  = source % l2
      target % l3  = source % l3
      target % l4  = source % l4
      target % np1 = source % np1

      target % k1 = source % k1
      target % k2 = source % k2
      target % k3 = source % k3
      target % k4 = source % k4
      target % k5 = source % k5

      target % DR    = source % DR
      target % lap   = source % lap
      target % coefs = source % coefs
   end subroutine copy_field
!=============================================================================
   pure subroutine shift_time_step(f)
      type(field), intent(inout) :: f 

      f % n  = f % np1
      f % k1 = f % k5
   end subroutine shift_time_step
!=============================================================================
end module mod_field
