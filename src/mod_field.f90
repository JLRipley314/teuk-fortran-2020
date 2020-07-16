!
! Provides the Field class and its methods.
!
module mod_field
!=============================================================================
   use mod_prec

   use mod_params, only: nx, ny
   implicit none
!=============================================================================
   private
   public :: Field, shift_time_step
!=============================================================================
   type :: Field

   character(:), allocatable :: name
   ! always have two levels: n, np1 and intermediate levels for RK4 time
   ! evolution. k1, etc. are the time derivatives (note 1==n, 5==np1).
   ! These functions are complex as we work in NP formalism
   complex(rp) :: &
   n( nx,ny),l2(nx,ny),l3(nx,ny),l4(nx,ny),np1(nx,ny), &
   k1(nx,ny),k2(nx,ny),k3(nx,ny),k4(nx,ny),k5(nx,ny), &
   DR(nx,ny), lap(nx,ny)

   end type Field
!=============================================================================
   interface Field
      module procedure :: field_constructor
   end interface Field
!=============================================================================
contains
!=============================================================================
   type(Field) function field_constructor(name) result(self)
      character(*), intent(in) :: name ! field name
      self % name = name 
      self % n   = 0
      self % l2  = 0
      self % l3  = 0
      self % l4  = 0
      self % np1 = 0
      self % k1  = 0
      self % k2  = 0
      self % k3  = 0
      self % k4  = 0
      self % k5  = 0

      self % DR  = 0
      self % lap = 0
   end function field_constructor
!=============================================================================
   pure subroutine from_field(target, source)
      ! Initializes Field instance target using components
      ! from Field instance source. Used to initialize a 
      ! Field from another Field without invoking the 
      ! assignment operator.
      type(Field), intent(in out) :: target
      type(Field), intent(in) :: source
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

      target % DR  = source % DR
      target % lap = source % lap
   end subroutine from_field
!=============================================================================
   pure subroutine shift_time_step(f)
      type(Field), intent(in out) :: f 

      f % n  = f % np1
      f % k1 = f % k5
   end subroutine shift_time_step
!=============================================================================
end module mod_field
