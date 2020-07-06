!
! Provides the Field class and its methods.
!
module mod_field
!-----------------------------------------------------------------------------
  use mod_def_prec

  implicit none
!-----------------------------------------------------------------------------
  private
  public :: Field 
!-----------------------------------------------------------------------------
  type :: Field

  character(:), allocatable :: name
  integer(ip)               :: dims(2)
  ! always have two levels: n, np1 and intermediate levels for RK4 time
  ! evolution. k1, etc. are the time derivatives (note 1==n, 5==np1).
  ! These functions are complex as we work in NP formalism
  complex(rp), allocatable :: &
    n( :,:),l2(:,:),l3(:,:),l4(:,:),np1(:,:), &
    k1(:,:),k2(:,:),k3(:,:),k4(:,:),k5(:,:)

  end type Field
!-----------------------------------------------------------------------------
  interface Field
    module procedure :: field_constructor
  end interface Field
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
  type(Field) function field_constructor(name, dims) result(self)
    character(*), intent(in) :: name ! field name
    integer(ip),  intent(in) :: dims(2) ! domain size in x and y
    self % name = name 
    self % dims = dims
    allocate(self % n(  dims(1),dims(2)))
    allocate(self % l2( dims(1),dims(2)))
    allocate(self % l3( dims(1),dims(2)))
    allocate(self % l4( dims(1),dims(2)))
    allocate(self % np1(dims(1),dims(2)))
    allocate(self % k1( dims(1),dims(2)))
    allocate(self % k2( dims(1),dims(2)))
    allocate(self % k3( dims(1),dims(2)))
    allocate(self % k4( dims(1),dims(2)))
    allocate(self % k5( dims(1),dims(2)))
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
  end function field_constructor
!-----------------------------------------------------------------------------
  pure subroutine from_field(target, source)
    ! Initializes Field instance target using components
    ! from Field instance source. Used to initialize a 
    ! Field from another Field without invoking the 
    ! assignment operator.
    type(Field), intent(in out) :: target
    type(Field), intent(in) :: source
    target % name = source % name
    target % dims = source % dims

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
  end subroutine from_field
!-----------------------------------------------------------------------------
end module mod_field
