!
! Teukolsky field solver
!
module mod_teuk
!-----------------------------------------------------------------------------
  use mod_prec

  use mod_field

  implicit none
!-----------------------------------------------------------------------------
  private
  public :: Teuk 
!-----------------------------------------------------------------------------
  type :: Teuk
  
  integer(ip) :: dims(2)
  complex(rp), allocatable :: &
    A_pp(:,:), A_pq(:,:), A_pf(:,:), &
    A_qp(:,:), A_qq(:,:), A_qf(:,:), &
    A_fp(:,:), A_fq(:,:), A_ff(:,:), &
    B_pp(:,:), B_pq(:,:), B_pf(:,:), &
    B_qp(:,:), B_qq(:,:), B_qf(:,:), &
    B_fp(:,:), B_fq(:,:), B_ff(:,:)

  end type Teuk
!-----------------------------------------------------------------------------
  interface Teuk 
    module procedure :: teuk_constructor
  end interface Teuk 
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
  type(Teuk) function teuk_constructor(dims) result(self)
    integer(ip), intent(in) :: dims(2) ! domain size in x and y

    integer(ip) :: i, j

    allocate(self % A_pp(dims(1),dims(2)))
    allocate(self % A_pq(dims(1),dims(2)))
    allocate(self % A_pf(dims(1),dims(2)))

    allocate(self % A_qp(dims(1),dims(2)))
    allocate(self % A_qq(dims(1),dims(2)))
    allocate(self % A_qf(dims(1),dims(2)))

    allocate(self % A_fp(dims(1),dims(2)))
    allocate(self % A_fq(dims(1),dims(2)))
    allocate(self % A_ff(dims(1),dims(2)))

    allocate(self % B_pp(dims(1),dims(2)))
    allocate(self % B_pq(dims(1),dims(2)))
    allocate(self % B_pf(dims(1),dims(2)))

    allocate(self % B_qp(dims(1),dims(2)))
    allocate(self % B_qq(dims(1),dims(2)))
    allocate(self % B_qf(dims(1),dims(2)))

    allocate(self % B_fp(dims(1),dims(2)))
    allocate(self % B_fq(dims(1),dims(2)))
    allocate(self % B_ff(dims(1),dims(2)))

    y_loop: do j=1,dims(2)
      x_loop: do i=1,dims(1)
        self % A_pp(i,j) = (0,0)
        self % A_pq(i,j) = (0,0)
        self % A_pf(i,j) = (0,0)

        self % A_qp(i,j) = (0,0)
        self % A_qq(i,j) = (0,0)
        self % A_qf(i,j) = (0,0)

        self % A_fp(i,j) = (0,0)
        self % A_fq(i,j) = (0,0)
        self % A_ff(i,j) = (0,0)

        self % B_pp(i,j) = (0,0)
        self % B_pq(i,j) = (0,0)
        self % B_pf(i,j) = (0,0)

        self % B_qp(i,j) = (0,0)
        self % B_qq(i,j) = (0,0)
        self % B_qf(i,j) = (0,0)

        self % B_fp(i,j) = (0,0)
        self % B_fq(i,j) = (0,0)
        self % B_ff(i,j) = (0,0)
      end do x_loop
    end do y_loop

  end function teuk_constructor
!-----------------------------------------------------------------------------
end module mod_teuk
