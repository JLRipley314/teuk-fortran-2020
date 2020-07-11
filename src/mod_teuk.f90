!
! Teukolsky field solver
!
module mod_teuk
!-----------------------------------------------------------------------------
   use mod_prec
   use mod_field
   use mod_sim_params

   implicit none
!-----------------------------------------------------------------------------
   private
   public :: Teuk 
!-----------------------------------------------------------------------------
   type :: Teuk

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
   type(Teuk) function teuk_constructor() result(self)
      integer(ip) :: i, j

      allocate(self % A_pp(nx,ny))
      allocate(self % A_pq(nx,ny))
      allocate(self % A_pf(nx,ny))

      allocate(self % A_qp(nx,ny))
      allocate(self % A_qq(nx,ny))
      allocate(self % A_qf(nx,ny))

      allocate(self % A_fp(nx,ny))
      allocate(self % A_fq(nx,ny))
      allocate(self % A_ff(nx,ny))

      allocate(self % B_pp(nx,ny))
      allocate(self % B_pq(nx,ny))
      allocate(self % B_pf(nx,ny))

      allocate(self % B_qp(nx,ny))
      allocate(self % B_qq(nx,ny))
      allocate(self % B_qf(nx,ny))

      allocate(self % B_fp(nx,ny))
      allocate(self % B_fq(nx,ny))
      allocate(self % B_ff(nx,ny))

      y_loop: do j=1,ny
      x_loop: do i=1,nx
         self % A_pp(i,j) = cmplx(0,0,rp)
         self % A_pq(i,j) = cmplx(0,0,rp)
         self % A_pf(i,j) = cmplx(0,0,rp)

         self % A_qp(i,j) = cmplx(0,0,rp)
         self % A_qq(i,j) = cmplx(0,0,rp)
         self % A_qf(i,j) = cmplx(0,0,rp)

         self % A_fp(i,j) = cmplx(0,0,rp)
         self % A_fq(i,j) = cmplx(0,0,rp)
         self % A_ff(i,j) = cmplx(0,0,rp)

         self % B_pp(i,j) = cmplx(0,0,rp)
         self % B_pq(i,j) = cmplx(0,0,rp)
         self % B_pf(i,j) = cmplx(0,0,rp)

         self % B_qp(i,j) = cmplx(0,0,rp)
         self % B_qq(i,j) = cmplx(0,0,rp)
         self % B_qf(i,j) = cmplx(0,0,rp)

         self % B_fp(i,j) = cmplx(0,0,rp)
         self % B_fq(i,j) = cmplx(0,0,rp)
         self % B_ff(i,j) = cmplx(0,0,rp)
      end do x_loop
      end do y_loop

   end function teuk_constructor
!-----------------------------------------------------------------------------
end module mod_teuk
