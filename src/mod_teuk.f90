!
! Teukolsky field solver
!
module mod_teuk
!-----------------------------------------------------------------------------
   use mod_prec
   use mod_field
   use mod_params, only: &
      nx, ny, &
      spin, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin

   use mod_cheb, only: Rvec=>R
   use mod_swal, only: Yvec=>Y

   implicit none
!-----------------------------------------------------------------------------
   private
   public :: Teuk, teuk_constructor 
!-----------------------------------------------------------------------------
   type :: Teuk

   private

   real(rp), allocatable :: &
      A_pp(:,:), A_pq(:,:), A_pf(:,:), &
      A_qp(:,:), A_qq(:,:), A_qf(:,:), &
      A_fp(:,:), A_fq(:,:), A_ff(:,:)

   complex(rp), allocatable :: &
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
   type(Teuk) function teuk_constructor(m_ang) result(self)
      integer(ip), intent(in) :: m_ang
      integer(ip) :: i, j
      real(rp) :: r, y

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
         r = Rvec(i)
         y = Yvec(j)
      !----------------------------
         self % A_pp(i,j) = 0.0_rp

         self % A_pq(i,j) = ((r**2)*((cl**4) - 2*(cl**2)*bhm*r + (bhs**2)*(r**2)))/(cl**4) 

         self % A_pf(i,j) = 0.0_rp
      !----------------------------
         self % A_qp(i,j) = cl**4/(&
            16*cl**2*bhm**2*(cl**2 + 2*bhm*r) &
         +  bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
         )

         self % A_qq(i,j) = (&
            2*(cl**6 + cl**2*(bhs**2 - 8*bhm**2)*r**2 + 4*bhs**2*bhm*r**3) &
         )/(&
            16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
         )

         self % A_qf(i,j) = 0.0_rp
      !----------------------------
         self % A_fp(i,j) = 0.0_rp

         self % A_fq(i,j) = 0.0_rp

         self % A_ff(i,j) = 0.0_rp
      !----------------------------
         self % B_pp(i,j) = 0.0_rp 

         self % B_pq(i,j) = cmplx(&
            -2*r*(-1 - spin + (r*(-2*bhs**2*r + cl**2*bhm*(3 + spin)))/cl**4) &
         ,&
            (-2*bhs*m_ang*r**2)/cl**2 &
         ,rp)

         self % B_pf(i,j) = cmplx(&
            (-2*r*(-(bhs**2*r) + cl**2*bhm*(1 + spin)))/cl**4 &
         ,&
            (-2*bhs*m_ang*r)/cl**2 &
         ,rp)
      !----------------------------
         self % B_qp(i,j) = cmplx( &
            -( &
            (32*cl**6*bhm**3 - 8*bhs**2*cl**4*bhm*(cl**2 + 4*bhm*r)) &
         /  (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2))**2) &
         , &
            0.0_rp &
         ,rp)

         self % B_qq(i,j) = cmplx( &
            -( &
            ( &
               64*cl**4*bhm**3*(12*cl**2*bhm*r - cl**4*(-1 + spin) + 4*bhm**2*r**2*(4 + spin)) &
            +  2*bhs**4*r*((cl**2 + 4*bhm*r)**2*(3*cl**2 + 10*bhm*r) - 3*cl**4*(cl**2 + 6*bhm*r)*y**2) &
            -  4*bhs**2*cl**2*bhm*( &
                  240*cl**2*bhm**2*r**2 + 32*bhm**3*r**3*(9 + spin) &
               -  2*cl**4*bhm*r*(-26 + 3*spin + (6 + spin)*y**2) + cl**6*(4 + spin*(-1 + y**2)) &
            ) &
            )/( &
            16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2))**2 &
            ) &
         , &
            ( &
               -2*bhs*cl**2*(4*m_ang*bhm*r + cl**2*(m_ang - spin*y)) &
            )/( &
               16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
            ) &
         ,rp)

         self % B_qf(i,j) = cmplx(&
            ( &
            -  2*cl**2*( &
                  128*cl**4*bhm**4*(1 + spin) &
               +  bhs**4*(32*bhm**2*r**2 - cl**4*(-1 + y**2) -  12*cl**2*bhm*r*(-1 + y**2)) &
               +  4*bhs**2*bhm**2*( &
                     16*bhm**2*r**2*(-1 + spin) - 16*cl**2*bhm*r*(3 + spin) + cl**4*(-6 - 5*spin + (2 + spin)*y**2) &
                  ) &
            ) &
            )/( &
               16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
            )**2 &
         ,&
            ( &
            -  8*bhs*cl**2*bhm*( &
                  8*cl**4*bhm**2*(m_ang + spin*y) + bhs**2*(m_ang*(cl**2 + 4*bhm*r)**2 & 
               -  2*cl**2*(cl**2 + 4*bhm*r)*spin*y + cl**4*m_ang*y**2) &
               ) &
            )/( &
               16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
            )**2 &
         ,rp)
      !----------------------------
         self % B_fp(i,j) = cmplx( &
            cl**4 &
         / &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2)) &
         , &
            0.0_rp &
         ,rp)

         self % B_fq(i,j) = cmplx( &
            (2*(cl**6 + cl**2*(bhs**2 - 8*bhm**2)*r**2 + 4*bhs**2*bhm*r**3)) &
         / &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2)) &
         , &
            0.0_rp &
         ,rp)

         self % B_ff(i,j) = cmplx( &
            -( & 
               (2*bhs**2*r*(cl**2 + 6*bhm*r) + 4*cl**2*bhm*(cl**2*spin - 2*bhm*r*(2 + spin))) &
         / &
               (-16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*((cl**2 + 4*bhm*r)**2 - cl**4*y**2)) &
            ) &
         , &
            (-2*bhs*cl**2*(4*m_ang*bhm*r + cl**2*(m_ang - spin*y))) &
         /  &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2)) &
         ,rp)
      end do x_loop
      end do y_loop

   end function teuk_constructor
!-----------------------------------------------------------------------------
end module mod_teuk
