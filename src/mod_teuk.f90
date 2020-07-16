!
! Teukolsky field solver
!
module mod_teuk
!-----------------------------------------------------------------------------
   use mod_prec
   use mod_field, only: Field
   use mod_params, only: &
      dt, nx, ny, &
      spin, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin

   use mod_cheb, only: Rvec=>r, set_DR
   use mod_swal, only: Yvec=>y, swal_laplacian

   implicit none
!-----------------------------------------------------------------------------
   private
   public :: Teuk, teuk_constructor 
!-----------------------------------------------------------------------------
   type :: Teuk

   private

   integer(ip) :: m_ang

   real(rp) :: &
      A_pp(nx,ny), A_pq(nx,ny), A_pf(nx,ny), &
      A_qp(nx,ny), A_qq(nx,ny), A_qf(nx,ny), &
      A_fp(nx,ny), A_fq(nx,ny), A_ff(nx,ny)

   complex(rp) :: &
      B_pp(nx,ny), B_pq(nx,ny), B_pf(nx,ny), &
      B_qp(nx,ny), B_qq(nx,ny), B_qf(nx,ny), &
      B_fp(nx,ny), B_fq(nx,ny), B_ff(nx,ny)

   complex(rp) :: &
      p_DR(nx,ny), q_DR(nx,ny), f_DR(nx,ny), f_laplacian(nx,ny)

   contains 

   procedure :: set_k 
   procedure, public :: time_step

   end type Teuk
!-----------------------------------------------------------------------------
   interface Teuk 
      module procedure :: teuk_constructor
   end interface Teuk 
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
   type(Teuk) function teuk_constructor(input_m_ang) result(self)
      integer(ip), intent(in) :: input_m_ang
      integer(ip) :: i, j
      real(rp) :: r, y

      self%m_ang = input_m_ang

      y_loop: do j=1,ny
      x_loop: do i=1,nx
         r = Rvec(i)
         y = Yvec(j)
      !----------------------------
         self%A_pp(i,j) = 0.0_rp

         self%A_pq(i,j) = ((r**2)*((cl**4) - 2*(cl**2)*bhm*r + (bhs**2)*(r**2)))/(cl**4) 

         self%A_pf(i,j) = 0.0_rp
      !----------------------------
         self%A_qp(i,j) = cl**4/(&
            16*cl**2*bhm**2*(cl**2 + 2*bhm*r) &
         +  bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
         )

         self%A_qq(i,j) = (&
            2*(cl**6 + cl**2*(bhs**2 - 8*bhm**2)*r**2 + 4*bhs**2*bhm*r**3) &
         )/(&
            16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
         )

         self%A_qf(i,j) = 0.0_rp
      !----------------------------
         self%A_fp(i,j) = 0.0_rp

         self%A_fq(i,j) = 0.0_rp

         self%A_ff(i,j) = 0.0_rp
      !----------------------------
         self%B_pp(i,j) = 0.0_rp 

         self%B_pq(i,j) = cmplx(&
            -2*r*(-1 - spin + (r*(-2*bhs**2*r + cl**2*bhm*(3 + spin)))/cl**4) &
         ,&
            (-2*bhs*input_m_ang*r**2)/cl**2 &
         ,rp)

         self%B_pf(i,j) = cmplx(&
            (-2*r*(-(bhs**2*r) + cl**2*bhm*(1 + spin)))/cl**4 &
         ,&
            (-2*bhs*input_m_ang*r)/cl**2 &
         ,rp)
      !----------------------------
         self%B_qp(i,j) = cmplx( &
            -( &
            (32*cl**6*bhm**3 - 8*bhs**2*cl**4*bhm*(cl**2 + 4*bhm*r)) &
         /  (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2))**2) &
         , &
            0.0_rp &
         ,rp)

         self%B_qq(i,j) = cmplx( &
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
               -2*bhs*cl**2*(4*input_m_ang*bhm*r + cl**2*(input_m_ang - spin*y)) &
            )/( &
               16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
            ) &
         ,rp)

         self%B_qf(i,j) = cmplx(&
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
                  8*cl**4*bhm**2*(input_m_ang + spin*y) + bhs**2*(input_m_ang*(cl**2 + 4*bhm*r)**2 & 
               -  2*cl**2*(cl**2 + 4*bhm*r)*spin*y + cl**4*input_m_ang*y**2) &
               ) &
            )/( &
               16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
            )**2 &
         ,rp)
      !----------------------------
         self%B_fp(i,j) = cmplx( &
            cl**4 &
         / &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2)) &
         , &
            0.0_rp &
         ,rp)

         self%B_fq(i,j) = cmplx( &
            (2*(cl**6 + cl**2*(bhs**2 - 8*bhm**2)*r**2 + 4*bhs**2*bhm*r**3)) &
         / &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2)) &
         , &
            0.0_rp &
         ,rp)

         self%B_ff(i,j) = cmplx( &
            -( & 
               (2*bhs**2*r*(cl**2 + 6*bhm*r) + 4*cl**2*bhm*(cl**2*spin - 2*bhm*r*(2 + spin))) &
         / &
               (-16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*((cl**2 + 4*bhm*r)**2 - cl**4*y**2)) &
            ) &
         , &
            (-2*bhs*cl**2*(4*input_m_ang*bhm*r + cl**2*(input_m_ang - spin*y))) &
         /  &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2)) &
         ,rp)
      end do x_loop
      end do y_loop

   end function teuk_constructor
!-----------------------------------------------------------------------------
   pure subroutine set_k(self, p, q, f, kp, kq, kf) 
      class(Teuk), intent(inout) :: self
      complex(rp), dimension(nx,ny),  intent(in) :: p, q, f 
      complex(rp), dimension(nx,ny), intent(out) :: kp, kq, kf 

      integer(ip) :: i, j

      call set_DR(p,self%p_DR)
      call set_DR(q,self%q_DR)
      call set_DR(f,self%f_DR)

      call swal_laplacian(spin,self%m_ang,f,self%f_laplacian)

      y_loop: do j=1,ny
      x_loop: do i=1,nx
         kp(i,j) = &
            self%A_pp(i,j) * self%p_DR(i,j) &
         +  self%A_pq(i,j) * self%q_DR(i,j) &
         +  self%A_pf(i,j) * self%f_DR(i,j) &
         +  self%B_pp(i,j) * p(i,j) &
         +  self%B_pq(i,j) * q(i,j) &
         +  self%B_pf(i,j) * f(i,j) &
         +  self%f_laplacian(i,j) 

         kq(i,j) = &
            self%A_qp(i,j) * self%p_DR(i,j) &
         +  self%A_qq(i,j) * self%q_DR(i,j) &
         +  self%A_qf(i,j) * self%f_DR(i,j) &
         +  self%B_qp(i,j) * p(i,j) &
         +  self%B_qq(i,j) * q(i,j) &
         +  self%B_qf(i,j) * f(i,j) 

         kf(i,j) = &
            self%A_fp(i,j) * self%p_DR(i,j) &
         +  self%A_fq(i,j) * self%q_DR(i,j) &
         +  self%A_ff(i,j) * self%f_DR(i,j) &
         +  self%B_fp(i,j) * p(i,j) &
         +  self%B_fq(i,j) * q(i,j) &
         +  self%B_ff(i,j) * f(i,j) 
      end do x_loop
      end do y_loop

   end subroutine set_k
!-----------------------------------------------------------------------------
! RK4 time integrator
!-----------------------------------------------------------------------------
   subroutine time_step(self, p, q, f) 
      class(Teuk), intent(inout) :: self
      type(Field), intent(inout) :: p, q, f 

      integer(ip) :: i, j
   !--------------------------------------------------------
      call self%set_k(p%n, q%n, f%n, p%k1, q%k1, f%k1) 
      do j=1,ny
      do i=1,nx
         p%l2(i,j)= p%n(i,j)+0.5_rp*dt*p%k1(i,j)
         q%l2(i,j)= q%n(i,j)+0.5_rp*dt*q%k1(i,j)
         f%l2(i,j)= f%n(i,j)+0.5_rp*dt*f%k1(i,j)
      end do
      end do
   !--------------------------------------------------------
      call self%set_k(p%l2, q%l2, f%l2, p%k2, q%k2, f%k2) 
      do j=1,ny
      do i=1,nx
         p%l3(i,j)= p%l2(i,j)+0.5_rp*dt*p%k2(i,j)
         q%l3(i,j)= q%l2(i,j)+0.5_rp*dt*q%k2(i,j)
         f%l3(i,j)= f%l2(i,j)+0.5_rp*dt*f%k2(i,j)
      end do
      end do
   !--------------------------------------------------------
      call self%set_k(p%l3, q%l3, f%l3, p%k3, q%k3, f%k3) 
      do j=1,ny
      do i=1,nx
         p%l4(i,j)= p%l3(i,j)+dt*p%k3(i,j)
         q%l4(i,j)= q%l3(i,j)+dt*q%k3(i,j)
         f%l4(i,j)= f%l3(i,j)+dt*f%k3(i,j)
      end do
      end do
   !--------------------------------------------------------
      call self%set_k(p%l4, q%l4, f%l4, p%k4, q%k4, f%k4) 
      do j=1,ny
      do i=1,nx
         p%np1(i,j)= p%n(i,j)+(dt/6.0_rp)*(p%k1(i,j)+2.0_rp*p%k2(i,j)+2.0_rp*p%k3(i,j)+p%k4(i,j))
         q%np1(i,j)= q%n(i,j)+(dt/6.0_rp)*(q%k1(i,j)+2.0_rp*q%k2(i,j)+2.0_rp*q%k3(i,j)+q%k4(i,j))
         f%np1(i,j)= f%n(i,j)+(dt/6.0_rp)*(f%k1(i,j)+2.0_rp*f%k2(i,j)+2.0_rp*f%k3(i,j)+f%k4(i,j))
      end do
      end do
   end subroutine time_step
!-----------------------------------------------------------------------------
end module mod_teuk
