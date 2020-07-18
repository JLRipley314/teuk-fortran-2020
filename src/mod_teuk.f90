!
! Teukolsky field solver
!
module mod_teuk
!=============================================================================
   use mod_prec
   use mod_field, only: field
   use mod_params, only: &
      dt, nx, ny, lmax, &
      spin, min_m, max_m, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin

   use mod_cheb, only: Rvec=>r, compute_DR
   use mod_swal, only: Yvec=>y, swal_laplacian
!=============================================================================
   implicit none
   private
   public :: teuk_init, teuk_time_step, compute_q_indep_res 

   type(field), public :: psi4_f, psi4_p, psi4_q
!=============================================================================
   real(rp) :: &
      A_pp(nx,ny,min_m:max_m), A_pq(nx,ny,min_m:max_m), A_pf(nx,ny,min_m:max_m), &
      A_qp(nx,ny,min_m:max_m), A_qq(nx,ny,min_m:max_m), A_qf(nx,ny,min_m:max_m), &
      A_fp(nx,ny,min_m:max_m), A_fq(nx,ny,min_m:max_m), A_ff(nx,ny,min_m:max_m)

   complex(rp) :: &
      B_pp(nx,ny,min_m:max_m), B_pq(nx,ny,min_m:max_m), B_pf(nx,ny,min_m:max_m), &
      B_qp(nx,ny,min_m:max_m), B_qq(nx,ny,min_m:max_m), B_qf(nx,ny,min_m:max_m), &
      B_fp(nx,ny,min_m:max_m), B_fq(nx,ny,min_m:max_m), B_ff(nx,ny,min_m:max_m)
!=============================================================================
contains
!=============================================================================
   subroutine teuk_init()
      integer(ip) :: i, j, m_ang
      real(rp) :: r, y

      m_loop: do m_ang=min_m,max_m
      y_loop: do j=1,ny
      x_loop: do i=1,nx
         r = Rvec(i)
         y = Yvec(j)
      !----------------------------
         A_pp(i,j,m_ang) = 0.0_rp

         A_pq(i,j,m_ang) = ((r**2)*((cl**4) - 2*(cl**2)*bhm*r + (bhs**2)*(r**2)))/(cl**4) 

         A_pf(i,j,m_ang) = 0.0_rp
      !----------------------------
         A_qp(i,j,m_ang) = cl**4/(&
            16*cl**2*bhm**2*(cl**2 + 2*bhm*r) &
         +  bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
         )

         A_qq(i,j,m_ang) = (&
            2*(cl**6 + cl**2*(bhs**2 - 8*bhm**2)*r**2 + 4*bhs**2*bhm*r**3) &
         )/(&
            16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2) &
         )

         A_qf(i,j,m_ang) = 0.0_rp
      !----------------------------
         A_fp(i,j,m_ang) = 0.0_rp

         A_fq(i,j,m_ang) = 0.0_rp

         A_ff(i,j,m_ang) = 0.0_rp
      !----------------------------
         B_pp(i,j,m_ang) = 0.0_rp 

         B_pq(i,j,m_ang) = cmplx(&
            -2*r*(-1 - spin + (r*(-2*bhs**2*r + cl**2*bhm*(3 + spin)))/cl**4) &
         ,&
            (-2*bhs*m_ang*r**2)/cl**2 &
         ,rp)

         B_pf(i,j,m_ang) = cmplx(&
            (-2*r*(-(bhs**2*r) + cl**2*bhm*(1 + spin)))/cl**4 &
         ,&
            (-2*bhs*m_ang*r)/cl**2 &
         ,rp)
      !----------------------------
         B_qp(i,j,m_ang) = cmplx( &
            -( &
            (32*cl**6*bhm**3 - 8*bhs**2*cl**4*bhm*(cl**2 + 4*bhm*r)) &
         /  (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2))**2) &
         , &
            0.0_rp &
         ,rp)

         B_qq(i,j,m_ang) = cmplx( &
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

         B_qf(i,j,m_ang) = cmplx(&
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
         B_fp(i,j,m_ang) = cmplx( &
            cl**4 &
         / &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2)) &
         , &
            0.0_rp &
         ,rp)

         B_fq(i,j,m_ang) = cmplx( &
            (2*(cl**6 + cl**2*(bhs**2 - 8*bhm**2)*r**2 + 4*bhs**2*bhm*r**3)) &
         / &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*r) + bhs**2*(-(cl**2 + 4*bhm*r)**2 + cl**4*y**2)) &
         , &
            0.0_rp &
         ,rp)

         B_ff(i,j,m_ang) = cmplx( &
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
      end do m_loop

   end subroutine teuk_init
!=============================================================================
   pure subroutine set_k(m_ang, & 
         p, q, f, &
         p_DR, q_DR, f_DR, f_coefs, f_laplacian, & 
         kp, kq, kf) 

      integer(ip), intent(in) :: m_ang 
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)    :: p, q, f 
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: &
         p_DR, q_DR, f_DR, f_laplacian, &
         kp, kq, kf 

      complex(rp), dimension(nx,0:lmax,min_m:max_m), intent(inout) :: f_coefs

      integer(ip) :: i, j

      call compute_DR(m_ang, p, p_DR)
      call compute_DR(m_ang, q, q_DR)
      call compute_DR(m_ang, f ,f_DR)

      call swal_laplacian(spin,m_ang,f,f_coefs,f_laplacian)

      y_loop: do j=1,ny
      x_loop: do i=1,nx
         kp(i,j,m_ang) = &
            A_pp(i,j,m_ang) * p_DR(i,j,m_ang) &
         +  A_pq(i,j,m_ang) * q_DR(i,j,m_ang) &
         +  A_pf(i,j,m_ang) * f_DR(i,j,m_ang) &
         +  B_pp(i,j,m_ang) * p(i,j,m_ang) &
         +  B_pq(i,j,m_ang) * q(i,j,m_ang) &
         +  B_pf(i,j,m_ang) * f(i,j,m_ang) &
         +  f_laplacian(i,j,m_ang) 

         kq(i,j,m_ang) = &
            A_qp(i,j,m_ang) * p_DR(i,j,m_ang) &
         +  A_qq(i,j,m_ang) * q_DR(i,j,m_ang) &
         +  A_qf(i,j,m_ang) * f_DR(i,j,m_ang) &
         +  B_qp(i,j,m_ang) * p(i,j,m_ang) &
         +  B_qq(i,j,m_ang) * q(i,j,m_ang) &
         +  B_qf(i,j,m_ang) * f(i,j,m_ang) 

         kf(i,j,m_ang) = &
            A_fp(i,j,m_ang) * p_DR(i,j,m_ang) &
         +  A_fq(i,j,m_ang) * q_DR(i,j,m_ang) &
         +  A_ff(i,j,m_ang) * f_DR(i,j,m_ang) &
         +  B_fp(i,j,m_ang) * p(i,j,m_ang) &
         +  B_fq(i,j,m_ang) * q(i,j,m_ang) &
         +  B_ff(i,j,m_ang) * f(i,j,m_ang) 
      end do x_loop
      end do y_loop

   end subroutine set_k
!=============================================================================
! RK4 time integrator
!=============================================================================
   pure subroutine teuk_time_step(m_ang, p, q, f) 
      integer(ip), intent(in) :: m_ang
      type(field), intent(inout) :: p, q, f 

      integer(ip) :: i, j
   !--------------------------------------------------------
      call set_k(m_ang, &
         p%n,  q%n,  f%n, & 
         p%DR, q%DR, f%DR, f%coefs, f%lap, &
         p%k1, q%k1, f%k1) 

      do j=1,ny
      do i=1,nx
         p%l2(i,j,m_ang)= p%n(i,j,m_ang)+0.5_rp*dt*p%k1(i,j,m_ang)
         q%l2(i,j,m_ang)= q%n(i,j,m_ang)+0.5_rp*dt*q%k1(i,j,m_ang)
         f%l2(i,j,m_ang)= f%n(i,j,m_ang)+0.5_rp*dt*f%k1(i,j,m_ang)
      end do
      end do
   !--------------------------------------------------------
      call set_k(m_ang, &
         p%l2, q%l2, f%l2, & 
         p%DR, q%DR, f%DR, f%coefs, f%lap, &
         p%k2, q%k2, f%k2) 

      do j=1,ny
      do i=1,nx
         p%l3(i,j,m_ang)= p%l2(i,j,m_ang)+0.5_rp*dt*p%k2(i,j,m_ang)
         q%l3(i,j,m_ang)= q%l2(i,j,m_ang)+0.5_rp*dt*q%k2(i,j,m_ang)
         f%l3(i,j,m_ang)= f%l2(i,j,m_ang)+0.5_rp*dt*f%k2(i,j,m_ang)
      end do
      end do
   !--------------------------------------------------------
      call set_k(m_ang, &
         p%l3, q%l3, f%l3, & 
         p%DR, q%DR, f%DR, f%coefs, f%lap, &
         p%k3, q%k3, f%k3) 

      do j=1,ny
      do i=1,nx
         p%l4(i,j,m_ang)= p%l3(i,j,m_ang)+dt*p%k3(i,j,m_ang)
         q%l4(i,j,m_ang)= q%l3(i,j,m_ang)+dt*q%k3(i,j,m_ang)
         f%l4(i,j,m_ang)= f%l3(i,j,m_ang)+dt*f%k3(i,j,m_ang)
      end do
      end do
   !--------------------------------------------------------
      call set_k(m_ang, &
         p%l4, q%l4, f%l4, & 
         p%DR, q%DR, f%DR, f%coefs, f%lap, &
         p%k4, q%k4, f%k4) 

      do j=1,ny
      do i=1,nx
         p%np1(i,j,m_ang)= p%n(i,j,m_ang) &
         +  (dt/6.0_rp)*(p%k1(i,j,m_ang)+2.0_rp*p%k2(i,j,m_ang)+2.0_rp*p%k3(i,j,m_ang)+p%k4(i,j,m_ang))

         q%np1(i,j,m_ang)= q%n(i,j,m_ang) &
         +  (dt/6.0_rp)*(q%k1(i,j,m_ang)+2.0_rp*q%k2(i,j,m_ang)+2.0_rp*q%k3(i,j,m_ang)+q%k4(i,j,m_ang))

         f%np1(i,j,m_ang)= f%n(i,j,m_ang) &
         +  (dt/6.0_rp)*(f%k1(i,j,m_ang)+2.0_rp*f%k2(i,j,m_ang)+2.0_rp*f%k3(i,j,m_ang)+f%k4(i,j,m_ang))
      end do
      end do
   end subroutine teuk_time_step
!=============================================================================
! independent residula: q - \partial_R f
!=============================================================================
   pure subroutine compute_q_indep_res(q, f, res) 
      type(field), intent(in)    :: q 
      type(field), intent(inout) :: f
      type(field), intent(out)   :: res

      integer(ip) :: m_ang

      do m_ang=min_m,max_m
         call compute_DR(m_ang,f%np1,f%DR)
      end do
      
      res%np1 = f%DR - q%np1

   end subroutine compute_q_indep_res
!=============================================================================
end module mod_teuk
