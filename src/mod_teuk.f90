!
! Teukolsky field solver
!
module mod_teuk
!=============================================================================
   use mod_prec
   use mod_field, only: field
   use mod_params, only: &
      dt, nx, ny, max_l, &
      min_m, max_m, &
      spin=>psi_spin, & 
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin

   use mod_cheb, only: R=>Rarr, compute_DR
   use mod_swal, only: Y=>Yarr, swal_laplacian

   use mod_scd_order_source, only: scd_order_source 
!=============================================================================
   implicit none
   private
   public :: teuk_init, teuk_time_step, compute_res_q 

   interface teuk_time_step 
      module procedure teuk_lin_time_step, teuk_scd_time_step
   end interface
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
      integer(ip) :: m_ang
      !----------------------------
      m_loop: do m_ang=min_m,max_m
      !----------------------------
         A_pp(:,:,m_ang) = 0.0_rp

         A_pq(:,:,m_ang) = ((R**2)*((cl**4) - 2*(cl**2)*bhm*R + (bhs**2)*(R**2)))/(cl**4) 

         A_pf(:,:,m_ang) = 0.0_rp
      !----------------------------
         A_qp(:,:,m_ang) = cl**4/(&
            16*cl**2*bhm**2*(cl**2 + 2*bhm*R) &
         +  bhs**2*(-(cl**2 + 4*bhm*R)**2 + cl**4*Y**2) &
         )

         A_qq(:,:,m_ang) = (&
            2*(cl**6 + cl**2*(bhs**2 - 8*bhm**2)*R**2 + 4*bhs**2*bhm*R**3) &
         )/(&
            16*cl**2*bhm**2*(cl**2 + 2*bhm*R) + bhs**2*(-(cl**2 + 4*bhm*R)**2 + cl**4*Y**2) &
         )

         A_qf(:,:,m_ang) = 0.0_rp
      !----------------------------
         A_fp(:,:,m_ang) = 0.0_rp

         A_fq(:,:,m_ang) = 0.0_rp

         A_ff(:,:,m_ang) = 0.0_rp
      !----------------------------
      !----------------------------
         B_pp(:,:,m_ang) = 0.0_rp 

         B_pq(:,:,m_ang) = cmplx(&
            -2*R*(-1 - spin + (R*(-2*bhs**2*R + cl**2*bhm*(3 + spin)))/cl**4) &
         ,&
            (-2*bhs*m_ang*R**2)/cl**2 &
         ,kind=rp)

         B_pf(:,:,m_ang) = cmplx(&
            (-2*R*(-(bhs**2*R) + cl**2*bhm*(1 + spin)))/cl**4 &
         ,&
            (-2*bhs*m_ang*R)/cl**2 &
         ,kind=rp)
      !----------------------------
         B_qp(:,:,m_ang) = cmplx( &
            -( &
            (32*cl**6*bhm**3 - 8*bhs**2*cl**4*bhm*(cl**2 + 4*bhm*R)) &
         /  (16*cl**2*bhm**2*(cl**2 + 2*bhm*R) + bhs**2*(-(cl**2 + 4*bhm*R)**2 + cl**4*Y**2))**2) &
         , &
            0.0_rp &
         ,kind=rp)

         B_qq(:,:,m_ang) = cmplx( &
            -( &
            ( &
               64*cl**4*bhm**3*(12*cl**2*bhm*R - cl**4*(-1 + spin) + 4*bhm**2*R**2*(4 + spin)) &
            +  2*bhs**4*R*((cl**2 + 4*bhm*R)**2*(3*cl**2 + 10*bhm*R) - 3*cl**4*(cl**2 + 6*bhm*R)*Y**2) &
            -  4*bhs**2*cl**2*bhm*( &
                  240*cl**2*bhm**2*R**2 + 32*bhm**3*R**3*(9 + spin) &
               -  2*cl**4*bhm*R*(-26 + 3*spin + (6 + spin)*Y**2) + cl**6*(4 + spin*(-1 + Y**2)) &
            ) &
            )/( &
            16*cl**2*bhm**2*(cl**2 + 2*bhm*R) + bhs**2*(-(cl**2 + 4*bhm*R)**2 + cl**4*Y**2))**2 &
            ) &
         , &
            ( &
               -2*bhs*cl**2*(4*m_ang*bhm*R + cl**2*(m_ang - spin*Y)) &
            )/( &
               16*cl**2*bhm**2*(cl**2 + 2*bhm*R) + bhs**2*(-(cl**2 + 4*bhm*R)**2 + cl**4*Y**2) &
            ) &
         ,kind=rp)

         B_qf(:,:,m_ang) = cmplx(&
            ( &
            -  2*cl**2*( &
                  128*cl**4*bhm**4*(1 + spin) &
               +  bhs**4*(32*bhm**2*R**2 - cl**4*(-1 + Y**2) -  12*cl**2*bhm*R*(-1 + Y**2)) &
               +  4*bhs**2*bhm**2*( &
                     16*bhm**2*R**2*(-1 + spin) - 16*cl**2*bhm*R*(3 + spin) + cl**4*(-6 - 5*spin + (2 + spin)*Y**2) &
                  ) &
            ) &
            )/( &
               16*cl**2*bhm**2*(cl**2 + 2*bhm*R) + bhs**2*(-(cl**2 + 4*bhm*R)**2 + cl**4*Y**2) &
            )**2 &
         ,&
            ( &
            -  8*bhs*cl**2*bhm*( &
                  8*cl**4*bhm**2*(m_ang + spin*Y) + bhs**2*(m_ang*(cl**2 + 4*bhm*R)**2 & 
               -  2*cl**2*(cl**2 + 4*bhm*R)*spin*Y + cl**4*m_ang*Y**2) &
               ) &
            )/( &
               16*cl**2*bhm**2*(cl**2 + 2*bhm*R) + bhs**2*(-(cl**2 + 4*bhm*R)**2 + cl**4*Y**2) &
            )**2 &
         ,kind=rp)
      !----------------------------
         B_fp(:,:,m_ang) = cmplx( &
            cl**4 &
         / &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*R) + bhs**2*(-(cl**2 + 4*bhm*R)**2 + cl**4*Y**2)) &
         , &
            0.0_rp &
         ,kind=rp)

         B_fq(:,:,m_ang) = cmplx( &
            (2*(cl**6 + cl**2*(bhs**2 - 8*bhm**2)*R**2 + 4*bhs**2*bhm*R**3)) &
         / &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*R) + bhs**2*(-(cl**2 + 4*bhm*R)**2 + cl**4*Y**2)) &
         , &
            0.0_rp &
         ,kind=rp)

         B_ff(:,:,m_ang) = cmplx( &
            -( & 
               (2*bhs**2*R*(cl**2 + 6*bhm*R) + 4*cl**2*bhm*(cl**2*spin - 2*bhm*R*(2 + spin))) &
         / &
               (-16*cl**2*bhm**2*(cl**2 + 2*bhm*R) + bhs**2*((cl**2 + 4*bhm*R)**2 - cl**4*Y**2)) &
            ) &
         , &
            (-2*bhs*cl**2*(4*m_ang*bhm*R + cl**2*(m_ang - spin*Y))) &
         /  &
            (16*cl**2*bhm**2*(cl**2 + 2*bhm*R) + bhs**2*(-(cl**2 + 4*bhm*R)**2 + cl**4*Y**2)) &
         ,kind=rp)
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
      complex(rp), dimension(nx,0:max_l,min_m:max_m), intent(inout) :: f_coefs

      call compute_DR(m_ang, p, p_DR)
      call compute_DR(m_ang, q, q_DR)
      call compute_DR(m_ang, f ,f_DR)

      call swal_laplacian(spin,m_ang,f,f_coefs,f_laplacian)

      !-------------------------------------
      kp(:,:,m_ang) = &
         A_pp(:,:,m_ang) * p_DR(:,:,m_ang) &
      +  A_pq(:,:,m_ang) * q_DR(:,:,m_ang) &
      +  A_pf(:,:,m_ang) * f_DR(:,:,m_ang) &
      +  B_pp(:,:,m_ang) * p(:,:,m_ang) &
      +  B_pq(:,:,m_ang) * q(:,:,m_ang) &
      +  B_pf(:,:,m_ang) * f(:,:,m_ang) &
      +  f_laplacian(:,:,m_ang) 
      !-------------------------------------
      kq(:,:,m_ang) = &
         A_qp(:,:,m_ang) * p_DR(:,:,m_ang) &
      +  A_qq(:,:,m_ang) * q_DR(:,:,m_ang) &
      +  A_qf(:,:,m_ang) * f_DR(:,:,m_ang) &
      +  B_qp(:,:,m_ang) * p(:,:,m_ang) &
      +  B_qq(:,:,m_ang) * q(:,:,m_ang) &
      +  B_qf(:,:,m_ang) * f(:,:,m_ang) 
      !-------------------------------------
      kf(:,:,m_ang) = &
         A_fp(:,:,m_ang) * p_DR(:,:,m_ang) &
      +  A_fq(:,:,m_ang) * q_DR(:,:,m_ang) &
      +  A_ff(:,:,m_ang) * f_DR(:,:,m_ang) &
      +  B_fp(:,:,m_ang) * p(:,:,m_ang) &
      +  B_fq(:,:,m_ang) * q(:,:,m_ang) &
      +  B_ff(:,:,m_ang) * f(:,:,m_ang) 

   end subroutine set_k
!=============================================================================
! RK4 time integrator: linear teukolsky wave
!=============================================================================
   pure subroutine teuk_lin_time_step(m_ang, p, q, f) 
      integer(ip), intent(in)    :: m_ang
      type(field), intent(inout) :: p, q, f 
   !--------------------------------------------------------
   ! if first time then k1 has not been set from k5
   !--------------------------------------------------------
      if (f%first_time) then
         call set_k(m_ang, &
            p%n,  q%n,  f%n, & 
            p%DR, q%DR, f%DR, f%coefs_swal, f%lap, &
            p%k1, q%k1, f%k1) 

         p%first_time = .false.
         q%first_time = .false.
         f%first_time = .false.
      end if

      p%l2(:,:,m_ang)= p%n(:,:,m_ang)+0.5_rp*dt*p%k1(:,:,m_ang)
      q%l2(:,:,m_ang)= q%n(:,:,m_ang)+0.5_rp*dt*q%k1(:,:,m_ang)
      f%l2(:,:,m_ang)= f%n(:,:,m_ang)+0.5_rp*dt*f%k1(:,:,m_ang)
   !--------------------------------------------------------
      call set_k(m_ang, &
         p%l2, q%l2, f%l2, & 
         p%DR, q%DR, f%DR, f%coefs_swal, f%lap, &
         p%k2, q%k2, f%k2) 

      p%l3(:,:,m_ang)= p%n(:,:,m_ang)+0.5_rp*dt*p%k2(:,:,m_ang)
      q%l3(:,:,m_ang)= q%n(:,:,m_ang)+0.5_rp*dt*q%k2(:,:,m_ang)
      f%l3(:,:,m_ang)= f%n(:,:,m_ang)+0.5_rp*dt*f%k2(:,:,m_ang)
   !--------------------------------------------------------
      call set_k(m_ang, &
         p%l3, q%l3, f%l3, & 
         p%DR, q%DR, f%DR, f%coefs_swal, f%lap, &
         p%k3, q%k3, f%k3) 

      p%l4(:,:,m_ang)= p%n(:,:,m_ang)+dt*p%k3(:,:,m_ang)
      q%l4(:,:,m_ang)= q%n(:,:,m_ang)+dt*q%k3(:,:,m_ang)
      f%l4(:,:,m_ang)= f%n(:,:,m_ang)+dt*f%k3(:,:,m_ang)
   !--------------------------------------------------------
      call set_k(m_ang, &
         p%l4, q%l4, f%l4, & 
         p%DR, q%DR, f%DR, f%coefs_swal, f%lap, &
         p%k4, q%k4, f%k4) 

      p%np1(:,:,m_ang)= p%n(:,:,m_ang) &
      +  (dt/6.0_rp)*(p%k1(:,:,m_ang)+2.0_rp*p%k2(:,:,m_ang)+2.0_rp*p%k3(:,:,m_ang)+p%k4(:,:,m_ang))

      q%np1(:,:,m_ang)= q%n(:,:,m_ang) &
      +  (dt/6.0_rp)*(q%k1(:,:,m_ang)+2.0_rp*q%k2(:,:,m_ang)+2.0_rp*q%k3(:,:,m_ang)+q%k4(:,:,m_ang))

      f%np1(:,:,m_ang)= f%n(:,:,m_ang) &
      +  (dt/6.0_rp)*(f%k1(:,:,m_ang)+2.0_rp*f%k2(:,:,m_ang)+2.0_rp*f%k3(:,:,m_ang)+f%k4(:,:,m_ang))
   !------------------------------------------------------------
   ! want k5 for computing source term and independent residuals
   !------------------------------------------------------------
      call set_k(m_ang, &
         p%np1, q%np1, f%np1, & 
         p%DR,  q%DR,  f%DR, f%coefs_swal, f%lap, &
         p%k5,  q%k5,  f%k5) 

   end subroutine teuk_lin_time_step
!=============================================================================
! RK4 time integrator: linear teukolsky wave with second order source term
!=============================================================================
   pure subroutine teuk_scd_time_step(m_ang, src, p, q, f) 
      integer(ip),            intent(in)    :: m_ang
      type(scd_order_source), intent(in)    :: src
      type(field),            intent(inout) :: p, q, f 
   !--------------------------------------------------------
   ! if first time then k1 has not been set from k5
   !--------------------------------------------------------
      if (f%first_time) then
         call set_k(m_ang, &
            p%n,  q%n,  f%n, & 
            p%DR, q%DR, f%DR, f%coefs_swal, f%lap, &
            p%k1, q%k1, f%k1) 

         p%k1(:,:,m_ang) = p%k1(:,:,m_ang) + src%n(:,:,m_ang)

         p%first_time = .false.
         q%first_time = .false.
         f%first_time = .false.
      end if

      p%l2(:,:,m_ang)= p%n(:,:,m_ang)+0.5_rp*dt*p%k1(:,:,m_ang)
      q%l2(:,:,m_ang)= q%n(:,:,m_ang)+0.5_rp*dt*q%k1(:,:,m_ang)
      f%l2(:,:,m_ang)= f%n(:,:,m_ang)+0.5_rp*dt*f%k1(:,:,m_ang)
   !--------------------------------------------------------
      call set_k(m_ang, &
         p%l2, q%l2, f%l2, & 
         p%DR, q%DR, f%DR, f%coefs_swal, f%lap, &
         p%k2, q%k2, f%k2) 

      p%k2(:,:,m_ang) = p%k2(:,:,m_ang) + src%n1h(:,:,m_ang)

      p%l3(:,:,m_ang)= p%n(:,:,m_ang)+0.5_rp*dt*p%k2(:,:,m_ang)
      q%l3(:,:,m_ang)= q%n(:,:,m_ang)+0.5_rp*dt*q%k2(:,:,m_ang)
      f%l3(:,:,m_ang)= f%n(:,:,m_ang)+0.5_rp*dt*f%k2(:,:,m_ang)
   !--------------------------------------------------------
      call set_k(m_ang, &
         p%l3, q%l3, f%l3, & 
         p%DR, q%DR, f%DR, f%coefs_swal, f%lap, &
         p%k3, q%k3, f%k3) 

      p%k3(:,:,m_ang) = p%k3(:,:,m_ang) + src%n1h(:,:,m_ang)

      p%l4(:,:,m_ang)= p%n(:,:,m_ang)+dt*p%k3(:,:,m_ang)
      q%l4(:,:,m_ang)= q%n(:,:,m_ang)+dt*q%k3(:,:,m_ang)
      f%l4(:,:,m_ang)= f%n(:,:,m_ang)+dt*f%k3(:,:,m_ang)
   !--------------------------------------------------------
      call set_k(m_ang, &
         p%l4, q%l4, f%l4, & 
         p%DR, q%DR, f%DR, f%coefs_swal, f%lap, &
         p%k4, q%k4, f%k4) 

      p%k4(:,:,m_ang) = p%k4(:,:,m_ang) + src%np1(:,:,m_ang)

      p%np1(:,:,m_ang)= p%n(:,:,m_ang) &
      +  (dt/6.0_rp)*(p%k1(:,:,m_ang)+2.0_rp*p%k2(:,:,m_ang)+2.0_rp*p%k3(:,:,m_ang)+p%k4(:,:,m_ang))

      q%np1(:,:,m_ang)= q%n(:,:,m_ang) &
      +  (dt/6.0_rp)*(q%k1(:,:,m_ang)+2.0_rp*q%k2(:,:,m_ang)+2.0_rp*q%k3(:,:,m_ang)+q%k4(:,:,m_ang))

      f%np1(:,:,m_ang)= f%n(:,:,m_ang) &
      +  (dt/6.0_rp)*(f%k1(:,:,m_ang)+2.0_rp*f%k2(:,:,m_ang)+2.0_rp*f%k3(:,:,m_ang)+f%k4(:,:,m_ang))
   !------------------------------------------------------------
   ! want k5 for computing source term and independent residuals
   !------------------------------------------------------------
      call set_k(m_ang, &
         p%np1, q%np1, f%np1, & 
         p%DR,  q%DR,  f%DR, f%coefs_swal, f%lap, &
         p%k5,  q%k5,  f%k5) 

      p%k5(:,:,m_ang) = p%k1(:,:,m_ang) + src%np1(:,:,m_ang)
   end subroutine teuk_scd_time_step
!=============================================================================
! independent residula: q - \partial_R f
!=============================================================================
   pure subroutine compute_res_q(q, f, res) 
      type(field), intent(in)    :: q 
      type(field), intent(inout) :: f
      type(field), intent(out)   :: res

      integer(ip) :: m_ang

      do m_ang=min_m,max_m
         call compute_DR(m_ang,f%np1,f%DR)
      end do
      
      res%np1 = f%DR - q%np1

   end subroutine compute_res_q
!=============================================================================
end module mod_teuk
