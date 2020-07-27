!
!  Compute source term for \Psi_4^{(2)} evolution
!
!=============================================================================
module mod_scnd_order_source
!=============================================================================
   use mod_prec

   use mod_cheb,     only: R, compute_DR
   use mod_field,    only: field, set_field, set_level
   use mod_ghp,      only: set_edth, set_edth_prime, set_thorn, set_thorn_prime
   use mod_bkgrd_np, only: mu_0, ta_0, pi_0, rh_0, psi2_0
   use mod_teuk,     only: psi4_lin_f
   use mod_metric_recon, &
                     only: psi3, psi2, la, pi, muhll, hlmb, hmbmb

   use mod_params,   only: &
      dt, nx, ny, min_m, max_m, max_l, &
      cl=>compactification_length, &
      bhm=>black_hole_mass

!=============================================================================
   implicit none
   private

   public :: scnd_order_source, set_scnd_order_source, &
      scnd_order_source_zero, &
      scnd_order_source_m1_plus_m2,  &
      scnd_order_source_compute
!=============================================================================
! scnd_order_source type:
! have 5 levels so can take time derivatives
!-----------------------------------------------------------------------------
   type :: scnd_order_source

   character(:), allocatable :: name
   character(:), allocatable :: error

   integer(ip) :: &
      pre_edth_spin,  pre_edth_boost,  pre_edth_falloff, &
      pre_thorn_spin, pre_thorn_boost, pre_thorn_falloff

   complex(rp) :: &
      n1h(nx,ny,min_m:max_m), &
      n(  nx,ny,min_m:max_m), &
      np1(nx,ny,min_m:max_m), &

      level(  nx,ny,min_m:max_m), &
      DT(     nx,ny,min_m:max_m), &
      DR(     nx,ny,min_m:max_m), &
      raised( nx,ny,min_m:max_m), &
      lowered(nx,ny,min_m:max_m), & 

      coefs_swal(nx,0:max_l,min_m:max_m), &
      coefs_cheb(nx,     ny,min_m:max_m), &
      coefs_both(nx,0:max_l,min_m:max_m), &

       edth_prime(nx,ny,min_m:max_m), &
      thorn_prime(nx,ny,min_m:max_m), &

      pre_edth_prime_np1( nx,ny,min_m:max_m), &
      pre_edth_prime_n(   nx,ny,min_m:max_m), &
      pre_edth_prime_nm1( nx,ny,min_m:max_m), &
      pre_edth_prime_nm2( nx,ny,min_m:max_m), &
      pre_edth_prime_nm3( nx,ny,min_m:max_m), &

      pre_thorn_prime_np1(nx,ny,min_m:max_m), &
      pre_thorn_prime_n(  nx,ny,min_m:max_m), &
      pre_thorn_prime_nm1(nx,ny,min_m:max_m), &
      pre_thorn_prime_nm2(nx,ny,min_m:max_m), &
      pre_thorn_prime_nm3(nx,ny,min_m:max_m)

   end type scnd_order_source
!=============================================================================
   contains
!=============================================================================
   pure subroutine set_scnd_order_source(name, sf)
      character(*),            intent(in)  :: name ! field name
      type(scnd_order_source), intent(out) :: sf

      sf % name = name

      ! make empty string long enough to hold error message
      sf % error = "                                                              "

      sf % pre_edth_spin  = 0_ip
      sf % pre_thorn_spin = 0_ip

      sf % pre_edth_boost  = 0_ip
      sf % pre_thorn_boost = 0_ip

      sf % pre_edth_falloff  = 0_ip
      sf % pre_thorn_falloff = 0_ip

      sf % n1h = 0.0_rp
      sf % n   = 0.0_rp
      sf % np1 = 0.0_rp

      sf % level = 0.0_rp
      sf % DT    = 0.0_rp
      sf % DR    = 0.0_rp
      sf % raised  = 0.0_rp
      sf % lowered = 0.0_rp 

      sf % coefs_swal = 0.0_rp
      sf % coefs_cheb = 0.0_rp
      sf % coefs_both = 0.0_rp

      sf %  edth_prime = 0.0_rp
      sf % thorn_prime = 0.0_rp

      sf % pre_edth_prime_np1 = 0.0_rp
      sf % pre_edth_prime_n   = 0.0_rp
      sf % pre_edth_prime_nm1 = 0.0_rp
      sf % pre_edth_prime_nm2 = 0.0_rp
      sf % pre_edth_prime_nm3 = 0.0_rp

      sf % pre_thorn_prime_np1 = 0.0_rp
      sf % pre_thorn_prime_n   = 0.0_rp
      sf % pre_thorn_prime_nm1 = 0.0_rp
      sf % pre_thorn_prime_nm2 = 0.0_rp
      sf % pre_thorn_prime_nm3 = 0.0_rp

   end subroutine set_scnd_order_source
!=============================================================================
   pure subroutine compute_DT(pre_type, m_ang, sf)
      character(*),            intent(in)    :: pre_type
      integer(ip),             intent(in)    :: m_ang
      type(scnd_order_source), intent(inout) :: sf

      select case (pre_type)
         case ("pre_edth_prime")
            sf%DT(:,:,m_ang) = (1.0_rp/dt) * (&
               (25.0_rp/12.0_rp)*(sf%pre_edth_prime_np1(:,:,m_ang)) &
            -  (4.0_rp         )*(sf%pre_edth_prime_n(  :,:,m_ang)) &
            +  (3.0_rp         )*(sf%pre_edth_prime_nm1(:,:,m_ang)) &
            -  (4.0_rp/3.0_rp  )*(sf%pre_edth_prime_nm2(:,:,m_ang)) &
            +  (1.0_rp/4.0_rp  )*(sf%pre_edth_prime_nm3(:,:,m_ang)) &
            )
         case ("pre_thorn_prime")
            sf%DT(:,:,m_ang) = (1.0_rp/dt) * (&
               (25.0_rp/12.0_rp)*(sf%pre_thorn_prime_np1(:,:,m_ang)) &
            -  (4.0_rp         )*(sf%pre_thorn_prime_n(  :,:,m_ang)) &
            +  (3.0_rp         )*(sf%pre_thorn_prime_nm1(:,:,m_ang)) &
            -  (4.0_rp/3.0_rp  )*(sf%pre_thorn_prime_nm2(:,:,m_ang)) &
            +  (1.0_rp/4.0_rp  )*(sf%pre_thorn_prime_nm3(:,:,m_ang)) &
            )
         case default
            write (sf%error,*) "ERROR(compute_DT): inter_type= ", pre_type 
      end select
   end subroutine compute_DT
!=============================================================================
   pure subroutine scnd_order_source_zero(m_ang,sf)
      integer(ip),             intent(in)    :: m_ang
      type(scnd_order_source), intent(inout) :: sf

      sf % pre_edth_prime_np1(:,:,m_ang)  = 0.0_rp
      sf % pre_thorn_prime_np1(:,:,m_ang) = 0.0_rp
      
   end subroutine scnd_order_source_zero
!=============================================================================
! add together m1 and m2
! look at section II D
!=============================================================================
   subroutine scnd_order_source_m1_plus_m2(m1_ang, m2_ang, sf)
      integer(ip), intent(in)                :: m1_ang, m2_ang
      type(scnd_order_source), intent(inout) :: sf

      integer(ip), parameter :: step = 5_ip
      integer(ip) :: i,j, mt_ang

      mt_ang = m1_ang + m2_ang
      !-----------------------------------------------------------------------
      call set_level(step,m1_ang,psi4_lin_f)
      call set_level(step,m1_ang,psi3)

      call set_level(step,m1_ang,la)
      call set_level(step,m1_ang,pi)

      call set_level(step,m1_ang,hlmb)
      call set_level(step,m1_ang,hmbmb)
      call set_level(step,m1_ang,muhll)

      call set_level(step,-m1_ang,hlmb)
      call set_level(step,-m1_ang,hmbmb)
      !-----------------------------------------------------------------------
      call set_level(step,m2_ang,psi4_lin_f)
      call set_level(step,m2_ang,psi3)
      call set_level(step,m2_ang,psi2)
      call set_level(step,m2_ang,hlmb)

      call set_thorn_prime(step,m2_ang,psi4_lin_f)
      call set_thorn_prime(step,m2_ang,psi3)
      call set_thorn_prime(step,m2_ang,muhll)
      call set_thorn_prime(step,m2_ang,hlmb)

      call set_edth(step,m2_ang,psi3)
      call set_edth(step,m2_ang,hmbmb)
      call set_edth(step,m2_ang,hlmb)

      call set_level(step,-m2_ang,pi)
      call set_level(step,-m2_ang,hlmb)
      call set_level(step,-m2_ang,hmbmb)

      call set_thorn_prime(step,-m2_ang,hlmb)

      call set_edth(step,-m2_ang,psi4_lin_f)
      call set_edth(step,-m2_ang,hlmb)
      call set_edth(step,-m2_ang,hmbmb)
      !-----------------------------------------------------------------------
      y_loop: do j=1,ny
      x_loop: do i=1,nx
         !--------------------------------------------------------------------
         sf % pre_thorn_prime_np1(i,j,mt_ang) = &
         sf % pre_thorn_prime_np1(i,j,mt_ang) &

         + 0.5_rp*(muhll%level(i,j,m1_ang)/mu_0(i,j))*( &
               psi4_lin_f%thorn_prime(i,j,m1_ang) &

            +  mu_0(i,j)*psi4_lin_f%level(i,j,m1_ang) &
            ) &

         + psi4_lin_f%level(i,j,m1_ang)*( &
               0.5_rp*hlmb%edth(i,j,m2_ang) &

            +  ( &
                  conjg(pi_0(i,j)) &
               +  2.0_rp*ta_0(i,j) &
               )*hlmb%level(i,j,m2_ang) &

            +  (muhll%thorn_prime(i,j,m2_ang)/mu_0(i,j)) &

            +  conjg(mu_0(i,j))*(muhll%level(i,j,m2_ang)/mu_0(i,j)) &

            -  0.5_rp*conjg(hlmb%edth(i,j,-m2_ang)) &

            +  0.5_rp*( &
                  5.0_rp*pi_0(i,j) &
               +  4.0_rp*conjg(ta_0(i,j)) &
               )*conjg(hlmb%level(i,j,-m2_ang)) &
            ) &

         -  0.5_rp*psi3%level(i,j,m1_ang)*( &
               hmbmb%edth(i,j,m2_ang) &

            +  ( &
                  conjg(pi_0(i,j)) &
               +  ta_0(i,j) &
               )*hmbmb%level(i,j,m2_ang) &

            +  hlmb%thorn_prime(i,j,m2_ang) &

            + ( &
               -  2.0_rp*mu_0(i,j) &
               +  conjg(mu_0(i,j)) &
               )*hlmb%level(i,j,m2_ang)  &
            ) &

            -  hlmb%level(i,j,m1_ang)*psi3%thorn_prime(i,j,m2_ang) &

            +  0.5_rp*hmbmb%level(i,j,m1_ang)*psi3%level(i,j,m2_ang) &

            +  4.0_rp*pi%level(i,j,m1_ang)*psi3%level(i,j,m2_ang) &

            -  3.0_rp*la%level(i,j,m1_ang)*psi2%level(i,j,m2_ang)
         !--------------------------------------------------------------------
         sf % pre_edth_prime_np1(i,j,mt_ang) = &
         sf % pre_edth_prime_np1(i,j,mt_ang) &

         -  conjg(hlmb%level(i,j,-m1_ang))*psi4_lin_f%thorn_prime(i,j,m2_ang) &

         -  (mu_0(i,j)+2.0_rp*conjg(mu_0(i,j)))*psi4_lin_f%level(i,j,m2_ang) &

         +  0.5_rp*conjg(hmbmb%level(i,j,-m1_ang))*conjg(psi4_lin_f%edth(i,j,-m2_ang)) &

         +  psi4_lin_f%level(i,j,m1_ang)*( &
               conjg(pi%level(i,j,-m2_ang)) &

            -  conjg(hlmb%thorn_prime(i,j,-m2_ang)) &

            +  conjg(hmbmb%thorn_prime(i,j,-m2_ang)) &

            -  0.5_rp*(pi_0(i,j)+conjg(ta_0(i,j)))*conjg(hmbmb%level(i,j,-m2_ang)) &
         )
         !--------------------------------------------------------------------
      end do x_loop
      end do y_loop
   end subroutine scnd_order_source_m1_plus_m2
!=============================================================================
   subroutine scnd_order_source_compute(m_ang, sf)
      integer(ip),             intent(in)    :: m_ang
      type(scnd_order_source), intent(inout) :: sf
 
      integer(ip) :: m1_ang, m2_ang 
 
      sf % n(:,:,m_ang) = sf % np1(:,:,m_ang)
       
      call scnd_order_source_zero(m_ang,sf)

      do m1_ang=min_m,max_m
      do m2_ang=min_m,max_m

         if (m1_ang+m2_ang==m_ang) then
            call scnd_order_source_m1_plus_m2( m1_ang, m2_ang, sf)

         else if (m1_ang-m2_ang==m_ang) then
            call scnd_order_source_m1_plus_m2( m1_ang,-m2_ang, sf)

         else if (-m1_ang+m2_ang==m_ang) then
            call scnd_order_source_m1_plus_m2(-m1_ang, m2_ang, sf)

         else
            continue
         end if
      end do
      end do



   end subroutine scnd_order_source_compute
!=============================================================================
end module mod_scnd_order_source
