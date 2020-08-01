!
!  Compute source term for \Psi_4^{(2)} evolution
!
!=============================================================================
module mod_scnd_order_source
!=============================================================================
   use mod_prec

   use mod_cheb,     only: R=>Rarr, compute_DR
   use mod_swal,     only: swal_lower, cy=>cyarr
   use mod_field,    only: field, set_level
   use mod_ghp,      only: set_edth, set_edth_prime, set_thorn, set_thorn_prime
   use mod_bkgrd_np, only: mu_0, ta_0, pi_0, rh_0, psi2_0

   use mod_fields_list, only: psi4_lin_f, psi3, psi2, la, pi, muhll, hlmb, hmbmb 

   use mod_params,   only: &
      dt, nx, ny, min_m, max_m, max_l, &
      cl=>compactification_length, &
      bhs=>black_hole_spin

!=============================================================================
   implicit none
   private

   public :: scnd_order_source, &
      scnd_order_source_init, scnd_order_source_compute
!=============================================================================
! scnd_order_source type:
! have 5 levels so can take time derivatives
!-----------------------------------------------------------------------------
   type :: scnd_order_source

   character(:), allocatable :: name
   character(:), allocatable :: error

   integer(ip) :: &
      pre_edth_prime_spin,  pre_edth_prime_boost,  pre_edth_prime_falloff, &
      pre_thorn_prime_spin, pre_thorn_prime_boost, pre_thorn_prime_falloff

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
   pure subroutine scnd_order_source_init(name, sf)
      character(*),            intent(in)  :: name ! field name
      type(scnd_order_source), intent(out) :: sf

      sf % name = name

      ! make empty string long enough to hold error message
      sf % error = "                                                              "

      sf % pre_edth_prime_spin  = -2_ip
      sf % pre_thorn_prime_spin = -1_ip

      sf % pre_edth_prime_spin   = -1_ip
      sf % pre_thorn_prime_boost = -2_ip

      sf % pre_edth_prime_falloff  = 0_ip
      sf % pre_thorn_prime_falloff = 0_ip

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

   end subroutine scnd_order_source_init
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
      integer(ip) :: mt_ang

      mt_ang = m1_ang + m2_ang
      !-----------------------------------------------------------------------
      call set_level(step,m1_ang,psi4_lin_f)
      call set_level(step,m1_ang,psi3)
      call set_level(step,m1_ang,la)
      call set_level(step,m1_ang,pi)
      call set_level(step,m1_ang,hmbmb)
      call set_level(step,m1_ang,hlmb)
      call set_level(step,m1_ang,muhll)

      call set_level(step,-m1_ang,hlmb)
      call set_level(step,-m1_ang,hmbmb)
      !-----------------------------------------------------------------------
      call set_level(step,m2_ang,psi4_lin_f)
      call set_level(step,m2_ang,psi3)
      call set_level(step,m2_ang,psi2)
      call set_level(step,m2_ang,hmbmb)
      call set_level(step,m2_ang,hlmb)
      call set_level(step,m2_ang,muhll)

      call set_thorn_prime(step,m2_ang,psi4_lin_f)
      call set_thorn_prime(step,m2_ang,psi3)
      call set_thorn_prime(step,m2_ang,muhll)
      call set_thorn_prime(step,m2_ang,hlmb)

      call set_edth(step,m2_ang,psi3)
      call set_edth(step,m2_ang,hmbmb)
      call set_edth(step,m2_ang,hlmb)

      call set_edth_prime(step,m2_ang,psi4_lin_f)

      call set_level(step,-m2_ang,pi)
      call set_level(step,-m2_ang,hlmb)
      call set_level(step,-m2_ang,hmbmb)

      call set_thorn_prime(step,-m2_ang,hlmb)

      call set_edth(step,-m2_ang,hlmb)
      call set_edth(step,-m2_ang,hmbmb)
      !--------------------------------------------------------------------
      sf % pre_thorn_prime_np1(:,:,mt_ang) = &
      sf % pre_thorn_prime_np1(:,:,mt_ang) &
      !-------------------------------------------------
      + 0.5_rp*(muhll%level(:,:,m1_ang)/mu_0)*( &
            psi4_lin_f%thorn_prime(:,:,m1_ang) &

         +  mu_0*psi4_lin_f%level(:,:,m1_ang) &
         ) &
      !-------------------------------------------------
      + psi4_lin_f%level(:,:,m1_ang)*( &
            0.5_rp*hlmb%edth(:,:,m2_ang) &
         +  ( &
               conjg(pi_0) &
            +  2.0_rp*ta_0 &
            )*conjg(hlmb%level(:,:,-m2_ang)) &

         +  (muhll%thorn_prime(:,:,m2_ang)/mu_0) &
         +  conjg(mu_0)*(muhll%level(:,:,m2_ang)/mu_0) &

         -  0.5_rp*conjg(hlmb%edth(:,:,-m2_ang)) &
         +  0.5_rp*( &
               5.0_rp*pi_0 &
            +  4.0_rp*conjg(ta_0) &
            )*conjg(hlmb%level(:,:,-m2_ang)) &

         ) &
      !-------------------------------------------------
      -  0.5_rp*psi3%level(:,:,m1_ang)*( &
            hmbmb%edth(:,:,m2_ang) &
         +  ( &
               conjg(pi_0) &
            +  ta_0 &
            )*hmbmb%level(:,:,m2_ang) &

         +  hlmb%thorn_prime(:,:,m2_ang) &
         + ( &
            -  2.0_rp*mu_0 &
            +  conjg(mu_0) &
            )*hlmb%level(:,:,m2_ang)  &
         ) &
      !-------------------------------------------------
         -  hlmb%level(:,:,m1_ang)*psi3%thorn_prime(:,:,m2_ang) &

         +  0.5_rp*hmbmb%level(:,:,m1_ang)*psi3%edth(:,:,m2_ang) &

         +  4.0_rp*pi%level(:,:,m1_ang)*psi3%level(:,:,m2_ang) &

         -  3.0_rp*la%level(:,:,m1_ang)*psi2%level(:,:,m2_ang)
      !--------------------------------------------------------------------
      sf % pre_edth_prime_np1(:,:,mt_ang) = &
      sf % pre_edth_prime_np1(:,:,mt_ang) &
      !-------------------------------------------------
      -  conjg(hlmb%level(:,:,-m1_ang))*psi4_lin_f%thorn_prime(:,:,m2_ang) &

      -  (mu_0+2.0_rp*conjg(mu_0))*psi4_lin_f%level(:,:,m2_ang) &

      +  0.5_rp*conjg(hmbmb%level(:,:,-m1_ang))*psi4_lin_f%edth_prime(:,:,m2_ang) &

      !-------------------------------------------------
      +  psi4_lin_f%level(:,:,m1_ang)*( &
            conjg(pi%level(:,:,-m2_ang)) &

         -  conjg(hlmb%thorn_prime(:,:,-m2_ang)) &

         +  conjg(hmbmb%edth(:,:,-m2_ang)) &

         -  0.5_rp*(pi_0+conjg(ta_0))*conjg(hmbmb%level(:,:,-m2_ang)) &
      )
      !--------------------------------------------------------------------
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

      call compute_DT("pre_edth_prime",m_ang,sf)
      call compute_DT("pre_edth_prime",m_ang,sf)
      call swal_lower( &
         sf%pre_edth_prime_spin, &
         m_ang, &
         sf%pre_edth_prime_np1, &
         sf%coefs_swal, &
         sf%lowered)

      call set_edth_prime( &
         m_ang, &
         sf%pre_edth_prime_spin, &
         sf%pre_edth_prime_boost, &
         sf%pre_edth_prime_np1, &
         sf%DT, &
         sf%lowered, &
         sf%edth_prime)

      call compute_DT("pre_thorn_prime",m_ang,sf)
      call compute_DR(m_ang,sf%pre_thorn_prime_np1,sf%DR)

      call set_thorn_prime( &
         m_ang, &
         sf%pre_thorn_prime_falloff, &
         sf%pre_thorn_prime_np1, &
         sf%DT, &
         sf%DR, &
         sf%thorn_prime)

      !-----------------------------------------------------------------------
      sf%np1(:,:,m_ang) = &
         sf%thorn_prime(:,:,m_ang) &
      +  (4.0_rp*mu_0+conjg(mu_0))*sf%pre_thorn_prime_np1(:,:,m_ang) &

      +  sf%edth_prime(:,:,m_ang) &
      +  (4.0_rp*pi_0-conjg(ta_0))*sf%pre_edth_prime_np1(:,:,m_ang) 
      !-----------------------------------------------------
      ! multiply by prefactor
      sf%np1(:,:,m_ang) = &
         (cl**4 + (bhs*R*cy)**2) * sf%np1(:,:,m_ang) 
      !-----------------------------------------------------
      ! estimate halfstep value 
      sf%n1h(:,:,m_ang) = 0.5_rp*(sf%np1(:,:,m_ang) + sf%n(:,:,m_ang))
      !-----------------------------------------------------------------------
   end subroutine scnd_order_source_compute
!=============================================================================
end module mod_scnd_order_source
