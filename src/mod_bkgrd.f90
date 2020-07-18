!
! Newman-Penrose fields for the background manifold
!
module mod_bkgrd_np
!-----------------------------------------------------------------------------
   use mod_prec
   use mod_params, only: &
      nx, ny, & 
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin

   use mod_cheb, only: Rvec=>R
   use mod_swal, only: cyvec=>cy, syvec=>sy

   implicit none
!-----------------------------------------------------------------------------
! all variables are public
!-----------------------------------------------------------------------------
   complex(rp), dimension(nx,ny), protected :: &
      mu_0, ta_0, pi_0, rh_0, ep_0, &
      thorn_prime_ta_0, &
      psi2_0
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
   subroutine bkgrd_np_init()

      integer(ip) :: i, j
      real(rp) :: r, cy, sy

      y_loop: do j=1,ny
      x_loop: do i=1,nx
         r = Rvec(i)
         cy = cyvec(j)
         sy = syvec(j)
      !----------------------------
         mu_0(i,j) = cmplx( &
            - (cl**2) &
            / ((cl**4) + (bhs**2)*(cy**2)*(r**2)) &
            , &
            - (bhs*cy*r) & 
            / ((cl**4) + (bhs**2)*(cy**2)*(r**2)) &
            , & 
            kind=rp)
      !----------------------------
         ta_0(i,j) = cmplx( &
            - (sqrt(2.0_rp)*(bhs**2)*cy*(cl**2)*r*sy) &
            / (((cl**4) + (bhs**2)*(cy**2)*(r**2))**2) &
            , &
            (bhs*sy*(4.0_rp*(cl**4) + (bhs**2)*(r**2)*(-1.0_rp - 3.0_rp*(cy**2) + (sy**2)))) &
            / (sqrt(2.0_rp)*((2*(cl**4) + (bhs**2)*(r**2)*(1.0_rp + (cy**2) - (sy**2)))**2)) &
            ,&
            kind=rp)
      !----------------------------
         pi_0(i,j) = cmplx( &
            0.0_rp &
            , &
            - (bhs*sy) &
            / (sqrt(2.0_rp)*(cl**4) + (bhs**2)*(cy**2)*(r**2)) &
            , &
            kind=rp)
      !----------------------------
         rh_0(i,j) = cmplx( &
            (-2.0_rp*((cl**6) - 2.0_rp*(cl**4)*bhm*r + (bhs**2)*(cl**2)*(r**2))) &
            / ((2.0_rp*(cl**4) + (bhs**2)*(r**2)*(1.0_rp + (cy**2) - (sy**2)))**2) &
            , & 
            (-2.0_rp*bhs*cy*r*((cl**4) - 2.0_rp*(cl**2)*bhm*r + (bhs**2)*(r**2))) &
            / ((2.0_rp*(cl**4) + (bhs**2)*(r**2)*(1.0_rp + (cy**2) - (sy**2)))**2) &
            , &
            kind=rp)
      !----------------------------
         ep_0(i,j) = cmplx( &
            ( &
               2*cl**4*bhm + bhs**2*cl**2*r*(-1 + cy**2 - sy**2) + bhs**2*bhm*r**2*(-1 - cy**2 + sy**2) &
            )/( &
               2*cl**4 + bhs**2*r**2*(1 + cy**2 - sy**2) &
            )**2 &
            , &
            ( &
               -2*bhs*cy*(cl**4 - 2*cl**2*bhm*r + bhs**2*r**2) &
            )/( &
               2*cl**4 + bhs**2*r**2*(1 + cy**2 - sy**2) &
            )**2 &
            , &
            kind=rp)
      !----------------------------
         thorn_prime_ta_0(i,j) = cmplx( &
            (4.0_rp*sqrt(2.0_rp)*(bhs**2)*cy*r*sy*(-6.0_rp*(cl**4) &
            + (bhs**2)*(r**2)*(1.0_rp + (cy**2) - (sy**2))) &
            ) &
            / ((2.0_rp*(cl**4) + (bhs**2)*(r**2)*(1.0_rp + (cy**2) - (sy**2)))**3) &
            , &
            (2.0_rp*sqrt(2.0_rp)*bhs*(cl**2)*sy*( &
              4.0_rp*(cl**4) &
            + 3.0_rp*(bhs**2)*(r**2)*(-1.0_rp - 3.0_rp*(cy**2) + (sy**2)) &
            )) &
            / ((2.0_rp*(cl**4) + (bhs**2)*(r**2)*(1.0_rp + (cy**2) - (sy**2)))**3) &
            , &
            kind=rp)
      !----------------------------
         psi2_0(i,j) = cmplx( &
            (4.0_rp*(cl**2)*bhm*( &
            - 2.0_rp*(cl**4) &
            + 3.0_rp*(bhs**2)*(r**2)*(1.0_rp + (cy**2) - (sy**2)) &
            )) &
            / ((2.0_rp*(cl**4) + (bhs**2)*(r**2)*(1.0_rp + (cy**2) - (sy**2)))**3) &
            , &
            (2.0_rp*bhs*cy*bhm*r*( &
            - 12.0_rp*(cl**4) &
            + (bhs**2)*(r**2)*(3.0_rp + (cy**2) - 3.0_rp*(sy**2)) &
            )) &
            / ((2.0_rp*(cl**4) + (bhs**2)*(r**2)*(1.0_rp + (cy**2) - (sy**2)))**3) &
            , &
            kind=rp)
      end do x_loop
      end do y_loop

   end subroutine bkgrd_np_init
!-----------------------------------------------------------------------------
end module mod_bkgrd_np
