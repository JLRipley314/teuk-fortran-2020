module mod_initial_data
!=============================================================================
   use mod_prec
   use mod_field, only: field
   use mod_cheb, only: Rvec=>R
   use mod_swal, only: Yvec=>Y, swal
   use mod_params, only: &
      dt, nx, ny, &
      spin, min_m, max_m, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin, &
      R_max, pm_ang, &
      initial_data_direction, initial_data_type, &
      amp_pm, rl_pm, ru_pm, l_ang_pm
!=============================================================================
   implicit none
   private

   public :: set_initial_data
!=============================================================================
contains
!=============================================================================

!=============================================================================
   subroutine set_initial_data(p,q,f)
      type(field), intent(inout) :: p, q, f

      integer(ip) :: i, j
      real(rp) :: max_val, bump, r
      real(rp) :: width 

      bump    = 0.0_rp
      max_val = 0.0_rp
      width = ru_pm-rl_pm
!-----------------------------------------------------------------------------
      y_loop: do j=1,ny
      x_loop: do i=1,nx-1 ! last index is at 'r=infinity'
         r = (cl**2) / Rvec(i)

         if ((r<ru_pm).and.(r>rl_pm)) then
            bump = exp(-1.0_rp*width/(r-rl_pm))*exp(-2.0_rp*width/(ru_pm-r))
         else
            bump = 0.0_rp
         end if

         f%n(i,j) = (((r-rl_pm)/width)**2) * (((ru_pm-r)/width)**2) * bump
         
         q%n(i,j) = ((2.0_rp*(((r-rl_pm)/width)   )*(((ru_pm-r)/width)**2)) &
                  -  (2.0_rp*(((r-rl_pm)/width)**2)*( (ru_pm-r)/width    )) &
                  +  (1.0_rp*(1.0_rp              )*(((ru_pm-r)/width)**2)) &
                  -  (2.0_rp*(((r-rl_pm)/width)**2)*(1.0_rp              )) &
                  )*bump/width
      ! rescale q as q = \partial_R f = -(r/cl)^2 partial_r f
         q%n(i,j) = q%n(i,j)*(-(r**2)/(cl**2))

         p%n(i,j) = 0.0_rp
                             
         f%n(i,j) = f%n(i,j) * swal(j,l_ang_pm,pm_ang,spin)
         q%n(i,j) = q%n(i,j) * swal(j,l_ang_pm,pm_ang,spin)

         if (abs(f%n(i,j)) > max_val) then
            max_val = abs(f%n(i,j))
         end if

      end do x_loop
      end do y_loop
!-----------------------------------------------------------------------------
! rescale to make max val = 'amp'
!-----------------------------------------------------------------------------
      f%n = f%n * (amp_pm / max_val)
      q%n = q%n * (amp_pm / max_val)
!-----------------------------------------------------------------------------
! copy to np1 so can be saved
!-----------------------------------------------------------------------------
      p%np1 = p%n
      q%np1 = q%n
      f%np1 = f%n

   end subroutine set_initial_data
!=============================================================================
end module mod_initial_data
