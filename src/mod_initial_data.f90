module mod_initial_data
!=============================================================================
   use mod_prec
   use mod_field, only: field
   use mod_cheb, only: Rvec=>R
   use mod_swal, only: Yvec=>Y, swal
   use mod_params, only: &
      dt, nx, ny, &
      min_m, max_m, &
      spin=>psi_spin, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin, &
      R_max, pm_ang, &
      initial_data_direction, initial_data_type, &
      amp_nm, rl_nm, ru_nm, l_ang_nm, &
      amp_pm, rl_pm, ru_pm, l_ang_pm
!=============================================================================
   implicit none
   private

   public :: set_initial_data
!=============================================================================
contains
!=============================================================================
   subroutine set_initial_data(m_ang,p,q,f)
      integer(ip), intent(in) :: m_ang
      type(field), intent(inout) :: p, q, f

      integer(ip) :: i, j
      integer(ip) :: l_ang
      real(rp) :: max_val, bump, r
      real(rp) :: width, amp, rl, ru 

      bump    = 0.0_rp
      max_val = 0.0_rp
      amp     = 0.0_rp
      rl      = 0.0_rp
      ru      = 0.0_rp
      width   = 0.0_rp

      if (m_ang==pm_ang) then
         amp   = amp_pm
         ru    = ru_pm
         rl    = rl_pm
         l_ang = l_ang_pm

      else if (m_ang==-pm_ang) then
         amp   = amp_nm
         ru    = ru_nm
         rl    = rl_nm

      else
         write (*,*) "ERROR(set_initial_data): m_ang = ", m_ang
         stop
      end if

      width = ru-rl
!-----------------------------------------------------------------------------
      y_loop: do j=1,ny
      x_loop: do i=1,nx-1 ! index 'i=nx' is at 'r=infinity'
         r = (cl**2) / Rvec(i)

         if ((r<ru).and.(r>rl)) then
            bump = exp(-1.0_rp*width/(r-rl))*exp(-2.0_rp*width/(ru-r))
         else
            bump = 0.0_rp
         end if

         f%n(i,j,m_ang) = (((r-rl)/width)**2) * (((ru-r)/width)**2) * bump
         
         q%n(i,j,m_ang) = ( &
            (2.0_rp*(((r-rl)/width)   )*(((ru-r)/width)**2)) &
         -  (2.0_rp*(((r-rl)/width)**2)*( (ru-r)/width    )) &
         +  (1.0_rp*(1.0_rp           )*(((ru-r)/width)**2)) &
         -  (2.0_rp*(((r-rl)/width)**2)*(1.0_rp              )) &
         )*bump/width
         !--------------------------------------------------------------------
         ! rescale q as q = \partial_R f = -(r/cl)^2 partial_r f
         !--------------------------------------------------------------------
         q%n(i,j,m_ang) = q%n(i,j,m_ang)*(-(r**2)/(cl**2))

         p%n(i,j,m_ang) = 0.0_rp
                             
         f%n(i,j,m_ang) = f%n(i,j,m_ang) * swal(j,l_ang,m_ang,spin)
         q%n(i,j,m_ang) = q%n(i,j,m_ang) * swal(j,l_ang,m_ang,spin)

         if (abs(f%n(i,j,m_ang)) > max_val) then
            max_val = abs(f%n(i,j,m_ang))
         end if

      end do x_loop
      end do y_loop
!-----------------------------------------------------------------------------
! rescale to make max val = 'amp'
!-----------------------------------------------------------------------------
      f%n(:,:,m_ang) = f%n(:,:,m_ang) * (amp / max_val)
      q%n(:,:,m_ang) = q%n(:,:,m_ang) * (amp / max_val)
!-----------------------------------------------------------------------------
! copy to np1 so can be saved
!-----------------------------------------------------------------------------
      p%np1(:,:,m_ang) = p%n(:,:,m_ang) 
      q%np1(:,:,m_ang) = q%n(:,:,m_ang) 
      f%np1(:,:,m_ang) = f%n(:,:,m_ang) 

   end subroutine set_initial_data
!=============================================================================
end module mod_initial_data
