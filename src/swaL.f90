module mod_swal 
!-----------------------------------------------------------------------------
  use mod_def_prec

  use mod_sim_params, only: nx, ny, nl, min_m, max_m, min_s, max_s

  implicit none
!-----------------------------------------------------------------------------
  private

  character(:), allocatable :: dir

  ! weights for Gaussian integration 
  real(rp), dimension(ny) :: weights

  ! p2 = 2, n1 = -1, etc. Refers to the spin weight 
  real(rp), dimension(ny, nl, min_m:max_m, min_s:max_s) :: swaL = 0
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
  subroutine set_vals(fn, arr)
    character(*), allocatable,  intent(in)  :: fn
    real(rp), intent(out) :: arr

    character(:), allocatable :: rn
    integer(ip) :: ierror
    integer(ip) :: uf = 3
    ! set the file name to read from
    rn = dir // fn

    ! Note: here we ASSUME the input file is correctly formatted
    open(unit=uf,file=rn,status='old',action='read',iostat=ierror)
      if (ierror/=0) then
        write (*,*) "Error: ierror=", ierror
        write (*,*) "file = ", rn
        stop
      end if
      read (uf,*,iostat=ierror) arr
    close(uf)
  end subroutine
!-----------------------------------------------------------------------------
  subroutine to_swaL_coef_space(spin,m_ang,vals,coefs)
    integer(ip), intent(in) :: spin
    integer(ip), intent(in) :: m_ang
    complex(rp), dimension(nx,ny), intent(in)  :: vals
    complex(rp), dimension(nx,nl), intent(out) :: coefs
    integer(ip) :: i, j, k

    coefs = 0


    if (m_ang<min_m .or. m_ang>max_m) then
      write (*,*) "ERROR(to_swaL_coef): m_ang vale"
      stop
    end if
    if (spin<min_s .or. spin>max_s) then
      write (*,*) "ERROR(to_swaL_coef): m_ang vale"
      stop
    end if

    ! Gaussian quadrature
      do k=1,nl
        do j=1,ny
          do i=1,nx
            coefs(i,k) = coefs(i,k) + (vals(i,j) * weights(j) * swaL(j,k,m_ang,spin))
          end do
        end do
      end do

  end subroutine
!-----------------------------------------------------------------------------
end module mod_swal
