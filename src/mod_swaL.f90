module mod_swal 
!-----------------------------------------------------------------------------
  use mod_prec

  use mod_sim_params, only: nx, ny, lmax, min_m, max_m, min_s, max_s

  implicit none
!-----------------------------------------------------------------------------
  private

  character(:), allocatable :: dir

  ! weights for Gaussian integration 
  real(rp), dimension(ny) :: weights

  ! p2 = 2, n1 = -1, etc. Refers to the spin weight 
  real(rp), dimension(ny, 0:lmax, min_m:max_m, min_s:max_s) ::  swal = 0
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
  end subroutine set_vals
!-----------------------------------------------------------------------------
  subroutine swal_real_to_coef(spin,m_ang,vals,coefs)
    integer(ip), intent(in) :: spin
    integer(ip), intent(in) :: m_ang
    complex(rp), dimension(nx,ny),     intent(in)  :: vals
    complex(rp), dimension(nx,0:lmax), intent(out) :: coefs
    integer(ip) :: i, j, k

    if (m_ang<min_m .or. m_ang>max_m) then
      write (*,*) "ERROR(to_swal_coef): m_ang = ", m_ang
      stop
    end if
    if (spin<min_s .or. spin>max_s) then
      write (*,*) "ERROR(to_swal_coef): spin = ", spin
      stop
    end if

    ! Gaussian quadrature
    coefs = 0
    do k=0,lmax
      do j=1,ny
        do i=1,nx
          coefs(i,k) = coefs(i,k) + (vals(i,j) * weights(j) * swal(j,k,m_ang,spin))
        end do
      end do
    end do

  end subroutine swal_real_to_coef
!-----------------------------------------------------------------------------
  subroutine swal_coef_to_real(spin,m_ang,coefs,vals)
    integer(ip), intent(in) :: spin
    integer(ip), intent(in) :: m_ang
    complex(rp), dimension(nx,0:lmax), intent(in)  :: coefs
    complex(rp), dimension(nx,ny),     intent(out) :: vals
    integer(ip) :: i, j, k

    if (m_ang<min_m .or. m_ang>max_m) then
      write (*,*) "ERROR(to_swal_coef): m_ang = ", m_ang
      stop
    end if
    if (spin<min_s .or. spin>max_s) then
      write (*,*) "ERROR(to_swal_coef): spin = ", spin
      stop
    end if

    ! synthesis
    vals = 0
    do k=0,lmax
      do j=1,ny
        do i=1,nx
          vals(i,j) = vals(i,j) + (coefs(i,k) * swal(j,k,m_ang,spin))
        end do
      end do
    end do

  end subroutine swal_coef_to_real
!-----------------------------------------------------------------------------
  subroutine swal_lower(spin,m_ang,coefs,vals)
    integer(ip), intent(in) :: spin
    integer(ip), intent(in) :: m_ang
    complex(rp), dimension(nx,0:lmax), intent(inout) :: coefs
    complex(rp), dimension(nx,ny),     intent(inout) :: vals

    real(rp)    :: pre
    integer(ip) :: i, k

    call swal_real_to_coef(spin,m_ang,vals,coefs) 

    do k=0,lmax
      do i=1,nx
        pre = -sqrt((real(k,rp)+real(spin,rp))*(real(k,rp)-real(spin,rp)+1.0_rp)) 
        coefs(i,k) = pre*coefs(i,k)
      end do
    end do

    call swal_coef_to_real(spin-1,m_ang,coefs,vals) 

  end subroutine swal_lower
!-----------------------------------------------------------------------------
  subroutine swal_raise(spin,m_ang,coefs,vals)
    integer(ip), intent(in) :: spin
    integer(ip), intent(in) :: m_ang
    complex(rp), dimension(nx,0:lmax), intent(inout) :: coefs
    complex(rp), dimension(nx,ny),     intent(inout) :: vals

    real(rp)    :: pre
    integer(ip) :: i, k

    call swal_real_to_coef(spin,m_ang,vals,coefs) 

    do k=0,lmax
      do i=1,nx
        pre = sqrt((real(k,rp)-real(spin,rp))*(real(k,rp)+real(spin,rp)+1.0_rp)) 
        coefs(i,k) = pre*coefs(i,k)
      end do
    end do

    call swal_coef_to_real(spin+1,m_ang,coefs,vals) 

  end subroutine swal_raise
!-----------------------------------------------------------------------------
end module mod_swal
