module mod_cheb 
!-----------------------------------------------------------------------------
  use mod_prec

  use mod_sim_params, only: nx, ny, R 

  implicit none
!-----------------------------------------------------------------------------
  private

  character(:), allocatable :: dir

  ! Chebyshev matrix and
  ! Chebyshev differentiation matrix  
  real(rp), dimension(nx,nx) ::    cheb = 0
  real(rp), dimension(nx,nx) ::  D_cheb = 0
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
  subroutine cheb_set(fn, arr)
    character(*), allocatable,  intent(in)  :: fn
    real(rp), dimension(nx,nx), intent(out) :: arr

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
  end subroutine cheb_set
!-----------------------------------------------------------------------------
  subroutine cheb_der(vals,D_vals)
    complex(rp), dimension(nx,nx), intent(in)  :: vals
    complex(rp), dimension(nx,nx), intent(out) :: D_vals 
    integer(ip) :: i, j, k

    D_vals = 0
    do k=0,nx
      do j=1,ny
        do i=1,nx
          D_vals(i,j) = D_vals(i,j) + (D_cheb(i,k) * vals(k,j))
        end do
      end do
    end do

  end subroutine cheb_der 
!-----------------------------------------------------------------------------
end module mod_cheb

