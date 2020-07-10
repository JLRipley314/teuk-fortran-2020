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
  subroutine cheb_init()
    integer(ip) :: i, j

    D_cheb(1,1) = (2.0_rp * (real(nx-1,rp)**2) + 1.0_rp) / 6.0_rp
    D_cheb(nx,nx) = - D_cheb(1,1)
    D_cheb(1,nx) = 0.5_rp * ((-1.0_rp) ** (nx-1))
    D_cheb(nx,1) = 0.5_rp * ((-1.0_rp) ** (nx-1))

    do j=1,nx
      do i=1,nx
        if (i/=j) then 
          if (i==1 .or. j==1) then
            D_cheb(i,j) = 2.0_rp*((-1.0_rp)**(i+j)) / (R(i)-R(j))
          else
            D_cheb(i,j) = 1.0_rp*((-1.0_rp)**(i+j)) / (R(i)-R(j))
          end if
        end if
      end do
    end do
    ! diagnoal 
    do i=2,nx-1
      cheb(i,i) = 0
      do j=2,nx-1
        if (j /= i) then
          cheb(i,i) = cheb(i,i) - cheb(i,j) 
        end if
      end do
    end do
  end subroutine cheb_init
!-----------------------------------------------------------------------------
  subroutine cheb_transform(vals,cheb_vals)
    complex(rp), dimension(nx,nx), intent(in)  :: vals
    complex(rp), dimension(nx,nx), intent(out) :: cheb_vals 
    integer(ip) :: i, j, k

    cheb_vals = 0
    do k=0,nx
      do j=1,ny
        do i=1,nx
          cheb_vals(i,j) = cheb_vals(i,j) + (cheb(i,k) * vals(k,j))
        end do
      end do
    end do

  end subroutine cheb_transform 
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

