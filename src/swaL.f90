module mod_swal 
!-----------------------------------------------------------------------------
  use mod_def_prec

  use mod_sim_params, only: ny, nl

  implicit none
!-----------------------------------------------------------------------------
  private

  character(:), allocatable :: dir
!-----------------------------------------------------------------------------
contains
  subroutine read_swaL(fn, swaL_arr)
    character(*), allocatable,  intent(in)  :: fn
    real(rp), dimension(ny,nl), intent(out) :: swaL_arr

    character(:), allocatable :: rn
    integer(ip) :: ierror
    integer(ip) :: uf = 3

    rn = dir // fn

    open(unit=uf,file=rn,status='old',action='read',iostat=ierror)
      if (ierror/=0) then
        write (*,*) "Error: ierror=", ierror
        write (*,*) "file = ", rn
        stop
      end if
      ! Note: here we ASSUME the input file is correctly formatted
      read (uf,*,iostat=ierror) swaL_arr
    close(uf)
  end subroutine
end module mod_swal
