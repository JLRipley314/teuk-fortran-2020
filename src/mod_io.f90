module mod_io
   use mod_prec
   use mod_params, only: dir_tables 

   implicit none
!-----------------------------------------------------------------------------
   private
   public :: set_arr
!-----------------------------------------------------------------------------
   interface set_arr
      module procedure set_arr_1d, set_arr_2d, set_arr_3d
   end interface
!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
   subroutine set_arr_1d(fn, n1, arr)
      character(*), intent(in)  :: fn
      integer(ip),  intent(in)  :: n1
      real(rp),     intent(out) :: arr(n1)

      character(:), allocatable :: rn
      integer(ip) :: ierror
      integer(ip) :: uf = 3
      ! set the file name to read from
      rn = dir_tables // fn

      ! Note: here we ASSUME the input file is correctly formatted
      open(unit=uf,file=fn,status='old',action='read',iostat=ierror)
         if (ierror/=0) then
            write (*,*) "Error(read_arr): ierror=", ierror
            write (*,*) "file = ", fn
            stop
         end if
         read (uf,*,iostat=ierror) arr
      close(uf)
   end subroutine set_arr_1d
!-----------------------------------------------------------------------------
   subroutine set_arr_2d(fn, n1, n2, arr)
      character(*), intent(in)  :: fn
      integer(ip),  intent(in)  :: n1, n2
      real(rp),     intent(out) :: arr(n1,n2)

      character(:), allocatable :: rn
      integer(ip) :: ierror
      integer(ip) :: uf = 3
      ! set the file name to read from
      rn = dir_tables // fn

      ! Note: here we ASSUME the input file is correctly formatted
      open(unit=uf,file=fn,status='old',action='read',iostat=ierror)
         if (ierror/=0) then
            write (*,*) "Error(read_arr): ierror=", ierror
            write (*,*) "file = ", fn
            stop
         end if
         read (uf,*,iostat=ierror) arr
      close(uf)
   end subroutine set_arr_2d
!-----------------------------------------------------------------------------
   subroutine set_arr_3d(fn, n1, n2, n3, arr)
      character(*), intent(in)  :: fn
      integer(ip),  intent(in)  :: n1, n2, n3
      real(rp),     intent(out) :: arr(n1,n2,n3)

      character(:), allocatable :: rn
      integer(ip) :: ierror
      integer(ip) :: uf = 3
      ! set the file name to read from
      rn = dir_tables // fn

      ! Note: here we ASSUME the input file is correctly formatted
      open(unit=uf,file=fn,status='old',action='read',iostat=ierror)
         if (ierror/=0) then
            write (*,*) "Error(read_arr): ierror=", ierror
            write (*,*) "file = ", fn
            stop
         end if
         read (uf,*,iostat=ierror) arr
      close(uf)
   end subroutine set_arr_3d
!-----------------------------------------------------------------------------
end module mod_io
