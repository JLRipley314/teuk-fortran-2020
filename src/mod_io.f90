module mod_io
   use mod_prec
   use mod_params, only: tables_dir

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
      rn = tables_dir // '/' // fn

      ! Note: here we ASSUME the input file is correctly formatted
      open(unit=uf,file=rn,status='old',action='read',iostat=ierror)
         if (ierror/=0) then
            write (*,*) "Error(read_arr): ierror=", ierror
            write (*,*) "file = ", rn
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
      integer(ip) :: i1, ierror
      integer(ip) :: uf = 3
      ! set the file name to read from
      rn = tables_dir // '/' // fn

      ! Note: here we ASSUME the input file is correctly formatted
      open(unit=uf,file=rn,status='old',action='read',iostat=ierror)
         if (ierror/=0) then
            write (*,*) "Error(read_arr): ierror=", ierror
            write (*,*) "file = ", rn
            stop
         end if
         do i1=1,n1
            read (uf,*,iostat=ierror) arr(i1,:)
         end do
      close(uf)
   end subroutine set_arr_2d
!-----------------------------------------------------------------------------
   subroutine set_arr_3d(fn, n1, n2, n3, arr)
      character(*), intent(in)  :: fn
      integer(ip),  intent(in)  :: n1, n2, n3
      real(rp),     intent(out) :: arr(n1,n2,n3)

      character(:), allocatable :: rn
      integer(ip) :: i1, i2, ierror
      integer(ip) :: uf = 3
      ! set the file name to read from
      rn = tables_dir // '/' // fn

      ! Note: here we ASSUME the input file is correctly formatted
      open(unit=uf,file=rn,status='old',action='read',iostat=ierror)
         if (ierror/=0) then
            write (*,*) "Error(read_arr): ierror=", ierror
            write (*,*) "file = ", rn
            stop
         end if
         do i1=1,n1
         do i2=1,n2
            read (uf,*,iostat=ierror) arr(i1,i2,:)
         end do   
         end do
      close(uf)
   end subroutine set_arr_3d
!-----------------------------------------------------------------------------
end module mod_io
