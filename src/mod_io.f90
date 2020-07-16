module mod_io
   use mod_prec
   use mod_params, only: nx, ny, tables_dir, output_dir
   use mod_field, only: field

   implicit none
!=============================================================================
   private
   public :: set_arr, write_csv
!=============================================================================
   interface set_arr
      module procedure set_arr_1d, set_arr_2d, set_arr_3d
   end interface

   interface write_csv
      module procedure write_field_csv, write_array_csv
   end interface
!=============================================================================
contains
!=============================================================================
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
!=============================================================================
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
!=============================================================================
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
!=============================================================================
! writes to one line, row by row
   subroutine write_field_csv(f)
      type(field),  intent(in) :: f

      character(:), allocatable :: rn
      logical :: exists
      integer(ip) :: i, j, ierror = 0
      integer(ip) :: uf = 3
      ! set the file name to read from
      rn = output_dir // '/' // f%name // '.csv'

      inquire(file=rn,exist=exists)
      if (exists) then
         open(unit=uf,file=rn,status='old',position='append',action='write',iostat=ierror)
      else
         open(unit=uf,file=rn,status='new',action='write',iostat=ierror) 
      end if
         do i=1,nx
            do j=1,ny
               write (uf,'(e12.6,a1)',advance='no',iostat=ierror) f%np1(i,j), ','
            end do
            ! line break
            write (uf,*) 
         end do
      close(uf)

      if (ierror/=0) then
         write (*,*) "Error(read_arr): ierror=", ierror
         write (*,*) "file = ", rn
         stop
      end if
   end subroutine write_field_csv
!=============================================================================
! writes to one line, row by row
   subroutine write_array_csv(fn, arr)
      character(*), intent(in) :: fn
      complex(rp),  dimension(nx,ny), intent(in) :: arr

      character(:), allocatable :: rn
      logical :: exists
      integer(ip) :: i, j, ierror = 0
      integer(ip) :: uf = 3
      ! set the file name to read from
      rn = output_dir // '/' // fn // '.csv'

      inquire(file=rn,exist=exists)
      if (exists) then
         open(unit=uf,file=rn,status='old',position='append',action='write',iostat=ierror)
      else
         open(unit=uf,file=rn,status='new',action='write',iostat=ierror) 
      end if
         do i=1,nx
            do j=1,ny
               write (uf,'(e12.6,a1)',advance='no',iostat=ierror) arr(i,j), ','
            end do
            write (uf,*) ! line break
         end do
      close(uf)

      if (ierror/=0) then
         write (*,*) "Error(read_arr): ierror=", ierror
         write (*,*) "file = ", rn
         stop
      end if
   end subroutine write_array_csv
!=============================================================================
end module mod_io
