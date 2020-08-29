module mod_io
   use mod_prec
   use mod_params, only: nx, ny, nl, tables_dir, output_dir
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
      module procedure write_field_csv, write_array_1d_csv, write_array_2d_csv
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
      integer(ip) :: uf
      ! set the file fname to read from
      rn = tables_dir // '/' // fn

      ! Note: here we ASSUME the input file is correctly formatted
      open(newunit=uf,file=rn,status='old',action='read',iostat=ierror)
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
      integer(ip) :: uf
      ! set the file fname to read from
      rn = tables_dir // '/' // fn

      ! Note: here we ASSUME the input file is correctly formatted
      open(newunit=uf,file=rn,status='old',action='read',iostat=ierror)
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
      integer(ip) :: uf
      ! set the file fname to read from
      rn = tables_dir // '/' // fn

      ! Note: here we ASSUME the input file is correctly formatted
      open(newunit=uf,file=rn,status='old',action='read',iostat=ierror)
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
   subroutine write_array_1d_csv(fn, time, m_ang, arr)
      character(*),              intent(in) :: fn
      real(rp),                  intent(in) :: time
      integer(ip),               intent(in) :: m_ang
      complex(rp), dimension(:), intent(in) :: arr

      character(:), allocatable  :: mstr, fn_re, fn_im
      integer(ip) :: ubx, lbx 
      logical :: exists
      integer(ip) :: i, ierror = 0
      integer(ip) :: uf

      lbx = lbound(arr,1)
      ubx = ubound(arr,1)

      ! inelegant int to str conversion
      mstr = '     '
      write (mstr,'(i5)') m_ang
      mstr = trim(adjustl(mstr))
      ! set the file fname to read from
      fn_re = output_dir // '/' // fn // '_m' // mstr // '_re.csv'
      fn_im = output_dir // '/' // fn // '_m' // mstr // '_im.csv'
      !----------------------------------------------------------------------
      ! save real part 
      !----------------------------------------------------------------------
      inquire(file=fn_re,exist=exists)
      if (exists) then
         open(newunit=uf,file=fn_re,status='old',position='append',action='write',iostat=ierror)
      else
         open(newunit=uf,file=fn_re,status='new',action='write',iostat=ierror) 
      end if

      write (uf,'(e14.6,a1,i3,a1)',advance='no',iostat=ierror) &
         time, ',', ubx-lbx+1, ','

      do i=lbx,ubx
         write (uf,'(e14.6,a1)',advance='no',iostat=ierror) &
            real(arr(i),kind=rp), ','
      end do
      ! line break 
      write (uf,*)

      close(uf)

      if (ierror/=0) then
         write (*,*) "Error(read_arr): ierror=", ierror
         write (*,*) "file = ", fn_re
         stop
      end if
      !----------------------------------------------------------------------
      ! save imaginary part 
      !----------------------------------------------------------------------
      inquire(file=fn_im,exist=exists)
      if (exists) then
         open(newunit=uf,file=fn_im,status='old',position='append',action='write',iostat=ierror)
      else
         open(newunit=uf,file=fn_im,status='new',action='write',iostat=ierror) 
      end if

      write (uf,'(e14.6,a1,i3,a1)',advance='no',iostat=ierror) &
         time, ',', ubx-lbx+1, ','

      do i=lbx,ubx
         write (uf,'(e14.6,a1)',advance='no',iostat=ierror) &
            aimag(arr(i)), ','
      end do
      ! line break      
      write (uf,*) 

      close(uf)

      if (ierror/=0) then
         write (*,*) "Error(read_arr): ierror=", ierror
         write (*,*) "file = ", fn_im
         stop
      end if
   end subroutine write_array_1d_csv
!=============================================================================
! writes to one line, row by row
   subroutine write_array_2d_csv(fn, time, m_ang, arr)
      character(*),                intent(in) :: fn
      real(rp),                    intent(in) :: time
      integer(ip),                 intent(in) :: m_ang
      complex(rp), dimension(:,:), intent(in) :: arr

      character(:), allocatable  :: mstr, fn_re, fn_im
      integer(ip) :: ubx, uby, lbx, lby 
      logical :: exists
      integer(ip) :: i, j, ierror = 0
      integer(ip) :: uf

      lbx = lbound(arr,1)
      lby = lbound(arr,2)

      ubx = ubound(arr,1)
      uby = ubound(arr,2)

      ! inelegant int to str conversion
      mstr = '     '
      write (mstr,'(i5)') m_ang
      mstr = trim(adjustl(mstr))
      ! set the file fname to read from
      fn_re = output_dir // '/' // fn // '_m' // mstr // '_re.csv'
      fn_im = output_dir // '/' // fn // '_m' // mstr // '_im.csv'
      !----------------------------------------------------------------------
      ! save real part 
      !----------------------------------------------------------------------
      inquire(file=fn_re,exist=exists)
      if (exists) then
         open(newunit=uf,file=fn_re,status='old',position='append',action='write',iostat=ierror)
      else
         open(newunit=uf,file=fn_re,status='new',action='write',iostat=ierror) 
      end if

      write (uf,'(e14.6,a1,i3,a1,i3,a1)',advance='no',iostat=ierror) &
         time, ',', ubx-lbx+1, ',', uby-lby+1, ','

      do i=lbx,ubx
      do j=lby,uby
         write (uf,'(e14.6,a1)',advance='no',iostat=ierror) &
            real(arr(i,j),kind=rp), ','
      end do
      end do
      ! line break 
      write (uf,*)

      close(uf)

      if (ierror/=0) then
         write (*,*) "Error(read_arr): ierror=", ierror
         write (*,*) "file = ", fn_re
         stop
      end if
      !----------------------------------------------------------------------
      ! save imaginary part 
      !----------------------------------------------------------------------
      inquire(file=fn_im,exist=exists)
      if (exists) then
         open(newunit=uf,file=fn_im,status='old',position='append',action='write',iostat=ierror)
      else
         open(newunit=uf,file=fn_im,status='new',action='write',iostat=ierror) 
      end if

      write (uf,'(e14.6,a1,i3,a1,i3,a1)',advance='no',iostat=ierror) &
         time, ',', ubx-lbx+1, ',', uby-lby+1, ','

      do i=lbx,ubx
      do j=lby,uby
         write (uf,'(e14.6,a1)',advance='no',iostat=ierror) &
            aimag(arr(i,j)), ','
      end do
      end do
      ! line break      
      write (uf,*) 

      close(uf)

      if (ierror/=0) then
         write (*,*) "Error(read_arr): ierror=", ierror
         write (*,*) "file = ", fn_im
         stop
      end if
   end subroutine write_array_2d_csv
!=============================================================================
! writes to one line, row by row
   subroutine write_field_csv(time,m_ang,f)
      real(rp),    intent(in) :: time
      integer(ip), intent(in) :: m_ang
      type(field), intent(in) :: f

      call write_array_2d_csv(f%fname, time, m_ang, f%np1(:,:,m_ang))

   end subroutine write_field_csv
!=============================================================================
end module mod_io
