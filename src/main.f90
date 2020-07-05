program main
!-----------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: prec=>real64

  use mod_field
  use mod_sim_params
  use mod_teuk
  use mod_bkgrd

  implicit none
!-----------------------------------------------------------------------------
  print '(A)', 'starting'

  call init_bkgrd()
!-----------------------------------------------------------------------------
end program main
