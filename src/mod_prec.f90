!
! define precision for int, real (and complex) numbers
!
module mod_prec
!   if using fftw use these:
   use, intrinsic :: iso_c_binding, ip=>c_int, rp=>c_double
   
!   use, intrinsic :: iso_fortran_env, only: ip=>int64, rp=>real64
end module mod_prec
