program main
    use, intrinsic :: iso_fortran_env, only: prec=>real64

    implicit none

    complex(prec) :: a = (0,1)
    complex(prec) :: b = (0,0)

    print '(A)', 'starting'
    print *, a*a 
end program main
