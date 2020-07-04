program main
    use, intrinsic :: iso_fortran_env, only: prec=>real64

    use teuk

    implicit none

    complex(prec) :: a = (0,1)
    complex(prec) :: b = (2,1)

    print '(A)', 'starting'
    print *, a*a 
    print *, b+a 
end program main
