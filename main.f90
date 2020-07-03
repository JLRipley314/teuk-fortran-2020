program main
    use iso_fortran_env, only: double=>real64

    implicit none

    complex(double) :: a = (0,1)
    complex(double) :: b = (0,0)

    print '(A)', 'starting'
    print *, a*a 
end program main
