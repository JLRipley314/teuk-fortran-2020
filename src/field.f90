module fields
    use, intrinsic :: iso_fortran_env, only: prec=>real64
    implicit none

    type field
        complex(prec), allocatable :: n(:,:)
        complex(prec), allocatable :: np1(:,:)
        complex(prec), allocatable :: l2(:,:)
        complex(prec), allocatable :: l3(:,:)
        complex(prec), allocatable :: l4(:,:)

        complex(prec), allocatable :: k1(:,:)
        complex(prec), allocatable :: k2(:,:)
        complex(prec), allocatable :: k3(:,:)
        complex(prec), allocatable :: k4(:,:)
        complex(prec), allocatable :: k5(:,:)
    end type field

end module fields
