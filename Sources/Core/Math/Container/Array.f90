module ArrayBase
    use, intrinsic :: iso_fortran_env
    use :: Vector
    implicit none

    type, abstract :: ArrayClass

    contains
        !Todo : interface or generic arrays
    end type ArrayClass

    type, extends(ArrayClass) :: Array2
        integer :: m_size(2) ! a vector to represent array dimension
        real(real64), allocatable :: m_data(:, :)
    end type Array2

    type, extends(ArrayClass) :: Array3
        integer :: m_size(3) ! a vector to represent array dimension
        real(real64), allocatable :: m_data(:, :, :)
    end type Array3

    type, extends(ArrayClass) :: Array4
        integer :: m_size(4) ! a vector to represent array dimension
        real(real64), allocatable :: m_data(:, :, :, :)
    end type Array4

end module ArrayBase
