module ArrayBase
    use, intrinsic :: iso_fortran_env
    use :: Vector
    implicit none

    !todo, implement generic finalizing interface for array2~4

    type :: Array2
        integer :: m_size(2) ! a vector to represent array dimension
        real(real64), allocatable :: m_data(:, :)
        contains
    end type Array2

    type :: Array3
        integer :: m_size(3) ! a vector to represent array dimension
        real(real64), allocatable :: m_data(:, :, :)
        contains
        procedure :: get_data => get_data3
    end type Array3

    ! todo : make constructor of same shape
    type :: Array4
        integer :: m_size(4) ! a vector to represent array dimension
        real(real64), allocatable :: m_data(:, :, :, :)
    end type Array4


    interface move
        procedure :: move_array2
        procedure :: move_array3
        procedure :: move_array4
    end interface

    interface get_data
        procedure :: get_data3
    end interface

    interface Array 
        procedure :: make_array2
        procedure :: make_array3
        procedure :: make_array4
        procedure :: make_array4_size
    end interface Array

contains
    pure function make_array2(from) result(to)
        class(Array2), intent(in) :: from
        class(Array2), allocatable :: to
        allocate(to)
        to%m_size = from%m_size
        allocate(to%m_data, mold=from%m_data)
    end function make_array2 
    
    pure function make_array3(from) result(to)
        class(Array3), intent(in) :: from
        class(Array3), allocatable :: to
        allocate(to)
        to%m_size = from%m_size
        allocate(to%m_data, mold=from%m_data)
    end function make_array3

    pure function make_array4(from) result(to)
        class(Array4), intent(in) :: from
        class(Array4), allocatable :: to
        allocate(to)
        to%m_size = from%m_size
        allocate(to%m_data, mold=from%m_data)
    end function make_array4

    pure function make_array4_size(i, j, k, l) result(arr)
        integer, intent(in) :: i, j, k, l
        class(Array4), allocatable :: arr
        allocate(arr)
        arr%m_size = [i, j, k, l]
        allocate(arr%m_data(i, j, k, l))
    end function     

    pure subroutine move_array2(to, from) 
        class(Array2), allocatable, intent(inout) :: from
        type(Array2), intent(inout) :: to
        to%m_size = from%m_size
        if (allocated(to%m_data)) then
            deallocate(to%m_data)
        end if 
        call move_alloc(from%m_data, to%m_data)
        if (allocated(from)) then
            deallocate(from)
        end if 
    end subroutine move_array2
    pure subroutine move_array3(to, from) 
        class(Array3), allocatable, intent(inout) :: from
        type(Array3), intent(inout) :: to
        to%m_size = from%m_size
        if (allocated(to%m_data)) then
            deallocate(to%m_data)
        end if 
        call move_alloc(from%m_data, to%m_data)
        if (allocated(from)) then
            deallocate(from)
        end if 
    end subroutine move_array3

    pure subroutine move_array4(to, from) 
        type(Array4), allocatable, intent(inout) :: from
        type(Array4), intent(inout) :: to
        to%m_size = from%m_size
        if (allocated(to%m_data)) then
            deallocate(to%m_data)
        end if 

        call move_alloc(from%m_data, to%m_data)

        if (allocated(from)) then
            deallocate(from)
        end if 
    end subroutine move_array4

    pure function get_data3(self, i, j, k) result(val)
        class(Array3), intent(in) :: self
        integer, intent(in) :: i, j, k
        real(real64) :: val
        val = self%m_data(i, j, k)
    end function
end module ArrayBase
