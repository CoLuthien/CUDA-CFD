module Debug
    use, intrinsic :: iso_fortran_env
    implicit none

    contains 
    pure function check(input) result(res)
        class(*), allocatable, intent(in) :: input
        logical :: res
        res = allocated(input)
    end function    
end module