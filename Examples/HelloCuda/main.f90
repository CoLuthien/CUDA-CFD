program testSaxpy
    use :: Species
    use :: ArrayBase
    implicit none
    integer, parameter :: N = 40000
    real :: x(N), y(N), a
    real :: x_d(N), y_d(N)

    call test()

contains 
    subroutine test()
    implicit none
    type(Array4), allocatable :: arr1, arr3 
    type(Array4) :: arr2, arr4, arr5
    arr1 = Array(100, 100, 100, 200)
    arr2 = Array(100, 100, 100, 200)


    call move(arr1, arr4)
    call move(arr2, arr5)
    pause
    
    end subroutine test
    

end program testSaxpy
