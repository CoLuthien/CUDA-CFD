module mathOps
    implicit none
contains
    subroutine saxpy(x, y, a)
        real :: x(:), y(:)
        real, value :: a
        integer :: i, n
        n = size(x)
        i = blockDim%x*(blockIdx%x - 1) + threadIdx%x
        if (i <= n) y(i) = y(i) + a*x(i)
    end subroutine saxpy
end module mathOps
