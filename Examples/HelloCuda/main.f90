program testSaxpy
    use mathOps
    use cudafor
    implicit none
    integer, parameter :: N = 40000
    real :: x(N), y(N), a
    real, device :: x_d(N), y_d(N)
    type(dim3) :: grid, tBlock

    tBlock = dim3(256, 1, 1)
    grid = dim3(ceiling(real(N)/tBlock%x), 1, 1)

    x = 1.0; y = 2.0; a = 2.0
    x_d = x
    y_d = y
    call saxpy <<<grid, tBlock>>>(x_d, y_d, a)
    y = y_d
    write (*, *) 'Max error: ', maxval(abs(y - 4.0))
end program testSaxpy