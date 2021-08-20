module Constants
    use, intrinsic :: iso_fortran_env
    implicit none

    real(real64), parameter :: pi = 4.d0 * atan(1.d0)

    type :: InitialCondition(n_spc)
        integer, len :: n_spc
        real(real64) :: u, v, w, temp, tk, tw, tv
        real(real64) :: spcs_density(n_spc)
    end type

end module  