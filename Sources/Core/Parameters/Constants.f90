module Constants
    use, intrinsic :: iso_fortran_env
    implicit none
    public 
    integer, parameter :: n_face = 6
    real(real64), parameter :: pi = 4.d0 * atan(1.d0)

    real(real64) :: prandtl_number = 0.9d0, schmidt_number = 0.9d0

    type :: InitialCondition
        real(real64) :: u, v, w, temp, tk, tw, tv
        real(real64), allocatable :: spcs_density(:)
    end type


end module
