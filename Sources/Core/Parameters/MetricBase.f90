module Metric
    use, intrinsic :: iso_fortran_env
    implicit none

    !  we will use MKS system

    type :: ReferenceState
        real(real64) :: ref_pressure, ref_temperature, ref_mach
        real(real64) :: ref_length, ref_gamma
        real(real64) :: ref_kinematic_viscosity
        real(real64), allocatable :: ref_mole_fraction(:)
    end type

contains

    subroutine set_metrics()
    end subroutine

    pure elemental function cm2meter(val) result(v)
        real(real64), parameter :: factor = 1.d0 / 100.d0
        real(real64), intent(in) :: val
        real(real64) :: v
        v = val * factor
    end function

    pure elemental function meter2cm(val) result(v)
        real(real64), parameter :: factor = 100.d0
        real(real64), intent(in) :: val
        real(real64) :: v
        v = val * factor
    end function

end module Metric
