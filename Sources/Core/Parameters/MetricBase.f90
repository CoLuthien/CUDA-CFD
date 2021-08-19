module MetricBase
    use, intrinsic :: iso_fortran_env
    implicit none

    !  we will use MKS system
    type :: Metric
        real(real64) :: ref_pressure, ref_temperature, ref_mach
        real(real64) :: ref_length
    contains
    end type

contains

    pure elemental function cm2meter(val) result(v)
        real(real64), parameter :: factor = 1.d0/100.d0
        real(real64), intent(in) :: val
        real(real64) :: v
        v = val*factor
    end function

    pure elemental function meter2cm(val) result(v)
        real(real64), parameter :: factor = 100.d0
        real(real64), intent(in) :: val
        real(real64) :: v
        v = val*factor
    end function

end module MetricBase
