module Species
    use, intrinsic :: iso_fortran_env
    implicit none(type, external)
    ! tot
    type, abstract :: SpecieBase
        real(real64) :: molar_weight
        real(real64) :: molar_hof, mass_hof ! heat of formation
        real(real64) :: molar_diam !molecular diameter
        real(real64) :: chrt_temp ! charateristic temperature
    contains
        procedure(entalphy), deferred, pass :: H
    end type SpecieBase
    type, extends(SpecieBase) :: SpecieNasa7
        real(real64), dimension(7, 2) :: poly_cp, poly_h, poly_s
    contains
        procedure :: H => entalphy7
    end type SpecieNasa7

    interface
        pure elemental function entalphy(self, temp) result(H)
            import SpecieBase, real64
            class(SpecieBase), intent(in) :: self
            real(real64), intent(in) :: temp
            real(real64) :: H
        end function
    end interface
contains

    pure elemental function entalphy7(self, temp) result(H)
        class(SpecieNasa7), intent(in) :: self
        real(real64), intent(in) :: temp
        real(real64) :: H
        H = 0.d0
    end function entalphy7

end module Species
