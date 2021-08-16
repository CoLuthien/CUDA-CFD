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
        procedure(kinematic_viscosity), deferred, pass :: Eu
        procedure(kinematic_conductivity), deferred, pass :: Cd
     end type SpecieBase

    type, extends(SpecieBase) :: SpecieNasa7
        real(real64), dimension(7, 2) :: poly_cp, poly_h, poly_s
        real(real64), dimension(5) :: poly_vis, poly_cd
    contains
        procedure :: H => entalphy7
        procedure :: Eu => viscosity5
        procedure :: Cd => conductivity5
    end type SpecieNasa7

    interface
        pure elemental function entalphy(self, temp, idx) result(H)
            import SpecieBase, real64
            class(SpecieBase), intent(in) :: self
            real(real64), intent(in) :: temp
            integer, intent(in) :: idx
            real(real64) :: H
        end function

    end interface

    interface 
        pure elemental function kinematic_viscosity(self, temp) result(H)
            import SpecieBase, real64
            class(SpecieBase), intent(in) :: self
            real(real64), intent(in) :: temp
            real(real64) :: H
        end function
        pure elemental function kinematic_conductivity(self, temp) result(H)
            import SpecieBase, real64
            class(SpecieBase), intent(in) :: self
            real(real64), intent(in) :: temp
            real(real64) :: H
        end function

    end interface 
contains

    pure elemental function entalphy7(self, temp, idx) result(H)
        class(SpecieNasa7), intent(in) :: self
        real(real64), intent(in) :: temp
        integer, intent(in) :: idx
        real(real64) :: H
        H = 0.d0
    end function entalphy7

    pure elemental function viscosity5(self, temp) result(Eu)
        class(SpecieNasa7), intent(in) :: self
        real(real64), intent(in) :: temp
        real(real64) :: Eu
        Eu = 0.d0
    end function viscosity5

    pure elemental function conductivity5(self, temp) result(Eu)
        class(SpecieNasa7), intent(in) :: self
        real(real64), intent(in) :: temp
        real(real64) :: Eu
        Eu = 0.d0
    end function 

end module Species
