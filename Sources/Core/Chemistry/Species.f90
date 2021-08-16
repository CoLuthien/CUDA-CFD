module SpecieBase
    use, intrinsic :: iso_fortran_env
    implicit none(type, external)
    type, abstract :: Specie
        real(real64) :: molar_weight
        real(real64) :: molar_hof, mass_hof ! heat of formation
        real(real64) :: molar_diam !molecular diameter
        real(real64) :: chrt_temp ! charateristic temperature
        real(real64), allocatable :: f1(:), f2(:) ! by doi:10.2514/6.2013-303, calculating scaling factor for ...
    contains
        procedure(entalphy), deferred, pass :: H
        procedure(kinematic_viscosity), deferred, pass :: Eu
        procedure(kinematic_conductivity), deferred, pass :: Cd
        procedure, pass, private :: calc_mixture_scale_coef
        procedure, pass :: get_scale_factor
    end type Specie

    type, extends(Specie) :: SpecieNasa7
        real(real64), dimension(7, 2) :: poly_cp, poly_h, poly_s
        real(real64), dimension(5) :: poly_vis, poly_cd
    contains
        procedure :: H => entalphy7
        procedure :: Eu => viscosity5
        procedure :: Cd => conductivity5
    end type SpecieNasa7

    interface
        pure elemental function entalphy(self, temp, idx) result(H)
            import Specie, real64
            class(Specie), intent(in) :: self
            real(real64), intent(in) :: temp
            integer, intent(in) :: idx
            real(real64) :: H
        end function

    end interface

    interface
        pure elemental function kinematic_viscosity(self, temp) result(H)
            import Specie, real64
            class(Specie), intent(in) :: self
            real(real64), intent(in) :: temp
            real(real64) :: H
        end function
        pure elemental function kinematic_conductivity(self, temp) result(H)
            import Specie, real64
            class(Specie), intent(in) :: self
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

    ! for convinient develop, spcs is in order, and include my self
    ! see third page or DOI: 10.2514/6.2013-303
    !
    subroutine calc_mixture_scale_coef(self, spcs)
        class(Specie), intent(inout) :: self
        class(Specie), intent(in) :: spcs(:)
        real(real64), dimension(size(spcs)) :: factor1, factor2, alpha, beta
        ! 1.0 / sqrt(8 * (1 + mw(mine) / mw(others...))
        alpha(:) = 1.d0 &
                   / sqrt(8.d0 * (1.d0 + self%molar_weight / spcs(:)%molar_weight))

        ! mw(others ...) / mw(mine)
        beta(:) = spcs(:)%molar_weight / self%molar_weight

        factor1(:) = alpha * (beta**2.d0)
        factor2(:) = (2.d0 * beta + 1.d0) * alpha

        self%f1 = factor1
        self%f2 = factor2
    end subroutine

    ! see third page or DOI: 10.2514/6.2013-303
    pure function get_scale_factor(self, ratio, molar_concent) result(weight)
        class(Specie), intent(in) :: self
        real(real64), intent(in) :: ratio(:), molar_concent(:)
        real(real64) :: weight
        ! now m_r * rest of term (in eqn 2)
        weight = sum(molar_concent(1:) & 
            * (ratio(1:) * self%f1(1:) + sqrt(ratio) * self%f2(1)))
    end function

end module SpecieBase
