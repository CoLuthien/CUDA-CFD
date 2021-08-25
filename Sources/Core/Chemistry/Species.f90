module SpecieBase
    use, intrinsic :: iso_fortran_env
    use :: Metric
    implicit none
    type, abstract :: Specie
        real(real64) :: molar_weight, molar_weight_nd
        real(real64) :: molar_hof, mass_hof ! heat of formation
        real(real64) :: molar_diam !molecular diameter
        real(real64) :: char_temp ! charateristic temperature
        real(real64), allocatable, private :: f1(:), f2(:) ! by doi:10.2514/6.2013-303, pre-calculated weight for viscosity coefficient calculation
        real(real64), allocatable :: teab(:) ! unknown
        real(real64), allocatable :: diff_coef(:) ! precalculated part of diffusion coefficient eqn
    contains
        ! some  eqns are based on mass fraction, rather than mole fraction (ex. Fick's 1st law of diffusion)
        ! converting such eqns into mole fraction base are tooooo cumbersome,
        ! so we will provide a way of getting mass base property
        ! Todo: add entrophy, etc
        procedure(entalphy), deferred, pass :: H
        procedure(dynamic_viscosity), deferred, pass :: Eu
        procedure(conductivity), deferred, pass :: Cd
        !procedure(specific_heat), deferred, pass :: Cp

        procedure(entalphy), deferred, pass :: H_mass
        procedure(dynamic_viscosity), deferred, pass :: Eu_mass
        procedure(conductivity), deferred, pass :: Cd_mass
        procedure(specific_heat), deferred, pass :: Cp_mass

        procedure, pass, private :: calc_mixture_viscosity_scale_coef
        procedure, pass, private :: calc_mixture_diffusion_coef
        procedure, pass :: get_scale_factor

        procedure :: set_spc_data
        procedure :: init_spc_derived_data
    end type Specie

    type, extends(Specie) :: SpecieNasa7
        real(real64), dimension(7, 2) :: poly_cp, poly_h, poly_s ! mole based non-dimensinal poly-coefficients
        real(real64), dimension(5) :: poly_vis, poly_cd ! mole base non-dimensional poly

        real(real64), dimension(7, 2) :: poly_cp_m, poly_h_m, poly_s_m ! mole based non-dimensinal poly-coefficients
        real(real64), dimension(5) :: poly_vis_m, poly_cd_m ! mole base non-dimensional poly
    contains
        procedure :: H => entalphy7
        procedure :: Eu => viscosity5
        procedure :: Cd => conductivity5

        procedure :: H_mass => entalphy7
        procedure :: Eu_mass => viscosity5
        procedure :: Cd_mass => conductivity5
        procedure :: Cp_mass => cp_mass5
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
        pure elemental function specific_heat(self, temp, idx) result(C)
            import Specie, real64
            class(Specie), intent(in) :: self
            real(real64), intent(in) :: temp
            integer, intent(in) :: idx
            real(real64) :: C
        end function
        pure elemental function dynamic_viscosity(self, temp) result(H)
            import Specie, real64
            class(Specie), intent(in) :: self
            real(real64), intent(in) :: temp
            real(real64) :: H
        end function
        pure elemental function conductivity(self, temp) result(H)
            import Specie, real64
            class(Specie), intent(in) :: self
            real(real64), intent(in) :: temp
            real(real64) :: H
        end function

    end interface
contains

    pure function check_temp_range(temp) result(idx)
        real(real64), intent(in) :: temp
        integer :: idx
        idx = 1
        if (temp >= 1000.d0) then
            idx = 2
        end if
    end function

    pure elemental function entalphy7(self, temp, idx) result(H)
        class(SpecieNasa7), intent(in) :: self
        real(real64), intent(in) :: temp
        integer, intent(in) :: idx
        real(real64) :: H
        H = 0.d0
    end function entalphy7

    pure elemental function viscosity5(self, temp) result(Eu)
        !$omp declare simd(viscosity5) uniform(self, temp)
        class(SpecieNasa7), intent(in) :: self
        real(real64), intent(in) :: temp
        real(real64) :: Eu
        Eu = 0.d0
    end function viscosity5

    pure elemental function conductivity5(self, temp) result(Eu)
        !$omp declare simd(conductivity5) uniform(self, temp)
        class(SpecieNasa7), intent(in) :: self
        real(real64), intent(in) :: temp
        real(real64) :: Eu
        Eu = 0.d0
    end function

    pure elemental function cp_mass5(self, temp, idx) result(C)
        class(SpecieNasa7), intent(in) :: self
        real(real64), intent(in) :: temp
        integer, intent(in) :: idx
        real(real64) :: C
        C = 0.d0
    end function

    subroutine set_spc_data(self, n_spc)
        class(Specie), intent(inout) :: self
        integer, intent(in):: n_spc

        allocate (self%f1(n_spc), self%f2(n_spc), self%diff_coef(n_spc), self%teab(n_spc))
    end subroutine

    subroutine init_spc_derived_data(self, spcs, ref_state)
        class(Specie) :: self
        class(Specie), intent(in)::  spcs(:)
        type(ReferenceState), intent(in) :: ref_state
        call self%calc_mixture_viscosity_scale_coef(spcs)
        call self%calc_mixture_diffusion_coef(spcs, ref_state)
    end subroutine

    ! for convinient develop, spcs is in order, and include my self
    ! see third page or DOI: 10.2514/6.2013-303
    !
    subroutine calc_mixture_viscosity_scale_coef(self, spcs)
        class(Specie), intent(inout) :: self
        class(Specie), intent(in) :: spcs(:)
        real(real64), dimension(size(spcs)) :: factor1, factor2, alpha, beta
        ! 1.0 / sqrt(8 * (1 + mw(mine) / mw(others...))
        alpha(1:) = 1.d0 &
                    / sqrt(8.d0 * (1.d0 + self%molar_weight / spcs(:)%molar_weight)) ! ok non-d

        ! mw(others ...) / mw(mine)
        beta(1:) = spcs(:)%molar_weight / self%molar_weight ! ok non-d

        ! inherently non-d
        factor1(1:) = alpha * (beta**2.d0)
        factor2(1:) = (2.d0 * beta + 1.d0) * alpha

        self%f1 = factor1
        self%f2 = factor2
    end subroutine

    ! see third page or DOI: 10.2514/6.2013-303
    pure function get_scale_factor(self, ratio, mole_fraction) result(weight)
        !$omp declare simd(get_scale_factor) uniform(ratio, mole_fraction)
        class(Specie), intent(in) :: self
        real(real64), intent(in) :: ratio(1:), mole_fraction(1:)
        real(real64) :: weight
        ! now m_r * rest of term (in eqn 2)
        weight = sum(mole_fraction(1:) &
                     * (ratio(1:) * self%f1(1:) + sqrt(ratio(1:)) * self%f2(1:)))
    end function

    subroutine calc_mixture_diffusion_coef(self, spcs, ref_state)
        class(Specie) :: self
        class(Specie), intent(in) :: spcs(:)
        type(ReferenceState), intent(in) :: ref_state
        real(real64), dimension(size(spcs)) :: sigma_ab, teab, dcoef, euk
        real(real64) :: gt1, t3p, gp1
        real(real64) :: ref_gamma, ref_temp, ref_pressure, ref_kine_vis
        ! ref mixture is mole fraction
        ref_gamma = ref_state%ref_gamma
        ref_temp = ref_state%ref_temperature
        ref_pressure = ref_state%ref_pressure
        ref_kine_vis = ref_state%ref_kinematic_viscosity

        gt1 = ref_gamma * ref_temp
        gp1 = ref_gamma * ref_pressure
        t3p = (gt1**3) / (gp1**2)

        sigma_ab(1:) = 0.5d0 * (self%molar_diam * spcs(:)%molar_diam)
        teab(1:) = gt1 / sqrt(self%char_temp * spcs(1:)%char_temp) ! non-d
        dcoef(1:) = 0.0018583d0 &
                    * sqrt(t3p &
                           * ((1.d0 / self%molar_weight) + (1.d0 / spcs(1:)%molar_weight)))

        ! now non-d
        dcoef(1:) = (dcoef(1:) / sigma_ab(1:)) &
                    / ref_kine_vis

        self%diff_coef(1:) = dcoef(1:)
        self%teab = teab
    end subroutine

end module SpecieBase
