module LegacyTransport
    use, intrinsic :: iso_fortran_env
    use :: SpecieBase
    implicit none

contains

    pure subroutine mixture_properties(spcs, mole_fraction, temperature& ! intent(in)
                                   , viscosity, thermal_conductivity, n_spc & ! intent(out)
                                   )
        use :: Constants, only:prandtl_number, schmidt_number
        class(Specie), intent(in) :: spcs(:)
        real(real64), intent(in) :: mole_fraction(:)
        real(real64), intent(in) :: temperature
        real(real64), intent(out) :: viscosity, thermal_conductivity
        integer, intent(in) :: n_spc
        real(real64) :: cd, eu, xsm, eu_k(n_spc), cd_k(n_spc), ratio(n_spc), weight
        integer :: m
        ! loop for calculating laminar property
        ! this loop quite difficult to understand,
        ! calculatie mole fraction, non-dimensional

        eu_k(1:n_spc) = spcs(1:n_spc)%Eu(temperature)! calculate single component, laminar viscosity coeff
        cd_k(1:n_spc) = spcs(1:n_spc)%Cd(temperature) ! calculate single component, laminar thermal conductivity coeff

        do m = 1, n_spc
            ratio(1:n_spc) = eu_k(m) / eu_k(1:n_spc) ! get viscosity coeff ratio for this spc
            weight = spcs(m)%get_scale_factor(ratio(1:n_spc), mole_fraction(1:n_spc)) ! get weighting factor for this spc
            xsm = mole_fraction(m) / weight
            eu = eu + xsm * eu_k(m)
            cd = cd + xsm * cd_k(m)
        end do

        viscosity = eu
        thermal_conductivity = cd
    end subroutine

    pure subroutine mixture_diffusivity(spcs, mole_fraction, mass_fraction, temperature, pressure,  &
                                            diffusion_coefficient, n_spc)
        use :: Constants, only:prandtl_number, schmidt_number
        class(Specie), intent(in) :: spcs(:)
        real(real64), intent(in) :: mole_fraction(n_spc), mass_fraction(n_spc), temperature, pressure
        real(real64), intent(out) :: diffusion_coefficient(:)
        integer, intent(in) :: n_spc
        real(real64) :: total, sum_cmk, t3p, td(n_spc), omg(n_spc), diffusion_coef(n_spc, n_spc)
        real(real64) :: cpt, cdt, cp_k(n_spc), difft
        integer :: m, idx

        t3p = sqrt((temperature**3) / (pressure**2))

        ! Diffusion Coef. calculation based on Chapmann-Enskog theory
        do m = 1, n_spc
            !DIR$ NOUNROLL
            td(1:n_spc) = temperature * spcs(m)%teab(1:n_spc)
            omg(1:n_spc) = 1.d0 / (td(1:n_spc)**0.145d0) &
                           + 1.d0 / ((td(1:n_spc) + 0.5d0)**2)
            diffusion_coef(m, 1:n_spc) = (spcs(m)%diff_coef(1:n_spc) * t3p) / omg(1:n_spc)
        end do

        !---- Diff. Coef. Correction According to SAND86-8246 by Kee et al.
        sum_cmk = sum(mole_fraction(:) * spcs(:)%molar_weight_nd)
        difft = (eddy_viscosity / density) * (1.d0 / schmidt_number)
        do m = 1, n_spc
            total = 0.d0

            total = sum(mole_fraction(1:n_spc) / diffusion_coef(m, 1:n_spc))
            total = total - mole_fraction(m) / diffusion_coef(m, m)

            diffusion_coefficient(m) = ((1.d0 - mass_fraction(m)) / total) + difft
        end do

    end subroutine 

    pure subroutine mixture_turbulent_properties(spcs, mole_fraction, mass_fraction, temperature, density, &
                                             eddy_viscosity, eddy_conductivity, diff_turbulent, n_spc)
        use :: Constants, only:prandtl_number, schmidt_number
        class(Specie), intent(in) :: spcs(:)
        real(real64), intent(in) :: mole_fraction(n_spc), mass_fraction(n_spc), temperature, density,eddy_viscosity
        real(real64), intent(out) :: eddy_conductivity, diff_turbulent
        integer, intent(in) :: n_spc
        real(real64) :: total, sum_cmk, t3p, td(n_spc), omg(n_spc), diffusion_coef(n_spc, n_spc)
        real(real64) :: cpt, cdt, cp_k(n_spc)
        integer :: m, idx

        cpt = 0.d0

        idx = check_temp_range(temperature)

        cp_k(1:n_spc) = spcs(1:n_spc)%Cp_mass(temperature, idx)

        cpt = sum(cp_k(1:n_spc) * mass_fraction(1:n_spc))

        cdt = eddy_viscosity * cpt * (1.d0 / prandtl_number)
        eddy_conductivity =  cdt



    end subroutine


end module 