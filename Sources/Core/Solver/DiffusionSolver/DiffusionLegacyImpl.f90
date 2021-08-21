module DiffusionImpl
    use, intrinsic :: iso_fortran_env
    use :: DiffusionBase, only:DiffusionSolver, TransportProperty
    use :: ArrayBase
    use :: SpecieBase
    implicit none

    type, extends(DiffusionSolver) :: DiffusionLegacy
    contains
        !procedure :: solve_diffusion => diffusion
        procedure :: diffusive_flux => legacy_diffusive_flux
    end type DiffusionLegacy

    type, extends(TransportProperty) :: TransportLegacy
    contains
        procedure :: transport
    end type TransportLegacy

    interface DiffusionLegacy
        module procedure :: init_diffusion_legacy
    end interface

contains
    ! todo: refine init diffusion
    function init_diffusion_legacy(n_spc) result(solver)
        integer :: n_spc
        type(DiffusionLegacy) :: solver
        type(TransportLegacy) :: prop

        solver%m_transport = prop
    end function

    subroutine diffusion(self)
        class(DiffusionLegacy), intent(in) :: self
        call self%diffusive_flux()
    end subroutine

    subroutine legacy_diffusive_flux(self)
        class(DiffusionLegacy), intent(in) :: self
    end subroutine legacy_diffusive_flux

    ! every data in or out for  this routine, is non-dimensional quantity
    pure subroutine transport(self, spcs, spcs_density, temperature, pressure, density, tv & ! intent(in)
                              , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient & ! intent(out)
                              )
        use :: Constants, only:rprt, rsct
        class(TransportLegacy), intent(in) :: self
        class(Specie), intent(in) :: spcs(:)
        real(real64), intent(in) :: spcs_density(size(spcs))
        real(real64), intent(in) :: temperature, pressure, density, tv
        real(real64), intent(out) :: diffusion_coefficient(size(spcs))
        real(real64), intent(out) :: viscosity, thermal_conductivity, turbulent_viscosity
        real(real64) :: temp, eu, cd, weight, p, t3p, xsm, total, sum_cmk, cpt, cdt, difft
        real(real64) :: td(size(spcs)), omg(size(spcs)), df(size(spcs), size(spcs))
        real(real64), dimension(size(spcs)) ::  eu_k, cd_k, c, ratio, cp_k
        integer :: i, j, k, l, m, n_spc, nx, ny, nz, cnt, idx
        n_spc = size(spcs)
        ! loop for calculating laminar property
        ! this loop quite difficult to understand,
        ! to do: change loop range into realistic one
        ! calculatie mole fraction, non-dimensional
        temp = temperature
        c(1:n_spc) = spcs_density(1:n_spc) * spcs(1:n_spc)%molar_weight_nd
        c(1:n_spc) = c(1:n_spc) / sum(c(1:n_spc))

        eu_k(1:n_spc) = spcs(1:n_spc)%Eu(temp)! calculate single component, laminar viscosity coeff
        cd_k(1:n_spc) = spcs(1:n_spc)%Cd(temp) ! calculate single component, laminar thermal conductivity coeff

        do m = 1, n_spc
            ratio(1:n_spc) = eu_k(m) / eu_k(1:n_spc) ! get viscosity coeff ratio for this spc
            weight = spcs(m)%get_scale_factor(ratio(:), c(1:)) ! get weighting factor for this spc
            xsm = c(m) / weight
            eu = eu + xsm * eu_k(m)
            cd = cd + xsm * cd_k(m)
        end do

        viscosity = eu
        thermal_conductivity = cd
        ! todo : calculate diffusion coefficient and turbulent viscosity

        p = pressure
        c(1:n_spc) = spcs_density(1:n_spc) * spcs(1:n_spc)%molar_weight
        t3p = sqrt((temp**3) / (p**2))
        do m = 1, n_spc
            td(1:n_spc) = temp * spcs(m)%teab(1:)
            omg(1:n_spc) = 1.d0 / (td(1:n_spc)**0.145d0) &
                           + 1.d0 / ((td(1:n_spc) + 0.5d0)**2)
            df(m, 1:n_spc) = (spcs(m)%diff_coef(1:n_spc) * t3p) / omg(1:n_spc)
        end do

!---- Diff. Coef. Correction According to SAND86-8246 by Kee et al.
        sum_cmk = sum(c(:) * spcs(:)%molar_weight)
        do m = 1, n_spc
            total = 0.d0
            total = sum(c(1:n_spc) / df(m, 1:n_spc))
            total = total - c(m) / df(m, m)
            diffusion_coefficient(m) = (sum_cmk - c(m) * spcs(m)%molar_weight) / total
        end do

        cpt = 0.d0

        idx = check_temp_range(temp)
        cp_k(1:n_spc) = spcs(1:n_spc)%Cp_mass(temp, idx)
        cpt = sum(cp_k(1:n_spc) * spcs_density(1:n_spc) / density)

        turbulent_viscosity = eu + tv
        cdt = tv * cpt * rprt
        difft = (tv / density) * rsct

        cd = cd + cdt

        diffusion_coefficient(1:n_spc) = diffusion_coefficient(1:) + difft

    end subroutine transport

end module DiffusionImpl
