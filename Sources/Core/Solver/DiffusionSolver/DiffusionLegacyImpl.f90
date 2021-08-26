module DiffusionImpl
    use, intrinsic :: iso_fortran_env
    use :: DiffusionBase
    use :: ArrayBase
    use :: SpecieBase
    implicit none

    type, extends(DiffusionSolver) :: DiffusionLegacy
    contains
        !procedure :: solve_diffusion => diffusion

        procedure, non_overridable :: diffusive_mass => diffusive_mass_impl
        procedure, non_overridable :: viscous_stress => viscous_stress_impl
    end type DiffusionLegacy

    type, extends(TransportProperty) :: TransportLegacy
    contains
        procedure :: transport => transport_impl
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

    pure subroutine viscous_stress_impl(self, du, dv, dw, turbulent_viscosity, i, j, k, stress)
        use :: ArrayBase
        use :: Vector
        class(DiffusionLegacy), intent(in) :: self
        type(Vector3), intent(in) :: du, dv, dw ! gradient representation of physical grid system
        type(ShearStress), intent(out) :: stress
        real(real64), intent(in) :: turbulent_viscosity
        integer, intent(in) :: i, j, k
        real(real64), parameter :: ort = 1.d0 / 3.d0

        stress%tau_xx = turbulent_viscosity * (4.d0 * du%x - 2.d0 * dv%y - 2.d0 * dw%z) * ort
        stress%tau_yy = turbulent_viscosity * (4.d0 * dv%y - 2.d0 * dw%z - 2.d0 * du%x) * ort
        stress%tau_zz = turbulent_viscosity * (4.d0 * dw%z - 2.d0 * du%x - 2.d0 * dv%y) * ort

        stress%tau_xy = turbulent_viscosity * (du%y + dv%x)
        stress%tau_yz = turbulent_viscosity * (dv%z + dw%y)
        stress%tau_xz = turbulent_viscosity * (dw%x + du%z)

    end subroutine

    pure subroutine diffusive_mass_impl(self, metrics, diffusion_coef, spcs_density, density, n_spc, i, j, k, diffused_mass)
        use :: Vector
        use :: ArrayBase
        class(DiffusionLegacy), intent(in) :: self
        class(CellMetricData3D), intent(in) :: metrics
        class(Array4), intent(in) :: spcs_density
        class(Array3), intent(in) :: density
        real(real64), intent(in) :: diffusion_coef(n_spc)
        integer, intent(in) :: n_spc, i, j, k
        type(Vector3) :: sec_x, sec_y, sec_z
        type(Vector3) :: s, e, c
        real(real64), dimension(n_spc) :: rear, front, left, right, top, bottom ! mass fraction of each cell
        real(real64), dimension(n_spc) :: dyk_ds, dyk_de, dyk_dc ! differenciated value for each computation grid direction
        real(real64), dimension(n_spc) :: diff_density
        real(real64), dimension(n_spc) :: flux_x, flux_y, flux_z
        real(real64), intent(out) :: diffused_mass(n_spc, 3) ! final result of calculated mass diffusion
        type(Vector3) :: sec_yk(n_spc), flux(n_spc)
        integer :: idx

        sec_x = [metrics%sx%m_data(i, j, k), metrics%ex%m_data(i, j, k), metrics%cx%m_data(i, j, k)]
        sec_y = [metrics%sy%m_data(i, j, k), metrics%ey%m_data(i, j, k), metrics%cy%m_data(i, j, k)]
        sec_z = [metrics%sz%m_data(i, j, k), metrics%ez%m_data(i, j, k), metrics%cz%m_data(i, j, k)]

        bottom(1:n_spc) = spcs_density%m_data(i, j, k - 1, 1:n_spc) / density%m_data(i, j, k - 1)
        top(1:n_spc) = spcs_density%m_data(i, j, k + 1, 1:n_spc) / density%m_data(i, j, k + 1)!&

        right(1:n_spc) = spcs_density%m_data(i, j - 1, k, 1:n_spc) / density%m_data(i, j - 1, k)
        left(1:n_spc) = spcs_density%m_data(i, j + 1, k, 1:n_spc) / density%m_data(i, j + 1, k)

        front(1:n_spc) = spcs_density%m_data(i - 1, j, k, 1:n_spc) / density%m_data(i - 1, j, k)
        rear(1:n_spc) = spcs_density%m_data(i + 1, j, k, 1:n_spc) / density%m_data(i + 1, j, k)

        dyk_ds(1:n_spc) = 0.5d0 * (rear - front)
        dyk_de(1:n_spc) = 0.5d0 * (left - right)
        dyk_dc(1:n_spc) = 0.5d0 * (top - bottom)

        ! diffusion of mass
        diff_density = -diffusion_coef(1:n_spc) * density%m_data(i, j, k)

        ! diffusion flux
        sec_yk = [dyk_ds, dyk_de, dyk_dc]
        !       ! scalar
        diffused_mass(1:n_spc, 1) = diff_density(1:n_spc) * (sec_x.dot.sec_yk(1:n_spc))! ruk
        diffused_mass(1:n_spc, 2) = diff_density(1:n_spc) * (sec_y.dot.sec_yk(1:n_spc))! rvk
        diffused_mass(1:n_spc, 3) = diff_density(1:n_spc) * (sec_z.dot.sec_yk(1:n_spc))! rwk


    end subroutine

    ! every data in or out for  this routine, is non-dimensional quantity
    pure subroutine transport_impl(self, spcs, spcs_density, temperature, pressure, density, tv & ! intent(in)
                                   , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient & ! intent(out)
                                   )
        use :: Constants, only:prandtl_number, schmidt_number
        class(TransportLegacy), intent(in) :: self
        class(Specie), intent(in) :: spcs(:)
        real(real64), intent(in) :: spcs_density(:)
        real(real64), intent(in) :: temperature, pressure, density, tv
        real(real64), intent(out) :: diffusion_coefficient(:)
        real(real64), intent(out) :: viscosity, thermal_conductivity, turbulent_viscosity
        real(real64) :: eu, cd, weight, t3p, xsm, total, sum_cmk, cpt, cdt, difft
        real(real64) :: td(size(spcs)), omg(size(spcs)), df(size(spcs), size(spcs))
        real(real64), dimension(size(spcs)) :: eu_k, cd_k, c, ratio, cp_k
        integer :: i, j, k, l, m, n_spc, nx, ny, nz, cnt, idx
        n_spc = size(spcs)
        ! loop for calculating laminar property
        ! this loop quite difficult to understand,
        ! to do: change loop range into realistic one
        ! calculatie mole fraction, non-dimensional
        c(1:n_spc) = spcs_density(1:n_spc) * spcs(1:n_spc)%molar_weight_nd
        c(1:n_spc) = c(1:n_spc) / sum(c(1:n_spc))

        eu_k(1:n_spc) = spcs(1:n_spc)%Eu(temperature)! calculate single component, laminar viscosity coeff
        cd_k(1:n_spc) = spcs(1:n_spc)%Cd(temperature) ! calculate single component, laminar thermal conductivity coeff

        do m = 1, n_spc
            ratio(1:n_spc) = eu_k(m) / eu_k(1:n_spc) ! get viscosity coeff ratio for this spc
            weight = spcs(m)%get_scale_factor(ratio(:), c(1:)) ! get weighting factor for this spc
            xsm = c(m) / weight
            eu = eu + xsm * eu_k(m)
            cd = cd + xsm * cd_k(m)
        end do

        viscosity = eu
        thermal_conductivity = cd

        c(1:n_spc) = spcs_density(1:n_spc) * spcs(1:n_spc)%molar_weight
        t3p = sqrt((temperature**3) / (pressure**2))
        do m = 1, n_spc
            td(1:n_spc) = temperature * spcs(m)%teab(1:)
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

        idx = check_temp_range(temperature)
        cp_k(1:n_spc) = spcs(1:n_spc)%Cp_mass(temperature, idx)
        cpt = sum(cp_k(1:n_spc) * spcs_density(1:n_spc) / density)

        turbulent_viscosity = eu + tv
        cdt = tv * cpt * (1.d0 / prandtl_number)
        difft = (tv / density) * (1.d0 / schmidt_number)

        cd = cd + cdt

        diffusion_coefficient(1:n_spc) = diffusion_coefficient(1:) + difft

    end subroutine transport_impl

end module DiffusionImpl
