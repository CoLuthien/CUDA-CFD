module DiffusionImpl
    use, intrinsic :: iso_fortran_env
    use :: DiffusionBase
    use :: ArrayBase
    use :: SpecieBase
    use :: Vector
    implicit none

    type, extends(DiffusionSolver) :: DiffusionLegacy
    contains
        !procedure :: solve_diffusion => diffusion
        procedure :: solve_diffusion_point => solve_diffusion_point_legacy
        procedure :: turbulent_property => calc_turbulent_property_impl
        procedure, non_overridable :: diffusive_mass => diffusive_mass_impl
        procedure:: viscous_stress => viscous_stress_impl
    end type DiffusionLegacy

    type, extends(TransportProperty) :: TransportLegacy
    contains
        procedure :: transport => transport_impl
    end type TransportLegacy

    interface DiffusionLegacy
        module procedure :: init_diffusion_legacy
    end interface

contains

    ! please refer https://doi.org/10.2514/3.12149
!    pure subroutine calc_turbulent_property_impl(self &
!                                                 , y_metric, viscosity, density &
!                                                 , tk, tw, wall_distance, volume &
!                                                 , dtk, dtw, du, dv, dw, property)
!        use :: Constants
!        class(DiffusionLegacy), intent(in) :: self
!        real(real64), intent(in) :: y_metric, viscosity, density, tk, tw, wall_distance, volume
!        type(Vector3), intent(in) :: dtk, dtw, du, dv, dw
!        type(TurbulentProperty), intent(out) :: property
!
    pure subroutine calc_turbulent_property_impl(self &
                                                 , y_metric, viscosity, density &
                                                 , tk, tw, wall_distance, volume &
                                                 , dtk, dtw, du, dv, dw, property)
        use :: Constants
        class(DiffusionLegacy), intent(in) :: self
        real(real64), intent(in) :: y_metric, viscosity
        real(real64), intent(in) :: density, tk, tw, wall_distance, volume
        type(Vector3), intent(in) :: dtk, dtw, du, dv, dw
        type(TurbulentProperty), intent(out) :: property

        real(real64) :: y_metric_sqr, f1_arg1, f1_arg2, f1_arg3
        real(real64) :: tv_arg1, tv_arg2, tv_arg3, tv_kinematic, tv_dynamic
        real(real64) :: cd_t, cd_kw, arg1, t_omega, arg2, correction
        real(real64) :: f1, f2, tk_production
        real(real64) :: tau_xx, tau_yy, tau_zz, tau_xz, tau_xy, tau_yz
        real(real64) :: sst_gamma, sst_sigma_k, sst_sigma_w, sst_beta
        type(ShearStress) :: reynolds_stress

        y_metric_sqr = y_metric**2
        cd_t = ((2.d0 * kw_sigma_w2 * density) * (dtk.dot.dtw)) / tw
        cd_kw = max(cd_t, 1d-10)

        ! Blending function F1
        f1_arg1 = sqrt(tk) / (kw_beta_star * tw * y_metric)
        f1_arg2 = (500 * viscosity) / (Re_inf_mach1 * density * tw * y_metric_sqr)
        f1_arg3 = (4.d0 * density * kw_sigma_w2 * tk) / (cd_kw * y_metric_sqr)

        arg1 = min(max(f1_arg1, f1_arg2), f1_arg3)
        f1 = tanh(arg1**4)

        !Turbulent Viscosity
        t_omega = size(curl(du, dv, dw)) !sqrt((dw%y - dv%z)**2 + (du%z - dw%x)**2 + (dv%x - du%y)**2) ! curl of velocity
        arg2 = max(2.d0 * f1_arg1, f1_arg2)
        f2 = tanh(arg2**2)

        tv_arg1 = kw_a1 * tw
        tv_arg2 = t_omega * f2
        tv_arg3 = tv_arg1 + tv_arg2 - 2.d0 * tv_arg1 * tv_arg2 / (tv_arg1 + tv_arg2)
        tv_kinematic = kw_a1 * tk / tv_arg3
        tv_dynamic = tv_kinematic * density
        property%SST%turbulent_viscosity = tv_dynamic * Re_inf_mach1

        ! Calculate SST constants
        sst_gamma = f1 * kw_gamma1 + (1 - f1) * kw_gamma2
        sst_sigma_k = f1 * kw_sigma_k1 + (1 - f1) * kw_sigma_k2
        sst_sigma_w = f1 * kw_sigma_w1 + (1 - f1) * kw_sigma_w2
        sst_beta = f1 * kw_beta1 + (1 - f1) * kw_beta2

        !calculate basic component of shear stress
        call self%viscous_stress(du, dv, dw, property%SST%turbulent_viscosity, reynolds_stress)
        ! Correct with (2/3) * rho * tk
        correction = (2.d0 / 3.d0) * density * tk
        tau_xx = reynolds_stress%tau_xx - correction
        tau_yy = reynolds_stress%tau_yy - correction
        tau_zz = reynolds_stress%tau_zz - correction

        tau_xy = reynolds_stress%tau_xy
        tau_yz = reynolds_stress%tau_yz
        tau_xz = reynolds_stress%tau_xz

        !turbulent kinetic energy production

        tk_production = tau_xx * du%x + tau_yy * dv%y + tau_zz * dw%z &
                        + tau_xy * (du%y + dv%x) + tau_yz * (dv%z + dw%y) + tau_xz * (dw%x + du%z)

        block ! code copy cat region
            real(real64) :: alkw, c_des, switch, f_sst, f_des
            alkw = sqrt(tk) / (sst_beta * tw)
            c_des = (1.d0 - f1) * 0.61d0 + f1 * 0.78d0
            switch = alkw / (c_des * wall_distance)
            f_sst = f2  ! Fsst = 0.d0, ff1 or ff2 ?
            f_des = dmax1(switch * (1.d0 - f_sst), 1.d0)

            property%SST%tk_source = (tk_production - kw_beta_star * tw * density * tk * f_des) * volume
            property%SST%tk_source = volume * &
                                     (sst_gamma / tv_kinematic * tk_production - sst_beta * density * tw * tw + (1.d0 - f1) * cd_t)

        end block
    end subroutine

    ! todo: refine init diffusion
    function init_diffusion_legacy(n_spc) result(solver)
        integer :: n_spc
        type(DiffusionLegacy) :: solver
        type(TransportLegacy) :: prop

        solver%m_transport = prop
    end function

    pure subroutine viscous_stress_impl(self, du, dv, dw, turbulent_viscosity, stress)
        use :: ArrayBase
        use :: Vector
        class(DiffusionLegacy), intent(in) :: self
        type(Vector3), intent(in) :: du, dv, dw ! gradient representation of physical grid system
        type(ShearStress), intent(out) :: stress
        real(real64), intent(in) :: turbulent_viscosity
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
    pure subroutine transport_impl(self, spcs, spcs_density, temperature, pressure, density, eddy_viscosity & ! intent(in)
                                   , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient & ! intent(out)
                                   )
        use :: Constants, only:prandtl_number, schmidt_number
        class(TransportLegacy), intent(in) :: self
        class(Specie), intent(in) :: spcs(:)
        real(real64), intent(in) :: spcs_density(:)
        real(real64), intent(in) :: temperature, pressure, density, eddy_viscosity
        real(real64), intent(out) :: diffusion_coefficient(:)
        real(real64), intent(out) :: viscosity, thermal_conductivity, turbulent_viscosity
        real(real64) :: eu, cd, weight, t3p, xsm, total, sum_cmk, cpt, cdt, difft
        real(real64) :: td(size(spcs)), omg(size(spcs)), df(size(spcs), size(spcs))
        real(real64), dimension(size(spcs)) :: eu_k, cd_k, c, ratio, cp_k
        integer :: i, j, k, l, m, n_spc, nx, ny, nz, cnt, idx
        n_spc = size(spcs)
        ! loop for calculating laminar property
        ! this loop quite difficult to understand,
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

        turbulent_viscosity = eu + eddy_viscosity
        cdt = eddy_viscosity * cpt * (1.d0 / prandtl_number)
        difft = (eddy_viscosity / density) * (1.d0 / schmidt_number)

        cd = cd + cdt

        diffusion_coefficient(1:n_spc) = diffusion_coefficient(1:) + difft

    end subroutine transport_impl

    pure subroutine solve_diffusion_point_legacy(self, prim, residual, metrics, spcs, i, j, k)
        use :: ArrayBase
        use :: SpecieBase, only:Specie
        class(DiffusionLegacy), intent(in) :: self
        class(CellMetricData3D), intent(in) :: metrics
        class(PrimitiveData3D), intent(in), allocatable :: prim
        class(ConservedData3D), intent(inout) :: residual
        class(Specie), intent(in) :: spcs(:)
        real(real64), dimension(size(spcs)) :: diffusion_coefficient, enthalpy
        real(real64), dimension(size(spcs)) :: spcs_density
        real(real64) :: viscosity, thermal_conductivity, turbulent_viscosity
        real(real64) :: diffused_energy_x, diffused_energy_y, diffused_energy_z
        real(real64) :: diffused_mass_xyz(size(spcs), 3)
        real(real64) :: volume, u, v, w, qx, qy, qz
        real(real64) :: rho_tk, rho_tw, density, tk, tw, dxx, dyy, dzz, wall_distance
        type(Vector3) :: diffused_mass(size(spcs)), viscous_work
        type(Vector3) :: du, dv, dw, dtk, dtw, dt ! gradient value in physical grid system, ex) du == [dudx, dudy, dudz]
        type(ShearStress) :: stress
        type(TurbulentProperty) :: turb_prop
        integer, intent(in) :: i, j, k
        integer :: idx, n_spc

        n_spc = size(spcs)
        idx = merge(1, 2, prim%t%m_data(i, j, k) >= self%state_nd%ref_temperature)
        u = prim%u%m_data(i, j, k); density = prim%rho%m_data(i, j, k)
        v = prim%v%m_data(i, j, k); tk = prim%tk%m_data(i, j, k)
        w = prim%w%m_data(i, j, k); tw = prim%tw%m_data(i, j, k)


        enthalpy(1:n_spc) = spcs(:)%H_mass(prim%t%m_data(i, j, k), idx)

        call self%m_transport%transport(spcs, prim%rhok%m_data(i, j, k, 1:n_spc), prim%t%m_data(i, j, k) &
                                        , prim%p%m_data(i, j, k), prim%rho%m_data(i, j, k), prim%tv%m_data(i, j, k) &
                                        , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient)

        ! calculate diffusion of mass for physical direction (x, y, z)
        call self%velocity_gradient(metrics, prim%u, prim%v, prim%w, i, j, k, du, dv, dw)
        call self%temperature_gradient(metrics, prim%t, i, j, k, dt)
        call self%turbulent_gradient(metrics, prim%tk, prim%tw, i, j, k, dtk, dtw)

        call self%diffusive_mass(metrics, diffusion_coefficient, prim%rhok, prim%rho, n_spc, i, j, k, diffused_mass_xyz)
        call self%viscous_stress(du, dv, dw, turbulent_viscosity, stress)

        diffused_mass(1:n_spc) = [diffused_mass_xyz(:, 1), diffused_mass_xyz(:, 2), diffused_mass_xyz(:, 3)]

        ! todo : hide calculation of viscous work
        diffused_energy_x = sum(enthalpy(1:n_spc) * diffused_mass(1:n_spc)%x)  ! diffuesed energy of each species in x dir
        diffused_energy_y = sum(enthalpy(1:n_spc) * diffused_mass(1:n_spc)%y)  ! diffuesed energy of each species in y dir
        diffused_energy_z = sum(enthalpy(1:n_spc) * diffused_mass(1:n_spc)%z)  ! diffuesed energy of each species in z dir

        qx = -thermal_conductivity * dt%x + diffused_energy_x
        qy = -thermal_conductivity * dt%y + diffused_energy_y
        qz = -thermal_conductivity * dt%z + diffused_energy_z

        !q_dot
        viscous_work%x = (u * stress%tau_xx + v * stress%tau_xy + w * stress%tau_xz) - qx
        viscous_work%y = (u * stress%tau_xy + v * stress%tau_yy + w * stress%tau_yz) - qy
        viscous_work%z = (u * stress%tau_xz + v * stress%tau_yz + w * stress%tau_zz) - qz
        ! end todo

        ! calculating non-turbulent residuals..! for now, we are using global residuals..
        call self%calc_residual_RANS(metrics, residual, diffused_mass, stress, viscous_work, i, j, k, n_spc)
        volume = 1.d0 / metrics%aj%m_data(i, j, k)
        dxx = volume * (metrics%ey%m_data(i, j, k) * metrics%cz%m_data(i, j, k) &
                        - metrics%ez%m_data(i, j, k) * metrics%cy%m_data(i, j, k))

        dyy = volume * (metrics%sx%m_data(i, j, k) * metrics%cz%m_data(i, j, k) &
                        - metrics%sz%m_data(i, j, k) * metrics%cx%m_data(i, j, k))

        dzz = volume * (metrics%sx%m_data(i, j, k) * metrics%ey%m_data(i, j, k) &
                        - metrics%sy%m_data(i, j, k) * metrics%ex%m_data(i, j, k))

        wall_distance = max(dxx, dyy, dzz)

        call self%turbulent_property(metrics%dely%m_data(i, j, k), viscosity, density, tk, tw, wall_distance, volume &
                                     , dtk, dtw, du, dv, dw, turb_prop)

    end subroutine

end module DiffusionImpl
