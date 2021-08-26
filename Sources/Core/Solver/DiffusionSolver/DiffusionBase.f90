module DiffusionBase
    use, intrinsic :: iso_fortran_env
    use :: GridBase, only:Grid3D, PrimitiveData3D, ConservedData3D, CellMetricData3D
    use :: Vector
    use :: Debug
    use :: Metric
    implicit none

    type, abstract :: TransportProperty
    contains
        procedure(calc_transport_property), pass, deferred :: transport
    end type TransportProperty
    type, abstract :: DiffusionSolver
        type(ReferenceState) :: state_d, state_nd ! dimensional, non-dimensional reference state of flow field
        class(TransportProperty), allocatable :: m_transport
    contains
        procedure, pass :: solve_diffusion_point
        procedure, nopass :: calc_residual_RANS => calc_residual_RANS_base_impl
        procedure, nopass :: velocity_gradient => calc_velocity_gradient
        procedure, nopass :: temperature_gradient => calc_temperature_gradient
        procedure, nopass :: turbulent_gradient => calc_turbulent_gradient
        procedure(calc_mass_diffusion), pass, deferred :: diffusive_mass
        procedure(calc_viscous_stress), pass, deferred :: viscous_stress
    end type DiffusionSolver

    type :: ShearStress
        real(real64) :: tau_xx, tau_yy, tau_zz, tau_xz, tau_xy, tau_yz
    end type

    interface
        pure subroutine calc_transport_property(self, spcs, spcs_density, temperature, pressure, density, tv & ! intent(in)
                                                , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient &
                                                ) ! intent(out)
            use :: SpecieBase, only:Specie
            use :: ArrayBase
            import TransportProperty
            implicit none
            class(TransportProperty), intent(in) :: self
            class(Specie), intent(in) :: spcs(:)
            real(real64), intent(in) :: spcs_density(:)
            real(real64), intent(in) :: temperature, pressure, density, tv
            real(real64), intent(out) :: diffusion_coefficient(:)
            real(real64), intent(out) :: viscosity, thermal_conductivity, turbulent_viscosity
        end subroutine calc_transport_property

        pure subroutine calc_mass_diffusion(self, metrics, diffusion_coef, spcs_density, density, n_spc, i, j, k, diffused_mass)
            use :: ArrayBase
            import DiffusionSolver, real64
            import CellMetricData3D
            implicit none
            class(DiffusionSolver), intent(in) :: self
            class(CellMetricData3D), intent(in) :: metrics
            real(real64), intent(in) :: diffusion_coef(n_spc)
            class(Array4), intent(in) :: spcs_density
            class(Array3), intent(in) :: density
            integer, intent(in) :: n_spc, i, j, k
            real(real64), intent(out) :: diffused_mass(n_spc, 3)
        end subroutine

        pure subroutine calc_viscous_stress(self, du, dv, dw, turbulent_viscosity, i, j, k, stress)
            use :: ArrayBase
            use :: DataGrid
            import DiffusionSolver, ShearStress
            implicit none
            class(DiffusionSolver), intent(in) :: self
            type(Vector3), intent(in) :: du, dv, dw
            real(real64), intent(in) :: turbulent_viscosity
            type(ShearStress), intent(out) :: stress
            integer, intent(in) :: i, j, k
        end subroutine

        pure subroutine calc_turbulent_quantity(self, metrics)
            use :: DataGrid
            import DiffusionSolver
            class(DiffusionSolver), intent(in) :: self
            class(CellMetricData3D), intent(in) :: metrics
        end subroutine

    end interface

contains

    pure subroutine calc_temperature_gradient(metrics, t, i, j, k, dt)
        use :: ArrayBase, only:Array3
        class(CellMetricData3D), intent(in) :: metrics
        class(Array3), intent(in) :: t
        type(Vector3), intent(out) :: dt
        type(Vector3) :: dt_c
        type(Vector3) :: sec_x, sec_y, sec_z
        integer, intent(in) :: i, j, k
        sec_x = [metrics%sx%m_data(i, j, k), metrics%ex%m_data(i, j, k), metrics%cx%m_data(i, j, k)]
        sec_y = [metrics%sy%m_data(i, j, k), metrics%ey%m_data(i, j, k), metrics%cy%m_data(i, j, k)]
        sec_z = [metrics%sz%m_data(i, j, k), metrics%ez%m_data(i, j, k), metrics%cz%m_data(i, j, k)]
        associate (t => t%m_data)
            dt_c = 0.5d0*[t(i + 1, j, k) - t(i - 1, j, k), &
                          t(i, j + 1, k) - t(i, j - 1, k), &
                          t(i, j, k + 1) - t(i, j, k - 1)]
        end associate

        dt = [dt_c .dot. sec_x, & !&
              dt_c .dot. sec_y, & !&
              dt_c .dot. sec_z]

    end subroutine
    pure subroutine calc_turbulent_gradient(metrics, tk, tw, i, j, k, dtk, dtw)
        use :: ArrayBase, only:Array3
        class(CellMetricData3D), intent(in) :: metrics
        class(Array3), intent(in) :: tk, tw
        type(Vector3), intent(out) :: dtk, dtw
        type(Vector3) :: dtk_c, dtw_c
        type(Vector3) :: sec_x, sec_y, sec_z
        integer, intent(in) :: i, j, k
        sec_x = [metrics%sx%m_data(i, j, k), metrics%ex%m_data(i, j, k), metrics%cx%m_data(i, j, k)]
        sec_y = [metrics%sy%m_data(i, j, k), metrics%ey%m_data(i, j, k), metrics%cy%m_data(i, j, k)]
        sec_z = [metrics%sz%m_data(i, j, k), metrics%ez%m_data(i, j, k), metrics%cz%m_data(i, j, k)]

        associate (tk => tk%m_data, tw => tw%m_data)
            dtk_c = 0.5d0*[tk(i + 1, j, k) - tk(i - 1, j, k), &
                           tk(i, j + 1, k) - tk(i, j - 1, k), &
                           tk(i, j, k + 1) - tk(i, j, k - 1)]

            dtw_c = 0.5d0*[tw(i + 1, j, k) - tw(i - 1, j, k), &
                           tw(i, j + 1, k) - tw(i, j - 1, k), &
                           tw(i, j, k + 1) - tw(i, j, k - 1)]

        end associate

        dtk = [dtk_c .dot. sec_x, & !&
               dtk_c .dot. sec_y, & !&
               dtk_c .dot. sec_z]

        dtw = [dtw_c .dot. sec_x, & !&
               dtw_c .dot. sec_y, & !&
               dtw_c .dot. sec_z]
    end subroutine

    pure subroutine calc_velocity_gradient(metrics, u, v, w, i, j, k, du, dv, dw)
        use :: ArrayBase, only:Array3
        class(CellMetricData3D), intent(in) :: metrics
        class(Array3), intent(in) :: u, v, w
        type(Vector3), intent(out) :: du, dv, dw
        type(Vector3) :: du_c, dv_c, dw_c
        type(Vector3) :: sec_x, sec_y, sec_z
        integer, intent(in) :: i, j, k

        sec_x = [metrics%sx%m_data(i, j, k), metrics%ex%m_data(i, j, k), metrics%cx%m_data(i, j, k)]
        sec_y = [metrics%sy%m_data(i, j, k), metrics%ey%m_data(i, j, k), metrics%cy%m_data(i, j, k)]
        sec_z = [metrics%sz%m_data(i, j, k), metrics%ez%m_data(i, j, k), metrics%cz%m_data(i, j, k)]
        associate (u => u%m_data, v => v%m_data, w => w%m_data)
            !du_ds, du_de, du_dc
            du_c = 0.5d0*[u(i + 1, j, k) - u(i - 1, j, k), &
                          u(i, j + 1, k) - u(i, j - 1, k), &
                          u(i, j, k + 1) - u(i, j, k - 1)]

            dv_c = 0.5d0*[v(i + 1, j, k) - v(i - 1, j, k), &
                          v(i, j + 1, k) - v(i, j - 1, k), &
                          v(i, j, k + 1) - v(i, j, k - 1)]

            dw_c = 0.5d0*[w(i + 1, j, k) - w(i - 1, j, k), &
                          w(i, j + 1, k) - w(i, j - 1, k), &
                          w(i, j, k + 1) - w(i, j, k - 1)]
        end associate

        du = [du_c .dot. sec_x, & !&
              du_c .dot. sec_y, & !&
              du_c .dot. sec_z]

        dv = [dv_c .dot. sec_x, & !&
              dv_c .dot. sec_y, & !&
              dv_c .dot. sec_z]

        dw = [dw_c .dot. sec_x, & !&
              dw_c .dot. sec_y, & !&
              dw_c .dot. sec_z]

    end subroutine

    pure subroutine calc_residual_RANS_base_impl(metrics, residual, mass, shear_stress, viscous_work, i, j, k, n_spc)
        class(CellMetricData3D), intent(in) :: metrics
        class(ConservedData3D), intent(inout) :: residual
        type(Vector3), intent(in) :: mass(n_spc), viscous_work
        type(ShearStress), intent(in) :: shear_stress
        integer, intent(in) :: n_spc, i, j, k
        type(Vector3) :: sec_x, sec_y, sec_z ! vector for converting sec system to xyz system
        type(Vector3) :: xyz_s, xyz_e, xyz_c ! vector for converting xyz system to sec system
        type(Vector3) :: tx, ty, tz ! shear stress matrix in order of row 1, 2, 3
        real(real64) :: diffused_momentum_s, diffused_momentum_e, diffused_momentum_c
        real(real64), dimension(n_spc) :: diffused_mass_s, diffused_mass_e, diffused_mass_c
        real(real64) :: volume
        volume = 1.d0 / metrics%aj%m_data(i, j, k)
        sec_x = [metrics%sx%m_data(i, j, k), metrics%ex%m_data(i, j, k), metrics%cx%m_data(i, j, k)]
        sec_y = [metrics%sy%m_data(i, j, k), metrics%ey%m_data(i, j, k), metrics%cy%m_data(i, j, k)]
        sec_z = [metrics%sz%m_data(i, j, k), metrics%ez%m_data(i, j, k), metrics%cz%m_data(i, j, k)]

        xyz_s = [metrics%sx%m_data(i, j, k), metrics%sy%m_data(i, j, k), metrics%sz%m_data(i, j, k)]
        xyz_e = [metrics%ex%m_data(i, j, k), metrics%ey%m_data(i, j, k), metrics%ez%m_data(i, j, k)]
        xyz_c = [metrics%cx%m_data(i, j, k), metrics%cy%m_data(i, j, k), metrics%cz%m_data(i, j, k)]

        tx = [shear_stress%tau_xx, shear_stress%tau_xy, shear_stress%tau_xz]
        ty = [shear_stress%tau_xy, shear_stress%tau_yy, shear_stress%tau_yz]
        tz = [shear_stress%tau_xz, shear_stress%tau_yz, shear_stress%tau_zz]

        diffused_mass_s(1:n_spc) = (xyz_s .dot. mass(1:n_spc)) * volume!& diffused mass for s direction
        diffused_mass_e(1:n_spc) = (xyz_e .dot. mass(1:n_spc)) * volume!& diffused mass for e direction
        diffused_mass_c(1:n_spc) = (xyz_c .dot. mass(1:n_spc)) * volume!& diffused mass for c direction

        residual%rhok%m_data(i + 1, j, k, 1:n_spc) = residual%rhok%m_data(i + 1, j, k, 1:n_spc) - diffused_mass_s(1:n_spc)
        residual%rhok%m_data(i - 1, j, k, 1:n_spc) = residual%rhok%m_data(i - 1, j, k, 1:n_spc) + diffused_mass_s(1:n_spc)

        residual%rhok%m_data(i, j + 1, k, 1:n_spc) = residual%rhok%m_data(i, j + 1, k, 1:n_spc) - diffused_mass_e(1:n_spc)
        residual%rhok%m_data(i, j - 1, k, 1:n_spc) = residual%rhok%m_data(i, j - 1, k, 1:n_spc) + diffused_mass_e(1:n_spc)

        residual%rhok%m_data(i, j, k + 1, 1:n_spc) = residual%rhok%m_data(i, j, k + 1, 1:n_spc) - diffused_mass_c(1:n_spc)
        residual%rhok%m_data(i, j, k - 1, 1:n_spc) = residual%rhok%m_data(i, j, k - 1, 1:n_spc) + diffused_mass_c(1:n_spc)

        !  u momentum to computational grid system
        diffused_momentum_s =  (xyz_s .dot. tx) * volume !&
        diffused_momentum_e =  (xyz_e .dot. tx) * volume !&
        diffused_momentum_c =  (xyz_c .dot. tx) * volume !&

        residual%u_momentum%m_data(i + 1, j, k) = residual%u_momentum%m_data(i + 1, j, k) - diffused_momentum_s
        residual%u_momentum%m_data(i - 1, j, k) = residual%u_momentum%m_data(i - 1, j, k) + diffused_momentum_s

        residual%u_momentum%m_data(i, j + 1, k) = residual%u_momentum%m_data(i, j + 1, k) - diffused_momentum_e
        residual%u_momentum%m_data(i, j - 1, k) = residual%u_momentum%m_data(i, j - 1, k) + diffused_momentum_e

        residual%u_momentum%m_data(i, j, k + 1) = residual%u_momentum%m_data(i, j, k + 1) - diffused_momentum_c
        residual%u_momentum%m_data(i, j, k - 1) = residual%u_momentum%m_data(i, j, k - 1) + diffused_momentum_c
        !  v momentum to computational grid system
        diffused_momentum_s =  (xyz_s .dot. ty) * volume !&
        diffused_momentum_e =  (xyz_e .dot. ty) * volume !&
        diffused_momentum_c =  (xyz_c .dot. ty) * volume !&

        residual%v_momentum%m_data(i + 1, j, k) = residual%v_momentum%m_data(i + 1, j, k) - diffused_momentum_s
        residual%v_momentum%m_data(i - 1, j, k) = residual%v_momentum%m_data(i - 1, j, k) + diffused_momentum_s

        residual%v_momentum%m_data(i, j + 1, k) = residual%v_momentum%m_data(i, j + 1, k) - diffused_momentum_e
        residual%v_momentum%m_data(i, j - 1, k) = residual%v_momentum%m_data(i, j - 1, k) + diffused_momentum_e

        residual%v_momentum%m_data(i, j, k + 1) = residual%v_momentum%m_data(i, j, k + 1) - diffused_momentum_c
        residual%v_momentum%m_data(i, j, k - 1) = residual%v_momentum%m_data(i, j, k - 1) + diffused_momentum_c

        !  w momentum to computational grid system
        diffused_momentum_s =  (xyz_s .dot. tz) * volume !&
        diffused_momentum_e =  (xyz_e .dot. tz) * volume !&
        diffused_momentum_c =  (xyz_c .dot. tz) * volume !&

        residual%w_momentum%m_data(i + 1, j, k) = residual%w_momentum%m_data(i + 1, j, k) - diffused_momentum_s
        residual%w_momentum%m_data(i - 1, j, k) = residual%w_momentum%m_data(i - 1, j, k) + diffused_momentum_s

        residual%w_momentum%m_data(i, j + 1, k) = residual%w_momentum%m_data(i, j + 1, k) - diffused_momentum_e
        residual%w_momentum%m_data(i, j - 1, k) = residual%w_momentum%m_data(i, j - 1, k) + diffused_momentum_e

        residual%w_momentum%m_data(i, j, k + 1) = residual%w_momentum%m_data(i, j, k + 1) - diffused_momentum_c
        residual%w_momentum%m_data(i, j, k - 1) = residual%w_momentum%m_data(i, j, k - 1) + diffused_momentum_c


    end subroutine

    pure subroutine solve_diffusion_point(self, prim, residual, metrics, spcs, i, j, k)
        use :: ArrayBase
        use :: SpecieBase, only:Specie
        class(DiffusionSolver), intent(in) :: self
        class(CellMetricData3D), intent(in) :: metrics
        class(PrimitiveData3D), intent(in), allocatable :: prim
        class(ConservedData3D), intent(inout) :: residual
        class(Specie), intent(in) :: spcs(:)
        real(real64), dimension(size(spcs)) :: diffusion_coefficient, enthalpy
        real(real64), dimension(size(spcs)) :: spcs_density
        real(real64) :: viscosity, thermal_conductivity, turbulent_viscosity
        real(real64) :: diffused_energy_x, diffused_energy_y, diffused_energy_z
        real(real64) :: diffused_momentum_s, diffused_momentum_e, diffused_momentum_c
        real(real64) :: diffused_mass_xyz(size(spcs), 3)
        real(real64) :: volume, u, v, w, qx, qy, qz
        type(Vector3) :: sec_x, sec_y, sec_z ! vector for converting sec system to xyz system
        type(Vector3) :: xyz_s, xyz_e, xyz_c ! vector for converting xyz system to sec system
        type(Vector3) :: diffused_mass(size(spcs)), viscous_work
        type(Vector3) :: du, dv, dw, dtk, dtw, dt ! gradient value in physical grid system, ex) du == [dudx, dudy, dudz]
        type(ShearStress) :: stress
        integer, intent(in) :: i, j, k
        integer :: idx, n_spc

        n_spc = size(spcs)
        idx = merge(1, 2, prim%t%m_data(i, j, k) >= self%state_nd%ref_temperature)
        enthalpy(1:n_spc) = spcs(:)%H_mass(prim%t%m_data(i, j, k), idx)

        call self%m_transport%transport(spcs, prim%rhok%m_data(i, j, k, 1:n_spc), prim%t%m_data(i, j, k) &
                                        , prim%p%m_data(i, j, k), prim%rho%m_data(i, j, k), prim%tv%m_data(i, j, k) &
                                        , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient)

        ! calculate diffusion of mass for physical direction (x, y, z)
        call self%diffusive_mass(metrics, diffusion_coefficient, prim%rhok, prim%rho, n_spc, i, j, k, diffused_mass_xyz)

        diffused_mass(1:n_spc) = [diffused_mass_xyz(:, 1), diffused_mass_xyz(:, 2), diffused_mass_xyz(:, 3)]

        call self%velocity_gradient(metrics, prim%u, prim%v, prim%w, i, j, k, du, dv, dw)
        call self%temperature_gradient(metrics, prim%t, i, j, k, dt)
        call self%turbulent_gradient(metrics, prim%tk, prim%tw, i, j, k, dtk, dtw)
        call self%viscous_stress(du, dv, dw, turbulent_viscosity, i, j, k, stress)

        ! todo : hide calculation of q_dot

        u = prim%u%m_data(i, j, k)
        v = prim%v%m_data(i, j, k)
        w = prim%w%m_data(i, j, k)
        diffused_energy_x = sum(enthalpy(1:n_spc) * diffused_mass(1:n_spc)%x)  ! diffuesed energy of each species in x dir
        diffused_energy_y = sum(enthalpy(1:n_spc) * diffused_mass(1:n_spc)%y)  ! diffuesed energy of each species in y dir
        diffused_energy_z = sum(enthalpy(1:n_spc) * diffused_mass(1:n_spc)%z)  ! diffuesed energy of each species in z dir

        qx = -thermal_conductivity * dt%x + diffused_energy_x
        qy = -thermal_conductivity * dt%y + diffused_energy_y
        qz = -thermal_conductivity * dt%z + diffused_energy_z

        viscous_work%x = (u * stress%tau_xx + v * stress%tau_xy + w * stress%tau_xz) - qx
        viscous_work%y = (u * stress%tau_xy + v * stress%tau_yy + w * stress%tau_yz) - qy
        viscous_work%z = (u * stress%tau_xz + v * stress%tau_yz + w * stress%tau_zz) - qz
        ! end todo

        call self%calc_residual_RANS(metrics, residual, diffused_mass, stress, viscous_work, i, j, k, n_spc)

        ! calculating residuals..! for now, we are using global residuals..

    end subroutine

end module DiffusionBase
