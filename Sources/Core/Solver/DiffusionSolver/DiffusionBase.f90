module DiffusionBase
    use, intrinsic :: iso_fortran_env
    use :: GridBase, only:Grid3D, PrimitiveData3D, ConservedData3D, CellMetricData3D
    use :: ArrayBase
    use :: Vector
    use :: Metric
    implicit none

    type, abstract :: TransportProperty
    contains
        procedure(calc_transport_property), pass, deferred :: transport
    end type TransportProperty

    type :: SST_DES
        real(real64) :: turbulent_viscosity
        real(real64) :: tk_source, tw_source
    end type

    type :: TurbulentProperty
        type(SST_DES) :: SST
    end type

    type, abstract :: DiffusionSolver
        type(ReferenceState) :: state_d, state_nd ! dimensional, non-dimensional reference state of flow field
        class(TransportProperty), allocatable :: m_transport
    contains
        procedure(solve_diffusion_point_interface), deferred, pass :: solve_diffusion_point
        procedure(calc_mass_diffusion_interface), pass, deferred :: diffusive_mass
        procedure(calc_viscous_stress_interface), pass, deferred :: viscous_stress
        procedure(calc_turbulent_property_interface), pass, deferred :: turbulent_property

        procedure :: calc_residual_RANS => calc_residual_RANS_base_impl
        procedure, nopass :: momentum_residual => calc_residual_momentum_impl
        procedure, nopass :: velocity_gradient => calc_velocity_gradient
        procedure, nopass :: temperature_gradient => calc_temperature_gradient
        procedure, nopass :: turbulent_gradient => calc_turbulent_gradient
    end type DiffusionSolver

    type :: ShearStress
        real(real64) :: tau_xx, tau_yy, tau_zz, tau_xz, tau_xy, tau_yz
    end type

    interface
        pure subroutine calc_turbulent_property_interface(self &
                                                          , y_metric, viscosity, density &
                                                          , tk, tw, wall_distance, volume &
                                                          , dtk, dtw, du, dv, dw, property)
            import DiffusionSolver, CellMetricData3D, Vector3, Array3, TurbulentProperty
            import real64
            class(DiffusionSolver), intent(in) :: self
            real(real64), intent(in) :: y_metric, viscosity
            real(real64), intent(in) :: density, tk, tw, wall_distance, volume
            type(Vector3), intent(in) :: dtk, dtw, du, dv, dw
            type(TurbulentProperty), intent(out) :: property
        end subroutine
        pure subroutine calc_transport_property(self, spcs, spcs_density, temperature, pressure, density, eddy_viscosity & ! intent(in)
                                                , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient &
                                                ) ! intent(out)
            use :: SpecieBase, only:Specie
            use :: ArrayBase
            import TransportProperty
            implicit none
            class(TransportProperty), intent(in) :: self
            class(Specie), intent(in) :: spcs(:)
            real(real64), intent(in) :: spcs_density(:)
            real(real64), intent(in) :: temperature, pressure, density, eddy_viscosity
            real(real64), intent(out) :: diffusion_coefficient(:)
            real(real64), intent(out) :: viscosity, thermal_conductivity, turbulent_viscosity
        end subroutine calc_transport_property

        pure subroutine calc_mass_diffusion_interface(self, metrics, diffusion_coef &
                                                      , spcs_density, density, n_spc, i, j, k, diffused_mass)
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

        pure subroutine calc_viscous_stress_interface(self, du, dv, dw, turbulent_viscosity, stress)
            use :: ArrayBase
            use :: DataGrid
            import DiffusionSolver, ShearStress
            implicit none
            class(DiffusionSolver), intent(in) :: self
            type(Vector3), intent(in) :: du, dv, dw
            real(real64), intent(in) :: turbulent_viscosity
            type(ShearStress), intent(out) :: stress
        end subroutine
        pure subroutine solve_diffusion_point_interface(self, prim, residual, metrics, spcs, i, j, k)
            use :: ArrayBase
            use :: GridBase
            use :: SpecieBase, only:Specie
            import DiffusionSolver
            class(DiffusionSolver), intent(in) :: self
            class(CellMetricData3D), intent(in) :: metrics
            class(PrimitiveData3D), intent(in), allocatable :: prim
            class(ConservedData3D), intent(inout) :: residual
            class(Specie), intent(in) :: spcs(:)
            integer, intent(in) :: i, j, k
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

    pure subroutine calc_residual_momentum_impl(target_component, volume, xyz_s, xyz_e, xyz_c, shear_comp, i, j, k)
        class(Array3), intent(inout) :: target_component
        type(Vector3), intent(in) :: xyz_s, xyz_e, xyz_c
        type(Vector3), intent(in) :: shear_comp
        real(real64), intent(in) :: volume
        integer, intent(in) :: i, j, k
        real(real64) :: diffused_momentum_s, diffused_momentum_e, diffused_momentum_c

        !  momentum to computational grid system
        diffused_momentum_s =  (xyz_s .dot. shear_comp) * volume !&
        diffused_momentum_e =  (xyz_e .dot. shear_comp) * volume !&
        diffused_momentum_c =  (xyz_c .dot. shear_comp) * volume !&

        target_component%m_data(i + 1, j, k) = target_component%m_data(i + 1, j, k) - diffused_momentum_s
        target_component%m_data(i - 1, j, k) = target_component%m_data(i - 1, j, k) + diffused_momentum_s

        target_component%m_data(i, j + 1, k) = target_component%m_data(i, j + 1, k) - diffused_momentum_e
        target_component%m_data(i, j - 1, k) = target_component%m_data(i, j - 1, k) + diffused_momentum_e

        target_component%m_data(i, j, k + 1) = target_component%m_data(i, j, k + 1) - diffused_momentum_c
        target_component%m_data(i, j, k - 1) = target_component%m_data(i, j, k - 1) + diffused_momentum_c
    end subroutine

    pure subroutine calc_residual_RANS_base_impl(self, metrics, residual, mass, shear_stress, viscous_work, i, j, k, n_spc)
        use :: Constants, only:Re_inf_mach1
        class(DiffusionSolver), intent(in) :: self
        class(CellMetricData3D), intent(in) :: metrics
        class(ConservedData3D), intent(inout) :: residual
        type(Vector3), intent(in) :: mass(n_spc), viscous_work
        type(ShearStress), intent(in) :: shear_stress
        integer, intent(in) :: n_spc, i, j, k
        type(Vector3) :: sec_x, sec_y, sec_z ! vector for converting sec system to xyz system
        type(Vector3) :: xyz_s, xyz_e, xyz_c ! vector for converting xyz system to sec system
        type(Vector3) :: tx, ty, tz ! shear stress matrix in order of row 1, 2, 3
        real(real64), dimension(n_spc) :: diffused_mass_s, diffused_mass_e, diffused_mass_c
        real(real64) :: volume, work_s, work_e, work_c
        volume = 1.d0 / (metrics%aj%m_data(i, j, k) * Re_inf_mach1)

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

        ! diffused mass residual
        residual%rhok%m_data(i + 1, j, k, 1:n_spc) = residual%rhok%m_data(i + 1, j, k, 1:n_spc) - diffused_mass_s(1:n_spc)
        residual%rhok%m_data(i - 1, j, k, 1:n_spc) = residual%rhok%m_data(i - 1, j, k, 1:n_spc) + diffused_mass_s(1:n_spc)

        residual%rhok%m_data(i, j + 1, k, 1:n_spc) = residual%rhok%m_data(i, j + 1, k, 1:n_spc) - diffused_mass_e(1:n_spc)
        residual%rhok%m_data(i, j - 1, k, 1:n_spc) = residual%rhok%m_data(i, j - 1, k, 1:n_spc) + diffused_mass_e(1:n_spc)

        residual%rhok%m_data(i, j, k + 1, 1:n_spc) = residual%rhok%m_data(i, j, k + 1, 1:n_spc) - diffused_mass_c(1:n_spc)
        residual%rhok%m_data(i, j, k - 1, 1:n_spc) = residual%rhok%m_data(i, j, k - 1, 1:n_spc) + diffused_mass_c(1:n_spc)

        call self%momentum_residual(residual%u_momentum, volume, xyz_s, xyz_e, xyz_c, tx, i, j, k)
        call self%momentum_residual(residual%v_momentum, volume, xyz_s, xyz_e, xyz_c, ty, i, j, k)
        call self%momentum_residual(residual%w_momentum, volume, xyz_s, xyz_e, xyz_c, tz, i, j, k)

        work_s = (viscous_work .dot. xyz_s) * volume!&
        work_e = (viscous_work .dot. xyz_e) * volume!&
        work_c = (viscous_work .dot. xyz_c) * volume!&

        residual%e%m_data(i + 1, j, k) = residual%e%m_data(i + 1, j, k) - work_s
        residual%e%m_data(i - 1, j, k) = residual%e%m_data(i - 1, j, k) + work_s

        residual%e%m_data(i, j + 1, k) = residual%e%m_data(i, j + 1, k) - work_e
        residual%e%m_data(i, j - 1, k) = residual%e%m_data(i, j - 1, k) + work_e

        residual%e%m_data(i, j, k + 1) = residual%e%m_data(i, j, k + 1) - work_c
        residual%e%m_data(i, j, k - 1) = residual%e%m_data(i, j, k - 1) + work_c

    end subroutine

end module DiffusionBase
