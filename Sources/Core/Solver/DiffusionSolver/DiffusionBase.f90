module DiffusionBase
    use, intrinsic :: iso_fortran_env
    use :: GridBase, only:Grid3D, PrimitiveData3D, ConservedData3D, CellMetricData3D
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
        procedure(calc_mass_diffusion), pass, deferred :: diffusive_mass
    end type DiffusionSolver

    interface
        subroutine calc_diffusive_flux(self)
            import DiffusionSolver
            class(DiffusionSolver), intent(in) :: self
        end subroutine
        pure subroutine calc_transport_property(self, spcs, spcs_density, temperature, pressure, density, tv & ! intent(in)
                                                , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient &
                                                ) ! intent(out)
            use :: SpecieBase, only:Specie
            use :: ArrayBase
            import TransportProperty
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
            class(DiffusionSolver), intent(in) :: self
            class(CellMetricData3D), intent(in) :: metrics
            real(real64), intent(in) ::  diffusion_coef(n_spc)
            class(Array4), intent(in) :: spcs_density
            class(Array3), intent(in):: density
            integer, intent(in) :: n_spc, i, j, k
            real(real64), intent(out) :: diffused_mass(n_spc, 3)
        end subroutine

    end interface

contains

    pure subroutine solve_diffusion_point(self, prim, conserv, metrics, spcs, i, j, k)
        use :: ArrayBase
        use :: SpecieBase, only:Specie
        class(DiffusionSolver), intent(in) :: self
        class(CellMetricData3D), intent(in) :: metrics
        class(PrimitiveData3D), intent(in), allocatable :: prim
        class(ConservedData3D), intent(inout) :: conserv
        class(Specie), intent(in) :: spcs(:)
        real(real64), dimension(size(spcs)) :: diffusion_coefficient, enthalpy
        real(real64), dimension(size(spcs)) :: spcs_density
        real(real64) :: viscosity, thermal_conductivity, turbulent_viscosity
        real(real64) :: diffused_mass(size(spcs), 3)
        integer, intent(in) ::  i, j, k
        integer :: idx, n_spc

        n_spc = size(spcs)
        idx = merge(1, 2, prim%t%m_data(i, j, k) >= self%state_nd%ref_temperature)
        ! value for i, j, k current cell
        !spcs_density = prim%rhok%m_data(i, j, k, 1:n_spc) 
        !total_density = prim%rho%m_data(i, j, k)
        enthalpy = spcs(:)%H_mass(prim%t%m_data(i, j, k), idx)
        call self%m_transport%transport(spcs, prim%rhok%m_data(i, j, k, 1:n_spc), prim%t%m_data(i, j, k) &
                                        , prim%p%m_data(i, j, k), prim%rho%m_data(i, j, k), prim%tv%m_data(i, j, k) &
                                        , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient)

        ! result of this subroutine does not multiplied by volume of each cell
        call self%diffusive_mass(metrics, diffusion_coefficient, prim%rhok, prim%rho, n_spc, i, j, k, diffused_mass)
    end subroutine

end module DiffusionBase
