module DiffusionBase
    use, intrinsic :: iso_fortran_env
    implicit none

    type, abstract :: DiffusionSolver
        class(TransportProperty), pointer :: m_transport
    contains
        procedure, pass :: solve_diffusion
        procedure(calc_diffusive_flux), pass, deferred :: diffusive_flux
    end type DiffusionSolver

    type, abstract :: TransportProperty
    contains
        procedure(calc_transport_property), pass, deferred :: transport
        !procedure(calc_viscosity), private, deferred :: viscosity
    end type TransportProperty

    interface
        subroutine calc_diffusive_flux(self)
            import DiffusionSolver
            class(DiffusionSolver), intent(in) :: self
        end subroutine
        pure subroutine calc_transport_property(self, spcs, spcs_density, temperature, pressure, density & ! intent(in)
                                                , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient) ! intent(out)
            use :: SpecieBase, only:Specie
            use :: ArrayBase
            import TransportProperty
            class(TransportProperty), intent(in) :: self
            class(Specie), intent(in) :: spcs(:)
            class(Array4), intent(in) :: spcs_density
            class(Array3), intent(in) :: temperature, pressure, density
            type(Array4), intent(out) :: diffusion_coefficient
            type(Array3), intent(out) :: viscosity, thermal_conductivity, turbulent_viscosity
        end subroutine calc_transport_property
    end interface
contains

    subroutine solve_diffusion(self, prim, conserv, metrics, spcs)
        use :: ArrayBase
        use :: GridBase, only:Grid3D, PrimitiveData3D, ConservedData3D, MetricData3D
        use :: SpecieBase, only:Specie
        class(DiffusionSolver), intent(in) :: self
        class(MetricData3D), intent(in) :: metrics
        class(PrimitiveData3D), intent(in) :: prim
        class(ConservedData3D), intent(inout) :: conserv
        class(Specie), intent(in) :: spcs(:)
        type(Array4) :: diffusion_coefficient
        type(Array3) :: viscosity, thermal_conductivity, turbulent_viscosity

        viscosity = Array(prim%t)
        thermal_conductivity = Array(prim%t)
        turbulent_viscosity = Array(prim%t)
        diffusion_coefficient = Array(prim%rhok)

        call self%m_transport%transport(spcs, prim%rhok, prim%t, prim%p, prim%rho &
                                        , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient)
        !call self%diffusive_flux
    end subroutine

end module DiffusionBase
