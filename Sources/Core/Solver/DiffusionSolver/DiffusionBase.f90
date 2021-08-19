module DiffusionBase
    use, intrinsic :: iso_fortran_env
    implicit none

    type, abstract :: TransportProperty
    contains
        procedure(calc_transport_property), pass, deferred :: transport
    end type TransportProperty
    type, abstract :: DiffusionSolver(n_spc)
        integer, len :: n_spc
        class(TransportProperty), allocatable :: m_transport
    contains
        procedure, pass :: solve_diffusion
        procedure(calc_diffusive_flux), pass, deferred :: diffusive_flux
    end type DiffusionSolver


    interface
        subroutine calc_diffusive_flux(self)
            import DiffusionSolver
            class(DiffusionSolver(*)), intent(in) :: self
        end subroutine
         subroutine calc_transport_property(self, spcs, spcs_density, temperature, pressure, density & ! intent(in)
                                                , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient&
                                                , resolution) ! intent(out)
            use :: SpecieBase, only:Specie
            use :: ArrayBase
            import TransportProperty
            class(TransportProperty), intent(in) :: self
            class(Specie), intent(in) :: spcs(:)
            class(Array4), intent(in) :: spcs_density
            class(Array3), intent(in) :: temperature, pressure, density
            type(Array4), intent(inout) :: diffusion_coefficient
            type(Array3), intent(inout) :: viscosity, thermal_conductivity, turbulent_viscosity
            integer, intent(in) :: resolution(3)
        end subroutine calc_transport_property
    end interface

contains

    subroutine solve_diffusion(self, prim, conserv, metrics, spcs)
        use :: ArrayBase
        use :: GridBase, only:Grid3D, PrimitiveData3D, ConservedData3D, CellMetricData3D
        use :: SpecieBase, only:Specie
        class(DiffusionSolver(*)), intent(in) :: self
        class(CellMetricData3D), intent(in) :: metrics
        class(PrimitiveData3D), intent(in) :: prim
        class(ConservedData3D), intent(inout) :: conserv
        class(Specie), intent(in) :: spcs(self%n_spc)
        type(Array4) :: diffusion_coefficient
        type(Array3) :: viscosity, thermal_conductivity, turbulent_viscosity

        viscosity = Array3(prim%t)
        thermal_conductivity = Array3(prim%t)
        turbulent_viscosity = Array3(prim%t)
        diffusion_coefficient = Array4(prim%rhok)

        call self%m_transport%transport(spcs, prim%rhok, prim%t, prim%p, prim%rho &
                                        , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient, prim%m_resolution)
        !call self%diffusive_flux
    end subroutine

end module DiffusionBase
