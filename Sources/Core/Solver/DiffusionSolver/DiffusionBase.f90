module DiffusionBase
    use, intrinsic :: iso_fortran_env
    implicit none

    type, abstract :: DiffusionSolver
        class(TransportProperty), pointer :: m_transport
    contains
        procedure, pass :: solve_diffusion
    end type DiffusionSolver

    type, abstract :: TransportProperty
    contains
        procedure(calc_transport_property), pass, deferred :: transport
        !procedure(calc_viscosity), private, deferred :: viscosity
    end type TransportProperty

    interface 
    pure subroutine calc_transport_property(self, spcs, spcs_density, mole_fraction, temperature, pressure, density & ! intent(in)
                              , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_quantity) ! intent(out)
        use :: SpecieBase, only:Specie
        use :: ArrayBase
        import TransportProperty
        class(TransportProperty), intent(in) :: self
        class(Specie), intent(in) :: spcs(:)
        class(Array4), intent(in) :: spcs_density, mole_fraction
        class(Array3), intent(in) :: temperature, pressure, density
        class(Array4), intent(out) :: diffusion_quantity
        class(Array3), intent(out) :: viscosity, thermal_conductivity, turbulent_viscosity
    end subroutine calc_transport_property
    end interface
contains
    subroutine solve_diffusion(self)
        implicit none
        class(DiffusionSolver),intent(in) :: self

    end subroutine 

end module DiffusionBase
