module DiffusionBase
    use, intrinsic :: iso_fortran_env
    implicit none

    type, abstract :: DiffusionSolver
        class(TransportProperty), pointer :: m_transport
    contains
    end type DiffusionSolver

    type, abstract :: TransportProperty
        real(real64), allocatable, dimension(:, :) :: fa, fb ! todo=>change this damn things name
    contains
        procedure(calc_transport_property), pass, deferred :: transport
        !procedure(calc_viscosity), private, deferred :: viscosity
    end type TransportProperty

    interface 
    pure subroutine calc_transport_property(self, spcs, spcs_density, mole_fraction, temperature, pressure, density & ! intent(in)
                              , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_quantity) ! intent(out)
        use :: Species, only:SpecieBase
        use :: ArrayBase
        import TransportProperty
        class(TransportProperty), intent(in) :: self
        class(SpecieBase), intent(in) :: spcs(:)
        class(Array4), intent(in) :: spcs_density, mole_fraction
        class(Array3), intent(in) :: temperature, pressure, density
        class(Array4), intent(out) :: diffusion_quantity
        class(Array3), intent(out) :: viscosity, thermal_conductivity, turbulent_viscosity
    end subroutine calc_transport_property
    end interface
    !interface
    !    pure subroutine calc_viscosity(spcs)
    !        import TransportProperty
    !        use :: Species, only:SpecieBase
    !        implicit none(type, external)
    !        class(SpecieBase), intent(in) :: spcs(:)
    !    end subroutine calc_viscosity
    !end interface
contains

end module DiffusionBase
