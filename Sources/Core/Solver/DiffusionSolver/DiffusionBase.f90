module DiffusionBase
    use, intrinsic :: iso_fortran_env
    implicit none

    type, abstract :: DiffusionSolver
        class(TransportProperty), pointer :: m_transport
    contains
    end type DiffusionSolver

    type, abstract :: TransportProperty
    contains
        procedure(calc_transport_property), pass, deferred :: transport
        !procedure(calc_viscosity), private, deferred :: viscosity
    end type TransportProperty

    interface 
    pure subroutine calc_transport_property(self &
                              , spcs, spcs_density, temperature, pressure, density &
                              , viscosity, thermal_conductivity)
        use :: Species, only:SpecieBase
        use :: ArrayBase
        import TransportProperty
        class(TransportProperty), intent(in) :: self
        class(SpecieBase), intent(in) :: spcs(:)
        type(Array4), intent(in) :: spcs_density
        type(Array3), intent(in) :: temperature, pressure, density
        type(Array3), allocatable, intent(out):: viscosity, thermal_conductivity
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
