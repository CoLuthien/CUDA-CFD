module DiffusionImpl
    use :: DiffusionBase
    implicit none(type, external)

    type, extends(DiffusionSolver) :: DiffusionLegacy
    end type DiffusionLegacy

    type, extends(TransportProperty) :: TransportLegacy
    contains
        procedure :: transport
        !procedure :: viscosity
    end type TransportLegacy

contains
    pure subroutine transport(self &
                              , spcs, spcs_density, temperature, pressure, density &
                              , viscosity, thermal_conductivity)
        use :: Species, only:SpecieBase
        use :: ArrayBase
        class(TransportLegacy), intent(in) :: self
        class(SpecieBase), intent(in) :: spcs(:)
        type(Array4), intent(in) :: spcs_density
        type(Array3), intent(in) :: temperature, pressure, density
        type(Array3), allocatable, intent(out):: viscosity, thermal_conductivity
    end subroutine transport

end module DiffusionImpl
