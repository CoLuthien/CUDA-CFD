module DiffusionBase
    use, intrinsic :: iso_fortran_env
    use :: Debug
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
        pure subroutine calc_transport_property(self, spcs, spcs_density, temperature, pressure, density, tv & ! intent(in)
                                           , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient &
                                           ) ! intent(out)
            use :: SpecieBase, only:Specie
            use :: ArrayBase
            import TransportProperty
            class(TransportProperty), intent(in) :: self
        class(Specie), intent(in) :: spcs(:)
        real(real64), intent(in) :: spcs_density(size(spcs))
        real(real64), intent(in) :: temperature, pressure, density, tv
        real(real64), intent(out) :: diffusion_coefficient(size(spcs))
        real(real64), intent(out) :: viscosity, thermal_conductivity, turbulent_viscosity
        end subroutine calc_transport_property
    end interface

    contains

    subroutine solve_diffusion(self, prim, conserv, metrics, spcs)
        use :: ArrayBase
        use :: GridBase, only:Grid3D, PrimitiveData3D, ConservedData3D, CellMetricData3D
        use :: SpecieBase, only:Specie
        class(DiffusionSolver(*)), intent(in) :: self
        class(CellMetricData3D), intent(in) :: metrics
        class(PrimitiveData3D), intent(in), allocatable :: prim
        class(ConservedData3D), intent(inout) :: conserv
        class(Specie), intent(in) :: spcs(self%n_spc)
        real(real64) :: diffusion_coefficient(size(spcs))
        real(real64) :: viscosity, thermal_conductivity, turbulent_viscosity
        integer :: nx, ny, nz, i, j, k
        nx = prim%m_resolution(1)
        ny = prim%m_resolution(2)
        nz = prim%m_resolution(3)
        do concurrent(i=1:nx, j=1:ny, k=1:nz)
            call self%m_transport%transport(spcs, prim%rhok%m_data(i, j, k, 1:), prim%t%m_data(i, j, k),&
             prim%p%m_data(i, j, k), prim%rho%m_data(i, j, k), prim%tv%m_data(i, j, k) &
                                   , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient)
        end do  
        !call self%diffusive_flux
    end subroutine

end module DiffusionBase
