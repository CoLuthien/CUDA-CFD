module SolverBase
    use, intrinsic :: iso_fortran_env
    use :: SpecieBase
    use :: AdvectionBase
    use :: DiffusionBase
    use :: ChemistryBase
    use :: GridBase
    use :: LegacyTransport
    implicit none

    type :: Solver
        integer :: n_spc
        class(AdvectionSolver), allocatable :: m_advection
        class(DiffusionSolver), allocatable :: m_diffusion
        class(ChemistrySolver), allocatable :: m_chemistry
        class(Grid3D), allocatable :: m_grid
    contains
        procedure :: solve
        procedure :: calc_dt
        procedure :: check_allocation
    end type Solver

contains

    subroutine check_allocation(self)
        class(Solver), intent(in) :: self
        if (self%n_spc /= self%m_grid%n_spc) then
            print *, "allocation status different, solver has n_spc:", self%n_spc, "while grid has n_spc:", self%m_grid%n_spc
            error stop
        end if
    end subroutine check_allocation

    subroutine solve(self, n_spc)
        class(Solver), intent(inout) :: self
        integer, intent(in) :: n_spc
        integer :: i, j, k, nx, ny, nz
        real(real64), dimension(n_spc) :: diffusion_coefficient, enthalpy
        real(real64), dimension(n_spc) :: spcs_density, mole_fraction, mass_fraction

        ! todo : add things to solve
        associate (grid => self%m_grid, spcs => self%m_chemistry%spcs, prim => self%m_grid%m_primitives, &
            consrv => grid%m_conservatives, metrics => grid%m_metrics)
        nx = grid%m_metrics%m_resolution(1)
        ny = grid%m_metrics%m_resolution(2)
        nz = grid%m_metrics%m_resolution(3)

        do concurrent(i=1:nx, j=1:ny, k=1:nz)
            !DIR$ NOUNROLL
            spcs_density(1:n_spc) = prim%rhok(1:n_spc)%m_data(i, j, k)
            mole_fraction(1:n_spc) = spcs(1:n_spc)%get_molar_concentration(spcs_density(1:n_spc)) 
            mole_fraction(1:n_spc) = mole_fraction(1:n_spc) / sum(mole_fraction(1:n_spc))
            mass_fraction(1:n_spc) = mole_fraction(1:n_spc) * spcs(1:n_spc)%molar_weight

            call mixture_properties(spcs, mole_fraction, temp, viscosity, thermal_conductivity, n_spc)
            call mixture_diffusivity(spcs, mole_fraction, mass_fraction, temp, pressure,diffusion_coefficient, n_spc)
            call mixture_turbulent_properties(spcs, mole_fraction, mass_fraction, temp, density, tv, turbulent_conductivity, difft, n_spc)

            turbulent_viscosity = viscosity + tv
            thermal_conductivity = thermal_conductivity + turbulent_conductivity
            diffusion_coefficient(1:n_spc) = diffusion_coefficient(1:n_spc) + difft

            !call self%m_diffusion%diffusive_mass(grid%m_metrics, diffusion_coefficient, )

            call self%m_diffusion%solve_diffusion_point(prim, consrv, grid%m_metrics, spcs, i, j, k, 13)
        end do 
        end associate
        !call self%m_advection%solve_advection

    end subroutine solve

    subroutine calc_dt(self, grid, dt)
        class(Solver), intent(in) :: self
        class(Grid3D), intent(in) :: grid
        real(real64), intent(out) :: dt

    end subroutine calc_dt
end module
