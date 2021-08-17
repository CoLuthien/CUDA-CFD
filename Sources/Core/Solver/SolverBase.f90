module SolverBase
    use, intrinsic :: iso_fortran_env
    use :: GridBase
    use :: SpecieBase
    use :: AdvectionBase
    use :: DiffusionBase
    use :: ChemistryBase
    implicit none

    type :: Solver
        class(AdvectionSolver), allocatable :: m_advection
        class(DiffusionSolver), allocatable :: m_diffusion
        class(ChemistrySolver), allocatable :: m_chemistry
        class(Grid3D), allocatable :: m_grid
    contains
        procedure :: set_grid
        procedure :: solve
        procedure :: calc_dt
    end type Solver

contains

    subroutine set_grid(self)
        class(Solver), intent(inout) :: self
    end subroutine
    subroutine solve(self)
        class(Solver), intent(in) :: self

        ! todo : add things to solve
        !call self%m_diffusion%solve_diffusion
        !call self%m_advection%solve_advection

    end subroutine solve

    subroutine calc_dt(self, grid, dt)
        class(Solver), intent(in) :: self
        class(Grid3D), intent(in) :: grid
        real(real64), intent(out) :: dt

    end subroutine calc_dt
end module
