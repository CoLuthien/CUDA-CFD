module SolverBase
    use, intrinsic :: iso_fortran_env
    use :: SpecieBase
    use :: AdvectionBase
    use :: DiffusionBase
    use :: ChemistryBase
    use :: GridBase
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

    subroutine solve(self)
        class(Solver), intent(inout) :: self
        integer :: i, j, k, nx, ny, nz

        ! todo : add things to solve
        associate (m_grid => self%m_grid, spcs => self%m_chemistry%spcs)
        nx = m_grid%m_metrics%m_resolution(1)
        ny = m_grid%m_metrics%m_resolution(2)
        nz = m_grid%m_metrics%m_resolution(3)
        do concurrent(i=1:nx, j=1:ny, k=1:nz)
            call self%m_diffusion%solve_diffusion_point(m_grid%m_primitives, m_grid%m_conservatives, m_grid%m_metrics, spcs, i, j, k)
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
