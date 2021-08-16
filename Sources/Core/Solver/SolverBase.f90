module SolverBase
    use, intrinsic :: iso_fortran_env
    use :: GridBase
    use :: SpecieBase
    use :: AdvectionBase
    use :: DiffusionBase
    implicit none

    type, abstract :: Solver
        class(AdvectionSolver), allocatable :: m_advection
        class(DiffusionSolver), allocatable :: m_diffusion
        class(Grids), allocatable :: m_grid
    contains
    procedure :: solve
    end type Solver
    contains
    subroutine solve(self)
        implicit none
        class(Solver),intent(in) :: self
    
    end subroutine solve
end module
