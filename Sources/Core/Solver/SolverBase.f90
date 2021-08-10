module SolverBase
    use, intrinsic :: iso_fortran_env
    use :: GridBase
    implicit none

    type, abstract :: Solver
    class(Grids), allocatable :: m_grid
    contains
    end type Solver
end module
