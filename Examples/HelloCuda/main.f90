
module chemtest
    use, intrinsic :: iso_fortran_env
    use :: ChemistryBase

    type, extends(ChemistrySolver) :: CustomChem
    end type
end module
program testSaxpy
    use, intrinsic :: iso_fortran_env
    use :: SpecieBase
    use :: ArrayBase
    use :: SolverBase
    use :: FCGrid
    use :: DiffusionImpl
    use :: chemtest
    implicit none
    integer, parameter :: n_spc = 13
    real(real64), dimension(:, :, :), allocatable :: x, y, z
    real(real64), allocatable :: k(:), j(:)
    type(Solver(n_spc)) :: fluid_solver
    type(InitialCondition(n_spc)) :: cond
    type(CustomChem(n_spc)) :: chemsolver
    type(DiffusionLegacy(n_spc)) :: diff_solver
    type(FCGrid3D(n_spc)) :: grid
    integer :: res(3)
    integer, parameter :: nx=100, ny=102, nz=300

    res= [nx, ny, nz]

    ! read from file
    allocate (x(-2:nx + 3, -2:ny+3, -2:nz+3))
    allocate (y, z, mold=x)
    x(1:9, 1:9, 1:9) = 2
    y = x
    z = x
    print*, size(x)
    fluid_solver%m_grid = FCGrid3D(x, y, z, cond, res, n_spc)
    fluid_solver%m_chemistry = chemsolver
    fluid_solver%m_diffusion = DiffusionLegacy(n_spc)

    allocate (SpecieNasa7::chemsolver%spcs(n_spc))

    call fluid_solver%check_allocation()

    call fluid_solver%solve()

    print*, fluid_solver%m_grid%n_spc

contains
    subroutine test()
        type(Array4), allocatable :: arr1, arr3
        type(Array4) :: arr2, arr4, arr5
        !arr2 = Array(100, 10, 10, 20)

    end subroutine test

end program testSaxpy

