
module chemtest
    use, intrinsic :: iso_fortran_env
    use :: ChemistryBase


    type, extends(ChemistrySolver) :: CustomChem
    end type
end module
program testSaxpy
    use, intrinsic :: iso_fortran_env
    use :: omp_lib
    use :: SpecieBase
    use :: ArrayBase
    use :: SolverBase
    use :: DiffusionImpl
    use :: chemtest
    use :: FCGrid
    implicit none
    integer, parameter :: n_spc = 13
    real(real64), dimension(:, :, :), allocatable :: x, y, z
    real(real64), allocatable :: k(:), j(:)
    type(Solver) :: fluid_solver
    type(InitialCondition) :: cond
    class(CustomChem), allocatable :: chemsolver
    class(DiffusionLegacy), ALLOCATABLE :: diff_solver
    type(ReferenceState) :: state
    integer :: res(3), i
    integer, parameter :: nx=32*2, ny=17 * 2, nz=27 * 2
    real(real64) :: start, fin

    fluid_solver%n_spc = n_spc

    res= [nx, ny, nz]

    ! read from file
    allocate (x(0:nx + 1, 0:ny+1, 0:nz+1))
    allocate (y, z, source=x)
    x(1:9, 1:9, 1:9) = 2
    y = x
    z = x
    print*, lbound(x), ubound(x)
    allocate(CustomChem ::chemsolver)

    allocate (SpecieNasa7::chemsolver%spcs(n_spc))
    do i=1, n_spc
        call chemsolver%spcs(i)%set_spc_data(n_spc)
    end do 

    do i=1, n_spc
        call chemsolver%spcs(i)%init_spc_derived_data(chemsolver%spcs(1:n_spc), state)
    end do 

    fluid_solver%m_chemistry = chemsolver

    fluid_solver%m_diffusion = DiffusionLegacy(n_spc)
    fluid_solver%m_grid = FCGrid3D(x, y, z, res, n_spc, cond)



    call fluid_solver%check_allocation()
    print*, size(fluid_solver%m_chemistry%spcs), lbound(fluid_solver%m_grid%m_primitives%rho%m_data)
    do i=1, 100
        start = omp_get_wtime()
        call fluid_solver%solve()
        fin = omp_get_wtime()
        print*, fin -start
    end  do


end program testSaxpy

