
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
    use :: FCGrid
    use :: DiffusionImpl
    use :: chemtest
    implicit none
    integer, parameter :: n_spc = 13
    real(real64), dimension(:, :, :), allocatable :: x, y, z
    real(real64), allocatable :: k(:), j(:)
    type(Solver(n_spc)) :: fluid_solver
    type(InitialCondition(n_spc)) :: cond
    type(CustomChem(n_spc)), allocatable :: chemsolver
    type(DiffusionLegacy(n_spc)) :: diff_solver
    type(FCGrid3D(n_spc)) :: grid
    type(ReferenceState) :: state
    integer :: res(3), i
    integer, parameter :: nx=32*1, ny=17 * 1, nz=27 * 1
    real(real64) :: start, fin

    res= [nx, ny, nz]

    ! read from file
    allocate (x(0:nx + 2, 0:ny+2, 0:nz+2))
    allocate (y, z, mold=x)
    x(1:9, 1:9, 1:9) = 2
    y = x
    z = x
    allocate(chemsolver)


    allocate (SpecieNasa7::chemsolver%spcs(n_spc))
    call chemsolver%spcs(:)%set_spc_data(n_spc)
    do i=1, n_spc
        call chemsolver%spcs(i)%init_spc_derived_data(chemsolver%spcs(:), state)
    end do 

    fluid_solver%m_grid = FCGrid3D(x, y, z, cond, res, n_spc)
    fluid_solver%m_chemistry = chemsolver
    fluid_solver%m_diffusion = DiffusionLegacy(n_spc)



    call fluid_solver%check_allocation()
    print*, size(fluid_solver%m_chemistry%spcs)
    do i=1, 100
        start = omp_get_wtime()
        call fluid_solver%solve()
        fin = omp_get_wtime()
        print*, fin -start
    end  do

    print*, fluid_solver%m_grid%n_spc

contains
    subroutine test()
        type(Array4), allocatable :: arr1, arr3
        type(Array4) :: arr2, arr4, arr5
        !arr2 = Array(100, 10, 10, 20)

    end subroutine test

end program testSaxpy

