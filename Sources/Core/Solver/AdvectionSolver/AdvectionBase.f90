module AdvectionBase
    use, intrinsic :: iso_fortran_env
    implicit none

    type, abstract :: AdvectionSolver(n_spc)
        integer, len :: n_spc
    end type 
end module AdvectionBase