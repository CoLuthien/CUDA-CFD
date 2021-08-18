module ChemistryBase
    use, intrinsic :: iso_fortran_env
    use :: SpecieBase, only : Specie
    implicit none

    type, abstract :: ChemistrySolver(n_spc)
        integer, len :: n_spc
        class(Specie), allocatable :: spcs(:)
    end type ChemistrySolver
end module
