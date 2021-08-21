module ChemistryBase
    use, intrinsic :: iso_fortran_env
    use :: SpecieBase, only:Specie
    implicit none

    type, abstract :: ChemistrySolver
        class(Specie), allocatable :: spcs(:)
    end type ChemistrySolver
end module
