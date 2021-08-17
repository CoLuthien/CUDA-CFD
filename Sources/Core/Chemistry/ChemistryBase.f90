module ChemistryBase
    use, intrinsic :: iso_fortran_env
    use :: SpecieBase, only : Specie
    implicit none

    ! todo : writing down a interface 
    type, abstract :: ChemistrySolver
        class(Specie), allocatable :: spcs(:)
    end type ChemistrySolver
end module
