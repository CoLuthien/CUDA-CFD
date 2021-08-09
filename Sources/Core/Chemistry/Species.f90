module Species
    use, intrinsic :: iso_fortran_env
    use :: ProblemSet
    implicit none

   type Specie
      integer :: idx, comp(4)
      character(len=20) :: name ! for debug
      real(real64) :: mole_weight, &
                      molar_hof, mass_hof, & ! heat of formation
                      mole_diam, & !molecular diameter
                      chrt_temp, & ! charateristic temperature
                      poly_cp(7, 2), & ! polynomial for cp fit
                      poly_h(7, 2), & !polynomial fitting for enthalpy
                      poly_s(7, 2) ! polynomial coefficients for entropy
   contains
   end type Specie
end module Species