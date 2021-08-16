module DiffusionImpl
    use, intrinsic :: iso_fortran_env
    use :: DiffusionBase
    implicit none(type, external)

    type, extends(DiffusionSolver) :: DiffusionLegacy
    end type DiffusionLegacy

    type, extends(TransportProperty) :: TransportLegacy
    contains
        procedure :: transport
        !procedure :: viscosity
    end type TransportLegacy

contains
    ! every data in or out for  this routine, is non-dimensional quantity
    pure subroutine transport(self, spcs, spcs_density, mole_fraction, temperature, pressure, density & ! intent(in)
                              , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_quantity) ! intent(out)
        use :: Species, only:SpecieBase
        use :: ArrayBase
        class(TransportLegacy), intent(in) :: self
        class(SpecieBase), intent(in) :: spcs(:)
        class(Array4), intent(in) :: spcs_density, mole_fraction
        class(Array3), intent(in) :: temperature, pressure, density
        class(Array4), intent(out) :: diffusion_quantity
        class(Array3), intent(out) :: viscosity, thermal_conductivity, turbulent_viscosity
        real(real64) :: temp
        real(real64) :: eu_k(size(spcs)), cd_k(size(spcs)), fi(size(spcs), size(spcs))
        integer :: i, j, k, l, m


        do concurrent(i=1:100, j=1:100, k=1:100) local(temp, eu_k, cd_k, fi)
            temp = temperature%get_data(i, j, k)
            eu_k = abs(spcs(:)%Eu(temp))
            cd_k = spcs(:)%Cd(temp)
            do l = 1, size(spcs)
                fi(:, l) = (self%fa(:, l) + self%fb(:, l)*sqrt(eu_k(:)/eu_k(l)))**2
            end do

            ! calculate cd, eu 
        end do  


    end subroutine transport

end module DiffusionImpl
