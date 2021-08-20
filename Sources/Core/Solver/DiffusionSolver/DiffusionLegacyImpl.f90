module DiffusionImpl
    use, intrinsic :: iso_fortran_env
    use :: DiffusionBase, only:DiffusionSolver, TransportProperty
    use :: ArrayBase
    use :: SpecieBase
    implicit none

    type, extends(DiffusionSolver) :: DiffusionLegacy
    contains
        !procedure :: solve_diffusion => diffusion
        procedure :: diffusive_flux => legacy_diffusive_flux
    end type DiffusionLegacy

    type, extends(TransportProperty) :: TransportLegacy
    contains
        procedure :: transport
    end type TransportLegacy

    interface DiffusionLegacy
        module procedure :: init_diffusion_legacy
    end interface

contains
    ! todo: refine init diffusion
    function init_diffusion_legacy(n_spc) result(solver)
        integer :: n_spc
        type(DiffusionLegacy(n_spc)) :: solver
        type(TransportLegacy) :: prop

        solver%m_transport = prop
    end function

    subroutine diffusion(self)
        class(DiffusionLegacy(*)), intent(in) :: self
        call self%diffusive_flux()
    end subroutine

    subroutine legacy_diffusive_flux(self)
        class(DiffusionLegacy(*)), intent(in) :: self
    end subroutine legacy_diffusive_flux

    ! every data in or out for  this routine, is non-dimensional quantity
    subroutine transport(self, spcs, spcs_density, temperature, pressure, density & ! intent(in)
                         , viscosity, turbulent_viscosity, thermal_conductivity, diffusion_coefficient & ! intent(out)
                         , resolution)
        class(TransportLegacy), intent(in) :: self
        class(Specie), intent(in) :: spcs(:)
        class(Array4), intent(in) :: spcs_density
        class(Array3), intent(in) :: temperature, pressure, density
        integer, intent(in) :: resolution(3)
        type(Array4), intent(inout) :: diffusion_coefficient
        type(Array3), intent(inout) :: viscosity, thermal_conductivity, turbulent_viscosity
        real(real64) :: temp, xsm, eu, cd, weight, p, t3p
        real(real64) :: td(size(spcs)), omg(size(spcs)), df(size(spcs), size(spcs))
        real(real64), dimension(size(spcs)) :: ratio, eu_k, cd_k, c
        integer :: i, j, k, l, m, n_spc, nx, ny, nz, cnt, idx
        n_spc = size(spcs)
        cnt = 0
        nx = resolution(1); ny = resolution(2); nz = resolution(3)
        ! loop for calculating laminar property
        ! this loop quite difficult to understand,
        ! to do: change loop range into realistic one
        do concurrent(i=1:nx, j=1:ny, k=1:nz) !local_init(temp, eu_k, cd_k, weight, eu, cd, xsm, ratio, c)
            ! calculatie mole fraction, non-dimensional
            c(1:n_spc) = spcs_density%m_data(i, j, k, 1:n_spc)*spcs(1:n_spc)%molar_weight
            c(1:n_spc) = c(1:n_spc)/sum(c(1:n_spc))

            temp = temperature%get_data(i, j, k) ! get flow field temperature
            do idx=1, n_spc
                eu_k(idx) = spcs(idx)%Eu(temp) ! calculate single component, laminar viscosity coeff
            end do
            cd_k = spcs(:)%Cd(temp) ! calculate single component, laminar thermal conductivity coeff
            do m = 1, n_spc
                ratio(1:n_spc) = eu_k(m)/eu_k(1:n_spc) ! get viscosity coeff ratio for this spc
                weight = spcs(m)%get_scale_factor(ratio(1:), c(1:)) ! get weighting factor for this spc
                xsm = c(m)/weight
                eu = eu + xsm*eu_k(m)
                cd = cd + xsm*cd_k(m)
            end do

            viscosity%m_data(i, j, k) = eu
            thermal_conductivity%m_data(i, j, k) = cd
        end do
        ! todo : calculate diffusion coefficient and turbulent viscosity

        do concurrent(i=1:nx, j=1:ny, k=1:nz) local(temp, t3p, p)
            temp = temperature%get_data(i, j, k)
            p = pressure%get_data(i, j, k)
            t3p = sqrt((temp**3)/(p**2))
            do m = 1, n_spc
                td(1:n_spc) = temp*spcs(m)%teab(1:)
                omg(1:n_spc) = 1.d0/(td(1:)**0.145d0) &
                          + 1.d0/ ((td + 0.5d0) ** 2)
                df(m, 1:n_spc) = (spcs(m)%diff_coef(1:n_spc) * t3p) / omg(1:)
            end do

        end do

    end subroutine transport

end module DiffusionImpl
