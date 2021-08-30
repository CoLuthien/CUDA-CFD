module Constants
    use, intrinsic :: iso_fortran_env
    implicit none
    public
    integer, parameter :: n_face = 6
    real(real64), parameter :: pi = 4.d0 * atan(1.d0)

    real(real64) :: Re_inf_mach1 = 0.d0
    real(real64) :: prandtl_number = 0.9d0, schmidt_number = 0.9d0

    real(real64) :: ref_kinematic_viscosity

    ! SST k-omega model constants, 1 for original k-w model, 2 for transformed k-e model 
    real(real64), parameter :: kw_sigma_k1 = .85d0, kw_sigma_w1 = .5d0, &
                               kw_sigma_k2 = 1.d0, kw_sigma_w2 = .856d0, &
                               kw_beta_star = .09d0, kw_beta1 = 3.d0 / 40.d0, kw_beta2 = .0828d0, &
                               kw_kappa = .41d0, kw_a1 = 0.31d0, &
                               kw_gamma1 = kw_beta1 / kw_beta_star - ((kw_sigma_w1 * (kw_kappa**2)) / sqrt(kw_beta_star)), &
                               kw_gamma2 = kw_beta2 / kw_beta_star - ((kw_sigma_w2 * (kw_kappa**2)) / sqrt(kw_beta_star))

    type :: InitialCondition
        real(real64) :: u, v, w, temp, tk, tw, tv
        real(real64), allocatable :: spcs_density(:)
    end type

end module
