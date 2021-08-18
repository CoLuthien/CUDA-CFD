module FCGrid
    use, intrinsic :: iso_fortran_env
    use :: GridBase
    implicit none

    type, extends(Grid3D) :: FCGrid3D

    end type FCGrid3D

    interface FCGrid3D
        procedure :: init_fcgrid3d
    end interface

contains
    function init_fcgrid3d(x, y, z, data) result(self)
        real(real64), allocatable, dimension(:, :, :) :: x, y, z
        type(FCGrid3D(:)), allocatable :: self
        class(InitialCondition(*)), intent(in) :: data

        allocate (FCGrid3D(data%n_spc)::self)
        call self%set_geometry(x, y, z)
        call self%set_data(data)
        ! todo => calculate cell metrics..

    end function
end module FCGrid
