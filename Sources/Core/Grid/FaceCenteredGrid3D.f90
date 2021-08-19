module FCGrid
    use, intrinsic :: iso_fortran_env
    use :: GridBase
    implicit none

    type, extends(Grid3D) :: FCGrid3D

        contains 
    end type FCGrid3D

    interface FCGrid3D
        procedure :: init_fcgrid3d
    end interface

contains
    function init_fcgrid3d(x, y, z, data, resolution) result(self)
        real(real64), allocatable, dimension(:, :, :) :: x, y, z
        type(FCGrid3D(n_spc=:)), allocatable :: self
        class(InitialCondition(*)), intent(in) :: data
        integer, intent(in) :: resolution(3)

        print*, "initializing FCGrid3D..."
        allocate (FCGrid3D(data%n_spc)::self)
        self%m_resolution = resolution

        call self%set_geometry(x, y, z)
        call self%set_data(data)
        ! todo => calculate cell metrics..
        associate (x => self%x, y => self%y, z => self%z)
            self%m_metrics = CellMetricData3D(x, y, z, resolution)
        end associate

        print*, "FCGrid3D initialized..."
    end function
end module FCGrid
