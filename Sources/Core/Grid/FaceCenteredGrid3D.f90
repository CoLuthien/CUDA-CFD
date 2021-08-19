module FCGrid
    use, intrinsic :: iso_fortran_env
    use :: GridBase
    use :: Constants
    implicit none

    type, extends(Grid3D) :: FCGrid3D

    end type FCGrid3D

    interface FCGrid3D
        module procedure :: init_fcgrid3d
    end interface

contains
    function init_fcgrid3d(x, y, z, cond, resolution, n_spc) result(self)
        integer, intent(in) :: n_spc
        real(real64), allocatable, dimension(:, :, :) :: x, y, z
        type(FCGrid3D(n_spc)), allocatable :: self
        type(InitialCondition(n_spc)), intent(in) :: cond
        integer, intent(in) :: resolution(3)

        print *, "initializing FCGrid3D..."
        allocate (FCGrid3D(n_spc)::self)

        self%m_resolution = resolution

        print *, "Calculating ghost points..."
        call set_geometry(self, x, y, z)
        print *, "Ghost point calculation done..."

        print *, "Calculating Cell metrics..."
        associate (x => self%x, y => self%y, z => self%z)
            self%m_metrics = CellMetricData3D(x, y, z, resolution)
        end associate
        print *, "Cell metrics calculation done..."

        print *, "Allocating flow field memories..."
        self%m_primitives = PrimitiveData3D(cond, self%x, resolution)
        self%m_conservatives = ConservedData3D(cond, self%x, resolution)

        print *, "Flow field allocation done..."

        print *, "Initializing flow field..."
        call self%set_data(cond)
        print *, "Flow field initialized"

        print *, "FCGrid3D initialized..."
    end function
end module FCGrid
