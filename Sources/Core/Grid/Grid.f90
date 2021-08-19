module GridBase
    use, intrinsic :: iso_fortran_env
    use :: DataGrid
    use :: ArrayBase
    use :: Constants
    implicit none

    type, abstract :: Grid3D(n_spc)
        integer, len :: n_spc=-1
        integer:: m_resolution(3)
        real(real64) :: m_origin(3)
        class(Array3), allocatable :: x, y, z
        class(CellMetricData3D), allocatable :: m_metrics
        class(PrimitiveData3D), allocatable :: m_primitives
        class(ConservedData3D), allocatable :: m_conservatives
    contains
        procedure, public :: set_geometry
        procedure :: set_data
        procedure :: calculate_ghost_point
    end type Grid3D

    private :: set_geometry
contains

    subroutine calculate_ghost_point(self) 
        class(Grid3D(*)) :: self
    end subroutine

    subroutine set_geometry(self, x, y, z)
        class(Grid3D(*)) :: self
        real(real64), dimension(:, :, :), allocatable :: x, y, z

        self%x = Array(x)
        self%y = Array(y)
        self%z = Array(z)
        call self%calculate_ghost_point()

    end subroutine set_geometry

    subroutine set_data(self, data)
        class(Grid3D(*)) :: self
        class(InitialCondition(*)) :: data
    end subroutine
end module GridBase
