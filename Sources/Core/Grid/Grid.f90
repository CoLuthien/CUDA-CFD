module GridBase
    use :: ArrayBase
    use :: Constants
    use :: DataGrid
    use, intrinsic :: iso_fortran_env
    implicit none

    type, abstract :: Grid3D
        integer :: n_spc
        integer :: m_resolution(3), lb(3), ub(3)
        real(real64) :: m_origin(3)
        type(Array3), allocatable :: x, y, z
        class(CellMetricData3D), allocatable :: m_metrics
        class(PrimitiveData3D), allocatable :: m_primitives
        class(ConservedData3D), allocatable :: m_conservatives
    end type

contains

    subroutine set_geometry(self, x, y, z)
        class(Grid3D) :: self
        real(real64), allocatable, dimension(:, :, :) :: x, y, z
        self%x = Array3D(x)
        self%y = Array3D(y)
        self%z = Array3D(z)
    end subroutine set_geometry

end module GridBase
