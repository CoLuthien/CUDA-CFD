module GridBase
    use, intrinsic :: iso_fortran_env
    use :: DataGrid
    use :: ArrayBase
    implicit none

    type, abstract :: Grid3D
        integer, private :: m_resolution(3)
        real(real64), private :: m_origin(3)
        class(Array3), allocatable :: x, y, z
        class(MetricData3D), allocatable :: m_metrics
        class(PrimitiveData3D), allocatable :: m_primitives
        class(ConservedData3D), allocatable :: m_conservatives
    contains
        procedure(read_mesh), deferred, pass :: read
        procedure(interpolate_points), deferred, pass :: interpolate
        procedure(make_vol), deferred, pass :: make_volume
    end type Grid3D

    

    interface
        subroutine read_mesh(self, file_name)
            import Grid3D
            class(Grid3D), intent(in) :: self
            character(len=*), intent(in) :: file_name
        end subroutine read_mesh
        subroutine make_vol(self)
            import Grid3D
            class(Grid3D), intent(in) :: self
        end subroutine make_vol
        subroutine interpolate_points(self)
            import Grid3D
            class(Grid3D), intent(in) :: self
        end subroutine interpolate_points
    end interface
contains
end module GridBase
