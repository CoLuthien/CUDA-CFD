module GridBase
    use, intrinsic :: iso_fortran_env
    implicit none(type, external)

    type, abstract :: Grids
        integer, private :: m_resolution(3)
        real(real64), private :: m_origin(3)

    contains
        procedure(read_mesh), deferred, pass :: read
        procedure(interpolate_points), deferred, pass :: interpolate
        procedure(make_vol), deferred, pass :: make_volume

    end type Grids

    interface
        subroutine read_mesh(self, file_name)
            import Grids
            class(Grids), intent(in) :: self
            character(len=*), intent(in) :: file_name
        end subroutine read_mesh
        subroutine make_vol(self)
            import Grids
            class(Grids), intent(in) :: self
        end subroutine make_vol
        subroutine interpolate_points(self)
            import Grids
            class(Grids), intent(in) :: self
        end subroutine interpolate_points
    end interface
contains
end module GridBase
