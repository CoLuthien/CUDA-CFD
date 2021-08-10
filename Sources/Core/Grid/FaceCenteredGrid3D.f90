module FCGrid
    use, intrinsic :: iso_fortran_env
    use :: GridBase
    implicit none

    type, extends(Grids) :: FCGrid3D
    
    real(real64), allocatable :: m_point(:, :, :)
    contains
        procedure, pass :: read => read_grid
        procedure, pass:: interpolate => interpolate
        procedure, pass :: make_volume => calc_volume
    end type FCGrid3D

contains

    subroutine read_grid(self, file_name)
        class(FCGrid3D), intent(in) :: self
        character(len=*), intent(in) :: file_name
    end subroutine read_grid

    subroutine calc_volume(self)
        class(FCGrid3D), intent(in) :: self
    end subroutine calc_volume
    subroutine interpolate(self)
        class(FCGrid3D), intent(in) :: self
    end subroutine interpolate
end module FCGrid
