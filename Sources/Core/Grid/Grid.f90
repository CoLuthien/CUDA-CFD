module GridBase
    use, intrinsic :: iso_fortran_env
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

    type :: PrimitiveData3D
        class(Array4), allocatable :: rhok ! density for each spc
        class(Array3), allocatable :: u, v, w
        class(Array3), allocatable :: e, t, p, a, rho
        class(Array3), allocatable :: tk, tw, tv ! turbulent kinetic energy, dissipation rate, viscosity
    end type

    type :: ConservedData3D 
        class(Array4), allocatable :: rhok
        class(Array3), allocatable :: u_momentum, v_momentum, w_momentum
        class(Array3), allocatable :: e
        class(Array3), allocatable :: tk, tw
    end type

    type :: MetricData3D
    end type
    

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
