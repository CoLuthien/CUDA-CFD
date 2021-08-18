module DataGrid
    use :: ArrayBase
    implicit none
    public
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

    type :: CellMetricData3D
        class(Array3), allocatable :: sx, sy, sz
        class(Array3), allocatable :: ex, ey, ez
        class(Array3), allocatable :: cx, cy, cz
        class(Array3), allocatable :: sxc, syc, szc
        class(Array3), allocatable :: exc, eyc, ezc
        class(Array3), allocatable :: cxc, cyc, czc, dely
    contains
        procedure :: init_cell_metrics
    end type

contains
    subroutine init_cell_metrics(self, x, y, z)
        class(CellMetricData3D) :: self
        class(Array3), intent(in) :: x, y, z
    end subroutine
end module DataGrid
