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

    type :: MetricData3D
    end type
end module DataGrid