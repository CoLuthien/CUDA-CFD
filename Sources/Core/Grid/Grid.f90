module GridBase
    use, intrinsic :: iso_fortran_env
    use :: DataGrid
    use :: ArrayBase
    use :: Constants
    implicit none

    type, abstract :: Grid3D(n_spc)
    integer, len :: n_spc = -1
    integer:: m_resolution(3)
    real(real64) :: m_origin(3)
    class(Array3), allocatable :: x, y, z
    class(CellMetricData3D), allocatable :: m_metrics
    class(PrimitiveData3D), allocatable :: m_primitives
    class(ConservedData3D), allocatable :: m_conservatives
contains
    procedure, public :: set_geometry
    procedure :: set_data
    procedure ::  calculate_ghost_points
    end type Grid3D

    private :: set_geometry
    contains

    subroutine calculate_ghost_points(self)
        class(Grid3D(*)) :: self
        type(Vector3), allocatable :: pts(:, :, :)
        type(Vector3) :: point
        integer :: nx, ny, nz, i, j, k
        integer :: ib, jb, kb, ie, je, ke
        integer :: ibm, jbm, kbm, ibp, jbp, kbp
        integer :: iem, jem, kem, iep, jep, kep
        integer :: lb(3), ub(3)

        nx = self%m_resolution(1)
        ny = self%m_resolution(2)
        nz = self%m_resolution(3)

        ib = 1
        jb = 1
        kb = 1
        ie = nx
        je = ny
        ke = nz
        ibm = ib - 1
        ibp = ib + 1
        iem = ie - 1
        iep = ie + 1
        jbm = jb - 1
        jbp = jb + 1
        jem = je - 1
        jep = je + 1
        kbm = kb - 1
        kbp = kb + 1
        kem = ke - 1
        kep = ke + 1
        lb = lbound(self%x%m_data)
        ub = ubound(self%x%m_data)

        allocate (pts(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
        print *, ub
        !allocate (pts, mold=reshape([Vector3::], &
        !                            self%m_resolution))

        ! 1: for reducing un-necessary operation
        pts(1:, 1:, 1:) = Vector3(self%x%m_data(1:, 1:, 1:), &
                                  self%y%m_data(1:, 1:, 1:), &
                                  self%z%m_data(1:, 1:, 1:))

        !---- Ghost Mesh Points of Bottom/Top Surface
        do concurrent(i=1:nx, j=1:ny) local(point)
            pts(i, j, kbm) = 2.d0*pts(i, j, ib) - pts(i, j, ibp)

            pts(i, j, kep) = 2.d0*pts(i, j, ke) - pts(i, j, kem)
        end do

!---- Ghost Mesh Points of Right/Left Surface
        do concurrent(i=1:nx, k=1:nz) local(point)
            pts(i, jbm, k) = 2.d0*pts(i, jb, k) - pts(i, jbp, k)

            pts(i, jep, k) = 2.d0*pts(i, je, k) - pts(i, jem, k)
        end do
        !---- Ghost Mesh Points of Front/Rear Surface
        do concurrent(j=jb:je, k=kb:ke)
            pts(ibm, j, k) = 2.d0*pts(ib, j, k) - pts(ibp, j, k)

            pts(iep, j, k) = 2.d0*pts(ie, j, k) - pts(iem, j, k)

        end do

!---- Ghost Mesh Points of i-direction edge

        do concurrent(i=ib:ie)
            pts(i, jbm, kbm) = (2.d0 &
                                *(2.d0*pts(i, jb, kb) - pts(i, jbp, kbp)) &
                                + 2.d0*pts(i, jbm, kb) - pts(i, jbm, kbp) &
                                + 2.d0*pts(i, jb, kbm) - pts(i, jbp, kbm))/4.d0

            pts(i, jbm, kep) = (2.d0 &
                                *(2.d0*pts(i, jb, ke) - pts(i, jbp, kem)) &
                                + 2.d0*pts(i, jbm, ke) - pts(i, jbm, kem) &
                                + 2.d0*pts(i, jb, kep) - pts(i, jbp, kep))/4.d0

            pts(i, jep, kep) = (2.d0 &
                                *(2.d0*pts(i, je, ke) - pts(i, jem, kem)) &
                                + 2.d0*pts(i, jep, ke) - pts(i, jep, kem) &
                                + 2.d0*pts(i, je, kep) - pts(i, jem, kep))/4.d0

            pts(i, jep, kbm) = (2.d0 &
                                *(2.d0*pts(i, je, kb) - pts(i, jem, kbp)) &
                                + 2.d0*pts(i, jep, kb) - pts(i, jep, kbp) &
                                + 2.d0*pts(i, je, kbm) - pts(i, jem, kbm))/4.d0
        end do
        !---- Ghost Mesh Points of j-direction edge
        do concurrent(j=jb:je)
            pts(ibm, j, kbm) = (2.d0 &
                                *(2.d0*pts(ib, j, kb) - pts(ibp, j, kbp)) &
                                + 2.d0*pts(ibm, j, kb) - pts(ibm, j, kbp) &
                                + 2.d0*pts(ib, j, kbm) - pts(ibp, j, kbm))/4.d0

            pts(iep, j, kbm) = (2.d0 &
                                *(2.d0*pts(ie, j, kb) - pts(iem, j, kbp)) &
                                + 2.d0*pts(iep, j, kb) - pts(iep, j, kbp) &
                                + 2.d0*pts(ie, j, kbm) - pts(iem, j, kbm))/4.d0

            pts(iep, j, kep) = (2.d0 &
                                *(2.d0*pts(ie, j, ke) - pts(iem, j, kem)) &
                                + 2.d0*pts(iep, j, ke) - pts(iep, j, kem) &
                                + 2.d0*pts(ie, j, kep) - pts(iem, j, kep))/4.d0

            pts(ibm, j, kep) = (2.d0 &
                                *(2.d0*pts(ib, j, ke) - pts(ibp, j, kem)) &
                                + 2.d0*pts(ibm, j, ke) - pts(ibm, j, kem) &
                                + 2.d0*pts(ib, j, kep) - pts(ibp, j, kep))/4.d0
        end do
        !---- Ghost Mesh Points of k-direction edge
        do concurrent(k=kb:ke)
            pts(ibm, jbm, k) = (2.d0 &
                                *(2.d0*pts(ib, jb, k) - pts(ibp, jbp, k)) &
                                + 2.d0*pts(ib, jbm, k) - pts(ibp, jbm, k) &
                                + 2.d0*pts(ibm, jb, k) - pts(ibm, jbp, k))/4.d0

            pts(iep, jbm, k) = (2.d0 &
                                *(2.d0*pts(ie, jb, k) - pts(iem, jbp, k)) &
                                + 2.d0*pts(ie, jbm, k) - pts(iem, jbm, k) &
                                + 2.d0*pts(iep, jb, k) - pts(iep, jbp, k))/4.d0

            pts(iep, jep, k) = (2.d0 &
                                *(2.d0*pts(ie, je, k) - pts(iem, jem, k)) &
                                + 2.d0*pts(ie, jep, k) - pts(iem, jep, k) &
                                + 2.d0*pts(iep, je, k) - pts(iep, jem, k))/4.d0

            pts(ibm, jep, k) = (2.d0 &
                                *(2.d0*pts(ib, je, k) - pts(ibp, jem, k)) &
                                + 2.d0*pts(ib, jep, k) - pts(ibp, jep, k) &
                                + 2.d0*pts(ibm, je, k) - pts(ibm, jem, k))/4.d0
        end do
        !---- Ghost Mesh Points of vertex
        pts(ibm, jbm, kbm) = (3.d0*(2.d0*pts(ib, jb, kb) - pts(ibp, jbp, kbp)) &
                              + 2.d0*pts(ibm, jbm, kb) - pts(ibm, jbm, kbp) &
                              + 2.d0*pts(ib, jbm, kbm) - pts(ibp, jbm, kbm) &
                              + 2.d0*pts(ibm, jb, kbm) - pts(ibm, jbp, kbm))/6.d0

        pts(iep, jbm, kbm) = (3.d0*(2.d0*pts(ie, jb, kb) - pts(iem, jbp, kbp)) &
                              + 2.d0*pts(iep, jbm, kb) - pts(iep, jbm, kbp) &
                              + 2.d0*pts(ie, jbm, kbm) - pts(iem, jbm, kbm) &
                              + 2.d0*pts(iep, jb, kbm) - pts(iep, jbp, kbm))/6.d0

        pts(iep, jep, kbm) = (3.d0*(2.d0*pts(ie, je, kb) - pts(iem, jem, kbp)) &
                              + 2.d0*pts(iep, jep, kb) - pts(iep, jep, kbp) &
                              + 2.d0*pts(ie, jep, kbm) - pts(iem, jep, kbm) &
                              + 2.d0*pts(iep, je, kbm) - pts(iep, jem, kbm))/6.d0

        pts(ibm, jep, kbm) = (3.d0*(2.d0*pts(ib, je, kb) - pts(ibp, jem, kbp)) &
                              + 2.d0*pts(ibm, jep, kb) - pts(ibm, jep, kbp) &
                              + 2.d0*pts(ib, jep, kbm) - pts(ibp, jep, kbm) &
                              + 2.d0*pts(ibm, je, kbm) - pts(ibm, jem, kbm))/6.d0

        pts(ibm, jbm, kep) = (3.d0*(2.d0*pts(ib, jb, ke) - pts(ibp, jbp, kem)) &
                              + 2.d0*pts(ibm, jbm, ke) - pts(ibm, jbm, kem) &
                              + 2.d0*pts(ib, jbm, kep) - pts(ibp, jbm, kep) &
                              + 2.d0*pts(ibm, jb, kep) - pts(ibm, jbp, kep))/6.d0

        pts(iep, jbm, kep) = (3.d0*(2.d0*pts(ie, jb, ke) - pts(iem, jbp, kem)) &
                              + 2.d0*pts(iep, jbm, ke) - pts(iep, jbm, kem) &
                              + 2.d0*pts(ie, jbm, kep) - pts(iem, jbm, kep) &
                              + 2.d0*pts(iep, jb, kep) - pts(iep, jbp, kep))/6.d0

        pts(iep, jep, kep) = (3.d0*(2.d0*pts(ie, je, ke) - pts(iem, jem, kem)) &
                              + 2.d0*pts(iep, jep, ke) - pts(iep, jep, kem) &
                              + 2.d0*pts(ie, jep, kep) - pts(iem, jep, kep) &
                              + 2.d0*pts(iep, je, kep) - pts(iep, jem, kep))/6.d0

        pts(ibm, jep, kep) = (3.d0*(2.d0*pts(ib, je, ke) - pts(ibp, jem, kem)) &
                              + 2.d0*pts(ibm, jep, ke) - pts(ibm, jep, kem) &
                              + 2.d0*pts(ib, jep, kep) - pts(ibp, jep, kep) &
                              + 2.d0*pts(ibm, je, kep) - pts(ibm, jem, kep))/6.d0

        self%x%m_data(:, :, :) = pts(:, :, :)%x
        self%y%m_data(:, :, :) = pts(:, :, :)%y
        self%z%m_data(:, :, :) = pts(:, :, :)%z

        deallocate (pts)

    end subroutine calculate_ghost_points

    subroutine set_geometry(self, x, y, z)
        class(Grid3D(*)) :: self
        real(real64), dimension(:, :, :), allocatable :: x, y, z

        self%x = Array(x)
        self%y = Array(y)
        self%z = Array(z)
        call self%calculate_ghost_points()

    end subroutine set_geometry

    subroutine set_data(self, data)
        class(Grid3D(*)) :: self
        class(InitialCondition(*)) :: data
    end subroutine
end module GridBase
