module Vector
    use, intrinsic :: iso_fortran_env
    implicit none

    type :: IntVect3
        integer :: x, y, z
    end type IntVect3
    !Cartesian Vector type
    type :: Vector3 
        real(real64) :: x = 0., &
                        y = 0., &
                        z = 0.
    end type Vector3

    interface operator(*)
        module procedure scalar_vec_mult
        module procedure vec_scalar_mult
        module procedure iscalar_vec_mult
    end interface operator(*)

    interface operator(/)
        module procedure scalar_vec_div
        module procedure iscalar_vec_div
    end interface operator(/)

    interface operator(+)
        module procedure vec_vec_add
        module procedure scalar_vec_add
        module procedure iscalar_vec_add
    end interface operator(+)

    interface operator(-)
        module procedure vec_vec_sub
        module procedure scalar_vec_sub
        module procedure iscalar_vec_sub
    end interface operator(-)

    interface operator(.dot.)
        module procedure dotproduct
    end interface operator(.dot.)

    interface operator(.cross.)
        module procedure crossproduct
    end interface operator(.cross.)
    public:: operator(*)
    public:: operator(/)
    public:: operator(+)
    public:: operator(-)
    ! public:: operator(/=)
    ! public:: operator(<)
    ! public:: operator(<=)
    ! public:: operator(==)
    ! public:: operator(>=)
    ! public:: operator(>)
    interface assignment(=)
        module procedure assign_vec
        module procedure assign_scalar
        module procedure assign_3p
    end interface assignment(=)

    ! interface size
    !     module procedure magnitude
    ! end interface size

    interface Vector3
        procedure :: assign_3pt
    end interface

contains

    pure elemental function magnitude(self) result(m)
        type(Vector3), intent(in) :: self
        real(real64) :: m
        m = self.dot.self
    end function magnitude

    pure elemental subroutine assign_vec(vec1, vec2)
        type(Vector3), intent(out) :: vec1
        type(Vector3), intent(in) :: vec2
        vec1%x = vec2%x
        vec1%y = vec2%y
        vec1%z = vec2%z
    end subroutine assign_vec

    pure elemental function assign_3pt(x, y, z) result(vec)
        type(Vector3) :: vec
        real(real64), intent(in) :: x, y, z
        vec%x = x
        vec%y = y
        vec%z = z
    end function

    pure subroutine assign_3p(vec, pt)
        type(Vector3), intent(inout) :: vec
        real(real64), intent(in) :: pt(3)
        vec%x = pt(1)
        vec%y = pt(2)
        vec%z = pt(3)
    end subroutine assign_3p

    pure elemental subroutine assign_scalar(vec, scalar)
        type(Vector3), intent(inout) :: vec
        real(real64), intent(in) :: scalar
        vec%x = scalar
        vec%y = scalar
        vec%z = scalar
    end subroutine assign_scalar

    pure elemental function scalar_vec_mult(scalar, vec) result(mul)
        type(Vector3), intent(in) :: vec
        real(real64), intent(in) :: scalar
        type(Vector3) :: mul
        mul%x = vec%x * scalar
        mul%y = vec%y * scalar
        mul%z = vec%z * scalar
    end function scalar_vec_mult

    pure elemental function vec_scalar_mult(vec, scalar) result(mul)
        type(Vector3), intent(in) :: vec
        real(real64), intent(in) :: scalar
        type(Vector3) :: mul
        mul%x = vec%x * scalar
        mul%y = vec%y * scalar
        mul%z = vec%z * scalar
    end function vec_scalar_mult

    pure elemental function iscalar_vec_mult(vec, scalar) result(mul)
        type(Vector3), intent(in) :: vec
        integer, intent(in) :: scalar
        type(Vector3) :: mul
        mul%x = vec%x * real(scalar, real64)
        mul%y = vec%y * real(scalar, real64)
        mul%z = vec%z * real(scalar, real64)
    end function iscalar_vec_mult

    pure elemental function scalar_vec_div(vec, scalar) result(div)
        type(Vector3), intent(in) :: vec
        real(real64), intent(in) :: scalar
        type(Vector3) :: div
        div%x = vec%x / scalar
        div%y = vec%y / scalar
        div%z = vec%z / scalar
    end function scalar_vec_div

    pure elemental function iscalar_vec_div(vec, scalar) result(div)
        type(Vector3), intent(in) :: vec
        integer, intent(in) :: scalar
        type(Vector3) :: div
        div%x = vec%x / real(scalar, real64)
        div%y = vec%y / real(scalar, real64)
        div%z = vec%z / real(scalar, real64)
    end function iscalar_vec_div

    pure elemental function vec_vec_add(vec1, vec2) result(add)
        !$omp declare simd(vec_vec_add) uniform(vec1, vec2)
        type(Vector3), intent(in) :: vec1, vec2
        type(Vector3) :: add
        add%x = vec1%x + vec2%x
        add%y = vec1%y + vec2%y
        add%z = vec1%z + vec2%z
    end function vec_vec_add

    pure elemental function scalar_vec_add(vec, scalar) result(add)
        type(Vector3), intent(in) :: vec
        real(real64), intent(in) :: scalar
        type(Vector3) :: add
        add%x = vec%x + scalar
        add%y = vec%y + scalar
        add%z = vec%z + scalar
    end function scalar_vec_add

    pure elemental function iscalar_vec_add(vec, scalar) result(add)
        type(Vector3), intent(in) :: vec
        integer, intent(in) :: scalar
        type(Vector3) :: add
        add%x = vec%x + real(scalar, real64)
        add%y = vec%y + real(scalar, real64)
        add%z = vec%z + real(scalar, real64)
    end function iscalar_vec_add

    pure elemental function vec_vec_sub(vec1, vec2) result(add)
        !$omp declare simd(vec_vec_sub) uniform(vec1, vec2)
        type(Vector3), intent(in) :: vec1, vec2
        type(Vector3) :: add
        add%x = vec1%x - vec2%x
        add%y = vec1%y - vec2%y
        add%z = vec1%z - vec2%z
    end function vec_vec_sub

    pure elemental function scalar_vec_sub(vec, scalar) result(add)
        type(Vector3), intent(in) :: vec
        real(real64), intent(in) :: scalar
        type(Vector3) :: add
        add%x = vec%x - scalar
        add%y = vec%y - scalar
        add%z = vec%z - scalar
    end function scalar_vec_sub

    pure elemental function iscalar_vec_sub(vec, scalar) result(add)
        type(Vector3), intent(in) :: vec
        integer, intent(in) :: scalar
        type(Vector3) :: add
        add%x = vec%x - real(scalar, real64)
        add%y = vec%y - real(scalar, real64)
        add%z = vec%z - real(scalar, real64)
    end function iscalar_vec_sub

    pure elemental function dotproduct(vec1, vec2) result(dot)
        !$omp declare simd(dotproduct) uniform(vec1, vec2)
        type(Vector3), intent(in) :: vec1, vec2
        real(real64) :: dot
        dot = vec1%x * vec2%x + &
              vec1%y * vec2%y + &
              vec1%z * vec2%z
    end function dotproduct

    !DIR$ ATTRIBUTES INLINE :: crossproduct
    pure elemental function crossproduct(vec1, vec2) result(cross)
        !$omp declare simd(crossproduct) uniform(vec1, vec2)
        type(Vector3), intent(in) :: vec1, vec2
        type(Vector3) :: cross
        cross%x = (vec1%y * vec2%z) - (vec1%z * vec2%y)
        cross%y = (vec1%z * vec2%x) - (vec1%x * vec2%z)
        cross%z = (vec1%x * vec2%y) - (vec1%y * vec2%x)
    end function crossproduct
end module Vector
