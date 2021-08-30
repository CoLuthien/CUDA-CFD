module DifferenceBase
    use, intrinsic :: iso_fortran_env
    use :: ArrayBase
    use :: Vector
    implicit none

contains

    pure elemental function central_diff_2nd(target, i, j, k) result(vec)
        class(Array3), intent(in) :: target
        integer, intent(in) :: i, j, k
        type(Vector3) :: vec

        vec = 0.5d0*[target%m_data(i + 1, j, k) - target%m_data(i - 1, j, k), &
                     target%m_data(i, j + 1, k) - target%m_data(i, j - 1, k), &
                     target%m_data(i, j, k + 1) - target%m_data(i, j, k - 1)]
    end function central_diff_2nd


end module
