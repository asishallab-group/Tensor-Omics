#include <src/macros.h>

!> summary: Implementation fixture for the compile-and-run test
!| author: A Developer
!| One implementation exercising a documented minimum, a work array and a permutation, so that both
!| generated wrappers are meaningful. Compiled and run by `test_end_to_end_fortran.py`.
module fx_ranks_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private

    public :: rank_scores_impl

contains

    !> summary: Shifts the scores by `min_score` and returns them in ascending order
    !| author: A Developer
    pure subroutine rank_scores_impl(n_scores, scores, scores_perm, tmp_shifted, min_score, ranked)
        integer(int32), intent(in) :: n_scores
            !! number of elements in `scores`
        real(real64), dimension(n_scores), intent(in) :: scores
            !! the values to rank
        integer(int32), dimension(n_scores), intent(in) :: scores_perm
            !! permutation that sorts `scores` ascending
        real(real64), dimension(n_scores), intent(out) :: tmp_shifted
            !! work array holding `scores - min_score`
        real(real64), intent(in) :: min_score
            !! the value subtracted from every score
            !! DM_MIN(0.0_real64)
        real(real64), dimension(n_scores), intent(out) :: ranked
            !! the shifted scores, ascending

        integer(int32) :: i_score

        do i_score = 1, n_scores
            tmp_shifted(i_score) = scores(i_score) - min_score
        end do
        do i_score = 1, n_scores
            ranked(i_score) = tmp_shifted(scores_perm(i_score))
        end do
    end subroutine rank_scores_impl

end module fx_ranks_impl
