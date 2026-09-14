!> Shared data for the tox_normalization suites: a matrix whose genes lie exactly on a linear
!| mean-sd trend. LOESS reproduces a straight line exactly, so normalizing by the fitted sd has a
!| closed form, and the LOESS-based cases can compare against hand-derived values.
module mod_test_tox_normalization_fixtures
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    private
    public :: fill_linear_trend, linear_trend_offset, linear_trend_sd

contains

    !> Offset d_j of replicate `i_replicate` from its gene's mean, in units of the gene index:
    !| j - (n + 1)/2, so the n offsets sum to zero.
    pure real(real64) function linear_trend_offset(i_replicate, n_replicates)
        integer(int32), intent(in) :: i_replicate
            !! Replicate index, 1..n_replicates
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene

        linear_trend_offset = real(i_replicate, real64) - 0.5_real64*real(n_replicates + 1, real64)
    end function linear_trend_offset

    !> The population sd of the offsets, sqrt(sum(d**2)/n) = sqrt((n**2 - 1)/12). Gene i of
    !| `fill_linear_trend` has sd i times this.
    pure real(real64) function linear_trend_sd(n_replicates)
        integer(int32), intent(in) :: n_replicates
            !! Number of replicates per gene

        linear_trend_sd = sqrt(real(n_replicates**2 - 1, real64)/12.0_real64)
    end function linear_trend_sd

    !> expr(j, i) = i*(10 + d_j): gene i has mean 10*i and sd i*linear_trend_sd(n), so every gene
    !| lies on the line sd = mean*linear_trend_sd(n)/10, and normalizing gene i by its fitted sd
    !| gives (10 + d_j)/linear_trend_sd(n), the same for every gene.
    pure subroutine fill_linear_trend(expr)
        real(real64), intent(out) :: expr(:, :)
            !! (n_replicates, n_genes) matrix to fill

        integer(int32) :: i_replicate, i_gene, n_replicates

        n_replicates = size(expr, 1, kind=int32)
        do i_gene = 1, size(expr, 2, kind=int32)
            do i_replicate = 1, n_replicates
                expr(i_replicate, i_gene) = real(i_gene, real64)* &
                                            (10.0_real64 + linear_trend_offset(i_replicate, n_replicates))
            end do
        end do
    end subroutine fill_linear_trend

end module mod_test_tox_normalization_fixtures
