!> Unit test suite for the multidimensional noise model (`noise_model_md`).
!|
!| Covers the ten checks of `md_noise_model_plan_independent_axes.md` §10:
!|
!|   1. `d = 1` on the ENUMERATED path reproduces `noise_model_exact`'s
!|      `pvalues_own` EXACTLY (both normalization branches). Both sides are
!|      deterministic there, so this is an exact-match test, not a tolerance test.
!|      It is the load-bearing one: the whole module is justified by reducing to
!|      the validated scalar path. It also pins the PER-SIDE `1/sqrt(n)` scaling,
!|      because the fixture is deliberately asymmetric (`n_case = 4`,
!|      `n_ctrl = 5`) — a single per-axis factor would fail it outright.
!|   2. `compute_beta_md` at `d = 1` reproduces the observed statistic the scalar
!|      pipelines are fed, on the linear AND the log2 (Frechet) branch.
!|   3. Both norms, `D` and `D_std`, against a hand-computed fixture.
!|   4. The closed-form `Var_null` against the mean of the enumerated null. Under
!|      full enumeration `mean(D^2_null) = sum_a [Var_null(a) + (mu_case,a -
!|      mu_ctrl,a)^2]`, and the pool means are exactly 0 (a pool is a union of
!|      whole genes, each of whose residuals sum to 0), so the two must agree.
!|   5. Ragged `n_case = [4,6,7,6,5]` against `n_ctrl = [3,5,5,9,3]`: pools are
!|      whole multiples of that axis's replicate count, and `Var_null(a)` tracks
!|      `1/n_case(a) + 1/n_ctrl(a)` across axes.
!|   6. `enum_max_product` routes correctly and `method_used` reports it.
!|   7. Systematic and RNG sampling agree within Monte-Carlo error.
!|   8. `beta_centre = 1` subtracts the per-axis median; a constant offset on one
!|      axis moves `delta_hat` by exactly that offset and leaves the centred
!|      statistic unchanged.
!|   9. Degenerate input: a permanently thin axis, a zero-variance axis, `D = 0`,
!|      a non-finite supplied `beta`, and the argument-validation error codes.
!|  10. `beta_mode = 1` with an off-scale `beta` -> correlation ~1, large MAD, run
!|      completes (the units-problem signature).
module mod_test_noise_model_md
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use noise_model, only: NULL_METHOD_POOLED
    use noise_model_exact, only: pipeline_exact => compute_noise_pvalue_pipeline
    use noise_model_md
    use tox_errors
    use test_suite, only: test_case
    implicit none
    public

    real(real64), parameter :: TOL = 1d-12
    real(real64), parameter :: EXACT = 0.0_real64

    ! Neighbourhood settings shared by every test that runs the full pipeline.
    integer(int32), parameter :: K_START = 10_int32
    integer(int32), parameter :: K_STEP = 2_int32
    integer(int32), parameter :: K_MAX = 30_int32
    real(real64), parameter :: TAU = 0.1_real64
    integer(int32), parameter :: MAX_POOL = 500_int32
    integer(int32), parameter :: SEED = 42_int32

    ! Big enough that any d = 1 fixture here enumerates; small enough that the
    ! d >= 2 fixtures sample unless they are deliberately tiny.
    integer(int64), parameter :: ENUM_BIG = 100000000_int64
    integer(int64), parameter :: ENUM_OFF = 0_int64

contains

    !> Get array of all available tests.
    function get_all_tests_noise_model_md() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate(all_tests(10))
        all_tests(1)  = test_case("test_md_reduces_to_exact_d1",     test_md_reduces_to_exact_d1)
        all_tests(2)  = test_case("test_compute_beta_md_d1",         test_compute_beta_md_d1)
        all_tests(3)  = test_case("test_both_norms_fixture",         test_both_norms_fixture)
        all_tests(4)  = test_case("test_var_null_closed_form",       test_var_null_closed_form)
        all_tests(5)  = test_case("test_ragged_axes",                test_ragged_axes)
        all_tests(6)  = test_case("test_enum_routing",               test_enum_routing)
        all_tests(7)  = test_case("test_systematic_vs_rng",          test_systematic_vs_rng)
        all_tests(8)  = test_case("test_beta_centre_median",         test_beta_centre_median)
        all_tests(9)  = test_case("test_degenerate_input",           test_degenerate_input)
        all_tests(10) = test_case("test_beta_mode_supplied_offscale", test_beta_mode_supplied_offscale)
    end function get_all_tests_noise_model_md

    ! =========================================================================
    ! helpers
    ! =========================================================================

    !> Integer to string, for assertion messages.
    function str(i) result(s)
        integer(int32), intent(in) :: i
        character(len=20) :: buffer
        character(len=:), allocatable :: s
        write(buffer, '(I0)') i
        s = trim(buffer)
    end function str

    !> Total length of a packed replicate buffer.
    pure function packed_len(n_genes, n_axes, n_rep) result(n)
        integer(int32), intent(in) :: n_genes, n_axes
        integer(int32), intent(in) :: n_rep(n_axes)
        integer(int32) :: n, a
        n = 0
        do a = 1, n_axes
            n = n + n_rep(a) * n_genes
        end do
    end function packed_len

    !> Build one group's packed replicates and matching per-axis means.
    !|
    !| Deterministic, replicate- and axis-dependent noise (no RNG, so the fixture is
    !| identical on every platform). `shift` displaces the whole group, which is how
    !| a case/control difference is introduced. `means` is set to the ACTUAL
    !| arithmetic replicate mean, matching how the R caller builds it.
    subroutine make_group(n_genes, n_axes, n_rep, phase, shift, noise_scale, means, packed)
        integer(int32), intent(in) :: n_genes, n_axes
        integer(int32), intent(in) :: n_rep(n_axes)
        real(real64), intent(in) :: phase, shift, noise_scale
        real(real64), intent(out) :: means(n_genes, n_axes)
        real(real64), intent(out) :: packed(:)

        integer(int32) :: a, g, r, off, base
        real(real64) :: mu, acc

        off = 0
        do a = 1, n_axes
            do g = 1, n_genes
                mu = 5.0_real64 + 10.0_real64 * real(g - 1, real64) / real(max(n_genes - 1, 1), real64) &
                     + 0.5_real64 * real(a - 1, real64) + shift
                base = off + (g - 1) * n_rep(a)
                acc = 0.0_real64
                do r = 1, n_rep(a)
                    packed(base + r) = mu + noise_scale * &
                        sin(phase + 1.7_real64 * real(g, real64) &
                            + 0.9_real64 * real(r, real64) + 2.3_real64 * real(a, real64))
                    acc = acc + packed(base + r)
                end do
                means(g, a) = acc / real(n_rep(a), real64)
            end do
            off = off + n_rep(a) * n_genes
        end do
    end subroutine make_group

    !> Build a group whose per-gene residual variance is EXACTLY known.
    !|
    !| Each gene's replicates are `mu + c*(2*(r-1)/(n-1) - 1)`: `n` points evenly
    !| spaced over `[mu-c, mu+c]`, so the gene mean is exactly `mu` and the
    !| population variance of the deviations is `c^2 (n+1) / (3(n-1))`.
    !| `prepare_sorted_data_helper` then applies the Bessel factor `sqrt(n/(n-1))`,
    !| so a pool built from any number of whole genes has variance exactly
    !| `c^2 n(n+1) / (3(n-1)^2)` and mean exactly 0 — which makes `Var_null(a)`
    !| predictable in closed form and turns the per-side `1/sqrt(n)` scaling into an
    !| exact assertion rather than a tolerance.
    subroutine make_group_ramp(n_genes, n_axes, n_rep, c, means, packed)
        integer(int32), intent(in) :: n_genes, n_axes
        integer(int32), intent(in) :: n_rep(n_axes)
        real(real64), intent(in) :: c
        real(real64), intent(out) :: means(n_genes, n_axes)
        real(real64), intent(out) :: packed(:)

        integer(int32) :: a, g, r, off, base
        real(real64) :: mu

        off = 0
        do a = 1, n_axes
            do g = 1, n_genes
                mu = 5.0_real64 + 10.0_real64 * real(g - 1, real64) / real(max(n_genes - 1, 1), real64) &
                     + 0.5_real64 * real(a - 1, real64)
                base = off + (g - 1) * n_rep(a)
                do r = 1, n_rep(a)
                    packed(base + r) = mu + c * (2.0_real64 * real(r - 1, real64) &
                                                 / real(n_rep(a) - 1, real64) - 1.0_real64)
                end do
                means(g, a) = mu
            end do
            off = off + n_rep(a) * n_genes
        end do
    end subroutine make_group_ramp

    !> Variance of one side\'s SCALED pool for the `make_group_ramp` fixture:
    !| `c^2 n(n+1)/(3(n-1)^2)` divided by `n`, i.e. `c^2 (n+1)/(3(n-1)^2)`.
    pure function ramp_scaled_var(n_rep, c) result(v)
        integer(int32), intent(in) :: n_rep
        real(real64), intent(in) :: c
        real(real64) :: v
        v = c * c * real(n_rep + 1, real64) &
            / (3.0_real64 * real(n_rep - 1, real64) ** 2)
    end function ramp_scaled_var

    !> Run the multidimensional pipeline with the suite's shared settings.
    !|
    !| Thin wrapper so each test states only what it actually varies.
    subroutine run_md(n_genes, n_axes, n_rep_case, n_rep_control, &
                      means_case, packed_case, means_control, packed_control, &
                      norm_method, beta_mode, beta_centre, &
                      sampling_mode, n_draws_max, n_exceed_target, enum_max_product, &
                      beta_obs, pvalues, d_obs, d_std_obs, d_sq_null_mean, &
                      method_used, n_draws_used, nb_case, nb_control, var_null, &
                      beta_cor, beta_mad, delta_hat, n_with_p, ierr)
        integer(int32), intent(in) :: n_genes, n_axes
        integer(int32), intent(in) :: n_rep_case(n_axes), n_rep_control(n_axes)
        real(real64), intent(in) :: means_case(n_genes, n_axes), means_control(n_genes, n_axes)
        real(real64), intent(in) :: packed_case(:), packed_control(:)
        integer(int32), intent(in) :: norm_method, beta_mode, beta_centre
        integer(int32), intent(in) :: sampling_mode, n_draws_max, n_exceed_target
        integer(int64), intent(in) :: enum_max_product
        real(real64), intent(inout) :: beta_obs(n_axes, n_genes)
        real(real64), intent(out) :: pvalues(n_genes), d_obs(n_genes), d_std_obs(n_genes)
        real(real64), intent(out) :: d_sq_null_mean(n_genes)
        integer(int32), intent(out) :: method_used(n_genes), n_draws_used(n_genes)
        integer(int32), intent(out) :: nb_case(n_axes, n_genes), nb_control(n_axes, n_genes)
        real(real64), intent(out) :: var_null(n_axes, n_genes)
        real(real64), intent(out) :: beta_cor(n_axes), beta_mad(n_axes), delta_hat(n_axes)
        integer(int32), intent(out) :: n_with_p, ierr

        integer(int32) :: valid(n_genes)

        valid = 1
        call compute_noise_pvalue_pipeline_md( &
            means_case, packed_case, n_rep_case, &
            means_control, packed_control, n_rep_control, &
            beta_obs, valid, beta_mode, beta_centre, &
            n_genes, n_axes, norm_method, K_START, K_STEP, K_MAX, TAU, 0.0_real64, &
            NULL_METHOD_POOLED, sampling_mode, n_draws_max, n_exceed_target, enum_max_product, &
            SEED, MAX_POOL, &
            pvalues, d_obs, d_std_obs, d_sq_null_mean, method_used, n_draws_used, &
            nb_case, nb_control, var_null, beta_cor, beta_mad, delta_hat, n_with_p, ierr)
    end subroutine run_md

    ! =========================================================================
    ! 1. d = 1, enumerated, reduces to noise_model_exact exactly
    ! =========================================================================
    subroutine test_md_reduces_to_exact_d1()
        integer(int32), parameter :: n_genes = 60, n_axes = 1
        integer(int32) :: n_rep_case(n_axes), n_rep_control(n_axes)
        real(real64) :: means_case(n_genes, n_axes), means_control(n_genes, n_axes)
        real(real64), allocatable :: packed_case(:), packed_control(:)
        real(real64) :: beta_obs(n_axes, n_genes)
        real(real64) :: p_md(n_genes), d_obs(n_genes), d_std(n_genes), d_sq_null(n_genes)
        integer(int32) :: method_used(n_genes), n_draws_used(n_genes)
        integer(int32) :: nb_case(n_axes, n_genes), nb_control(n_axes, n_genes)
        real(real64) :: var_null(n_axes, n_genes)
        real(real64) :: beta_cor(n_axes), beta_mad(n_axes), delta_hat(n_axes)
        integer(int32) :: n_with_p, ierr

        real(real64), allocatable :: rep_case(:, :), rep_control(:, :)
        real(real64) :: p_exact(n_genes), obs_own(n_genes)
        integer(int32) :: valid(n_genes)
        integer(int32) :: nb_own_case(n_genes), nb_own_control(n_genes), nb_scalar(n_genes)
        integer(int32) :: n_with_p_exact, ierr_exact, i_norm, norm_method, g

        ! Deliberately asymmetric: a single per-axis 1/sqrt(n) factor instead of the
        ! required per-SIDE one would change every null value and fail this test.
        n_rep_case = [4]
        n_rep_control = [5]
        allocate(packed_case(packed_len(n_genes, n_axes, n_rep_case)))
        allocate(packed_control(packed_len(n_genes, n_axes, n_rep_control)))
        allocate(rep_case(n_rep_case(1), n_genes))
        allocate(rep_control(n_rep_control(1), n_genes))

        call make_group(n_genes, n_axes, n_rep_case, 0.0_real64, 1.2_real64, 0.30_real64, &
                        means_case, packed_case)
        call make_group(n_genes, n_axes, n_rep_control, 1.1_real64, 0.0_real64, 0.25_real64, &
                        means_control, packed_control)

        ! At d = 1 the packed buffer IS the (n_rep, n_genes) matrix the scalar
        ! pipeline takes -- same element order, no repacking.
        rep_case = reshape(packed_case, [n_rep_case(1), n_genes])
        rep_control = reshape(packed_control, [n_rep_control(1), n_genes])

        valid = 1

        do i_norm = 0, 1
            norm_method = i_norm

            ! The exact pipeline is fed exactly the statistic compute_beta_md
            ! produces, so the only remaining difference would be the null itself.
            call compute_beta_md(packed_case, n_rep_case, packed_control, n_rep_control, &
                                 n_genes, n_axes, norm_method, beta_obs, ierr)
            call assert_equal_int(ierr, ERR_OK, "d1: compute_beta_md ierr (norm "//str(norm_method)//")")
            obs_own = beta_obs(1, :)

            call pipeline_exact( &
                means_case(:, 1), rep_case, n_genes, n_rep_case(1), &
                means_control(:, 1), rep_control, n_genes, n_rep_control(1), &
                obs_own, valid, n_genes, norm_method, K_START, K_STEP, K_MAX, TAU, &
                0.0_real64, p_exact, n_with_p_exact, MAX_POOL, &
                nb_own_case, nb_own_control, nb_scalar, ierr_exact)
            call assert_equal_int(ierr_exact, ERR_OK, &
                                  "d1: exact pipeline ierr (norm "//str(norm_method)//")")

            call run_md(n_genes, n_axes, n_rep_case, n_rep_control, &
                        means_case, packed_case, means_control, packed_control, &
                        norm_method, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                        SAMPLING_SYSTEMATIC, 1000, 0, ENUM_BIG, &
                        beta_obs, p_md, d_obs, d_std, d_sq_null, method_used, n_draws_used, &
                        nb_case, nb_control, var_null, beta_cor, beta_mad, delta_hat, n_with_p, ierr)
            call assert_equal_int(ierr, ERR_OK, "d1: md pipeline ierr (norm "//str(norm_method)//")")

            call assert_equal_int(n_with_p, n_with_p_exact, &
                                  "d1: gene count with p-value differs (norm "//str(norm_method)//")")
            call assert_equal_array_real(p_md, p_exact, n_genes, EXACT, &
                                         "d1: p-values must match noise_model_exact EXACTLY (norm " &
                                         //str(norm_method)//")")
            call assert_equal_array_int(nb_case(1, :), nb_own_case, n_genes, &
                                        "d1: case pool sizes differ (norm "//str(norm_method)//")")
            call assert_equal_array_int(nb_control(1, :), nb_own_control, n_genes, &
                                        "d1: control pool sizes differ (norm "//str(norm_method)//")")

            do g = 1, n_genes
                if (p_md(g) < 0.0_real64) cycle
                call assert_equal_int(method_used(g), METHOD_ENUMERATED, &
                                      "d1: gene "//str(g)//" should have enumerated")
                ! At d = 1 the unweighted norm is |beta|.
                call assert_equal_real(d_obs(g), abs(beta_obs(1, g)), EXACT, &
                                       "d1: D must equal |beta| at gene "//str(g))
                ! ... and the enumerated walk visits every case/control pair exactly once.
                call assert_equal_int(n_draws_used(g), nb_case(1, g) * nb_control(1, g), &
                                      "d1: enumerated atom count wrong at gene "//str(g))
            end do
        end do

        call assert_true(n_with_p > 0, "d1: no gene received a p-value")

        deallocate(packed_case, packed_control, rep_case, rep_control)
    end subroutine test_md_reduces_to_exact_d1

    ! =========================================================================
    ! 2. compute_beta_md at d = 1, both normalization branches
    ! =========================================================================
    subroutine test_compute_beta_md_d1()
        integer(int32), parameter :: n_genes = 12, n_axes = 1
        integer(int32) :: n_rep_case(n_axes), n_rep_control(n_axes)
        real(real64) :: means_case(n_genes, n_axes), means_control(n_genes, n_axes)
        real(real64), allocatable :: packed_case(:), packed_control(:)
        real(real64) :: beta_obs(n_axes, n_genes)
        real(real64) :: expected, sum_c, sum_t, log2_factor
        integer(int32) :: ierr, g, r, bc, bt

        n_rep_case = [3]
        n_rep_control = [4]
        allocate(packed_case(packed_len(n_genes, n_axes, n_rep_case)))
        allocate(packed_control(packed_len(n_genes, n_axes, n_rep_control)))
        call make_group(n_genes, n_axes, n_rep_case, 0.4_real64, 0.7_real64, 0.4_real64, &
                        means_case, packed_case)
        call make_group(n_genes, n_axes, n_rep_control, 2.0_real64, 0.0_real64, 0.4_real64, &
                        means_control, packed_control)

        log2_factor = 1.0_real64 / log(2.0_real64)

        ! --- linear branch: plain difference of arithmetic replicate means -------
        call compute_beta_md(packed_case, n_rep_case, packed_control, n_rep_control, &
                             n_genes, n_axes, 0, beta_obs, ierr)
        call assert_equal_int(ierr, ERR_OK, "beta d1 linear: ierr")
        do g = 1, n_genes
            bc = (g - 1) * n_rep_case(1)
            bt = (g - 1) * n_rep_control(1)
            sum_c = 0.0_real64
            do r = 1, n_rep_case(1)
                sum_c = sum_c + packed_case(bc + r)
            end do
            sum_t = 0.0_real64
            do r = 1, n_rep_control(1)
                sum_t = sum_t + packed_control(bt + r)
            end do
            expected = sum_c / real(n_rep_case(1), real64) - sum_t / real(n_rep_control(1), real64)
            call assert_equal_real(beta_obs(1, g), expected, EXACT, &
                                   "beta d1 linear: gene "//str(g))
            ! `means` was built as the arithmetic replicate mean, so it must agree.
            call assert_equal_real(beta_obs(1, g), means_case(g, 1) - means_control(g, 1), TOL, &
                                   "beta d1 linear: must equal the mean difference at gene "//str(g))
        end do

        ! --- log2 branch: difference of FRECHET means ---------------------------
        ! mean_i log2(x_i + c), NOT log2(mean_i x_i + c). This is the term the
        ! residual centring in prepare_sorted_data_helper uses, and the reason this
        ! routine takes the replicates rather than the plan's `means_*`.
        call compute_beta_md(packed_case, n_rep_case, packed_control, n_rep_control, &
                             n_genes, n_axes, 1, beta_obs, ierr)
        call assert_equal_int(ierr, ERR_OK, "beta d1 log2: ierr")
        do g = 1, n_genes
            bc = (g - 1) * n_rep_case(1)
            bt = (g - 1) * n_rep_control(1)
            sum_c = 0.0_real64
            do r = 1, n_rep_case(1)
                sum_c = sum_c + log(max(packed_case(bc + r), 0.0_real64) + 1.0_real64)
            end do
            sum_t = 0.0_real64
            do r = 1, n_rep_control(1)
                sum_t = sum_t + log(max(packed_control(bt + r), 0.0_real64) + 1.0_real64)
            end do
            expected = sum_c * log2_factor / real(n_rep_case(1), real64) &
                     - sum_t * log2_factor / real(n_rep_control(1), real64)
            call assert_equal_real(beta_obs(1, g), expected, EXACT, &
                                   "beta d1 log2: gene "//str(g))
        end do

        deallocate(packed_case, packed_control)
    end subroutine test_compute_beta_md_d1

    ! =========================================================================
    ! 3. both norms against a hand-computed fixture
    ! =========================================================================
    !| `beta_hat` on a small ragged fixture, then the two norms of §2 built from the
    !| `var_null` the pipeline itself reports — so the test pins the relationship
    !| `D = sqrt(sum beta^2)`, `D_std = sqrt(sum beta^2/Var_null)` independently of
    !| how `Var_null` happens to come out.
    subroutine test_both_norms_fixture()
        integer(int32), parameter :: n_genes_f = 2, n_axes = 2
        integer(int32) :: n_rep_case(n_axes), n_rep_control(n_axes)
        real(real64) :: packed_case(2 * 2 + 3 * 2), packed_control(3 * 2 + 2 * 2)
        real(real64) :: beta_obs(n_axes, n_genes_f)
        real(real64) :: log2_factor, e11, e12, e21, e22
        integer(int32) :: ierr

        integer(int32), parameter :: n_genes = 40
        integer(int32) :: nrc(n_axes), nrt(n_axes)
        real(real64) :: mc(n_genes, n_axes), mt(n_genes, n_axes)
        real(real64), allocatable :: pc(:), pt(:)
        real(real64) :: b(n_axes, n_genes)
        real(real64) :: p(n_genes), d_obs(n_genes), d_std(n_genes), s(n_genes)
        integer(int32) :: meth(n_genes), nd(n_genes)
        integer(int32) :: nbc(n_axes, n_genes), nbt(n_axes, n_genes)
        real(real64) :: vn(n_axes, n_genes), cc(n_axes), mm(n_axes), dh(n_axes)
        integer(int32) :: nwp, g, a
        real(real64) :: exp_d_sq, exp_d_std_sq

        n_rep_case = [2, 3]
        n_rep_control = [3, 2]

        ! axis 1 case: gene1 = (10, 14) -> 12 ; gene2 = (20, 30) -> 25
        ! axis 2 case: gene1 = (1, 2, 3) -> 2  ; gene2 = (4, 5, 9) -> 6
        packed_case = [10.0_real64, 14.0_real64, 20.0_real64, 30.0_real64, &
                        1.0_real64,  2.0_real64,  3.0_real64, &
                        4.0_real64,  5.0_real64,  9.0_real64]
        ! axis 1 ctrl: gene1 = (7, 8, 9) -> 8   ; gene2 = (1, 2, 3) -> 2
        ! axis 2 ctrl: gene1 = (0, 4) -> 2      ; gene2 = (6, 10) -> 8
        packed_control = [7.0_real64, 8.0_real64, 9.0_real64, &
                          1.0_real64, 2.0_real64, 3.0_real64, &
                          0.0_real64, 4.0_real64, &
                          6.0_real64, 10.0_real64]

        call compute_beta_md(packed_case, n_rep_case, packed_control, n_rep_control, &
                             n_genes_f, n_axes, 0, beta_obs, ierr)
        call assert_equal_int(ierr, ERR_OK, "norms fixture linear: ierr")
        call assert_equal_real(beta_obs(1, 1), 12.0_real64 - 8.0_real64, TOL, "norms fixture: axis 1 gene 1")
        call assert_equal_real(beta_obs(1, 2), 25.0_real64 - 2.0_real64, TOL, "norms fixture: axis 1 gene 2")
        call assert_equal_real(beta_obs(2, 1),  2.0_real64 - 2.0_real64, TOL, "norms fixture: axis 2 gene 1")
        call assert_equal_real(beta_obs(2, 2),  6.0_real64 - 8.0_real64, TOL, "norms fixture: axis 2 gene 2")

        ! log2 branch, hand-written Frechet means with c = NOISE_LOG_OFFSET = 1.
        log2_factor = 1.0_real64 / log(2.0_real64)
        e11 = (log(11.0_real64) + log(15.0_real64)) * log2_factor / 2.0_real64 &
            - (log(8.0_real64) + log(9.0_real64) + log(10.0_real64)) * log2_factor / 3.0_real64
        e12 = (log(21.0_real64) + log(31.0_real64)) * log2_factor / 2.0_real64 &
            - (log(2.0_real64) + log(3.0_real64) + log(4.0_real64)) * log2_factor / 3.0_real64
        e21 = (log(2.0_real64) + log(3.0_real64) + log(4.0_real64)) * log2_factor / 3.0_real64 &
            - (log(1.0_real64) + log(5.0_real64)) * log2_factor / 2.0_real64
        e22 = (log(5.0_real64) + log(6.0_real64) + log(10.0_real64)) * log2_factor / 3.0_real64 &
            - (log(7.0_real64) + log(11.0_real64)) * log2_factor / 2.0_real64

        call compute_beta_md(packed_case, n_rep_case, packed_control, n_rep_control, &
                             n_genes_f, n_axes, 1, beta_obs, ierr)
        call assert_equal_int(ierr, ERR_OK, "norms fixture log2: ierr")
        call assert_equal_real(beta_obs(1, 1), e11, TOL, "norms fixture log2: axis 1 gene 1")
        call assert_equal_real(beta_obs(1, 2), e12, TOL, "norms fixture log2: axis 1 gene 2")
        call assert_equal_real(beta_obs(2, 1), e21, TOL, "norms fixture log2: axis 2 gene 1")
        call assert_equal_real(beta_obs(2, 2), e22, TOL, "norms fixture log2: axis 2 gene 2")

        ! --- the two norms, on a fixture big enough to actually run --------------
        nrc = [4, 3]
        nrt = [3, 4]
        allocate(pc(packed_len(n_genes, n_axes, nrc)))
        allocate(pt(packed_len(n_genes, n_axes, nrt)))
        call make_group(n_genes, n_axes, nrc, 0.5_real64, 1.1_real64, 0.30_real64, mc, pc)
        call make_group(n_genes, n_axes, nrt, 1.9_real64, 0.0_real64, 0.30_real64, mt, pt)

        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, 2000, 0, ENUM_OFF, &
                    b, p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "norms: ierr")
        call assert_true(nwp > 0, "norms: no gene received a p-value")

        do g = 1, n_genes
            if (p(g) < 0.0_real64) cycle
            exp_d_sq = 0.0_real64
            exp_d_std_sq = 0.0_real64
            do a = 1, n_axes
                exp_d_sq = exp_d_sq + b(a, g) * b(a, g)
                call assert_true(vn(a, g) > 0.0_real64, &
                                 "norms: var_null must be positive, axis "//str(a)//" gene "//str(g))
                exp_d_std_sq = exp_d_std_sq + b(a, g) * b(a, g) / vn(a, g)
            end do
            call assert_equal_real(d_obs(g), sqrt(exp_d_sq), TOL, &
                                   "norms: D is not the unweighted norm at gene "//str(g))
            call assert_equal_real(d_std(g), sqrt(exp_d_std_sq), TOL, &
                                   "norms: D_std is not the Var_null-standardised norm at gene "//str(g))
            ! The two are genuinely different statistics -- if they coincided the
            ! standardisation would be doing nothing.
            call assert_true(abs(d_obs(g) - d_std(g)) > 1.0e-9_real64, &
                             "norms: D and D_std coincide at gene "//str(g))
        end do

        deallocate(pc, pt)
    end subroutine test_both_norms_fixture

    ! =========================================================================
    ! 4. Var_null closed form vs the enumerated null's own second moment
    ! =========================================================================
    !| Under FULL enumeration the mean of `D^2_null` is exactly
    !| `sum_a [ Var(pool_case,a) + Var(pool_ctrl,a) + (mu_case,a - mu_ctrl,a)^2 ]`
    !| on the scaled pools, i.e. `sum_a Var_null(a)` plus the squared pool-mean gap.
    !| That gap is zero here: a pool is a union of whole genes and each gene's
    !| residuals sum to exactly zero by construction, so both pool means vanish.
    !| The identity is therefore exact, and it checks the closed form against the
    !| null the module actually walks.
    subroutine test_var_null_closed_form()
        integer(int32), parameter :: n_genes = 30, n_axes = 2
        integer(int32) :: nrc(n_axes), nrt(n_axes)
        real(real64) :: mc(n_genes, n_axes), mt(n_genes, n_axes)
        real(real64), allocatable :: pc(:), pt(:)
        real(real64) :: b(n_axes, n_genes)
        real(real64) :: p(n_genes), d_obs(n_genes), d_std(n_genes), s(n_genes)
        integer(int32) :: meth(n_genes), nd(n_genes)
        integer(int32) :: nbc(n_axes, n_genes), nbt(n_axes, n_genes)
        real(real64) :: vn(n_axes, n_genes), cc(n_axes), mm(n_axes), dh(n_axes)
        integer(int32) :: nwp, ierr, g, a, n_checked
        real(real64) :: sum_vn

        ! Small pools (k_max is capped by the gene count here) so the joint space
        ! stays inside ENUM_BIG and the enumerated path is actually taken.
        nrc = [3, 3]
        nrt = [3, 3]
        allocate(pc(packed_len(n_genes, n_axes, nrc)))
        allocate(pt(packed_len(n_genes, n_axes, nrt)))
        call make_group(n_genes, n_axes, nrc, 0.3_real64, 0.8_real64, 0.30_real64, mc, pc)
        call make_group(n_genes, n_axes, nrt, 1.4_real64, 0.0_real64, 0.30_real64, mt, pt)

        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, 1000, 0, ENUM_BIG, &
                    b, p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "var_null: ierr")
        call assert_true(nwp > 0, "var_null: no gene received a p-value")

        n_checked = 0
        do g = 1, n_genes
            if (p(g) < 0.0_real64) cycle
            call assert_equal_int(meth(g), METHOD_ENUMERATED, &
                                  "var_null: gene "//str(g)//" must enumerate for this identity")
            sum_vn = 0.0_real64
            do a = 1, n_axes
                sum_vn = sum_vn + vn(a, g)
            end do
            ! Relative tolerance: the two are the same sum accumulated in a different
            ! order, so only rounding separates them.
            call assert_true(abs(s(g) - sum_vn) <= 1.0e-9_real64 * max(1.0_real64, abs(sum_vn)), &
                             "var_null: mean(D^2_null) /= sum_a Var_null(a) at gene "//str(g))
            n_checked = n_checked + 1
        end do
        call assert_true(n_checked > 0, "var_null: nothing checked")

        deallocate(pc, pt)
    end subroutine test_var_null_closed_form

    ! =========================================================================
    ! 5. ragged replicate counts, both sides, five axes
    ! =========================================================================
    subroutine test_ragged_axes()
        integer(int32), parameter :: n_genes = 50, n_axes = 5
        integer(int32) :: nrc(n_axes), nrt(n_axes)
        real(real64) :: mc(n_genes, n_axes), mt(n_genes, n_axes)
        real(real64), allocatable :: pc(:), pt(:)
        real(real64) :: b(n_axes, n_genes)
        real(real64) :: p(n_genes), d_obs(n_genes), d_std(n_genes), s(n_genes)
        integer(int32) :: meth(n_genes), nd(n_genes)
        integer(int32) :: nbc(n_axes, n_genes), nbt(n_axes, n_genes)
        real(real64) :: vn(n_axes, n_genes), cc(n_axes), mm(n_axes), dh(n_axes)
        integer(int32) :: nwp, ierr, g, a, n_checked
        real(real64) :: expect_vn(n_axes)
        real(real64), parameter :: C_DEV = 0.4_real64

        ! The plan's fixture. Asymmetric per side AND ragged across axes: a single
        ! per-axis factor would be wrong on every axis here.
        nrc = [4, 6, 7, 6, 5]
        nrt = [3, 5, 5, 9, 3]
        allocate(pc(packed_len(n_genes, n_axes, nrc)))
        allocate(pt(packed_len(n_genes, n_axes, nrt)))
        ! Ramp fixture: identical deviation amplitude everywhere, and a residual
        ! variance that is known in closed form for every replicate count. The only
        ! thing that can then move Var_null between axes is the replicate counts.
        call make_group_ramp(n_genes, n_axes, nrc, C_DEV, mc, pc)
        call make_group_ramp(n_genes, n_axes, nrt, C_DEV, mt, pt)

        ! Var_null(a) = Var(pool_case,a)/n_case(a) + Var(pool_ctrl,a)/n_ctrl(a),
        ! each side by its OWN replicate count. Using a single per-axis factor, or
        ! dropping the Bessel correction, changes these numbers.
        do a = 1, n_axes
            expect_vn(a) = ramp_scaled_var(nrc(a), C_DEV) + ramp_scaled_var(nrt(a), C_DEV)
        end do

        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, 2000, 0, ENUM_OFF, &
                    b, p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)

        call assert_equal_int(ierr, ERR_OK, "ragged: ierr")
        call assert_true(nwp > 0, "ragged: no gene received a p-value")

        n_checked = 0
        do g = 1, n_genes
            if (p(g) < 0.0_real64) cycle
            n_checked = n_checked + 1
            call assert_in_range_real(p(g), 1.0_real64 / 2001.0_real64, 1.0_real64, &
                                      "ragged: p out of range at gene "//str(g))
            do a = 1, n_axes
                ! Every neighbour gene contributes exactly n_rep residuals, so a pool
                ! that is not a whole multiple means the gather or the packing is off.
                call assert_true(nbc(a, g) >= 10, &
                                 "ragged: thin case pool, axis "//str(a)//" gene "//str(g))
                call assert_equal_int(mod(nbc(a, g), nrc(a)), 0, &
                                      "ragged: case pool not a multiple of n_rep, axis "//str(a))
                call assert_equal_int(mod(nbt(a, g), nrt(a)), 0, &
                                      "ragged: control pool not a multiple of n_rep, axis "//str(a))
                call assert_true(vn(a, g) > 0.0_real64, &
                                 "ragged: var_null must be positive, axis "//str(a))
            end do

            ! The exact per-axis null variance, each side scaled by its own sqrt(n).
            do a = 1, n_axes
                call assert_true(abs(vn(a, g) - expect_vn(a)) <= 1.0e-9_real64 * expect_vn(a), &
                                 "ragged: Var_null does not match the per-side 1/n closed form, "// &
                                 "axis "//str(a)//" gene "//str(g))
            end do
        end do
        call assert_true(n_checked > 0, "ragged: nothing checked")

        deallocate(pc, pt)
    end subroutine test_ragged_axes

    ! =========================================================================
    ! 6. enum_max_product routing, reported through method_used
    ! =========================================================================
    subroutine test_enum_routing()
        integer(int32), parameter :: n_genes = 24, n_axes = 2
        integer(int32) :: nrc(n_axes), nrt(n_axes)
        real(real64) :: mc(n_genes, n_axes), mt(n_genes, n_axes)
        real(real64), allocatable :: pc(:), pt(:)
        real(real64) :: b(n_axes, n_genes)
        real(real64) :: p(n_genes), d_obs(n_genes), d_std(n_genes), s(n_genes)
        integer(int32) :: meth(n_genes), nd(n_genes)
        integer(int32) :: nbc(n_axes, n_genes), nbt(n_axes, n_genes)
        real(real64) :: vn(n_axes, n_genes), cc(n_axes), mm(n_axes), dh(n_axes)
        integer(int32) :: nwp, ierr, g, n_enum, n_sys, n_rng

        nrc = [3, 3]
        nrt = [3, 3]
        allocate(pc(packed_len(n_genes, n_axes, nrc)))
        allocate(pt(packed_len(n_genes, n_axes, nrt)))
        call make_group(n_genes, n_axes, nrc, 0.1_real64, 0.6_real64, 0.30_real64, mc, pc)
        call make_group(n_genes, n_axes, nrt, 1.2_real64, 0.0_real64, 0.30_real64, mt, pt)

        ! (a) a cap above the joint space -> every gene enumerates
        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, 500, 0, ENUM_BIG, &
                    b, p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "routing: enumerated ierr")
        n_enum = 0
        do g = 1, n_genes
            if (p(g) < 0.0_real64) cycle
            n_enum = n_enum + 1
            call assert_equal_int(meth(g), METHOD_ENUMERATED, &
                                  "routing: gene "//str(g)//" should have enumerated")
            ! The enumerated walk visits prod_a S_a atoms, and reports it.
            call assert_equal_int(nd(g), (nbc(1, g) * nbt(1, g)) * (nbc(2, g) * nbt(2, g)), &
                                  "routing: enumerated atom count wrong at gene "//str(g))
        end do
        call assert_true(n_enum > 0, "routing: nothing enumerated")

        ! (b) cap 0 forces sampling, and the mode picks which sampler
        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, 500, 0, ENUM_OFF, &
                    b, p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "routing: systematic ierr")
        n_sys = 0
        do g = 1, n_genes
            if (p(g) < 0.0_real64) cycle
            n_sys = n_sys + 1
            call assert_equal_int(meth(g), METHOD_SYSTEMATIC, &
                                  "routing: gene "//str(g)//" should have sampled systematically")
            call assert_equal_int(nd(g), 500, "routing: fixed-M run must use all draws")
        end do
        call assert_equal_int(n_sys, n_enum, "routing: gene count changed with the path")

        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_RNG, 500, 0, ENUM_OFF, &
                    b, p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "routing: rng ierr")
        n_rng = 0
        do g = 1, n_genes
            if (p(g) < 0.0_real64) cycle
            n_rng = n_rng + 1
            call assert_equal_int(meth(g), METHOD_RNG, &
                                  "routing: gene "//str(g)//" should have sampled by RNG")
        end do
        call assert_equal_int(n_rng, n_enum, "routing: gene count changed with the sampler")

        ! (c) Besag-Clifford stops early and reports the reduced draw count.
        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, 500, 20, ENUM_OFF, &
                    b, p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "routing: sequential ierr")
        do g = 1, n_genes
            if (p(g) < 0.0_real64) cycle
            call assert_true(nd(g) <= 500, "routing: sequential run exceeded M at gene "//str(g))
            call assert_in_range_real(p(g), 1.0_real64 / 501.0_real64, 1.0_real64, &
                                      "routing: sequential p out of range at gene "//str(g))
        end do

        deallocate(pc, pt)
    end subroutine test_enum_routing

    ! =========================================================================
    ! 7. systematic and RNG sampling agree within Monte-Carlo error
    ! =========================================================================
    subroutine test_systematic_vs_rng()
        integer(int32), parameter :: n_genes = 40, n_axes = 2
        integer(int32), parameter :: M = 20000
        integer(int32) :: nrc(n_axes), nrt(n_axes)
        real(real64) :: mc(n_genes, n_axes), mt(n_genes, n_axes)
        real(real64), allocatable :: pc(:), pt(:)
        real(real64) :: b1(n_axes, n_genes), b2(n_axes, n_genes)
        real(real64) :: p1(n_genes), d1(n_genes), ds1(n_genes), s1(n_genes)
        real(real64) :: p2(n_genes), d2(n_genes), ds2(n_genes), s2(n_genes)
        integer(int32) :: m1(n_genes), nd1(n_genes), m2(n_genes), nd2(n_genes)
        integer(int32) :: nbc(n_axes, n_genes), nbt(n_axes, n_genes)
        real(real64) :: vn(n_axes, n_genes), cc(n_axes), mm(n_axes), dh(n_axes)
        integer(int32) :: nwp, ierr, g, n_cmp
        real(real64) :: sum_abs, max_abs

        nrc = [4, 4]
        nrt = [4, 4]
        allocate(pc(packed_len(n_genes, n_axes, nrc)))
        allocate(pt(packed_len(n_genes, n_axes, nrt)))
        call make_group(n_genes, n_axes, nrc, 0.0_real64, 0.9_real64, 0.30_real64, mc, pc)
        call make_group(n_genes, n_axes, nrt, 1.3_real64, 0.0_real64, 0.30_real64, mt, pt)

        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, M, 0, ENUM_OFF, &
                    b1, p1, d1, ds1, s1, m1, nd1, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "sys-vs-rng: systematic ierr")

        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_RNG, M, 0, ENUM_OFF, &
                    b2, p2, d2, ds2, s2, m2, nd2, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "sys-vs-rng: rng ierr")

        sum_abs = 0.0_real64
        max_abs = 0.0_real64
        n_cmp = 0
        do g = 1, n_genes
            if (p1(g) < 0.0_real64 .or. p2(g) < 0.0_real64) cycle
            n_cmp = n_cmp + 1
            ! The observed side is deterministic and must be identical.
            call assert_equal_real(d1(g), d2(g), EXACT, &
                                   "sys-vs-rng: D differs between samplers at gene "//str(g))
            call assert_equal_real(ds1(g), ds2(g), EXACT, &
                                   "sys-vs-rng: D_std differs between samplers at gene "//str(g))
            sum_abs = sum_abs + abs(p1(g) - p2(g))
            max_abs = max(max_abs, abs(p1(g) - p2(g)))
        end do
        call assert_true(n_cmp > 0, "sys-vs-rng: nothing compared")
        ! Monte-Carlo SE at M = 20,000 is <= 0.0035; 0.02 mean / 0.08 max is a loose
        ! but non-vacuous band around that.
        call assert_true(sum_abs / real(n_cmp, real64) < 0.02_real64, &
                         "sys-vs-rng: mean |p_sys - p_rng| exceeds Monte-Carlo tolerance")
        call assert_true(max_abs < 0.08_real64, &
                         "sys-vs-rng: max |p_sys - p_rng| exceeds Monte-Carlo tolerance")

        deallocate(pc, pt)
    end subroutine test_systematic_vs_rng

    ! =========================================================================
    ! 8. beta_centre = 1 subtracts the per-axis median
    ! =========================================================================
    subroutine test_beta_centre_median()
        integer(int32), parameter :: n_genes = 40, n_axes = 2
        real(real64), parameter :: OFFSET = 3.75_real64
        integer(int32) :: nrc(n_axes), nrt(n_axes)
        real(real64) :: mc(n_genes, n_axes), mt(n_genes, n_axes)
        real(real64), allocatable :: pc(:), pt(:)
        real(real64) :: beta_plain(n_axes, n_genes), beta_c(n_axes, n_genes), beta_off(n_axes, n_genes)
        real(real64) :: p0(n_genes), d0(n_genes), ds0(n_genes), s0(n_genes)
        real(real64) :: pcv(n_genes), dc(n_genes), dsc(n_genes), sc(n_genes)
        real(real64) :: po(n_genes), dof(n_genes), dso(n_genes), so(n_genes)
        integer(int32) :: m0(n_genes), n0(n_genes), mcm(n_genes), nc(n_genes), mo(n_genes), no(n_genes)
        integer(int32) :: nbc(n_axes, n_genes), nbt(n_axes, n_genes)
        real(real64) :: vn(n_axes, n_genes)
        real(real64) :: cor0(n_axes), mad0(n_axes), delta0(n_axes)
        real(real64) :: corc(n_axes), madc(n_axes), deltac(n_axes)
        real(real64) :: coro(n_axes), mado(n_axes), deltao(n_axes)
        integer(int32) :: nwp, ierr, g, a

        nrc = [4, 3]
        nrt = [3, 4]
        allocate(pc(packed_len(n_genes, n_axes, nrc)))
        allocate(pt(packed_len(n_genes, n_axes, nrt)))
        call make_group(n_genes, n_axes, nrc, 0.6_real64, 1.4_real64, 0.30_real64, mc, pc)
        call make_group(n_genes, n_axes, nrt, 2.2_real64, 0.0_real64, 0.30_real64, mt, pt)

        ! Run 0: no centring -- gives the uncentred beta and the shift estimate.
        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, 500, 0, ENUM_OFF, &
                    beta_plain, p0, d0, ds0, s0, m0, n0, nbc, nbt, vn, cor0, mad0, delta0, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "centre: uncentred ierr")
        call assert_true(nwp > 0, "centre: no gene received a p-value")

        ! delta_hat is the per-axis median of beta across genes, reported even with
        ! centring off.
        do a = 1, n_axes
            call assert_equal_real(delta0(a), median_ref(beta_plain(a, :), n_genes), TOL, &
                                   "centre: delta_hat is not the per-axis median, axis "//str(a))
        end do

        ! Run C: supply the uncentred beta and centre it.
        beta_c = beta_plain
        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_SUPPLIED, BETA_CENTRE_MEDIAN, &
                    SAMPLING_SYSTEMATIC, 500, 0, ENUM_OFF, &
                    beta_c, pcv, dc, dsc, sc, mcm, nc, nbc, nbt, vn, corc, madc, deltac, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "centre: centred ierr")
        do a = 1, n_axes
            call assert_equal_real(deltac(a), delta0(a), TOL, &
                                   "centre: delta_hat changed, axis "//str(a))
            do g = 1, n_genes
                call assert_equal_real(beta_c(a, g), beta_plain(a, g) - delta0(a), EXACT, &
                                       "centre: beta not shifted by the median, axis "//str(a) &
                                       //" gene "//str(g))
            end do
        end do

        ! Run O: the SAME beta with a constant offset added to axis 1 only. The
        ! median moves by exactly the offset, so the centred statistic is unchanged --
        ! which is the property that makes centring a composition correction rather
        ! than a rescaling.
        beta_off = beta_plain
        beta_off(1, :) = beta_off(1, :) + OFFSET
        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_SUPPLIED, BETA_CENTRE_MEDIAN, &
                    SAMPLING_SYSTEMATIC, 500, 0, ENUM_OFF, &
                    beta_off, po, dof, dso, so, mo, no, nbc, nbt, vn, coro, mado, deltao, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "centre: offset ierr")

        call assert_equal_real(deltao(1) - deltac(1), OFFSET, TOL, &
                               "centre: delta_hat did not recover the injected offset")
        call assert_equal_real(deltao(2), deltac(2), TOL, &
                               "centre: axis 2 delta_hat moved when only axis 1 was offset")

        do g = 1, n_genes
            call assert_equal_real(beta_off(1, g), beta_c(1, g), TOL, &
                                   "centre: centred beta changed under the offset, gene "//str(g))
            if (dc(g) >= 0.0_real64) then
                call assert_equal_real(dof(g), dc(g), TOL, &
                                       "centre: centred D changed under the offset, gene "//str(g))
                call assert_equal_real(dso(g), dsc(g), TOL, &
                                       "centre: centred D_std changed under the offset, gene "//str(g))
            end if
        end do

        deallocate(pc, pt)
    end subroutine test_beta_centre_median

    !> Independent median, by insertion sort of a copy -- deliberately not the
    !| module's own routine, so the two can disagree.
    function median_ref(values, n) result(med)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: values(n)
        real(real64) :: med
        real(real64) :: v(n), tmp
        integer(int32) :: i, j

        v = values
        do i = 2, n
            tmp = v(i)
            j = i - 1
            do while (j >= 1)
                if (v(j) <= tmp) exit
                v(j + 1) = v(j)
                j = j - 1
            end do
            v(j + 1) = tmp
        end do
        if (mod(n, 2) == 1) then
            med = v((n + 1) / 2)
        else
            med = 0.5_real64 * (v(n / 2) + v(n / 2 + 1))
        end if
    end function median_ref

    ! =========================================================================
    ! 9. degenerate input
    ! =========================================================================
    subroutine test_degenerate_input()
        integer(int32), parameter :: n_genes = 40, n_axes = 2
        integer(int32) :: nrc(n_axes), nrt(n_axes)
        real(real64) :: mc(n_genes, n_axes), mt(n_genes, n_axes)
        real(real64), allocatable :: pc(:), pt(:)
        real(real64) :: b(n_axes, n_genes)
        real(real64) :: p(n_genes), d_obs(n_genes), d_std(n_genes), s(n_genes)
        integer(int32) :: meth(n_genes), nd(n_genes)
        integer(int32) :: nbc(n_axes, n_genes), nbt(n_axes, n_genes)
        real(real64) :: vn(n_axes, n_genes), cc(n_axes), mm(n_axes), dh(n_axes)
        integer(int32) :: nwp, ierr, g, off, r, base

        ! --- thin axis: too few genes on axis 2 to ever reach 10 residuals -------
        block
            integer(int32), parameter :: n_thin = 4
            integer(int32) :: trc(2), trt(2)
            integer(int32) :: tm(2, n_thin), tn(2, n_thin)
            real(real64) :: tmc(n_thin, 2), tmt(n_thin, 2)
            real(real64), allocatable :: tpc(:), tpt(:)
            real(real64) :: tb(2, n_thin)
            real(real64) :: tp(n_thin), td(n_thin), tds(n_thin), ts(n_thin)
            integer(int32) :: tmeth(n_thin), tnd(n_thin)
            real(real64) :: tvn(2, n_thin), tcc(2), tmm(2), tdh(2)
            integer(int32) :: tnwp, terr

            ! Axis 2 holds 2 replicates over 4 genes -> at most 8 residuals, forever
            ! under the MIN_POOL_RESIDUALS gate.
            trc = [4, 2]
            trt = [4, 2]
            allocate(tpc(packed_len(n_thin, 2, trc)))
            allocate(tpt(packed_len(n_thin, 2, trt)))
            call make_group(n_thin, 2, trc, 0.0_real64, 1.0_real64, 0.3_real64, tmc, tpc)
            call make_group(n_thin, 2, trt, 1.0_real64, 0.0_real64, 0.3_real64, tmt, tpt)

            call run_md(n_thin, 2, trc, trt, tmc, tpc, tmt, tpt, &
                        0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                        SAMPLING_SYSTEMATIC, 200, 0, ENUM_OFF, &
                        tb, tp, td, tds, ts, tmeth, tnd, tm, tn, tvn, tcc, tmm, tdh, tnwp, terr)

            call assert_equal_int(terr, ERR_OK, "degenerate: a thin axis is data, not an error")
            call assert_equal_int(tnwp, 0, "degenerate: thin axis must drop every gene")
            do g = 1, n_thin
                call assert_equal_real(tp(g), -1.0_real64, EXACT, &
                                       "degenerate: thin-axis gene "//str(g)//" must have no p-value")
                call assert_equal_int(tmeth(g), -1, &
                                      "degenerate: untested gene must report method_used = -1")
                ! The failing axis stays identifiable: its pool size is reported as
                ! gathered (0 here, since the gather gives up below 10).
                call assert_true(tm(2, g) < 10, &
                                 "degenerate: the thin axis must be visible at gene "//str(g))
            end do
            deallocate(tpc, tpt)
        end block

        ! --- zero-variance axis, and therefore D = 0 ----------------------------
        nrc = [4, 4]
        nrt = [4, 4]
        allocate(pc(packed_len(n_genes, n_axes, nrc)))
        allocate(pt(packed_len(n_genes, n_axes, nrt)))
        call make_group(n_genes, n_axes, nrc, 0.0_real64, 0.0_real64, 0.3_real64, mc, pc)
        call make_group(n_genes, n_axes, nrt, 1.0_real64, 0.0_real64, 0.3_real64, mt, pt)

        ! Flatten axis 2 on BOTH sides to the identical constant per gene: zero
        ! residual variance, zero beta. Var_null(2) is then 0, so axis 2 must get
        ! weight 0 and drop out of D_std rather than producing an infinity.
        off = nrc(1) * n_genes
        do g = 1, n_genes
            base = off + (g - 1) * nrc(2)
            do r = 1, nrc(2)
                pc(base + r) = 20.0_real64
            end do
            mc(g, 2) = 20.0_real64
            base = off + (g - 1) * nrt(2)
            do r = 1, nrt(2)
                pt(base + r) = 20.0_real64
            end do
            mt(g, 2) = 20.0_real64
        end do

        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    0, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, 500, 0, ENUM_OFF, &
                    b, p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "degenerate: zero-variance axis ierr")
        call assert_true(nwp > 0, "degenerate: zero-variance axis dropped everything")
        do g = 1, n_genes
            if (p(g) < 0.0_real64) cycle
            call assert_equal_real(b(2, g), 0.0_real64, EXACT, &
                                   "degenerate: flat axis must give beta = 0 at gene "//str(g))
            call assert_equal_real(vn(2, g), 0.0_real64, EXACT, &
                                   "degenerate: flat axis must give Var_null = 0 at gene "//str(g))
            call assert_equal_real(d_obs(g), abs(b(1, g)), TOL, &
                                   "degenerate: flat axis must not enter D at gene "//str(g))
            call assert_true(d_std(g) == d_std(g), &
                             "degenerate: D_std must not be NaN at gene "//str(g))
            call assert_in_range_real(p(g), 1.0_real64 / 501.0_real64, 1.0_real64, &
                                      "degenerate: p out of range at gene "//str(g))
        end do

        ! --- D = 0: a gene with no difference on any axis is maximally unremarkable
        block
            real(real64) :: b0(n_axes, n_genes)
            real(real64) :: pp(n_genes), dd(n_genes), dds(n_genes), ss(n_genes)
            integer(int32) :: mth(n_genes), ndr(n_genes), nwp2, err2
            real(real64) :: v2(n_axes, n_genes), c2(n_axes), m2(n_axes), d2(n_axes)

            b0 = b
            b0(:, 3) = 0.0_real64
            call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                        0, BETA_MODE_SUPPLIED, BETA_CENTRE_NONE, &
                        SAMPLING_SYSTEMATIC, 500, 0, ENUM_OFF, &
                        b0, pp, dd, dds, ss, mth, ndr, nbc, nbt, v2, c2, m2, d2, nwp2, err2)
            call assert_equal_int(err2, ERR_OK, "degenerate: D = 0 ierr")
            call assert_equal_real(dd(3), 0.0_real64, EXACT, "degenerate: D must be 0 at gene 3")
            call assert_equal_real(dds(3), 0.0_real64, EXACT, "degenerate: D_std must be 0 at gene 3")
            ! Every null draw satisfies D_std_null >= 0, so p is pinned at 1.
            call assert_equal_real(pp(3), 1.0_real64, TOL, &
                                   "degenerate: D = 0 must give p = 1")
        end block

        ! --- non-finite supplied beta: that gene is dropped, the run continues ---
        block
            real(real64) :: bn(n_axes, n_genes)
            real(real64) :: pp(n_genes), dd(n_genes), dds(n_genes), ss(n_genes)
            integer(int32) :: mth(n_genes), ndr(n_genes), nwp2, err2
            real(real64) :: v2(n_axes, n_genes), c2(n_axes), m2(n_axes), d2(n_axes)

            bn = b
            bn(1, 5) = ieee_value(1.0_real64, ieee_quiet_nan)
            call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                        0, BETA_MODE_SUPPLIED, BETA_CENTRE_NONE, &
                        SAMPLING_SYSTEMATIC, 500, 0, ENUM_OFF, &
                        bn, pp, dd, dds, ss, mth, ndr, nbc, nbt, v2, c2, m2, d2, nwp2, err2)
            call assert_equal_int(err2, ERR_OK, "degenerate: NaN beta is data, not an error")
            call assert_equal_real(pp(5), -1.0_real64, EXACT, &
                                   "degenerate: a NaN beta must drop exactly its own gene")
            call assert_true(nwp2 > 0, "degenerate: a single NaN beta must not drop the run")
        end block

        ! --- argument validation -------------------------------------------------
        call assert_equal_int(bad_run(n_genes, n_axes, [1, 4], nrt, mc, pc, mt, pt, &
                                      BETA_MODE_INTERNAL, BETA_CENTRE_NONE, NULL_METHOD_POOLED, &
                                      SAMPLING_SYSTEMATIC, 500, ENUM_OFF), &
                              ERR_INVALID_INPUT, "degenerate: n_rep < 2 must be ERR_INVALID_INPUT")
        call assert_equal_int(bad_run(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                                      7, BETA_CENTRE_NONE, NULL_METHOD_POOLED, &
                                      SAMPLING_SYSTEMATIC, 500, ENUM_OFF), &
                              ERR_INVALID_INPUT, "degenerate: bad beta_mode must be ERR_INVALID_INPUT")
        call assert_equal_int(bad_run(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                                      BETA_MODE_INTERNAL, 7, NULL_METHOD_POOLED, &
                                      SAMPLING_SYSTEMATIC, 500, ENUM_OFF), &
                              ERR_INVALID_INPUT, "degenerate: bad beta_centre must be ERR_INVALID_INPUT")
        ! Gene-blocking is not a distinct null here, so it is rejected rather than
        ! silently ignored.
        call assert_equal_int(bad_run(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                                      BETA_MODE_INTERNAL, BETA_CENTRE_NONE, 1, &
                                      SAMPLING_SYSTEMATIC, 500, ENUM_OFF), &
                              ERR_INVALID_INPUT, "degenerate: null_method /= 0 must be ERR_INVALID_INPUT")
        call assert_equal_int(bad_run(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                                      BETA_MODE_INTERNAL, BETA_CENTRE_NONE, NULL_METHOD_POOLED, &
                                      9, 500, ENUM_OFF), &
                              ERR_INVALID_INPUT, "degenerate: bad sampling_mode must be ERR_INVALID_INPUT")
        call assert_equal_int(bad_run(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                                      BETA_MODE_INTERNAL, BETA_CENTRE_NONE, NULL_METHOD_POOLED, &
                                      SAMPLING_SYSTEMATIC, 0, ENUM_OFF), &
                              ERR_EMPTY_INPUT, "degenerate: n_draws_max = 0 must be ERR_EMPTY_INPUT")
        call assert_equal_int(bad_run(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                                      BETA_MODE_INTERNAL, BETA_CENTRE_NONE, NULL_METHOD_POOLED, &
                                      SAMPLING_SYSTEMATIC, 500, -1_int64), &
                              ERR_INVALID_INPUT, "degenerate: negative enum cap must be ERR_INVALID_INPUT")

        deallocate(pc, pt)
    end subroutine test_degenerate_input

    !> Run the pipeline purely for its error code, discarding every output.
    function bad_run(n_genes, n_axes, n_rep_case, n_rep_control, &
                     means_case, packed_case, means_control, packed_control, &
                     beta_mode, beta_centre, null_method, sampling_mode, &
                     n_draws_max, enum_max_product) result(ierr)
        integer(int32), intent(in) :: n_genes, n_axes
        integer(int32), intent(in) :: n_rep_case(n_axes), n_rep_control(n_axes)
        real(real64), intent(in) :: means_case(n_genes, n_axes), means_control(n_genes, n_axes)
        real(real64), intent(in) :: packed_case(:), packed_control(:)
        integer(int32), intent(in) :: beta_mode, beta_centre, null_method, sampling_mode, n_draws_max
        integer(int64), intent(in) :: enum_max_product
        integer(int32) :: ierr

        real(real64) :: beta_obs(n_axes, n_genes)
        real(real64) :: p(n_genes), d_obs(n_genes), d_std(n_genes), s(n_genes)
        integer(int32) :: meth(n_genes), nd(n_genes), valid(n_genes)
        integer(int32) :: nbc(n_axes, n_genes), nbt(n_axes, n_genes)
        real(real64) :: vn(n_axes, n_genes), cc(n_axes), mm(n_axes), dh(n_axes)
        integer(int32) :: nwp

        beta_obs = 0.0_real64
        valid = 1
        call compute_noise_pvalue_pipeline_md( &
            means_case, packed_case, n_rep_case, &
            means_control, packed_control, n_rep_control, &
            beta_obs, valid, beta_mode, beta_centre, &
            n_genes, n_axes, 0, K_START, K_STEP, K_MAX, TAU, 0.0_real64, &
            null_method, sampling_mode, n_draws_max, 0, enum_max_product, SEED, MAX_POOL, &
            p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
    end function bad_run

    ! =========================================================================
    ! 10. beta_mode = 1 with an off-scale beta
    ! =========================================================================
    !| The concrete failure mode from the plan's Q4: edgeR's `coefficients` are on
    !| the NATURAL log scale, a factor of ln 2 against log2-space residuals. The
    !| diagnostic must separate that from shrinkage -- correlation stays ~1 while the
    !| MAD blows up -- and the run must still complete rather than refuse.
    subroutine test_beta_mode_supplied_offscale()
        integer(int32), parameter :: n_genes = 40, n_axes = 2
        integer(int32) :: nrc(n_axes), nrt(n_axes)
        real(real64) :: mc(n_genes, n_axes), mt(n_genes, n_axes)
        real(real64), allocatable :: pc(:), pt(:)
        real(real64) :: beta_int(n_axes, n_genes), beta_bad(n_axes, n_genes)
        real(real64) :: p(n_genes), d_obs(n_genes), d_std(n_genes), s(n_genes)
        integer(int32) :: meth(n_genes), nd(n_genes)
        integer(int32) :: nbc(n_axes, n_genes), nbt(n_axes, n_genes)
        real(real64) :: vn(n_axes, n_genes), cc(n_axes), mm(n_axes), dh(n_axes)
        integer(int32) :: nwp, ierr, a

        nrc = [4, 4]
        nrt = [4, 4]
        allocate(pc(packed_len(n_genes, n_axes, nrc)))
        allocate(pt(packed_len(n_genes, n_axes, nrt)))
        call make_group(n_genes, n_axes, nrc, 0.0_real64, 1.5_real64, 0.30_real64, mc, pc)
        call make_group(n_genes, n_axes, nrt, 1.0_real64, 0.0_real64, 0.30_real64, mt, pt)

        ! Internal mode: the check is exact by construction.
        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    1, BETA_MODE_INTERNAL, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, 500, 0, ENUM_OFF, &
                    beta_int, p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)
        call assert_equal_int(ierr, ERR_OK, "offscale: internal ierr")
        do a = 1, n_axes
            call assert_equal_real(cc(a), 1.0_real64, EXACT, &
                                   "offscale: internal mode must report cor = 1, axis "//str(a))
            call assert_equal_real(mm(a), 0.0_real64, EXACT, &
                                   "offscale: internal mode must report MAD = 0, axis "//str(a))
        end do

        ! Supplied, on the wrong log base: a pure rescaling by ln 2.
        beta_bad = beta_int * log(2.0_real64)
        call run_md(n_genes, n_axes, nrc, nrt, mc, pc, mt, pt, &
                    1, BETA_MODE_SUPPLIED, BETA_CENTRE_NONE, &
                    SAMPLING_SYSTEMATIC, 500, 0, ENUM_OFF, &
                    beta_bad, p, d_obs, d_std, s, meth, nd, nbc, nbt, vn, cc, mm, dh, nwp, ierr)

        call assert_equal_int(ierr, ERR_OK, "offscale: supplied beta must not be refused")
        call assert_true(nwp > 0, "offscale: the run must complete and produce p-values")
        do a = 1, n_axes
            ! A pure positive rescaling leaves the correlation at 1 (up to rounding)
            ! -- the units signature.
            call assert_true(cc(a) > 0.999_real64, &
                             "offscale: correlation should stay ~1 under a rescaling, axis "//str(a))
            call assert_true(mm(a) > 1.0e-6_real64, &
                             "offscale: MAD should be large under a rescaling, axis "//str(a))
        end do

        deallocate(pc, pt)
    end subroutine test_beta_mode_supplied_offscale

end module mod_test_noise_model_md
