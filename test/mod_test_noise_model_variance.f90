!> Unit tests for the variance-aware noise model (`noise_model_variance`,
!| implementation notes A and F).
!|
!|   1. Trigamma against closed forms and note A's tau rule `sqrt(psi1(d/2))`.
!|   2. Hand-built band gather: which genes enter, in which order, scan count,
!|      `scan_max`, short pools, non-finite candidates.
!|   3. Regression (note F test 1): band, raw, no scaling reproduces the pre-note-F
!|      outputs bit-identically (fixture `test/test_files/variance_band_regression.txt`,
!|      recorded with the old module on 28/10/2026).
!|   4. Parity: mean mode == `noise_model_exact` (adaptive kNN); band with an
!|      all-accepting band == exact at fixed k; cross-fit / same-group scaling with
!|      lambda = 1 and a constant target == exact; adaptive with neutral settings ==
!|      band (note F §5).
!|   5. Adaptive gather, hand-built: every stop reason.
!|   6. Shrinkage: the contrast `Vhat` is unbiased (normal, skewed, t10; n = 3, 5, 10;
!|      unbalanced); the v1 estimator equals `s2^2/4` at n = 3 (documentation test);
!|      lambda near 1 without heterogeneity, clearly below with it; relative ==
!|      absolute weights for a constant target.
!|   7. Calibration regression (n = 3, heterogeneous variances, null): cross-fitted
!|      JS tail `p <= 0.001` at most 2x nominal, same-group clearly above.
!|   8. Scale invariance (linear branch), the floor path, validation.
!|   9. 2-D ellipse gather: equals a brute-force search for the k nearest in the
!|      joint distance; its stop reasons; a pipeline run.
module mod_test_noise_model_variance
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
    use noise_model_exact, only: sorted_data_t, prepare_sorted_data, &
                                 pipeline_exact => compute_noise_pvalue_pipeline, BETA_CENTRE_NONE
    use noise_model_variance
    use tox_errors
    use test_suite, only: test_case
    implicit none
    public

    real(real64), parameter :: TOL = 1d-12
    real(real64), parameter :: EXACT = 0.0_real64
    real(real64), parameter :: PI = 3.14159265358979323846_real64

    !> Pipeline settings, defaulting to the benchmark's (mean kNN k20_50 / tau 0.1).
    type :: vcfg
        integer(int32) :: mode = VARIANCE_MODE_MEAN, match = MATCH_RAW, scaling = SCALING_NONE
        integer(int32) :: target = SHRINK_TARGET_MEDIAN, lscale = LAMBDA_RELATIVE, norm = 1
        integer(int32) :: k_start = 20, k_step = 1, k_max = 50, k = 50
        integer(int32) :: b = 5, k_min = 20, w = 2, robust = 1, srobust = 0
        real(real64) :: tau = 0.1_real64, c_tau = 1.5_real64, c_js = 1.0_real64
        real(real64) :: r_max = 0.5_real64, tau_stab = 0.1_real64, a_min = 0.2_real64
    end type vcfg

    !> Every pipeline output.
    type :: vout
        real(real64), allocatable :: p(:), obs(:), l2(:), u(:), s2s(:), radc(:), radt(:)
        integer(int32), allocatable :: nbc(:), nbt(:), nsc(:), nst(:), nac(:), nat(:), stc(:), stt(:)
        real(real64) :: dh, tau_var, scatter, diag(N_SHRINK_DIAG)
        integer(int32) :: nwp, nf, ierr
    end type vout

contains

    function get_all_tests_noise_model_variance() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate(all_tests(15))
        all_tests(14) = test_case("test_ellipse_brute_force",         test_ellipse_brute_force)
        all_tests(15) = test_case("test_ellipse_stops_and_pipeline",  test_ellipse_stops_and_pipeline)
        all_tests(1)  = test_case("test_trigamma_and_tau_rule",       test_trigamma_and_tau_rule)
        all_tests(2)  = test_case("test_band_gather_hand_built",      test_band_gather_hand_built)
        all_tests(3)  = test_case("test_band_regression_fixture",     test_band_regression_fixture)
        all_tests(4)  = test_case("test_mean_parity_with_exact",      test_mean_parity_with_exact)
        all_tests(5)  = test_case("test_band_parity_with_exact",      test_band_parity_with_exact)
        all_tests(6)  = test_case("test_scaling_parity_lambda_one",   test_scaling_parity_lambda_one)
        all_tests(7)  = test_case("test_adaptive_parity_with_band",   test_adaptive_parity_with_band)
        all_tests(8)  = test_case("test_adaptive_stop_reasons",       test_adaptive_stop_reasons)
        all_tests(9)  = test_case("test_vhat_unbiased",               test_vhat_unbiased)
        all_tests(10) = test_case("test_lambda_behaviour",            test_lambda_behaviour)
        all_tests(11) = test_case("test_xfit_calibration_regression", test_xfit_calibration_regression)
        all_tests(12) = test_case("test_scale_invariance_and_floor",  test_scale_invariance_and_floor)
        all_tests(13) = test_case("test_variance_validation",         test_variance_validation)
    end function get_all_tests_noise_model_variance

    function str(i) result(s)
        integer(int32), intent(in) :: i
        character(len=20) :: buffer
        character(len=:), allocatable :: s
        write(buffer, '(I0)') i
        s = trim(buffer)
    end function str

    !> Heterogeneous fixture: `n_genes` genes over a wide mean range, per-gene noise
    !| amplitude varying smoothly AND gene-to-gene, deterministic.
    subroutine make_hetero(n_genes, n_rep, phase, shift, reps)
        integer(int32), intent(in) :: n_genes, n_rep
        real(real64), intent(in) :: phase, shift
        real(real64), intent(out) :: reps(n_rep, n_genes)
        integer(int32) :: g, r
        real(real64) :: mu, amp
        do g = 1, n_genes
            mu = 10.0_real64 + 990.0_real64 * real(g - 1, real64) / real(n_genes - 1, real64) + shift
            amp = mu * (0.05_real64 + 0.25_real64 * (0.5_real64 + 0.5_real64 * sin(2.7_real64 * real(g, real64))))
            do r = 1, n_rep
                reps(r, g) = mu + amp * sin(phase + 1.3_real64 * real(g, real64) + 2.1_real64 * real(r, real64))
            end do
        end do
    end subroutine make_hetero

    !> Standard normal draws (Box-Muller on `random_number`).
    subroutine rnorm(x)
        real(real64), intent(out) :: x(:)
        real(real64) :: u1(size(x)), u2(size(x))
        call random_number(u1)
        call random_number(u2)
        x = sqrt(-2.0_real64 * log(1.0_real64 - u1)) * cos(2.0_real64 * PI * u2)
    end subroutine rnorm

    !> Run the variance pipeline with `cfg`; allocates every output in `o`.
    subroutine run_v(cfg, G, nc, nt, rc, rt, o, lambda_fixed, target_fixed)
        type(vcfg), intent(in) :: cfg
        integer(int32), intent(in) :: G, nc, nt
        real(real64), intent(in) :: rc(nc, G), rt(nt, G)
        type(vout), intent(out) :: o
        real(real64), intent(in), optional :: lambda_fixed, target_fixed
        integer(int32) :: valid(G)
        valid = 1
        allocate(o%p(G), o%obs(G), o%l2(G), o%u(G), o%s2s(G), o%radc(G), o%radt(G))
        allocate(o%nbc(G), o%nbt(G), o%nsc(G), o%nst(G), o%nac(G), o%nat(G), o%stc(G), o%stt(G))
        call compute_noise_pvalue_pipeline_variance(rc, G, nc, rt, G, nt, valid, BETA_CENTRE_NONE, G, cfg%norm, &
            cfg%mode, cfg%match, cfg%scaling, cfg%target, cfg%lscale, &
            cfg%k_start, cfg%k_step, cfg%k_max, cfg%tau, cfg%k, cfg%c_tau, cfg%c_js, &
            cfg%b, cfg%k_min, cfg%r_max, cfg%tau_stab, cfg%w, cfg%a_min, &
            0.3_real64, 1, cfg%robust, cfg%srobust, 5000, &
            o%p, o%obs, o%l2, o%dh, o%nwp, o%nbc, o%nbt, o%nsc, o%nst, o%nac, o%nat, &
            o%radc, o%radt, o%stc, o%stt, o%u, o%s2s, o%tau_var, o%scatter, o%diag, o%nf, o%ierr, &
            lambda_fixed, target_fixed)
    end subroutine run_v

    !> The plain exact model.
    subroutine run_exact(G, nc, nt, rc, rt, k_start, k_step, k_max, tau, norm, p, nbc, nbt, nwp, ierr)
        integer(int32), intent(in) :: G, nc, nt, k_start, k_step, k_max, norm
        real(real64), intent(in) :: rc(nc, G), rt(nt, G), tau
        real(real64), intent(out) :: p(G)
        integer(int32), intent(out) :: nbc(G), nbt(G), nwp, ierr
        integer(int32) :: valid(G), nb3(G)
        real(real64) :: obs(G), l2(G), dh
        valid = 1
        call pipeline_exact(rc, G, nc, rt, G, nt, valid, BETA_CENTRE_NONE, G, norm, k_start, k_step, k_max, tau, &
                            p, obs, l2, dh, nwp, 5000, nbc, nbt, nb3, ierr)
    end subroutine run_exact

    ! =========================================================================
    ! 1. trigamma and the tau rule
    ! =========================================================================
    subroutine test_trigamma_and_tau_rule()
        integer(int32), parameter :: G = 300
        integer(int32) :: n_list(3), i, n, nf, ierr
        real(real64) :: expect(3)
        real(real64), allocatable :: rc(:, :), rt(:, :)
        real(real64) :: u(G), tau, sc

        call assert_equal_real(trigamma_helper(1.0_real64), PI**2 / 6.0_real64, 1d-10, "trigamma(1)")
        call assert_equal_real(trigamma_helper(0.5_real64), PI**2 / 2.0_real64, 1d-10, "trigamma(1/2)")
        call assert_equal_real(trigamma_helper(2.0_real64), PI**2 / 6.0_real64 - 1.0_real64, 1d-10, "trigamma(2)")

        n_list = [3, 5, 10]
        expect = [0.80_real64, 0.53_real64, 0.34_real64]
        do i = 1, 3
            n = n_list(i)
            allocate(rc(n, G), rt(n, G))
            call make_hetero(G, n, 0.0_real64, 0.0_real64, rc)
            call make_hetero(G, n, 1.0_real64, 0.0_real64, rt)
            call compute_variance_match(rc, n, rt, n, G, 1.0_real64, 0.3_real64, 1, 1, u, tau, sc, nf, ierr)
            call assert_equal_int(ierr, ERR_OK, "tau rule: ierr at n = "//str(n))
            call assert_equal_real(tau, sqrt(trigamma_helper(real(n - 1, real64))), TOL, &
                                   "tau rule: tau is c*sqrt(trigamma(d/2)) at n = "//str(n))
            call assert_equal_real(tau, expect(i), 0.006_real64, "tau rule: note-A value at n = "//str(n))
            deallocate(rc, rt)
        end do
    end subroutine test_trigamma_and_tau_rule

    ! =========================================================================
    ! 2. hand-built band gather
    ! =========================================================================
    !| 12 genes with means 1..12 (3 replicates, residuals -1/0/+1 scaled by the
    !| gene index so pools are traceable). Target = gene 6 (mean 6, u = 0). Visit
    !| order from slot 6 under `choose_index` (ties go left): 6, 5, 7, 4, 8, 3, 9,
    !| 2, 10, 1, 11, 12.
    subroutine hand_built(sd)
        type(sorted_data_t), intent(out) :: sd
        integer(int32), parameter :: G = 12, NR = 3
        real(real64) :: means(G), reps(NR, G)
        integer(int32) :: ig, ierr
        do ig = 1, G
            means(ig) = real(ig, real64)
            reps(:, ig) = means(ig) + real(ig, real64) * [-1.0_real64, 0.0_real64, 1.0_real64]
        end do
        call prepare_sorted_data(means, reps, NR, G, 0, sd, ierr)
        call assert_equal_int(ierr, ERR_OK, "hand-built: prepare_sorted_data")
    end subroutine hand_built

    subroutine test_band_gather_hand_built()
        integer(int32), parameter :: G = 12
        real(real64) :: u(G), pool(100), nan
        type(sorted_data_t) :: sd
        integer(int32) :: ig, n_pooled, n_gp, n_scan
        integer(int32), parameter :: expect_genes(4) = [6, 7, 4, 9]
        real(real64) :: expect_pool(12)

        nan = ieee_value(1.0_real64, ieee_quiet_nan)
        call hand_built(sd)
        u = 5.0_real64
        u(6) = 0.0_real64; u(5) = 5.0_real64; u(7) = 0.2_real64; u(4) = -0.4_real64
        u(8) = nan; u(3) = -0.6_real64; u(9) = 0.5_real64

        call gather_residuals_band_helper(6.0_real64, 0.0_real64, sd, u, 4, 0.5_real64, 1000, &
                                          pool, n_pooled, n_gp, n_scan, 100)
        call assert_equal_int(n_gp, 4, "hand-built: 4 genes accepted")
        call assert_equal_int(n_scan, 7, "hand-built: visits 6,5,7,4,8,3,9")
        call assert_equal_int(n_pooled, 12, "hand-built: 4 genes x 3 residuals")
        do ig = 1, 4
            expect_pool(3 * ig - 2 : 3 * ig) = real(expect_genes(ig), real64) * sqrt(1.5_real64) &
                                             * [-1.0_real64, 0.0_real64, 1.0_real64]
        end do
        call assert_equal_array_real(pool(1:12), expect_pool, 12, TOL, &
                                     "hand-built: pool = genes 6, 7, 4, 9 in visit order")

        call gather_residuals_band_helper(6.0_real64, 0.0_real64, sd, u, 4, 0.5_real64, 3, &
                                          pool, n_pooled, n_gp, n_scan, 100)
        call assert_equal_int(n_scan, 3, "hand-built: scan_max = 3 stops after 3 slots")
        call assert_equal_int(n_gp, 2, "hand-built: 2 genes accepted within scan_max")

        call gather_residuals_band_helper(6.0_real64, 0.0_real64, sd, u, 50, 0.5_real64, 1000, &
                                          pool, n_pooled, n_gp, n_scan, 100)
        call assert_equal_int(n_gp, 4, "hand-built: short pool keeps the 4 eligible genes")
        call assert_equal_int(n_scan, G, "hand-built: short pool scans every slot")

        call gather_residuals_band_helper(6.0_real64, 0.0_real64, sd, u, 4, 0.5_real64, 1000, &
                                          pool, n_pooled, n_gp, n_scan, 5)
        call assert_equal_int(n_pooled, 5, "hand-built: max_pool_size caps the residuals")
        call assert_equal_int(n_gp, 2, "hand-built: stops once the pool is full")
    end subroutine test_band_gather_hand_built

    ! =========================================================================
    ! 3. regression against the pre-note-F band
    ! =========================================================================
    subroutine test_band_regression_fixture()
        integer(int32), parameter :: G = 300
        integer(int32) :: u_in, i, n, nwp_ref, ig, ios
        real(real64) :: tau_ref, p_ref(G), obs_ref(G), u_ref(G)
        integer(int32) :: nbc_ref(G), nbt_ref(G), nsc_ref(G), nst_ref(G)
        real(real64), allocatable :: rc(:, :), rt(:, :)
        type(vcfg) :: cfg
        type(vout) :: o

        cfg%mode = VARIANCE_MODE_BAND
        cfg%c_tau = 3.0_real64
        open(newunit=u_in, file="test/test_files/variance_band_regression.txt", status="old", action="read", &
             iostat=ios)
        call assert_equal_int(ios, 0, "regression: fixture file readable")
        if (ios /= 0) return
        do i = 1, 2
            read(u_in, *) n, nwp_ref, tau_ref
            do ig = 1, G
                read(u_in, *) p_ref(ig), obs_ref(ig), u_ref(ig), nbc_ref(ig), nbt_ref(ig), nsc_ref(ig), nst_ref(ig)
            end do
            allocate(rc(n, G), rt(n, G))
            call make_hetero(G, n, 0.0_real64, 30.0_real64, rc)
            call make_hetero(G, n, 1.7_real64, 0.0_real64, rt)
            call run_v(cfg, G, n, n, rc, rt, o)
            call assert_equal_int(o%ierr, ERR_OK, "regression: ierr at n = "//str(n))
            call assert_equal_int(o%nwp, nwp_ref, "regression: genes tested at n = "//str(n))
            call assert_equal_real(o%tau_var, tau_ref, EXACT, "regression: tau_var at n = "//str(n))
            call assert_equal_array_real(o%p, p_ref, G, EXACT, "regression: p-values bit-identical at n = "//str(n))
            call assert_equal_array_real(o%obs, obs_ref, G, EXACT, "regression: statistic at n = "//str(n))
            call assert_equal_array_real(o%u, u_ref, G, EXACT, "regression: u at n = "//str(n))
            call assert_equal_array_int(o%nbc, nbc_ref, G, "regression: case pools at n = "//str(n))
            call assert_equal_array_int(o%nbt, nbt_ref, G, "regression: control pools at n = "//str(n))
            call assert_equal_array_int(o%nsc, nsc_ref, G, "regression: case scans at n = "//str(n))
            call assert_equal_array_int(o%nst, nst_ref, G, "regression: control scans at n = "//str(n))
            deallocate(rc, rt)
        end do
        close(u_in)
    end subroutine test_band_regression_fixture

    ! =========================================================================
    ! 4. parities
    ! =========================================================================
    subroutine test_mean_parity_with_exact()
        integer(int32), parameter :: G = 400, N = 3
        real(real64) :: rc(N, G), rt(N, G), pe(G)
        integer(int32) :: nbce(G), nbte(G), nwe, ierre
        type(vcfg) :: cfg
        type(vout) :: o

        call make_hetero(G, N, 0.0_real64, 30.0_real64, rc)
        call make_hetero(G, N, 1.7_real64, 0.0_real64, rt)
        call run_v(cfg, G, N, N, rc, rt, o)
        call run_exact(G, N, N, rc, rt, 20, 1, 50, 0.1_real64, 1, pe, nbce, nbte, nwe, ierre)
        call assert_equal_int(o%ierr, ERR_OK, "mean parity: ierr")
        call assert_equal_int(ierre, ERR_OK, "mean parity: exact ierr")
        call assert_true(o%nwp > 0, "mean parity: p-values computed")
        call assert_equal_int(o%nwp, nwe, "mean parity: same tested genes")
        call assert_equal_array_real(o%p, pe, G, EXACT, "mean parity: p-values identical to the exact model")
        call assert_equal_array_int(o%nbc, nbce, G, "mean parity: case pools identical")
        call assert_equal_array_int(o%nbt, nbte, G, "mean parity: control pools identical")
        call assert_true(all(o%nsc == -1) .and. all(o%stc == -1), "mean parity: band diagnostics -1")
    end subroutine test_mean_parity_with_exact

    subroutine test_band_parity_with_exact()
        integer(int32), parameter :: G = 200, N = 4, K = 30
        real(real64) :: rc(N, G), rt(N, G), pe(G)
        integer(int32) :: nbce(G), nbte(G), nwe, ierre
        type(vcfg) :: cfg
        type(vout) :: o

        call make_hetero(G, N, 0.0_real64, 30.0_real64, rc)
        call make_hetero(G, N, 1.7_real64, 0.0_real64, rt)
        cfg%mode = VARIANCE_MODE_BAND
        cfg%k = K
        cfg%k_min = K
        cfg%c_tau = 1.0e12_real64
        call run_v(cfg, G, N, N, rc, rt, o)
        call run_exact(G, N, N, rc, rt, K, K, K, 0.1_real64, 1, pe, nbce, nbte, nwe, ierre)
        call assert_equal_int(o%ierr, ERR_OK, "band parity: ierr")
        call assert_equal_array_real(o%p, pe, G, EXACT, "band parity: p-values identical to exact fixed-k")
        call assert_equal_array_int(o%nbc, nbce, G, "band parity: case pools identical")
        call assert_true(all(o%nsc == K) .and. all(o%nac == K), "band parity: exactly k slots, k accepted")
        call assert_true(maxval(o%u) - minval(o%u) > 0.1_real64, "band parity: fixture heterogeneous in u")
    end subroutine test_band_parity_with_exact

    !| With lambda = 1 and a constant target T, every gene's variance is T: the
    !| scaling multiplies the statistic and every residual by 1/sqrt(T). For T = 1
    !| and T = 4 that factor is a power of two, so the p-values are bit-identical.
    subroutine test_scaling_parity_lambda_one()
        integer(int32), parameter :: G = 300, N = 3
        real(real64) :: rc(N, G), rt(N, G), pe(G)
        integer(int32) :: nbce(G), nbte(G), nwe, ierre, s, it
        real(real64) :: targets(2)
        integer(int32) :: scalings(2)
        type(vcfg) :: cfg
        type(vout) :: o

        call make_hetero(G, N, 0.0_real64, 30.0_real64, rc)
        call make_hetero(G, N, 1.7_real64, 0.0_real64, rt)
        call run_exact(G, N, N, rc, rt, 20, 1, 50, 0.1_real64, 1, pe, nbce, nbte, nwe, ierre)
        targets = [1.0_real64, 4.0_real64]
        scalings = [SCALING_XFIT_JS, SCALING_SAME_GROUP]
        do s = 1, 2
            do it = 1, 2
                cfg%scaling = scalings(s)
                call run_v(cfg, G, N, N, rc, rt, o, lambda_fixed=1.0_real64, target_fixed=targets(it))
                call assert_equal_int(o%ierr, ERR_OK, "scaling parity: ierr")
                call assert_equal_array_real(o%p, pe, G, EXACT, "scaling parity: lambda = 1, T = " // &
                                             str(int(targets(it))) // ", scaling " // str(scalings(s)) // &
                                             " == exact model")
                call assert_true(all(o%s2s == targets(it)), "scaling parity: s2s = T")
            end do
        end do
    end subroutine test_scaling_parity_lambda_one

    subroutine test_adaptive_parity_with_band()
        integer(int32), parameter :: G = 200, N = 4, K = 30
        real(real64) :: rc(N, G), rt(N, G), inf
        type(vcfg) :: cfg
        type(vout) :: ob, oa

        inf = ieee_value(1.0_real64, ieee_positive_inf)
        call make_hetero(G, N, 0.0_real64, 30.0_real64, rc)
        call make_hetero(G, N, 1.7_real64, 0.0_real64, rt)
        cfg%mode = VARIANCE_MODE_BAND
        cfg%k = K
        cfg%k_min = K
        cfg%c_tau = 1.5_real64
        call run_v(cfg, G, N, N, rc, rt, ob)
        cfg%mode = VARIANCE_MODE_ADAPTIVE
        cfg%b = 1
        cfg%r_max = inf
        cfg%tau_stab = inf
        cfg%a_min = 0.0_real64
        call run_v(cfg, G, N, N, rc, rt, oa)
        call assert_equal_int(oa%ierr, ERR_OK, "adaptive parity: ierr")
        call assert_true(ob%nwp > 0, "adaptive parity: band tests genes")
        call assert_equal_array_real(oa%p, ob%p, G, EXACT, "adaptive parity: neutral adaptive == band")
        call assert_equal_array_int(oa%nsc, ob%nsc, G, "adaptive parity: same scans")
        call assert_equal_array_int(oa%nac, ob%nac, G, "adaptive parity: same accepted")
        call assert_true(all(oa%stc == STOP_K_MAX .or. oa%stc == STOP_SCAN_CAP .or. oa%stc == STOP_EXHAUSTED), &
                         "adaptive parity: only k / scan-cap / exhausted stops")
    end subroutine test_adaptive_parity_with_band

    ! =========================================================================
    ! 5. adaptive stop reasons (hand-built structure, target gene 6)
    ! =========================================================================
    subroutine test_adaptive_stop_reasons()
        integer(int32), parameter :: G = 12
        real(real64) :: u(G), pool(100), inf, rad
        type(sorted_data_t) :: sd
        integer(int32) :: n_pooled, n_gp, n_scan, stp

        inf = ieee_value(1.0_real64, ieee_positive_inf)
        call hand_built(sd)
        u = 0.0_real64

        ! k reached: 3 of 3 accepted
        call gather_residuals_adaptive_helper(6.0_real64, 0.0_real64, sd, u, 1.0_real64, 1, 1, 3, inf, inf, &
                                              2, 0.0_real64, 1000, pool, n_pooled, n_gp, n_scan, rad, stp, 100)
        call assert_equal_int(stp, STOP_K_MAX, "adaptive: k reached")
        call assert_equal_int(n_gp, 3, "adaptive: 3 genes at k = 3")
        ! radius: |log2(6) - log2(7)| = 0.222 > 0.21, while gene 7 (0.193) is inside
        call gather_residuals_adaptive_helper(6.0_real64, 0.0_real64, sd, u, 1.0_real64, 1, 1, 10, 0.21_real64, &
                                              inf, 2, 0.0_real64, 1000, pool, n_pooled, n_gp, n_scan, rad, stp, 100)
        call assert_equal_int(stp, STOP_RADIUS, "adaptive: radius stop")
        call assert_equal_int(n_gp, 1, "adaptive: radius stop before gene 5")
        ! stability: residual scale grows with the gene index; after genes 6, 5 the
        ! third (gene 7) moves mean|r| by more than 1 %
        call gather_residuals_adaptive_helper(6.0_real64, 0.0_real64, sd, u, 1.0_real64, 1, 2, 10, inf, &
                                              0.01_real64, 2, 0.0_real64, 1000, pool, n_pooled, n_gp, n_scan, &
                                              rad, stp, 100)
        call assert_equal_int(stp, STOP_STABILITY, "adaptive: stability stop")
        call assert_equal_int(n_gp, 2, "adaptive: the destabilising block is discarded")
        call assert_equal_int(n_pooled, 6, "adaptive: discarded residuals removed")
        ! acceptance rate: only the target passes; blocks 2 and 3 accept 0 of 2
        u = 5.0_real64
        u(6) = 0.0_real64
        call gather_residuals_adaptive_helper(6.0_real64, 0.0_real64, sd, u, 0.5_real64, 1, 1, 10, inf, inf, &
                                              2, 0.5_real64, 1000, pool, n_pooled, n_gp, n_scan, rad, stp, 100)
        call assert_equal_int(stp, STOP_ACCEPTANCE, "adaptive: acceptance-rate stop")
        call assert_equal_int(n_scan, 3, "adaptive: stops after the third block")
        ! scan cap
        call gather_residuals_adaptive_helper(6.0_real64, 0.0_real64, sd, u, 0.5_real64, 1, 1, 10, inf, inf, &
                                              2, 0.0_real64, 4, pool, n_pooled, n_gp, n_scan, rad, stp, 100)
        call assert_equal_int(stp, STOP_SCAN_CAP, "adaptive: scan cap")
        call assert_equal_int(n_scan, 4, "adaptive: 4 slots at scan cap 4")
        ! exhausted: everything accepted, k larger than the gene count
        u = 0.0_real64
        call gather_residuals_adaptive_helper(6.0_real64, 0.0_real64, sd, u, 1.0_real64, 5, 1, 50, inf, inf, &
                                              2, 0.0_real64, 1000, pool, n_pooled, n_gp, n_scan, rad, stp, 100)
        call assert_equal_int(stp, STOP_EXHAUSTED, "adaptive: candidates exhausted")
        call assert_equal_int(n_gp, G, "adaptive: every gene accepted")
        call assert_equal_real(rad, log(7.0_real64 / 2.0_real64) / log(2.0_real64), 1d-12, &
                               "adaptive: radius = log2 distance of the farthest accepted gene (gene 1)")
    end subroutine test_adaptive_stop_reasons

    ! =========================================================================
    ! 6. shrinkage
    ! =========================================================================
    !| Linear branch, every gene with true variance 1. With the target fixed at the
    !| truth, E[sum Vhat] = sum Var(s2) = E[sum (s2 - 1)^2], so the unclipped ratios
    !| (shrink_diag 4:6) are about 1 iff the code's `Vhat` is unbiased.
    subroutine test_vhat_unbiased()
        integer(int32), parameter :: G = 200000
        integer(int32) :: nlist(3), i, dist, nc, nt, nf, ierr, seed_size, jg
        integer(int32), allocatable :: seed(:)
        real(real64), allocatable :: rc(:, :), rt(:, :), z(:), c5(:)
        real(real64), allocatable :: s2c(:), s2t(:), s2(:), m(:), tgt(:), sc(:), st(:), sp(:)
        real(real64) :: diag(N_SHRINK_DIAG), tol, y(3), v1, s2_3
        character(len=8) :: names(3)

        allocate(s2c(G), s2t(G), s2(G), m(G), tgt(G), sc(G), st(G), sp(G))
        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        seed = 20261028
        call random_seed(put=seed)
        names = [character(len=8) :: "normal", "skewed", "t10"]
        nlist = [3, 5, 10]
        do dist = 1, 3
            do i = 1, 4
                if (i <= 3) then
                    nc = nlist(i); nt = nlist(i)
                else
                    if (dist /= 1) cycle
                    nc = 3; nt = 5                         ! unbalanced: exact for normal data
                end if
                allocate(rc(nc, G), rt(nt, G), z(max(nc, nt) * G), c5(max(nc, nt) * G))
                call draw(rc, nc)
                call draw(rt, nt)
                call compute_group_variances_helper(rc, nc, rt, nt, G, 0, s2c, s2t, s2, m)
                call compute_moment_shrinkage(s2c, s2t, s2, m, G, nc, nt, SHRINK_TARGET_MEDIAN, LAMBDA_ABSOLUTE, &
                                              0.3_real64, 1, 0, tgt, sc, st, sp, diag, nf, ierr, &
                                              target_fixed=1.0_real64)
                call assert_equal_int(ierr, ERR_OK, "vhat: ierr")
                tol = merge(0.08_real64, 0.04_real64, dist == 3)
                call assert_equal_real(diag(4), 1.0_real64, tol, "vhat: pooled ratio ~ 1, "//trim(names(dist)) &
                                       //" n = "//str(nc)//"/"//str(nt))
                call assert_equal_real(diag(5), 1.0_real64, tol, "vhat: case ratio ~ 1, "//trim(names(dist)))
                call assert_equal_real(diag(6), 1.0_real64, tol, "vhat: control ratio ~ 1, "//trim(names(dist)))
                deallocate(rc, rt, z, c5)
            end do
        end do

        ! documentation: at n = 3 the v1 estimator n/(n-1)^3 sum (w - wbar)^2 equals s2^2/4
        y = [0.3_real64, 1.9_real64, -0.7_real64]
        s2_3 = sum((y - sum(y) / 3.0_real64)**2) / 2.0_real64
        v1 = 3.0_real64 / 8.0_real64 * sum(((y - sum(y) / 3.0_real64)**2 &
             - sum((y - sum(y) / 3.0_real64)**2) / 3.0_real64)**2)
        call assert_equal_real(v1, s2_3**2 / 4.0_real64, 1d-12, "vhat: v1 estimator = s2^2/4 at n = 3")
    contains
        !> Unit-variance draws: normal, centred exponential (skewed), t10 / sqrt(10/8).
        !| Heavy tail is t10, not note F's t5: the MC check needs a finite 8th moment
        !| (the variance of the sums of delta^2 and (s2 - T)^2); for t5 it is
        !| infinite and the ratio does not settle at any feasible gene count.
        subroutine draw(x, n)
            integer(int32), intent(in) :: n
            real(real64), intent(out) :: x(n, G)
            real(real64), allocatable :: u(:)
            real(real64) :: w(10)
            integer(int32) :: j, kk
            select case (dist)
            case (1)
                call rnorm(z(1:n * G))
            case (2)
                allocate(u(n * G))
                call random_number(u)
                z(1:n * G) = -log(1.0_real64 - u) - 1.0_real64
            case default
                call rnorm(z(1:n * G))
                do j = 1, n * G
                    call rnorm(w)
                    c5(j) = sum(w**2)
                end do
                z(1:n * G) = z(1:n * G) / sqrt(c5(1:n * G) / 10.0_real64) / sqrt(10.0_real64 / 8.0_real64)
            end select
            do jg = 1, G
                do kk = 1, n
                    x(kk, jg) = 100.0_real64 + z((jg - 1) * n + kk)
                end do
            end do
        end subroutine draw
    end subroutine test_vhat_unbiased

    subroutine test_lambda_behaviour()
        integer(int32), parameter :: G = 5000
        integer(int32) :: nf, ierr, jg, n
        real(real64), allocatable :: rc(:, :), rt(:, :), z(:)
        real(real64) :: sd(G), s2c(G), s2t(G), s2(G), m(G), tgt(G), sc(G), st(G), sp(G)
        real(real64) :: drel(N_SHRINK_DIAG), dabs(N_SHRINK_DIAG)

        ! homogeneous (all variances 1) at n = 3: lambda near 1
        n = 3
        allocate(rc(n, G), rt(n, G), z(2 * n * G))
        call rnorm(z)
        do jg = 1, G
            rc(:, jg) = 50.0_real64 + z((jg - 1) * 2 * n + 1 : (jg - 1) * 2 * n + n)
            rt(:, jg) = 50.0_real64 + z((jg - 1) * 2 * n + n + 1 : jg * 2 * n)
        end do
        call compute_group_variances_helper(rc, n, rt, n, G, 0, s2c, s2t, s2, m)
        call compute_moment_shrinkage(s2c, s2t, s2, m, G, n, n, SHRINK_TARGET_MEDIAN, LAMBDA_RELATIVE, &
                                      0.3_real64, 1, 0, tgt, sc, st, sp, drel, nf, ierr)
        call assert_equal_int(ierr, ERR_OK, "lambda: ierr")
        call assert_true(drel(1) > 0.85_real64, "lambda: near 1 without heterogeneity")
        call assert_true(drel(7) > 2.0_real64, "lambda: Bock condition p_hat > 2")
        ! relative == absolute for a constant target
        call compute_moment_shrinkage(s2c, s2t, s2, m, G, n, n, SHRINK_TARGET_MEDIAN, LAMBDA_ABSOLUTE, &
                                      0.3_real64, 1, 0, tgt, sc, st, sp, dabs, nf, ierr)
        call assert_equal_real(drel(1), dabs(1), 1d-12, "lambda: relative == absolute for a constant target")
        deallocate(rc, rt, z)

        ! heterogeneous (log-variance scatter 1.05) at n = 10: clearly below 1
        n = 10
        allocate(rc(n, G), rt(n, G), z(2 * n * G))
        call rnorm(z)
        call rnorm(sd)
        sd = exp(0.525_real64 * sd)
        do jg = 1, G
            rc(:, jg) = 50.0_real64 + sd(jg) * z((jg - 1) * 2 * n + 1 : (jg - 1) * 2 * n + n)
            rt(:, jg) = 50.0_real64 + sd(jg) * z((jg - 1) * 2 * n + n + 1 : jg * 2 * n)
        end do
        call compute_group_variances_helper(rc, n, rt, n, G, 0, s2c, s2t, s2, m)
        call compute_moment_shrinkage(s2c, s2t, s2, m, G, n, n, SHRINK_TARGET_MEDIAN, LAMBDA_RELATIVE, &
                                      0.3_real64, 1, 0, tgt, sc, st, sp, drel, nf, ierr)
        call assert_true(drel(1) < 0.5_real64, "lambda: clearly below 1 with heterogeneity at n = 10")
        call assert_true(drel(1) >= 0.0_real64, "lambda: clipped at 0")
        call assert_true(all(sp >= min(tgt(1), minval(s2)) - 1d-12) .and. &
                         all(sp <= max(tgt(1), maxval(s2)) + 1d-12), &
                         "lambda: shrunk variances between the raw values and the target")
    end subroutine test_lambda_behaviour

    ! =========================================================================
    ! 7. calibration regression (note F test 5)
    ! =========================================================================
    !| Null data, n = 3, linear branch, normal noise, a mean trend and gene-wise SDs
    !| scattering around it (log-variance SD 1.05). Mean-matched neighbours then
    !| differ in noise, which is what cross-fitting has to handle.
    !|
    !| Fixed k = 300 isolates the internal-studentisation mechanism this test guards
    !| from the finite-pool tail effect. Measured on this fixture (28/10/2026), rate
    !| of p <= 0.001 for none / xfit_js / xfit_raw: adaptive k20_50 0.0053 / 0.0033 /
    !| 0.0015; fixed 50 0.0041 / 0.0027 / 0.0005; fixed 100 0.0034 / 0.0024 / 0.0003;
    !| fixed 300 0.0026 / 0.0017 / 0.0001. The excess shrinks with the pool for every
    !| scaling (the empirical null under-represents its own extreme tail), so it is
    !| not a cross-fitting defect; same_group is 0.024 at k20_50.
    subroutine test_xfit_calibration_regression()
        integer(int32), parameter :: G = 20000, N = 3
        real(real64), allocatable :: rc(:, :), rt(:, :), z(:), sd(:)
        real(real64) :: rate_x, rate_s
        integer(int32) :: jg, seed_size
        integer(int32), allocatable :: seed(:)
        type(vcfg) :: cfg
        type(vout) :: o

        allocate(rc(N, G), rt(N, G), z(2 * N * G), sd(G))
        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        seed = 4242
        call random_seed(put=seed)
        call rnorm(z)
        call rnorm(sd)
        sd = exp(0.525_real64 * sd)
        do jg = 1, G
            rc(:, jg) = 200.0_real64 + 0.05_real64 * real(jg, real64) + sd(jg) * z((jg - 1) * 2 * N + 1 : (jg - 1) * 2 * N + N)
            rt(:, jg) = 200.0_real64 + 0.05_real64 * real(jg, real64) + sd(jg) * z((jg - 1) * 2 * N + N + 1 : jg * 2 * N)
        end do
        cfg%norm = 0
        cfg%k_start = 300
        cfg%k_step = 300
        cfg%k_max = 300
        cfg%scaling = SCALING_XFIT_JS
        call run_v(cfg, G, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_OK, "calibration: xfit ierr")
        rate_x = real(count(o%p >= 0.0_real64 .and. o%p <= 0.001_real64), real64) &
                 / real(count(o%p >= 0.0_real64), real64)
        cfg%scaling = SCALING_SAME_GROUP
        call run_v(cfg, G, N, N, rc, rt, o)
        rate_s = real(count(o%p >= 0.0_real64 .and. o%p <= 0.001_real64), real64) &
                 / real(count(o%p >= 0.0_real64), real64)
        print '(A,F8.5,A,F8.5)', "    calibration: p <= 0.001 rate, cross-fit ", rate_x, ", same-group ", rate_s
        call assert_true(rate_x <= 0.002_real64, "calibration: cross-fit tail p <= 0.001 at most 2x nominal")
        call assert_true(rate_s > 0.004_real64, "calibration: same-group (internally studentised) tail too heavy")
    end subroutine test_xfit_calibration_regression

    ! =========================================================================
    ! 8. scale invariance, floor, validation
    ! =========================================================================
    subroutine test_scale_invariance_and_floor()
        integer(int32), parameter :: G = 300, N = 3
        real(real64) :: rc(N, G), rt(N, G), p1(G)
        type(vcfg) :: cfg
        type(vout) :: o

        call make_hetero(G, N, 0.0_real64, 30.0_real64, rc)
        call make_hetero(G, N, 1.7_real64, 0.0_real64, rt)
        cfg%norm = 0
        cfg%scaling = SCALING_XFIT_JS
        call run_v(cfg, G, N, N, rc, rt, o)
        p1 = o%p
        call run_v(cfg, G, N, N, 4.0_real64 * rc, 4.0_real64 * rt, o)
        call assert_equal_array_real(o%p, p1, G, EXACT, "scale invariance: data x 4 gives the same p-values")

        ! floor: a gene constant in the control group under xfit_raw
        rt(:, 150) = rt(1, 150)
        cfg%scaling = SCALING_XFIT_RAW
        call run_v(cfg, G, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_OK, "floor: ierr")
        call assert_true(o%nf >= 1, "floor: the constant gene is counted")
        call assert_true(all(o%p == -1.0_real64 .or. (o%p >= 0.0_real64 .and. o%p <= 1.0_real64)), &
                         "floor: every p-value valid")
    end subroutine test_scale_invariance_and_floor

    subroutine test_variance_validation()
        integer(int32), parameter :: G = 50, N = 3
        real(real64) :: rc(N, G), rt(N, G), rbad(G, N)
        type(vcfg) :: cfg, bad
        type(vout) :: o
        integer(int32) :: valid(G), nwp, nb(G, 8), nf, ierr
        real(real64) :: p(G), obs(G), l2(G), u(G), s2s(G), rad(G, 2), dh, tv, sc, dg(N_SHRINK_DIAG)

        call make_hetero(G, N, 0.0_real64, 0.0_real64, rc)
        call make_hetero(G, N, 1.0_real64, 0.0_real64, rt)

        bad = cfg; bad%mode = 2
        call run_v(bad, G, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_INVALID_INPUT, "validation: removed mode 2 (note B) rejected")
        bad = cfg; bad%match = 2
        call run_v(bad, G, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_INVALID_INPUT, "validation: unknown match_variable")
        bad = cfg; bad%scaling = 4
        call run_v(bad, G, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_INVALID_INPUT, "validation: unknown scaling")
        bad = cfg; bad%c_js = 0.0_real64
        call run_v(bad, G, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_INVALID_INPUT, "validation: c_js must be > 0")
        bad = cfg; bad%k_min = 60
        call run_v(bad, G, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_INVALID_INPUT, "validation: k_min > k rejected")
        bad = cfg; bad%a_min = 1.5_real64
        call run_v(bad, G, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_INVALID_INPUT, "validation: a_min outside [0, 1]")
        bad = cfg; bad%srobust = 2
        call run_v(bad, G, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_INVALID_INPUT, "validation: shrink_loess_robust must be 0/1")
        bad = cfg; bad%k = 0
        call run_v(bad, G, N, N, rc, rt, o)
        call assert_true(o%ierr /= ERR_OK, "validation: k = 0 rejected")

        rbad = transpose(rt)
        valid = 1
        call compute_noise_pvalue_pipeline_variance(rc, G, N, rbad, N, G, valid, BETA_CENTRE_NONE, G, 1, &
            0, 0, 0, 0, 0, 20, 1, 50, 0.1_real64, 50, 1.5_real64, 1.0_real64, 5, 20, 0.5_real64, 0.1_real64, &
            2, 0.2_real64, 0.3_real64, 1, 1, 0, 5000, p, obs, l2, dh, nwp, nb(:, 1), nb(:, 2), nb(:, 3), &
            nb(:, 4), nb(:, 5), nb(:, 6), rad(:, 1), rad(:, 2), nb(:, 7), nb(:, 8), u, s2s, tv, sc, dg, nf, ierr)
        call assert_equal_int(ierr, ERR_DIM_MISMATCH, "validation: transposed control matrix")
    end subroutine test_variance_validation

    ! =========================================================================
    ! 9. 2-D ellipse
    ! =========================================================================
    !> Random structure: G genes, linear means over a wide range, 3 replicates.
    subroutine random_structure(G, sd, u, seedv)
        integer(int32), intent(in) :: G, seedv
        type(sorted_data_t), intent(out) :: sd
        real(real64), intent(out) :: u(G)
        real(real64) :: means(G), reps(3, G), z(3 * G)
        integer(int32) :: ig, ierr, seed_size
        integer(int32), allocatable :: seed(:)
        call random_seed(size=seed_size); allocate(seed(seed_size)); seed = seedv; call random_seed(put=seed)
        call rnorm(z)
        call random_number(means)
        means = 2.0_real64**(12.0_real64 * means) - 1.0_real64
        do ig = 1, G
            reps(:, ig) = means(ig) + (1.0_real64 + 0.1_real64 * means(ig)) * z(3 * ig - 2 : 3 * ig)
            means(ig) = sum(reps(:, ig)) / 3.0_real64
        end do
        call rnorm(u)
        call prepare_sorted_data(means, reps, 3, G, 0, sd, ierr)
        call assert_equal_int(ierr, ERR_OK, "ellipse: prepare_sorted_data")
    end subroutine random_structure

    subroutine test_ellipse_brute_force()
        integer(int32), parameter :: G = 600
        type(sorted_data_t) :: sd
        real(real64) :: u(G), pool(5000), expect(5000), d(G), dm(G), rad, r_list(3), t_list(3), lt
        integer(int32) :: it, ir, ik, j, m, cnt, kk, n_pooled, n_expect, n_gp, n_scan, stp, n_in
        integer(int32) :: order(G), k_list(2)

        call random_structure(G, sd, u, 77)
        r_list = [0.5_real64, 1.5_real64, 1.0e6_real64]
        t_list = [0.3_real64, 1.0_real64, 1.0e12_real64]
        k_list = [10, 50]
        do it = 1, 40
            j = 1 + mod(37 * it, G)                       ! target: sorted slot j
            lt = log(max(sd%means_sorted(j), 0.0_real64) + 1.0_real64) / log(2.0_real64)
            dm = abs(log(max(sd%means_sorted, 0.0_real64) + 1.0_real64) / log(2.0_real64) - lt)
            do ir = 1, 3
                d = sqrt((dm / r_list(ir))**2 + ((u(sd%original_indices) - u(sd%original_indices(j))) &
                         / t_list(ir))**2)
                n_in = count(d <= 1.0_real64)
                order = [(m, m = 1, G)]
                call sort_by(d, order, G)                 ! brute force: slots by ascending D
                do ik = 1, 2
                    kk = k_list(ik)
                    cnt = 0
                    do m = 1, G
                        if (d(order(m)) > 1.0_real64 .or. cnt >= kk) exit
                        cnt = cnt + 1
                        expect(3 * cnt - 2 : 3 * cnt) = sd%residuals_packed(:, order(m))
                    end do
                    n_expect = 3 * cnt
                    call gather_residuals_ellipse_helper(sd%means_sorted(j), u(sd%original_indices(j)), sd, u, &
                                                         t_list(ir), r_list(ir), kk, 100000, pool, n_pooled, n_gp, &
                                                         n_scan, rad, stp, 5000)
                    call assert_equal_int(n_pooled, n_expect, "ellipse: pool size = brute force, target "//str(it))
                    call assert_equal_array_real(pool(1:n_pooled), expect(1:n_expect), n_expect, EXACT, &
                                                 "ellipse: pool = brute-force k nearest inside the ellipse")
                    call assert_equal_int(n_gp, min(n_in, kk), "ellipse: genes pooled")
                    if (ir == 1) call assert_true(n_scan < G, "ellipse: the walk stops early inside a small ellipse")
                end do
            end do
        end do
    end subroutine test_ellipse_brute_force

    !> Stable insertion sort of `idx` by `key(idx)` (test helper; small n).
    subroutine sort_by(key, idx, n)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: key(n)
        integer(int32), intent(inout) :: idx(n)
        integer(int32) :: i, j, t
        do i = 2, n
            t = idx(i); j = i - 1
            do while (j >= 1)
                if (key(idx(j)) <= key(t)) exit
                idx(j + 1) = idx(j); j = j - 1
            end do
            idx(j + 1) = t
        end do
    end subroutine sort_by

    subroutine test_ellipse_stops_and_pipeline()
        integer(int32), parameter :: G = 12, GP = 200, N = 4
        type(sorted_data_t) :: sd
        real(real64) :: u(G), pool(100), rad, inf, rc(N, GP), rt(N, GP)
        integer(int32) :: n_pooled, n_gp, n_scan, stp
        type(vcfg) :: cfg
        type(vout) :: o

        inf = ieee_value(1.0_real64, ieee_positive_inf)
        call hand_built(sd)
        u = 0.0_real64
        call gather_residuals_ellipse_helper(6.0_real64, 0.0_real64, sd, u, 1.0_real64, 10.0_real64, 3, 1000, &
                                             pool, n_pooled, n_gp, n_scan, rad, stp, 100)
        call assert_equal_int(stp, STOP_K_MAX, "ellipse: k nearest found")
        call assert_equal_int(n_gp, 3, "ellipse: 3 genes")
        ! |log2(6) - log2(7)| = 0.222 (gene 5), 0.193 (gene 7): semi-axis 0.21 keeps only genes 6, 7
        call gather_residuals_ellipse_helper(6.0_real64, 0.0_real64, sd, u, 1.0_real64, 0.21_real64, 10, 1000, &
                                             pool, n_pooled, n_gp, n_scan, rad, stp, 100)
        call assert_equal_int(stp, STOP_RADIUS, "ellipse: nothing left inside, short pool")
        call assert_equal_int(n_gp, 2, "ellipse: genes 6 and 7 inside")
        call gather_residuals_ellipse_helper(6.0_real64, 0.0_real64, sd, u, 1.0_real64, inf, 10, 4, &
                                             pool, n_pooled, n_gp, n_scan, rad, stp, 100)
        call assert_equal_int(stp, STOP_SCAN_CAP, "ellipse: scan cap")
        call gather_residuals_ellipse_helper(6.0_real64, 0.0_real64, sd, u, 1.0_real64, inf, 50, 1000, &
                                             pool, n_pooled, n_gp, n_scan, rad, stp, 100)
        call assert_equal_int(stp, STOP_EXHAUSTED, "ellipse: exhausted")
        call assert_equal_int(n_gp, G, "ellipse: all genes")

        call make_hetero(GP, N, 0.0_real64, 30.0_real64, rc)
        call make_hetero(GP, N, 1.7_real64, 0.0_real64, rt)
        cfg%mode = VARIANCE_MODE_ELLIPSE
        cfg%k = 30
        cfg%k_min = 30
        cfg%r_max = 1.0_real64
        cfg%c_tau = 3.0_real64
        call run_v(cfg, GP, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_OK, "ellipse pipeline: ierr")
        call assert_true(o%nwp > 0, "ellipse pipeline: p-values computed")
        call assert_true(all(o%nac(1:GP) <= 30), "ellipse pipeline: at most k genes")
        call assert_true(all(o%p == -1.0_real64 .or. (o%p >= 0.0_real64 .and. o%p <= 1.0_real64)), &
                         "ellipse pipeline: valid p-values")
        cfg%scaling = SCALING_XFIT_JS
        cfg%match = MATCH_JS
        call run_v(cfg, GP, N, N, rc, rt, o)
        call assert_equal_int(o%ierr, ERR_OK, "ellipse pipeline: JS + cross-fit ierr")
    end subroutine test_ellipse_stops_and_pipeline

end module mod_test_noise_model_variance
