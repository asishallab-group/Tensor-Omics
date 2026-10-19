!> Unit tests for the variance-matched noise model (`noise_model_variance`,
!| implementation notes A and B).
!|
!|   1. Polygamma helpers against closed forms, and the note-A tau rule
!|      `sqrt(trigamma(d/2))` = 0.80 / 0.53 / 0.34 at n = 3 / 5 / 10 per group.
!|   2. Hand-built band gather: exactly which genes enter, in which order, the scan
!|      count, `scan_max`, short pools, non-finite candidates.
!|   3. Parity (note A §4.2): mode A with an all-accepting band reproduces
!|      `noise_model_exact` with `k_start = k_step = k_max = k` EXACTLY.
!|   4. Degradation (note B §5.1): mode B on data with identical per-gene noise has
!|      `d0 = Inf`, `w = 0`, `u = 0`, and reproduces the same plain fixed-k model.
!|   5. A band so strict that no pool reaches `BAND_MIN_GENES` drops every gene and
!|      reports the full scan.
!|   6. Argument validation.
module mod_test_noise_model_variance
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
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

contains

    function get_all_tests_noise_model_variance() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate(all_tests(6))
        all_tests(1) = test_case("test_polygamma_and_tau_rule",  test_polygamma_and_tau_rule)
        all_tests(2) = test_case("test_band_gather_hand_built",  test_band_gather_hand_built)
        all_tests(3) = test_case("test_mode_a_parity_with_exact", test_mode_a_parity_with_exact)
        all_tests(4) = test_case("test_mode_b_degrades_to_exact", test_mode_b_degrades_to_exact)
        all_tests(5) = test_case("test_strict_band_drops_genes",  test_strict_band_drops_genes)
        all_tests(6) = test_case("test_variance_validation",      test_variance_validation)
    end function get_all_tests_noise_model_variance

    function str(i) result(s)
        integer(int32), intent(in) :: i
        character(len=20) :: buffer
        character(len=:), allocatable :: s
        write(buffer, '(I0)') i
        s = trim(buffer)
    end function str

    !> Heterogeneous fixture: `n_genes` genes over a wide mean range, per-gene noise
    !| amplitude varying smoothly AND gene-to-gene (so the variance band has
    !| something to select on), deterministic.
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

    !> Homogeneous fixture: on L = log2(x + 1) every gene has the SAME deviation
    !| pattern around its own level, so every gene's pooled log2 variance is equal.
    subroutine make_homog(n_genes, n_rep, offset, level_shift, reps)
        integer(int32), intent(in) :: n_genes, n_rep
        real(real64), intent(in) :: offset, level_shift
        real(real64), intent(out) :: reps(n_rep, n_genes)
        integer(int32) :: g, r
        real(real64) :: lvl
        do g = 1, n_genes
            lvl = 3.0_real64 + 7.0_real64 * real(g - 1, real64) / real(n_genes - 1, real64) + level_shift
            do r = 1, n_rep
                reps(r, g) = 2.0_real64**(lvl + 0.3_real64 * sin(offset + 2.1_real64 * real(r, real64))) - 1.0_real64
            end do
        end do
    end subroutine make_homog

    ! =========================================================================
    ! 1. polygamma helpers and the tau rule
    ! =========================================================================
    subroutine test_polygamma_and_tau_rule()
        integer(int32), parameter :: G = 300
        integer(int32) :: n_list(3), i, n, nf, ierr
        real(real64) :: expect(3), y
        real(real64), allocatable :: rc(:, :), rt(:, :)
        real(real64) :: u(G), tau, d0, w, sc

        call assert_equal_real(trigamma_helper(1.0_real64), PI**2 / 6.0_real64, 1d-10, "trigamma(1)")
        call assert_equal_real(trigamma_helper(0.5_real64), PI**2 / 2.0_real64, 1d-10, "trigamma(1/2)")
        call assert_equal_real(trigamma_helper(2.0_real64), PI**2 / 6.0_real64 - 1.0_real64, 1d-10, "trigamma(2)")
        call assert_equal_real(digamma_helper(1.0_real64), -0.5772156649015329_real64, 1d-10, "digamma(1)")
        call assert_equal_real(digamma_helper(0.5_real64), -0.5772156649015329_real64 - 2.0_real64 * log(2.0_real64), &
                               1d-10, "digamma(1/2)")
        do i = 1, 4
            y = 10.0_real64**(i - 2)                 ! 0.1, 1, 10, 100
            call assert_true(abs(trigamma_inverse_helper(trigamma_helper(y)) - y) <= 1d-9 * y, &
                             "trigamma inverse round trip at y = 10^"//str(i - 2))
        end do

        ! tau rule of note A: c * sqrt(trigamma(d/2)), d = 2n - 2.
        n_list = [3, 5, 10]
        expect = [0.80_real64, 0.53_real64, 0.34_real64]
        do i = 1, 3
            n = n_list(i)
            allocate(rc(n, G), rt(n, G))
            call make_hetero(G, n, 0.0_real64, 0.0_real64, rc)
            call make_hetero(G, n, 1.0_real64, 0.0_real64, rt)
            call compute_variance_match(rc, n, rt, n, G, VARIANCE_MODE_BAND, 1.0_real64, 0.3_real64, 1, 1, 1, &
                                        u, tau, d0, w, sc, nf, ierr)
            call assert_equal_int(ierr, ERR_OK, "tau rule: ierr at n = "//str(n))
            call assert_equal_real(tau, sqrt(trigamma_helper(real(n - 1, real64))), TOL, &
                                   "tau rule: tau is c*sqrt(trigamma(d/2)) at n = "//str(n))
            call assert_equal_real(tau, expect(i), 0.006_real64, "tau rule: note-A value at n = "//str(n))
            call assert_equal_real(d0, -1.0_real64, EXACT, "mode A reports d0 = -1")
            deallocate(rc, rt)
        end do
    end subroutine test_polygamma_and_tau_rule

    ! =========================================================================
    ! 2. hand-built band gather
    ! =========================================================================
    !| 12 genes with means 1..12 (3 replicates, residuals -1/0/+1 scaled by the
    !| gene index so pools are traceable). Target = gene 6 (mean 6, u = 0). Visit
    !| order from slot 6 under `choose_index` (ties go left): 6, 5, 7, 4, 8, 3, 9,
    !| 2, 10, 1, 11, 12.
    subroutine test_band_gather_hand_built()
        integer(int32), parameter :: G = 12, NR = 3
        real(real64) :: means(G), reps(NR, G), u(G), pool(100), nan
        type(sorted_data_t) :: sd
        integer(int32) :: ig, ierr, n_pooled, n_gp, n_scan
        integer(int32), parameter :: expect_genes(4) = [6, 7, 4, 9]
        real(real64) :: expect_pool(12)

        nan = ieee_value(1.0_real64, ieee_quiet_nan)
        do ig = 1, G
            means(ig) = real(ig, real64)
            reps(:, ig) = means(ig) + real(ig, real64) * [-1.0_real64, 0.0_real64, 1.0_real64]
        end do
        call prepare_sorted_data(means, reps, NR, G, 0, sd, ierr)
        call assert_equal_int(ierr, ERR_OK, "hand-built: prepare_sorted_data")

        ! in band (|u| <= 0.5): 6, 7, 4, 9 (9 exactly on the boundary); out: 5 (far),
        ! 8 (NaN, never a neighbour), 3 (-0.6), all others far.
        u = 5.0_real64
        u(6) = 0.0_real64; u(5) = 5.0_real64; u(7) = 0.2_real64; u(4) = -0.4_real64
        u(8) = nan; u(3) = -0.6_real64; u(9) = 0.5_real64

        call gather_residuals_band_helper(6.0_real64, 0.0_real64, sd, u, 4, 0.5_real64, 1000, &
                                          pool, n_pooled, n_gp, n_scan, 100)
        call assert_equal_int(n_gp, 4, "hand-built: 4 genes accepted")
        call assert_equal_int(n_scan, 7, "hand-built: visits 6,5,7,4,8,3,9")
        call assert_equal_int(n_pooled, 12, "hand-built: 4 genes x 3 residuals")
        ! residuals are Bessel-scaled deviations: g * [-1,0,1] * sqrt(3/2)
        do ig = 1, 4
            expect_pool(3 * ig - 2 : 3 * ig) = real(expect_genes(ig), real64) * sqrt(1.5_real64) &
                                             * [-1.0_real64, 0.0_real64, 1.0_real64]
        end do
        call assert_equal_array_real(pool(1:12), expect_pool, 12, TOL, &
                                     "hand-built: pool = genes 6, 7, 4, 9 in visit order")

        ! scan_max stops the walk: 3 slots -> genes 6 (in), 5 (out), 7 (in)
        call gather_residuals_band_helper(6.0_real64, 0.0_real64, sd, u, 4, 0.5_real64, 3, &
                                          pool, n_pooled, n_gp, n_scan, 100)
        call assert_equal_int(n_scan, 3, "hand-built: scan_max = 3 stops after 3 slots")
        call assert_equal_int(n_gp, 2, "hand-built: 2 genes accepted within scan_max")

        ! short pool: asking for more than are eligible exhausts both sides and keeps
        ! the 4 eligible genes (the band is never relaxed)
        call gather_residuals_band_helper(6.0_real64, 0.0_real64, sd, u, 50, 0.5_real64, 1000, &
                                          pool, n_pooled, n_gp, n_scan, 100)
        call assert_equal_int(n_gp, 4, "hand-built: short pool keeps the 4 eligible genes")
        call assert_equal_int(n_scan, G, "hand-built: short pool scans every slot")

        ! max_pool_size caps the pool and stops the walk
        call gather_residuals_band_helper(6.0_real64, 0.0_real64, sd, u, 4, 0.5_real64, 1000, &
                                          pool, n_pooled, n_gp, n_scan, 5)
        call assert_equal_int(n_pooled, 5, "hand-built: max_pool_size caps the residuals")
        call assert_equal_int(n_gp, 2, "hand-built: stops once the pool is full")
    end subroutine test_band_gather_hand_built

    !> Run the variance pipeline with the suite's defaults.
    subroutine run_var(G, n, rc, rt, mode, k, c_tau, p, obs, l2, nwp, nbc, nbt, nsc, nst, u, tau, d0, w, ierr)
        integer(int32), intent(in) :: G, n, mode, k
        real(real64), intent(in) :: rc(n, G), rt(n, G), c_tau
        real(real64), intent(out) :: p(G), obs(G), l2(G), u(G), tau, d0, w
        integer(int32), intent(out) :: nwp, nbc(G), nbt(G), nsc(G), nst(G), ierr
        integer(int32) :: valid(G), nf
        real(real64) :: dh, sc
        valid = 1
        call compute_noise_pvalue_pipeline_variance(rc, G, n, rt, G, n, valid, BETA_CENTRE_NONE, &
            G, 1, mode, k, c_tau, 0.3_real64, 1, 1, 1, 5000, &
            p, obs, l2, dh, nwp, nbc, nbt, nsc, nst, u, tau, d0, w, sc, nf, ierr)
    end subroutine run_var

    !> Run the plain exact model at a fixed k (k_start = k_step = k_max = k).
    subroutine run_exact_fixed_k(G, n, rc, rt, k, p, obs, l2, nwp, nbc, nbt, ierr)
        integer(int32), intent(in) :: G, n, k
        real(real64), intent(in) :: rc(n, G), rt(n, G)
        real(real64), intent(out) :: p(G), obs(G), l2(G)
        integer(int32), intent(out) :: nwp, nbc(G), nbt(G), ierr
        integer(int32) :: valid(G), nb3(G)
        real(real64) :: dh
        valid = 1
        call pipeline_exact(rc, G, n, rt, G, n, valid, BETA_CENTRE_NONE, G, 1, k, k, k, 0.1_real64, &
                            p, obs, l2, dh, nwp, 5000, nbc, nbt, nb3, ierr)
    end subroutine run_exact_fixed_k

    ! =========================================================================
    ! 3. parity: mode A, all-accepting band == exact at fixed k
    ! =========================================================================
    subroutine test_mode_a_parity_with_exact()
        integer(int32), parameter :: G = 200, N = 4, K = 30
        real(real64) :: rc(N, G), rt(N, G)
        real(real64) :: pv(G), ov(G), lv(G), u(G), tau, d0, w
        real(real64) :: pe(G), oe(G), le(G)
        integer(int32) :: nwv, nbcv(G), nbtv(G), nsc(G), nst(G), ierr
        integer(int32) :: nwe, nbce(G), nbte(G), ierre

        call make_hetero(G, N, 0.0_real64, 30.0_real64, rc)
        call make_hetero(G, N, 1.7_real64, 0.0_real64, rt)

        call run_var(G, N, rc, rt, VARIANCE_MODE_BAND, K, 1.0e12_real64, pv, ov, lv, nwv, nbcv, nbtv, nsc, nst, &
                     u, tau, d0, w, ierr)
        call run_exact_fixed_k(G, N, rc, rt, K, pe, oe, le, nwe, nbce, nbte, ierre)
        call assert_equal_int(ierr, ERR_OK, "parity: variance ierr")
        call assert_equal_int(ierre, ERR_OK, "parity: exact ierr")
        call assert_true(nwv > 0, "parity: no p-values")
        call assert_equal_int(nwv, nwe, "parity: same number of tested genes")
        call assert_equal_array_real(pv, pe, G, EXACT, "parity: p-values identical to exact fixed-k")
        call assert_equal_array_real(ov, oe, G, EXACT, "parity: statistic identical")
        call assert_equal_array_real(lv, le, G, EXACT, "parity: log2fc identical")
        call assert_equal_array_int(nbcv, nbce, G, "parity: case pool sizes identical")
        call assert_equal_array_int(nbtv, nbte, G, "parity: control pool sizes identical")
        ! all accepted -> exactly K slots scanned per side
        call assert_true(all(nsc == K) .and. all(nst == K), "parity: an all-accepting band scans exactly k slots")
        ! u is the detrended log-variance: centred near 0 and not constant
        call assert_true(maxval(u) - minval(u) > 0.1_real64, "parity: fixture must be heterogeneous in u")
    end subroutine test_mode_a_parity_with_exact

    ! =========================================================================
    ! 4. degradation: mode B without heterogeneity == exact at fixed k
    ! =========================================================================
    subroutine test_mode_b_degrades_to_exact()
        integer(int32), parameter :: G = 200, N = 4, K = 30
        real(real64) :: rc(N, G), rt(N, G)
        real(real64) :: pv(G), ov(G), lv(G), u(G), tau, d0, w
        real(real64) :: pe(G), oe(G), le(G)
        integer(int32) :: nwv, nbcv(G), nbtv(G), nsc(G), nst(G), ierr
        integer(int32) :: nwe, nbce(G), nbte(G), ierre

        call make_homog(G, N, 0.0_real64, 0.4_real64, rc)
        call make_homog(G, N, 1.1_real64, 0.0_real64, rt)

        call run_var(G, N, rc, rt, VARIANCE_MODE_SHRUNK, K, 1.0_real64, pv, ov, lv, nwv, nbcv, nbtv, nsc, nst, &
                     u, tau, d0, w, ierr)
        call assert_equal_int(ierr, ERR_OK, "degrade: ierr")
        call assert_true(.not. ieee_is_finite(d0), "degrade: d0 must be Inf without heterogeneity")
        call assert_equal_real(w, 0.0_real64, EXACT, "degrade: w = 0")
        call assert_equal_real(tau, 1.0_real64, EXACT, "degrade: tau_var = 1 (accept everything)")
        call assert_true(all(u == 0.0_real64), "degrade: u identically 0")

        call run_exact_fixed_k(G, N, rc, rt, K, pe, oe, le, nwe, nbce, nbte, ierre)
        call assert_equal_int(ierre, ERR_OK, "degrade: exact ierr")
        call assert_equal_array_real(pv, pe, G, EXACT, "degrade: mode B reproduces plain fixed-k TOX exactly")
        call assert_equal_array_int(nbcv, nbce, G, "degrade: case pool sizes identical")
    end subroutine test_mode_b_degrades_to_exact

    ! =========================================================================
    ! 5. a band nothing passes but the target itself
    ! =========================================================================
    subroutine test_strict_band_drops_genes()
        integer(int32), parameter :: G = 200, N = 4, K = 30
        real(real64) :: rc(N, G), rt(N, G)
        real(real64) :: pv(G), ov(G), lv(G), u(G), tau, d0, w
        integer(int32) :: nwv, nbcv(G), nbtv(G), nsc(G), nst(G), ierr

        call make_hetero(G, N, 0.0_real64, 30.0_real64, rc)
        call make_hetero(G, N, 1.7_real64, 0.0_real64, rt)
        call run_var(G, N, rc, rt, VARIANCE_MODE_BAND, K, 1.0e-14_real64, pv, ov, lv, nwv, nbcv, nbtv, nsc, nst, &
                     u, tau, d0, w, ierr)
        call assert_equal_int(ierr, ERR_OK, "strict: ierr")
        call assert_equal_int(nwv, 0, "strict: no pool reaches BAND_MIN_GENES, every gene dropped")
        call assert_true(all(pv == -1.0_real64), "strict: every p is -1")
        ! the case walk runs to scan_max (20k = 600 > G) i.e. exhausts all G slots;
        ! the control side is never gathered once the case side failed
        call assert_true(all(nsc == G), "strict: the case walk scans every slot")
        call assert_true(all(nst == -1), "strict: control side not gathered after the case side failed")
    end subroutine test_strict_band_drops_genes

    ! =========================================================================
    ! 6. validation
    ! =========================================================================
    subroutine test_variance_validation()
        integer(int32), parameter :: G = 50, N = 3
        real(real64) :: rc(N, G), rt(N, G), rbad(G, N)
        real(real64) :: pv(G), ov(G), lv(G), u(G), tau, d0, w, dh, sc
        integer(int32) :: nwv, nbcv(G), nbtv(G), nsc(G), nst(G), ierr, valid(G), nf

        call make_hetero(G, N, 0.0_real64, 0.0_real64, rc)
        call make_hetero(G, N, 1.0_real64, 0.0_real64, rt)
        valid = 1

        call run_var(G, N, rc, rt, 3, 10, 1.0_real64, pv, ov, lv, nwv, nbcv, nbtv, nsc, nst, u, tau, d0, w, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "validation: unknown variance_mode")
        call run_var(G, N, rc, rt, VARIANCE_MODE_BAND, 10, 0.0_real64, pv, ov, lv, nwv, nbcv, nbtv, nsc, nst, &
                     u, tau, d0, w, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "validation: c_tau must be > 0")
        call run_var(G, N, rc, rt, VARIANCE_MODE_BAND, 0, 1.0_real64, pv, ov, lv, nwv, nbcv, nbtv, nsc, nst, &
                     u, tau, d0, w, ierr)
        call assert_true(ierr /= ERR_OK, "validation: k = 0 rejected")

        call compute_noise_pvalue_pipeline_variance(rc, G, N, rt, G, N, valid, BETA_CENTRE_NONE, &
            G, 1, VARIANCE_MODE_BAND, 10, 1.0_real64, 0.3_real64, 1, 2, 1, 5000, &
            pv, ov, lv, dh, nwv, nbcv, nbtv, nsc, nst, u, tau, d0, w, sc, nf, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "validation: loess_robust must be 0/1")

        rbad = transpose(rt)
        call compute_noise_pvalue_pipeline_variance(rc, G, N, rbad, N, G, valid, BETA_CENTRE_NONE, &
            G, 1, VARIANCE_MODE_BAND, 10, 1.0_real64, 0.3_real64, 1, 1, 1, 5000, &
            pv, ov, lv, dh, nwv, nbcv, nbtv, nsc, nst, u, tau, d0, w, sc, nf, ierr)
        call assert_equal_int(ierr, ERR_DIM_MISMATCH, "validation: transposed control matrix")
    end subroutine test_variance_validation

end module mod_test_noise_model_variance
