module tox_ejscomptest
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    private

    integer(int32), parameter, public :: EJS_OK = 0_int32
    integer(int32), parameter, public :: EJS_ERR_INVALID_INPUT = 1_int32
    integer(int32), parameter, public :: EJS_ERR_ALLOC = 2_int32

    type, private :: study_values_t
        real(real64), allocatable :: v(:)
    end type study_values_t

    type, private :: bin_spec_t
        logical :: use_log_mode = .true.
        logical :: has_positive = .false.
        integer(int32) :: n_bins = 0_int32
        real(real64) :: xmin = 0.0_real64
        real(real64) :: xmax = 0.0_real64
        real(real64) :: delta = 0.0_real64
        real(real64) :: log_xmin = 0.0_real64
        real(real64) :: inv_log_rho = 0.0_real64
    end type bin_spec_t

    public :: ejscomptest_run

contains

    subroutine ejscomptest_run(values, study_offsets, n_studies, use_log_mode, n_bins, n_boot, n_perm, &
                               jsd_obs, ci_low, ci_high, p_value, n_i, ierr, &
                               min_value_keep, random_seed_value, boot_jsd_samples, perm_jsd_samples)
        real(real64), intent(in) :: values(:)
        integer(int32), intent(in) :: study_offsets(:)
        integer(int32), intent(in) :: n_studies
        logical, intent(in) :: use_log_mode
        integer(int32), intent(in) :: n_bins
        integer(int32), intent(in) :: n_boot
        integer(int32), intent(in) :: n_perm
        real(real64), intent(out) :: jsd_obs(n_studies)
        real(real64), intent(out) :: ci_low(n_studies)
        real(real64), intent(out) :: ci_high(n_studies)
        real(real64), intent(out) :: p_value(n_studies)
        integer(int32), intent(out) :: n_i(n_studies)
        integer(int32), intent(out) :: ierr
        real(real64), intent(in), optional :: min_value_keep
        integer(int32), intent(in), optional :: random_seed_value
        real(real64), intent(out), optional :: boot_jsd_samples(n_studies, n_boot)
        real(real64), intent(out), optional :: perm_jsd_samples(n_studies, n_perm)

        type(study_values_t), allocatable :: studies(:)
        real(real64), allocatable :: mix(:)
        real(real64), allocatable :: boot_work(:)
        real(real64), allocatable :: jsd_boot(:)
        real(real64), allocatable :: mix_boot(:)
        real(real64), allocatable :: perm_pool(:)
        real(real64), allocatable :: jsd_tmp(:)
        real(real64), allocatable :: perm_counter(:)
        real(real64), allocatable :: p_mix(:)
        integer(int32), allocatable :: idx_perm(:)
        integer(int32), allocatable :: n_i_perm(:)
        integer(int32) :: i, b, istat
        real(real64) :: min_keep
        type(bin_spec_t) :: spec

        ierr = EJS_OK

        if (.not. validate_inputs(values, study_offsets, n_studies, n_bins, n_boot, n_perm)) then
            ierr = EJS_ERR_INVALID_INPUT
            return
        end if

        min_keep = -1.0_real64
        if (present(min_value_keep)) min_keep = min_value_keep

        if (present(random_seed_value)) call init_random_seed(random_seed_value)

        allocate(studies(n_studies), stat=istat)
        if (istat /= 0) then
            ierr = EJS_ERR_ALLOC
            return
        end if

        call build_studies(values, study_offsets, n_studies, use_log_mode, min_keep, studies, mix, n_i, ierr)
        if (ierr /= EJS_OK) return

        if (size(mix) == 0) then
            jsd_obs = 0.0_real64
            ci_low = 0.0_real64
            ci_high = 0.0_real64
            p_value = 1.0_real64
            return
        end if

        call compute_jsd_all(studies, mix, use_log_mode, n_bins, jsd_obs, ierr)
        if (ierr /= EJS_OK) return

        allocate(jsd_tmp(n_studies), stat=istat)
        if (istat /= 0) then
            ierr = EJS_ERR_ALLOC
            return
        end if

        if (n_boot > 0) then
            allocate(mix_boot(size(mix)), jsd_boot(n_boot), stat=istat)
            if (istat /= 0) then
                ierr = EJS_ERR_ALLOC
                return
            end if

            do i = 1, n_studies
                if (size(studies(i)%v) == 0) then
                    ci_low(i) = 0.0_real64
                    ci_high(i) = 0.0_real64
                    if (present(boot_jsd_samples)) boot_jsd_samples(i, :) = 0.0_real64
                    cycle
                end if

                allocate(boot_work(size(studies(i)%v)), stat=istat)
                if (istat /= 0) then
                    ierr = EJS_ERR_ALLOC
                    return
                end if

                do b = 1, n_boot
                    call bootstrap_sample(studies(i)%v, boot_work)
                    call bootstrap_sample(mix, mix_boot)
                    call compute_jsd_single(boot_work, mix_boot, use_log_mode, n_bins, jsd_tmp(i), ierr)
                    if (ierr /= EJS_OK) return
                    jsd_boot(b) = jsd_tmp(i)
                    if (present(boot_jsd_samples)) boot_jsd_samples(i, b) = jsd_tmp(i)
                end do

                ci_low(i) = quantile_linear(jsd_boot, 0.025_real64)
                ci_high(i) = quantile_linear(jsd_boot, 0.975_real64)
                deallocate(boot_work)
            end do
        else
            ci_low = jsd_obs
            ci_high = jsd_obs
        end if

        if (n_perm > 0) then
            allocate(perm_counter(n_studies), perm_pool(size(mix)), idx_perm(size(mix)), &
                     n_i_perm(n_studies), stat=istat)
            if (istat /= 0) then
                ierr = EJS_ERR_ALLOC
                return
            end if

            perm_counter = 0.0_real64
            perm_pool = mix
            n_i_perm = n_i

            call build_bins_and_mix_pmf(mix, use_log_mode, n_bins, spec, p_mix, ierr)
            if (ierr /= EJS_OK) return

            do b = 1, n_perm
                call random_permutation(size(perm_pool), idx_perm)
                call evaluate_permutation(idx_perm, perm_pool, n_i_perm, n_studies, spec, p_mix, jsd_tmp, ierr)
                if (ierr /= EJS_OK) return

                do i = 1, n_studies
                    if (jsd_tmp(i) >= jsd_obs(i)) perm_counter(i) = perm_counter(i) + 1.0_real64
                    if (present(perm_jsd_samples)) perm_jsd_samples(i, b) = jsd_tmp(i)
                end do
            end do

            do i = 1, n_studies
                p_value(i) = (1.0_real64 + perm_counter(i)) / real(n_perm + 1_int32, real64)
            end do
        else
            p_value = 1.0_real64
        end if
    end subroutine ejscomptest_run

    logical function validate_inputs(values, study_offsets, n_studies, n_bins, n_boot, n_perm)
        real(real64), intent(in) :: values(:)
        integer(int32), intent(in) :: study_offsets(:)
        integer(int32), intent(in) :: n_studies, n_bins, n_boot, n_perm
        integer(int32) :: i

        validate_inputs = .false.

        if (n_studies <= 0 .or. n_bins <= 0 .or. n_boot < 0 .or. n_perm < 0) return
        if (size(study_offsets) /= n_studies + 1) return
        if (study_offsets(1) /= 1_int32) return
        if (study_offsets(n_studies + 1) /= size(values) + 1) return

        do i = 1, n_studies
            if (study_offsets(i) > study_offsets(i + 1)) return
        end do

        validate_inputs = .true.
    end function validate_inputs

    subroutine build_studies(values, study_offsets, n_studies, use_log_mode, min_keep, studies, mix, n_i, ierr)
        real(real64), intent(in) :: values(:)
        integer(int32), intent(in) :: study_offsets(:)
        integer(int32), intent(in) :: n_studies
        logical, intent(in) :: use_log_mode
        real(real64), intent(in) :: min_keep
        type(study_values_t), intent(inout) :: studies(:)
        real(real64), allocatable, intent(out) :: mix(:)
        integer(int32), intent(out) :: n_i(:)
        integer(int32), intent(out) :: ierr

        integer(int32) :: i, j, lo, hi, c, total_n, istat

        ierr = EJS_OK
        total_n = 0_int32

        do i = 1, n_studies
            lo = study_offsets(i)
            hi = study_offsets(i + 1) - 1
            c = 0_int32
            do j = lo, hi
                if (ieee_is_nan(values(j))) cycle
                if (min_keep >= 0.0_real64) then
                    if (values(j) < min_keep) cycle
                end if
                c = c + 1
            end do

            n_i(i) = c
            total_n = total_n + c
            allocate(studies(i)%v(c), stat=istat)
            if (istat /= 0) then
                ierr = EJS_ERR_ALLOC
                return
            end if

            c = 0_int32
            do j = lo, hi
                if (ieee_is_nan(values(j))) cycle
                if (min_keep >= 0.0_real64) then
                    if (values(j) < min_keep) cycle
                end if
                c = c + 1
                if (use_log_mode) then
                    studies(i)%v(c) = log(1.0_real64 + max(values(j), 0.0_real64))
                else
                    studies(i)%v(c) = max(values(j), 0.0_real64)
                end if
            end do
        end do

        allocate(mix(total_n), stat=istat)
        if (istat /= 0) then
            ierr = EJS_ERR_ALLOC
            return
        end if

        c = 0_int32
        do i = 1, n_studies
            if (size(studies(i)%v) > 0) then
                mix(c + 1:c + size(studies(i)%v)) = studies(i)%v
                c = c + size(studies(i)%v)
            end if
        end do
    end subroutine build_studies

    subroutine compute_jsd_all(studies, mix, use_log_mode, n_bins, jsd_obs, ierr)
        type(study_values_t), intent(in) :: studies(:)
        real(real64), intent(in) :: mix(:)
        logical, intent(in) :: use_log_mode
        integer(int32), intent(in) :: n_bins
        real(real64), intent(out) :: jsd_obs(:)
        integer(int32), intent(out) :: ierr

        type(bin_spec_t) :: spec
        real(real64), allocatable :: p_mix(:), p_study(:)
        integer(int32) :: i, istat

        ierr = EJS_OK
        call build_bins_and_mix_pmf(mix, use_log_mode, n_bins, spec, p_mix, ierr)
        if (ierr /= EJS_OK) return

        allocate(p_study(n_bins), stat=istat)
        if (istat /= 0) then
            ierr = EJS_ERR_ALLOC
            return
        end if

        do i = 1, size(studies)
            call histogram_pmf(studies(i)%v, spec, p_study)
            jsd_obs(i) = jsd_from_pmfs(p_study, p_mix)
        end do
    end subroutine compute_jsd_all

    subroutine compute_jsd_single(study_vals, mix_vals, use_log_mode, n_bins, jsd, ierr)
        real(real64), intent(in) :: study_vals(:)
        real(real64), intent(in) :: mix_vals(:)
        logical, intent(in) :: use_log_mode
        integer(int32), intent(in) :: n_bins
        real(real64), intent(out) :: jsd
        integer(int32), intent(out) :: ierr

        type(bin_spec_t) :: spec
        real(real64), allocatable :: p_mix(:), p_study(:)
        integer(int32) :: istat

        ierr = EJS_OK
        call build_bins_and_mix_pmf(mix_vals, use_log_mode, n_bins, spec, p_mix, ierr)
        if (ierr /= EJS_OK) return

        allocate(p_study(n_bins), stat=istat)
        if (istat /= 0) then
            ierr = EJS_ERR_ALLOC
            return
        end if

        call histogram_pmf(study_vals, spec, p_study)
        jsd = jsd_from_pmfs(p_study, p_mix)
    end subroutine compute_jsd_single

    subroutine build_bins_and_mix_pmf(mix_vals, use_log_mode, n_bins, spec, p_mix, ierr)
        real(real64), intent(in) :: mix_vals(:)
        logical, intent(in) :: use_log_mode
        integer(int32), intent(in) :: n_bins
        type(bin_spec_t), intent(out) :: spec
        real(real64), allocatable, intent(out) :: p_mix(:)
        integer(int32), intent(out) :: ierr

        integer(int32) :: n_pos, istat

        ierr = EJS_OK
        spec%use_log_mode = use_log_mode
        spec%n_bins = n_bins

        allocate(p_mix(n_bins), stat=istat)
        if (istat /= 0) then
            ierr = EJS_ERR_ALLOC
            return
        end if

        if (size(mix_vals) == 0) then
            p_mix = 0.0_real64
            return
        end if

        if (use_log_mode) then
            spec%xmin = minval(mix_vals)
            spec%xmax = maxval(mix_vals)
            spec%delta = (spec%xmax - spec%xmin) / real(n_bins, real64)
            spec%has_positive = .false.
        else
            n_pos = count(mix_vals > 0.0_real64)
            if (n_pos > 0) then
                spec%has_positive = .true.
                spec%xmin = minval(mix_vals, mask=mix_vals > 0.0_real64)
                spec%xmax = maxval(mix_vals, mask=mix_vals > 0.0_real64)
                if (spec%xmax > spec%xmin) then
                    spec%log_xmin = log(spec%xmin)
                    spec%inv_log_rho = real(n_bins - 1_int32, real64) / log(spec%xmax / spec%xmin)
                else
                    spec%log_xmin = log(spec%xmin)
                    spec%inv_log_rho = 0.0_real64
                end if
            else
                spec%has_positive = .false.
                spec%xmin = 0.0_real64
                spec%xmax = 0.0_real64
            end if
        end if

        call histogram_pmf(mix_vals, spec, p_mix)

        if (sum(p_mix) > 0.0_real64) then
            p_mix = p_mix / sum(p_mix)
        end if
    end subroutine build_bins_and_mix_pmf

    subroutine histogram_pmf(x, spec, p)
        real(real64), intent(in) :: x(:)
        type(bin_spec_t), intent(in) :: spec
        real(real64), intent(out) :: p(:)

        integer(int32) :: i, b

        p = 0.0_real64
        if (size(x) == 0) return

        do i = 1, size(x)
            b = locate_bin(x(i), spec)
            p(b) = p(b) + 1.0_real64
        end do

        p = p / real(size(x), real64)
    end subroutine histogram_pmf

    integer(int32) function locate_bin(x, spec)
        real(real64), intent(in) :: x
        type(bin_spec_t), intent(in) :: spec

        real(real64) :: rel
        integer(int32) :: idx_pos

        if (spec%use_log_mode) then
            if (spec%delta <= 0.0_real64) then
                locate_bin = spec%n_bins
                return
            end if
            rel = (x - spec%xmin) / spec%delta
            locate_bin = min(spec%n_bins, max(1_int32, int(rel) + 1_int32))
            return
        end if

        if (x <= 0.0_real64) then
            locate_bin = 1_int32
            return
        end if

        if (.not. spec%has_positive) then
            locate_bin = 1_int32
            return
        end if

        if (spec%inv_log_rho <= 0.0_real64) then
            locate_bin = 2_int32
            return
        end if

        idx_pos = int((log(x) - spec%log_xmin) * spec%inv_log_rho) + 1_int32
        idx_pos = min(spec%n_bins - 1_int32, max(1_int32, idx_pos))
        locate_bin = idx_pos + 1_int32
    end function locate_bin

    real(real64) function jsd_from_pmfs(p, q)
        real(real64), intent(in) :: p(:)
        real(real64), intent(in) :: q(:)

        real(real64) :: m, kl_pm, kl_qm
        integer(int32) :: i

        kl_pm = 0.0_real64
        kl_qm = 0.0_real64

        do i = 1, size(p)
            m = 0.5_real64 * (p(i) + q(i))
            if (p(i) > 0.0_real64 .and. m > 0.0_real64) kl_pm = kl_pm + p(i) * log(p(i) / m)
            if (q(i) > 0.0_real64 .and. m > 0.0_real64) kl_qm = kl_qm + q(i) * log(q(i) / m)
        end do

        jsd_from_pmfs = 0.5_real64 * (kl_pm + kl_qm)
    end function jsd_from_pmfs

    subroutine bootstrap_sample(pool, out)
        real(real64), intent(in) :: pool(:)
        real(real64), intent(out) :: out(:)

        integer(int32) :: i, idx
        real(real64) :: u

        do i = 1, size(out)
            call random_number(u)
            idx = min(size(pool), int(u * real(size(pool), real64)) + 1_int32)
            out(i) = pool(idx)
        end do
    end subroutine bootstrap_sample

    subroutine random_permutation(n, perm)
        integer(int32), intent(in) :: n
        integer(int32), intent(out) :: perm(n)

        integer(int32) :: i, j, tmp
        real(real64) :: u

        do i = 1, n
            perm(i) = i
        end do

        do i = n, 2, -1
            call random_number(u)
            j = int(u * real(i, real64)) + 1_int32
            tmp = perm(i)
            perm(i) = perm(j)
            perm(j) = tmp
        end do
    end subroutine random_permutation

    subroutine evaluate_permutation(perm, pool, n_i, n_studies, spec, p_mix, jsd_out, ierr)
        integer(int32), intent(in) :: perm(:)
        real(real64), intent(in) :: pool(:)
        integer(int32), intent(in) :: n_i(:)
        integer(int32), intent(in) :: n_studies
        type(bin_spec_t), intent(in) :: spec
        real(real64), intent(in) :: p_mix(:)
        real(real64), intent(out) :: jsd_out(n_studies)
        integer(int32), intent(out) :: ierr

        integer(int32) :: i, lo, hi, istat
        real(real64), allocatable :: tmp(:), p_study(:)

        ierr = EJS_OK
        allocate(p_study(spec%n_bins), stat=istat)
        if (istat /= 0) then
            ierr = EJS_ERR_ALLOC
            return
        end if

        lo = 1_int32
        do i = 1, n_studies
            hi = lo + n_i(i) - 1_int32
            if (n_i(i) <= 0) then
                jsd_out(i) = 0.0_real64
            else
                allocate(tmp(n_i(i)), stat=istat)
                if (istat /= 0) then
                    ierr = EJS_ERR_ALLOC
                    return
                end if
                tmp = pool(perm(lo:hi))
                call histogram_pmf(tmp, spec, p_study)
                jsd_out(i) = jsd_from_pmfs(p_study, p_mix)
                deallocate(tmp)
            end if
            lo = hi + 1_int32
        end do
    end subroutine evaluate_permutation

    real(real64) function quantile_linear(x, q)
        real(real64), intent(in) :: x(:)
        real(real64), intent(in) :: q

        real(real64), allocatable :: tmp(:)
        real(real64) :: pos, frac
        integer(int32) :: n, lo, hi, istat

        if (size(x) == 0) then
            quantile_linear = 0.0_real64
            return
        end if

        n = size(x)
        allocate(tmp(n), stat=istat)
        if (istat /= 0) then
            quantile_linear = 0.0_real64
            return
        end if
        tmp = x
        call sort_real(tmp)

        pos = q * real(n - 1_int32, real64) + 1.0_real64
        lo = max(1_int32, min(n, int(floor(pos))))
        hi = max(1_int32, min(n, lo + 1_int32))
        frac = pos - real(lo, real64)

        quantile_linear = (1.0_real64 - frac) * tmp(lo) + frac * tmp(hi)
    end function quantile_linear

    subroutine sort_real(a)
        real(real64), intent(inout) :: a(:)
        call quicksort_real(a, 1_int32, size(a, kind=int32))
    end subroutine sort_real

    recursive subroutine quicksort_real(a, left, right)
        real(real64), intent(inout) :: a(:)
        integer(int32), intent(in) :: left, right

        integer(int32) :: i, j
        real(real64) :: pivot, t

        if (left >= right) return

        i = left
        j = right
        pivot = a((left + right) / 2_int32)

        do
            do while (a(i) < pivot)
                i = i + 1_int32
            end do
            do while (a(j) > pivot)
                j = j - 1_int32
            end do
            if (i <= j) then
                t = a(i)
                a(i) = a(j)
                a(j) = t
                i = i + 1_int32
                j = j - 1_int32
            end if
            if (i > j) exit
        end do

        if (left < j) call quicksort_real(a, left, j)
        if (i < right) call quicksort_real(a, i, right)
    end subroutine quicksort_real

    subroutine init_random_seed(seed_value)
        integer(int32), intent(in) :: seed_value

        integer(int32) :: n, i
        integer(int32), allocatable :: seed(:)

        call random_seed(size=n)
        allocate(seed(n))

        do i = 1, n
            seed(i) = modulo(seed_value + 7919_int32 * i, huge(1_int32) - 1_int32)
            if (seed(i) <= 0) seed(i) = i
        end do
        call random_seed(put=seed)
    end subroutine init_random_seed

end module tox_ejscomptest
