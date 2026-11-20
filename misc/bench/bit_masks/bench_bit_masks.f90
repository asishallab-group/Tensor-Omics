! Why f42_bit_masks is a core of procedures on plain int32 words and not a derived type.
!
! Every operation below is timed in each form the two designs offer (see bench_kernels.f90),
! with logical(c_bool) arrays -- today's masks -- as the reference. Results are ns per bit of
! the mask; the merge is ms per full merge.
!
! Build and run with ./run_bench.sh. By hand -- three units, and NOT with -flto / -ipo:
!     gfortran -O3 -c mask_candidates.f90 bench_kernels.f90
!     gfortran -O3 mask_candidates.o bench_kernels.o bench_bit_masks.f90 -o bench
! Under ifx, `pack:assignment` at 10^7 bits needs `ulimit -s unlimited`
! (probes/probe_assignment_stack.f90 shows why).
program bench_bit_masks
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use mask_candidates
    use bench_kernels
    implicit none

    integer(int32), parameter :: SIZES(3) = [1000_int32, 100000_int32, 10000000_int32]
    !: each measurement repeats its kernel until it runs this long, then keeps the best trial
    real(real64), parameter :: MIN_SECONDS = 0.1_real64
    integer(int32), parameter :: TRIALS = 3_int32
    integer(int32), parameter :: N_KERNELS = 14_int32
    character(len=20), parameter :: KERNEL_NAMES(N_KERNELS) = [character(len=20) :: &
        'and_count:flags', 'and_count:core', 'and_count:eager', 'and_count:lazy_bits', &
        'and_count:lazy_words', &
        'and:flags', 'and:core', 'and:eager', 'and:view', &
        'copy:flags', 'copy:core', 'copy:owned', &
        'pack:core', 'pack:assignment']
    !: Shatter's merge: every pair of MERGE_MASKS columns of MERGE_BITS bits
    integer(int32), parameter :: MERGE_BITS = 100000_int32, MERGE_MASKS = 100_int32

    logical(c_bool), allocatable :: flags_left(:), flags_right(:), flags_result(:)
    integer(int32), allocatable, target :: words_left(:), words_right(:), words_result(:)
    type(owned_mask), target :: owned_left, owned_right, owned_result
    type(view_mask) :: view_left, view_right, view_result
    integer(int32) :: n_bits, n_words, i_size, i_kernel
    integer(int64) :: sink = 0_int64
    real(real64) :: ns_per_bit(N_KERNELS, size(SIZES))

    do i_size = 1, size(SIZES, kind=int32)
        call set_up(SIZES(i_size))
        do i_kernel = 1, N_KERNELS
            ns_per_bit(i_kernel, i_size) = best_ns_per_bit(i_kernel)
        end do
    end do

    print '(a)', '# ns per bit, best of 3 trials of >= 0.1 s each'
    print '(a20,3i10)', 'n_bits', SIZES
    do i_kernel = 1, N_KERNELS
        print '(a20,3f10.4)', KERNEL_NAMES(i_kernel), ns_per_bit(i_kernel, :)
    end do
    call bench_merge()
    print '(a,i0)', '# checksum (ignore, defeats dead-code elimination): ', sink

contains

    !> Random masks of density 0.5 in every representation, all holding the same bits, and a
    !| check that every form of an operation agrees with the logical(c_bool) reference.
    subroutine set_up(n)
        integer(int32), intent(in) :: n
        real(real64), allocatable :: uniform(:)
        integer(int32) :: seed_size
        integer(int32), allocatable :: seed(:)

        call random_seed(size=seed_size)
        allocate (seed(seed_size), source=12345_int32)
        call random_seed(put=seed)
        n_bits = n
        n_words = words_for_bits(n_bits)
        if (allocated(flags_left)) deallocate (flags_left, flags_right, flags_result, words_left, &
                                               words_right, words_result)
        allocate (flags_left(n_bits), flags_right(n_bits), flags_result(n_bits), uniform(n_bits))
        allocate (words_left(n_words), words_right(n_words), words_result(n_words), source=0_int32)
        call random_number(uniform)
        flags_left = uniform < 0.5_real64
        call random_number(uniform)
        flags_right = uniform < 0.5_real64
        flags_result = .false._c_bool
        call core_pack(n_bits, n_words, flags_left, words_left)
        call core_pack(n_bits, n_words, flags_right, words_right)
        owned_left = flags_left
        owned_right = flags_right
        owned_result = flags_result
        view_left%words => words_left
        view_right%words => words_right
        view_result%words => words_result

        if (and_count_core(n_words, words_left, words_right) /= count(flags_left .and. flags_right) .or. &
            and_count_eager(owned_left, owned_right) /= count(flags_left .and. flags_right) .or. &
            and_count_lazy_bits(owned_left, owned_right) /= count(flags_left .and. flags_right) .or. &
            and_count_lazy_words(owned_left, owned_right) /= count(flags_left .and. flags_right) .or. &
            any(owned_left%words /= words_left)) then
            error stop 'the forms of an operation disagree'
        end if
    end subroutine set_up

    !> Flips one bit of every left operand, so no repetition can be hoisted out of the loop.
    subroutine perturb(i_rep)
        integer(int64), intent(in) :: i_rep
        integer(int32) :: i_bit, i_word

        i_bit = int(mod(i_rep, int(n_bits, int64)), int32) + 1_int32
        i_word = (i_bit - 1_int32)/BITS_PER_WORD + 1_int32
        flags_left(i_bit) = .not. flags_left(i_bit)
        words_left(i_word) = ieor(words_left(i_word), 1_int32)
        owned_left%words(i_word) = ieor(owned_left%words(i_word), 1_int32)
    end subroutine perturb

    subroutine run_kernel(kernel)
        integer(int32), intent(in) :: kernel

        select case (kernel)
        case (1); sink = sink + and_count_flags(n_bits, flags_left, flags_right)
        case (2); sink = sink + and_count_core(n_words, words_left, words_right)
        case (3); sink = sink + and_count_eager(owned_left, owned_right)
        case (4); sink = sink + and_count_lazy_bits(owned_left, owned_right)
        case (5); sink = sink + and_count_lazy_words(owned_left, owned_right)
        case (6); call and_flags(n_bits, flags_left, flags_right, flags_result)
        case (7); call and_core(n_words, words_left, words_right, words_result)
        case (8); call and_eager(owned_left, owned_right, owned_result)
        case (9); call and_view(view_left, view_right, view_result)
        case (10); call copy_flags(n_bits, flags_left, flags_result)
        case (11); call copy_core(n_words, words_left, words_result)
        case (12); call copy_owned(owned_left, owned_result)
        case (13); call pack_core(n_bits, n_words, flags_left, words_result)
        case (14); call pack_assignment(flags_left, owned_result)
        end select
        sink = sink + merge(1_int64, 0_int64, logical(flags_result(1))) + words_result(1) + owned_result%words(1)
    end subroutine run_kernel

    !> `cpu_time`, not `system_clock`: ifx's system_clock is wall-clock based and jumps
    !| backwards under WSL2 (see ../README.md).
    real(real64) function best_ns_per_bit(kernel)
        integer(int32), intent(in) :: kernel
        integer(int64) :: n_reps, i_rep
        integer(int32) :: i_trial
        real(real64) :: started, finished, elapsed, best

        call run_kernel(kernel)   ! warm-up
        n_reps = 1_int64
        do   ! calibrate: double the repetitions until one trial takes MIN_SECONDS
            call cpu_time(started)
            do i_rep = 1, n_reps
                call perturb(i_rep)
                call run_kernel(kernel)
            end do
            call cpu_time(finished)
            elapsed = finished - started
            if (elapsed >= MIN_SECONDS) exit
            n_reps = 2_int64*n_reps
        end do
        best = elapsed
        do i_trial = 2, TRIALS
            call cpu_time(started)
            do i_rep = 1, n_reps
                call perturb(i_rep)
                call run_kernel(kernel)
            end do
            call cpu_time(finished)
            best = min(best, finished - started)
        end do
        best_ns_per_bit = best*1.0e9_real64/(real(n_reps, real64)*real(n_bits, real64))
    end function best_ns_per_bit

    subroutine bench_merge()
        logical(c_bool), allocatable :: mask_flags(:, :)
        integer(int32), allocatable :: mask_words(:, :)
        real(real64), allocatable :: uniform(:, :)
        integer(int32) :: merge_words, i_mask, i_trial
        integer(int64) :: expected
        real(real64) :: started, finished, best_flags, best_core

        merge_words = words_for_bits(MERGE_BITS)
        allocate (mask_flags(MERGE_BITS, MERGE_MASKS), mask_words(merge_words, MERGE_MASKS), &
                  uniform(MERGE_BITS, MERGE_MASKS))
        call random_number(uniform)
        mask_flags = uniform < 0.3_real64
        do i_mask = 1, MERGE_MASKS
            call core_pack(MERGE_BITS, merge_words, mask_flags(:, i_mask), mask_words(:, i_mask))
        end do
        expected = merge_flags(MERGE_BITS, MERGE_MASKS, mask_flags)
        if (merge_core(merge_words, MERGE_MASKS, mask_words) /= expected) error stop 'merge disagrees'

        best_flags = huge(1.0_real64)
        best_core = huge(1.0_real64)
        do i_trial = 1, TRIALS
            call cpu_time(started)
            sink = sink + merge_flags(MERGE_BITS, MERGE_MASKS, mask_flags)
            call cpu_time(finished)
            best_flags = min(best_flags, finished - started)
            call cpu_time(started)
            sink = sink + merge_core(merge_words, MERGE_MASKS, mask_words)
            call cpu_time(finished)
            best_core = min(best_core, finished - started)
        end do
        print '(a,i0,a,i0,a)', '# merge: all pairs of ', MERGE_MASKS, ' masks of ', MERGE_BITS, &
            ' bits, ms per full merge'
        print '(a20,f10.1)', 'merge:flags', 1.0e3_real64*best_flags
        print '(a20,f10.1)', 'merge:core', 1.0e3_real64*best_core
    end subroutine bench_merge

end program bench_bit_masks
