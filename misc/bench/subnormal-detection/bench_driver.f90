!> Times each kernel over a matrix of ordinary vectors (no subnormals), best of 5 trials, at three
!| vector lengths. The pass count is calibrated per length so that a trial of today's kernel takes
!| at least 0.15 s, and then held fixed for all seven kernels.
!|
!| The clock is the process CPU time from `clock_gettime`, not `cpu_time`: nvfortran's `cpu_time`
!| is `gettimeofday` underneath, a wall clock that WSL2 steps backwards by about a second every
!| half minute -- the trap misc/bench/README.md describes for ifx's `system_clock`. On gfortran
!| and ifx the two clocks measure the same thing. For a threaded build it sums the threads.
program bench_driver
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: iso_c_binding, only: c_int, c_long
    use bench_kernels
    implicit none
    type, bind(c) :: timespec
        integer(c_long) :: seconds, nanoseconds
    end type timespec
    interface
        integer(c_int) function clock_gettime(clock_id, time) bind(c, name="clock_gettime")
            import :: c_int, timespec
            integer(c_int), value :: clock_id
            type(timespec), intent(out) :: time
        end function clock_gettime
    end interface
    integer(c_int), parameter :: CLOCK_PROCESS_CPUTIME_ID = 2       ! Linux
    integer(int32), parameter :: TOTAL = 2**18, N_TRIALS = 5, N_KERNELS = 7
    integer(int32), parameter :: DIMS(3) = [8, 64, 1024]
    real(real64), allocatable :: pristine(:, :), work(:, :)
    real(real64) :: best(N_KERNELS), checksum(N_KERNELS), t0, t1, elapsed
    integer(int32) :: i_size, n_dims, n_vecs, passes, i_kernel, i_trial, i_pass
    character(len=*), parameter :: fmt_row = '(i6, i8, 7f12.3, 3x, 6f7.0)'

    print '(a)', "  dims  passes       plain    exp-call    cls-call     exp-inl     cls-inl   bits-call    bits-inl" // &
                 "   overhead %: exp-c  cls-c  exp-i  cls-i  bit-c  bit-i"
    do i_size = 1, size(DIMS)
        n_dims = DIMS(i_size)
        n_vecs = TOTAL/n_dims
        allocate (pristine(n_dims, n_vecs), work(n_dims, n_vecs))
        call random_seed(put=[(42, i_kernel = 1, 64)])
        call random_number(pristine)
        pristine = 0.5_real64 + 1.5_real64*pristine

        passes = 1
        do
            work = pristine
            t0 = process_seconds()
            call run(1, passes, checksum(1))
            t1 = process_seconds()
            if (t1 - t0 >= 0.15_real64) exit
            passes = passes*2
        end do

        do i_kernel = 1, N_KERNELS
            best(i_kernel) = huge(1.0_real64)
            do i_trial = 1, N_TRIALS
                work = pristine
                checksum(i_kernel) = 0.0_real64
                t0 = process_seconds()
                call run(i_kernel, passes, checksum(i_kernel))
                t1 = process_seconds()
                elapsed = t1 - t0
                if (elapsed > 0.0_real64) best(i_kernel) = min(best(i_kernel), elapsed)
            end do
        end do
        best = best/(real(passes, real64)*real(TOTAL, real64))*1.0e9_real64
        print fmt_row, n_dims, passes, best, 100*(best(2:) - best(1))/best(1)
        if (any(checksum /= checksum(1))) then
            print '(a, 7es24.16)', "  CHECKSUM MISMATCH: ", checksum
        end if
        deallocate (pristine, work)
    end do
contains
    real(real64) function process_seconds()
        type(timespec) :: now
        if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, now) /= 0) error stop "clock_gettime failed"
        process_seconds = real(now%seconds, real64) + 1.0e-9_real64*real(now%nanoseconds, real64)
    end function process_seconds

    subroutine run(i_kernel, passes, checksum)
        integer(int32), intent(in) :: i_kernel, passes
        real(real64), intent(inout) :: checksum
        real(real64) :: scaled_norm
        integer(int32) :: i_vec
        do i_pass = 1, passes
            work(1, 1) = work(1, 1) + 1.0e-12_real64          ! vary the input every pass
            do i_vec = 1, n_vecs
                select case (i_kernel)
                case (1); call length_plain(n_dims, work(:, i_vec), scaled_norm)
                case (2); call length_exponent_call(n_dims, work(:, i_vec), scaled_norm)
                case (3); call length_class_call(n_dims, work(:, i_vec), scaled_norm)
                case (4); call length_exponent_inline(n_dims, work(:, i_vec), scaled_norm)
                case (5); call length_class_inline(n_dims, work(:, i_vec), scaled_norm)
                case (6); call length_bits_call(n_dims, work(:, i_vec), scaled_norm)
                case (7); call length_bits_inline(n_dims, work(:, i_vec), scaled_norm)
                end select
                checksum = checksum + scaled_norm
            end do
        end do
    end subroutine run
end program bench_driver
