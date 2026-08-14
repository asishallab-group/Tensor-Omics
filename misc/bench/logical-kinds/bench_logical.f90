! What the c_bool marshalling copy costs.
!
! Every generated C wrapper that takes or returns a logical array declares an automatic
! array of default `logical` and copies elementwise, because `c_bool` is one byte and a
! default `logical` is four:
!
!     logical, dimension(n_vectors) :: vectors_selection_mask_f
!     vectors_selection_mask_f = vectors_selection_mask
!
! Declaring the implementation's dummy `logical(c_bool)` would remove both the temporary
! and the copy. This measures what that is worth, and what it would cost inside the
! kernels, where every intrinsic yielding a default `logical` would then need converting.
!
! Build with ./run_bench.sh, which does every configuration. By hand -- two units, and NOT
! with -flto / -ipo, or the copies are optimised away (see ../README.md):
!     gfortran -O0 -c bench_kernels.f90 && gfortran -O0 bench_kernels.o bench_logical.f90 -o bench
!     gfortran -O3 -march=native -mtune=native -funroll-loops -ftree-vectorize ...
!     ifx -O0 -heap-arrays ...          (-heap-arrays: see the note on marshal_automatic)
!     ifx -O3 -xHost -heap-arrays ...
program bench_logical
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: iso_c_binding, only: c_bool
    use bench_kernels
    implicit none

    integer(int32), parameter :: SIZES(4) = [10000, 100000, 1000000, 10000000]
    !: total element-visits per measurement, so every size takes roughly the same wall time
    integer(int64), parameter :: WORK = 200000000_int64
    !: each measurement is repeated and the fastest kept. Partly the usual noise floor, but
    !| also necessary: ifx's `system_clock` is wall-clock based (its count is the Unix epoch
    !| in microseconds, where gfortran's is monotonic uptime), so an NTP step lands as a
    !| backwards jump mid-measurement. Discarding non-positive results is what catches that.
    integer(int32), parameter :: TRIALS = 3
    integer(int32) :: i
    integer(int64) :: sink = 0

    print '(a)', "# ns per element. 'in' and 'out' are the two marshalling directions,"
    print '(a)', "# 'auto' includes the automatic-array allocation the wrapper really pays,"
    print '(a)', "# 'mk' is computing a mask, 'rd' is reading one (count)."
    print '(a)', ""
    print '(a10,6a12)', "n", "in", "out", "auto", "mk:logical", "mk:c_bool", "rd:logical"
    do i = 1, size(SIZES)
        call one_size(SIZES(i), sink)
    end do
    print '(a)', ""
    print '(a,i0)', "# checksum (ignore, defeats dead-code elimination): ", sink

contains

    !> Both endpoints are read the same way. `system_clock`'s tick rate depends on the kind
    !| of its COUNT argument, so mixing a one-argument and a three-argument call can silently
    !| pair ticks from one clock with the rate of another.
    integer(int64) function ticks_now()
        integer(int64) :: rate, cmax
        call system_clock(ticks_now, rate, cmax)
    end function ticks_now

    real(real64) function elapsed_ns(t0, n, reps)
        integer(int64), intent(in) :: t0
        integer(int32), intent(in) :: n, reps
        integer(int64) :: t1, rate, cmax

        call system_clock(t1, rate, cmax)
        ! No wraparound correction: an int64 counter at these rates runs for millennia, and
        ! a correction that adds `count_max` overflows int64 and reports 1e300 instead.
        elapsed_ns = real(t1 - t0, real64)*1.0e9_real64 &
                     /(real(rate, real64)*real(n, real64)*real(reps, real64))
    end function elapsed_ns

    subroutine one_size(n, sink)
        integer(int32), intent(in) :: n
        integer(int64), intent(inout) :: sink

        logical(c_bool), allocatable :: src(:), out_c(:)
        logical, allocatable :: dst(:), kern(:)
        real(real64), allocatable :: x(:)
        integer(int32) :: rep, reps, trial
        integer(int64) :: t0
        real(real64) :: ns(6), this(6)

        reps = int(max(3_int64, WORK/int(n, int64)), int32)
        allocate (src(n), dst(n), kern(n), out_c(n), x(n))
        call random_number(x)
        src = x > 0.5_real64
        kern = x > 0.5_real64
        ns = huge(1.0_real64)

        do trial = 1, TRIALS
        t0 = ticks_now()
        do rep = 1, reps
            src(1 + mod(rep, n)) = .not. src(1 + mod(rep, n))   ! defeat hoisting
            call marshal_in(n, src, dst, sink)
        end do
        this(1) = elapsed_ns(t0, n, reps)

        t0 = ticks_now()
        do rep = 1, reps
            kern(1 + mod(rep, n)) = .not. kern(1 + mod(rep, n))
            call marshal_out(n, kern, out_c, sink)
        end do
        this(2) = elapsed_ns(t0, n, reps)

        t0 = ticks_now()
        do rep = 1, reps
            src(1 + mod(rep, n)) = .not. src(1 + mod(rep, n))
            call marshal_automatic(n, src, sink)
        end do
        this(3) = elapsed_ns(t0, n, reps)

        t0 = ticks_now()
        do rep = 1, reps
            call make_default(n, x, real(mod(rep, 100), real64)/100.0_real64, kern, sink)
        end do
        this(4) = elapsed_ns(t0, n, reps)

        t0 = ticks_now()
        do rep = 1, reps
            call make_cbool(n, x, real(mod(rep, 100), real64)/100.0_real64, src, sink)
        end do
        this(5) = elapsed_ns(t0, n, reps)

        t0 = ticks_now()
        do rep = 1, reps
            kern(1 + mod(rep, n)) = .not. kern(1 + mod(rep, n))
            call read_mask(n, kern, sink)
        end do
        this(6) = elapsed_ns(t0, n, reps)

        where (this > 0.0_real64) ns = min(ns, this)
        end do

        where (ns == huge(1.0_real64)) ns = -1.0_real64   ! every trial was disturbed
        print '(i10,6f12.3)', n, ns
        deallocate (src, dst, kern, out_c, x)
    end subroutine one_size

end program bench_logical
