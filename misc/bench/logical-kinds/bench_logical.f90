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
! implementations, where every intrinsic yielding a default `logical` would then need
! converting.
!
! `logical(c_bool)` was adopted, so the marshalling columns are now history: no wrapper
! copies a mask any more. The live question is what the implementations pay for holding
! their masks in c_bool, and that is answered by comparing the two kinds column-wise --
! `mk:` for writing a mask, `rd:`/`br:` for the two ways an implementation reads one back.
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
    !: each measurement is repeated and the fastest kept -- the usual noise floor.
    integer(int32), parameter :: TRIALS = 3
    integer(int32) :: i
    integer(int64) :: sink = 0

    print '(a)', "# ns per element. 'in' and 'out' are the two marshalling directions,"
    print '(a)', "# 'auto' includes the automatic-array allocation the wrapper really pays,"
    print '(a)', "# 'mk' is computing a mask, 'rd' is reading one whole (count), 'br' is"
    print '(a)', "# reading one element at a time to guard work -- what the implementations do."
    print '(a)', ""
    print '(a10,9a12)', "n", "in", "out", "auto", "mk:logical", "mk:c_bool", &
        "rd:logical", "rd:c_bool", "br:logical", "br:c_bool"
    do i = 1, size(SIZES)
        call one_size(SIZES(i), sink)
    end do
    print '(a)', ""
    print '(a,i0)', "# checksum (ignore, defeats dead-code elimination): ", sink

contains

    !> `cpu_time`, not `system_clock`. ifx's `system_clock` is wall-clock based -- its count
    !| is the Unix epoch in microseconds, where gfortran's is monotonic uptime -- so an NTP
    !| step under WSL2 lands as a backwards jump mid-measurement. That does not merely add
    !| noise: it makes a trial look *faster*, and best-of-N then keeps precisely the corrupted
    !| trial. It showed as single cells an order of magnitude below their neighbours, in a
    !| different cell on every run. `cpu_time` is monotonic on both compilers, this benchmark
    !| is single-threaded so its CPU time is its wall time, and every measurement runs for
    !| 0.1 s or more -- far above the clock's resolution.
    real(real64) function seconds_now()
        call cpu_time(seconds_now)
    end function seconds_now

    real(real64) function elapsed_ns(t0, n, reps)
        real(real64), intent(in) :: t0
        integer(int32), intent(in) :: n, reps
        real(real64) :: t1

        call cpu_time(t1)
        elapsed_ns = (t1 - t0)*1.0e9_real64/(real(n, real64)*real(reps, real64))
    end function elapsed_ns

    subroutine one_size(n, sink)
        integer(int32), intent(in) :: n
        integer(int64), intent(inout) :: sink

        logical(c_bool), allocatable :: src(:), out_c(:)
        logical, allocatable :: dst(:), kern(:)
        real(real64), allocatable :: x(:)
        integer(int32) :: rep, reps, trial
        real(real64) :: t0
        real(real64) :: ns(9), this(9)

        reps = int(max(3_int64, WORK/int(n, int64)), int32)
        allocate (src(n), dst(n), kern(n), out_c(n), x(n))
        call random_number(x)
        src = x > 0.5_real64
        kern = x > 0.5_real64
        ns = huge(1.0_real64)

        do trial = 1, TRIALS
        t0 = seconds_now()
        do rep = 1, reps
            src(1 + mod(rep, n)) = .not. src(1 + mod(rep, n))   ! defeat hoisting
            call marshal_in(n, src, dst, sink)
        end do
        this(1) = elapsed_ns(t0, n, reps)

        t0 = seconds_now()
        do rep = 1, reps
            kern(1 + mod(rep, n)) = .not. kern(1 + mod(rep, n))
            call marshal_out(n, kern, out_c, sink)
        end do
        this(2) = elapsed_ns(t0, n, reps)

        t0 = seconds_now()
        do rep = 1, reps
            src(1 + mod(rep, n)) = .not. src(1 + mod(rep, n))
            call marshal_automatic(n, src, sink)
        end do
        this(3) = elapsed_ns(t0, n, reps)

        t0 = seconds_now()
        do rep = 1, reps
            call make_default(n, x, real(mod(rep, 100), real64)/100.0_real64, kern, sink)
        end do
        this(4) = elapsed_ns(t0, n, reps)

        t0 = seconds_now()
        do rep = 1, reps
            call make_cbool(n, x, real(mod(rep, 100), real64)/100.0_real64, src, sink)
        end do
        this(5) = elapsed_ns(t0, n, reps)

        t0 = seconds_now()
        do rep = 1, reps
            kern(1 + mod(rep, n)) = .not. kern(1 + mod(rep, n))
            call read_mask(n, kern, sink)
        end do
        this(6) = elapsed_ns(t0, n, reps)

        t0 = seconds_now()
        do rep = 1, reps
            src(1 + mod(rep, n)) = .not. src(1 + mod(rep, n))
            call read_mask_cbool(n, src, sink)
        end do
        this(7) = elapsed_ns(t0, n, reps)

        t0 = seconds_now()
        do rep = 1, reps
            kern(1 + mod(rep, n)) = .not. kern(1 + mod(rep, n))
            call branch_default(n, kern, x, sink)
        end do
        this(8) = elapsed_ns(t0, n, reps)

        t0 = seconds_now()
        do rep = 1, reps
            src(1 + mod(rep, n)) = .not. src(1 + mod(rep, n))
            call branch_cbool(n, src, x, sink)
        end do
        this(9) = elapsed_ns(t0, n, reps)

        where (this > 0.0_real64) ns = min(ns, this)
        end do

        where (ns == huge(1.0_real64)) ns = -1.0_real64   ! no trial measured above zero
        print '(i10,9f12.3)', n, ns
        deallocate (src, dst, kern, out_c, x)
    end subroutine one_size

end program bench_logical
