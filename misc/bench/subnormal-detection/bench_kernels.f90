!> `scaled_length`'s two loops (#220), seven ways: today's, the three spellings called out of line
!| (what fpm builds: no cross-file inlining without -flto/-ipo), and the three spellings written
!| inline, which separates the cost of the test from the cost of the call.
module bench_kernels
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_class, ieee_class_type, ieee_positive_denormal, &
                                             ieee_negative_denormal, operator(==)
    use bench_helpers, only: saz_exponent, saz_class, saz_bits
    implicit none
    private
    public :: length_plain, length_exponent_call, length_class_call, length_exponent_inline, length_class_inline
    public :: length_bits_call, length_bits_inline
    integer(int64), parameter :: MAGNITUDE_MASK = huge(0_int64)
    integer(int64), parameter :: SMALLEST_NORMAL_BITS = transfer(tiny(1.0_real64), 0_int64)
    real(real64), parameter :: LARGEST_UNSCALED = 2.0_real64**500, SMALLEST_UNSCALED = 2.0_real64**(-500)
    real(real64), parameter :: SCALE_DOWN = 2.0_real64**(-600), SCALE_UP = 2.0_real64**600
contains
    pure subroutine pick_factor(largest, factor)
        real(real64), intent(in) :: largest
        real(real64), intent(out) :: factor
        factor = 1.0_real64
        if (largest > LARGEST_UNSCALED) factor = SCALE_DOWN
        if (largest < SMALLEST_UNSCALED) factor = SCALE_UP
    end subroutine pick_factor

    pure subroutine length_plain(n_dims, vector, scaled_norm)
        integer(int32), intent(in) :: n_dims
        real(real64), intent(in) :: vector(n_dims)
        real(real64), intent(out) :: scaled_norm
        real(real64) :: largest, factor
        integer(int32) :: i_dim
        largest = 0.0_real64
        do i_dim = 1, n_dims
            largest = max(largest, abs(vector(i_dim)))
        end do
        scaled_norm = 0.0_real64
        if (largest < tiny(1.0_real64)) return
        call pick_factor(largest, factor)
        do concurrent (i_dim = 1:n_dims) shared(vector, factor) reduce(+:scaled_norm)
            scaled_norm = scaled_norm + (vector(i_dim)*factor)**2
        end do
        scaled_norm = sqrt(scaled_norm)
    end subroutine length_plain

    pure subroutine length_exponent_call(n_dims, vector, scaled_norm)
        integer(int32), intent(in) :: n_dims
        real(real64), intent(in) :: vector(n_dims)
        real(real64), intent(out) :: scaled_norm
        real(real64) :: largest, factor
        integer(int32) :: i_dim
        largest = 0.0_real64
        do i_dim = 1, n_dims
            largest = max(largest, abs(saz_exponent(vector(i_dim))))
        end do
        scaled_norm = 0.0_real64
        if (largest <= 0.0_real64) return
        call pick_factor(largest, factor)
        do concurrent (i_dim = 1:n_dims) shared(vector, factor) reduce(+:scaled_norm)
            scaled_norm = scaled_norm + (saz_exponent(vector(i_dim))*factor)**2
        end do
        scaled_norm = sqrt(scaled_norm)
    end subroutine length_exponent_call

    pure subroutine length_class_call(n_dims, vector, scaled_norm)
        integer(int32), intent(in) :: n_dims
        real(real64), intent(in) :: vector(n_dims)
        real(real64), intent(out) :: scaled_norm
        real(real64) :: largest, factor
        integer(int32) :: i_dim
        largest = 0.0_real64
        do i_dim = 1, n_dims
            largest = max(largest, abs(saz_class(vector(i_dim))))
        end do
        scaled_norm = 0.0_real64
        if (largest <= 0.0_real64) return
        call pick_factor(largest, factor)
        do concurrent (i_dim = 1:n_dims) shared(vector, factor) reduce(+:scaled_norm)
            scaled_norm = scaled_norm + (saz_class(vector(i_dim))*factor)**2
        end do
        scaled_norm = sqrt(scaled_norm)
    end subroutine length_class_call

    pure subroutine length_exponent_inline(n_dims, vector, scaled_norm)
        integer(int32), intent(in) :: n_dims
        real(real64), intent(in) :: vector(n_dims)
        real(real64), intent(out) :: scaled_norm
        real(real64) :: largest, factor
        integer(int32) :: i_dim
        largest = 0.0_real64
        do i_dim = 1, n_dims
            largest = max(largest, abs(merge(0.0_real64, vector(i_dim), &
                                             exponent(vector(i_dim)) < minexponent(1.0_real64))))
        end do
        scaled_norm = 0.0_real64
        if (largest <= 0.0_real64) return
        call pick_factor(largest, factor)
        do concurrent (i_dim = 1:n_dims) shared(vector, factor) reduce(+:scaled_norm)
            scaled_norm = scaled_norm + (merge(0.0_real64, vector(i_dim), &
                                               exponent(vector(i_dim)) < minexponent(1.0_real64))*factor)**2
        end do
        scaled_norm = sqrt(scaled_norm)
    end subroutine length_exponent_inline

    pure subroutine length_class_inline(n_dims, vector, scaled_norm)
        integer(int32), intent(in) :: n_dims
        real(real64), intent(in) :: vector(n_dims)
        real(real64), intent(out) :: scaled_norm
        real(real64) :: largest, factor
        integer(int32) :: i_dim
        largest = 0.0_real64
        do i_dim = 1, n_dims
            largest = max(largest, abs(merge(0.0_real64, vector(i_dim), &
                ieee_class(vector(i_dim)) == ieee_positive_denormal .or. &
                ieee_class(vector(i_dim)) == ieee_negative_denormal)))
        end do
        scaled_norm = 0.0_real64
        if (largest <= 0.0_real64) return
        call pick_factor(largest, factor)
        do concurrent (i_dim = 1:n_dims) shared(vector, factor) reduce(+:scaled_norm)
            scaled_norm = scaled_norm + (merge(0.0_real64, vector(i_dim), &
                ieee_class(vector(i_dim)) == ieee_positive_denormal .or. &
                ieee_class(vector(i_dim)) == ieee_negative_denormal)*factor)**2
        end do
        scaled_norm = sqrt(scaled_norm)
    end subroutine length_class_inline
    pure subroutine length_bits_call(n_dims, vector, scaled_norm)
        integer(int32), intent(in) :: n_dims
        real(real64), intent(in) :: vector(n_dims)
        real(real64), intent(out) :: scaled_norm
        real(real64) :: largest, factor
        integer(int32) :: i_dim
        largest = 0.0_real64
        do i_dim = 1, n_dims
            largest = max(largest, abs(saz_bits(vector(i_dim))))
        end do
        scaled_norm = 0.0_real64
        if (largest <= 0.0_real64) return
        call pick_factor(largest, factor)
        do concurrent (i_dim = 1:n_dims) shared(vector, factor) reduce(+:scaled_norm)
            scaled_norm = scaled_norm + (saz_bits(vector(i_dim))*factor)**2
        end do
        scaled_norm = sqrt(scaled_norm)
    end subroutine length_bits_call

    pure subroutine length_bits_inline(n_dims, vector, scaled_norm)
        integer(int32), intent(in) :: n_dims
        real(real64), intent(in) :: vector(n_dims)
        real(real64), intent(out) :: scaled_norm
        real(real64) :: largest, factor, component
        integer(int64) :: magnitude
        integer(int32) :: i_dim
        largest = 0.0_real64
        do i_dim = 1, n_dims
            magnitude = iand(transfer(vector(i_dim), 0_int64), MAGNITUDE_MASK)
            component = merge(0.0_real64, vector(i_dim), magnitude > 0_int64 .and. magnitude < SMALLEST_NORMAL_BITS)
            largest = max(largest, abs(component))
        end do
        scaled_norm = 0.0_real64
        if (largest <= 0.0_real64) return
        call pick_factor(largest, factor)
        do concurrent (i_dim = 1:n_dims) shared(vector, factor) local(magnitude, component) reduce(+:scaled_norm)
            magnitude = iand(transfer(vector(i_dim), 0_int64), MAGNITUDE_MASK)
            component = merge(0.0_real64, vector(i_dim), magnitude > 0_int64 .and. magnitude < SMALLEST_NORMAL_BITS)
            scaled_norm = scaled_norm + (component*factor)**2
        end do
        scaled_norm = sqrt(scaled_norm)
    end subroutine length_bits_inline
end module bench_kernels
