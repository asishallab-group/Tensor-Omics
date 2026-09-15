#include <src/macros.h>

!> Random number generation backed by the GNU Scientific Library (GSL).
!| AUTHOR_LASZLO_LANG
!|
!| Wraps GSL's RNG (`gsl_rng`) and its hypergeometric/binomial samplers behind the opaque
!| [[f42_random_gsl(module):rng_t(type)]] handle. For selecting the RNG algorithm, use the
!| environment variable `GSL_RNG_TYPE`, as described in
!| [the GSL documentation](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-environment-variables).
!|
!| Hand-written, Path II infrastructure -- **not** an `_impl` module. `rng_t` carries a bare
!| `c_ptr` into GSL's own allocation with no `DM_*`-expressible contract (the generator refuses
!| a derived-type argument outright), so nothing here is `M_EXPORT_C` and nothing generates a
!| wrapper -- the same category as [[f42_safeguard(module)]] and [[f42_config(module)]]. Every
!| procedure therefore validates its own inputs by hand and reports failure through its own
!| `ierr`, using [[tox_errors(module)]] exactly as a generated wrapper would.
!|
!| Callers reach this only from other `_impl` modules (it sits on
!| `Conventions.impl_import_whitelist`) or from Fortran test code; it is never called from C,
!| Python or R directly.
module f42_random_gsl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_double, c_int, c_long
    use tox_errors, only: set_ok, set_err_once, is_err
    use tox_errors, only: validate_dimension_size, validate_in_range_int, validate_in_range_real
    use tox_errors, only: ERR_POINTER_NULL, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: rng_t, create_rng, destroy_rng, reset_rng
    public :: random_uniform, random_hypergeom, random_binomial, rand_range
    public :: random_multinomial, random_multiv_hypergeom

    interface gsl_get_rng_type
        !> Returns algorithm for the RNG, set by the environment variable `GSL_RNG_TYPE`
        function gsl_get_rng_type() result(rng_type) bind(C, name="gsl_rng_env_setup")
            import
            type(c_ptr) :: rng_type
        end function gsl_get_rng_type
    end interface gsl_get_rng_type

    interface gsl_alloc_rng
        !> Creates a RNG structure to generate random numbers with, uses seed set by environment variable `GSL_RNG_SEED`
        function gsl_alloc_rng(rng_type) result(rng) bind(C, name="gsl_rng_alloc")
            import
            type(c_ptr), intent(in), value :: rng_type
                !! RNG type returned by [[f42_random_gsl(module):gsl_get_rng_type(interface)]]
            type(c_ptr) :: rng
                !! Created RNG object, `c_null_ptr` on allocation failure
        end function gsl_alloc_rng
    end interface gsl_alloc_rng

    interface gsl_free_rng
        !> Deallocates all memory related to `rng`
        subroutine gsl_free_rng(rng) bind(C, name="gsl_rng_free")
            import
            type(c_ptr), intent(in), value :: rng
                !! RNG object, originally created by [[f42_random_gsl(module):gsl_alloc_rng(interface)]]
        end subroutine gsl_free_rng
    end interface gsl_free_rng

    interface gsl_set_rng_seed
        !> Sets the seed for a RNG structure to a specific one. Will reset the RNG stream if it is the same seed as already used.
        subroutine gsl_set_rng_seed(rng, seed) bind(C, name="gsl_rng_set")
            import
            type(c_ptr), intent(in), value :: rng
                !! RNG object, originally created by [[f42_random_gsl(module):gsl_alloc_rng(interface)]]
            integer(c_long), intent(in), value :: seed
                !! Seed to initialize `rng` with
        end subroutine gsl_set_rng_seed
    end interface gsl_set_rng_seed

    interface gsl_random_uniform
        !> Produces a uniform random number
        function gsl_random_uniform(rng) result(rand) bind(C, name="gsl_rng_uniform")
            import
            type(c_ptr), intent(in), value :: rng
                !! RNG object, originally created by [[f42_random_gsl(module):gsl_alloc_rng(interface)]]
            real(c_double) :: rand
                !! Random uniform number
        end function gsl_random_uniform
    end interface gsl_random_uniform

    interface gsl_random_hypergeom
        !> Performs a hypergeometric draw
        function gsl_random_hypergeom(rng, n_population1, n_population2, n_samples) result(n_drawn) bind(C, name="gsl_ran_hypergeometric")
            import
            type(c_ptr), intent(in), value :: rng
                !! RNG object, originally created by [[f42_random_gsl(module):gsl_alloc_rng(interface)]]
            integer(c_int), intent(in), value :: n_population1
                !! Size of success population
            integer(c_int), intent(in), value :: n_population2
                !! Size of failure population
            integer(c_int), intent(in), value :: n_samples
                !! Number of draws to perform
            integer(c_int) :: n_drawn
                !! Number of drawn elements from success population
        end function gsl_random_hypergeom
    end interface gsl_random_hypergeom

    interface gsl_random_binomial
        !> Performs a binomial draw
        function gsl_random_binomial(rng, p, n_samples) result(n_drawn) bind(C, name="gsl_ran_binomial")
            import
            type(c_ptr), intent(in), value :: rng
                !! RNG object, originally created by [[f42_random_gsl(module):gsl_alloc_rng(interface)]]
            real(c_double), intent(in), value :: p
                !! Probability to draw a sample from success population
            integer(c_int), intent(in), value :: n_samples
                !! Number of draws to perform
            integer(c_int) :: n_drawn
                !! Number of drawn elements from success population
        end function gsl_random_binomial
    end interface gsl_random_binomial

    !> Opaque handle onto a GSL random number generator. Created by
    !| [[f42_random_gsl(module):create_rng(function)]]; every other procedure in this module
    !| takes one and reports `ERR_POINTER_NULL` (at argument position 1) if it was never
    !| initialized that way, instead of dereferencing a null `c_ptr`.
    type :: rng_t
        type(c_ptr), private :: rng = c_null_ptr
    end type rng_t
contains

    !> Checks that `rng` was produced by `create_rng` and not just default-initialized.
    !|
    !| Private. Every public procedure below calls this first, at argument position 1 (`rng`
    !| is always the first dummy here), and skips its computation -- returning a zeroed result
    !| -- when it fails. This replaces what the upstream GSL wrapper this module was ported
    !| from did with `error stop`: that aborted the whole process, which made the failure
    !| untestable and took every other running test down with it.
    pure subroutine check_rng_initialized(rng, ierr)
        type(rng_t), intent(in) :: rng
            !! RNG object to check
        integer(int32), intent(inout) :: ierr
            !! Error accumulator; unchanged if `rng` is initialized

        if (.not. c_associated(rng%rng)) call set_err_once(ierr, ERR_POINTER_NULL, arg_pos=1_int32)
    end subroutine check_rng_initialized

    !> Creates a random number generator to be used with the RNG procedures of this module
    type(rng_t) function create_rng(ierr, seed) result(rng)
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_ALLOC_FAIL` if GSL could not allocate the generator
        integer(int32), intent(in), optional :: seed
            !! Seed to initialize the `rng` with; defaults to 42 (see `reset_rng`)

        call set_ok(ierr)

        rng%rng = gsl_alloc_rng(gsl_get_rng_type())
        if (.not. c_associated(rng%rng)) then
            call set_err_once(ierr, ERR_ALLOC_FAIL)
            return
        end if

        call reset_rng(rng, ierr, seed)
    end function create_rng

    !> Deallocates all memory related to `rng`. A never-initialized or already-destroyed `rng` is left untouched.
    subroutine destroy_rng(rng)
        type(rng_t), intent(inout) :: rng
            !! RNG object, originally created by [[f42_random_gsl(module):create_rng(function)]]

        if (c_associated(rng%rng)) then
            call gsl_free_rng(rng%rng)
            rng%rng = c_null_ptr
        end if
    end subroutine destroy_rng

    !> Resets `rng` to a given seed (or 42 by default), restarting its stream from the beginning.
    subroutine reset_rng(rng, ierr, seed)
        type(rng_t), intent(in) :: rng
            !! RNG object, originally created by [[f42_random_gsl(module):create_rng(function)]]
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_POINTER_NULL` if `rng` was never initialized
        integer(int32), intent(in), optional :: seed
            !! Seed to reset `rng` to; defaults to 42

        call set_ok(ierr)
        call check_rng_initialized(rng, ierr)
        if (is_err(ierr)) return

        if (present(seed)) then
            call gsl_set_rng_seed(rng%rng, int(seed, kind=c_long))
        else
            call gsl_set_rng_seed(rng%rng, 42_c_long)
        end if
    end subroutine reset_rng

    !> Produces a uniform random number in `[0, 1)`
    real(real64) function random_uniform(rng, ierr) result(rand)
        type(rng_t), intent(in) :: rng
            !! RNG object, originally created by [[f42_random_gsl(module):create_rng(function)]]
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_POINTER_NULL` if `rng` was never initialized

        call set_ok(ierr)
        call check_rng_initialized(rng, ierr)
        if (is_err(ierr)) then
            rand = 0.0_real64
            return
        end if

        rand = gsl_random_uniform(rng%rng)
    end function random_uniform

    !> Performs a hypergeometric draw
    function random_hypergeom(rng, n_population1, n_population2, n_samples, ierr) result(n_drawn)
        type(rng_t), intent(in) :: rng
            !! RNG object, originally created by [[f42_random_gsl(module):create_rng(function)]]
        integer(int32), intent(in) :: n_population1
            !! Size of success population; must be non-negative
        integer(int32), intent(in) :: n_population2
            !! Size of failure population; must be non-negative
        integer(int32), intent(in) :: n_samples
            !! Number of draws to perform; must be non-negative
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_POINTER_NULL` if `rng` was never initialized, `ERR_INVALID_INPUT` for a negative argument
        integer(int32) :: n_drawn
            !! Number of drawn elements from success population

        call set_ok(ierr)
        call check_rng_initialized(rng, ierr)
        call validate_in_range_int(n_population1, ierr, arg_pos=2_int32, min=0_int32)
        call validate_in_range_int(n_population2, ierr, arg_pos=3_int32, min=0_int32)
        call validate_in_range_int(n_samples, ierr, arg_pos=4_int32, min=0_int32)
        if (is_err(ierr)) then
            n_drawn = 0_int32
            return
        end if

        n_drawn = gsl_random_hypergeom(rng%rng, n_population1, n_population2, n_samples)
    end function random_hypergeom

    !> Performs a binomial draw
    function random_binomial(rng, p, n_samples, ierr) result(n_drawn)
        type(rng_t), intent(in) :: rng
            !! RNG object, originally created by [[f42_random_gsl(module):create_rng(function)]]
        real(real64), intent(in) :: p
            !! Probability to draw a sample from success population, in `[0, 1]`
        integer(int32), intent(in) :: n_samples
            !! Number of draws to perform; must be non-negative
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_POINTER_NULL` if `rng` was never initialized, `ERR_INVALID_INPUT` if `p` or `n_samples` is out of range
        integer(int32) :: n_drawn
            !! Number of drawn elements from success population

        call set_ok(ierr)
        call check_rng_initialized(rng, ierr)
        call validate_in_range_real(p, ierr, arg_pos=2_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_int(n_samples, ierr, arg_pos=3_int32, min=0_int32)
        if (is_err(ierr)) then
            n_drawn = 0_int32
            return
        end if

        n_drawn = gsl_random_binomial(rng%rng, p, n_samples)
    end function random_binomial

    !> Returns a uniform random real number `min <= rand_num < max`. If `min > max`, it will be `max <= rand_num < min`. If `min == max`, it will be `min`.
    real(real64) function rand_range(rng, min, max, ierr) result(res)
        type(rng_t), intent(in) :: rng
            !! The rng, originally created by [[f42_random_gsl(module):create_rng(function)]]
        real(real64), intent(in) :: min
            !! Lower bound
        real(real64), intent(in) :: max
            !! Upper bound
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_POINTER_NULL` if `rng` was never initialized

        call set_ok(ierr)
        call check_rng_initialized(rng, ierr)
        call validate_in_range_real(min, ierr, arg_pos=2_int32)
        call validate_in_range_real(max, ierr, arg_pos=3_int32)
        if (is_err(ierr)) then
            res = 0.0_real64
            return
        end if

        res = min + gsl_random_uniform(rng%rng) * (max - min)
    end function rand_range

    !> Performs a binomial multivariate draw
    subroutine random_multinomial(rng, n_populations, population_sizes, total_population, n_to_draw, drawn, ierr)
        type(rng_t), intent(in) :: rng
            !! The rng, originally created by [[f42_random_gsl(module):create_rng(function)]]
        integer(int32), intent(in) :: n_populations
            !! Number of variates/subpopulations; must be at least 1
        integer(int32), dimension(n_populations), intent(in) :: population_sizes
            !! Sizes of subpopulations; each must be non-negative and sum to `total_population`
        integer(int32), intent(in) :: total_population
            !! The total population size -> `sum(population_sizes)`
        integer(int32), intent(in) :: n_to_draw
            !! Number of elements to draw in total; must be non-negative
        integer(int32), dimension(n_populations), intent(out) :: drawn
            !! Drawn sample for `population_sizes` with `sum(drawn) == n_to_draw`
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_POINTER_NULL` if `rng` was never initialized, `ERR_INVALID_INPUT`/`ERR_EMPTY_INPUT` for a bad size or count

        integer(int32) :: remaining_population, i_population, remaining_draws
        real(real64) :: p

        call set_ok(ierr)
        call check_rng_initialized(rng, ierr)
        call validate_dimension_size(n_populations, ierr, arg_pos=2_int32)
        call validate_in_range_int(total_population, ierr, arg_pos=4_int32, min=0_int32)
        call validate_in_range_int(n_to_draw, ierr, arg_pos=5_int32, min=0_int32)
        if (.not. is_err(ierr)) then
            if (n_populations >= 1_int32) then
                if (any(population_sizes < 0_int32)) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=3_int32)
                ! Cross-check, in the spirit of a mask's `n_selected_<arg>` (codegen_guide.md
                ! §5.3): `total_population` is documented as `sum(population_sizes)`, so a
                ! mismatch blames the derived value's own position, not the array it was
                ! supposed to summarize.
                if (sum(population_sizes) /= total_population) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
            end if
        end if
        if (is_err(ierr)) then
            drawn = 0_int32
            return
        end if

        remaining_population = total_population
        remaining_draws = n_to_draw

        ! Draw for each population, last one gets rest
        do i_population = 1, n_populations - 1
            associate (&
                current_population => population_sizes(i_population),&
                n_drawn => drawn(i_population)&
            )

                if (remaining_population > 0) then
                    ! Draw elements of current population from pool with remaining elements
                    p = real(current_population, real64) / real(remaining_population, real64)
                    n_drawn = min(remaining_draws, gsl_random_binomial(rng%rng, p, remaining_draws))

                    ! For further draws, the current subpopulation is not taken into account anymore
                    remaining_population = remaining_population - current_population
                    remaining_draws = remaining_draws - n_drawn
                else
                    n_drawn = 0_int32
                end if

            end associate
        end do

        drawn(n_populations) = remaining_draws
    end subroutine random_multinomial

    !> Performs a hypergeometric multivariate draw
    subroutine random_multiv_hypergeom(rng, n_populations, population_sizes, total_population, n_to_draw, drawn, ierr)
        type(rng_t), intent(in) :: rng
            !! The rng, originally created by [[f42_random_gsl(module):create_rng(function)]]
        integer(int32), intent(in) :: n_populations
            !! Number of variates/subpopulations; must be at least 1
        integer(int32), dimension(n_populations), intent(inout) :: population_sizes
            !! Sizes of subpopulations (will be reduced by the number of drawn elements per population -> will be the remaining pool);
            !! each must be non-negative and sum to `total_population`.
            !! Known limitation, ported as-is from the upstream algorithm: only elements
            !! `1..n_populations-1` are actually reduced -- the loop below never revisits
            !! `population_sizes(n_populations)`, so that one element is returned unchanged
            !! regardless of how many of it were drawn. Callers that need the true remaining
            !! pool for the last subpopulation must compute it themselves as
            !! `population_sizes(n_populations) - drawn(n_populations)`.
        integer(int32), intent(in) :: total_population
            !! The total population size -> `sum(population_sizes)`
        integer(int32), intent(in) :: n_to_draw
            !! Number of elements to draw in total; must be non-negative
        integer(int32), dimension(n_populations), intent(out) :: drawn
            !! Drawn sample for `population_sizes` with `sum(drawn) == n_to_draw`
        integer(int32), intent(out) :: ierr
            !! Error code; `ERR_POINTER_NULL` if `rng` was never initialized, `ERR_INVALID_INPUT`/`ERR_EMPTY_INPUT` for a bad size or count

        integer(int32) :: remaining_population, i_population, remaining_draws

        call set_ok(ierr)
        call check_rng_initialized(rng, ierr)
        call validate_dimension_size(n_populations, ierr, arg_pos=2_int32)
        call validate_in_range_int(total_population, ierr, arg_pos=4_int32, min=0_int32)
        call validate_in_range_int(n_to_draw, ierr, arg_pos=5_int32, min=0_int32)
        if (.not. is_err(ierr)) then
            if (n_populations >= 1_int32) then
                if (any(population_sizes < 0_int32)) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=3_int32)
                ! Cross-check, in the spirit of a mask's `n_selected_<arg>` (codegen_guide.md
                ! §5.3): `total_population` is documented as `sum(population_sizes)`, so a
                ! mismatch blames the derived value's own position, not the array it was
                ! supposed to summarize.
                if (sum(population_sizes) /= total_population) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
            end if
        end if
        if (is_err(ierr)) then
            drawn = 0_int32
            return
        end if

        remaining_population = total_population
        remaining_draws = min(total_population, n_to_draw)
        drawn = 0_int32

        ! Draw for each population, last one gets rest
        do i_population = 1, n_populations - 1
            associate (&
                current_population => population_sizes(i_population),&
                n_drawn => drawn(i_population)&
            )

                if (remaining_population > 0) then
                        ! Draw elements of current population from pool with remaining elements
                        remaining_population = remaining_population - current_population
                        n_drawn = gsl_random_hypergeom(rng%rng, current_population, remaining_population, remaining_draws)

                        ! For further draws, the current subpopulation is not taken into account anymore
                        current_population = current_population - n_drawn
                        remaining_draws = remaining_draws - n_drawn
                else
                    n_drawn = 0_int32
                end if

            end associate
        end do

        drawn(n_populations) = remaining_draws
    end subroutine random_multiv_hypergeom
end module f42_random_gsl
