!> This module interfaces with the random number generation from the GNU Scientific Library.
!| Thus, for selecting the RNG algorithm, use the environment variable `GSL_RNG_TYPE` as defined [in the documentation](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-environment-variables):
module f42_random_gsl
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64, int64
    use, intrinsic :: iso_c_binding, only: c_ptr, c_associated, c_double, c_int, c_long
    implicit none

    interface get_rng_type
        !> Returns algorithm for the RNG, set by the environment variable `GSL_RNG_TYPE`
        function get_rng_type() result(rng_type) bind(C, name="gsl_rng_env_setup")
            import
            type(c_ptr) :: rng_type
        end function get_rng_type
    end interface get_rng_type

    interface alloc_rng
        !> Creates a RNG structure to generate random numbers with, uses seed set by environment variable `GSL_RNG_SEED`
        function alloc_rng(rng_type) result(rng) bind(C, name="gsl_rng_alloc")
            import
            type(c_ptr), intent(in), value :: rng_type
                !! RNG type returned by [[f42_random_gsl(module):get_rng_type(interface)]]
            type(c_ptr) :: rng
                !! Created RNG object
        end function alloc_rng
    end interface alloc_rng

    interface free_rng
        !> Deallocates all memory related to `rng`
        subroutine free_rng(rng) bind(C, name="gsl_rng_free")
            import
            type(c_ptr), intent(in), value :: rng
                !! RNG object, originally created by [[f42_random_gsl(module):alloc_rng(interface)]]
        end subroutine free_rng
    end interface free_rng

    interface set_rng_seed
        !> Sets the seed for a RNG structure to a specific one. Will reset the RNG stream if it is the same seed as already used.
        subroutine set_rng_seed(rng, seed) bind(C, name="gsl_rng_set")
            import
            type(c_ptr), intent(in), value :: rng
                !! RNG object, originally created by [[f42_random_gsl(module):alloc_rng(interface)]]
            integer(c_long), intent(in) :: seed
                !! Seed to initialize `rng` with
        end subroutine set_rng_seed
    end interface set_rng_seed

    interface random_uniform
        function random_uniform(rng) result(rand) bind(C, name="gsl_rng_uniform")
            import
            type(c_ptr), intent(in), value :: rng
            real(c_double) :: rand
        end function random_uniform
    end interface random_uniform

    interface random_hypergeom
        !> Performs a hypergeometric draw
        function random_hypergeom(rng, n_population1, n_population2, n_samples) result(n_drawn) bind(C, name="gsl_ran_hypergeometric")
            import
            type(c_ptr), intent(in), value :: rng
                !! RNG object, originally created by [[f42_random_gsl(module):alloc_rng(interface)]]
            integer(c_int), intent(in), value :: n_population1
                !! Size of success population
            integer(c_int), intent(in), value :: n_population2
                !! Size of failure population
            integer(c_int), intent(in), value :: n_samples
                !! Number of draws to perform
            integer(c_int) :: n_drawn
                !! Number of drawn elements from success population
        end function random_hypergeom
    end interface random_hypergeom

    interface random_binomial
        !> Performs a binomial draw
        function random_binomial(rng, p, n_samples) result(n_drawn) bind(C, name="gsl_ran_binomial")
            import
            type(c_ptr), intent(in), value :: rng
                !! RNG object, originally created by [[f42_random_gsl(module):alloc_rng(interface)]]
            real(c_double), intent(in), value :: p
                !! Propability to draw a sample from success population
            integer(c_int), intent(in), value :: n_samples
                !! Number of draws to perform
            integer(c_int) :: n_drawn
                !! Number of drawn elements from success population
        end function random_binomial
    end interface random_binomial
contains

    !> Creates a random number generator to be used with RNG functions/subroutines
    type(c_ptr) function create_rng(seed) result(rng)
        integer(int32), intent(in), optional :: seed
            !! Seed to initialize the `rng` with

        rng = alloc_rng(get_rng_type())
        call reset_rng(rng, seed)
    end function create_rng

    !> Creates a random number generator to be used with RNG functions/subroutines
    subroutine reset_rng(rng, seed)
        type(c_ptr), intent(in) :: rng
            !! RNG object, originally created by [[f42_random_gsl(module):create_rng(function)]]
        integer(int32), intent(in), optional :: seed
            !! Seed to reset `rng` to

        if (present(seed)) then
            call set_rng_seed(rng, int(seed, kind=c_long))
        else
            call set_rng_seed(rng, 42_c_long)
        end if
    end subroutine reset_rng


    !> Returns a uniform random real number `min <= rand_num < max`. If `min > max`, it will be `max <= rand_num < min`. If `min == max`, it will be `min`.
    real(real64) function rand_range(rng, min, max) result(res)
        type(c_ptr), intent(in) :: rng
            !! The rng, originally created by [[f42_random_gsl(module):create_rng(function)]]
        real(real64), intent(in) :: min
            !! Lower bound
        real(real64), intent(in) :: max
            !! Upper bound

        res = min + random_uniform(rng) * (max - min)
    end function rand_range

    !> Performs a binomial multivariate draw
    subroutine random_multinomial(rng, population_sizes, n_populations, total_population, n_to_draw, drawn)
        integer(int32), intent(in) :: n_populations
            !! Number of variates/subpopulations
        type(c_ptr), intent(in) :: rng
            !! The rng, originally created by [[f42_random_gsl(module):create_rng(function)]]
        integer(int32), dimension(n_populations), intent(in) :: population_sizes
            !! Sizes of subpopulations
        integer(int32), intent(in) :: total_population
            !! The total population size -> `sum(population_sizes)`
        integer(int32), intent(in) :: n_to_draw
            !! Number of elements to draw in total
        integer(int32), dimension(n_populations), intent(out) :: drawn
            !! Drawn sample for `population_sizes` with `sum(drawn) == n_to_draw`

        integer(int32) :: remaining_population, i_population, remaining_draws
        real(real64) :: p

        remaining_population = total_population
        remaining_draws = n_to_draw
        drawn = 0_int32

        ! Draw for each population, last one gets rest
        do i_population = 1, n_populations - 1
            associate (&
                current_population => population_sizes(i_population),&
                n_drawn => drawn(i_population)&
            )

                if (remaining_population > 0) then
                    ! Draw elements of current population from pool with remaining elements
                    p = real(current_population, real64) / real(remaining_population, real64)
                    n_drawn = min(remaining_draws, random_binomial(rng, p, remaining_draws))

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
    subroutine random_multiv_hypergeom(rng, population_sizes, n_populations, total_population, n_to_draw, drawn)
        integer(int32), intent(in) :: n_populations
            !! Number of variates/subpopulations
        type(c_ptr), intent(in) :: rng
            !! The rng, originally created by [[f42_random_gsl(module):create_rng(function)]]
        integer(int32), dimension(n_populations), intent(inout) :: population_sizes
            !! Sizes of subpopulations (will be reduced by the number of drawn elements per population -> will be the remaining pool)
        integer(int32), intent(in) :: total_population
            !! The total population size -> `sum(population_sizes)`
        integer(int32), intent(in) :: n_to_draw
            !! Number of elements to draw in total
        integer(int32), dimension(n_populations), intent(out) :: drawn
            !! Drawn sample for `population_sizes` with `sum(drawn) == n_to_draw`

        integer(int32) :: remaining_population, i_population, remaining_draws

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
                        n_drawn = random_hypergeom(rng, current_population, remaining_population, remaining_draws)

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