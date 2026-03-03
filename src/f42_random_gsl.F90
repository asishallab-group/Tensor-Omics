!> This module interfaces with the random number generation from the GNU Scientific Library.
!| Thus, for selecting the RNG algorithm or seed, use the environment variables `GSL_RNG_TYPE` and `GSL_RNG_SEED` as defined [in the documentation](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-environment-variables):
module f42_random_gsl
    use, intrinsic :: iso_fortran_env, only: int32, real64, int64
    use, intrinsic :: iso_c_binding, only: c_ptr, c_associated, c_double, c_int
    implicit none

    interface
        function get_rng_type() result(rng_type) bind(C, name="gsl_rng_env_setup")
            import
            type(c_ptr) :: rng_type
        end function get_rng_type

        function alloc_rng(rng_type) result(rng) bind(C, name="gsl_rng_alloc")
            import
            type(c_ptr), intent(in), value :: rng_type
            type(c_ptr) :: rng
        end function alloc_rng

        subroutine free_rng(rng) bind(C, name="gsl_rng_free")
            import
            type(c_ptr), intent(in), value :: rng
        end subroutine free_rng

        function random_uniform(rng) result(rand) bind(C, name="gsl_rng_uniform")
            import
            type(c_ptr), intent(in), value :: rng
            real(c_double) :: rand
        end function random_uniform

        function random_hypergeom(rng, n_population1, n_population2, n_samples) result(n_drawn) bind(C, name="gsl_ran_hypergeometric")
            import
            type(c_ptr), intent(in), value :: rng
            integer(c_int), intent(in), value :: n_population1
            integer(c_int), intent(in), value :: n_population2
            integer(c_int), intent(in), value :: n_samples
            integer(c_int) :: n_drawn
        end function random_hypergeom

        function random_binomial(rng, p, n_samples) result(n_drawn) bind(C, name="gsl_ran_binomial")
            import
            type(c_ptr), intent(in), value :: rng
            real(real64), intent(in), value :: p
            integer(c_int), intent(in), value :: n_samples
            integer(c_int) :: n_drawn
        end function random_binomial
    end interface

contains
    type(c_ptr) function create_rng() result(rng)
        rng = alloc_rng(get_rng_type())
    end function create_rng

    real(real64) function rand_range(min, max, rng) result(res)
        type(c_ptr), intent(in) :: rng
        real(real64), intent(in) :: min
        real(real64), intent(in) :: max

        res = min + random_uniform(rng) * (max - min)
    end function rand_range
end module f42_random_gsl