!> @file mod_test_empirical_pvalue.f90
!> Unit test suite for EDF (Empirical Distribution Function)
!> Contains dedicated tests for compute_edf from f42_utils.

module mod_test_empirical_pvalue
  use asserts
  use f42_utils, only: lower_bound_ge, compute_empirical_p_values
  use tox_errors
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use test_suite
  implicit none
  public

contains

  !> Get array of all available EPV tests.
  function get_all_tests_empirical_pvalue() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(9))

    all_tests(1) = test_case("test_lower_bound_ge_singleton", test_lower_bound_ge_singleton)
    all_tests(2) = test_case("test_lower_bound_ge_increasing", test_lower_bound_ge_increasing)
    all_tests(3) = test_case("test_lower_bound_ge_duplicates", test_lower_bound_ge_duplicates)
    all_tests(4) = test_case("test_empirical_pvalues_denom_le_zero", test_empirical_pvalues_denom_le_zero)
    all_tests(5) = test_case("test_empirical_pvalues_negative_is_one", test_empirical_pvalues_negative_is_one)
    all_tests(6) = test_case("test_empirical_pvalues_extremes", test_empirical_pvalues_extremes)
    all_tests(7) = test_case("test_empirical_pvalues_duplicates", test_empirical_pvalues_duplicates)
    all_tests(8) = test_case("test_empirical_pvalues_matches_naive_all", test_empirical_pvalues_matches_naive_all)
    all_tests(9) = test_case("test_empirical_pvalues_monotonicity", test_empirical_pvalues_monotonicity)
  end function get_all_tests_empirical_pvalue

        ! ==========================================================================
    ! TESTS: lower_bound_ge + compute_empirical_p_values
    ! ==========================================================================

    !> Helper: naive count of elements >= d in the perm-sorted distribution
    pure integer(int32) function naive_count_ge(sorted_rdi, perm, n, d) result(k)
      use, intrinsic :: iso_fortran_env, only: int32, real64
      implicit none
      real(real64), intent(in) :: sorted_rdi(n)
      integer(int32), intent(in) :: perm(n)
      integer(int32), intent(in) :: n
      real(real64), intent(in) :: d
      integer(int32) :: i
      k = 0_int32
      do i = 1, n
        if (sorted_rdi(perm(i)) >= d) k = k + 1_int32
      end do
    end function naive_count_ge


    ! ==========================================================================
    ! LOWER_BOUND_ge
    ! (First position pos in [1..n] such that sorted_rdi(perm(pos)) >= x; returns n+1 if none)
    ! ==========================================================================

    !> Edge case: n=1, test x below/eq/above
    subroutine test_lower_bound_ge_singleton()
      use, intrinsic :: iso_fortran_env, only: int32, real64
      implicit none
      integer(int32), parameter :: n = 1
      real(real64) :: s(n)
      integer(int32) :: perm(n)
      integer(int32) :: pos

      s(1) = 2.0_real64
      perm(1) = 1_int32

      pos = lower_bound_ge(s, perm, n, 1.0_real64)
      call assert_equal_int(pos, 1_int32, "lower_bound perm singleton: x<min -> pos=1")

      pos = lower_bound_ge(s, perm, n, 2.0_real64)
      call assert_equal_int(pos, 1_int32, "lower_bound perm singleton: x==min -> pos=1")

      pos = lower_bound_ge(s, perm, n, 3.0_real64)
      call assert_equal_int(pos, 2_int32, "lower_bound perm singleton: x>max -> pos=n+1")
    end subroutine test_lower_bound_ge_singleton


    !> Strictly increasing: verify typical positions
    subroutine test_lower_bound_ge_increasing()
      use, intrinsic :: iso_fortran_env, only: int32, real64
      implicit none
      integer(int32), parameter :: n = 5
      real(real64) :: s(n)
      integer(int32) :: perm(n)
      integer(int32) :: pos

      ! We'll store values in-place but keep perm identity (already sorted)
      s = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
      perm = [1_int32, 2_int32, 3_int32, 4_int32, 5_int32]

      pos = lower_bound_ge(s, perm, n, 0.5_real64)
      call assert_equal_int(pos, 1_int32, "lower_bound perm: x<min -> 1")

      pos = lower_bound_ge(s, perm, n, 1.0_real64)
      call assert_equal_int(pos, 1_int32, "lower_bound perm: x==min -> 1")

      pos = lower_bound_ge(s, perm, n, 2.5_real64)
      call assert_equal_int(pos, 3_int32, "lower_bound perm: x=2.5 -> first >=3 at pos=3")

      pos = lower_bound_ge(s, perm, n, 5.0_real64)
      call assert_equal_int(pos, 5_int32, "lower_bound perm: x==max -> pos=n")

      pos = lower_bound_ge(s, perm, n, 6.0_real64)
      call assert_equal_int(pos, 6_int32, "lower_bound perm: x>max -> n+1")
    end subroutine test_lower_bound_ge_increasing


    !> Duplicates: must return FIRST occurrence (lower bound)
    subroutine test_lower_bound_ge_duplicates()
      use, intrinsic :: iso_fortran_env, only: int32, real64
      implicit none
      integer(int32), parameter :: n = 6
      real(real64) :: s(n)
      integer(int32) :: perm(n)
      integer(int32) :: pos

      ! Values are [0, 2, 2, 2, 5, 9] but stored scrambled in s
      s = [2.0_real64, 0.0_real64, 9.0_real64, 2.0_real64, 5.0_real64, 2.0_real64]
      ! Indices sorted by s: 2(0), 1(2), 4(2), 6(2), 5(5), 3(9)
      perm = [2_int32, 1_int32, 4_int32, 6_int32, 5_int32, 3_int32]

      pos = lower_bound_ge(s, perm, n, 2.0_real64)
      call assert_equal_int(pos, 2_int32, "lower_bound perm duplicates: x=2 -> first 2 at pos=2")

      pos = lower_bound_ge(s, perm, n, 3.0_real64)
      call assert_equal_int(pos, 5_int32, "lower_bound perm duplicates: x=3 -> first >=5 at pos=5")

      pos = lower_bound_ge(s, perm, n, 9.0_real64)
      call assert_equal_int(pos, 6_int32, "lower_bound perm duplicates: x=9 -> pos=6")
    end subroutine test_lower_bound_ge_duplicates


    ! ==========================================================================
    ! compute_empirical_p_values
    ! ==========================================================================

    !> Edge: denom <= 0 => p_values = 1
    subroutine test_empirical_pvalues_denom_le_zero()
      use, intrinsic :: iso_fortran_env, only: int32, real64
      implicit none
      integer(int32), parameter :: n = 4
      real(real64) :: rdi(n), s(n), p(n)
      integer(int32) :: perm(n)

      ! distribution values stored but perm sorts them already
      s = [2.0_real64, 0.0_real64, 3.0_real64, 1.0_real64]          ! sorted would be [0,1,2,3]
      perm = [2_int32, 4_int32, 1_int32, 3_int32]

      rdi = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64]

      ! choose c so denom = n + c <= 0
      call compute_empirical_p_values(n, rdi, s, perm, p, -4.0_real64)

      call assert_equal_real(p(1), 1.0_real64, 1d-12, "denom<=0 -> p=1 (1)")
      call assert_equal_real(p(2), 1.0_real64, 1d-12, "denom<=0 -> p=1 (2)")
      call assert_equal_real(p(3), 1.0_real64, 1d-12, "denom<=0 -> p=1 (3)")
      call assert_equal_real(p(4), 1.0_real64, 1d-12, "denom<=0 -> p=1 (4)")
    end subroutine test_empirical_pvalues_denom_le_zero


    !> Negative rdi => p=1, regardless of distribution
    subroutine test_empirical_pvalues_negative_is_one()
      use, intrinsic :: iso_fortran_env, only: int32, real64
      implicit none
      integer(int32), parameter :: n = 5
      real(real64) :: rdi(n), s(n), p(n)
      integer(int32) :: perm(n)
      real(real64), parameter :: c = 1.0_real64

      s    = [4.0_real64, 0.0_real64, 2.0_real64, 1.0_real64, 3.0_real64]  ! sorted [0,1,2,3,4]
      perm = [2_int32, 4_int32, 3_int32, 5_int32, 1_int32]

      rdi = [-1.0_real64, 0.0_real64, 2.0_real64, -5.0_real64, 4.0_real64]

      call compute_empirical_p_values(n, rdi, s, perm, p, c)

      call assert_equal_real(p(1), 1.0_real64, 1d-12, "negative rdi -> p=1")
      call assert_equal_real(p(4), 1.0_real64, 1d-12, "negative rdi -> p=1")
    end subroutine test_empirical_pvalues_negative_is_one


    !> Extremes: d below min => p=1; d above max => p = c/(n+c)
    subroutine test_empirical_pvalues_extremes()
      use, intrinsic :: iso_fortran_env, only: int32, real64
      implicit none
      integer(int32), parameter :: n = 4
      real(real64) :: rdi(n), s(n), p(n)
      integer(int32) :: perm(n)
      real(real64), parameter :: c = 1.0_real64
      real(real64) :: denom, expected

      ! Distribution is [0,1,2,3] but stored scrambled
      s = [2.0_real64, 0.0_real64, 3.0_real64, 1.0_real64]
      perm = [2_int32, 4_int32, 1_int32, 3_int32]

      denom = real(n, real64) + c

      ! gene values: negative, equal min, above max, equal max
      rdi = [-0.5_real64, 0.0_real64, 10.0_real64, 3.0_real64]

      call compute_empirical_p_values(n, rdi, s, perm, p, c)

      call assert_equal_real(p(1), 1.0_real64, 1d-12, "negative -> p=1")

      ! d=0: all >=0 => count=n => (n+c)/(n+c)=1
      call assert_equal_real(p(2), 1.0_real64, 1d-12, "d<=min -> p=1")

      ! d=10: count=0 => c/(n+c)
      expected = c / denom
      call assert_equal_real(p(3), expected, 1d-12, "d>max -> p=c/(n+c)")

      ! d=3: count=1 => (1+c)/(n+c)
      expected = (1.0_real64 + c) / denom
      call assert_equal_real(p(4), expected, 1d-12, "d==max -> p=(1+c)/(n+c)")
    end subroutine test_empirical_pvalues_extremes


    !> Duplicates: all genes with same d get same p-value; counts include equals (>=)
    subroutine test_empirical_pvalues_duplicates()
      use, intrinsic :: iso_fortran_env, only: int32, real64
      implicit none
      integer(int32), parameter :: n = 6
      real(real64) :: rdi(n), s(n), p(n)
      integer(int32) :: perm(n)
      real(real64), parameter :: c = 1.0_real64
      real(real64) :: denom, expected_for_2

      ! distribution values [0,2,2,2,5,9] stored scrambled
      s = [2.0_real64, 0.0_real64, 9.0_real64, 2.0_real64, 5.0_real64, 2.0_real64]
      perm = [2_int32, 1_int32, 4_int32, 6_int32, 5_int32, 3_int32]

      denom = real(n, real64) + c

      ! genes include several 2's
      rdi = [2.0_real64, 9.0_real64, 1.0_real64, 2.0_real64, 5.0_real64, 2.0_real64]

      call compute_empirical_p_values(n, rdi, s, perm, p, c)

      ! for d=2: values >=2 are [2,2,2,5,9] => count=5 => (5+c)/(n+c)
      expected_for_2 = (5.0_real64 + c) / denom

      call assert_equal_real(p(1), expected_for_2, 1d-12, "d=2 -> (5+c)/(n+c)")
      call assert_equal_real(p(4), expected_for_2, 1d-12, "d=2 -> same p-value")
      call assert_equal_real(p(6), expected_for_2, 1d-12, "d=2 -> same p-value")
    end subroutine test_empirical_pvalues_duplicates


    !> Exhaustive correctness vs naive counting for all genes (fixed values)
    subroutine test_empirical_pvalues_matches_naive_all()
      use, intrinsic :: iso_fortran_env, only: int32, real64
      implicit none
      integer(int32), parameter :: n = 7
      real(real64) :: rdi(n), s(n), p(n)
      integer(int32) :: perm(n)
      real(real64), parameter :: c = 1.0_real64
      real(real64) :: denom, expected
      integer(int32) :: i, k

      ! distribution sorted would be: [0.0, 0.5, 1.0, 1.0, 2.0, 4.0, 10.0]
      ! store scrambled:
      s = [1.0_real64, 10.0_real64, 0.0_real64, 1.0_real64, 4.0_real64, 0.5_real64, 2.0_real64]
      ! indices sorted by s: 3(0), 6(0.5), 1(1.0), 4(1.0), 7(2.0), 5(4.0), 2(10.0)
      perm = [3_int32, 6_int32, 1_int32, 4_int32, 7_int32, 5_int32, 2_int32]

      denom = real(n, real64) + c

      rdi = [-1.0_real64, 0.0_real64, 1.0_real64, 1.7_real64, 2.0_real64, 10.0_real64, 99.0_real64]

      call compute_empirical_p_values(n, rdi, s, perm, p, c)

      do i = 1, n
        if (rdi(i) < 0.0_real64) then
          expected = 1.0_real64
        else
          k = naive_count_ge(s, perm, n, rdi(i))
          expected = (real(k, real64) + c) / denom
        end if
        call assert_equal_real(p(i), expected, 1d-12, "p-value must match naive counting for each gene")
      end do
    end subroutine test_empirical_pvalues_matches_naive_all


    !> Monotonicity sanity: if d1 < d2 then p(d1) >= p(d2) for non-negative d
    subroutine test_empirical_pvalues_monotonicity()
      use, intrinsic :: iso_fortran_env, only: int32, real64
      implicit none
      integer(int32), parameter :: n = 5
      real(real64) :: s(n), p(n), rdi(n)
      integer(int32) :: perm(n)
      real(real64), parameter :: c = 1.0_real64

      ! distribution [0,1,2,3,4] stored scrambled
      s = [2.0_real64, 0.0_real64, 4.0_real64, 1.0_real64, 3.0_real64]
      perm = [2_int32, 4_int32, 1_int32, 5_int32, 3_int32]

      rdi = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]

      call compute_empirical_p_values(n, rdi, s, perm, p, c)

      call assert_true(p(1) >= p(2), "monotonic: p(0) >= p(1)")
      call assert_true(p(2) >= p(3), "monotonic: p(1) >= p(2)")
      call assert_true(p(3) >= p(4), "monotonic: p(2) >= p(3)")
      call assert_true(p(4) >= p(5), "monotonic: p(3) >= p(4)")
    end subroutine test_empirical_pvalues_monotonicity



end module mod_test_empirical_pvalue
