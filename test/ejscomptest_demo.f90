program ejscomptest_demo
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_ejscomptest, only: ejscomptest_run, EJS_OK
    implicit none

    integer(int32), parameter :: n_studies = 3_int32
    integer(int32), parameter :: n_values = 18_int32
    integer(int32), parameter :: n_bins = 20_int32
    integer(int32), parameter :: n_boot = 100_int32
    integer(int32), parameter :: n_perm = 250_int32

    real(real64) :: values(n_values)
    integer(int32) :: study_offsets(n_studies + 1)
    real(real64) :: jsd(n_studies), ci_low(n_studies), ci_high(n_studies), pval(n_studies)
    integer(int32) :: n_i(n_studies)
    integer(int32) :: ierr, i

    ! Flattened per-study values using offsets below.
    values = [ &
        0.0_real64, 0.0_real64, 1.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, & ! study 1
        0.0_real64, 1.2_real64, 2.8_real64, 2.9_real64, 5.1_real64, 5.2_real64, & ! study 2
        0.0_real64, 0.0_real64, 0.1_real64, 0.2_real64, 9.0_real64, 10.0_real64  & ! study 3
    ]

    study_offsets = [1_int32, 7_int32, 13_int32, 19_int32]

    call ejscomptest_run(values=values, study_offsets=study_offsets, n_studies=n_studies, &
                         use_log_mode=.true., n_bins=n_bins, n_boot=n_boot, n_perm=n_perm, &
                         jsd_obs=jsd, ci_low=ci_low, ci_high=ci_high, p_value=pval, n_i=n_i, ierr=ierr, &
                         min_value_keep=0.0_real64, random_seed_value=42_int32)

    if (ierr /= EJS_OK) then
        print *, "eJSCompTest failed with ierr=", ierr
        stop 1
    end if

    print *, "Study   N_i   JSD_obs      CI_low       CI_high      p_value"
    do i = 1, n_studies
        write(*,'(I3,2X,I5,2X,F10.6,2X,F10.6,2X,F10.6,2X,F10.6)') i, n_i(i), jsd(i), ci_low(i), ci_high(i), pval(i)
    end do
end program ejscomptest_demo
