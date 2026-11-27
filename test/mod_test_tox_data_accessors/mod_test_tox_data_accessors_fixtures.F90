!> What the children of the `tox_data_accessors` suite share: a data set of three genes in two
!| families with two samples, small enough that every value it holds can be written out, and the
!| value an output holds before a call.
!|
!| Each value says where it is: expression `10*gene + sample`, centroid `100*family + sample`,
!| and the shift vector of gene j holds its family's centroid in rows 1-2 (`100*family + row`)
!| and `-(10*j + row)` in rows 3-4. So a value copied from the wrong row or column shows.
module mod_test_tox_data_accessors_fixtures
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    private
    public :: N_SAMPLES, N_GENES, N_FAMILIES, UNWRITTEN
    public :: FIXTURE_GENE_IDS, FIXTURE_FAMILY_IDS, FIXTURE_GENE_TO_FAM, expression_vectors, family_centroids, shift_vectors

    integer(int32), parameter :: N_SAMPLES = 2, N_GENES = 3, N_FAMILIES = 2

    !> What an output holds before a call, so that an element the procedure never writes shows.
    !| No value of the data set equals it.
    real(real64), parameter :: UNWRITTEN = -777.0_real64

    character(len=8), parameter :: FIXTURE_GENE_IDS(N_GENES) = ["g1", "g2", "g3"]
    character(len=8), parameter :: FIXTURE_FAMILY_IDS(N_FAMILIES) = ["OG1", "OG2"]
    !> g1 in family 2, g2 in none (M_GENE_TO_FAM_SENTINEL), g3 in family 1.
    integer(int32), parameter :: FIXTURE_GENE_TO_FAM(N_GENES) = [2, 0, 1]

contains

    !> (samples, genes): gene j's value of sample i is 10*j + i.
    function expression_vectors() result(values)
        real(real64) :: values(N_SAMPLES, N_GENES)

        values = reshape([11.0_real64, 12.0_real64, 21.0_real64, 22.0_real64, 31.0_real64, 32.0_real64], &
                         [N_SAMPLES, N_GENES])
    end function expression_vectors

    !> (samples, families): family k's value of sample i is 100*k + i.
    function family_centroids() result(values)
        real(real64) :: values(N_SAMPLES, N_FAMILIES)

        values = reshape([101.0_real64, 102.0_real64, 201.0_real64, 202.0_real64], [N_SAMPLES, N_FAMILIES])
    end function family_centroids

    !> (2*samples, genes): gene j's start in rows 1-2 (its family's centroid, 0 for g2), its
    !| shift -(10*j + row) in rows 3-4.
    function shift_vectors() result(values)
        real(real64) :: values(2*N_SAMPLES, N_GENES)

        values = reshape([201.0_real64, 202.0_real64, -13.0_real64, -14.0_real64, &
                          0.0_real64, 0.0_real64, -23.0_real64, -24.0_real64, &
                          101.0_real64, 102.0_real64, -33.0_real64, -34.0_real64], [2*N_SAMPLES, N_GENES])
    end function shift_vectors
end module mod_test_tox_data_accessors_fixtures
