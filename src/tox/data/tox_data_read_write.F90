#include <src/macros.h>

!> Wrappers for serialization/deserialization of arrays
!|
!| Each `save_*` writes its array to a file in the layout of
!| [[f42_serde_arrays_utils(module)]], and each `load_*` reads one back into an array of the
!| file's type, rank and shape. They store and read every value as it is, NaN and Inf included.
!| After an error, a loaded array is undefined.
module tox_data_read_write
    use iso_fortran_env, only: int32, real64
    use f42_serde_arrays_serialize, only: serialize_int_1d, serialize_real_2D, serialize_char_1D
    use f42_serde_arrays_deserialize, only: deserialize_int_1D, deserialize_real_2D, deserialize_char_1D
    M_IMPLICIT_NONE
    private

    public :: save_gene_ids, save_expression_vectors, save_gene_to_family, save_family_ids
    public :: save_family_centroids, save_shift_vectors
    public :: load_gene_ids, load_expression_vectors, load_gene_to_family, load_family_ids
    public :: load_family_centroids, load_shift_vectors

contains
    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for serialize_char_1D function
    subroutine save_gene_ids(gene_ids, filename, ierr)
        character(len=*), intent(in), contiguous :: gene_ids(:)
            !! Gene IDs to write
        character(len=*), intent(in) :: filename
            !! Filename to write to
        integer(int32), intent(out) :: ierr
            !! Error code (0 if successful)

        call serialize_char_1D(gene_ids, filename, ierr)
    end subroutine

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for serialize_real_2D function
    subroutine save_expression_vectors(expression_vectors, filename, ierr)
        character(len=*), intent(in) :: filename
            !! Filename to write to
        real(real64), intent(in), contiguous :: expression_vectors(:, :)
            !! Expression vectors to write
        integer(int32), intent(out) :: ierr
            !! Error code (0 if successful)

        call serialize_real_2D(expression_vectors, filename, ierr)
    end subroutine

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for serialize_int_1D function
    subroutine save_gene_to_family(gene_to_fam, filename, ierr)
        character(len=*), intent(in) :: filename
            !! Filename to write to
        integer(int32), intent(in), contiguous :: gene_to_fam(:)
            !! Gene to family mapping to write
        integer(int32), intent(out) :: ierr
            !! Error code (0 if successful)

        call serialize_int_1D(gene_to_fam, filename, ierr)
    end subroutine

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for serialize_char_1D function
    subroutine save_family_ids(family_ids, filename, ierr)
        character(len=*), intent(in) :: filename
            !! Filename to write to
        character(len=*), intent(in), contiguous :: family_ids(:)
            !! Family IDs to write
        integer(int32), intent(out) :: ierr
            !! Error code (0 if successful)

        call serialize_char_1D(family_ids, filename, ierr)
    end subroutine

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for serialize_real_2D function
    subroutine save_family_centroids(family_centroids, filename, ierr)
        character(len=*), intent(in) :: filename
            !! Filename to write to
        real(real64), intent(in), contiguous :: family_centroids(:, :)
            !! Family centroids to write
        integer(int32), intent(out) :: ierr
            !! Error code (0 if successful)

        call serialize_real_2D(family_centroids, filename, ierr)
    end subroutine

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for serialize_real_2D function
    subroutine save_shift_vectors(shift_vectors, filename, ierr)
        character(len=*), intent(in) :: filename
            !! Filename to write to
        real(real64), intent(in), contiguous :: shift_vectors(:, :)
            !! Shift vectors to write
        integer(int32), intent(out) :: ierr
            !! Error code (0 if successful)

        call serialize_real_2D(shift_vectors, filename, ierr)
    end subroutine

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for deserialize_char_1D function
    subroutine load_gene_ids(gene_ids, filename, ierr)
        character(len=*), intent(in)  :: filename
            !! Filename to read from
        character(len=*), intent(out), contiguous :: gene_ids(:)
            !! Gene IDs to read
        integer(int32), intent(out)   :: ierr
            !! Error code (0 if successful)

        call deserialize_char_1D(gene_ids, filename, ierr)
    end subroutine load_gene_ids

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for deserialize_real_2D function
    subroutine load_expression_vectors(expression_vectors, filename, ierr)
        character(len=*), intent(in) :: filename
            !! Filename to read from
        real(real64), intent(out), contiguous    :: expression_vectors(:, :)
            !! Expression vectors to read
        integer(int32), intent(out)  :: ierr
            !! Error code (0 if successful)

        call deserialize_real_2D(expression_vectors, filename, ierr)
    end subroutine load_expression_vectors

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for deserialize_int_1D function
    subroutine load_gene_to_family(gene_to_fam, filename, ierr)
        character(len=*), intent(in) :: filename
            !! Filename to read from
        integer(int32), intent(out), contiguous :: gene_to_fam(:)
            !! Gene to family mapping to read
        integer(int32), intent(out) :: ierr
            !! Error code (0 if successful)

        call deserialize_int_1D(gene_to_fam, filename, ierr)
    end subroutine

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for deserialize_char_1D function
    subroutine load_family_ids(family_ids, filename, ierr)
        character(len=*), intent(in) :: filename
            !! Filename to read from
        CHARACTER(len=*), intent(out), contiguous :: family_ids(:)
            !! Family IDs to read
        integer(int32), intent(out) :: ierr
            !! Error code (0 if successful)

        call deserialize_char_1D(family_ids, filename, ierr)
    end subroutine

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for deserialize_real_2D function
    subroutine load_family_centroids(family_centroids, filename, ierr)
        character(len=*), intent(in) :: filename
            !! Filename to read from
        real(real64), intent(out), contiguous    :: family_centroids(:, :)
            !! Family centroids to read
        integer(int32), intent(out)  :: ierr
            !! Error code

        call deserialize_real_2D(family_centroids, filename, ierr)
    end subroutine

    !> AUTHOR_AARON_SCHROEDER
    !| Wrapper for deserialize_real_2D function
    subroutine load_shift_vectors(shift_vectors, filename, ierr)
        character(len=*), intent(in) :: filename
            !! Filename to read from
        real(real64), intent(out), contiguous    :: shift_vectors(:, :)
            !! Shift vectors to read
        integer(int32), intent(out)  :: ierr
            !! Error code

        call deserialize_real_2D(shift_vectors, filename, ierr)
    end subroutine

end module tox_data_read_write
