#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_archive(module)]]
!| Zip-archive backed persistence for TensorOmics data sets.
!|
!| Wraps libzip (via C bindings declared in the interface block below) to create/extract zip
!| archives whose members are [[tox_data_read_write(module)]]-serialized arrays, indexed by a
!| plain-text `manifest.txt` mapping logical keys (e.g. `gene_ids`, `expression`) to member
!| filenames. [[tox_data_archive(module):save_tox_data(subroutine)]] and
!| [[tox_data_archive(module):read_tox_data(subroutine)]] are the standard entry points for the
!| fixed TensorOmics data set schema; `create_zip_archive`/`extract_zip_archive` and the
!| `*_manifest*` routines below are the generic key/filename building blocks they are built on.
module tox_data_archive_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string, c_char_2d_as_string, string_as_c_char_2d
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE
    private

    public :: create_zip_archive_c
    public :: save_tox_data_c
    public :: get_tox_data_dims_c
    public :: read_tox_data_into_c

contains

    !> summary: C-wrapper for [[tox_data_archive(module):create_zip_archive(subroutine)]]
    subroutine create_zip_archive_c(&
            zip_filename,&
            zip_filename_strlen,&
            keys,&
            keys_strlen,&
            n_keys_elements,&
            filenames,&
            filenames_strlen,&
            n_filenames_elements,&
            ierr&
        ) bind(C, name="create_zip_archive_c")
        use tox_data_archive, only: create_zip_archive

        integer(c_int), intent(in), target :: zip_filename_strlen
            !! length of the strings in `zip_filename`
        integer(c_int), intent(in), target :: keys_strlen
            !! length of the strings in `keys`
        integer(c_int), intent(in), target :: n_keys_elements
            !! number of elements in `keys`
        integer(c_int), intent(in), target :: filenames_strlen
            !! length of the strings in `filenames`
        integer(c_int), intent(in), target :: n_filenames_elements
            !! number of elements in `filenames`
        character(len=1, kind=c_char), dimension(zip_filename_strlen), intent(in), target :: zip_filename
            !! Name of the zip file to create
        character(len=1, kind=c_char), dimension(keys_strlen, n_keys_elements), intent(in), target :: keys
            !! Array of keys for manifest entries
        character(len=1, kind=c_char), dimension(filenames_strlen, n_filenames_elements), intent(in), target :: filenames
            !! Array of filenames to add to zip
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable :: zip_filename_f
        character(len=:), allocatable, dimension(:) :: keys_f
        character(len=:), allocatable, dimension(:) :: filenames_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(zip_filename_strlen)
        M_CHECK_NON_NULL(keys_strlen)
        M_CHECK_NON_NULL(n_keys_elements)
        M_CHECK_NON_NULL(filenames_strlen)
        M_CHECK_NON_NULL(n_filenames_elements)
        M_CHECK_ARRAY_NON_NULL(zip_filename, zip_filename_strlen)
        M_CHECK_ARRAY_NON_NULL(keys, keys_strlen * n_keys_elements)
        M_CHECK_ARRAY_NON_NULL(filenames, filenames_strlen * n_filenames_elements)

        call c_char_1d_as_string(zip_filename, zip_filename_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(keys, keys_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(filenames, filenames_f, ierr)
        if (is_err(ierr)) return

        call create_zip_archive(&
            zip_filename = zip_filename_f,&
            keys = keys_f,&
            filenames = filenames_f,&
            ierr = ierr&
        )
    end subroutine create_zip_archive_c

    !> summary: C-wrapper for [[tox_data_archive(module):save_tox_data(subroutine)]]
    subroutine save_tox_data_c(&
            zip_filename,&
            zip_filename_strlen,&
            ierr,&
            gene_ids,&
            gene_ids_strlen,&
            n_gene_ids_elements,&
            gene_ids_file,&
            gene_ids_file_strlen,&
            expression,&
            n_expression_elements_dim_1,&
            n_expression_elements_dim_2,&
            expression_file,&
            expression_file_strlen,&
            gene_to_family,&
            n_gene_to_family_elements,&
            gene_to_family_file,&
            gene_to_family_file_strlen,&
            family_ids,&
            family_ids_strlen,&
            n_family_ids_elements,&
            family_ids_file,&
            family_ids_file_strlen,&
            family_centroids,&
            n_family_centroids_elements_dim_1,&
            n_family_centroids_elements_dim_2,&
            family_centroids_file,&
            family_centroids_file_strlen,&
            shift_vectors,&
            n_shift_vectors_elements_dim_1,&
            n_shift_vectors_elements_dim_2,&
            shift_vectors_file,&
            shift_vectors_file_strlen&
        ) bind(C, name="save_tox_data_c")
        use tox_data_archive, only: save_tox_data

        integer(c_int), intent(in), target :: zip_filename_strlen
            !! length of the strings in `zip_filename`
        integer(c_int), intent(in), target :: gene_ids_strlen
            !! length of the strings in `gene_ids`
        integer(c_int), intent(in), target :: n_gene_ids_elements
            !! number of elements in `gene_ids`
        integer(c_int), intent(in), target :: gene_ids_file_strlen
            !! length of the strings in `gene_ids_file`
        integer(c_int), intent(in), target :: n_expression_elements_dim_1
            !! 1. dimension of `expression`
        integer(c_int), intent(in), target :: n_expression_elements_dim_2
            !! 2. dimension of `expression`
        integer(c_int), intent(in), target :: expression_file_strlen
            !! length of the strings in `expression_file`
        integer(c_int), intent(in), target :: n_gene_to_family_elements
            !! number of elements in `gene_to_family`
        integer(c_int), intent(in), target :: gene_to_family_file_strlen
            !! length of the strings in `gene_to_family_file`
        integer(c_int), intent(in), target :: family_ids_strlen
            !! length of the strings in `family_ids`
        integer(c_int), intent(in), target :: n_family_ids_elements
            !! number of elements in `family_ids`
        integer(c_int), intent(in), target :: family_ids_file_strlen
            !! length of the strings in `family_ids_file`
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_1
            !! 1. dimension of `family_centroids`
        integer(c_int), intent(in), target :: n_family_centroids_elements_dim_2
            !! 2. dimension of `family_centroids`
        integer(c_int), intent(in), target :: family_centroids_file_strlen
            !! length of the strings in `family_centroids_file`
        integer(c_int), intent(in), target :: n_shift_vectors_elements_dim_1
            !! 1. dimension of `shift_vectors`
        integer(c_int), intent(in), target :: n_shift_vectors_elements_dim_2
            !! 2. dimension of `shift_vectors`
        integer(c_int), intent(in), target :: shift_vectors_file_strlen
            !! length of the strings in `shift_vectors_file`
        character(len=1, kind=c_char), dimension(zip_filename_strlen), intent(in), target :: zip_filename
            !! Zip filename
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=1, kind=c_char), dimension(gene_ids_strlen, n_gene_ids_elements), intent(in), optional :: gene_ids
            !! Gene ids array, will be saved if provided
        character(len=1, kind=c_char), dimension(gene_ids_file_strlen), intent(in), optional :: gene_ids_file
            !! Name of the gene ids file
        real(c_double), dimension(n_expression_elements_dim_1, n_expression_elements_dim_2), intent(in), optional :: expression
            !! Expression vectors array, will be saved if provided
        character(len=1, kind=c_char), dimension(expression_file_strlen), intent(in), optional :: expression_file
            !! Name of the expression file
        integer(c_int), dimension(n_gene_to_family_elements), intent(in), optional :: gene_to_family
            !! Gene to family mapping array, will be saved if provided
        character(len=1, kind=c_char), dimension(gene_to_family_file_strlen), intent(in), optional :: gene_to_family_file
            !! Name of the gene to family mapping file
        character(len=1, kind=c_char), dimension(family_ids_strlen, n_family_ids_elements), intent(in), optional :: family_ids
            !! Family ids array, will be saved if provided
        character(len=1, kind=c_char), dimension(family_ids_file_strlen), intent(in), optional :: family_ids_file
            !! Name of the family ids file
        real(c_double), dimension(n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2), intent(in), optional :: family_centroids
            !! Family centroids array, will be saved if provided
        character(len=1, kind=c_char), dimension(family_centroids_file_strlen), intent(in), optional :: family_centroids_file
            !! Name of the family centroids file
        real(c_double), dimension(n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2), intent(in), optional :: shift_vectors
            !! Shift vectors array, will be saved if provided
        character(len=1, kind=c_char), dimension(shift_vectors_file_strlen), intent(in), optional :: shift_vectors_file
            !! Name of the shift vectors file
        character(len=:), allocatable :: zip_filename_f
        character(len=:), allocatable, dimension(:) :: gene_ids_f
        character(len=:), allocatable :: gene_ids_file_f
        character(len=:), allocatable :: expression_file_f
        character(len=:), allocatable :: gene_to_family_file_f
        character(len=:), allocatable, dimension(:) :: family_ids_f
        character(len=:), allocatable :: family_ids_file_f
        character(len=:), allocatable :: family_centroids_file_f
        character(len=:), allocatable :: shift_vectors_file_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(zip_filename_strlen)
        M_CHECK_NON_NULL(gene_ids_strlen)
        M_CHECK_NON_NULL(n_gene_ids_elements)
        M_CHECK_NON_NULL(gene_ids_file_strlen)
        M_CHECK_NON_NULL(n_expression_elements_dim_1)
        M_CHECK_NON_NULL(n_expression_elements_dim_2)
        M_CHECK_NON_NULL(expression_file_strlen)
        M_CHECK_NON_NULL(n_gene_to_family_elements)
        M_CHECK_NON_NULL(gene_to_family_file_strlen)
        M_CHECK_NON_NULL(family_ids_strlen)
        M_CHECK_NON_NULL(n_family_ids_elements)
        M_CHECK_NON_NULL(family_ids_file_strlen)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_1)
        M_CHECK_NON_NULL(n_family_centroids_elements_dim_2)
        M_CHECK_NON_NULL(family_centroids_file_strlen)
        M_CHECK_NON_NULL(n_shift_vectors_elements_dim_1)
        M_CHECK_NON_NULL(n_shift_vectors_elements_dim_2)
        M_CHECK_NON_NULL(shift_vectors_file_strlen)
        M_CHECK_ARRAY_NON_NULL(zip_filename, zip_filename_strlen)

        call c_char_1d_as_string(zip_filename, zip_filename_f, ierr)
        if (is_err(ierr)) return
        if (present(gene_ids)) then
            call c_char_2d_as_string(gene_ids, gene_ids_f, ierr)
            if (is_err(ierr)) return
        end if
        if (present(gene_ids_file)) then
            call c_char_1d_as_string(gene_ids_file, gene_ids_file_f, ierr)
            if (is_err(ierr)) return
        end if
        if (present(expression_file)) then
            call c_char_1d_as_string(expression_file, expression_file_f, ierr)
            if (is_err(ierr)) return
        end if
        if (present(gene_to_family_file)) then
            call c_char_1d_as_string(gene_to_family_file, gene_to_family_file_f, ierr)
            if (is_err(ierr)) return
        end if
        if (present(family_ids)) then
            call c_char_2d_as_string(family_ids, family_ids_f, ierr)
            if (is_err(ierr)) return
        end if
        if (present(family_ids_file)) then
            call c_char_1d_as_string(family_ids_file, family_ids_file_f, ierr)
            if (is_err(ierr)) return
        end if
        if (present(family_centroids_file)) then
            call c_char_1d_as_string(family_centroids_file, family_centroids_file_f, ierr)
            if (is_err(ierr)) return
        end if
        if (present(shift_vectors_file)) then
            call c_char_1d_as_string(shift_vectors_file, shift_vectors_file_f, ierr)
            if (is_err(ierr)) return
        end if

        call save_tox_data(&
            zip_filename = zip_filename_f,&
            ierr = ierr,&
            gene_ids = gene_ids_f,&
            gene_ids_file = gene_ids_file_f,&
            expression = expression,&
            expression_file = expression_file_f,&
            gene_to_family = gene_to_family,&
            gene_to_family_file = gene_to_family_file_f,&
            family_ids = family_ids_f,&
            family_ids_file = family_ids_file_f,&
            family_centroids = family_centroids,&
            family_centroids_file = family_centroids_file_f,&
            shift_vectors = shift_vectors,&
            shift_vectors_file = shift_vectors_file_f&
        )
    end subroutine save_tox_data_c

    !> summary: C-wrapper for [[tox_data_archive(module):get_tox_data_dims(subroutine)]]
    !| Each count (and each string length) is 0 when the corresponding member is absent, so a
    !| caller can size all six output buffers up front. Character members report both an element
    !| count and a per-element string length. Pairs with
    !| [[tox_data_archive(module):read_tox_data_into(subroutine)]].
    subroutine get_tox_data_dims_c(&
            zip_filename,&
            zip_filename_strlen,&
            n_gene_ids,&
            gene_id_len,&
            n_expression_rows,&
            n_expression_cols,&
            n_gene_to_family,&
            n_family_ids,&
            family_id_len,&
            n_family_centroids_rows,&
            n_family_centroids_cols,&
            n_shift_vectors_rows,&
            n_shift_vectors_cols,&
            ierr&
        ) bind(C, name="get_tox_data_dims_c")
        use tox_data_archive, only: get_tox_data_dims

        integer(c_int), intent(in), target :: zip_filename_strlen
            !! length of the strings in `zip_filename`
        character(len=1, kind=c_char), dimension(zip_filename_strlen), intent(in), target :: zip_filename
            !! Name of the zip file
        integer(c_int), intent(out), target :: n_gene_ids
            !! Number of gene ids, 0 if absent
        integer(c_int), intent(out), target :: gene_id_len
            !! String length of each gene id, 0 if absent
        integer(c_int), intent(out), target :: n_expression_rows
            !! Rows (samples) of the expression matrix, 0 if absent
        integer(c_int), intent(out), target :: n_expression_cols
            !! Columns (genes) of the expression matrix, 0 if absent
        integer(c_int), intent(out), target :: n_gene_to_family
            !! Number of gene-to-family entries, 0 if absent
        integer(c_int), intent(out), target :: n_family_ids
            !! Number of family ids, 0 if absent
        integer(c_int), intent(out), target :: family_id_len
            !! String length of each family id, 0 if absent
        integer(c_int), intent(out), target :: n_family_centroids_rows
            !! Rows (samples) of the family centroids matrix, 0 if absent
        integer(c_int), intent(out), target :: n_family_centroids_cols
            !! Columns (families) of the family centroids matrix, 0 if absent
        integer(c_int), intent(out), target :: n_shift_vectors_rows
            !! Rows of the shift vectors matrix, 0 if absent
        integer(c_int), intent(out), target :: n_shift_vectors_cols
            !! Columns of the shift vectors matrix, 0 if absent
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable :: zip_filename_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(zip_filename_strlen)
        M_CHECK_NON_NULL(n_gene_ids)
        M_CHECK_NON_NULL(gene_id_len)
        M_CHECK_NON_NULL(n_expression_rows)
        M_CHECK_NON_NULL(n_expression_cols)
        M_CHECK_NON_NULL(n_gene_to_family)
        M_CHECK_NON_NULL(n_family_ids)
        M_CHECK_NON_NULL(family_id_len)
        M_CHECK_NON_NULL(n_family_centroids_rows)
        M_CHECK_NON_NULL(n_family_centroids_cols)
        M_CHECK_NON_NULL(n_shift_vectors_rows)
        M_CHECK_NON_NULL(n_shift_vectors_cols)
        M_CHECK_ARRAY_NON_NULL(zip_filename, zip_filename_strlen)

        call c_char_1d_as_string(zip_filename, zip_filename_f, ierr)
        if (is_err(ierr)) return

        call get_tox_data_dims(&
            zip_filename = zip_filename_f,&
            n_gene_ids = n_gene_ids,&
            gene_id_len = gene_id_len,&
            n_expression_rows = n_expression_rows,&
            n_expression_cols = n_expression_cols,&
            n_gene_to_family = n_gene_to_family,&
            n_family_ids = n_family_ids,&
            family_id_len = family_id_len,&
            n_family_centroids_rows = n_family_centroids_rows,&
            n_family_centroids_cols = n_family_centroids_cols,&
            n_shift_vectors_rows = n_shift_vectors_rows,&
            n_shift_vectors_cols = n_shift_vectors_cols,&
            ierr = ierr&
        )
    end subroutine get_tox_data_dims_c

    !> summary: C-wrapper for [[tox_data_archive(module):read_tox_data_into(subroutine)]]
    !| Fills every buffer from the archive; size them from
    !| [[tox_data_archive(module):get_tox_data_dims(subroutine)]] first. A member that is absent
    !| has a zero extent and is left untouched.
    subroutine read_tox_data_into_c(&
            zip_filename,&
            zip_filename_strlen,&
            n_gene_ids,&
            gene_id_len,&
            gene_ids,&
            n_expression_rows,&
            n_expression_cols,&
            expression,&
            n_gene_to_family,&
            gene_to_family,&
            n_family_ids,&
            family_id_len,&
            family_ids,&
            n_family_centroids_rows,&
            n_family_centroids_cols,&
            family_centroids,&
            n_shift_vectors_rows,&
            n_shift_vectors_cols,&
            shift_vectors,&
            ierr&
        ) bind(C, name="read_tox_data_into_c")
        use tox_data_archive, only: read_tox_data_into

        integer(c_int), intent(in), target :: zip_filename_strlen
            !! length of the strings in `zip_filename`
        integer(c_int), intent(in), target :: n_gene_ids
            !! Number of gene ids.
            !! It is *VERY IMPORTANT* to compute this argument from the `n_gene_ids` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        integer(c_int), intent(in), target :: gene_id_len
            !! String length of each gene id.
            !! It is *VERY IMPORTANT* to compute this argument from the `gene_id_len` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        integer(c_int), intent(in), target :: n_expression_rows
            !! Rows (samples) of the expression matrix.
            !! It is *VERY IMPORTANT* to compute this argument from the `n_expression_rows` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        integer(c_int), intent(in), target :: n_expression_cols
            !! Columns (genes) of the expression matrix.
            !! It is *VERY IMPORTANT* to compute this argument from the `n_expression_cols` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        integer(c_int), intent(in), target :: n_gene_to_family
            !! Number of gene-to-family entries.
            !! It is *VERY IMPORTANT* to compute this argument from the `n_gene_to_family` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        integer(c_int), intent(in), target :: n_family_ids
            !! Number of family ids.
            !! It is *VERY IMPORTANT* to compute this argument from the `n_family_ids` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        integer(c_int), intent(in), target :: family_id_len
            !! String length of each family id.
            !! It is *VERY IMPORTANT* to compute this argument from the `family_id_len` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        integer(c_int), intent(in), target :: n_family_centroids_rows
            !! Rows (samples) of the family centroids matrix.
            !! It is *VERY IMPORTANT* to compute this argument from the `n_family_centroids_rows` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        integer(c_int), intent(in), target :: n_family_centroids_cols
            !! Columns (families) of the family centroids matrix.
            !! It is *VERY IMPORTANT* to compute this argument from the `n_family_centroids_cols` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        integer(c_int), intent(in), target :: n_shift_vectors_rows
            !! Rows of the shift vectors matrix.
            !! It is *VERY IMPORTANT* to compute this argument from the `n_shift_vectors_rows` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        integer(c_int), intent(in), target :: n_shift_vectors_cols
            !! Columns of the shift vectors matrix.
            !! It is *VERY IMPORTANT* to compute this argument from the `n_shift_vectors_cols` output produced by [[tox_data_archive(module):get_tox_data_dims]].
        character(len=1, kind=c_char), dimension(zip_filename_strlen), intent(in), target :: zip_filename
            !! Name of the zip file
        character(len=1, kind=c_char), dimension(gene_id_len, n_gene_ids), intent(out), target :: gene_ids
            !! Gene ids
        real(c_double), dimension(n_expression_rows, n_expression_cols), intent(out), target :: expression
            !! Expression vectors
        integer(c_int), dimension(n_gene_to_family), intent(out), target :: gene_to_family
            !! Gene to family mapping
        character(len=1, kind=c_char), dimension(family_id_len, n_family_ids), intent(out), target :: family_ids
            !! Family ids
        real(c_double), dimension(n_family_centroids_rows, n_family_centroids_cols), intent(out), target :: family_centroids
            !! Family centroids
        real(c_double), dimension(n_shift_vectors_rows, n_shift_vectors_cols), intent(out), target :: shift_vectors
            !! Shift vectors
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable :: zip_filename_f
        character(len=gene_id_len), allocatable, dimension(:) :: gene_ids_f
        character(len=family_id_len), allocatable, dimension(:) :: family_ids_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(zip_filename_strlen)
        M_CHECK_NON_NULL(n_gene_ids)
        M_CHECK_NON_NULL(gene_id_len)
        M_CHECK_NON_NULL(n_expression_rows)
        M_CHECK_NON_NULL(n_expression_cols)
        M_CHECK_NON_NULL(n_gene_to_family)
        M_CHECK_NON_NULL(n_family_ids)
        M_CHECK_NON_NULL(family_id_len)
        M_CHECK_NON_NULL(n_family_centroids_rows)
        M_CHECK_NON_NULL(n_family_centroids_cols)
        M_CHECK_NON_NULL(n_shift_vectors_rows)
        M_CHECK_NON_NULL(n_shift_vectors_cols)
        M_CHECK_ARRAY_NON_NULL(zip_filename, zip_filename_strlen)
        M_CHECK_ARRAY_NON_NULL(gene_ids, gene_id_len * n_gene_ids)
        M_CHECK_ARRAY_NON_NULL(expression, n_expression_rows * n_expression_cols)
        M_CHECK_ARRAY_NON_NULL(gene_to_family, n_gene_to_family)
        M_CHECK_ARRAY_NON_NULL(family_ids, family_id_len * n_family_ids)
        M_CHECK_ARRAY_NON_NULL(family_centroids, n_family_centroids_rows * n_family_centroids_cols)
        M_CHECK_ARRAY_NON_NULL(shift_vectors, n_shift_vectors_rows * n_shift_vectors_cols)

        call c_char_1d_as_string(zip_filename, zip_filename_f, ierr)
        if (is_err(ierr)) return
        M_ALLOCATE(gene_ids_f(n_gene_ids))
        M_ALLOCATE(family_ids_f(n_family_ids))

        call read_tox_data_into(&
            zip_filename = zip_filename_f,&
            n_gene_ids = n_gene_ids,&
            gene_id_len = gene_id_len,&
            gene_ids = gene_ids_f,&
            n_expression_rows = n_expression_rows,&
            n_expression_cols = n_expression_cols,&
            expression = expression,&
            n_gene_to_family = n_gene_to_family,&
            gene_to_family = gene_to_family,&
            n_family_ids = n_family_ids,&
            family_id_len = family_id_len,&
            family_ids = family_ids_f,&
            n_family_centroids_rows = n_family_centroids_rows,&
            n_family_centroids_cols = n_family_centroids_cols,&
            family_centroids = family_centroids,&
            n_shift_vectors_rows = n_shift_vectors_rows,&
            n_shift_vectors_cols = n_shift_vectors_cols,&
            shift_vectors = shift_vectors,&
            ierr = ierr&
        )

        call string_as_c_char_2d(gene_ids_f, gene_ids)
        call string_as_c_char_2d(family_ids_f, family_ids)
    end subroutine read_tox_data_into_c

end module tox_data_archive_c
#endif
