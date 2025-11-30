!> Wrappers for serialization/deserialization of arrays
module tox_data_read_write
    use safeguard
    use iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err, is_err, is_ok
    use serialize_int, only: serialize_int_1d
    use int_deserialize_mod, only: deserialize_int_1D
    use serialize_real, only: serialize_real_2D
    use real_deserialize_mod, only: deserialize_real_2D
    use serialize_char, only: serialize_char_1D
    use char_deserialize_mod, only: deserialize_char_1D
    implicit none

contains
    !> Wrapper for serialize_char_1D function
subroutine save_gene_ids(gene_ids, filename, ierr)
    character(len=*), intent(in) :: gene_ids(:)
    !! Gene IDs to write
    character(len=*), intent(in) :: filename
    !! Filename to write to
    integer(int32), intent(out) :: ierr
    !! Error code (0 if successful)

    call serialize_char_1D(gene_ids, filename, ierr)
end subroutine

!> Wrapper for serialize_real_2D function
subroutine save_expression_vectors(expression_vectors, filename, ierr)
    character(len=*), intent(in) :: filename
    !! Filename to write to
    real(real64), intent(in) :: expression_vectors(:,:)
    !! Expression vectors to write
    integer(int32), intent(out) :: ierr
    !! Error code (0 if successful)

    call serialize_real_2D(expression_vectors, filename, ierr)
end subroutine

!> Wrapper for serialize_int_1D function
subroutine save_gene_to_family(gene_to_fam, filename, ierr)
    character(len=*), intent(in) :: filename
    !! Filename to write to
    integer(int32), intent(in) :: gene_to_fam(:)
    !! Gene to family mapping to write
    integer(int32), intent(out) :: ierr

    call serialize_int_1D(gene_to_fam, filename, ierr)
end subroutine

!> Wrapper for serialize_char_1D function
subroutine save_family_ids(family_ids, filename, ierr)
    character(len=*), intent(in) :: filename
    !! Filename to write to
    character(len=*), intent(in) :: family_ids(:)
    !! Family IDs to write
    integer(int32), intent(out) :: ierr
    !! Error code (0 if successful)

    call serialize_char_1D(family_ids, filename, ierr)
end subroutine

!> Wrapper for serialize_real_2D function
subroutine save_family_centroids(family_centroids, filename, ierr)
    character(len=*), intent(in) :: filename
    !! Filename to write to
    real(real64), intent(in) :: family_centroids(:,:)
    !! Family centroids to write
    integer(int32), intent(out) :: ierr
    !! Error code (0 if successful)

    call serialize_real_2D(family_centroids, filename, ierr)
end subroutine

!> Wrapper for serialize_real_2D function
subroutine save_shift_vectors(shift_vectors, filename, ierr)
    character(len=*), intent(in) :: filename
    !! Filename to write to
    real(real64), intent(in) :: shift_vectors(:,:)
    !! Shift vectors to write
    integer(int32), intent(out) :: ierr
    !! Error code (0 if successful)

    call serialize_real_2D(shift_vectors, filename, ierr)
end subroutine

!> Wrapper for deserialize_char_1D function
subroutine load_gene_ids(gene_ids, filename, ierr)
    character(len=*), intent(in)  :: filename
    !! Filename to read from
    character(len=*), intent(out) :: gene_ids(:)
    !! Gene IDs to read
    integer(int32), intent(out)   :: ierr
    !! Error code (0 if successful)

    call deserialize_char_1D(gene_ids, filename, ierr)
end subroutine load_gene_ids

!> Wrapper for deserialize_real_2D function
subroutine load_expression_vectors(expression_vectors, filename, ierr)
    character(len=*), intent(in) :: filename
    !! Filename to read from
    real(real64), intent(out)    :: expression_vectors(:,:)
    !! Expression vectors to read
    integer(int32), intent(out)  :: ierr
    !! Error code (0 if successful)

    call deserialize_real_2D(expression_vectors, filename, ierr)
end subroutine load_expression_vectors

!> Wrapper for deserialize_int_1D function
subroutine load_gene_to_family(gene_to_fam, filename, ierr)
    character(len=*), intent(in) :: filename
    !! Filename to read from
    integer(int32), intent(out) :: gene_to_fam(:)
    !! Gene to family mapping to read
    integer(int32), intent(out) :: ierr
    !! Error code (0 if successful)

    call deserialize_int_1D(gene_to_fam, filename, ierr)
end subroutine

!> Wrapper for deserialize_char_1D function
subroutine load_family_ids(family_ids, filename, ierr)
    character(len=*), intent(in) :: filename
    !! Filename to read from
    CHARACTER(len=*), intent(out) :: family_ids(:)
    !! Family IDs to read
    integer(int32), intent(out) :: ierr
    !! Error code (0 if successful)

    call deserialize_char_1D(family_ids, filename, ierr)
end subroutine

!> Wrapper for deserialize_real_2D function
subroutine load_family_centroids(family_centroids, filename, ierr)
    character(len=*), intent(in) :: filename
    !! Filename to read from
    real(real64), intent(out)    :: family_centroids(:,:)
    !! Family centroids to read
    integer(int32), intent(out)  :: ierr
    !! Error code

    call deserialize_real_2D(family_centroids, filename, ierr)
end subroutine

!> Wrapper for deserialize_real_2D function
subroutine load_shift_vectors(shift_vectors, filename, ierr)
    character(len=*), intent(in) :: filename
    !! Filename to read from
    real(real64), intent(out)    :: shift_vectors(:,:)
    !! Shift vectors to read
    integer(int32), intent(out)  :: ierr
    !! Error code

    call deserialize_real_2D(shift_vectors, filename, ierr)
end subroutine

end module tox_data_read_write