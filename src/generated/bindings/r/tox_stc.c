// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void calculate_density_radius_c(const double*, const int*, const int*, const double*, double*, int*);
void calculate_labels_as_density_c(const double*, const int*, const int*, const double*, const int*, const int*, const int*, double*, int*);
void identify_ensemble_seeds_c(const double*, const int*, const int*, const double*, const int*, const int*, const int*, int*, int*, unsigned char*, double*, const double*, int*);
void compute_ambient_density_stats_c(const double*, const int*, double*, double*, int*);
void merge_ensembles_c(const unsigned char*, const int*, const int*, const double*, unsigned char*, int*, int*);
void grow_ensemble_c(const double*, const int*, const int*, const unsigned char*, const double*, const int*, const int*, const double*, const double*, unsigned char*, int*);
void compute_ensemble_observable_c(const unsigned char*, const double*, const int*, const int*, const int*, double*, const int*, const int*, int*);
void accept_ensemble_c(const double*, const int*, const int*, const int*, const double*, unsigned char*, int*);
void obtain_ensembles_global_radius_c(const double*, const int*, const int*, const int*, const int*, const double*, const int*, const int*, const double*, const double*, const double*, const int*, const int*, unsigned char*, int*, double*, int*, int*);
void obtain_ensembles_seed_radius_c(const double*, const int*, const int*, const int*, const int*, const double*, const int*, const int*, const double*, const double*, const double*, const int*, const int*, unsigned char*, int*, double*, int*, int*);

SEXP calculate_density_radius_call(SEXP vectors, SEXP mean_to_other_vecs_dist_quant) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double _v_mean_to_other_vecs_dist_quant = Rf_asReal(mean_to_other_vecs_dist_quant);

    // outputs and work space
    double radius = 0;
    int ierr = 0;

    calculate_density_radius_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        &_v_mean_to_other_vecs_dist_quant,
        &radius,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(radius));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("radius"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP calculate_labels_as_density_call(SEXP vectors, SEXP r, SEXP dimension_order, SEXP kd_indices, SEXP n_tiles) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double _v_r = Rf_asReal(r);
    int _v_n_tiles = Rf_asInteger(n_tiles);

    // convert what Fortran cannot take from R directly
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");
    const int* _c_kd_indices = tox_int_in(kd_indices, "kd_indices");

    // outputs and work space
    SEXP label_densities = PROTECT(Rf_allocVector(REALSXP, n_vectors)); _nprot++;
    int ierr = 0;

    calculate_labels_as_density_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        &_v_r,
        _c_dimension_order,
        _c_kd_indices,
        &_v_n_tiles,
        REAL(label_densities),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, label_densities);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("label_densities"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP identify_ensemble_seeds_call(SEXP vectors, SEXP density_labels, SEXP dimension_order, SEXP kd_indices, SEXP k_seeding, SEXP coverage_quant) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_k_seeding = Rf_asInteger(k_seeding);
    double _v_coverage_quant = Rf_asReal(coverage_quant);

    // convert what Fortran cannot take from R directly
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");
    const int* _c_kd_indices = tox_int_in(kd_indices, "kd_indices");

    // outputs and work space
    SEXP sorted_perm = PROTECT(Rf_allocVector(INTSXP, n_vectors)); _nprot++;
    int n_seeds = 0;
    unsigned char* _c_seed_mask = tox_bool_alloc(n_vectors);
    SEXP seed_radii = PROTECT(Rf_allocVector(REALSXP, n_vectors)); _nprot++;
    int ierr = 0;

    identify_ensemble_seeds_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        REAL(density_labels),
        _c_dimension_order,
        _c_kd_indices,
        &_v_k_seeding,
        INTEGER(sorted_perm),
        &n_seeds,
        _c_seed_mask,
        REAL(seed_radii),
        &_v_coverage_quant,
        &ierr
    );

    // convert the outputs back
    SEXP seed_mask = PROTECT(tox_bool_out(_c_seed_mask, n_vectors)); _nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); _nprot++;
    SET_VECTOR_ELT(_out, 0, sorted_perm);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_seeds));
    SET_VECTOR_ELT(_out, 2, seed_mask);
    SET_VECTOR_ELT(_out, 3, seed_radii);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("sorted_perm"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_seeds"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("seed_mask"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("seed_radii"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_ambient_density_stats_call(SEXP density_labels) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_vectors = (int) Rf_length(density_labels);

    // outputs and work space
    double median_ambient = 0;
    double mad_ambient = 0;
    int ierr = 0;

    compute_ambient_density_stats_c(
        REAL(density_labels),
        &n_vectors,
        &median_ambient,
        &mad_ambient,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(median_ambient));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(mad_ambient));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("median_ambient"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("mad_ambient"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP merge_ensembles_call(SEXP raw_masks, SEXP min_overlap_coefficient) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_vectors = INTEGER(Rf_getAttrib(raw_masks, R_DimSymbol))[0];
    int n_seeds = INTEGER(Rf_getAttrib(raw_masks, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double _v_min_overlap_coefficient = Rf_asReal(min_overlap_coefficient);

    // convert what Fortran cannot take from R directly
    unsigned char* _c_raw_masks = tox_bool_in(raw_masks, "raw_masks");

    // outputs and work space
    unsigned char* _c_merged_masks = tox_bool_alloc(n_vectors * n_seeds);
    int n_ensembles = 0;
    int ierr = 0;

    merge_ensembles_c(
        _c_raw_masks,
        &n_vectors,
        &n_seeds,
        &_v_min_overlap_coefficient,
        _c_merged_masks,
        &n_ensembles,
        &ierr
    );

    // convert the outputs back
    SEXP merged_masks = PROTECT(tox_bool_out(_c_merged_masks, n_vectors * n_seeds)); _nprot++;
    { SEXP _dim_merged_masks = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_merged_masks)[0] = n_vectors; INTEGER(_dim_merged_masks)[1] = n_seeds; Rf_setAttrib(merged_masks, R_DimSymbol, _dim_merged_masks); UNPROTECT(1); }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, merged_masks);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_ensembles));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("merged_masks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_ensembles"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP grow_ensemble_call(SEXP vectors, SEXP ensemble_mask, SEXP r, SEXP dimension_order, SEXP kd_indices, SEXP density_labels, SEXP alpha_mad) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double _v_r = Rf_asReal(r);
    double _v_alpha_mad = Rf_asReal(alpha_mad);

    // convert what Fortran cannot take from R directly
    unsigned char* _c_ensemble_mask = tox_bool_in(ensemble_mask, "ensemble_mask");
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");
    const int* _c_kd_indices = tox_int_in(kd_indices, "kd_indices");

    // outputs and work space
    unsigned char* _c_grown_mask = tox_bool_alloc(n_vectors);
    int ierr = 0;

    grow_ensemble_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        _c_ensemble_mask,
        &_v_r,
        _c_dimension_order,
        _c_kd_indices,
        REAL(density_labels),
        &_v_alpha_mad,
        _c_grown_mask,
        &ierr
    );

    // convert the outputs back
    SEXP grown_mask = PROTECT(tox_bool_out(_c_grown_mask, n_vectors)); _nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, grown_mask);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("grown_mask"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_ensemble_observable_call(SEXP ensemble_mask, SEXP density_labels, SEXP candidate_count, SEXP current_iter, SEXP observables) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_vectors = (int) Rf_length(ensemble_mask);
    int n_history = INTEGER(Rf_getAttrib(observables, R_DimSymbol))[1];
    int n_observables = INTEGER(Rf_getAttrib(observables, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int _v_candidate_count = Rf_asInteger(candidate_count);
    int _v_current_iter = Rf_asInteger(current_iter);

    // copy what is modified in place, so the caller's stays intact
    SEXP _out_observables = PROTECT(Rf_duplicate(observables)); _nprot++;

    // convert what Fortran cannot take from R directly
    unsigned char* _c_ensemble_mask = tox_bool_in(ensemble_mask, "ensemble_mask");

    // outputs and work space
    int ierr = 0;

    compute_ensemble_observable_c(
        _c_ensemble_mask,
        REAL(density_labels),
        &n_vectors,
        &_v_candidate_count,
        &_v_current_iter,
        REAL(_out_observables),
        &n_history,
        &n_observables,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, _out_observables);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("observables"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP accept_ensemble_call(SEXP observables, SEXP current_iter, SEXP alpha_accept) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_history = INTEGER(Rf_getAttrib(observables, R_DimSymbol))[1];
    int n_observables = INTEGER(Rf_getAttrib(observables, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int _v_current_iter = Rf_asInteger(current_iter);
    double _v_alpha_accept = Rf_asReal(alpha_accept);

    // outputs and work space
    unsigned char is_accepted = 0;
    int ierr = 0;

    accept_ensemble_c(
        REAL(observables),
        &n_history,
        &n_observables,
        &_v_current_iter,
        &_v_alpha_accept,
        &is_accepted,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarLogical(is_accepted != 0));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_accepted"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP obtain_ensembles_global_radius_call(SEXP vectors, SEXP dimension_order, SEXP kd_indices, SEXP density_labels, SEXP seed_indices, SEXP r, SEXP alpha_mad, SEXP alpha_accept, SEXP t_observables, SEXP n_tiles) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];
    int n_seeds = (int) Rf_length(seed_indices);

    // scalar inputs, pulled from their length-1 vectors
    double _v_r = Rf_asReal(r);
    double _v_alpha_mad = Rf_asReal(alpha_mad);
    double _v_alpha_accept = Rf_asReal(alpha_accept);
    int _v_t_observables = Rf_asInteger(t_observables);
    int _v_n_tiles = Rf_asInteger(n_tiles);

    // convert what Fortran cannot take from R directly
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");
    const int* _c_kd_indices = tox_int_in(kd_indices, "kd_indices");
    const int* _c_seed_indices = tox_int_in(seed_indices, "seed_indices");

    // outputs and work space
    unsigned char* _c_ensemble_matrix = tox_bool_alloc(n_vectors * n_seeds);
    SEXP stop_reasons = PROTECT(Rf_allocVector(INTSXP, n_seeds)); _nprot++;
    double mad_ambient = 0;
    int n_ensembles = 0;
    int ierr = 0;

    obtain_ensembles_global_radius_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        _c_dimension_order,
        _c_kd_indices,
        REAL(density_labels),
        _c_seed_indices,
        &n_seeds,
        &_v_r,
        &_v_alpha_mad,
        &_v_alpha_accept,
        &_v_t_observables,
        &_v_n_tiles,
        _c_ensemble_matrix,
        INTEGER(stop_reasons),
        &mad_ambient,
        &n_ensembles,
        &ierr
    );

    // convert the outputs back
    SEXP ensemble_matrix = PROTECT(tox_bool_out(_c_ensemble_matrix, n_vectors * n_seeds)); _nprot++;
    { SEXP _dim_ensemble_matrix = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_ensemble_matrix)[0] = n_vectors; INTEGER(_dim_ensemble_matrix)[1] = n_seeds; Rf_setAttrib(ensemble_matrix, R_DimSymbol, _dim_ensemble_matrix); UNPROTECT(1); }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); _nprot++;
    SET_VECTOR_ELT(_out, 0, ensemble_matrix);
    SET_VECTOR_ELT(_out, 1, stop_reasons);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarReal(mad_ambient));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(n_ensembles));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ensemble_matrix"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("stop_reasons"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("mad_ambient"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("n_ensembles"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP obtain_ensembles_seed_radius_call(SEXP vectors, SEXP dimension_order, SEXP kd_indices, SEXP density_labels, SEXP seed_indices, SEXP seed_radii, SEXP alpha_mad, SEXP alpha_accept, SEXP t_observables, SEXP n_tiles) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];
    int n_seeds = (int) Rf_length(seed_indices);

    // scalar inputs, pulled from their length-1 vectors
    double _v_alpha_mad = Rf_asReal(alpha_mad);
    double _v_alpha_accept = Rf_asReal(alpha_accept);
    int _v_t_observables = Rf_asInteger(t_observables);
    int _v_n_tiles = Rf_asInteger(n_tiles);

    // convert what Fortran cannot take from R directly
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");
    const int* _c_kd_indices = tox_int_in(kd_indices, "kd_indices");
    const int* _c_seed_indices = tox_int_in(seed_indices, "seed_indices");

    // outputs and work space
    unsigned char* _c_ensemble_matrix = tox_bool_alloc(n_vectors * n_seeds);
    SEXP stop_reasons = PROTECT(Rf_allocVector(INTSXP, n_seeds)); _nprot++;
    double mad_ambient = 0;
    int n_ensembles = 0;
    int ierr = 0;

    obtain_ensembles_seed_radius_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        _c_dimension_order,
        _c_kd_indices,
        REAL(density_labels),
        _c_seed_indices,
        &n_seeds,
        REAL(seed_radii),
        &_v_alpha_mad,
        &_v_alpha_accept,
        &_v_t_observables,
        &_v_n_tiles,
        _c_ensemble_matrix,
        INTEGER(stop_reasons),
        &mad_ambient,
        &n_ensembles,
        &ierr
    );

    // convert the outputs back
    SEXP ensemble_matrix = PROTECT(tox_bool_out(_c_ensemble_matrix, n_vectors * n_seeds)); _nprot++;
    { SEXP _dim_ensemble_matrix = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_ensemble_matrix)[0] = n_vectors; INTEGER(_dim_ensemble_matrix)[1] = n_seeds; Rf_setAttrib(ensemble_matrix, R_DimSymbol, _dim_ensemble_matrix); UNPROTECT(1); }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); _nprot++;
    SET_VECTOR_ELT(_out, 0, ensemble_matrix);
    SET_VECTOR_ELT(_out, 1, stop_reasons);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarReal(mad_ambient));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(n_ensembles));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ensemble_matrix"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("stop_reasons"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("mad_ambient"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("n_ensembles"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
