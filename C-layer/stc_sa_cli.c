/* stc_sa_cli: simulated annealing over Shape Truthful Clustering (STC) parameters, scored
 * against known-manifold benchmark datasets. See misc/mod_STC.md, "Objective Function and
 * Simulated Annealing for Tangent-Space Shape Truthful Clustering" and "Simulated annealing
 * CLI in C". This file is wiring only -- it never reimplements any part of STC itself, only
 * calls the same generated `_c` bindings `stc_cli.c` calls, in the same order.
 *
 * The tuned parameter vector phi = (k_min, chordal_dist_max_as_prcnt_of_range, d_max, G_max,
 * RMSE_change_max, radius_percentile, min_stable_iterations) -- see mod_STC.md's "Simulated
 * annealing" section. `k_density` has been abolished as a separate parameter -- `k_min` is
 * used everywhere seeding, density estimation, and growth radius need a neighborhood size.
 * `estimate_stc_parameters`'s own `density_quantile` output has no STC input to feed back
 * into, so it is not part of phi; `RMSE_change_max`/`radius_percentile`/`min_stable_iterations`
 * have no estimator output at all, so the estimator-seeded chain starts them at the midpoint
 * of their own phi range instead. Every other STC/reconciliation parameter is fixed for the
 * whole run via CLI flags, exactly as in `stc_cli`.
 *
 * --store-iteration-shatter-results additionally writes, for every (chain, dataset, iteration),
 * the same JSON+CSV/TSV artifacts `stc_cli` writes for a standalone run (never the HTML report,
 * which would re-embed the multi-megabyte D3 bundle every iteration -- see mod_STC.md's
 * "Full per-iteration STC output") into
 * `--output-dir/{dataset_stem}_chain{chain_id}/iter_{iteration:06d}/`.
 */

#include <argp.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <omp.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <cjson/cJSON.h>

/* ==== the STC `_c` symbols this file calls =========================================
 * Same convention as `stc_cli.c`: no shared header exists, every consumer declares its own
 * prototypes matching the exact `_c.F90` signatures under src/generated/bindings/c/. */

extern void build_kd_index_c(const double *points, const int *n_dimensions, const int *n_points,
                             int *kd_indices, const int *dimension_order, int *ierr);

extern void seeds_c(const double *vectors, const int *n_dimensions, const int *n_vectors,
                    const int *kd_indices, const int *dimension_order, const int *k_min,
                    const double *bandwidth_percentile, const double *exclusion_radius_percentile,
                    unsigned char *is_seed_mask, int *ierr);

extern void estimate_stc_parameters_c(const double *vectors, const int *n_dimensions, const int *n_vectors,
                                      const int *kd_indices, const int *dimension_order,
                                      const int *k_min, const double *bandwidth_percentile,
                                      const int *n_anchors, const double *seed_max_set_size,
                                      const double *quantile_pairwise_ea_comparison,
                                      double *estimated_k_min,
                                      double *estimated_density_quantile,
                                      double *estimated_chordal_dist_max_as_prcnt_of_range,
                                      double *estimated_G_max, double *estimated_d_max, int *ierr);

extern void ensemble_identification_merged_c(
    const double *vectors, const int *n_dimensions, const int *n_vectors, const int *kd_indices,
    const int *dimension_order, const unsigned char *seed_selection_mask, const int *n_selected_seed,
    const int *k_min, const double *chordal_dist_max_as_prcnt_of_range, const int *d_max,
    const double *G_max, const double *RMSE_change_max, const double *f_max,
    const int *min_stable_iterations, const double *radius_percentile, const int *o,
    unsigned char *ensemble_masks, int *ensemble_stop_reason, double *ensemble_growth_radii,
    double *ensemble_U_history, double *ensemble_S_history, int *ensemble_d_history,
    double *ensemble_G_history, double *ensemble_mu_history, int *ensemble_k_history,
    unsigned char *ensemble_accepted_history, int *ensemble_member_added_at_step,
    unsigned char *ensemble_low_confidence_masks, double *ensemble_U_first, int *ensemble_d_first,
    int *ierr);

extern void ensemble_reconciliation_c(const unsigned char *ensemble_masks,
                                      const int *ensemble_stop_reason, const int *n_dimensions,
                                      const int *n_vectors, const int *n_ensembles,
                                      const double *ensemble_U_history, const int *ensemble_d_history,
                                      const double *ensemble_S_history, const double *ensemble_mu_history,
                                      const double *ensemble_G_history, const int *ensemble_k_history,
                                      const unsigned char *ensemble_accepted_history, const int *o,
                                      const char *mode,
                                      const double *min_overlap_coefficient,
                                      const unsigned char *report_overlap_coefficient,
                                      const unsigned char *allowed_stop_reasons,
                                      const int *d_min, const int *d_max, const double *var_explained_min,
                                      const int *max_group_size, int *super_ensembles,
                                      int *n_super_ensembles, double *super_ensembles_overlap_coefficient,
                                      unsigned char *eligible, unsigned char *eligible_by_stop_condition,
                                      unsigned char *eligible_by_dimension, unsigned char *eligible_by_var_explained,
                                      int *ierr);

/* The ensemble's real *final* accepted state (tangent basis, center, ...) -- distinct from
 * `ensemble_U_first`/`ensemble_d_first` above, which is the *first*-accepted, still-small
 * snapshot from early in growth. L_A must compare against the final tangent basis, the same
 * ensemble state L_D's center is computed from (see compute_losses) -- using U_first there
 * was a real bug: it let a bloated, structure-agnostic ensemble (huge k_min swallowing most of
 * a dataset) score a deceptively good L_A from its small early tangent estimate, before it
 * grew past the point where that estimate was still locally meaningful. */
extern void ensemble_final_observable_c(
    const int *n_dimensions, const int *o, const int *n_ensembles, const double *ensemble_U_history,
    const int *ensemble_d_history, const double *ensemble_S_history, const double *ensemble_mu_history,
    const double *ensemble_G_history, const int *ensemble_k_history,
    const unsigned char *ensemble_accepted_history, double *ensemble_U_final, int *ensemble_d_final,
    double *ensemble_S_final, double *ensemble_mu_final, double *ensemble_G_final, int *ensemble_k_final,
    unsigned char *ensemble_has_final, int *ensemble_final_index, int *ierr);

/* For --store-iteration-shatter-results: identical signatures to what `stc_cli.c` calls for a
 * standalone run (JSON + CSV/TSV companions only -- never the HTML report, see the header
 * comment above). */

extern void serialize_stc_results_as_json_c(
    const char *filename, const int *filename_strlen, const int *n_dimensions, const int *n_vectors,
    const int *n_selected_seed, const int *o, const int *max_group_size, const int *n_super_ensembles,
    const double *vectors, const char *dim_names, const int *dim_names_strlen,
    const unsigned char *seed_selection_mask, const unsigned char *ensemble_masks,
    const int *ensemble_stop_reason, const double *ensemble_growth_radii, const double *ensemble_U_history,
    const double *ensemble_S_history, const int *ensemble_d_history, const double *ensemble_G_history,
    const double *ensemble_mu_history, const int *ensemble_k_history,
    const unsigned char *ensemble_accepted_history, const int *ensemble_member_added_at_step,
    const unsigned char *ensemble_low_confidence_masks,
    const double *ensemble_U_first, const int *ensemble_d_first,
    const int *super_ensembles, const int *k_min,
    const double *chordal_dist_max_as_prcnt_of_range, const int *d_max,
    const double *G_max, const double *RMSE_change_max, const double *f_max,
    const int *min_stable_iterations, const double *radius_percentile,
    const double *exclusion_radius_percentile, const double *bandwidth_percentile,
    const char *reconciliation_mode, const double *min_overlap_coefficient,
    const unsigned char *allowed_stop_reasons, const int *filter_dim_min, const int *filter_dim_max,
    const double *filter_var_explained_min, const unsigned char *ensemble_eligible,
    const unsigned char *ensemble_eligible_by_stop_condition, const unsigned char *ensemble_eligible_by_dimension,
    const unsigned char *ensemble_eligible_by_var_explained, const int *estimated_k_min,
    const double *estimated_density_quantile,
    const double *estimated_chordal_dist_max_as_prcnt_of_range, const double *estimated_G_max,
    const int *estimated_d_max, int *ierr);

extern void serialize_stc_points_as_csv_c(const char *filename, const int *filename_strlen,
                                          const int *n_vectors, const int *n_selected_seed,
                                          const int *max_group_size, const int *n_super_ensembles,
                                          const unsigned char *seed_selection_mask,
                                          const unsigned char *ensemble_masks,
                                          const unsigned char *ensemble_low_confidence_masks,
                                          const int *super_ensembles, int *ierr);

extern void serialize_stc_ensemble_overlap_as_csv_c(const char *filename, const int *filename_strlen,
                                                    const int *n_vectors, const int *n_selected_seed,
                                                    const unsigned char *ensemble_masks, int *ierr);

extern void serialize_stc_super_ensembles_as_tsv_c(const char *filename, const int *filename_strlen,
                                                   const int *max_group_size, const int *n_super_ensembles,
                                                   const int *super_ensembles, int *ierr);

/* ==== small helpers ================================================================ */

static void die(const char *fmt, ...) {
    va_list ap;
    fprintf(stderr, "stc_sa_cli: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    exit(1);
}

static void check_ierr(int ierr, const char *what) {
    if (ierr != 0) die("%s failed (ierr=%d)", what, ierr);
}

static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (p == NULL && n > 0) die("out of memory (requested %zu bytes)", n);
    return p;
}

static void *xcalloc(size_t n, size_t sz) {
    void *p = calloc(n, sz);
    if (p == NULL && n > 0) die("out of memory (requested %zu x %zu bytes)", n, sz);
    return p;
}

/* Fills `dst[0..width)` with `src`'s bytes, null-padding the rest -- the layout
 * `tox_conversions::c_char_1d_as_string` (used by every generated `_c` wrapper that takes a
 * Fortran string) expects: it scans for a null byte within the slot to find the string's real
 * length, falling back to the whole slot if none is found. Same helper as `stc_cli.c`. */
static void fill_padded_string(char *dst, size_t width, const char *src) {
    size_t len = strlen(src);
    if (len > width) len = width;
    memset(dst, 0, width);
    memcpy(dst, src, len);
}

static char *join_path(const char *dir, const char *name) {
    size_t len = strlen(dir) + 1 + strlen(name) + 1;
    char *out = xmalloc(len);
    snprintf(out, len, "%s/%s", dir, name);
    return out;
}

/* Fortran `logical(c_bool)` layout used by every `_c` mask argument. */
typedef unsigned char c_bool_t;

/* ==== RNG: splitmix64-seeded xoshiro256** ==========================================
 * Self-contained (no libc rand()/drand48() dependency -- those are not guaranteed portable
 * or reentrant across platforms/threads), deterministic from an explicit seed, one
 * independent stream per chain (see rng_seed_for_chain below). Public-domain algorithm
 * (Blackman & Vigna). */

typedef struct { uint64_t s[4]; } rng_state;

static uint64_t splitmix64_next(uint64_t *state) {
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static void rng_seed(rng_state *r, uint64_t seed) {
    uint64_t sm = seed;
    r->s[0] = splitmix64_next(&sm);
    r->s[1] = splitmix64_next(&sm);
    r->s[2] = splitmix64_next(&sm);
    r->s[3] = splitmix64_next(&sm);
}

/* One independent stream per chain, deterministic from the CLI's --seed and the chain index. */
static void rng_seed_for_chain(rng_state *r, uint64_t base_seed, int chain_id) {
    rng_seed(r, base_seed + 0x9E3779B97F4A7C15ULL * (uint64_t)(chain_id + 1));
}

static inline uint64_t rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static uint64_t xoshiro_next(rng_state *r) {
    uint64_t *s = r->s;
    const uint64_t result = rotl64(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl64(s[3], 45);
    return result;
}

static double rng_uniform01(rng_state *r) {
    return (double)(xoshiro_next(r) >> 11) * (1.0 / 9007199254740992.0); /* 2^53 */
}

static double rng_uniform(rng_state *r, double lo, double hi) { return lo + (hi - lo) * rng_uniform01(r); }

static double rng_normal(rng_state *r) {
    /* Box-Muller, one draw per call (discards the paired value) -- simplicity over
     * throughput; at most p=5 draws are needed per SA iteration. */
    double u1 = rng_uniform01(r);
    if (u1 < 1e-300) u1 = 1e-300;
    double u2 = rng_uniform01(r);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ==== small linear algebra: principal angles between two small orthonormal bases ===
 * Both bases are given as d x D row-major arrays (d orthonormal row vectors of length D).
 * Principal angles' cosines are the singular values of A*B^T (d_a x d_b); we only need the
 * singular values, so it is enough to find the eigenvalues of the smaller of (A*B^T)*(A*B^T)^T
 * or (A*B^T)^T*(A*B^T), via a classic cyclic Jacobi eigenvalue solver -- exact and simple for
 * the small sizes (d <= a handful) every benchmark manifold here has. */

static void jacobi_eigenvalues_symmetric(double *a, int n, double *eigenvalues) {
    /* `a` is n x n row-major, overwritten as scratch. Off-diagonal Jacobi rotations until
     * converged (or a generous iteration cap -- always converges quickly for tiny n). */
    for (int iter = 0; iter < 100; iter++) {
        double off = 0.0;
        int p = 0, q = 1;
        double max_off = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double v = fabs(a[i * n + j]);
                off += v * v;
                if (v > max_off) {
                    max_off = v;
                    p = i;
                    q = j;
                }
            }
        }
        if (off < 1e-24 || n < 2) break;
        double app = a[p * n + p], aqq = a[q * n + q], apq = a[p * n + q];
        double theta = 0.5 * atan2(2.0 * apq, aqq - app);
        double c = cos(theta), s = sin(theta);
        for (int k = 0; k < n; k++) {
            double akp = a[k * n + p], akq = a[k * n + q];
            a[k * n + p] = c * akp - s * akq;
            a[k * n + q] = s * akp + c * akq;
        }
        for (int k = 0; k < n; k++) {
            double apk = a[p * n + k], aqk = a[q * n + k];
            a[p * n + k] = c * apk - s * aqk;
            a[q * n + k] = s * apk + c * aqk;
        }
    }
    for (int i = 0; i < n; i++) eigenvalues[i] = a[i * n + i];
}

/* Returns the average normalized angular loss (1/d) * sum(theta_i / (pi/2)) over the
 * "common" intrinsic dimension d = min(d_a, d_b) -- the same shared-rank convention
 * `estimate_stc_parameters`/`accept_ensemble` already use elsewhere in this codebase for
 * comparing tangent bases of possibly-differing rank. */
static double angular_loss(const double *basis_a, int d_a, const double *basis_b, int d_b, int D) {
    int d = d_a < d_b ? d_a : d_b;
    if (d <= 0) return 1.0; /* no meaningful tangent direction on either side: worst case */

    /* M = basis_a * basis_b^T, d_a x d_b */
    double *m = xmalloc(sizeof(double) * (size_t)d_a * (size_t)d_b);
    for (int i = 0; i < d_a; i++) {
        for (int j = 0; j < d_b; j++) {
            double dot = 0.0;
            for (int k = 0; k < D; k++) dot += basis_a[i * D + k] * basis_b[j * D + k];
            m[i * d_b + j] = dot;
        }
    }

    /* Eigenvalues of the smaller Gram matrix (M M^T or M^T M) -> squared singular values. */
    int small_dim = d_a < d_b ? d_a : d_b;
    double *gram = xmalloc(sizeof(double) * (size_t)small_dim * (size_t)small_dim);
    if (d_a <= d_b) {
        for (int i = 0; i < d_a; i++)
            for (int j = 0; j < d_a; j++) {
                double s = 0.0;
                for (int k = 0; k < d_b; k++) s += m[i * d_b + k] * m[j * d_b + k];
                gram[i * d_a + j] = s;
            }
    } else {
        for (int i = 0; i < d_b; i++)
            for (int j = 0; j < d_b; j++) {
                double s = 0.0;
                for (int k = 0; k < d_a; k++) s += m[k * d_b + i] * m[k * d_b + j];
                gram[i * d_b + j] = s;
            }
    }
    double *eig = xmalloc(sizeof(double) * (size_t)small_dim);
    jacobi_eigenvalues_symmetric(gram, small_dim, eig);

    /* Sort descending (selection sort -- small_dim is tiny) so the top `d` singular values
     * (largest cosines, smallest angles) are used, matching "the common intrinsic dimension". */
    for (int i = 0; i < small_dim; i++) {
        int best = i;
        for (int j = i + 1; j < small_dim; j++)
            if (eig[j] > eig[best]) best = j;
        double tmp = eig[i];
        eig[i] = eig[best];
        eig[best] = tmp;
    }

    double sum = 0.0;
    for (int i = 0; i < d; i++) {
        double cos_theta = i < small_dim ? eig[i] : 0.0;
        if (cos_theta < 0.0) cos_theta = 0.0;
        cos_theta = sqrt(cos_theta);
        if (cos_theta > 1.0) cos_theta = 1.0;
        double theta = acos(cos_theta);
        sum += theta / (M_PI / 2.0);
    }

    free(m);
    free(gram);
    free(eig);
    return sum / (double)d;
}

/* ==== benchmark dataset (JSON schema, see mod_STC.md "Benchmark dataset JSON schema") ==== */

typedef struct {
    int manifold_id;
    int intrinsic_dim;
    double noise_sd;
    int n_ref;
    double *ref_points;      /* n_ref x ambient_dim, row-major */
    double *ref_bases;       /* n_ref x intrinsic_dim x ambient_dim, row-major per ref point */
} manifold_t;

typedef struct {
    char *dataset_id;
    char *path;
    int ambient_dim;
    int n_manifolds;
    manifold_t *manifolds;
    int n_vectors;
    double *vectors;         /* ambient_dim x n_vectors, Fortran column-major (== point-major flat) */
    int *true_manifold_id;   /* n_vectors, values 0..n_manifolds-1 or -1 (sentinel: no true manifold) */
    int *kd_indices;         /* precomputed once, independent of phi */
    int *dimension_order;
    char *dim_names_buf;     /* "dim1".."dimD", padded per tox_conversions' string layout -- for
                              * --store-iteration-shatter-results, computed once at load time */
    int dim_names_strlen;
} dataset_t;

static double *parse_double_matrix(cJSON *arr2d, int *n_rows_out, int n_cols_expected) {
    int n_rows = cJSON_GetArraySize(arr2d);
    double *out = xmalloc(sizeof(double) * (size_t)n_rows * (size_t)n_cols_expected);
    int r = 0;
    cJSON *row;
    cJSON_ArrayForEach(row, arr2d) {
        int c = 0;
        cJSON *val;
        cJSON_ArrayForEach(val, row) {
            if (c >= n_cols_expected) die("dataset row has more than %d columns", n_cols_expected);
            out[r * n_cols_expected + c] = val->valuedouble;
            c++;
        }
        if (c != n_cols_expected) die("dataset row has %d columns, expected %d", c, n_cols_expected);
        r++;
    }
    *n_rows_out = n_rows;
    return out;
}

static dataset_t load_dataset(const char *path) {
    FILE *fh = fopen(path, "rb");
    if (fh == NULL) die("cannot open '%s'", path);
    fseek(fh, 0, SEEK_END);
    long len = ftell(fh);
    fseek(fh, 0, SEEK_SET);
    char *buf = xmalloc((size_t)len + 1);
    if (fread(buf, 1, (size_t)len, fh) != (size_t)len) die("short read on '%s'", path);
    buf[len] = '\0';
    fclose(fh);

    cJSON *root = cJSON_Parse(buf);
    free(buf);
    if (root == NULL) die("invalid JSON in '%s'", path);

    dataset_t d;
    memset(&d, 0, sizeof(d));
    d.path = strdup(path);

    cJSON *dataset_id = cJSON_GetObjectItemCaseSensitive(root, "dataset_id");
    d.dataset_id = strdup(dataset_id != NULL && cJSON_IsString(dataset_id) ? dataset_id->valuestring : path);

    cJSON *ambient_dim = cJSON_GetObjectItemCaseSensitive(root, "ambient_dim");
    if (ambient_dim == NULL) die("'%s': missing ambient_dim", path);
    d.ambient_dim = ambient_dim->valueint;

    cJSON *manifolds = cJSON_GetObjectItemCaseSensitive(root, "manifolds");
    if (manifolds == NULL || !cJSON_IsArray(manifolds)) die("'%s': missing manifolds[]", path);
    d.n_manifolds = cJSON_GetArraySize(manifolds);
    d.manifolds = xcalloc((size_t)d.n_manifolds, sizeof(manifold_t));

    int mi = 0;
    cJSON *m;
    cJSON_ArrayForEach(m, manifolds) {
        manifold_t *out = &d.manifolds[mi];
        out->manifold_id = cJSON_GetObjectItemCaseSensitive(m, "manifold_id")->valueint;
        out->intrinsic_dim = cJSON_GetObjectItemCaseSensitive(m, "intrinsic_dim")->valueint;
        out->noise_sd = cJSON_GetObjectItemCaseSensitive(m, "noise_sd")->valuedouble;
        if (out->noise_sd <= 0.0) die("'%s': manifold %d has noise_sd <= 0", path, out->manifold_id);

        cJSON *ref_points = cJSON_GetObjectItemCaseSensitive(m, "reference_points");
        out->ref_points = parse_double_matrix(ref_points, &out->n_ref, d.ambient_dim);

        cJSON *ref_bases = cJSON_GetObjectItemCaseSensitive(m, "reference_tangent_bases");
        int n_bases = cJSON_GetArraySize(ref_bases);
        if (n_bases != out->n_ref) die("'%s': manifold %d has mismatched reference_points/reference_tangent_bases counts", path, out->manifold_id);
        out->ref_bases = xmalloc(sizeof(double) * (size_t)out->n_ref * (size_t)out->intrinsic_dim * (size_t)d.ambient_dim);
        int ri = 0;
        cJSON *basis;
        cJSON_ArrayForEach(basis, ref_bases) {
            int n_rows;
            double *flat = parse_double_matrix(basis, &n_rows, d.ambient_dim);
            if (n_rows != out->intrinsic_dim) die("'%s': manifold %d ref point %d basis has %d rows, expected intrinsic_dim=%d", path, out->manifold_id, ri, n_rows, out->intrinsic_dim);
            memcpy(out->ref_bases + (size_t)ri * out->intrinsic_dim * d.ambient_dim, flat,
                   sizeof(double) * (size_t)out->intrinsic_dim * (size_t)d.ambient_dim);
            free(flat);
            ri++;
        }
        mi++;
    }

    cJSON *points = cJSON_GetObjectItemCaseSensitive(root, "points");
    if (points == NULL) die("'%s': missing points", path);
    cJSON *coordinates = cJSON_GetObjectItemCaseSensitive(points, "coordinates");
    /* Row-per-point (N x ambient_dim) JSON flattens identically to Fortran's column-major
     * vectors(ambient_dim, N) -- both conventions store one point's ambient_dim coordinates
     * contiguously, so no transpose is needed (see mod_STC.md's CLI section). */
    d.vectors = parse_double_matrix(coordinates, &d.n_vectors, d.ambient_dim);

    cJSON *true_ids = cJSON_GetObjectItemCaseSensitive(points, "true_manifold_id");
    if (cJSON_GetArraySize(true_ids) != d.n_vectors) die("'%s': true_manifold_id length mismatch", path);
    d.true_manifold_id = xmalloc(sizeof(int) * (size_t)d.n_vectors);
    int pi = 0;
    cJSON *idv;
    cJSON_ArrayForEach(idv, true_ids) {
        int v = idv->valueint;
        if (v != -1 && (v < 0 || v >= d.n_manifolds)) die("'%s': true_manifold_id[%d]=%d resolves to no declared manifold and is not the -1 sentinel", path, pi, v);
        d.true_manifold_id[pi] = v;
        pi++;
    }

    cJSON_Delete(root);

    /* k-d tree: independent of phi, computed once here and reused across every chain/iteration. */
    d.dimension_order = xmalloc(sizeof(int) * (size_t)d.ambient_dim);
    for (int i = 0; i < d.ambient_dim; i++) d.dimension_order[i] = i + 1;
    d.kd_indices = xmalloc(sizeof(int) * (size_t)d.n_vectors);
    int ierr = 0;
    build_kd_index_c(d.vectors, &d.ambient_dim, &d.n_vectors, d.kd_indices, d.dimension_order, &ierr);
    check_ierr(ierr, "build_kd_index");

    /* "dim1".."dimD", the same generated-name convention stc_cli.c uses for a header-less
     * CSV -- these JSON benchmark datasets never carry dimension names of their own. */
    d.dim_names_strlen = 1;
    for (int i = 0; i < d.ambient_dim; i++) {
        char name[32];
        snprintf(name, sizeof(name), "dim%d", i + 1);
        if ((int)strlen(name) > d.dim_names_strlen) d.dim_names_strlen = (int)strlen(name);
    }
    d.dim_names_buf = xmalloc((size_t)d.dim_names_strlen * (size_t)d.ambient_dim);
    for (int i = 0; i < d.ambient_dim; i++) {
        char name[32];
        snprintf(name, sizeof(name), "dim%d", i + 1);
        fill_padded_string(d.dim_names_buf + (size_t)i * (size_t)d.dim_names_strlen, (size_t)d.dim_names_strlen, name);
    }

    return d;
}

static char *dataset_stem(const char *path) {
    const char *base = strrchr(path, '/');
    base = base != NULL ? base + 1 : path;
    char *out = strdup(base);
    char *dot = strrchr(out, '.');
    if (dot != NULL) *dot = '\0';
    return out;
}

/* ==== fixed (non-phi) STC parameters, one shared value for the whole SA run ========= */

typedef struct {
    int o;
    double f_max;
    double bandwidth_percentile;
    double exclusion_radius_percentile;
    const char *reconciliation_mode;
    double min_overlap_coefficient;
    int n_anchors;
    double seed_max_set_size;
    double quantile_pairwise_ea_comparison;
} fixed_params_t;

/* ==== phi: the 7-parameter STC vector SA tunes ====================================== */

#define PHI_P 7
enum {
    PHI_K_MIN = 0,
    PHI_CHORDAL,
    PHI_D_MAX,
    PHI_G_MAX,
    PHI_RMSE_CHANGE_MAX,
    PHI_RADIUS_PERCENTILE,
    PHI_MIN_STABLE_ITERATIONS,
};

typedef struct {
    double lo[PHI_P];
    double hi[PHI_P];
} phi_bounds_t;

static void phi_reflect(double *phi, const phi_bounds_t *b) {
    for (int i = 0; i < PHI_P; i++) {
        double lo = b->lo[i], hi = b->hi[i];
        double v = phi[i];
        /* Reflect at the boundary (preferable to clipping, which accumulates mass at the
         * edge -- see mod_STC.md's "Parameter scales and proposal step sizes"). A single
         * reflection is enough in practice since proposal steps stay well inside the range. */
        if (v < lo) v = lo + (lo - v);
        if (v > hi) v = hi - (v - hi);
        if (v < lo) v = lo;
        if (v > hi) v = hi;
        phi[i] = v;
    }
}

/* ==== one STC run's raw output, and per-dataset objective computation ============== */

typedef struct {
    int n_selected_seed;
    unsigned char *is_seed_mask;
    unsigned char *ensemble_masks;       /* n_vectors x n_selected_seed */
    int *ensemble_stop_reason;
    double *ensemble_growth_radii;
    double *ensemble_U_history;
    double *ensemble_S_history;
    int *ensemble_d_history;
    double *ensemble_G_history;
    double *ensemble_mu_history;
    int *ensemble_k_history;
    unsigned char *ensemble_accepted_history;
    int *ensemble_member_added_at_step;
    unsigned char *ensemble_low_confidence_masks;
    double *ensemble_U_first;            /* ambient_dim x ambient_dim x n_selected_seed, columns 0..d_first-1 valid */
    int *ensemble_d_first;
    double *ensemble_U_final;            /* ambient_dim x ambient_dim x n_selected_seed, columns 0..d_final-1 valid -- what L_A uses */
    int *ensemble_d_final;
    int max_group_size;
    int n_super_ensembles;
    int *super_ensembles;                /* max_group_size x super_ensembles_capacity, 1-indexed ensemble ids, 0 = empty */
    unsigned char *eligible, *eligible_by_stop, *eligible_by_dim, *eligible_by_var;
} stc_run_t;

static void free_stc_run(stc_run_t *run) {
    free(run->is_seed_mask);
    free(run->ensemble_masks);
    free(run->ensemble_stop_reason);
    free(run->ensemble_growth_radii);
    free(run->ensemble_U_history);
    free(run->ensemble_S_history);
    free(run->ensemble_d_history);
    free(run->ensemble_G_history);
    free(run->ensemble_mu_history);
    free(run->ensemble_k_history);
    free(run->ensemble_accepted_history);
    free(run->ensemble_member_added_at_step);
    free(run->ensemble_low_confidence_masks);
    free(run->ensemble_U_first);
    free(run->ensemble_d_first);
    free(run->ensemble_U_final);
    free(run->ensemble_d_final);
    free(run->super_ensembles);
    free(run->eligible);
    free(run->eligible_by_stop);
    free(run->eligible_by_dim);
    free(run->eligible_by_var);
}

/* Runs seeds -> ensemble_identification_merged -> ensemble_reconciliation for one dataset at
 * one phi, exactly the sequence stc_cli.c uses. Every raw output buffer is kept (not just the
 * ones compute_losses reads) so --store-iteration-shatter-results can serialize the full
 * result without a second, duplicate run -- "plenty of RAM" per mod_STC.md's own sizing
 * discussion, and simpler than two different buffer-lifetime paths. `eligible`/`super_ensembles`
 * are allocated unconditionally (matching stc_cli.c's own defensive sizing) even when
 * reconciliation itself does not run (n_selected_seed < 2), so downstream consumers never have
 * to treat "reconciliation didn't run" as "these pointers might be NULL". */
static stc_run_t run_stc(const dataset_t *ds, const double *phi, const fixed_params_t *fx) {
    stc_run_t run;
    memset(&run, 0, sizeof(run));

    int k_min = (int)lround(phi[PHI_K_MIN]);
    double chordal = phi[PHI_CHORDAL];
    int d_max = (int)lround(phi[PHI_D_MAX]);
    double G_max = phi[PHI_G_MAX];
    double RMSE_change_max = phi[PHI_RMSE_CHANGE_MAX];
    double radius_percentile = phi[PHI_RADIUS_PERCENTILE];
    int min_stable_iterations = (int)lround(phi[PHI_MIN_STABLE_ITERATIONS]);

    run.is_seed_mask = xmalloc((size_t)ds->n_vectors);
    int ierr = 0;
    seeds_c(ds->vectors, &ds->ambient_dim, &ds->n_vectors, ds->kd_indices, ds->dimension_order, &k_min,
           &fx->bandwidth_percentile, &fx->exclusion_radius_percentile, run.is_seed_mask, &ierr);
    check_ierr(ierr, "seeds");

    int n_selected_seed = 0;
    for (int i = 0; i < ds->n_vectors; i++)
        if (run.is_seed_mask[i]) n_selected_seed++;
    run.n_selected_seed = n_selected_seed;

    if (n_selected_seed == 0) return run; /* nothing grew: caller treats this as worst-case losses */

    /* Pure Fortran-array-sizing safety cap, not a clustering-quality knob -- no CLI flag, see
     * mod_STC.md's own reasoning for dropping it from stc_cli's surface too. */
    int max_group_size = n_selected_seed < 1024 ? n_selected_seed : 1024;
    if (max_group_size < 2) max_group_size = 2;
    run.max_group_size = max_group_size;

    run.ensemble_masks = xmalloc((size_t)ds->n_vectors * (size_t)n_selected_seed);
    run.ensemble_stop_reason = xmalloc(sizeof(int) * (size_t)n_selected_seed);
    run.ensemble_growth_radii = xmalloc(sizeof(double) * (size_t)n_selected_seed);
    run.ensemble_U_history = xmalloc(sizeof(double) * (size_t)ds->ambient_dim * (size_t)ds->ambient_dim * (size_t)fx->o * (size_t)n_selected_seed);
    run.ensemble_S_history = xmalloc(sizeof(double) * (size_t)ds->ambient_dim * (size_t)fx->o * (size_t)n_selected_seed);
    run.ensemble_d_history = xmalloc(sizeof(int) * (size_t)fx->o * (size_t)n_selected_seed);
    run.ensemble_G_history = xmalloc(sizeof(double) * (size_t)fx->o * (size_t)n_selected_seed);
    run.ensemble_mu_history = xmalloc(sizeof(double) * (size_t)ds->ambient_dim * (size_t)fx->o * (size_t)n_selected_seed);
    run.ensemble_k_history = xmalloc(sizeof(int) * (size_t)fx->o * (size_t)n_selected_seed);
    run.ensemble_accepted_history = xmalloc((size_t)fx->o * (size_t)n_selected_seed);
    run.ensemble_member_added_at_step = xmalloc(sizeof(int) * (size_t)ds->n_vectors * (size_t)n_selected_seed);
    run.ensemble_low_confidence_masks = xmalloc((size_t)ds->n_vectors * (size_t)n_selected_seed);
    run.ensemble_U_first = xmalloc(sizeof(double) * (size_t)ds->ambient_dim * (size_t)ds->ambient_dim * (size_t)n_selected_seed);
    run.ensemble_d_first = xmalloc(sizeof(int) * (size_t)n_selected_seed);

    ensemble_identification_merged_c(
        ds->vectors, &ds->ambient_dim, &ds->n_vectors, ds->kd_indices, ds->dimension_order, run.is_seed_mask,
        &n_selected_seed, &k_min, &chordal, &d_max, &G_max, &RMSE_change_max, &fx->f_max,
        &min_stable_iterations, &radius_percentile, &fx->o,
        run.ensemble_masks, run.ensemble_stop_reason, run.ensemble_growth_radii, run.ensemble_U_history,
        run.ensemble_S_history, run.ensemble_d_history, run.ensemble_G_history, run.ensemble_mu_history,
        run.ensemble_k_history, run.ensemble_accepted_history, run.ensemble_member_added_at_step,
        run.ensemble_low_confidence_masks, run.ensemble_U_first, run.ensemble_d_first, &ierr);
    check_ierr(ierr, "ensemble_identification_merged");

    /* The ensemble's real final accepted tangent basis -- what L_A must compare against (see
     * the extern declaration's comment above). Only U_final/d_final are kept; the rest of this
     * kernel's outputs (S/mu/G/k_final, has_final, final_index) are scratch mod_STC.md's L_D
     * doesn't need, since L_D's center is the literal mean of final membership it documents,
     * not the growth algorithm's own internal running center estimate. */
    run.ensemble_U_final = xmalloc(sizeof(double) * (size_t)ds->ambient_dim * (size_t)ds->ambient_dim * (size_t)n_selected_seed);
    run.ensemble_d_final = xmalloc(sizeof(int) * (size_t)n_selected_seed);
    {
        double *scratch_S_final = xmalloc(sizeof(double) * (size_t)ds->ambient_dim * (size_t)n_selected_seed);
        double *scratch_mu_final = xmalloc(sizeof(double) * (size_t)ds->ambient_dim * (size_t)n_selected_seed);
        double *scratch_G_final = xmalloc(sizeof(double) * (size_t)n_selected_seed);
        int *scratch_k_final = xmalloc(sizeof(int) * (size_t)n_selected_seed);
        unsigned char *scratch_has_final = xmalloc((size_t)n_selected_seed);
        int *scratch_final_index = xmalloc(sizeof(int) * (size_t)n_selected_seed);

        ensemble_final_observable_c(&ds->ambient_dim, &fx->o, &n_selected_seed, run.ensemble_U_history,
                                    run.ensemble_d_history, run.ensemble_S_history, run.ensemble_mu_history,
                                    run.ensemble_G_history, run.ensemble_k_history, run.ensemble_accepted_history,
                                    run.ensemble_U_final, run.ensemble_d_final, scratch_S_final, scratch_mu_final,
                                    scratch_G_final, scratch_k_final, scratch_has_final, scratch_final_index, &ierr);
        check_ierr(ierr, "ensemble_final_observable");

        free(scratch_S_final);
        free(scratch_mu_final);
        free(scratch_G_final);
        free(scratch_k_final);
        free(scratch_has_final);
        free(scratch_final_index);
    }

    int n_super_ensembles_capacity = n_selected_seed * (n_selected_seed - 1);
    if (n_super_ensembles_capacity < 1) n_super_ensembles_capacity = 1;
    run.super_ensembles = xcalloc((size_t)max_group_size * (size_t)n_super_ensembles_capacity, sizeof(int));
    run.eligible = xmalloc((size_t)n_selected_seed);
    run.eligible_by_stop = xmalloc((size_t)n_selected_seed);
    run.eligible_by_dim = xmalloc((size_t)n_selected_seed);
    run.eligible_by_var = xmalloc((size_t)n_selected_seed);
    memset(run.eligible, 1, (size_t)n_selected_seed);
    memset(run.eligible_by_stop, 1, (size_t)n_selected_seed);
    memset(run.eligible_by_dim, 1, (size_t)n_selected_seed);
    memset(run.eligible_by_var, 1, (size_t)n_selected_seed);

    if (n_selected_seed >= 2) {
        double *super_ensembles_overlap_coefficient =
            xmalloc(sizeof(double) * (size_t)(max_group_size - 1) * (size_t)n_super_ensembles_capacity);

        char mode_buf[25];
        fill_padded_string(mode_buf, sizeof(mode_buf), fx->reconciliation_mode);
        unsigned char report_oc = 0;

        ensemble_reconciliation_c(run.ensemble_masks, run.ensemble_stop_reason, &ds->ambient_dim, &ds->n_vectors,
                                  &n_selected_seed, run.ensemble_U_history, run.ensemble_d_history,
                                  run.ensemble_S_history, run.ensemble_mu_history, run.ensemble_G_history,
                                  run.ensemble_k_history, run.ensemble_accepted_history, &fx->o, mode_buf,
                                  &fx->min_overlap_coefficient, &report_oc, NULL, NULL, NULL, NULL,
                                  &max_group_size, run.super_ensembles, &run.n_super_ensembles,
                                  super_ensembles_overlap_coefficient, run.eligible, run.eligible_by_stop,
                                  run.eligible_by_dim, run.eligible_by_var, &ierr);
        check_ierr(ierr, "ensemble_reconciliation");

        free(super_ensembles_overlap_coefficient);
    }
    return run;
}

/* Writes the same JSON+CSV/TSV artifacts stc_cli.c writes for a standalone run (never the HTML
 * report) into `output_subdir`, which the caller has already created. No per-iteration
 * "estimated_*" values are attached (estimate_stc_parameters only ever runs once, before the SA
 * loop, to seed chain 0's initial phi -- not once per iteration), so those six pointers are
 * always NULL here. Silently does nothing when nothing grew (n_selected_seed == 0): there is
 * nothing meaningful to report, and several buffers below stay unallocated in that case. */
static void write_shatter_iteration_output(const dataset_t *ds, const stc_run_t *run, const double *phi,
                                           const fixed_params_t *fx, const char *output_subdir) {
    if (run->n_selected_seed == 0) return;

    int k_min = (int)lround(phi[PHI_K_MIN]);
    double chordal = phi[PHI_CHORDAL];
    int d_max = (int)lround(phi[PHI_D_MAX]);
    double G_max = phi[PHI_G_MAX];
    double RMSE_change_max = phi[PHI_RMSE_CHANGE_MAX];
    double radius_percentile = phi[PHI_RADIUS_PERCENTILE];
    int min_stable_iterations = (int)lround(phi[PHI_MIN_STABLE_ITERATIONS]);
    int ierr = 0;

    char mode_buf[25];
    fill_padded_string(mode_buf, sizeof(mode_buf), fx->reconciliation_mode);

    char *json_path = join_path(output_subdir, "results.json");
    int json_path_len = (int)strlen(json_path);
    serialize_stc_results_as_json_c(
        json_path, &json_path_len, &ds->ambient_dim, &ds->n_vectors, &run->n_selected_seed, &fx->o,
        &run->max_group_size, &run->n_super_ensembles, ds->vectors, ds->dim_names_buf, &ds->dim_names_strlen,
        run->is_seed_mask, run->ensemble_masks, run->ensemble_stop_reason, run->ensemble_growth_radii,
        run->ensemble_U_history, run->ensemble_S_history, run->ensemble_d_history, run->ensemble_G_history,
        run->ensemble_mu_history, run->ensemble_k_history, run->ensemble_accepted_history,
        run->ensemble_member_added_at_step, run->ensemble_low_confidence_masks, run->ensemble_U_first,
        run->ensemble_d_first, run->super_ensembles, &k_min, &chordal, &d_max, &G_max,
        &RMSE_change_max, &fx->f_max, &min_stable_iterations, &radius_percentile,
        &fx->exclusion_radius_percentile, &fx->bandwidth_percentile,
        mode_buf, &fx->min_overlap_coefficient, NULL, NULL, NULL, NULL, run->eligible, run->eligible_by_stop,
        run->eligible_by_dim, run->eligible_by_var, NULL, NULL, NULL, NULL, NULL, &ierr);
    check_ierr(ierr, "serialize_stc_results_as_json (shatter)");
    free(json_path);

    char *points_path = join_path(output_subdir, "points.csv");
    int points_path_len = (int)strlen(points_path);
    serialize_stc_points_as_csv_c(points_path, &points_path_len, &ds->n_vectors, &run->n_selected_seed,
                                  &run->max_group_size, &run->n_super_ensembles, run->is_seed_mask,
                                  run->ensemble_masks, run->ensemble_low_confidence_masks, run->super_ensembles,
                                  &ierr);
    check_ierr(ierr, "serialize_stc_points_as_csv (shatter)");
    free(points_path);

    char *overlap_path = join_path(output_subdir, "ensemble_overlap_coefficients.csv");
    int overlap_path_len = (int)strlen(overlap_path);
    serialize_stc_ensemble_overlap_as_csv_c(overlap_path, &overlap_path_len, &ds->n_vectors,
                                            &run->n_selected_seed, run->ensemble_masks, &ierr);
    check_ierr(ierr, "serialize_stc_ensemble_overlap_as_csv (shatter)");
    free(overlap_path);

    char *super_path = join_path(output_subdir, "super_ensembles.tsv");
    int super_path_len = (int)strlen(super_path);
    serialize_stc_super_ensembles_as_tsv_c(super_path, &super_path_len, &run->max_group_size,
                                           &run->n_super_ensembles, run->super_ensembles, &ierr);
    check_ierr(ierr, "serialize_stc_super_ensembles_as_tsv (shatter)");
    free(super_path);
}

/* mkdir -p, minimal: only ever needs to create at most two missing levels
 * (--output-dir/{dataset_stem}_chain{N}/iter_{NNNNNN}), so no need for a general recursive walk. */
static void mkdir_or_die(const char *path) {
    if (mkdir(path, 0755) != 0 && errno != EEXIST) die("cannot create directory '%s': %s", path, strerror(errno));
}

/* ==== objective function: mod_STC.md's L_A, L_D, L_C, L_E, and reporting-only F ===== */

typedef struct {
    double L_A, L_D, L_C, L_E, F;
} losses_t;

/* Brute-force nearest reference point (linear scan -- reference clouds are at most tens of
 * thousands of points and this runs once per accepted ensemble per iteration, not per data
 * point, so a k-d tree here would be premature). Returns squared distance via *out_dist2. */
static int nearest_ref_point(const manifold_t *m, const double *query, int D, double *out_dist2) {
    int best = -1;
    double best_d2 = 1e300;
    for (int r = 0; r < m->n_ref; r++) {
        double d2 = 0.0;
        for (int k = 0; k < D; k++) {
            double diff = query[k] - m->ref_points[r * D + k];
            d2 += diff * diff;
        }
        if (d2 < best_d2) {
            best_d2 = d2;
            best = r;
        }
    }
    *out_dist2 = best_d2;
    return best;
}

/* Exact optimal one-to-one matching (manifolds -> super-ensemble groups) maximizing total
 * matched point count, via branch-and-bound permutation search -- exact for the small ensemble
 * counts these benchmarks have; falls back to a greedy approximation (documented, reporting-only
 * metric, not used by the SA objective) above a defensive size cap to avoid combinatorial blowup. */
static void match_search(const long *contingency, int rows, int cols, int *col_used, int depth,
                         long current_sum, long remaining_max_sum, long *best) {
    if (current_sum + remaining_max_sum <= *best) return; /* admissible upper-bound prune */
    if (depth == rows) {
        if (current_sum > *best) *best = current_sum;
        return;
    }
    long row_max = 0;
    for (int c = 0; c < cols; c++)
        if (!col_used[c] && contingency[depth * cols + c] > row_max) row_max = contingency[depth * cols + c];

    /* leave this manifold unmatched */
    match_search(contingency, rows, cols, col_used, depth + 1, current_sum, remaining_max_sum - row_max, best);
    for (int c = 0; c < cols; c++) {
        if (col_used[c]) continue;
        col_used[c] = 1;
        match_search(contingency, rows, cols, col_used, depth + 1, current_sum + contingency[depth * cols + c],
                    remaining_max_sum - row_max, best);
        col_used[c] = 0;
    }
}

static long best_matching(const long *contingency, int rows, int cols) {
    if (rows == 0 || cols == 0) return 0;
    long remaining_max = 0;
    for (int r = 0; r < rows; r++) {
        long row_max = 0;
        for (int c = 0; c < cols; c++)
            if (contingency[r * cols + c] > row_max) row_max = contingency[r * cols + c];
        remaining_max += row_max;
    }
    if (rows > 12 || cols > 12) {
        /* Greedy fallback: repeatedly take the largest remaining cell. Documented
         * approximation for pathologically large ensemble counts -- not expected in practice
         * for this benchmark suite (never more than a handful of manifolds). */
        long *work = xmalloc(sizeof(long) * (size_t)rows * (size_t)cols);
        memcpy(work, contingency, sizeof(long) * (size_t)rows * (size_t)cols);
        int *row_done = xcalloc((size_t)rows, sizeof(int));
        int *col_done = xcalloc((size_t)cols, sizeof(int));
        long total = 0;
        for (int step = 0; step < rows; step++) {
            long best_v = -1;
            int best_r = -1, best_c = -1;
            for (int r = 0; r < rows; r++) {
                if (row_done[r]) continue;
                for (int c = 0; c < cols; c++) {
                    if (col_done[c]) continue;
                    if (work[r * cols + c] > best_v) {
                        best_v = work[r * cols + c];
                        best_r = r;
                        best_c = c;
                    }
                }
            }
            if (best_r < 0) break;
            row_done[best_r] = 1;
            col_done[best_c] = 1;
            total += best_v;
        }
        free(work);
        free(row_done);
        free(col_done);
        return total;
    }
    int *col_used = xcalloc((size_t)cols, sizeof(int));
    long best = 0;
    match_search(contingency, rows, cols, col_used, 0, 0, remaining_max, &best);
    free(col_used);
    return best;
}

static double compute_ari(const long *contingency, const long *row_sum, const long *col_sum, int rows,
                          int cols, long n) {
    /* Degenerate cases where no disagreement between the two partitions is even possible --
     * conventionally (matching scikit-learn's adjusted_rand_score) defined as ARI=1 (perfect),
     * not 0: fewer than 2 covered points, only one true class and one predicted class present
     * (e.g. a single-manifold dataset correctly recovered as a single super-ensemble), or every
     * covered point is alone in both its true class and its predicted class. */
    if (n < 2) return 1.0;
    int rows_present = 0, cols_present = 0;
    for (int r = 0; r < rows; r++)
        if (row_sum[r] > 0) rows_present++;
    for (int c = 0; c < cols; c++)
        if (col_sum[c] > 0) cols_present++;
    if (rows_present <= 1 && cols_present <= 1) return 1.0;
    if (rows_present == (int)n && cols_present == (int)n) return 1.0;

    double sum_comb_c = 0.0;
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++) {
            long v = contingency[r * cols + c];
            sum_comb_c += (double)v * (double)(v - 1) / 2.0;
        }
    double sum_comb_rows = 0.0;
    for (int r = 0; r < rows; r++) sum_comb_rows += (double)row_sum[r] * (double)(row_sum[r] - 1) / 2.0;
    double sum_comb_cols = 0.0;
    for (int c = 0; c < cols; c++) sum_comb_cols += (double)col_sum[c] * (double)(col_sum[c] - 1) / 2.0;
    double comb_n = (double)n * (double)(n - 1) / 2.0;
    double expected = comb_n > 0.0 ? (sum_comb_rows * sum_comb_cols) / comb_n : 0.0;
    double max_index = 0.5 * (sum_comb_rows + sum_comb_cols);
    double denom = max_index - expected;
    if (fabs(denom) < 1e-12) return 1.0; /* remaining degenerate cases: no room to disagree */
    return (sum_comb_c - expected) / denom;
}

static losses_t compute_losses(const dataset_t *ds, const stc_run_t *run) {
    /* L_A/L_D default to worst-case (1.0) when nothing is accepted: unlike L_E below, these
     * measure the quality of what STC actually produced, so "produced nothing" is a real,
     * distinct failure worth penalizing on top of L_C, not one to leave neutral.
     * L_E defaults to 0.0 (perfect/no ARI penalty) when there are no covered points at all --
     * mod_STC.md's L_E section explicitly avoids double-penalizing lack of coverage (L_C
     * already does that); see compute_ari's own degenerate-case handling for the matching
     * n<2-covered-points convention when some but very few points are covered. */
    losses_t out = {1.0, 1.0, 1.0, 0.0, -1.0};

    /* per-point predicted group id (0-indexed), -1 if not covered by any accepted super-ensemble */
    int *point_group = xmalloc(sizeof(int) * (size_t)ds->n_vectors);
    for (int i = 0; i < ds->n_vectors; i++) point_group[i] = -1;

    if (run->n_super_ensembles > 0) {
        for (int g = 0; g < run->n_super_ensembles; g++) {
            for (int r = 0; r < run->max_group_size; r++) {
                int ensemble_id_1based = run->super_ensembles[r + g * run->max_group_size];
                if (ensemble_id_1based == 0) continue;
                int e = ensemble_id_1based - 1;
                for (int i = 0; i < ds->n_vectors; i++) {
                    if (run->ensemble_masks[i + (size_t)e * ds->n_vectors] && point_group[i] < 0) point_group[i] = g;
                }
            }
        }
    }

    /* ---- L_A, L_D: point-weighted average over ensembles that are members of an accepted
     * super-ensemble (an ensemble "belongs" to group g if any of its member points map to g
     * above -- since group assignment is per-point via ensemble_masks, re-derive per-ensemble
     * membership directly from the same super_ensembles listing). ---- */
    double sum_w_A = 0.0, sum_wL_A = 0.0;
    double sum_w_D = 0.0, sum_wL_D = 0.0;
    if (run->n_super_ensembles > 0) {
        for (int g = 0; g < run->n_super_ensembles; g++) {
            for (int r = 0; r < run->max_group_size; r++) {
                int ensemble_id_1based = run->super_ensembles[r + g * run->max_group_size];
                if (ensemble_id_1based == 0) continue;
                int e = ensemble_id_1based - 1;

                int size = 0;
                double *center = xcalloc((size_t)ds->ambient_dim, sizeof(double));
                for (int i = 0; i < ds->n_vectors; i++) {
                    if (run->ensemble_masks[i + (size_t)e * ds->n_vectors]) {
                        size++;
                        for (int k = 0; k < ds->ambient_dim; k++)
                            center[k] += ds->vectors[k + (size_t)i * ds->ambient_dim];
                    }
                }
                if (size == 0) {
                    free(center);
                    continue;
                }
                for (int k = 0; k < ds->ambient_dim; k++) center[k] /= (double)size;

                /* nearest true manifold to this ensemble's center */
                int best_manifold = -1;
                double best_d2 = 1e300;
                int best_ref = -1;
                for (int mi = 0; mi < ds->n_manifolds; mi++) {
                    double d2;
                    int ref = nearest_ref_point(&ds->manifolds[mi], center, ds->ambient_dim, &d2);
                    if (d2 < best_d2) {
                        best_d2 = d2;
                        best_manifold = mi;
                        best_ref = ref;
                    }
                }
                free(center);
                if (best_manifold < 0) continue; /* no manifolds at all: skip (shouldn't happen) */

                const manifold_t *mref = &ds->manifolds[best_manifold];
                double d_k = sqrt(best_d2);
                double L_D_k = d_k / mref->noise_sd;
                if (L_D_k > 1.0) L_D_k = 1.0;
                sum_w_D += (double)size;
                sum_wL_D += (double)size * L_D_k;

                int d_ens = run->ensemble_d_final[e];
                if (d_ens > 0) {
                    const double *ens_basis_cols = run->ensemble_U_final + (size_t)e * ds->ambient_dim * ds->ambient_dim;
                    /* ensemble_U_final is Fortran column-major (ambient_dim, ambient_dim, seed):
                     * basis vector j (0-indexed) is the contiguous block of ambient_dim doubles
                     * starting at column j -- already exactly the "row-major d x D" shape
                     * angular_loss expects if we treat consecutive columns as consecutive rows.
                     * Deliberately the ensemble's *final* accepted tangent basis, matching the
                     * final membership L_D's center above is computed from -- ensemble_U_first
                     * would be the early, still-small first-accepted snapshot, a real bug fixed
                     * after a huge k_min (bloating an ensemble structure-agnostically) turned
                     * out to score a deceptively good L_A from that stale early estimate. */
                    const double *ref_basis = mref->ref_bases + (size_t)best_ref * mref->intrinsic_dim * ds->ambient_dim;
                    double L_A_k = angular_loss(ens_basis_cols, d_ens, ref_basis, mref->intrinsic_dim, ds->ambient_dim);
                    sum_w_A += (double)size;
                    sum_wL_A += (double)size * L_A_k;
                }
            }
        }
    }
    if (sum_w_A > 0.0) out.L_A = sum_wL_A / sum_w_A;
    if (sum_w_D > 0.0) out.L_D = sum_wL_D / sum_w_D;

    /* ---- L_C: coverage over points with a real true manifold (sentinel -1 excluded from
     * both N and N_covered, per mod_STC.md's amended "Data-point coverage"). ---- */
    long N = 0, N_covered = 0;
    for (int i = 0; i < ds->n_vectors; i++) {
        if (ds->true_manifold_id[i] == -1) continue;
        N++;
        if (point_group[i] >= 0) N_covered++;
    }
    if (N > 0) out.L_C = 1.0 - (double)N_covered / (double)N;

    /* ---- L_E (ARI, over ALL covered points, sentinel included as its own class) and F
     * (reporting-only matched fraction, over covered points with a real true manifold only). ---- */
    if (run->n_super_ensembles > 0) {
        int true_classes = ds->n_manifolds + 1; /* index n_manifolds = sentinel -1 */
        int pred_classes = run->n_super_ensembles;
        long *contingency = xcalloc((size_t)true_classes * (size_t)pred_classes, sizeof(long));
        long *row_sum = xcalloc((size_t)true_classes, sizeof(long));
        long *col_sum = xcalloc((size_t)pred_classes, sizeof(long));
        long n_covered_any = 0;
        for (int i = 0; i < ds->n_vectors; i++) {
            if (point_group[i] < 0) continue;
            int t = ds->true_manifold_id[i] == -1 ? ds->n_manifolds : ds->true_manifold_id[i];
            int p = point_group[i];
            contingency[t * pred_classes + p]++;
            row_sum[t]++;
            col_sum[p]++;
            n_covered_any++;
        }
        double ari = compute_ari(contingency, row_sum, col_sum, true_classes, pred_classes, n_covered_any);
        out.L_E = 1.0 - (ari > 0.0 ? ari : 0.0);

        if (N_covered > 0) {
            /* matching restricted to real manifolds (rows 0..n_manifolds-1, excluding the
             * sentinel row) against super-ensemble groups */
            long matched = best_matching(contingency, ds->n_manifolds, pred_classes);
            out.F = (double)matched / (double)N_covered;
        }
        free(contingency);
        free(row_sum);
        free(col_sum);
    }

    free(point_group);
    return out;
}

static double combined_loss(const losses_t *l) { return (l->L_A + l->L_D + l->L_C + l->L_E) / 4.0; }

/* ==== CSV logging: one file per (dataset, chain), see mod_STC.md "Logging" ========== */

typedef struct {
    FILE **files; /* one per dataset, indexed the same as the datasets array */
} chain_logs_t;

static chain_logs_t open_chain_logs(const char *output_dir, const dataset_t *datasets, int n_datasets,
                                    int chain_id) {
    chain_logs_t logs;
    logs.files = xmalloc(sizeof(FILE *) * (size_t)n_datasets);
    for (int s = 0; s < n_datasets; s++) {
        char *stem = dataset_stem(datasets[s].path);
        char name[512];
        snprintf(name, sizeof(name), "%s_chain%02d.csv", stem, chain_id);
        free(stem);
        char *path = join_path(output_dir, name);
        logs.files[s] = fopen(path, "w");
        if (logs.files[s] == NULL) die("cannot create '%s'", path);
        free(path);
        fprintf(logs.files[s],
                "iteration,k_min,chordal_dist_max_as_prcnt_of_range,d_max,G_max,RMSE_change_max,"
                "radius_percentile,min_stable_iterations,"
                "L_A,L_D,L_C,L_E,F,L,T,accepted\n");
    }
    return logs;
}

static void close_chain_logs(chain_logs_t *logs, int n_datasets) {
    for (int s = 0; s < n_datasets; s++) fclose(logs->files[s]);
    free(logs->files);
}

static void log_iteration(chain_logs_t *logs, int dataset_idx, long iteration, const double *phi,
                          const losses_t *l, double L, double T, int accepted) {
    FILE *f = logs->files[dataset_idx];
    fprintf(f, "%ld,%d,%.10g,%d,%.10g,%.10g,%.10g,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%d\n", iteration,
            (int)lround(phi[PHI_K_MIN]), phi[PHI_CHORDAL], (int)lround(phi[PHI_D_MAX]), phi[PHI_G_MAX],
            phi[PHI_RMSE_CHANGE_MAX], phi[PHI_RADIUS_PERCENTILE], (int)lround(phi[PHI_MIN_STABLE_ITERATIONS]),
            l->L_A, l->L_D, l->L_C, l->L_E, l->F, L, T, accepted);
}

/* ==== simulated annealing loop, one chain ========================================== */

typedef struct {
    double T0;
    double alpha;
    double q_min;
    double q_max;
    long max_iterations;
    double min_temperature;
} sa_params_t;

/* One phi's full evaluation against every dataset in the suite: the raw STC run and the
 * resulting losses, kept per dataset (not just the aggregate) so a chain can cache its
 * *current* state across iterations and never re-run STC just to log or shatter-dump a state
 * that was already computed -- see run_chain below. */
typedef struct {
    stc_run_t *runs;     /* n_datasets */
    losses_t *losses;    /* n_datasets */
    double L;            /* mean combined loss, mod_STC.md "Aggregation across a benchmark suite" */
} phi_evaluation_t;

static phi_evaluation_t evaluate_phi_fresh(const double *phi, const dataset_t *datasets, int n_datasets,
                                           const fixed_params_t *fx) {
    phi_evaluation_t ev;
    ev.runs = xmalloc(sizeof(stc_run_t) * (size_t)n_datasets);
    ev.losses = xmalloc(sizeof(losses_t) * (size_t)n_datasets);
    double sum_L = 0.0;
    for (int s = 0; s < n_datasets; s++) {
        ev.runs[s] = run_stc(&datasets[s], phi, fx);
        ev.losses[s] = compute_losses(&datasets[s], &ev.runs[s]);
        sum_L += combined_loss(&ev.losses[s]);
    }
    ev.L = sum_L / (double)n_datasets;
    return ev;
}

static void free_phi_evaluation(phi_evaluation_t *ev, int n_datasets) {
    for (int s = 0; s < n_datasets; s++) free_stc_run(&ev->runs[s]);
    free(ev->runs);
    free(ev->losses);
}

/* Logs `ev` (an already-computed evaluation -- never recomputes anything) as the chain's state
 * at `iteration`, one CSV row per dataset, and -- if `store_shatter` -- the full per-iteration
 * JSON+CSV/TSV dump alongside it. Called once per iteration for whatever the chain's *current*
 * state actually is (the accepted proposal, or the retained previous state on a rejection) --
 * never for a proposal that gets rejected, which never becomes part of the logged trajectory. */
static void log_and_maybe_shatter(chain_logs_t *logs, const dataset_t *datasets, int n_datasets,
                                  long iteration, const double *phi, const phi_evaluation_t *ev, double T,
                                  int accepted, int store_shatter, const char *output_dir, int chain_id,
                                  const fixed_params_t *fx) {
    for (int s = 0; s < n_datasets; s++) {
        log_iteration(logs, s, iteration, phi, &ev->losses[s], ev->L, T, accepted);
        if (store_shatter) {
            char *stem = dataset_stem(datasets[s].path);
            char chain_dir_name[256];
            snprintf(chain_dir_name, sizeof(chain_dir_name), "%s_chain%02d", stem, chain_id);
            free(stem);
            char *chain_dir = join_path(output_dir, chain_dir_name);
            mkdir_or_die(chain_dir);
            char iter_dir_name[32];
            snprintf(iter_dir_name, sizeof(iter_dir_name), "iter_%06ld", iteration);
            char *iter_dir = join_path(chain_dir, iter_dir_name);
            mkdir_or_die(iter_dir);
            write_shatter_iteration_output(&datasets[s], &ev->runs[s], phi, fx, iter_dir);
            free(iter_dir);
            free(chain_dir);
        }
    }
}

typedef struct {
    double phi[PHI_P];
    double L;
} chain_result_t;

static chain_result_t run_chain(int chain_id, const double *initial_phi, const dataset_t *datasets,
                                int n_datasets, const fixed_params_t *fx, const phi_bounds_t *bounds,
                                const sa_params_t *sa, const char *output_dir, uint64_t seed,
                                int store_shatter) {
    rng_state rng;
    rng_seed_for_chain(&rng, seed, chain_id);

    chain_logs_t logs = open_chain_logs(output_dir, datasets, n_datasets, chain_id);

    double phi[PHI_P], best_phi[PHI_P];
    memcpy(phi, initial_phi, sizeof(phi));
    double T = sa->T0;

    /* `current` is the chain's own state, evaluated exactly once per phi it ever actually
     * holds -- never recomputed just to log or shatter-dump it. A rejected proposal's
     * evaluation is freed immediately below without ever being logged; an accepted one
     * replaces `current` outright, so the buffer that gets logged this iteration is always
     * exactly the buffer that was already computed to make the accept/reject decision. */
    phi_evaluation_t current = evaluate_phi_fresh(phi, datasets, n_datasets, fx);
    log_and_maybe_shatter(&logs, datasets, n_datasets, 0, phi, &current, T, 1, store_shatter, output_dir,
                          chain_id, fx);
    memcpy(best_phi, phi, sizeof(phi));
    double best_L = current.L;

    for (long iter = 1; iter < sa->max_iterations && T > sa->min_temperature; iter++) {
        double proposal[PHI_P];
        for (int i = 0; i < PHI_P; i++) {
            double range = bounds->hi[i] - bounds->lo[i];
            double sigma = range * (sa->q_min + (sa->q_max - sa->q_min) * (T / sa->T0));
            proposal[i] = phi[i] + sigma * rng_normal(&rng);
        }
        phi_reflect(proposal, bounds);

        phi_evaluation_t candidate = evaluate_phi_fresh(proposal, datasets, n_datasets, fx);
        int accept = candidate.L <= current.L;
        if (!accept) {
            double p_accept = exp(-(candidate.L - current.L) / T);
            accept = rng_uniform01(&rng) < p_accept;
        }

        if (accept) {
            memcpy(phi, proposal, sizeof(phi));
            free_phi_evaluation(&current, n_datasets);
            current = candidate;
        } else {
            free_phi_evaluation(&candidate, n_datasets);
        }

        log_and_maybe_shatter(&logs, datasets, n_datasets, iter, phi, &current, T, accept, store_shatter,
                              output_dir, chain_id, fx);

        if (current.L < best_L) {
            best_L = current.L;
            memcpy(best_phi, phi, sizeof(phi));
        }
        T *= sa->alpha;
    }

    free_phi_evaluation(&current, n_datasets);
    close_chain_logs(&logs, n_datasets);

    chain_result_t result;
    memcpy(result.phi, best_phi, sizeof(best_phi));
    result.L = best_L;
    return result;
}

/* ==== --input-dir: list *.json files, load, validate shared ambient_dim ============ */

static char **list_json_files(const char *dir, int *n_out) {
    DIR *d = opendir(dir);
    if (d == NULL) die("cannot open --input-dir '%s'", dir);
    char **paths = NULL;
    int n = 0, cap = 0;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        size_t len = strlen(ent->d_name);
        if (len < 6 || strcmp(ent->d_name + len - 5, ".json") != 0) continue;
        if (n == cap) {
            cap = cap == 0 ? 16 : cap * 2;
            paths = realloc(paths, sizeof(char *) * (size_t)cap);
        }
        paths[n++] = join_path(dir, ent->d_name);
    }
    closedir(d);
    if (n == 0) die("--input-dir '%s' contains no *.json files", dir);
    /* sort for reproducible chain-log ordering across runs/platforms */
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (strcmp(paths[i], paths[j]) > 0) {
                char *tmp = paths[i];
                paths[i] = paths[j];
                paths[j] = tmp;
            }
    *n_out = n;
    return paths;
}

/* ==== argument parsing ============================================================= */

enum {
    OPT_INPUT_DIR = 256,
    OPT_OUTPUT_DIR,
    OPT_CHAINS,
    OPT_SEED,
    OPT_T0,
    OPT_ALPHA,
    OPT_Q_MIN,
    OPT_Q_MAX,
    OPT_MAX_ITERATIONS,
    OPT_MIN_TEMPERATURE,
    OPT_PHI_MAX_D_MAX,
    OPT_PHI_MAX_G_MAX,
    OPT_PHI_MAX_RMSE_CHANGE_MAX,
    OPT_PHI_MIN_MIN_STABLE_ITERATIONS,
    OPT_PHI_MAX_MIN_STABLE_ITERATIONS,
    OPT_O,
    OPT_F_MAX,
    OPT_BANDWIDTH_PERCENTILE,
    OPT_EXCLUSION_RADIUS_PERCENTILE,
    OPT_RECONCILIATION_MODE,
    OPT_MIN_OVERLAP_COEFFICIENT,
    OPT_N_ANCHORS,
    OPT_SEED_MAX_SET_SIZE,
    OPT_QUANTILE_PAIRWISE_EA_COMPARISON,
    OPT_STORE_ITERATION_SHATTER_RESULTS,
};

struct arguments {
    const char *input_dir;
    const char *output_dir;
    int chains, have_chains;
    uint64_t seed;
    int have_seed;

    double T0, alpha, q_min, q_max, min_temperature;
    long max_iterations;

    double phi_max_d_max;
    int have_phi_max_d_max;
    double phi_max_g_max;
    double phi_max_rmse_change_max;
    double phi_min_min_stable_iterations;
    double phi_max_min_stable_iterations;

    fixed_params_t fx;
    int have_o;

    int store_shatter;
};

static struct argp_option options[] = {
    {"input-dir", OPT_INPUT_DIR, "DIR", 0, "Directory of benchmark dataset JSON files (required)", 0},
    {"output-dir", OPT_OUTPUT_DIR, "DIR", 0, "Directory to write per-(dataset,chain) CSV logs into (required)", 0},
    {"chains", OPT_CHAINS, "N", 0,
     "Randomly-initialized chains, in addition to the always-present estimator-seeded chain "
     "(default: omp_get_max_threads() - 1)", 0},
    {"seed", OPT_SEED, "N", 0, "RNG seed (required)", 0},

    {"t0", OPT_T0, "REAL", 0, "Initial temperature (default: 0.15)", 1},
    {"alpha", OPT_ALPHA, "REAL", 0, "Geometric cooling factor, 0<alpha<1 (default: 0.995)", 1},
    {"q-min", OPT_Q_MIN, "FRACTION", 0, "Minimum relative proposal size (default: 0.0)", 1},
    {"q-max", OPT_Q_MAX, "FRACTION", 0, "Initial relative proposal size (default: 0.2)", 1},
    {"max-iterations", OPT_MAX_ITERATIONS, "N", 0, "Maximum iterations per chain (default: 2000)", 1},
    {"min-temperature", OPT_MIN_TEMPERATURE, "REAL", 0, "Stop a chain once T falls below this (default: 0.001)", 1},

    {"phi-max-d-max", OPT_PHI_MAX_D_MAX, "N", 0,
     "Upper bound for the tuned d_max (default: the suite's shared ambient_dim)", 1},
    {"phi-max-g-max", OPT_PHI_MAX_G_MAX, "REAL", 0, "Upper bound for the tuned G_max (default: 5.0)", 1},
    {"phi-max-rmse-change-max", OPT_PHI_MAX_RMSE_CHANGE_MAX, "REAL", 0,
     "Upper bound for the tuned RMSE_change_max (default: 5.0)", 1},
    {"phi-min-min-stable-iterations", OPT_PHI_MIN_MIN_STABLE_ITERATIONS, "N", 0,
     "Lower bound for the tuned min_stable_iterations (default: 1)", 1},
    {"phi-max-min-stable-iterations", OPT_PHI_MAX_MIN_STABLE_ITERATIONS, "N", 0,
     "Upper bound for the tuned min_stable_iterations (default: 5)", 1},

    {"o", OPT_O, "N", 0, "Fixed (non-tuned) STC parameter: observable-history window depth (required)", 2},
    {"f-max", OPT_F_MAX, "FRACTION", 0, "Fixed (non-tuned) STC parameter (default: 1.0, i.e. disabled)", 2},
    {"bandwidth-percentile", OPT_BANDWIDTH_PERCENTILE, "PERCENTILE", 0, "Fixed (non-tuned) STC parameter (default: 68.27)", 2},
    {"exclusion-radius-percentile", OPT_EXCLUSION_RADIUS_PERCENTILE, "PERCENTILE", 0, "Fixed (non-tuned) STC parameter (default: 50.0)", 2},
    {"reconciliation-mode", OPT_RECONCILIATION_MODE, "MODE", 0, "Fixed (non-tuned) STC parameter (default: merge_overlap_coefficient)", 2},
    {"min-overlap-coefficient", OPT_MIN_OVERLAP_COEFFICIENT, "FRACTION", 0, "Fixed (non-tuned) STC parameter (default: 0.9)", 2},

    {"n-anchors", OPT_N_ANCHORS, "N", 0, "estimate_stc_parameters parameter, for the estimator-seeded chain only (default: 5)", 3},
    {"seed-max-set-size", OPT_SEED_MAX_SET_SIZE, "PERCENTILE", 0, "estimate_stc_parameters parameter (default: 5.0)", 3},
    {"quantile-pairwise-ea-comparison", OPT_QUANTILE_PAIRWISE_EA_COMPARISON, "PERCENTILE", 0, "estimate_stc_parameters parameter (default: 25.0)", 3},

    {"store-iteration-shatter-results", OPT_STORE_ITERATION_SHATTER_RESULTS, 0, 0,
     "Also write the full per-iteration JSON+CSV/TSV STC output (never the HTML report) into "
     "--output-dir/{dataset_stem}_chain{N}/iter_{NNNNNN}/ for every (chain, dataset, iteration) "
     "-- off by default; substantial extra disk space and I/O time for a long run (see "
     "mod_STC.md's \"Full per-iteration STC output\")", 4},

    {0},
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *args = state->input;
    switch (key) {
        case OPT_INPUT_DIR: args->input_dir = arg; break;
        case OPT_OUTPUT_DIR: args->output_dir = arg; break;
        case OPT_CHAINS: args->chains = atoi(arg); args->have_chains = 1; break;
        case OPT_SEED: args->seed = strtoull(arg, NULL, 10); args->have_seed = 1; break;

        case OPT_T0: args->T0 = atof(arg); break;
        case OPT_ALPHA: args->alpha = atof(arg); break;
        case OPT_Q_MIN: args->q_min = atof(arg); break;
        case OPT_Q_MAX: args->q_max = atof(arg); break;
        case OPT_MAX_ITERATIONS: args->max_iterations = atol(arg); break;
        case OPT_MIN_TEMPERATURE: args->min_temperature = atof(arg); break;

        case OPT_PHI_MAX_D_MAX: args->phi_max_d_max = atof(arg); args->have_phi_max_d_max = 1; break;
        case OPT_PHI_MAX_G_MAX: args->phi_max_g_max = atof(arg); break;
        case OPT_PHI_MAX_RMSE_CHANGE_MAX: args->phi_max_rmse_change_max = atof(arg); break;
        case OPT_PHI_MIN_MIN_STABLE_ITERATIONS: args->phi_min_min_stable_iterations = atof(arg); break;
        case OPT_PHI_MAX_MIN_STABLE_ITERATIONS: args->phi_max_min_stable_iterations = atof(arg); break;

        case OPT_O: args->fx.o = atoi(arg); args->have_o = 1; break;
        case OPT_F_MAX: args->fx.f_max = atof(arg); break;
        case OPT_BANDWIDTH_PERCENTILE: args->fx.bandwidth_percentile = atof(arg); break;
        case OPT_EXCLUSION_RADIUS_PERCENTILE: args->fx.exclusion_radius_percentile = atof(arg); break;
        case OPT_RECONCILIATION_MODE: args->fx.reconciliation_mode = arg; break;
        case OPT_MIN_OVERLAP_COEFFICIENT: args->fx.min_overlap_coefficient = atof(arg); break;

        case OPT_N_ANCHORS: args->fx.n_anchors = atoi(arg); break;
        case OPT_SEED_MAX_SET_SIZE: args->fx.seed_max_set_size = atof(arg); break;
        case OPT_QUANTILE_PAIRWISE_EA_COMPARISON: args->fx.quantile_pairwise_ea_comparison = atof(arg); break;
        case OPT_STORE_ITERATION_SHATTER_RESULTS: args->store_shatter = 1; break;

        case ARGP_KEY_END:
            if (args->input_dir == NULL) argp_error(state, "--input-dir is required");
            if (args->output_dir == NULL) argp_error(state, "--output-dir is required");
            if (!args->have_seed) argp_error(state, "--seed is required");
            if (!args->have_o) argp_error(state, "--o is required");
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, 0,
                           "Simulated annealing over Shape Truthful Clustering parameters, scored "
                           "against known-manifold benchmark datasets (misc/mod_STC.md).",
                           0, 0, 0};

/* ==== main ========================================================================== */

int main(int argc, char **argv) {
    struct arguments args;
    memset(&args, 0, sizeof(args));
    args.T0 = 0.15;
    args.alpha = 0.995;
    args.q_min = 0.0;
    args.q_max = 0.2;
    args.max_iterations = 2000;
    args.min_temperature = 0.001;
    args.phi_max_g_max = 5.0;
    args.phi_max_rmse_change_max = 5.0;
    args.phi_min_min_stable_iterations = 1;
    args.phi_max_min_stable_iterations = 5;
    args.fx.f_max = 1.0; /* disabled by default: a dataset may legitimately be one cluster */
    args.fx.bandwidth_percentile = 68.27;
    args.fx.exclusion_radius_percentile = 50.0;
    args.fx.reconciliation_mode = "merge_overlap_coefficient";
    args.fx.min_overlap_coefficient = 0.9;
    args.fx.n_anchors = 5;
    args.fx.seed_max_set_size = 5.0;
    args.fx.quantile_pairwise_ea_comparison = 25.0;

    argp_parse(&argp, argc, argv, 0, 0, &args);

    struct stat sb;
    if (stat(args.output_dir, &sb) != 0 || !S_ISDIR(sb.st_mode))
        die("--output-dir '%s' does not exist or is not a directory", args.output_dir);

    /* Deterministic regardless of OMP_NESTED/OMP_MAX_ACTIVE_LEVELS in the environment: the
     * outer parallel-for over chains already saturates every core one-to-one, so any nested
     * parallel region STC's own `_impl` kernels start must run single-threaded, not rely on
     * a default that a caller's environment could override -- see mod_STC.md "Parallelization". */
    omp_set_max_active_levels(1);

    if (!args.have_chains) args.chains = omp_get_max_threads() > 1 ? omp_get_max_threads() - 1 : 0;
    if (args.chains < 0) die("--chains must be >= 0");

    int n_datasets;
    char **paths = list_json_files(args.input_dir, &n_datasets);
    dataset_t *datasets = xmalloc(sizeof(dataset_t) * (size_t)n_datasets);
    for (int s = 0; s < n_datasets; s++) {
        datasets[s] = load_dataset(paths[s]);
        if (s > 0 && datasets[s].ambient_dim != datasets[0].ambient_dim) {
            die("--input-dir mixes ambient_dim=%d ('%s') with ambient_dim=%d ('%s') -- every "
                "dataset in one run must share the same ambient_dim (see mod_STC.md)",
                datasets[0].ambient_dim, paths[0], datasets[s].ambient_dim, paths[s]);
        }
    }
    printf("stc_sa_cli: loaded %d dataset(s) from '%s' (ambient_dim=%d)\n", n_datasets, args.input_dir,
           datasets[0].ambient_dim);

    /* ---- phi bounds, shared across the suite ---- */
    phi_bounds_t bounds;
    long min_n_vectors = datasets[0].n_vectors;
    for (int s = 1; s < n_datasets; s++)
        if (datasets[s].n_vectors < min_n_vectors) min_n_vectors = datasets[s].n_vectors;
    bounds.lo[PHI_K_MIN] = 1;
    bounds.hi[PHI_K_MIN] = (double)(min_n_vectors - 1);
    bounds.lo[PHI_CHORDAL] = 0.0;
    bounds.hi[PHI_CHORDAL] = 1.0;
    bounds.lo[PHI_D_MAX] = 0;
    bounds.hi[PHI_D_MAX] = args.have_phi_max_d_max ? args.phi_max_d_max : (double)datasets[0].ambient_dim;
    bounds.lo[PHI_G_MAX] = 0.0;
    bounds.hi[PHI_G_MAX] = args.phi_max_g_max;
    bounds.lo[PHI_RMSE_CHANGE_MAX] = 0.0;
    bounds.hi[PHI_RMSE_CHANGE_MAX] = args.phi_max_rmse_change_max;
    bounds.lo[PHI_RADIUS_PERCENTILE] = 0.0;
    bounds.hi[PHI_RADIUS_PERCENTILE] = 100.0;
    bounds.lo[PHI_MIN_STABLE_ITERATIONS] = args.phi_min_min_stable_iterations;
    bounds.hi[PHI_MIN_STABLE_ITERATIONS] = args.phi_max_min_stable_iterations;

    /* ---- estimator-seeded chain's initial phi: estimate_stc_parameters per dataset, averaged ----
     * A per-dataset failure here is non-fatal: estimate_stc_parameters is an explicitly heuristic,
     * best-effort starting point (mod_STC.md: "a reasonable value to start from and refine, not a
     * guarantee") -- with the default n_anchors=5, a high-noise or mixed-structure dataset can
     * occasionally leave fewer than 2 usable anchor clouds, or no anchor pair sharing a nonzero
     * tangent dimension (ERR_INTERNAL, a clean documented failure path, not a crash). One
     * dataset's estimator failing should not abort the whole suite's SA run -- skip it from the
     * average; if every dataset fails, fall back to the midpoint of each phi bound. */
    double estimator_phi[PHI_P] = {0};
    int n_estimated = 0;
    for (int s = 0; s < n_datasets; s++) {
        double est_k_min, est_density_quantile, est_chordal, est_G_max, est_d_max;
        int ierr = 0;
        estimate_stc_parameters_c(datasets[s].vectors, &datasets[s].ambient_dim, &datasets[s].n_vectors,
                                  datasets[s].kd_indices, datasets[s].dimension_order,
                                  NULL /* k_min: kernel default, matching stc_cli's own --estimate-parameters */,
                                  &args.fx.bandwidth_percentile, &args.fx.n_anchors, &args.fx.seed_max_set_size,
                                  &args.fx.quantile_pairwise_ea_comparison, &est_k_min,
                                  &est_density_quantile, &est_chordal, &est_G_max, &est_d_max, &ierr);
        if (ierr != 0) {
            fprintf(stderr, "stc_sa_cli: warning: estimate_stc_parameters failed for '%s' (ierr=%d) -- "
                            "excluded from the estimator-seeded chain's starting phi\n",
                    datasets[s].path, ierr);
            continue;
        }
        estimator_phi[PHI_K_MIN] += est_k_min;
        estimator_phi[PHI_CHORDAL] += est_chordal;
        estimator_phi[PHI_D_MAX] += est_d_max;
        estimator_phi[PHI_G_MAX] += est_G_max;
        n_estimated++;
    }
    if (n_estimated > 0) {
        estimator_phi[PHI_K_MIN] /= (double)n_estimated;
        estimator_phi[PHI_CHORDAL] /= (double)n_estimated;
        estimator_phi[PHI_D_MAX] /= (double)n_estimated;
        estimator_phi[PHI_G_MAX] /= (double)n_estimated;
    } else {
        fprintf(stderr, "stc_sa_cli: warning: estimate_stc_parameters failed for every dataset -- "
                        "starting the estimator-seeded chain at the midpoint of each phi range instead\n");
        estimator_phi[PHI_K_MIN] = 0.5 * (bounds.lo[PHI_K_MIN] + bounds.hi[PHI_K_MIN]);
        estimator_phi[PHI_CHORDAL] = 0.5 * (bounds.lo[PHI_CHORDAL] + bounds.hi[PHI_CHORDAL]);
        estimator_phi[PHI_D_MAX] = 0.5 * (bounds.lo[PHI_D_MAX] + bounds.hi[PHI_D_MAX]);
        estimator_phi[PHI_G_MAX] = 0.5 * (bounds.lo[PHI_G_MAX] + bounds.hi[PHI_G_MAX]);
    }
    /* No estimator output at all for these three -- always start at the midpoint of their own
     * phi range, see the header comment above. */
    estimator_phi[PHI_RMSE_CHANGE_MAX] = 0.5 * (bounds.lo[PHI_RMSE_CHANGE_MAX] + bounds.hi[PHI_RMSE_CHANGE_MAX]);
    estimator_phi[PHI_RADIUS_PERCENTILE] = 0.5 * (bounds.lo[PHI_RADIUS_PERCENTILE] + bounds.hi[PHI_RADIUS_PERCENTILE]);
    estimator_phi[PHI_MIN_STABLE_ITERATIONS] =
        0.5 * (bounds.lo[PHI_MIN_STABLE_ITERATIONS] + bounds.hi[PHI_MIN_STABLE_ITERATIONS]);
    phi_reflect(estimator_phi, &bounds);

    sa_params_t sa = {args.T0, args.alpha, args.q_min, args.q_max, args.max_iterations, args.min_temperature};

    int total_chains = 1 + args.chains;
    chain_result_t *results = xmalloc(sizeof(chain_result_t) * (size_t)total_chains);

    printf("stc_sa_cli: running %d chain(s) (1 estimator-seeded + %d random) for up to %ld iterations each\n",
           total_chains, args.chains, args.max_iterations);

#pragma omp parallel for schedule(dynamic)
    for (int c = 0; c < total_chains; c++) {
        double initial_phi[PHI_P];
        if (c == 0) {
            memcpy(initial_phi, estimator_phi, sizeof(initial_phi));
        } else {
            rng_state init_rng;
            rng_seed_for_chain(&init_rng, args.seed, c + total_chains /* distinct stream from the chain's own */);
            for (int i = 0; i < PHI_P; i++) initial_phi[i] = rng_uniform(&init_rng, bounds.lo[i], bounds.hi[i]);
        }
        results[c] = run_chain(c, initial_phi, datasets, n_datasets, &args.fx, &bounds, &sa, args.output_dir,
                               args.seed, args.store_shatter);
    }

    int best_chain = 0;
    for (int c = 1; c < total_chains; c++)
        if (results[c].L < results[best_chain].L) best_chain = c;

    printf("stc_sa_cli: best chain %d (%s), L=%.6g, phi=(k_min=%d, chordal=%.6g, d_max=%d, G_max=%.6g, "
           "RMSE_change_max=%.6g, radius_percentile=%.6g, min_stable_iterations=%d)\n",
           best_chain, best_chain == 0 ? "estimator-seeded" : "random", results[best_chain].L,
           (int)lround(results[best_chain].phi[PHI_K_MIN]),
           results[best_chain].phi[PHI_CHORDAL], (int)lround(results[best_chain].phi[PHI_D_MAX]),
           results[best_chain].phi[PHI_G_MAX], results[best_chain].phi[PHI_RMSE_CHANGE_MAX],
           results[best_chain].phi[PHI_RADIUS_PERCENTILE],
           (int)lround(results[best_chain].phi[PHI_MIN_STABLE_ITERATIONS]));
    printf("stc_sa_cli: estimator-seeded chain (chain 0) achieved L=%.6g -- compare against the best "
           "random chain to judge whether estimate_stc_parameters is already landing near the optimum\n",
           results[0].L);

    return 0;
}
