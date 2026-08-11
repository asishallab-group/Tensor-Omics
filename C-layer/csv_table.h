#ifndef TENSOR_OMICS_CSV_TABLE_H
#define TENSOR_OMICS_CSV_TABLE_H

/* Reads an all-real-number CSV into a flat buffer laid out exactly as Fortran's own
 * column-major `vectors(n_dimensions, n_records)` expects -- i.e. one input row's cells are
 * contiguous, in column order, so `data` can be handed to a `_c` binding with zero copy. See
 * misc/mod_STC.md, "Command line interface (CLI) in C", "Real 2D array from a CSV with
 * libcsv" for the algorithm this implements.
 */

typedef struct {
    double *data;
        /* n_dimensions * n_records doubles; data[record * n_dimensions + dim] is row
         * `record` (0-based), column `dim` (0-based) -- equivalently Fortran's
         * vectors(dim + 1, record + 1). Owned by this struct; free via csv_table_free. */
    char **dim_names;
        /* n_dimensions null-terminated strings taken from the header row, or NULL if `header`
         * was 0 in the csv_table_read call that produced this table. Owned by this struct. */
    int n_dimensions;
    int n_records;
} csv_table;

/* Reads `path` as a CSV of exclusively real-number cells into `out` (zeroed and then filled
 * on success; unspecified on failure -- do not call csv_table_free on a failed table).
 *
 * `separator` is the column delimiter. If `header` is nonzero, the first line is taken as
 * column names (populating `dim_names`), not data. `n_records` is the exact number of DATA
 * rows expected (excluding the header line, if any) -- libcsv has no cheap way to report a
 * file's line count up front, so, per misc/mod_STC.md, this is the caller's job; a mismatch
 * between `n_records` and what the file actually contains is reported as an error, not
 * silently truncated or overflowed.
 *
 * Returns 0 on success, nonzero on failure; on failure, csv_table_last_error() describes why. */
int csv_table_read(const char *path, char separator, int header, int n_records, csv_table *out);

/* Frees everything csv_table_read allocated into `t` (data, dim_names, and each dim_names
 * string), and zeroes `*t`. Safe to call on a zero-initialized (never successfully read)
 * table; a no-op then. */
void csv_table_free(csv_table *t);

/* The message belonging to the most recent failing csv_table_read call on this thread.
 * Overwritten by the next csv_table_read call; copy it out before making another call if it
 * needs to outlive that. */
const char *csv_table_last_error(void);

#endif
