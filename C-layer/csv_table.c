#include "csv_table.h"

#include <csv.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Generous for a numeric-cell CSV (no realistic real-number literal comes close); a fixed
 * stack buffer keeps the per-field hot path allocation-free. */
#define CSV_TABLE_MAX_FIELD 1024

static char g_last_error[512] = "";

const char *csv_table_last_error(void) {
    return g_last_error;
}

static void set_error(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(g_last_error, sizeof(g_last_error), fmt, ap);
    va_end(ap);
}

/* ==== Phase 1: how many columns does the first line have? ======================== */
/* Parsed as its own small, throwaway csv_parser instance over just the first line's bytes --
 * cheaper than instantiating the real output buffer before its shape is even known, and
 * exactly what misc/mod_STC.md's own algorithm description asks for. */

struct count_state {
    int n_fields;
};

static void count_field_cb(void *s, size_t len, void *data) {
    (void)s;
    (void)len;
    struct count_state *st = (struct count_state *)data;
    st->n_fields++;
}

static void count_row_cb(int c, void *data) {
    (void)c;
    (void)data;
    /* first-line-only parse never needs more than one row; nothing to do here */
}

static int count_columns(const char *buf, size_t buf_len, char separator, int *out_n_dimensions) {
    const char *newline = memchr(buf, '\n', buf_len);
    size_t first_line_len = newline ? (size_t)(newline - buf) : buf_len;

    struct csv_parser parser;
    if (csv_init(&parser, 0) != 0) {
        set_error("csv_table: csv_init failed");
        return -1;
    }
    csv_set_delim(&parser, (unsigned char)separator);

    struct count_state st = {0};
    size_t consumed = csv_parse(&parser, buf, first_line_len, count_field_cb, count_row_cb, &st);
    int failed = (consumed != first_line_len);
    if (!failed && csv_fini(&parser, count_field_cb, count_row_cb, &st) != 0) {
        failed = 1;
    }
    if (failed) {
        set_error("csv_table: %s", csv_strerror(csv_error(&parser)));
    }
    csv_free(&parser);
    if (failed) return -1;

    if (st.n_fields < 1) {
        set_error("csv_table: the first line has no columns");
        return -1;
    }
    *out_n_dimensions = st.n_fields;
    return 0;
}

/* ==== Phase 2: the real, whole-file parse ========================================= */

struct parse_state {
    int n_dimensions;
    int n_records;
    int header;
    double *data;
    char **dim_names;
    int cur_row;    /* 0-based; counts the header row too, if present */
    int cur_col;    /* 0-based, within the current row */
    int data_row;   /* 0-based index into `data`'s row dimension; only advances for data rows */
    int error;
    char field_buf[CSV_TABLE_MAX_FIELD];
};

static void parse_field_cb(void *s, size_t len, void *data) {
    struct parse_state *st = (struct parse_state *)data;
    if (st->error) return;

    if (len >= sizeof(st->field_buf)) {
        st->error = 1;
        set_error("csv_table: field at row %d, column %d is longer than %zu bytes -- not a "
                  "real number", st->cur_row + 1, st->cur_col + 1, sizeof(st->field_buf) - 1);
        return;
    }
    memcpy(st->field_buf, s, len);
    st->field_buf[len] = '\0';

    if (st->cur_row == 0 && st->header) {
        if (st->cur_col >= st->n_dimensions) {
            st->error = 1;
            set_error("csv_table: header row has more than %d column(s)", st->n_dimensions);
            return;
        }
        st->dim_names[st->cur_col] = strdup(st->field_buf);
        if (st->dim_names[st->cur_col] == NULL) {
            st->error = 1;
            set_error("csv_table: out of memory copying a header name");
            return;
        }
    } else {
        if (st->cur_col >= st->n_dimensions) {
            st->error = 1;
            set_error("csv_table: row %d has more than %d column(s)", st->cur_row + 1, st->n_dimensions);
            return;
        }
        if (st->data_row >= st->n_records) {
            st->error = 1;
            set_error("csv_table: found more data rows than n_records=%d allows", st->n_records);
            return;
        }

        char *end = NULL;
        errno = 0;
        double val = strtod(st->field_buf, &end);
        if (end == st->field_buf || *end != '\0' || errno == ERANGE) {
            st->error = 1;
            set_error("csv_table: row %d, column %d ('%s') is not a real number",
                      st->cur_row + 1, st->cur_col + 1, st->field_buf);
            return;
        }
        st->data[(size_t)st->data_row * (size_t)st->n_dimensions + (size_t)st->cur_col] = val;
    }

    st->cur_col++;
}

static void parse_row_cb(int c, void *data) {
    (void)c;
    struct parse_state *st = (struct parse_state *)data;
    if (st->error) return;

    if (st->cur_col != st->n_dimensions) {
        st->error = 1;
        set_error("csv_table: row %d has %d column(s), expected %d",
                  st->cur_row + 1, st->cur_col, st->n_dimensions);
        return;
    }

    if (!(st->cur_row == 0 && st->header)) {
        st->data_row++;
    }
    st->cur_row++;
    st->cur_col = 0;
}

static int read_whole_file(const char *path, char **out_buf, size_t *out_len) {
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) {
        set_error("csv_table: cannot open '%s': %s", path, strerror(errno));
        return -1;
    }
    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        set_error("csv_table: cannot seek '%s'", path);
        return -1;
    }
    long size = ftell(fp);
    if (size < 0) {
        fclose(fp);
        set_error("csv_table: cannot determine the size of '%s'", path);
        return -1;
    }
    rewind(fp);

    char *buf = NULL;
    if (size > 0) {
        buf = malloc((size_t)size);
        if (buf == NULL) {
            fclose(fp);
            set_error("csv_table: out of memory reading '%s' (%ld bytes)", path, size);
            return -1;
        }
    }
    size_t read_bytes = (size > 0) ? fread(buf, 1, (size_t)size, fp) : 0;
    fclose(fp);
    if (read_bytes != (size_t)size) {
        free(buf);
        set_error("csv_table: short read on '%s'", path);
        return -1;
    }
    if (read_bytes == 0) {
        free(buf);
        set_error("csv_table: '%s' is empty", path);
        return -1;
    }

    *out_buf = buf;
    *out_len = read_bytes;
    return 0;
}

int csv_table_read(const char *path, char separator, int header, int n_records, csv_table *out) {
    g_last_error[0] = '\0';
    memset(out, 0, sizeof(*out));

    if (n_records < 1) {
        set_error("csv_table: n_records must be at least 1, got %d", n_records);
        return -1;
    }

    char *buf = NULL;
    size_t buf_len = 0;
    if (read_whole_file(path, &buf, &buf_len) != 0) return -1;

    int n_dimensions = 0;
    if (count_columns(buf, buf_len, separator, &n_dimensions) != 0) {
        free(buf);
        return -1;
    }

    double *data = malloc(sizeof(double) * (size_t)n_dimensions * (size_t)n_records);
    if (data == NULL) {
        free(buf);
        set_error("csv_table: out of memory allocating a %d x %d table", n_dimensions, n_records);
        return -1;
    }
    char **dim_names = NULL;
    if (header) {
        dim_names = calloc((size_t)n_dimensions, sizeof(char *));
        if (dim_names == NULL) {
            free(buf);
            free(data);
            set_error("csv_table: out of memory allocating dim_names");
            return -1;
        }
    }

    struct csv_parser parser;
    if (csv_init(&parser, 0) != 0) {
        free(buf);
        free(data);
        free(dim_names);
        set_error("csv_table: csv_init failed");
        return -1;
    }
    csv_set_delim(&parser, (unsigned char)separator);

    struct parse_state st;
    memset(&st, 0, sizeof(st));
    st.n_dimensions = n_dimensions;
    st.n_records = n_records;
    st.header = header;
    st.data = data;
    st.dim_names = dim_names;

    size_t consumed = csv_parse(&parser, buf, buf_len, parse_field_cb, parse_row_cb, &st);
    if (!st.error && consumed != buf_len) {
        st.error = 1;
        set_error("csv_table: %s", csv_strerror(csv_error(&parser)));
    }
    if (!st.error && csv_fini(&parser, parse_field_cb, parse_row_cb, &st) != 0 && !st.error) {
        st.error = 1;
        set_error("csv_table: %s", csv_strerror(csv_error(&parser)));
    }
    csv_free(&parser);
    free(buf);

    if (!st.error && st.data_row != n_records) {
        st.error = 1;
        set_error("csv_table: found %d data row(s), but n_records=%d", st.data_row, n_records);
    }

    if (st.error) {
        free(data);
        if (dim_names != NULL) {
            for (int i = 0; i < n_dimensions; i++) free(dim_names[i]);
            free(dim_names);
        }
        return -1;
    }

    out->data = data;
    out->dim_names = dim_names;
    out->n_dimensions = n_dimensions;
    out->n_records = n_records;
    return 0;
}

void csv_table_free(csv_table *t) {
    if (t == NULL) return;
    free(t->data);
    if (t->dim_names != NULL) {
        for (int i = 0; i < t->n_dimensions; i++) free(t->dim_names[i]);
        free(t->dim_names);
    }
    memset(t, 0, sizeof(*t));
}
