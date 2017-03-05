/* CML 1.0.0 -- A header-only C matrix library.
 *
 * ABOUT: 
 *      CML is a (fully C++ compatible) header-only C matrix library designed with a focus on portability, simplicity, 
 *      and efficiency written under the POSIX standard. Several computationally-intense functions require/are enhanced 
 *      by BLAS/LAPACK. The user is able to define optional flags prior to including this library in source code to do 
 *      things such as removing BLAS/LAPACK as a dependency (while simultaneously removing several functions from API), 
 *   	or giving all functions static storage. CML operates on type 'MATRIX' for all processes. 'MATRIX' is column major 
 *   	so the user can easily interface with the most popular linear-algebra-related API's. 'errno' is used for error 
 *   	handling extensively. 'errno' is set only in the presence of an error. The user is free to disregard 'errno' if 
 *   	an error occurred, but the state of the resultant matrix is undefined. 
 *
 * USAGE:
 *      CML can be used in either header-only mode or implementation mode. The header-only mode is used by default when
 *      included and allows including this header in other headers and does not contain the actual implementation. The
 *      implementation mode requires the user to define the preprocessor macro 'CML_IMPLEMENTATION' in one .c/.cpp file
 *      before #including this file, e.g.:
 *          
 *          #define CML_IMPLEMENTATION
 *          #include "cml.h"
 *
 *      Also optionally define the flags listed in the section 'OPTIONAL DEFINES' below in header and implementation
 *      mode if you want to use additional functionality or need more control over the library. Everytime you include 
 *      "cml.h", you have to define the same optional flags. Not doing so will lead to compile-time errors.
 *
 * OPTIONAL DEFINES:
 *      CML_PRIVATE
 *          If defined, declares all functions as static, so they can only be accessed inside the file that contains the
 *          implementation.
 *
 *      CML_NO_DEPENDENCIES
 *          If defined, will remove all functions that fully depend on external API's. The user will then be able to
 *          compile and link without the external API's.
 */
#ifndef CML_H
#define CML_H


#ifdef __cplusplus
extern "C" {
#endif


#include <stddef.h>
#include <stdbool.h>


#ifdef CML_PRIVATE
    #define CML_API static inline
#else
    #define CML_API extern
#endif


/* ========================================================
 *
 *                          API
 *
 * ======================================================== */
typedef struct {
    size_t rows, cols;
    double *data;
} MATRIX;


/*    ALLOCATION    */
CML_API     MATRIX*     cml_new(size_t rows, size_t cols);
CML_API     MATRIX*     cml_dup(MATRIX *m);
CML_API     MATRIX*     cml_ones(size_t rows, size_t cols);
CML_API     MATRIX*     cml_identity(size_t dim);
CML_API     MATRIX*     cml_lower_tri(size_t dim);
CML_API     MATRIX*     cml_upper_tri(size_t dim);
CML_API     void        cml_free(MATRIX *m);


/*    GET/SET    */
CML_API     double      cml_get(MATRIX *m, size_t row, size_t col);
CML_API     void        cml_set(MATRIX *m, size_t row, size_t col, double data);
CML_API     void        cml_set_all(MATRIX *m, double data);
CML_API     void        cml_set_row(MATRIX *m, size_t row, double data);
CML_API     void        cml_set_col(MATRIX *m, size_t col, double data);


/*    PHYSICAL MANIPULATION    */
CML_API     void        cml_cpy(MATRIX *dest, MATRIX *src);
CML_API     void        cml_cpy_row(MATRIX *dest, MATRIX *src, size_t dest_row, size_t src_row);
CML_API     void        cml_cpy_col(MATRIX *dest, MATRIX *src, size_t dest_col, size_t src_col);
CML_API     void        cml_cpy_elem(MATRIX *dest, MATRIX *src, size_t dest_row, size_t dest_col, size_t src_row, size_t src_col);
CML_API     void        cml_cpy_self_row(MATRIX *m, size_t dest_row, size_t src_row);
CML_API     void        cml_cpy_self_col(MATRIX *m, size_t dest_col, size_t src_col);
CML_API     void        cml_cpy_self_elem(MATRIX *m, size_t dest_row, size_t dest_col, size_t src_row, size_t src_col);
CML_API     void        cml_swap(MATRIX *m1, MATRIX *m2);
CML_API     void        cml_swap_row(MATRIX *m1, MATRIX *m2, size_t m1_row, size_t m2_row);
CML_API     void        cml_swap_col(MATRIX *m1, MATRIX *m2, size_t m1_col, size_t m2_col);
CML_API     void        cml_swap_elem(MATRIX *m1, MATRIX *m2, size_t m1_row, size_t m1_col, size_t m2_row, size_t m2_col);
CML_API     void        cml_swap_self_row(MATRIX *m, size_t row_1, size_t row_2);
CML_API     void        cml_swap_self_col(MATRIX *m, size_t col_1, size_t col_2);
CML_API     void        cml_swap_self_elem(MATRIX *m, size_t row_1, size_t col_1, size_t row_2, size_t col_2);
CML_API     void        cml_ins_row(MATRIX *dest, MATRIX *src, size_t dest_row, size_t src_row);
CML_API     void        cml_ins_col(MATRIX *dest, MATRIX *src, size_t dest_col, size_t src_col);
CML_API     void        cml_ins_self_row(MATRIX *m, size_t dest_row, size_t src_row);
CML_API     void        cml_ins_self_col(MATRIX *m, size_t dest_col, size_t src_col);
CML_API     void        cml_del_all(MATRIX *m);
CML_API     void        cml_del_row(MATRIX *m, size_t row);
CML_API     void        cml_del_col(MATRIX *m, size_t col);
CML_API     void        cml_adjoin_top(MATRIX *dest, MATRIX *src);
CML_API     void        cml_adjoin_bottom(MATRIX *dest, MATRIX *src);
CML_API     void        cml_adjoin_left(MATRIX *dest, MATRIX *src);
CML_API     void        cml_adjoin_right(MATRIX *dest, MATRIX *src);


/*    ARITHMETIC    */
CML_API     void        cml_add_const(MATRIX *m, double data, MATRIX *opt);
CML_API     void        cml_mul_const(MATRIX *m, double data, MATRIX *opt);
CML_API     void        cml_add(MATRIX *m1, MATRIX *m2, MATRIX *opt);
CML_API     void        cml_sub(MATRIX *m1, MATRIX *m2, MATRIX *opt);
CML_API     void        cml_mul(MATRIX *m1, MATRIX *m2, MATRIX *opt);
CML_API     void        cml_mul_elem(MATRIX *m1, MATRIX *m2, MATRIX *opt);
CML_API     void        cml_div_elem(MATRIX *m1, MATRIX *m2, MATRIX *opt);


/*    OTHER OPERATIONS    */
CML_API     double      cml_min(MATRIX *m);
CML_API     double      cml_max(MATRIX *m);
CML_API     bool        cml_is_zero(MATRIX *m);
CML_API     bool        cml_is_pos(MATRIX *m);
CML_API     bool        cml_is_neg(MATRIX *m);
CML_API     bool        cml_is_nonneg(MATRIX *m);
CML_API     bool        cml_is_equal(MATRIX *m1, MATRIX *m2);
CML_API     void        cml_transpose(MATRIX *m, MATRIX *opt);
CML_API     void        cml_normalize(MATRIX *m, MATRIX *opt);


/*    OTHER OPERATIONS (BLAS/LAPACK REQUIRED)    */
#ifndef CML_NO_DEPENDENCIES
CML_API     bool        cml_inverse(MATRIX *m, MATRIX *opt);
CML_API     bool        cml_sys_equ(MATRIX *A, MATRIX *B, MATRIX *opt);
#endif


#ifdef __cplusplus
}
#endif


/* ========================================================
 *
 *                      IMPLEMENTATION
 *
 * ======================================================== */
#ifdef CML_IMPLEMENTATION


#ifdef __cplusplus
extern "C" {
#endif


#ifdef CML_PRIVATE
    #define VEC_PRIVATE
#endif


#include <stdlib.h>
#include <errno.h>
#include "vec.h" 


#ifndef CML_NO_DEPENDENCIES 
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
void dgetrf_(int*, int*, double*, int*, int*, int*);
void dgetri_(int*, double*, int*, int*, double*, int*, int*);
void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
#endif


#ifdef __cplusplus
}
#endif


CML_API MATRIX* cml_new(size_t rows, size_t cols) {
    if ((rows == 0 && cols > 0) || (cols == 0 && rows > 0)) {
        errno = EINVAL;
        return NULL;
    }

    int old_errno = errno;
    errno = 0;
    MATRIX *m = (MATRIX *) malloc(sizeof(MATRIX));

    if (m == NULL) {
        return NULL;
    }

    m->rows = rows;
    m->cols = cols;
    m->data = (double *) vec_new_cap(sizeof(double), rows * cols);

    if (m->data == NULL) {
        free(m);
        return NULL;
    } else {
        errno = old_errno;
    }

    for (size_t i = 0; i < rows * cols; ++i) {
        vec_push(m->data, 0);
    }

    return m;
}


CML_API MATRIX* cml_dup(MATRIX *m) {
    if (m == NULL) {
        errno = EINVAL;
        return NULL;
    }

    MATRIX *dup = cml_new(m->rows, m->cols);

    if (dup == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < dup->rows; ++i) {
        for (size_t j = 0; j < dup->cols; ++j) {
            cml_set(dup, i, j, cml_get(m, i, j));
        }
    }

    return dup;
}


CML_API MATRIX* cml_ones(size_t rows, size_t cols) {
    MATRIX *m = cml_new(rows, cols);

    if (m == NULL) {
        return NULL;
    }

    cml_set_all(m, 1);

    return m;
}


CML_API MATRIX* cml_identity(size_t dim) {
    MATRIX *m = cml_new(dim, dim);

    if (m == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < dim; ++i) {
        cml_set(m, i, i, 1);
    }

    return m;
}


CML_API MATRIX* cml_lower_tri(size_t dim) {
    MATRIX *m = cml_new(dim, dim);

    if (m == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            cml_set(m, i, j, 1);
        }
    }

    return m;
}


CML_API MATRIX* cml_upper_tri(size_t dim) {
    MATRIX *m = cml_new(dim, dim);

    if (m == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = i; j < dim; ++j) {
            cml_set(m, i, j, 1);
        }
    }

    return m;
}


CML_API void cml_free(MATRIX *m) {
    if (m == NULL) {
        return;
    }

    vec_free(m->data);
    free(m);
}


CML_API double cml_get(MATRIX *m, size_t row, size_t col) {
    if (m == NULL || row >= m->rows || col >= m->cols) {
        errno = EINVAL;
        return 0;
    }

    return *(m->data + m->rows * col + row);
}


CML_API void cml_set(MATRIX *m, size_t row, size_t col, double data) {
    if (m == NULL || row >= m->rows || col >= m->cols) {
        errno = EINVAL;
        return;
    }

    *(m->data + m->rows * col + row) = data;
}


CML_API void cml_set_all(MATRIX *m, double data) {
    if (m == NULL) {
        errno = EINVAL;
        return;
    }

    for (size_t i = 0; i < vec_len(m->data); ++i) {
        m->data[i] = data;
    }
}


CML_API void cml_set_row(MATRIX *m, size_t row, double data) {
    if (m == NULL || row >= m->rows) {
        errno = EINVAL;
        return;
    }

    for (size_t j = 0; j < m->cols; ++j) {
        cml_set(m, row, j, data);
    }
}


CML_API void cml_set_col(MATRIX *m, size_t col, double data) {
    if (m == NULL || col >= m->cols) {
        errno = EINVAL;
        return;
    }

    for (size_t i = 0; i < m->rows; ++i) {
        cml_set(m, i, col, data);
    }
}


CML_API void cml_cpy(MATRIX *dest, MATRIX *src) {
    if (dest == NULL || src == NULL || dest == src) {
        errno = EINVAL;
        return;
    }

    int old_errno = errno;
    errno = 0;
    double *copied_data = (double *) vec_dup(src->data);

    if (errno) {
        return;
    } else {
        errno = old_errno;
    }

    vec_free(dest->data);
    dest->rows = src->rows;
    dest->cols = src->cols;
    dest->data = copied_data;
}


CML_API void cml_cpy_row(MATRIX *dest, MATRIX *src, size_t dest_row, size_t src_row) {
    if (dest == NULL || src == NULL || dest == src || dest->cols != src->cols || dest_row >= dest->rows || src_row >= src->rows) {
        errno = EINVAL;
        return;
    }

    for (size_t j = 0; j < dest->cols; ++j) {
        cml_set(dest, dest_row, j, cml_get(src, src_row, j));
    }
}


CML_API void cml_cpy_col(MATRIX *dest, MATRIX *src, size_t dest_col, size_t src_col) {
    if (dest == NULL || src == NULL || dest == src || dest->rows != src->rows || dest_col >= dest->cols || src_col >= src->cols) {
        errno = EINVAL;
        return;
    }

    for (size_t i = 0; i < dest->rows; ++i) {
        cml_set(dest, i, dest_col, cml_get(src, i, src_col));
    }
}


CML_API void cml_cpy_elem(MATRIX *dest, MATRIX *src, size_t dest_row, size_t dest_col, size_t src_row, size_t src_col) {
    if (dest == NULL || src == NULL || dest == src || dest_row >= dest->rows || dest_col >= dest->cols || src_row >= src->rows || src_col >= src->cols) {
        errno = EINVAL;
        return;
    }

    cml_set(dest, dest_row, dest_col, cml_get(src, src_row, src_col));
}


CML_API void cml_cpy_self_row(MATRIX *m, size_t dest_row, size_t src_row) {
    if (m == NULL || dest_row >= m->rows || src_row >= m->rows) {
        errno = EINVAL;
        return;
    }

    for (size_t j = 0; j < m->cols; ++j) {
        cml_set(m, dest_row, j, cml_get(m, src_row, j));
    }
}


CML_API void cml_cpy_self_col(MATRIX *m, size_t dest_col, size_t src_col) {
    if (m == NULL || dest_col >= m->cols || src_col >= m->cols) {
        errno = EINVAL;
        return;
    }

    for (size_t i = 0; i < m->rows; ++i) {
        cml_set(m, i, dest_col, cml_get(m, i, src_col));
    }
}


CML_API void cml_cpy_self_elem(MATRIX *m, size_t dest_row, size_t dest_col, size_t src_row, size_t src_col) {
    if (m == NULL || dest_row >= m->rows || dest_col >= m->cols || src_row >= m->rows || src_col >= m->cols) {
        errno = EINVAL;
        return;
    }

    cml_set(m, dest_row, dest_col, cml_get(m, src_row, src_col));
}


CML_API void cml_swap(MATRIX *m1, MATRIX *m2) {
    if (m1 == NULL || m2 == NULL || m1 == m2) {
        errno = EINVAL;
        return;
    }

    MATRIX temp = { m1->rows, m1->cols, m1->data };

    m1->rows = m2->rows;
    m1->cols = m2->cols;
    m1->data = m2->data;

    m2->rows = temp.rows;
    m2->cols = temp.cols;
    m2->data = temp.data;
}


CML_API void cml_swap_row(MATRIX *m1, MATRIX *m2, size_t m1_row, size_t m2_row) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1->cols != m2->cols || m1_row >= m1->rows || m2_row >= m2->rows) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(1, m1->cols);

    if (temp == NULL) {
        return;
    }

    cml_cpy_row(temp, m1, 0, m1_row);
    cml_cpy_row(m1, m2, m1_row, m2_row);
    cml_cpy_row(m2, temp, m2_row, 0);
    cml_free(temp);
}


CML_API void cml_swap_col(MATRIX *m1, MATRIX *m2, size_t m1_col, size_t m2_col) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1->rows != m2->rows || m1_col >= m1->cols || m2_col >= m2->cols) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(m1->rows, 1);

    if (temp == NULL) {
        return;
    }

    cml_cpy_col(temp, m1, 0, m1_col);
    cml_cpy_col(m1, m2, m1_col, m2_col);
    cml_cpy_col(m2, temp, m2_col, 0);
    cml_free(temp);
}


CML_API void cml_swap_elem(MATRIX *m1, MATRIX *m2, size_t m1_row, size_t m1_col, size_t m2_row, size_t m2_col) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1_row >= m1->rows || m1_col >= m1->cols || m2_row >= m2->rows || m2_col >= m2->cols) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(1, 1);

    if (temp == NULL) {
        return;
    }

    cml_cpy_elem(temp, m1, 0, 0, m1_row, m1_col);
    cml_cpy_elem(m1, m2, m1_row, m1_col, m2_row, m2_col);
    cml_cpy_elem(m2, temp, m2_row, m2_col, 0, 0);
    cml_free(temp);
}


CML_API void cml_swap_self_row(MATRIX *m, size_t row_1, size_t row_2) {
    if (m == NULL || row_1 >= m->rows || row_2 >= m->rows) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(1, m->cols);

    if (temp == NULL) {
        return;
    }

    cml_cpy_row(temp, m, 0, row_1);
    cml_cpy_self_row(m, row_1, row_2);
    cml_cpy_row(m, temp, row_2, 0);
    cml_free(temp);
}


CML_API void cml_swap_self_col(MATRIX *m, size_t col_1, size_t col_2) {
    if (m == NULL || col_1 >= m->cols || col_2 >= m->cols) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(m->rows, 1);

    if (temp == NULL) {
        return;
    }

    cml_cpy_col(temp, m, 0, col_1);
    cml_cpy_self_col(m, col_1, col_2);
    cml_cpy_col(m, temp, col_2, 0);
    cml_free(temp);
}


CML_API void cml_swap_self_elem(MATRIX *m, size_t row_1, size_t col_1, size_t row_2, size_t col_2) {
    if (m == NULL || row_1 >= m->rows || col_1 >= m->cols || row_2 >= m->rows || col_2 >= m->cols) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(1, 1);

    if (temp == NULL) {
        return;
    }

    cml_cpy_elem(temp, m, 0, 0, row_1, col_1);
    cml_cpy_self_elem(m, row_1, col_1, row_2, col_2);
    cml_cpy_elem(m, temp, row_2, col_2, 0, 0);
    cml_free(temp);
}


CML_API void cml_ins_row(MATRIX *dest, MATRIX *src, size_t dest_row, size_t src_row) {
    if (dest == NULL || src == NULL || dest == src || dest->cols != src->cols || dest_row > dest->rows || src_row >= src->rows) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dest->data, vec_len(dest->data) + dest->cols);

    if (vec_cap(dest->data) < vec_len(dest->data) + dest->cols) {
        errno = ENOMEM;
        return;
    }

    for (size_t j = dest->cols; j-- > 0; ) {
        vec_insert(dest->data, cml_get(src, src_row, j), dest->rows * j + dest_row);
    }
    dest->rows += 1;
}


CML_API void cml_ins_col(MATRIX *dest, MATRIX *src, size_t dest_col, size_t src_col) {
    if (dest == NULL || src == NULL || dest == src || dest->rows != src->rows || dest_col > dest->cols || src_col >= src->cols) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dest->data, vec_len(dest->data) + dest->rows);
    
    if (vec_cap(dest->data) < vec_len(dest->data) + dest->rows) {
        errno = ENOMEM;
        return;
    }

    for (size_t i = dest->rows; i-- > 0; ) {
        vec_insert(dest->data, cml_get(src, i, src_col), dest->rows * dest_col);
    }
    dest->cols += 1;
}


CML_API void cml_ins_self_row(MATRIX *m, size_t dest_row, size_t src_row) {
    if (m == NULL || dest_row > m->rows || src_row >= m->rows) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(1, m->cols);

    if (temp == NULL) {
        return;
    }

    cml_cpy_row(temp, m, 0, src_row);
    cml_ins_row(m, temp, dest_row, 0);
    cml_free(temp);
}


CML_API void cml_ins_self_col(MATRIX *m, size_t dest_col, size_t src_col) {
    if (m == NULL || dest_col > m->cols || src_col >= m->cols) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(m->rows, 1);

    if (temp == NULL) {
        return;
    }

    cml_cpy_col(temp, m, 0, src_col);
    cml_ins_col(m, temp, dest_col, 0);
    cml_free(temp);
}


CML_API void cml_del_all(MATRIX *m) {
    if (m == NULL) {
        errno = EINVAL;
        return;
    }

    m->rows = 0;
    m->cols = 0;
    vec_clear(m->data);
    vec_shrink(m->data);
}


CML_API void cml_del_row(MATRIX *m, size_t row) {
    if (m == NULL || row >= m->rows) {
        errno = EINVAL;
        return;
    }

    for (size_t j = m->cols; j-- > 0; ) {
        vec_remove(m->data, m->rows * j + row);
    }
    m->rows -= 1;
    vec_shrink(m->data);
}


CML_API void cml_del_col(MATRIX *m, size_t col) {
    if (m == NULL || col >= m->cols) {
        errno = EINVAL;
        return;
    }

    for (size_t i = m->rows; i-- > 0; ) {
        vec_remove(m->data, m->rows * col);
    }
    m->cols -= 1;
    vec_shrink(m->data);
}


CML_API void cml_adjoin_top(MATRIX *dest, MATRIX *src) {
    if (dest == NULL || src == NULL || dest == src || dest->cols != src->cols) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dest->data, vec_len(dest->data) + vec_len(src->data));

    if (vec_cap(dest->data) < vec_len(dest->data) + vec_len(src->data)) {
        errno = ENOMEM;
        return;
    }
    
    for (size_t i = src->rows; i-- > 0; ) {
        cml_ins_row(dest, src, 0, i);
    }
}


CML_API void cml_adjoin_bottom(MATRIX *dest, MATRIX *src) {
    if (dest == NULL || src == NULL || dest == src || dest->cols != src->cols) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dest->data, vec_len(dest->data) + vec_len(src->data));

    if (vec_cap(dest->data) < vec_len(dest->data) + vec_len(src->data)) {
        errno = ENOMEM;
        return;
    }

    for (size_t i = 0; i < src->rows; ++i) {
        cml_ins_row(dest, src, dest->rows, i);
    }
}


CML_API void cml_adjoin_left(MATRIX *dest, MATRIX *src) {
    if (dest == NULL || src == NULL || dest == src || dest->rows != src->rows) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dest->data, vec_len(dest->data) + vec_len(src->data));

    if (vec_cap(dest->data) < vec_len(dest->data) + vec_len(src->data)) {
        errno = ENOMEM;
        return;
    }

    for (size_t j = src->cols; j-- > 0; ) {
        cml_ins_col(dest, src, 0, j);
    }
}


CML_API void cml_adjoin_right(MATRIX *dest, MATRIX *src) {
    if (dest == NULL || src == NULL || dest == src || dest->rows != src->rows) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dest->data, vec_len(dest->data) + vec_len(src->data));

    if (vec_cap(dest->data) < vec_len(dest->data) + vec_len(src->data)) {
        errno = ENOMEM;
        return;
    }
    
    for (size_t j = 0; j < src->cols; ++j) {
        cml_ins_col(dest, src, dest->cols, j);
    }
}


CML_API void cml_add_const(MATRIX *m, double data, MATRIX *opt) {
    if (m == NULL || m == opt || m->rows == 0) {
        errno = EINVAL;
        return;
    }

    if (opt != NULL) {
        cml_cpy(opt, m);
        m = opt;
    } 
    for (size_t i = 0; i < vec_len(m->data); ++i) {
        m->data[i] += data;
    }
}


CML_API void cml_mul_const(MATRIX *m, double data, MATRIX *opt) {
    if (m == NULL || m == opt || m->rows == 0) {
        errno = EINVAL;
        return;
    }

    if (opt != NULL) {
        cml_cpy(opt, m);
        m = opt;
    } 
    for (size_t i = 0; i < vec_len(m->data); ++i) {
        m->data[i] *= data;
    }
}


CML_API void cml_add(MATRIX *m1, MATRIX *m2, MATRIX *opt) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1 == opt || m2 == opt || m1->rows != m2->rows || m1->cols != m2->cols || m1->rows == 0) {
        errno = EINVAL;
        return;
    }

    if (opt != NULL) {
        cml_cpy(opt, m1);
        m1 = opt;
    }
    for (size_t i = 0; i < vec_len(m1->data); ++i) {
        m1->data[i] += m2->data[i];
    }
}


CML_API void cml_sub(MATRIX *m1, MATRIX *m2, MATRIX *opt) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1 == opt || m2 == opt || m1->rows != m2->rows || m1->cols != m2->cols || m1->rows == 0) {
        errno = EINVAL;
        return;
    }
    
    if (opt != NULL) {
        cml_cpy(opt, m1);
        m1 = opt;
    }
    for (size_t i = 0; i < vec_len(m1->data); ++i) {
        m1->data[i] -= m2->data[i];
    }
}


CML_API void cml_mul(MATRIX *m1, MATRIX *m2, MATRIX *opt) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1 == opt || m2 == opt || m1->cols != m2->rows || m1->rows == 0) {
        errno = EINVAL;
        return;
    }

    MATRIX *result = cml_new(m1->rows, m2->cols);

    if (result == NULL) {
        return;
    } 

    #ifndef CML_NO_DEPENDENCIES
    char notrans[2] = "N";
    int m = (int) m1->rows;
    int n = (int) m2->cols;
    int k = (int) m1->cols;
    double alpha = 1;
    double beta = 0;

    dgemm_(notrans, notrans, &m, &n, &k, &alpha, m1->data, &m, m2->data, &k, &beta, result->data, &m);
    #else
    for (size_t i = 0; i < m1->rows; ++i) {
        for (size_t j = 0; j < m2->cols; ++j) {
            double sum = 0;
            for (size_t k = 0; k < m1->cols; ++k) {
                sum += cml_get(m1, i, k) * cml_get(m2, k, j);
            }
            cml_set(result, i, j, sum);
        }
    }
    #endif
    
    if (opt != NULL) {
        cml_cpy(opt, result);
    } else {
        cml_cpy(m1, result);
    }
    cml_free(result);
}


CML_API void cml_mul_elem(MATRIX *m1, MATRIX *m2, MATRIX *opt) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1 == opt || m2 == opt || m1->rows != m2->rows || m1->cols != m2->cols || m1->rows == 0) {
        errno = EINVAL;
        return;
    }
    
    if (opt != NULL) {
        cml_cpy(opt, m1);
        m1 = opt;
    }
    for (size_t i = 0; i < vec_len(m1->data); ++i) {
        m1->data[i] *= m2->data[i];
    }
}


CML_API void cml_div_elem(MATRIX *m1, MATRIX *m2, MATRIX *opt) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1 == opt || m2 == opt || m1->rows != m2->rows || m1->cols != m2->cols || m1->rows == 0) {
        errno = EINVAL;
        return;
    }
    
    if (opt != NULL) {
        cml_cpy(opt, m1);
        m1 = opt;
    }
    for (size_t i = 0; i < vec_len(m1->data); ++i) {
        m1->data[i] /= m2->data[i];
    }
}


CML_API double cml_min(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return 0;
    }

    double min = cml_get(m, 0, 0);
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            if (min > cml_get(m, i, j)) {
                min = cml_get(m, i, j);
            }
        }
    }

    return min;
}


CML_API double cml_max(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return 0;
    }

    double max = cml_get(m, 0, 0);
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            if (max < cml_get(m, i, j)) {
                max = cml_get(m, i, j);
            }
        }
    }

    return max;
}


CML_API bool cml_is_zero(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return false;
    }

    for (size_t i = 0; i < vec_len(m->data); ++i) {
        if (m->data[i] != 0) {
            return false;
        }
    }

    return true;
}


CML_API bool cml_is_pos(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return false;
    }

    for (size_t i = 0; i < vec_len(m->data); ++i) {
        if (m->data[i] <= 0) {
            return false;
        }
    }

    return true;
}


CML_API bool cml_is_nonneg(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return false;
    }

    for (size_t i = 0; i < vec_len(m->data); ++i) {
        if (m->data[i] < 0) {
            return false;
        }
    }

    return true;
}


CML_API bool cml_is_equal(MATRIX *m1, MATRIX *m2) {
    if (m1 == NULL || m2 == NULL || m1 == m2) {
        errno = EINVAL;
        return false;
    }

    if (m1->rows != m2->rows || m1->cols != m2->cols) {
        return false;
    }
    for (size_t i = 0; i < vec_len(m1->data); ++i) {
        if (m1->data[i] != m2->data[i]) {
            return false;
        }
    }
    
    return true;
}


CML_API void cml_transpose(MATRIX *m, MATRIX *opt) {
    if (m == NULL || m == opt) {
        errno = EINVAL;
        return;
    }

    if (opt != NULL) {
        cml_cpy(opt, m);
        m = opt;
    }
    MATRIX *trans = cml_new(m->cols, m->rows);
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            cml_set(trans, j, i, cml_get(m, i, j));
        }
    }
    cml_cpy(m, trans);
    cml_free(trans);
}


CML_API void cml_normalize(MATRIX *m, MATRIX *opt) {
    if (m == NULL || m == opt || m->rows == 0) {
        errno = EINVAL;
        return;
    }

    if (opt != NULL) {
        cml_cpy(opt, m);
        m = opt;
    }
    double min = cml_min(m);
    double max = cml_max(m);
    for (size_t i = 0; i < vec_len(m->data); ++i) {
        m->data[i] = (m->data[i] - min) / (max - min);
    }
}


#ifndef CML_NO_DEPENDENCIES


/* Returns false if matrix m doesn't have inverse or error occurred. Otherwise returns true.
 */
CML_API bool cml_inverse(MATRIX *m, MATRIX *opt) {
    if (m == NULL || m == opt || m->rows != m->cols || m->rows == 0) {
        errno = EINVAL;
        return false;
    }

    int old_errno = errno;
    errno = 0;
    MATRIX *backup = cml_new(0, 0);
    cml_cpy(backup, m);
    if (opt != NULL) {
        cml_cpy(opt, m);
        m = opt;
    } 

    if (errno) {
        cml_free(backup);
        return false;
    } else {
        errno = old_errno;
    }

    int n = (int) m->rows;
    int ipiv[n];
    int lwork = n * n;
    double work[lwork];
    int info;
    dgetrf_(&n, &n, m->data, &n, ipiv, &info);
    dgetri_(&n, m->data, &n, ipiv, work, &lwork, &info);

    if (info < 0) {
        cml_cpy(m, backup);
        cml_free(backup);
        errno = EDOM;
        return false;
    } else if (info > 0) {
        cml_cpy(m, backup);
        cml_free(backup);
        return false;
    }

    cml_free(backup);

    return true;
}


/* Solves the system of equations 'AX = B' for X and stores it in B (NOT A, since B already has desired dimensions), or
 * 'opt' if 'opt != NULL'. Returns false if unsolvable or error occurred. Otherwise returns true.
 */
CML_API bool cml_sys_equ(MATRIX *A, MATRIX *B, MATRIX *opt) {
    if (A == NULL || B == NULL || A == B || A == opt || B == opt || A->rows != A->cols || A->cols != B->rows || A->rows == 0) {
        errno = EINVAL;
        return false;
    }

    int old_errno = errno;
    errno = 0;
    MATRIX *A_copy = cml_new(0, 0);
    cml_cpy(A_copy, A);
    A = A_copy;
    MATRIX *backup = cml_new(0, 0);
    cml_cpy(backup, B);
    if (opt != NULL) {
        cml_cpy(opt, B);
        B = opt;
    }

    if (errno) {
        cml_free(A);
        cml_free(backup);
        return false;
    } else {
        errno = old_errno;
    }

    int n = (int) A->rows;
    int nrhs = (int) B->cols;
    int ldb = (int) B->rows;
    int ipiv[n];
    int info;
    dgesv_(&n, &nrhs, A->data, &n, ipiv, B->data, &ldb, &info);

    if (info < 0) {
        cml_cpy(B, backup);
        cml_free(A);
        cml_free(backup);
        errno = EDOM;
        return false;
    } else if (info > 0) {
        cml_cpy(B, backup);
        cml_free(A);
        cml_free(backup);
        return false;
    }

    cml_free(A);
    cml_free(backup);

    return true;
}


#endif


#endif


#endif 
