#ifndef _JIMSMATRIX_P_H_
#define _JIMSMATRIX_P_H_

#include "jimsmatrix.h"
#include "jimsmatrixinterface.h"

typedef Matrix* Matrix_p;

// Constructors:
Matrix_p *rmatrix_p_new(index_t len);
Matrix_p *matrix_p_new(index_t len);

// Destructor:
void matrix_p_free(Matrix_p *mat_p, index_t len);

// Selectors and mutators:
void set_mat_p_ptr(Matrix_p *mat_p, index_t i, Matrix *mat);
Matrix *get_mat_p_ptr(Matrix_p *mat_p, index_t i);

// Printing tools:
void matrix_p_print_all(Matrix_p *mat_p, index_t len);
void matrix_p_print_subset(Matrix_p *mat_p, index_t lo, index_t hi, index_t rl,
                           index_t rh, index_t cl, index_t ch);

// Unpacking tools:
Matrix_p *matrix_p_unpack_new(SEXP rLst);
Matrix_p *rmatrix_p_unpack_new(SEXP rLst);

#endif
