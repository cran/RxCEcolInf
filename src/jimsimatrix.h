
#ifndef _JIMSIMATRIX_H
#define _JIMSIMATRIX_H

#ifndef _JIMSMATRIX_H
#include "jimsmatrix.h"
#endif

typedef int Matrix_int;

/****************************************************************/
/*******************  INT MATRIX ALGEBRA LIBRARY ****************/
/****************************************************************/

// R-allocation:
Matrix_int *
rmatrix_new_int(index_t nrow, index_t ncol);

// C-allocation:
Matrix_int *
matrix_new_int(index_t nrow, index_t ncol);

void
matrix_free_int(Matrix_int *xx);

index_t GQ_INLINE
numrows_int(Matrix_int *xx)
{
  return (index_t)xx[-1];
}

index_t GQ_INLINE
numcols_int(Matrix_int *xx)
{
  return (index_t)xx[-2];
}

int GQ_INLINE
matrix_get_int_element(Matrix_int *xx, index_t row, index_t col)
{
  return xx[row+(index_t)xx[-1]*col];
}

int GQ_INLINE
matrix_fast_get_int_element(Matrix_int *xx, index_t row, index_t col, index_t nrows)
{
  return xx[row+nrows*col];
}

void GQ_INLINE
matrix_set_int_element(Matrix_int *xx, index_t row, index_t col, int setval)
{
  xx[row+(index_t)xx[-1]*col] = setval;
  return;
}

void GQ_INLINE
matrix_fast_set_int_element(Matrix_int *xx, index_t row, index_t col, index_t nrows, int setval)
{
  xx[row+nrows*col] = setval;
  return;
}

void
sample_equ_pr_wo_replace_int(Matrix_int *bigvec_int,
			     Matrix_int *smallvec_int);
void
matrix_print_int_element(Matrix_int *xx_int, index_t row, index_t col);

void
matrix_print_all_int(Matrix_int *xx_int);

#endif // jimsimatrix.h
