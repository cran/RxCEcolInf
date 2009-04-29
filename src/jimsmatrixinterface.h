
#ifndef _JIMSMATRIXINTERFACE_H
#define _JIMSMATRIXINTERFACE_H

#include "jimsmatrix.h"
#include "jimsimatrix.h"

/****************************************************************************/
/*****************  C UNPACK AND REPACK FUNCTION DECLARATIONS  *************/
/****************************************************************************/

Matrix *
rmatrix_unpack_new(SEXP mat);

Matrix *
matrix_unpack_new(SEXP mat);

Matrix_int *
rmatrix_int_unpack_new(SEXP mat);

Matrix_int *
matrix_int_unpack_new(SEXP mat);

SEXP
matrix_repack_new(Matrix *xx);

Matrix *
rmatrix_vector_unpack_new(SEXP inVec);

Matrix *
matrix_vector_unpack_new(SEXP inVec);

SEXP
matrix_vector_repack_new(Matrix *inVec);

int
matrix_match_dims(Matrix *xx, Matrix *yy);

#endif /* jimsmatrixinterface.h */
