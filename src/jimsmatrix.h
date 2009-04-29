#ifndef _JIMSMATRIX_H
#define _JIMSMATRIX_H

// Added 12/08/08:
//#define _MULTINOMIAL_CHK_
 
//#define _DBG_
//#define _MATRIX_BOUNDS_CHECK_
//#define _EP_DBG_
//#define _STARTING_STATE_DBG_
//#define _DBG_SIGMA_DRAW_ 
//#define _DBG_ETA_DRAW_
//#define_NCHG_BOUNDS_CHECK_
//#define _USE_ADJUST_
//#define _MAT_MULT_DBG_
#define _LAPACK_INVERSE_
#define _LAPACK_DET_
//#define _LAPACK_MULTIPLY_
#define _LAPACK_CHOL_

#define GQ_INLINE static R_INLINE

#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>
#include<R_ext/Lapack.h>

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<float.h>

typedef unsigned index_t;
typedef double Matrix;

double GQ_INLINE
max(const double a, const double b)
{
  return ((a>b)?a:b);
}

double GQ_INLINE
min(const double a, const double b)
{
  return ((a<b)?a:b);
}

/******************************************************************************/
/**********************  MATRIX LIBRARY DECLARATIONS **************************/
/******************************************************************************/

// R-allocation:
Matrix *
rmatrix_new(index_t nrow, index_t ncol);

// C-allocation:
Matrix *
matrix_new(index_t nrow, index_t ncol);

// R-allocation:
Matrix *
rmatrix_diag_new(index_t nrow, index_t ncol, double diag);

// C-allocation:
Matrix *
matrix_diag_new(index_t nrow, index_t ncol, double diag);

// Regular function declarations:

index_t GQ_INLINE numrows(Matrix const * const xx)
{
  return (index_t)xx[-1];
}

index_t GQ_INLINE numcols(Matrix const * const xx)
{
  return (index_t)xx[-2];
}

double GQ_INLINE
matrix_get_element(Matrix const * const xx, index_t row, index_t col)
{
#ifdef _MATRIX_BOUNDS_CHECK_
if ((row<0)||(row>=numrows(xx))||(col<0)||(col>=numcols(xx)))
  error("Invalid call to matrix_get_element()");
#endif
  return xx[row+(index_t)xx[-1]*col];
}

double GQ_INLINE
matrix_fast_get_element(Matrix *xx, index_t row, index_t col, index_t nrows)
{
#ifdef _MATRIX_BOUNDS_CHECK_
if (nrows!=numrows(xx))
  error("Argument nrows in matrix_fast_get_element() does not match actual numrows: arg:%u, act:%u",
	nrows,numrows(xx));
if ((row<0)||(row>=nrows)||(col<0)||(col>=numcols(xx)))
  error("Invalid call to matrix_fast_get_element()\nRequested [%u,%u], Allowed: lo=[%u,%u], hi=[%u,%u]",
	row,col,0,0,numrows(xx)-1,numcols(xx)-1);
#endif
  return xx[row+nrows*col];
}
double GQ_INLINE
matrix_fast_get_element_const(Matrix const * const xx, const index_t row, const index_t col, const index_t nrows)
{
#ifdef _MATRIX_BOUNDS_CHECK_
if (nrows!=numrows(xx))
  error("Argument nrows in matrix_fast_get_element_const() does not match actual numrows: arg:%u, act:%u",
	nrows,numrows(xx));
if ((row<0)||(row>=nrows)||(col<0)||(col>=numcols(xx)))
  error("Invalid call to matrix_fast_get_element_const()\nRequested [%u,%u], Allowed: lo=[%u,%u], hi=[%u,%u]",
	row,col,0,0,numrows(xx)-1,numcols(xx)-1);
#endif
  return xx[row+nrows*col];
}

void GQ_INLINE
matrix_set_element(Matrix *xx, index_t row, index_t col, double setval)
{
#ifdef _MATRIX_BOUNDS_CHECK_
if ((row<0)||(row>=numrows(xx))||(col<0)||(col>=numcols(xx)))
  error("Invalid call to matrix_set_element()");
#endif
  xx[row+(index_t)xx[-1]*col] = setval;
  return;
}

void GQ_INLINE
matrix_fast_set_element(Matrix *xx, index_t row, index_t col, index_t nrows, double setval)
{
#ifdef _MATRIX_BOUNDS_CHECK_
if (nrows!=numrows(xx))
  error("Argument nrows in matrix_fast_set_element() does not match actual numrows: arg:%u, act:%u",
	nrows,numrows(xx));
if ((row<0)||(row>=nrows)||(col<0)||(col>=numcols(xx)))
  error("Invalid call to matrix_fast_set_element()\nRequested [%u,%u], Allowed: lo=[%u,%u], hi=[%u,%u]",
	row,col,0,0,numrows(xx)-1,numcols(xx)-1);
#endif
  xx[row+nrows*col] = setval;
  return;
}

void GQ_INLINE
matrix_increment_element(Matrix *xx, index_t row, index_t col, double byval)
{
#ifdef _MATRIX_BOUNDS_CHECK_
if ((row<0)||(row>=numrows(xx))||(col<0)||(col>=numcols(xx)))
  error("Invalid call to matrix_increment_element()");
#endif
  xx[row+(index_t)xx[-1]*col] += byval;
  return;
}

void GQ_INLINE
matrix_fast_increment_element(Matrix *xx, index_t row, index_t col, index_t nrows, double byval)
{
#ifdef _MATRIX_BOUNDS_CHECK_
if (nrows!=numrows(xx))
  error("Argument nrows in matrix_fast_increment_element() does not match actual numrows: arg:%u, act:%u",nrows,numrows(xx));
if ((row<0)||(row>=nrows)||(col<0)||(col>=numcols(xx)))
  error("Invalid call to matrix_fast_increment_element()");
#endif
  xx[row+nrows*col] += byval;
  return;
}

void GQ_INLINE matrix_free(Matrix *xx)
{ // Free a double matrix allocated by matrix_new()
  if (xx!=NULL){free(xx-2);}
  return;
}

void GQ_INLINE
matrix_fastset(Matrix *xx, int val)
{
  memset(xx,val,sizeof(double)*(numrows(xx)*numcols(xx)));
  return;
}

int
matrix_assert(Matrix *xx);

int
matrix_assert_column_vec(Matrix *xx);

int
matrix_assert_row_vec(Matrix *xx);

int
matrix_assert_vec(Matrix *xx);


void //GQ_INLINE
matrix_get_set_block(Matrix *yy, index_t r1_yy, index_t r2_yy, index_t c1_yy, index_t c2_yy,
		     Matrix *xx, index_t r1_xx, index_t r2_xx, index_t c1_xx, index_t c2_xx);
/*{
  index_t ii_xx, jj_xx, ii_yy, jj_yy;
  index_t r2_xx_p1 = r2_xx+1;
  index_t c2_xx_p1 = c2_xx+1;
  ii_yy = r1_yy;
  jj_yy = c1_yy;

  for (ii_xx=r1_xx; ii_xx<r2_xx_p1; ii_xx++){
    for (jj_xx=c1_xx; jj_xx<c2_xx_p1; jj_xx++){
      matrix_set_element(yy, ii_yy, jj_yy, matrix_get_element(xx, ii_xx, jj_xx));
      jj_yy++;
    }
    ii_yy++;
    jj_yy = c1_yy;
  }
  return;
}*/
/*{
  index_t ii, jj;
  for (ii=r1_yy; ii<=r2_yy; ii++)
    for (jj=c1_yy; ii<=c2_yy; ii++)
      matrix_set_element(yy,ii,jj, matrix_get_element(xx,ii-r1_yy+r1_xx,jj-c1_yy+c1_xx));
  return;
}*/

void
matrix_print_element(Matrix const * const xx, index_t row, index_t col);

void
matrix_print_all(Matrix const * const xx);

void
matrix_print_subset(Matrix const * const xx, index_t irmin, index_t irmax, index_t icmin, index_t icmax);

void
matrix_transpose(Matrix  *xx, Matrix *yy);

void
matrix_transpose_same(Matrix *xx);

double matrix_quadform(Matrix *x, Matrix *A, Matrix *y);

void matrix_multiply(Matrix *A, char tA, Matrix *B, char tB, Matrix *C);

void
matrix_scalar_multiply(Matrix *xx, double ss, Matrix *yy);

void
matrix_add(Matrix *xx, Matrix *yy, Matrix *zz);

//void matrix_get_row(Matrix *m, index_t i, Matrix *v);

void GQ_INLINE
matrix_get_row(Matrix *m, index_t i, Matrix *v)
{
  // Extracts the ith row of m to the column vector v
  index_t j;
  const index_t nc=numcols(m);
#ifdef _MATRIX_BOUNDS_CHECK_
  const index_t len=numrows(v);
  if (len!=nc)
    error("Incompatible dimensions in matrix_get_row()");
#endif
  for (j=0; j<nc; j++)
    matrix_set_element(v,j,0, matrix_get_element(m,i,j));
  return;
}

double GQ_INLINE
matrix_as_scalar(Matrix  *xx)
{
  return matrix_get_element(xx,0,0);
}

int
matrix_is_symmetric(Matrix *xx);

void
matrix_ADJUST(Matrix *xx, index_t kk);

void matrix_identity(Matrix *xx);
/*
void
matrix_inverse(Matrix *xx, Matrix *yy, Matrix *AUG_in_matrix_inverse);*/
void matrix_inverse(Matrix *X, Matrix *X_inverse, Matrix *pointlessArgument);
void lapack_det(Matrix *X, Matrix *v, unsigned int useLog);

void
matrix_DOOLITTLE(Matrix *xx, index_t kk);

/*double
matrix_determinant(Matrix *xx, Matrix *yy);*/
double 
matrix_determinant(Matrix *X, Matrix *deformableArg, unsigned useLog);

//void blas_matrix_cholesky(Matrix  *X, Matrix *Y);

void
matrix_cholesky(Matrix *xx, Matrix *yy);

void
matrix_copy(Matrix *ToBeCopied, Matrix *Copy);

void
matrix_subtract(Matrix *xx, Matrix *yy, Matrix *zz);

void
matrix_sum_xx_m_mu(Matrix *YY, Matrix *XX, Matrix *mu_vec);

void
matrix_sum_xx_m_mu_by_precinct(Matrix * const SS, 
                               Matrix * const OMEGAS, 
                               Matrix * const mu_mat_cu);

#endif /* jimsmatrix.h */
