
#include "jimsmatrixinterface.h"

/**********************************************************************************/
/********************  C UNPACK AND REPACK FUNCTION DEFINITIONS ********************/
/**********************************************************************************/

Matrix *
rmatrix_unpack_new(SEXP mat)
{
  //  Takes an R object of type SEXP matrix and returns a
  //  pointer to a Matrix object
  //  Extract dimensions and create the matrix:
  const index_t nrow = nrows(mat);
  const index_t ncol = ncols(mat);
  Matrix *retMat = rmatrix_new(nrow, ncol);
  double *aa = REAL(mat);
  memcpy(retMat,aa,sizeof(double)*nrow*ncol);
  return retMat;
}

Matrix *
matrix_unpack_new(SEXP mat)
{
  //  Takes an R object of type SEXP matrix and returns a
  //  pointer to a Matrix object
  //  Extract dimensions and create the matrix:
  const index_t nrow = nrows(mat);
  const index_t ncol = ncols(mat);
  Matrix *retMat = matrix_new(nrow, ncol);
  double *aa = REAL(mat);
  memcpy(retMat,aa,sizeof(double)*nrow*ncol);
  return retMat;
}

Matrix_int *
rmatrix_int_unpack_new(SEXP mat)
{
  int *aa = INTEGER(mat);
  //  Extract the dimensions.
  const index_t nrow = nrows(mat);
  const index_t ncol = ncols(mat);
  Matrix_int *retMat = rmatrix_new_int(nrow, ncol);
  index_t ii, jj;
  for (ii=0; ii<nrow; ii++)
    for (jj=0; jj<ncol; jj++)
      matrix_set_int_element(retMat, ii, jj, aa[ii+nrow*jj]);
  return retMat;
}

Matrix_int *
matrix_int_unpack_new(SEXP mat)
{
  int *aa = INTEGER(mat);
  //  Extract the dimensions.
  const index_t nrow = nrows(mat);
  const index_t ncol = ncols(mat);
  Matrix_int *retMat = matrix_new_int(nrow, ncol);
  index_t ii, jj;
  for (ii=0; ii<nrow; ii++)
    for (jj=0; jj<ncol; jj++)
      matrix_set_int_element(retMat, ii, jj, aa[ii+nrow*jj]);
  return retMat;
}


/****************************************************************************/

SEXP
matrix_repack_new(Matrix *xx)
{
  // Packages a matrix into SEXP form to return to R.  The SEXP will 
  // be in matrix form, not in the vec form of the matrix.

  const index_t nrow_xx=numrows(xx);
  const index_t ncol_xx=numcols(xx);
  SEXP retMat;
  PROTECT(retMat = allocMatrix(REALSXP,nrow_xx,ncol_xx));
  double *aa=REAL(retMat);
  memcpy(aa,xx,sizeof(double)*nrow_xx*ncol_xx);
  return retMat;
}


/******************************************************************************/

Matrix *
rmatrix_vector_unpack_new(SEXP inVec)
{
  // Returns a Greiner type matrix with 1 row (i.e., a row vector)
  // from an R SEXP object of type vector.
  const index_t ncol = length(inVec);
  double *aa = REAL(inVec);
  Matrix * const retVec = rmatrix_new(1, ncol);
  memcpy(retVec,aa,sizeof(double)*ncol);
  return retVec;
}

Matrix *
matrix_vector_unpack_new(SEXP inVec)
{
  // Returns a Greiner type matrix with 1 row (i.e., a row vector)
  // from an R SEXP object of type vector.
  const index_t ncol = length(inVec);
  double *aa = REAL(inVec);
  Matrix * const retVec = matrix_new(1, ncol);
  memcpy(retVec,aa,sizeof(double)*ncol);
  return retVec;
}

/******************************************************************************/

SEXP
matrix_vector_repack_new(Matrix *inVec)
{
  //  Takes either a row or a column vector and returns an SEXP vector

  index_t ii;
  const index_t nrow_inVec = numrows(inVec);
  const index_t ncol_inVec = numcols(inVec);

  SEXP retVec;
  PROTECT(retVec = allocVector(REALSXP, max(nrow_inVec,ncol_inVec)));
  double *aa = REAL(retVec);

  if (nrow_inVec == 1){  // if row vector
    for (ii=0; ii<ncol_inVec; ii++)
      aa[ii] = matrix_get_element(inVec, 0, ii);
  }
  else { // if column vector
    if (ncol_inVec!=1)
      error("inVec is not a valid vector in matrix_vector_repack_new()");
    for (ii=0; ii<nrow_inVec; ii++)
      aa[ii] = matrix_get_element(inVec, ii, 0);
  }
  return retVec;
}


int
matrix_match_dims(Matrix *xx, Matrix *yy)
{
  return (((numrows(xx) != numrows(yy)) || (numcols(xx) != numcols(yy))) ? 0 : 1);
}
