#include "jimsmatrix.h"

/****************************************************************/
/*******************  MATRIX ALGEBRA LIBRARY ********************/
/****************************************************************/

Matrix *
rmatrix_new(index_t nrow, index_t ncol)
{
  // R-allocation (has its own error check):
  double *mat = (double *) R_alloc((nrow*ncol)+2,sizeof(double));

  // Initialize to 0:
  memset(mat,0,sizeof(double)*((nrow*ncol)+2));

  // Return pointer to array of pointers to rows
  mat[0] = (double)ncol;
  mat[1] = (double)nrow;
  return mat+2;
}

Matrix *
matrix_new(index_t nrow, index_t ncol)
{
  // C-allocation:
  double *mat = (double *) malloc((size_t)((nrow*ncol)+2)*sizeof(double));
  if (mat==NULL) 
    error("allocation failure in matrix_new()");

  // Initialize to 0:
  memset(mat,0,sizeof(double)*((nrow*ncol)+2));

  // Return pointer to array of pointers to rows
  mat[0] = (double)ncol;
  mat[1] = (double)nrow;
  return mat+2;
}

/*****************************************************************************/

Matrix *
rmatrix_diag_new(index_t nrow, index_t ncol, double diag)
{
  // R-allocation:
  Matrix *mat = rmatrix_new(nrow,ncol);
 
  // Initialize diagonal elements to diag:
  index_t ii;
  for (ii=0; ii<ncol; ii++)
    mat[ii+nrow*ii] = diag;
  return mat;
}

Matrix *
matrix_diag_new(index_t nrow, index_t ncol, double diag)
{
  // C-allocation:
  Matrix *mat = matrix_new(nrow,ncol);
 
  // Initialize diagonal elements to diag:
  index_t ii;
  for (ii=0; ii<ncol; ii++)
    mat[ii+nrow*ii] = diag;
  return mat;
}

/*****************************************************************************/

int
matrix_assert(Matrix *xx)
{
  //  Runs several standard checks on a matrix. Returns 0
  //  if anything is amiss, 1 if all is well 
  int ret_val = 1;

  if (xx == NULL){
    Rprintf("Error:  Matrix that should not be NULL is NULL.\n");
    ret_val = 0;
  }
  if ((numrows(xx) <= 0) || (numcols(xx) <= 0)){
    Rprintf("Error:  Matrix has fewer than 1 row or fewer than 1 column.\n");
    ret_val = 0;
  }
  return ret_val;
}


/*****************************************************************************/

int
matrix_assert_column_vec(Matrix *xx)
{
  //  Runs several standard checks on a matrix that is supposed
  //  to be a column vector. Returns 0
  //  if anything is amiss, 1 if all is well.
  int ret_val = 1;

  if (xx == NULL)
    error("Error:  Column vector that should not be NULL is NULL.\n");
    
  if (numrows(xx) <= 0)
    error("Error:  Column vector has fewer than 1 row.\n");
  
  if (numcols(xx) != 1)
    error("Error:  Column vector has number of columns not equal to 1.\n");
  
  return ret_val;
}

/*****************************************************************************/
/*  We use the R notation that the call matrix_get_block(Mat1, 0, 3, 3, 5, Mat2, etc.)
    requires Mat2 to have space for FOUR rows and THREE columns, consisting
    of rows 0:3 and of columns 3:5 of Mat 1.  */
void
matrix_get_set_block(Matrix *yy, index_t r1_yy, index_t r2_yy, index_t c1_yy, index_t c2_yy,
		     Matrix *xx, index_t r1_xx, index_t r2_xx, index_t c1_xx, index_t c2_xx)
{
  index_t ii_xx, jj_xx, ii_yy, jj_yy;
  const index_t r2_xx_p1 = r2_xx+1;
  const index_t c2_xx_p1 = c2_xx+1;
  const index_t nr_yy = numrows(yy);
  const index_t nr_xx = numrows(xx);
  ii_yy = r1_yy;
  jj_yy = c1_yy;

  for (ii_xx=r1_xx; ii_xx<r2_xx_p1; ii_xx++){
    for (jj_xx=c1_xx; jj_xx<c2_xx_p1; jj_xx++){
      matrix_fast_set_element(yy, ii_yy, jj_yy, nr_yy, matrix_fast_get_element(xx, ii_xx, jj_xx, nr_xx));
      jj_yy++;
    }
    ii_yy++;
    jj_yy = c1_yy;
  }
  return;
}

/*****************************************************************************/

void
matrix_print_element(Matrix const * const xx, index_t row, index_t col)
{
  Rprintf("%g\t", matrix_get_element(xx, row, col));
  return;
}

/*****************************************************************************/

void
matrix_print_all(Matrix const * const xx)
{
  index_t ii, jj, nrow_xx = numrows(xx), ncol_xx = numcols(xx);
  for (ii=0; ii<nrow_xx; ii++){
    for (jj=0; jj<ncol_xx; jj++){
      matrix_print_element(xx, ii, jj);
    }
    Rprintf("\n");
  }
  return;
}

void
matrix_print_subset(Matrix const * const xx, index_t irmin, index_t irmax, index_t icmin, index_t icmax)
{
  index_t ii, jj;
  for (ii=irmin; ii<=irmax; ii++){
    for (jj=icmin; jj<=icmax; jj++){
      matrix_print_element(xx, ii, jj);
    }
    Rprintf("\n");
  }
  return;
}

/*****************************************************************************/

void
matrix_transpose(Matrix  *xx, Matrix *yy)
{ //  yy <- t(xx)
  index_t ii, jj;
  const index_t nrow_xx=numrows(xx);
  const index_t nrow_yy=numrows(yy);
  const index_t ncol_xx=numcols(xx);
  for (ii=0; ii<nrow_xx; ii++)
    for (jj=0; jj<ncol_xx; jj++)
      matrix_fast_set_element(yy, jj,ii,nrow_yy, matrix_fast_get_element(xx, ii,jj,nrow_xx));
  return;
}

/*****************************************************************************/
double matrix_quadform(Matrix *x, Matrix *A, Matrix *y)
{
  // Computes x^{T}Ay
  index_t i,j, nrowy=numrows(y), nrowx=numrows(x);
  double ret=0.0;
  if ((nrowy!=numcols(A)) || (nrowx!=numrows(A)))
    error("Incompatible dims in matrix_quadform()");
  for (i=0; i<nrowx; i++)
    for (j=0; j<nrowy; j++)
      ret += (matrix_get_element(x,i,0)*matrix_get_element(A,i,j)*matrix_get_element(y,j,0));
  return ret;
}

void
matrix_transpose_same(Matrix *xx)
{
  //  In R notation, xx<-t(xx), but ONLY FOR SQUARE xx.
  //  Second function written for speed.

  double aa;
  index_t ii, jj, nrow_xx = numrows(xx), ncol_xx = numcols(xx);

  for (ii=0; ii<nrow_xx; ii++){
    for (jj=(ii+1); jj<ncol_xx; jj++){
      aa = matrix_get_element(xx, ii, jj);
      matrix_set_element(xx, ii, jj, matrix_get_element(xx, jj, ii));
      matrix_set_element(xx, jj, ii, aa);
    }
  }
  return;
}

/*****************************************************************************/
#ifndef _LAPACK_MULTIPLY_
void
matrix_multiply(Matrix * xx, 
		char tX, 
		Matrix * yy, 
		char tY, 
		Matrix * zz)
{
  //  In R notation, this function does zz <- xx %*% yy.
  //  User is responsible for allocating memory for zz.

  index_t mm, nn, pp;
  double aa;
  const index_t nrow_xx = numrows(xx);
  const index_t ncol_xx = numcols(xx);
  const index_t ncol_yy = numcols(yy);
  const index_t nrow_yy = numrows(yy);
  const index_t nrow_zz = numrows(zz);

  if ((tX=='N')&&(tY=='N')){

    for (mm=0; mm<nrow_xx; mm++){
      for (pp=0; pp<ncol_yy; pp++){
        aa = 0.0;
        for (nn=0; nn<ncol_xx; nn++){
	  aa += matrix_fast_get_element(xx,mm,nn,nrow_xx)*matrix_fast_get_element(yy,nn,pp,nrow_yy);
        }
        matrix_fast_set_element(zz, mm,pp,nrow_zz, aa);
      }
    }
    return;

  } else if ((tX=='T')&&(tY=='N')){

    for (mm=0; mm<ncol_xx; mm++){
      for (pp=0; pp<ncol_yy; pp++){
        aa = 0.0;
        for (nn=0; nn<nrow_xx; nn++){
	  aa += matrix_fast_get_element(xx,nn,mm,nrow_xx)*matrix_fast_get_element(yy,nn,pp,nrow_yy);
        }
        matrix_fast_set_element(zz, mm,pp,nrow_zz, aa);
      }
    }
    return;

  } else if ((tX=='N')&&(tY=='T')){
    
    for (mm=0; mm<nrow_xx; mm++){
      for (pp=0; pp<nrow_yy; pp++){
        aa = 0.0;
        for (nn=0; nn<ncol_xx; nn++){
	  aa += matrix_fast_get_element(xx,mm,nn,nrow_xx)*matrix_fast_get_element(yy,pp,nn,nrow_yy);
        }
        matrix_fast_set_element(zz, mm,pp,nrow_zz, aa);
      }
    }
    return;

  } else if ((tX=='T')&&(tY=='T')){

    for (mm=0; mm<ncol_xx; mm++){
      for (pp=0; pp<nrow_yy; pp++){
        aa = 0.0;
        for (nn=0; nn<nrow_xx; nn++){
	  aa += matrix_fast_get_element(xx,nn,mm,nrow_xx)*matrix_fast_get_element(yy,pp,nn,nrow_yy);
        }
        matrix_fast_set_element(zz, mm,pp,nrow_zz, aa);
      }
    }
    return;

  } else 
   error("Invalid tX and tY arguments in matrix multiply");

/*
    Matrix *xx_copy, *yy_copy;
    if (tX=='T'){
      nrow_xx=numcols(xx), ncol_xx=numrows(xx);
      xx_copy = matrix_new(nrow_xx,ncol_xx);
      matrix_transpose(xx,xx_copy);
    } else {
      xx_copy = xx;
    }
    if (tY=='T'){
      nrow_yy=numcols(yy), ncol_yy=numrows(yy);
      yy_copy = matrix_new(nrow_yy,ncol_yy);
      matrix_transpose(yy,yy_copy);
    } else {
      yy_copy = yy;
    }
    const index_t nrow_xxcpy = numrows(xx_copy);
    const index_t nrow_yycpy = numrows(yy_copy);

    for (mm=0; mm<nrow_xx; mm++){
      for (pp=0; pp<ncol_yy; pp++){
        aa = 0.0;
        for (nn=0; nn<ncol_xx; nn++){
	  aa += matrix_fast_get_element(xx_copy, mm,nn,nrow_xxcpy)*
		matrix_fast_get_element(yy_copy, nn,pp,nrow_yycpy);
        }
        matrix_fast_set_element(zz, mm,pp,nrow_zz, aa);
      }
    }
    if (tX=='T')
      matrix_free(xx_copy);
    if (tY=='T')
      matrix_free(yy_copy);
  } 
*/ 
  return;
}
#endif
#ifdef _LAPACK_MULTIPLY_
void matrix_multiply(Matrix *A, char tA, Matrix *B, char tB, Matrix *C)
{
  char *transa, *transb;
  int lda = numrows(A);
  int ldb = numrows(B);
  int ldc = numrows(C);
  int m; // numrows(TRANSA(A)) == numrows(C)
  int n; // numcols(TRANSB(B)) == numcols(C)
  int k; // numcols(TRANSA(A)) == numrows(TRANSB(B))
  int kChk;

  if (tA=='T'){
    transa="T";
    m = numcols(A);
    k = numrows(A);
  } else {
    transa="N";
    m = numrows(A);
    k = numcols(A);
  }

  if (tB=='T'){
    transb="T";
    n    = numrows(B);
    kChk = numcols(B);
  } else {
    transb = "N";
    n    = numcols(B);
    kChk = numrows(B);
  }

  // Checks:
  if (k!=kChk)
    error("Incompatible dimensions in matrix_multiply()");
  if ( (m!=numrows(C)) || (n!=numcols(C)) )
    error("Incompatible output matrix in matrix_multiply()");

  double one=1.0, zero=0.0;

  // Compute: C := alpha*TRANSA(A)%*%TRANSB(B) + beta*C
  //
  //SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  //
  // M = numrows(TRANSA(A)) == numrows(C)
  // N = numcols(TRANSB(B)) == numcols(C)
  // K = numcols(TRANSA(A)) == numcols(TRANSB(B))
  // ALPHA = Scalar multiple of TRANSA(A)%*%TRANSB(B)
  //
  // A = Matrix A (LDA x ka) where 'N' => ka=K, 'T' => ka=M
  // LDA = First dim of TRANS(A) (i.e. 'N'>=max(1,M), 'T'>=max(1,K))
  // [NOTE: Only the first M (if 'N') or K (if 'T') rows are used]
  //
  // B = Matrix B (LDB x kb) where 'N' => kb=N, 'T' => kb=K
  // LDB = First dim of TRANS(A) (i.e. 'N'>=max(1,K), 'T'>=max(1,N))
  // [NOTE: Only the first N (if 'N') or K (if 'T') rows are used]
  //
  // BETA = Scalar multiple of C on RHS
  // C = Matrix C (output)
  // LDC = First dim of C (i.e. LDC>=max(1,M))
  //
  // C is the only input that is changed on exit
  
  if (m>0 && n>0 && k>0){
    F77_CALL(dgemm)(transa,transb,&m,&n,&k,&one,A,&lda,B,&ldb,&zero,C,&ldc FCONE FCONE);
  } else
    matrix_fastset(C,0);
  return;
}
#endif
/*****************************************************************************/

void
matrix_scalar_multiply(Matrix *xx, double ss, Matrix *yy)
{
  //  In R notation, executes yy <- ss * xx.  User responsible
  //  for allocating memory to yy.  yy can be the same as xx.

  const index_t nrtnc=numrows(xx)*numcols(xx);
  matrix_copy(xx,yy);
  index_t ii;
  for (ii=0; ii<nrtnc; ii++)
    yy[ii] *= ss;
  return;
}

/*****************************************************************************/

void
matrix_add(Matrix *xx, Matrix *yy, Matrix *zz)
{
  //  In R notation, implements zz <- xx + yy.  User responsible
  //  for allocating memory to zz.  zz could be the same as xx
  //  or yy.
  index_t ii, jj;
  const index_t nrow_xx = numrows(xx);
  const index_t ncol_xx = numcols(xx);
  const index_t nrow_yy = numrows(yy);
  const index_t nrow_zz = numrows(zz);
  for (ii=0; ii<nrow_xx; ii++)
    for (jj=0; jj<ncol_xx; jj++)
      matrix_fast_set_element(zz,ii,jj,nrow_zz, matrix_fast_get_element(xx, ii,jj,nrow_xx) + 
						matrix_fast_get_element(yy, ii,jj,nrow_yy));
  return;
}

/*****************************************************************************/
/*
void
matrix_get_row(Matrix *m, index_t i, Matrix *v)
{
  // Extracts the ith row of m to the column vector v
  index_t j;
  const index_t nc=numcols(m);
#ifdef _DBG_
 len=numrows(v), 
  if (len!=nc)
    error("Incompatible dimensions in matrix_get_row()");
#endif
  for (j=0; j<nc; j++)
    matrix_set_element(v,j,0, matrix_get_element(m,i,j));
  return;
}
*/
int
matrix_is_symmetric(Matrix *xx)
{
  //  Checks a matrix for symmetry, aindex_t with other basic checks.
  //  Returns 1 if the matrix is
  //  symmetric, 0 if not symmetric (if something else is wrong).

  int retval = 1;
  index_t ii, jj, nrow_xx = numrows(xx), ncol_xx = numcols(xx);

  Matrix *yy = matrix_new(nrow_xx, ncol_xx);

  matrix_transpose(xx, yy);
  matrix_scalar_multiply(yy, -1.0, yy);
  matrix_add(xx, yy, yy);
  for (ii=0; ii<nrow_xx; ii++)
    for (jj=0; jj<ncol_xx; jj++)
      if (matrix_get_element(yy, ii, jj) > DBL_EPSILON) 
	retval = 0;

  matrix_free(yy);
  return retval;
}

void
matrix_ADJUST(Matrix *xx, index_t kk)
{
  //  A function needed to sweep a matrix.  See Goodnight,
  //  33(3) Am. Stat. 149 (1979) for complete explanation.
  //  Function performs ADJUST(kk) on matrix xx.

  index_t ii, jj;
  double aa, aa_kk, aa_ik;
  index_t nrow_xx = numrows(xx);
  index_t ncol_xx = numcols(xx);

  //  Adjust row kk.
  aa_kk = matrix_get_element(xx, kk, kk);
  for (jj=(kk+1); jj<ncol_xx; jj++){
    aa = matrix_get_element(xx, kk, jj);
    matrix_set_element(xx, kk, jj, aa/aa_kk);
  }
  matrix_set_element(xx, kk, kk, 1.0);

  //  Adjust rows != kk.
  for (ii=0; ii<nrow_xx; ii++){
    if (ii == kk) 
	continue;
    aa_ik = matrix_get_element(xx, ii, kk);
    matrix_set_element(xx, ii, kk, 0.0);
    for (jj=(kk+1); jj<ncol_xx; jj++){
      aa = matrix_get_element(xx, ii, jj);
      matrix_set_element(xx, ii, jj, aa-(aa_ik * matrix_get_element(xx, kk, jj)));
    }
  }
}


/*****************************************************************************/

void matrix_identity(Matrix *xx)
{
  index_t ii, nr=numrows(xx), nc=numcols(xx);
  if (nr!=nc)
    error("Non-square matrix in matrix_identity()");
  matrix_fastset(xx,0);
  for (ii=0; ii<nr; ii++)
    matrix_set_element(xx,ii,ii, 1.0);
  return;
}
#ifdef _LAPACK_INVERSE_
void matrix_inverse(Matrix *X, Matrix *X_inverse, Matrix *Xsamedims)
{
  int n=numrows(X), e_code, ipiv[n];

  // Need to set X_inverse to the identity matrix on input:
  matrix_identity(X_inverse);

  // Copy X to Xsamedims (error check for dims inside matrix_copy):
  matrix_copy(X, Xsamedims);

  // Compute: Solution to a real system of linear equations: A * X = B
  // Where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
  // The LU decomposition with partial pivoting and row interchanges is
  // used to factor A as A = P * L * U,
  // where P is a permutation matrix, L is unit lower triangular, and U is
  // upper triangular.  The factored form of A is then used to solve the
  // system of equations A * X = B.
  //
  // N    = The number of linear equations, i.e., numrows(A)
  // NRHS = The number of right hand sides, i.e., numcols(B)
  //
  // A    = LDA-by-N matrix, the leading N-by-N matrix of A is the 
  //        coefficient matrix A. On exit, the factors L and U from the
  //        factorization. A = P*L*U
  // LDA = The leading dimension of the array A (LDA >= max(1,N))
  //
  // IPIV = N-vector containing the pivot indices that define P;
  //        row i of the matrix was interchanged with row IPIV(i)
  //
  // B    = LDB-by-NRHS matrix, the leading N-by-NRHS matrix of B is the
  //        right hand side matrix. On exit, the N-by-NRHS solution X.
  //
  // LDB = The leading dimension of the array B (LDB >= max(1,N))
  // INFO  =0 => Successful exit
  //       <0 => If INFO = -i, the i-th argument had an illegal value
  //       >0 => If INFO = i, U(i,i) is exactly zero.  The factorization
  //             has been completed, but the factor U is exactly
  //              singular, so the solution could not be computed.

//dgesv(n,n,Xsamedims,n,ipiv,X_inverse,n,&e_code);               // C version
  F77_CALL(dgesv)(&n,&n,Xsamedims,&n,ipiv,X_inverse,&n,&e_code); // R version

  if (!e_code)
    return;
  if (e_code<0)
    error("Singular value in mat_inverse.\n");
  else 
    error("Illegal value in mat_inverse.\n");
  return;
}
#endif
#ifdef _LAPACK_DET_
double matrix_determinant(Matrix *X, Matrix *Xsamedims, unsigned useLog)
{
  // Returns = log(det(X)) if useLog!=0, det(X) if useLog==0
  
  // Compute: An LU factorization of a general MxN matrix A
  // using partial pivoting with row interchanges.
  //
  //SUBROUTINE DGETRF(M,N,A,LDA,IPIV,INFO)
  //
  // M = numrows(A)
  // N = numcols(A)
  // A = LDA-by-N matrix, with the leading M-by-N submatrix to be factored
  //     On exit, the factors L and U from the factorization
  //     A = P*L*U; the unit diagonal elements of L are not stored.
  // LDA  = Leading dimension of the array A.  LDA >= max(1,M)
  // IPIV = (output) INTEGER array, dimension (min(M,N))
  //        The pivot indices; for 1 <= i <= min(M,N), row i of the
  //        matrix was interchanged with row IPIV(i).
  // INFO  =0 => successful exit
  //       <0 => If INFO = -i, the i-th argument had an illegal value
  //       >0 => If INFO = i, U(i,i) is exactly zero. The factorization
  //             has been completed, but the factor U is exactly
  //             singular, and division by zero will occur if it is used
  //             to solve a system of equations.
  
  int m = numrows(X), n = numcols(X); 

  // Copy to a deformable version (error check for dims in matrix_copy)
  matrix_copy(X,Xsamedims);

  int ii, info, sign, ipiv[(m<n)?m:n];
  double modulus;

//dgetrf(m,n,Xsamedims,m,ipiv,&info);              // C version
  F77_CALL(dgetrf)(&m,&n,Xsamedims,&m,ipiv,&info); // R version

  if (!info)
  {//LU-decomposition successful:
    sign = 1;
    for (ii=0; ii<n; ii++)
      if (ipiv[ii] != (ii+1))
	sign = -sign;
    if (useLog)
    {
      modulus = 0.0;
      for (ii=0; ii<n; ii++){
        double d_ii = Xsamedims[ii*(n+1)];
        modulus += log( d_ii<0 ? -d_ii : d_ii );
	if (d_ii<0)
	  sign = -sign;
      }
    } else {
      modulus = 1.0;
      for (ii=0; ii<n; ii++)
        modulus *= Xsamedims[ii*(n+1)];
      if (modulus<0){
        modulus = -modulus;
        sign = -sign;
      }
    }
    // Customized for Jim's application where SIGMA should be p.d.:
    if (sign<0)
      error("Matrix not positive definite in matrix_determinant()");

    return modulus;

  } else 
  {//LU-decomposition failed:
    if (info<0)
      error("Illegal value in matrix_determinant()");
    else
      return (useLog ? R_NegInf : 0.0 );
  }
  return R_NegInf; // Can never reach here :)
}
#endif
/*****************************************************************************/

void
matrix_DOOLITTLE(Matrix *xx, index_t kk)
{
  //  A function needed to find the determinant and to
  //  find the cholesky decomposition of a matrix.  See Goodnight,
  //  33(3) Am. Stat. 149 (1979) for complete explanation.
  //  Function performs DOOLITTLE(kk) on matrix xx.
  //  Function will only work on a square matrix.

  double aa_ij, aa_ik, aa_kk, aa_kj;
  index_t ii, jj, nrow_xx = numrows(xx), ncol_xx = numcols(xx);

  //  Adjust rows below kk
  aa_kk = matrix_get_element(xx, kk, kk);
  for (ii=(kk+1); ii<nrow_xx; ii++){
    aa_ik = matrix_get_element(xx, ii, kk);
    for(jj=(kk+1); jj<ncol_xx; jj++){
      aa_ij = matrix_get_element(xx, ii, jj);
      aa_kj = matrix_get_element(xx, kk, jj);
      matrix_set_element(xx, ii, jj, aa_ij - ((aa_ik/aa_kk)*aa_kj));
    }
    matrix_set_element(xx, ii, kk, 0.0);
  }
}

/*****************************************************************************/
#ifdef _LAPACK_CHOL_
void matrix_cholesky(Matrix  *X, Matrix *Y)
{
  //  Sets Y equal to the cholesky decomp of X.  Note per the definition,
  //  the cholesky decomp is an upper triangular matrix.

  int i,j, m=numrows(X),n=numcols(X);
  if (n!=m)
    error("Non-square matrix in matrix_cholesky()");

  // Copy X to Y (error check for dims inside matrix_copy)
  matrix_copy(X,Y);

  // Zero out lower triangle
  for (j=0; j<n; j++)
    for (i=j+1; i<n; i++)
      matrix_set_element(Y,i,j,0.0);

  // Compute: Cholesky factorization of Y (upper triangular)
  //
  //SUBROUTINE DPOTRF(UPLO,N,A,LDA,INFO)
  //
  // UPLO = 'U' => Upper triangle stored, 'L' => Lower triangle
  // N = numrows(A)
  // A = LDA-by-N matrix, leading N-by-N matrix to be factored
  //     On exit, the Cholesky factor
  // LDA = Leading dimension of A
  // INFO =0 => Successful exit
  //      >0 => if INFO = -i the ith argument had an illegal value
  //      <0 => if INFO = i the leading minor of order i is not p.d.
  //   
  F77_CALL(dpotrf)("Upper",&m,Y,&m,&i FCONE); 
//F77_CALL(chol)(X,&n,&n,Y,i);
  if (i!=0){
    if (i>0)
       error("Leading minor is not positive definite in matrix_cholesky()");
    error("Illegal value in matrix_cholesky()");
  }
  return;
}
#endif
#ifndef _LAPACK_CHOL_
void
matrix_cholesky(Matrix  *xx, Matrix *yy)
{
  //  Sets yy equal to the cholesky decomp of xx.  Note per the definition,
  //  the cholesky decomp is an upper triangular matrix.

  index_t kk, jj;
  double aa;

  matrix_get_set_block(yy, 0, numrows(yy)-1, 0, numcols(yy)-1, xx, 0, numrows(xx)-1, 0, numcols(xx)-1);
  for (kk=0; kk<(numrows(yy)-1); kk++){
    matrix_DOOLITTLE(yy, kk);
    aa = matrix_get_element(yy, kk, kk);
    for (jj=kk; jj<numcols(yy); jj++){
      matrix_set_element(yy, kk, jj, matrix_get_element(yy, kk, jj)/sqrt(aa));
    }
  }
  matrix_set_element(yy, numcols(yy)-1, numcols(yy)-1,
		     sqrt(matrix_get_element(yy, numcols(yy)-1, numcols(yy)-1)));
}
#endif
/*****************************************************************************/

void
matrix_copy(Matrix *ToBeCopied, Matrix *Copy)
{
  //  In R notation, Copy <-ToBeCopied
  const index_t nrow=numrows(Copy);
  const index_t ncol=numcols(Copy);
  if ((nrow!=numrows(ToBeCopied))||(ncol!=numcols(ToBeCopied)))
    error("Incompatible dims in matrix_copy()");
  // Much faster copying:
  memcpy(Copy,ToBeCopied,sizeof(double)*nrow*ncol);
  return;
}

void
matrix_subtract(Matrix *xx, Matrix *yy, Matrix *zz)
{
  //  In R notation, implements zz <- xx - yy
  index_t ii, jj, nrow_xx = numrows(xx), ncol_xx = numcols(xx);
  for (ii=0; ii<nrow_xx; ii++)
    for (jj=0; jj<ncol_xx; jj++)
      matrix_set_element(zz, ii, jj, matrix_get_element(xx, ii, jj) - matrix_get_element(yy, ii, jj));
  return;
}

int
matrix_assert_row_vec(Matrix *xx)
{
  //  Runs several standard checks on a matrix that is supposed
  //  to be a row vetor.  Returns 0 if anything is amiss, 1 if
  //  all is well.

  if (xx == NULL)
    error("Error:  Row vector that should not be NULL is NULL.\n");
  
  if (numcols(xx) <= 0)
    error("Error:  Row vector has fewer than 1 column.\n");
  
  if (numrows(xx) != 1)
    error("Error:  Row vector has number of rows not equal to 1.\n");
  
  return 1;
}

int
matrix_assert_vec(Matrix *xx)
{
  //  Runs several standard checks on a matrix that is supposed
  //  to be a vector (column or row).  Returns 0 if anything is amiss,
  //  1 if all is well.

  if (xx == NULL)
    error("Error:  Vector that should not be NULL is NULL.\n");
  
  index_t minn = min(numrows(xx), numcols(xx));
  index_t maxx = max(numrows(xx), numcols(xx));

  if (minn != 1)
    error("Error:  Vector has dimension less than 1.\n");
  if (maxx < 1)
    error("Error:  Vector has no room for elements.\n");
  
  return 1;
}


void
matrix_sum_xx_m_mu(Matrix *yy, Matrix *xx, Matrix *mu_vec)
{
  //  Subtracts the row vector mu_vec from each row of XX,
  //  then multiplies the transpose of the resulting matrix
  //  by itself, sets the result equal to YY.  In R notation,
  //  YY <- t(XX - t(mu_vec)) %*% (XX-t(mu_vec)).  User
  //  allocates all memory.

  // First, set y to zero matrix:
  matrix_fastset(yy,0);

  index_t ii, jj, kk;
  const index_t nrow_mu = numrows(mu_vec);
  const index_t nrow_yy = numrows(yy);
  const index_t nrow_xx = numrows(xx);
  const index_t ncol_xx = numcols(xx);

	//double tmp_jj_factor;

  //  Begin by filling in the diagonals and upper triangle.
  for (ii=0; ii<nrow_xx; ii++){
    for (jj=0; jj<ncol_xx; jj++){

			// Cache the unchanging elements (added 12/08/08) -- REMOVED, ACTUALLY SLOWER.
			// tmp_jj_factor = (matrix_fast_get_element(xx, ii,jj,nrow_xx) - matrix_fast_get_element(mu_vec, 0,jj,nrow_mu));

      for (kk=jj; kk<ncol_xx; kk++){
	matrix_fast_increment_element(yy, jj,kk,nrow_yy,
			(matrix_fast_get_element(xx, ii,jj,nrow_xx) - matrix_fast_get_element(mu_vec, 0,jj,nrow_mu))*
			(matrix_fast_get_element(xx, ii,kk,nrow_xx) - matrix_fast_get_element(mu_vec, 0,kk,nrow_mu)));

      }
    }
  }
  //  Now fill in the part below the diagonal.
  for (jj=1; jj<ncol_xx; jj++){
    for (kk=0; kk<jj; kk++){
      matrix_fast_set_element(yy,jj,kk,nrow_yy, matrix_fast_get_element(yy,kk,jj,nrow_yy));
    }
  }
  return;
}


void
matrix_sum_xx_m_mu_by_precinct(Matrix * const SS, 
                               Matrix * const OMEGAS, 
                               Matrix * const mu_mat_cu)
{
//
// Computes: \sum_{i=1}^{n} (\omega_{i}-\mu_{i})(\omega_{i}-\mu_{i})^{T}
// Implemented as: matrix_sum_xx_m_mu_by_precinct(SS, OMEGAS, mu_mat_cu);
//
  const index_t n = numrows(OMEGAS);
  if (numrows(mu_mat_cu)!=n)
    error("OMEGAS and mu_mat_cu must have the same number of rows");

  const index_t p = numcols(OMEGAS);
  if (numcols(mu_mat_cu)!=p)
    error("OMEGAS and mu_mat_cu must have the same number of cols");

  if ((numrows(SS)!=p) || (numcols(SS)!=p))
    error("SS must be p x p");

  double tmp_omega_i_minus_mu_i[p];
  index_t ii, jj, kk;

  // First, set SS to zero matrix:
  matrix_fastset(SS,0);

  // Precinct-by-precinct:
  for (ii=0; ii<n; ii++){
    // First, do the omega_i - mu_i bit:
    for (jj=0; jj<p; jj++){
      tmp_omega_i_minus_mu_i[jj] = matrix_fast_get_element(OMEGAS,ii,jj,n) - matrix_fast_get_element(mu_mat_cu,ii,jj,n);
    }
    // Now do the xx^{T} bit:
    for (jj=0; jj<p; jj++){
      for (kk=0; kk<p; kk++){
        matrix_fast_increment_element(SS,jj,kk,p, (tmp_omega_i_minus_mu_i[jj]*tmp_omega_i_minus_mu_i[kk]));
      }
    }
  }

#ifdef _DBG_
  Rprintf("Computing sum_{i=1}^{n}(omega_i-mu_i)(omega_i-mu_i)^{T}:\n\n");
  Rprintf("OMEGAS:\n");
  matrix_print_subset(OMEGAS,0,1,0,numcols(OMEGAS)-1);
  Rprintf("MUS:\n");
  matrix_print_subset(mu_mat_cu,0,1,0,numcols(mu_mat_cu)-1);
  Rprintf("Result:\n");
  matrix_print_all(SS);
#endif
  
  return;
}
