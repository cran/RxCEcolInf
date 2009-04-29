
#include "jimsimatrix.h"

/****************************************************************/
/*******************  INT MATRIX ALGEBRA LIBRARY ****************/
/****************************************************************/

Matrix_int *
rmatrix_new_int(index_t nrow, index_t ncol)
{
  // R-allocation (has its own error checks):
  int *mat = (int *) R_alloc((nrow*ncol)+2,sizeof(int));

  // Initialize to 0:
  memset(mat,0,sizeof(int)*((nrow*ncol)+2));

  // Return pointer to array of pointers to rows
  mat[0] = (int)ncol;
  mat[1] = (int)nrow;
  return mat+2;
}

Matrix_int *
matrix_new_int(index_t nrow, index_t ncol)
{
  // C-allocation (has its own error checks):
  int *mat = (int *) malloc((size_t)((nrow*ncol)+2)*sizeof(int));
  if (mat==NULL) 
    error("Allocation failure in matrix_new_int()");

  // Initialize to 0:
  memset(mat,0,sizeof(int)*((nrow*ncol)+2));

  // Return pointer to array of pointers to rows
  mat[0] = (int)ncol;
  mat[1] = (int)nrow;
  return mat+2;
}

void
matrix_free_int(Matrix_int *xx)
{
  if (xx!=NULL)
    free(xx-2);
  return;
}

void
sample_equ_pr_wo_replace_int(Matrix_int *bigvec_int,
			     Matrix_int *smallvec_int)
{
  // Sets smallvec_int a set of random draws without replacement
  // from the numbers 0 to numcols_int(bigvec_int) - 1.  smallvec_int
  // and bigvec_int must be row vectors of integers

  index_t uplim = numcols_int(bigvec_int), ii, ranint;
  double randou;

  for (ii=0; ii<uplim; ii++)
    matrix_set_int_element(bigvec_int, 0, ii, ii);
  
  index_t ncol_svi = numcols_int(smallvec_int);
  for (ii=0; ii<ncol_svi; ii++){

    do {
      randou = runif(0.0, 1.0);
    } while (randou >= 1.0);

    ranint = (index_t) (uplim * randou);
    matrix_set_int_element(smallvec_int, 0,ii,
			   matrix_get_int_element(bigvec_int, 0, ranint));
    --uplim;
    matrix_set_int_element(bigvec_int, 0,ranint,
			   matrix_get_int_element(bigvec_int, 0, uplim));
  }
  return;
}


void
matrix_print_int_element(Matrix_int *xx_int, index_t row, index_t col)
{
  Rprintf("%d\t", matrix_get_int_element(xx_int, row, col));
}

void
matrix_print_all_int(Matrix_int *xx_int)
{
  index_t ii, jj, nrow_xx=numrows_int(xx_int), ncol_xx=numcols_int(xx_int);

  for (ii=0; ii<nrow_xx; ii++){
    for (jj=0; jj<ncol_xx; jj++){
      matrix_print_int_element(xx_int, ii, jj);
    }
    Rprintf("\n");
  }
  return;
}
