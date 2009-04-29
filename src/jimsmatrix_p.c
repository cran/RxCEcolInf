#include "jimsmatrix_p.h"

// Constructors:

// Single pointer constructors:
/*
Matrix_p *rmatrix_p_new(index_t len);
Matrix_p *matrix_p_new(index_t len);
*/

Matrix_p *
rmatrix_p_new(index_t len)
{ // Allocates (R) storage for len Matrix pointers (not the matrices):
  Matrix_p *mat_p = (Matrix_p *)R_alloc(len,sizeof(Matrix_p));
  if (mat_p==NULL)
    error("Allocation failure in rmatrix_p_new()");

  // Set all pointers to null:
  index_t i;
  Matrix *null_mat_ptr = NULL;
  for (i=0; i<len; i++)
    set_mat_p_ptr(mat_p,i, null_mat_ptr);
  return mat_p;
}
Matrix_p *
matrix_p_new(index_t len)
{ // Allocates (C) storage for len Matrix pointers (not the matrices):
  Matrix_p *mat_p = (Matrix_p *)malloc((size_t)len*sizeof(Matrix_p));
  if (mat_p==NULL)
    error("Allocation failure in matrix_p_new()");

  // Set all pointers to null:
  index_t i;
  Matrix *null_mat_ptr = NULL;
  for (i=0; i<len; i++)
    set_mat_p_ptr(mat_p,i, null_mat_ptr);
  return mat_p;
}

// Destructors:

// Single pointer destructors:
/*
void matrix_p_free(Matrix_p *mat_p, index_t len);
*/

void matrix_p_free(Matrix_p *mat_p, index_t len)
{
  // Frees the matrices pointed to by each pointer, as well as the pointers:
  index_t i;
  Matrix *tmp;
  for (i=0; i<len; i++){
    tmp = get_mat_p_ptr(mat_p,i);
    matrix_free(tmp);
  }
  free(mat_p);
  return;
}

// Mutators:

// Single pointer mutators:
/*
void set_mat_p_ptr(Matrix_p *mat_p, index_t i, Matrix *mat);
*/

void set_mat_p_ptr(Matrix_p *mat_p, index_t i, Matrix *mat)
{
  mat_p[i] = mat;
  return;
}

// Selectors:

// Single pointer selectors:
/*
Matrix *get_mat_p_ptr(Matrix_p *mat_p, index_t i);
*/

Matrix *get_mat_p_ptr(Matrix_p *mat_p, index_t i)
{
  return mat_p[i];
}

// Methods:

// Single pointer printing methods:
/*
void matrix_p_print_all(Matrix_p *mat_p, index_t len);
void matrix_p_print_subset(Matrix_p *mat_p, index_t lo, index_t hi, index_t rl,
                           index_t rh, index_t cl, index_t ch);
*/

void matrix_p_print_all(Matrix_p *mat_p, index_t len)
{
  index_t i;
  for (i=0; i<len; i++)
    matrix_print_all(get_mat_p_ptr(mat_p,i));
  return;
}
void matrix_p_print_subset(Matrix_p *mat_p, index_t lo, index_t hi, index_t rl, index_t rh, index_t cl, index_t ch)
{
  index_t i;
  for (i=lo; i<=hi; i++)
    matrix_print_subset(get_mat_p_ptr(mat_p,i),rl,rh,cl,ch);
  return;
}

// Single pointer list unpacking:
/*
Matrix_p *matrix_p_unpack_new(SEXP rLst);
Matrix_p *rmatrix_p_unpack_new(SEXP rLst);
*/

Matrix_p *matrix_p_unpack_new(SEXP rLst)
{
  index_t i, len=length(rLst);
  Matrix *mat;
  Matrix_p *mat_p = matrix_p_new(len);
  for (i=0; i<len; i++){
    mat = matrix_unpack_new(VECTOR_ELT(rLst,i));
    set_mat_p_ptr(mat_p,i, mat);
  }
  return mat_p;
}
Matrix_p *rmatrix_p_unpack_new(SEXP rLst)
{
  index_t i, len=length(rLst);
  Matrix *mat;
  Matrix_p *mat_p = rmatrix_p_new(len);
  for (i=0; i<len; i++){
    mat = rmatrix_unpack_new(VECTOR_ELT(rLst,i));
    set_mat_p_ptr(mat_p,i, mat);
  }
  return mat_p;
}

