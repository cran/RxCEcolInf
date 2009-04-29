#include "jimsmisc.h"

SEXP
getListElement(SEXP list, char *str)
{
  /*  A most excellent function to extract something from a list
      passed from R to C. */
   SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
   int ii;
     
   for (ii = 0; ii < length(list); ii++)
     if (strcmp(CHAR(STRING_ELT(names, ii)), str) == 0) {
       elmt = VECTOR_ELT(list, ii);
       break;
     }
   return elmt;
}

Matrix *
create_log_factorial_lookup_table(const index_t imax)
{
  // Log-Factorial vector in the form:
  // [0]: log(0!), [1]: log(1!), ..., [imax]: log(imax!)

  Matrix * lfvec = rmatrix_new(1,imax+1); // Include room for 0!
  index_t i;
  double tmp = 0.0;

  // Avoid log(0):
  matrix_fast_set_element(lfvec,0,0,1, tmp); 

  // Compute the rest of the log-factorials iteratively:
  for (i=1; i<=imax; i++){
    tmp += log((double)i);
    matrix_fast_set_element(lfvec,0,i,1, tmp); 
  }
  return lfvec;
}

void 
store_internals(Matrix * const ToBeStored, 
                Matrix * const Store_mat, 
                long * const colnum_store)
{
  //  ToBeStored is a (number of precincts) by numcells_pt matrix holding
  //     the current values of either the THETAS or the internal cell counts.
  //     Store_mat is a ((number of precincts)*numcells_pt) by (number of draws
  //     to be saved) matrix.  colnum_store is a pointer to the integer
  //     column number where the values of ToBeStored will be placed in
  //     Store_mat.  In the colnum_store'th column, the first row will hold
  //     the value from the first precinct, row 1, column 1; the second row will
  //     hold the value from the first precinct, row 1, column 2; etc.

  long ii = 0, jj = 0, row_num = 0, col_num = *colnum_store;
  long nrow_TBS = numrows(ToBeStored);
  long ncol_TBS = numcols(ToBeStored);

  for (ii=0; ii<nrow_TBS; ii++){
    for (jj=0; jj<ncol_TBS; jj++){
      matrix_set_element(Store_mat, row_num, col_num,
			 matrix_get_element(ToBeStored, ii, jj));
      row_num++;
    }
  }

  *colnum_store = col_num + 1;
}

/****************************************************************************/
/*****************  TRANSFORMATION FUNCTION DECLARATIONS  *******************/
/****************************************************************************/

void
THETAS_to_OMEGAS(Matrix * const THETAS, 
                 Matrix * const OMEGAS, 
		             const index_t numrows_pt, 
                 const index_t numcols_pt)
{
  //  Takes matrix of THETAS where each row is a precinct and transforms
  //  them into OMEGAS using the precinct-table-row-by-row additive logistic
  //  transformation.  User allocates memory to OMEGAS.

  index_t ii, jj, kk;
  const index_t numcols_pt_m1 = numcols_pt-1;
  const index_t nrow_THETAS = numrows(THETAS);
  const index_t nr_OMEGAS = numrows(OMEGAS);
  double aa_ref;

  for (ii=0; ii<nrow_THETAS; ii++){
    for (jj=0; jj<numrows_pt; jj++){
      aa_ref = matrix_fast_get_element(THETAS, ii,((jj+1)*numcols_pt)-1,nrow_THETAS);
      for (kk=0; kk<numcols_pt_m1; kk++){
	matrix_fast_set_element(OMEGAS, ii,(jj*numcols_pt_m1)+kk,nr_OMEGAS,
	  log(matrix_fast_get_element(THETAS, ii,(jj*numcols_pt)+kk,nrow_THETAS)/aa_ref));
      }
    }
  }
  return;
}

void
check_ep_validity(Matrix * const NNs,
                  Matrix * const MMs,
                  Matrix * const KKs,
                  const index_t numcells_pt,
                  const index_t numrows_pt,
                  const index_t numcols_pt)
{
// Checks that the exit poll (KKs), and non-ep 'voters' (MMs), sum to
// the totals in NNs. This doesn't check the bounds condition, the 
// NNs state may be invalid, but that check is performed separately
// by bounds_check_all.

  const index_t n=numrows(NNs);
#ifdef _DBG_
  const index_t nr_MMs=numrows(MMs);
  const index_t nr_KKs=numrows(KKs);
  if (nr_MMs!=n)
    error("Number of rows in MMs (%u) doesn't match that of NNs (%u)\n",nr_MMs,n);
  if (nr_KKs!=n)
    error("Number of rows in KKs (%u) doesn't match that of NNs (%u)\n",nr_KKs,n);
#endif
  double tmp_NN, tmp_MM, tmp_KK;
  index_t ii, jj, kk, ll;

  // Precinct by precinct:
  for (ii=0; ii<n; ii++){

    // Check each cell is valid:
    for (jj=0; jj<numcells_pt; jj++){

      // Extract the cell counts:
      tmp_NN = matrix_fast_get_element(NNs,ii,jj,n);
      tmp_MM = matrix_fast_get_element(MMs,ii,jj,n);
      tmp_KK = matrix_fast_get_element(KKs,ii,jj,n);

      // Check for negatives and sum equality:
      if (tmp_NN<0.0 || tmp_MM<0.0 || tmp_KK<0.0 || (tmp_NN!=(tmp_MM+tmp_KK))){

        Rprintf("Failed validity check: invalid state in precinct %u, cell %u\n(NN_ij,MM_ij,KK_ij)=(%g,%g,%g)\n",
              ii,jj,tmp_NN,tmp_MM,tmp_KK);

        Rprintf("\nNNs for this precinct:\n");
        for (kk=0; kk<numrows_pt; kk++){
          for (ll=0; ll<numcols_pt; ll++)
            Rprintf("%g\t",matrix_get_element(NNs,ii,kk*numcols_pt+ll));
          Rprintf("\n");
        }

        Rprintf("\nMMs for this precinct:\n");
        for (kk=0; kk<numrows_pt; kk++){
          for (ll=0; ll<numcols_pt; ll++)
            Rprintf("%g\t",matrix_get_element(MMs,ii,kk*numcols_pt+ll));
          Rprintf("\n");
        }

        Rprintf("\nKKs for this precinct:\n");
        for (kk=0; kk<numrows_pt; kk++){
          for (ll=0; ll<numcols_pt; ll++)
            Rprintf("%g\t",matrix_get_element(KKs,ii,kk*numcols_pt+ll));
          Rprintf("\n");
        }

        error("Failed validity check: invalid state in precinct %u, cell %u\n(NN_ij,MM_ij,KK_ij)=(%g,%g,%g)\n",
              ii,jj,tmp_NN,tmp_MM,tmp_KK);
      }

    } // Finish cell jj check
  } // Finish precinct ii check
  return;
}

void 
adjust_rho_vec(Matrix * const rho_vec, 
               SEXP acc_THETAS_t_vec)
{
  // Adjusts the multiplier of SIGMA_cu used in the t proposal distribution 
  // for the THETAS M-H algorithm draws.
  double acc_cu;
  double *acc_t_p = REAL(acc_THETAS_t_vec);
  const index_t len_rho = numcols(rho_vec);
  const index_t nr_rho  = numrows(rho_vec);
  index_t ii;

  for (ii=0; ii<len_rho; ii++){
    acc_cu = acc_t_p[ii];
    //  First see if accepting too many draws
    if ((0.4<=acc_cu)&&(acc_cu<0.5)){
      matrix_fast_set_element(rho_vec,0,ii,nr_rho, matrix_fast_get_element(rho_vec,0,ii,nr_rho)*1.1);
      continue;
    }
    if ((0.5 <= acc_cu) && (acc_cu < 0.7)){
      matrix_fast_set_element(rho_vec,0,ii,nr_rho, matrix_fast_get_element(rho_vec,0,ii,nr_rho)*1.4);
      continue;
    }
    if (0.7 <= acc_cu){
      matrix_fast_set_element(rho_vec,0,ii,nr_rho, matrix_fast_get_element(rho_vec,0,ii,nr_rho)*1.7);
      continue;
    }
    //  Now see if accepting too few draws
    if((0.2 < acc_cu) && (acc_cu <= 0.3)){
      matrix_fast_set_element(rho_vec,0,ii,nr_rho, matrix_fast_get_element(rho_vec,0,ii,nr_rho)*0.9);
      continue;
    }
    if((0.1 < acc_cu) && (acc_cu <= 0.2)){
      matrix_fast_set_element(rho_vec,0,ii,nr_rho, matrix_fast_get_element(rho_vec,0,ii,nr_rho)*0.7);
      continue;
    }
    if (acc_cu<=0.1){
      matrix_fast_set_element(rho_vec,0,ii,nr_rho, matrix_fast_get_element(rho_vec,0,ii,nr_rho)*0.5);
      continue;
    }
  }
  return;
}

void 
adjust_acc_vector(SEXP acc_vec, 
                  Matrix * const count_use_vec)
{
  //  Turns counts of number of accepted proposals
  //  in acc_vec into fractions based on the number
  //  of times the proposal method was used.
  const index_t len = length(acc_vec);
  const index_t nr_cuv = numrows(count_use_vec);
  double *acc_p = REAL(acc_vec);
  index_t i;
  for (i=0; i<len; i++)
    acc_p[i] /= matrix_fast_get_element(count_use_vec,0,i,nr_cuv);
  return;
}

void 
adjust_acc_matrix(SEXP acc_mat, 
                  Matrix * const count_use_mat)
{
  //  Turns counts of number of accepted proposals
  //  in acc_mat into fractions based on the number
  //  of times the proposal method was used.
  const index_t nn = nrows(acc_mat);
  const index_t nr = ncols(acc_mat);
  const index_t nr_cum = numrows(count_use_mat);
  const index_t nc_cum = numcols(count_use_mat);
  if ((nn<nr_cum)||(nr<nc_cum))
		error("acc_mat too small to hold acceptance fractions");
  double *acc_p = REAL(acc_mat);
  index_t i, j;
  for (i=0; i<nr_cum; i++)
    for (j=0; j<nc_cum; j++)
	    acc_p[i+nn*j] /= matrix_fast_get_element(count_use_mat,i,j,nr_cum);
  return;
}

void
report_dry_run(SEXP acc_THETAS_t_vec, 
               SEXP acc_THETAS_Diri_vec,
	             const long num_run)
{
  // Prints the results of a dry run to the screen:
  Rprintf("\nRun %d resulted in the following:\n",num_run+1);
  Rprintf("\tFor the t proposal:\n");
  Rprintf("\t\tLowest acceptance fraction = %f\n",Rmatrix_min(acc_THETAS_t_vec));
  Rprintf("\t\tHighest acceptance fraction = %f\n",Rmatrix_max(acc_THETAS_t_vec));
  Rprintf("\t\tFraction of precints, acceptance fractions < .2 = %f\n",
	  Rmatrix_get_fraction_under_c(acc_THETAS_t_vec,0.2));
  Rprintf("\t\tFraction of precints, acceptance fractions > .5 = %f\n",
	  Rmatrix_get_fraction_over_c(acc_THETAS_t_vec,0.5));
  Rprintf("\n\tFor the Diri proposal:\n");
  Rprintf("\t\tLowest acceptance fraction = %f\n",Rmatrix_min(acc_THETAS_Diri_vec));
  Rprintf("\t\tHighest acceptance fraction = %f\n",Rmatrix_max(acc_THETAS_Diri_vec));
  Rprintf("\t\tFraction of precints, acceptance fractions < .2 = %f\n",
	  Rmatrix_get_fraction_under_c(acc_THETAS_Diri_vec,0.2));
  Rprintf("\t\tFraction of precints, acceptance fractions > .5 = %f\n",
	  Rmatrix_get_fraction_over_c(acc_THETAS_Diri_vec,0.5));
  return;
}

void
check_bounds_all(Matrix * const NNs, 
                 Matrix * const NNbounds, 
                 const index_t numcells_pt)
{
  //  Checks all values of NNs to make sure each is within
  //  its bounds.  Crashes the function if one wanders outside.
  //  User allocates all memory.

  index_t mm, jj;
  const index_t nrow_NNs = numrows(NNs);
  const index_t ncol_NNs = numcols(NNs);

  for (mm=0; mm<nrow_NNs; mm++){
    for (jj=0; jj<ncol_NNs; jj++){
      if (!check_bounds(matrix_fast_get_element(NNs,mm,jj,nrow_NNs), 
		NNbounds, mm, jj, numcells_pt)){
	Rprintf("Error:  Count value outside bounds in precinct %d, position %d.\n", mm, jj);
	Rprintf("\tCurrent count in this position:  %f\n",
	       matrix_get_element(NNs, mm, jj));
	Rprintf("\tCorresponding lower bound:  %f\n",
	       matrix_get_element(NNbounds, mm, jj));
	Rprintf("\tCorresponding upper bound:  %f\n",
	       matrix_get_element(NNbounds, mm, jj+numcells_pt));
	error("Exiting\n");
      }
    }
  }
  return;
}


void
store_draws(Matrix * const NNs, 
	    Matrix * const mu_vec_cu, 
	    Matrix * const SIGMA_cu,
	    Matrix * const draws, 
	    long * const rownum_store)
{
  //  Stores the draw of mu_vec_cu, then of SIGMA_cu, then
  //  of the sum (over precincts) of the NNs matrix, all
  //  in one large storage matrix.
  //  The SIGMA draw is stored as variances and correlations (not
  //  covariances).  First all variances are stored; then all
  //  correlations along the first row, then all non-duplicate
  //  correlations along the second row, etc.
  //  User allocates all memory.

  index_t ii, jj;
  long rownum = *rownum_store, place = 0;
  const index_t nr_draws=numrows(draws);
  double diag1_sqrt, diag2_sqrt, temp_val;

  // Store the mu_vec_cu draw:
  const index_t nr_mu=numrows(mu_vec_cu);
  const index_t len_mu=numcols(mu_vec_cu);
  for (jj=0; jj<len_mu; jj++){
    matrix_fast_set_element(draws,rownum,place,nr_draws, matrix_fast_get_element(mu_vec_cu,0,jj,nr_mu));
    place++;
  }

  // Transform the SIGMA, store the draw, (diagonal first)
  const index_t nr_SIG=numrows(SIGMA_cu);
  const index_t nc_SIG=numcols(SIGMA_cu);
  for (ii=0; ii<nr_SIG; ii++){
    matrix_set_element(draws, rownum, place, matrix_get_element(SIGMA_cu, ii, ii));
    place++;
  }
  // Now the off diagonal terms:
  for (ii=0; ii<(nr_SIG-1); ii++){
    diag1_sqrt = sqrt(matrix_fast_get_element(SIGMA_cu,ii,ii,nr_SIG));
    for (jj=(ii+1); jj<nc_SIG; jj++){
      diag2_sqrt = sqrt(matrix_fast_get_element(SIGMA_cu,jj,jj,nr_SIG)) * diag1_sqrt;
      matrix_set_element(draws, rownum, place, matrix_fast_get_element(SIGMA_cu,ii,jj,nr_SIG)/diag2_sqrt);
      place++;
    }
  }

  // Store the sum, across the precincts, of the NNs
  const index_t nr_NNs=numrows(NNs);
  const index_t nc_NNs=numcols(NNs);
  for (jj=0; jj<nc_NNs; jj++){
    temp_val = 0.0;
    for (ii=0; ii<nr_NNs; ii++){
      temp_val += matrix_fast_get_element(NNs,ii,jj,nr_NNs);
    }
    matrix_set_element(draws, rownum, place, temp_val);
    place++;
  }

  *rownum_store = *rownum_store + 1;
  return;
}

void
store_the_draws(Matrix * const mu,
		SEXP mu_draws,
		Matrix * const SIGMA,
		SEXP SIGMA_draws,
		Matrix * const NNs,
		SEXP NNs_draws,
		SEXP LAMBDA_draws,
		SEXP TURNOUT_draws,
		SEXP GAMMA_draws,
		SEXP BETA_draws,
		const index_t R,
		const index_t C,
		long * const rownum_store,
    const int ncol_LAMBDA,
    const int ncol_BETA)
{
  index_t ii, jj;
  const long rownum = *rownum_store;
  const index_t nr_draws = nrows(mu_draws);// = nrow(*_draws);
  double diag1_sqrt, diag2_sqrt, temp_val;

  // Store the mu_vec_cu or eta_vec_cu draw:
  const index_t nr_mu =numrows(mu);
  const index_t len_mu=numcols(mu);
  double *mu_p = REAL(mu_draws);
  double *mu_q = &(mu_p[rownum]);
  for (jj=0; jj<len_mu; jj++)
    mu_q[nr_draws*jj] = matrix_fast_get_element(mu,0,jj,nr_mu);

  // Transform the SIGMA, store the draw, (diagonal first)
  index_t place=0;
  const index_t nr_SIG=numrows(SIGMA);
  const index_t nc_SIG=numcols(SIGMA);
  double *SIGMA_p = REAL(SIGMA_draws);
  double *SIGMA_q = &(SIGMA_p[rownum]);
  for (ii=0; ii<nr_SIG; ii++)
    SIGMA_q[nr_draws*ii] = sqrt(matrix_fast_get_element(SIGMA,ii,ii,nr_SIG));

  // Increment the placeholder:
  place += nr_SIG;

  // Now the off diagonal terms:
  for (ii=0; ii<(nr_SIG-1); ii++){
    diag1_sqrt = sqrt(matrix_fast_get_element(SIGMA,ii,ii,nr_SIG));
    for (jj=(ii+1); jj<nc_SIG; jj++){
      diag2_sqrt = sqrt(matrix_fast_get_element(SIGMA,jj,jj,nr_SIG)) * diag1_sqrt;
      SIGMA_q[nr_draws*place] = matrix_fast_get_element(SIGMA,ii,jj,nr_SIG)/diag2_sqrt;
      place++;
    }
  }

  // Store the sum, across the precincts, of the NNs
  // Note that this is done in row-major form:
  // NNr1c1_sum, NNr1c2_sum, ..., NNr1cC_sum, NNr2c1_sum etc.
  // since NNs is arranged in row-major form 
  // This is just apply(NNs,2,sum) in R code.
  // We want to extract NNs_p[row_num+nr_draws*place] so for faster indexing,
  // just get rid of the row_num bit:
  double *NNs_p = REAL(NNs_draws);
  double *NNs_q = &(NNs_p[rownum]);
  const index_t num_precincts = numrows(NNs);
  const index_t R_times_C = R*C;
  place=0;
  for (jj=0; jj<R_times_C; jj++){
    temp_val = 0.0;
    for (ii=0; ii<num_precincts; ii++){
      temp_val += matrix_fast_get_element(NNs,ii,jj,num_precincts);
    }
    NNs_q[nr_draws*place] = temp_val;
    place++;
  }

  // Compute the R*ncol_LAMBDA LAMBDA's (the fractions wrt rowsums excluding non-voters and eschew-ers):
  double *LAMBDA_p = REAL(LAMBDA_draws);
  double *LAMBDA_q = &(LAMBDA_p[rownum]);
  double *TURNOUT_p = REAL(TURNOUT_draws);
  double *TURNOUT_q = &(TURNOUT_p[rownum]);
  double *GAMMA_p = REAL(GAMMA_draws);
  double *GAMMA_q = &(GAMMA_p[rownum]);
  double *BETA_p = REAL(BETA_draws);
  double *BETA_q = &(BETA_p[rownum]);
  const index_t Cm1 = (C-1); // This never changes even if eschew=TRUE
  const index_t Cm1or2 = ncol_LAMBDA;
  double tmp;
  double tmp_full_rowsums[R];   // Row sums including all columns
  double tmp_voter_rowsums[R];  // Excluding non-voters and eschew-ers (last 1 or 2 columns)
  double gamma_denominator=0.0; // All cells excluding the last 1 or 2 columns
  // Compute the gamma denominator and rowsums so we can do everything in one loop:
  place=0;
  for (ii=0; ii<R; ii++){

    // Reset for row ii:
    tmp = 0.0;

    // Depending on eschew:
    for (jj=0; jj<Cm1or2; jj++)
      tmp += NNs_q[nr_draws*(jj+ii*C)];
    
    tmp_voter_rowsums[ii] = tmp;
    gamma_denominator += tmp;

    // Tack on the remaining 1 or 2 columns:
    for (jj=Cm1or2; jj<C; jj++)
      tmp += NNs_q[nr_draws*(jj+ii*C)];

    // Recored the total row sum for row ii:
    tmp_full_rowsums[ii] = tmp;
  }

  // Now all of the fractions in one go:
  place=0;
  for (ii=0; ii<R; ii++){
    for (jj=0; jj<Cm1or2; jj++){
      LAMBDA_q[nr_draws*place] = NNs_q[nr_draws*(jj+ii*C)]/tmp_voter_rowsums[ii];
      place++;
    }
    TURNOUT_q[nr_draws*ii] = tmp_voter_rowsums[ii]/tmp_full_rowsums[ii];
    GAMMA_q[nr_draws*ii]   = tmp_voter_rowsums[ii]/gamma_denominator;
  }

  // Do BETAs separately to avoid problems with place:
  place=0;
  for (ii=0; ii<R; ii++){
    for (jj=0; jj<Cm1; jj++){
      BETA_q[nr_draws*place] = NNs_q[nr_draws*(jj+ii*C)]/tmp_full_rowsums[ii];
      place++;
    }
  }

  // Increment the number of draws so far:
  *rownum_store = rownum+1;
  return;
}


Matrix_int *get_which_rc(const index_t nr_pt, const index_t nc_pt)
{
  const index_t nr_rows_mat = (index_t)choose(nr_pt,2.0);
  const index_t nr_cols_mat = (index_t)choose(nc_pt,2.0);
  const index_t nr_pt_m1 = nr_pt-1;
  const index_t nc_pt_m1 = nc_pt-1;

  // Allocate the matrices:
  Matrix_int * const which_rc = rmatrix_new_int(nr_rows_mat*nr_cols_mat, 4); // R-Allocation
  Matrix_int * const rows_mat = matrix_new_int(nr_rows_mat, 2); // C-allocation
  Matrix_int * const cols_mat = matrix_new_int(nr_cols_mat, 2); // C-allocation

  index_t ii, jj, rnum;
  rnum=0;
  for (ii=0; ii<nr_pt_m1; ii++){
    for (jj=(ii+1); jj<nr_pt; jj++){
      matrix_set_int_element(rows_mat,rnum,0, ii);
      matrix_set_int_element(rows_mat,rnum,1, jj);
      rnum++;
    }
  }
  rnum=0;
  for (ii=0; ii<nc_pt_m1; ii++){
    for (jj=(ii+1); jj<nc_pt; jj++){
      matrix_set_int_element(cols_mat,rnum,0, ii);
      matrix_set_int_element(cols_mat,rnum,1, jj);
      rnum++;
    }
  }
  rnum=0;
  for (ii=0; ii<nr_rows_mat; ii++){
    for (jj=0; jj<nr_cols_mat; jj++){
      matrix_set_int_element(which_rc,rnum,0, matrix_get_int_element(rows_mat,ii,0));
      matrix_set_int_element(which_rc,rnum,1, matrix_get_int_element(rows_mat,ii,1));
      matrix_set_int_element(which_rc,rnum,2, matrix_get_int_element(cols_mat,jj,0));
      matrix_set_int_element(which_rc,rnum,3, matrix_get_int_element(cols_mat,jj,1));
      rnum++;
    }
  }
  // Tidy up:
  matrix_free_int(rows_mat);
  matrix_free_int(cols_mat);
  return which_rc;
}

double
Rmatrix_get_fraction_under_c(SEXP xx, double c)
{
  // Returns the fraction of elements in the matrix that are below c
  index_t ret_val=0;
  const index_t nrow_xx=nrows(xx);
  const index_t ncol_xx=ncols(xx);
  double *x_p = REAL(xx);
  index_t i, j;
  for (i=0; i<nrow_xx; i++)
    for (j=0; j<ncol_xx; j++)
      if (x_p[i+nrow_xx*j]<c) 
	ret_val++;
  return ((double)ret_val)/(nrow_xx*ncol_xx);
}

double
Rmatrix_get_fraction_over_c(SEXP xx, double c)
{
  // Returns the fraction of elements in the matrix that are above c
  index_t ret_val = 0;
  const index_t nrow_xx = nrows(xx);
  const index_t ncol_xx = ncols(xx);
  double *x_p = REAL(xx);
  index_t i, j;
  for (i=0; i<nrow_xx; i++)
    for (j=0; j<ncol_xx; j++)
      if (x_p[i+nrow_xx*j]>c) 
	ret_val++;
  return ((double)ret_val)/(nrow_xx*ncol_xx);
}

double
Rmatrix_min(SEXP xx)
{
  // Returns the smallest elt in the matrix (SEXP) xx:
  const index_t nr = nrows(xx);
  const index_t nc = ncols(xx);

  // Check for degenerate case:
  if ((nr==0)||(nc==0))
    return R_NegInf;

  double *x_p = REAL(xx);
  double tmin = x_p[0];
  index_t i, j;
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      if (x_p[i+nr*j]<tmin)
	tmin = x_p[i+nr*j];
  return tmin;
}

double
Rmatrix_max(SEXP xx)
{
  // Returns the largest elt in the matrix (SEXP) xx:
  const index_t nr = nrows(xx);
  const index_t nc = ncols(xx);
  double *x_p = REAL(xx);
  double tmax = x_p[0];
  index_t i, j;
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      if (x_p[i+nr*j]>tmax)
	tmax = x_p[i+nr*j];
  return tmax;
}

void
initialize_KKtots_and_MMtots(Matrix * const KKtots,
                             Matrix * const MMtots,
                             Matrix * const NNtots,
                             Matrix * const KKs,
                             const index_t numrows_pt,
                             const index_t numcols_pt)
{
  // NOTE: All of the *tots matrices have rowsums in the first numrows_pt
  // columns, and colsums in the next numcols_pt columns.

  index_t ii, jj, kk;
  const index_t n=numrows(NNtots);
  double tmp_ep_count, tmp_rowsum, tmp_colsum;
#ifdef _DBG_
  const index_t nr_MMtots = numrows(MMtots);
  const index_t nr_KKtots = numrows(KKtots);
  const index_t nr_KKs = numrows(KKs);
  if (nr_MMtots!=n)
    error("Number of rows of MMtots (%u) doesn't make that of NNtots (%u)\n",nr_MMtots,n);
  if (nr_KKtots!=n)
    error("Number of rows of KKtots (%u) doesn't make that of NNtots (%u)\n",nr_KKtots,n);
  if (nr_KKs!=n)
    error("Number of rows of KKs (%u) doesn't make that of NNtots (%u)\n",nr_KKs,n);
#endif
  // Precinct-by-precinct:
  for (ii=0; ii<n; ii++){

    // Compute the row totals first:
    for (jj=0; jj<numrows_pt; jj++){

      // Reset the counter to zero, find NNs row total:
      tmp_ep_count = 0.0;
      tmp_rowsum = matrix_fast_get_element(NNtots,ii,jj,n);

      // Need to sum over the columns:
      for (kk=0; kk<numcols_pt; kk++)
        tmp_ep_count += matrix_fast_get_element(KKs,ii,(jj*numcols_pt) + kk,n);
      
      // Another error check:
      if (tmp_ep_count>tmp_rowsum){
        // Switch to R indexing for the error report:
        error("Invalid exit poll totals in precinct %u, EP row %u sum = %g, NN row %u sum = %g\n",
              ii+1,jj+1,tmp_ep_count,jj+1,tmp_rowsum);
      }

      // Put the KKtots and MMtots=NNtots-KKtots in place:
      matrix_fast_set_element(KKtots,ii,jj,n, tmp_ep_count);
      matrix_fast_set_element(MMtots,ii,jj,n, tmp_rowsum-tmp_ep_count);

    } // Finished row jj

    // Next, compute the column totals:
    for (jj=0; jj<numcols_pt; jj++){

      // Reset the counter to zero, find NNs column total:
      tmp_ep_count = 0.0;
      tmp_colsum = matrix_fast_get_element(NNtots,ii,numrows_pt+jj,n);

      // Need to sum over the rows:
      for (kk=0; kk<numrows_pt; kk++)
        tmp_ep_count += matrix_fast_get_element(KKs,ii,(kk*numcols_pt)+jj,n);
      
      // Another error check:
      if (tmp_ep_count>tmp_colsum){
        // Switch to R indexing for the error report:
        error("Invalid exit poll totals in precinct %u, EP col %u sum = %g, NN col %u sum = %g\n",
              ii+1,jj+1,tmp_ep_count,jj+1,tmp_colsum);
      }

      // Put the KKtots and MMtots=NNtots-KKtots in place:
      matrix_fast_set_element(KKtots,ii,numrows_pt+jj,n, tmp_ep_count);
      matrix_fast_set_element(MMtots,ii,numrows_pt+jj,n, tmp_colsum-tmp_ep_count);

    } // Finished col jj

  } // Finished precinct ii

#ifdef _DBG_
  Rprintf("Finished initialization of KKtots and MMtots\n");
#endif

  return;
}

void
multiply_list_of_X_by_eta(Matrix * const mu_mat_cu,
                          Matrix_p * const Xmatrix_list,
                          Matrix * const eta_vec_cu)
{
  // This function assumes all dimensions are correct, not error checked:
  const index_t n = numrows(mu_mat_cu);
  const index_t p = numcols(mu_mat_cu);
  const index_t eta_len = numcols(eta_vec_cu);
  index_t ii, jj, kk;
  double tmp = 0.0;

  // Temporary Xmatrix pointer:
  Matrix *Xmatrix = NULL;

  for (ii=0; ii<n; ii++){

    // Extract precinct ii's design matrix:
    Xmatrix = get_mat_p_ptr(Xmatrix_list,ii);

    // Now multiply Xmatrix and eta_vec_cu and store in row ii of mu_mat_cu:
    for (jj=0; jj<p; jj++){

      // Reset counter:
      tmp = 0.0;

      for (kk=0; kk<eta_len; kk++){
        tmp += matrix_fast_get_element(Xmatrix,jj,kk,p)*
               matrix_fast_get_element(eta_vec_cu,0,kk,1);
      }

      // Put the computed result into mu_mat_cu:
      matrix_fast_set_element(mu_mat_cu,ii,jj,n, tmp);

    } // end loop over columns of mu_mat_cu
  } // end loop over precincts

#ifdef _DBG_
  if (0){ // All checks out, no need for verbosity...
  for (ii=0; ii<n; ii++){
    Xmatrix = get_mat_p_ptr(Xmatrix_list,ii);
    Rprintf("Multiply X %*% eta, where:\n\n");
    Rprintf("X = \n");
    matrix_print_all(Xmatrix);
    Rprintf("eta = \n");
    matrix_print_all(eta_vec_cu);
    Rprintf("Result is:\n");
    for (jj=0; jj<p; jj++){
      Rprintf("%g\t",matrix_fast_get_element(mu_mat_cu,ii,jj,n));
    }
    Rprintf("\n#############################\n");
  }
  }
#endif

  return;
}

