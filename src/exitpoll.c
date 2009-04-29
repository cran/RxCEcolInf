#include "gqjrssa.h"

SEXP TuneWithExitPoll(SEXP args)
{
  long nProtected = 0, ii, jj, kk;
  double m12_log_SIGMA_det_cu;
  SEXP list_ret, names_ret;
  SEXP acc_THETAS_t_vec, acc_THETAS_t_matrix;
  SEXP acc_THETAS_Diri_vec, acc_THETAS_Diri_matrix;
  SEXP rho_vec_ret;
  SEXP acc_NNs_multinom_vec, vld_NNs_multinom_vec;
  SEXP acc_NNs_multinom_matrix, vld_NNs_multinom_matrix;
  const unsigned int dbg = INTEGER(getListElement(args,"dbg"))[0];

  GetRNGstate();

  //  Print a notice
  if (dbg){
    Rprintf("***********************************************\n");
    Rprintf("Currently running the TuneWithExitPoll function\n");
    Rprintf("***********************************************\n\n");
  }
  // Unpack the NNtots & NNbounds matrices:
  Matrix * const NNtots   = rmatrix_unpack_new(getListElement(args,"NNtots"));
  Matrix * const NNbounds = rmatrix_unpack_new(getListElement(args,"NNbounds"));
  Matrix * const KKs = rmatrix_unpack_new(getListElement(args,"ExitPoll"));
  Matrix * const MMbounds = rmatrix_unpack_new(getListElement(args,"MMbounds"));

  // The frequency and number of repititions for the MH proposal:
  Matrix * const sr_probs = rmatrix_unpack_new(getListElement(args,"sr_probs"));
  Matrix_int * const sr_reps  = rmatrix_int_unpack_new(getListElement(args,"sr_reps"));

  //  Unpack the scalar elements and report them
  const long num_iters      = INTEGER(getListElement(args, "num_iters"))[0];
  const index_t numrows_pt  = INTEGER(getListElement(args, "numrows_pt"))[0];
  const index_t numcols_pt  = INTEGER(getListElement(args, "numcols_pt"))[0];
  const index_t num_runs    = INTEGER(getListElement(args, "num_runs"))[0];
  const index_t num_scans   = INTEGER(getListElement(args, "numscans"))[0];
  const double psi_0        = REAL(getListElement(args, "psi_0"))[0];
  const double kappa_0      = REAL(getListElement(args, "kappa_0"))[0];
  const double nu_0         = REAL(getListElement(args, "nu_0"))[0];
  const double nolocalmode  = REAL(getListElement(args, "nolocalmode"))[0];
  const double dof          = REAL(getListElement(args, "dof"))[0];
  
  const index_t numcells_pt = numrows_pt*numcols_pt;
  const index_t p   = (numcols_pt-1)*numrows_pt;

	// Error-checking:
  const index_t n = numrows(NNtots);
  const index_t nc_NNbounds = numcols(NNbounds);
  if (numrows(NNbounds)!=n || numrows(MMbounds)!=n)
    error("Bounds and totals must be matrices of with the same number of rows");
  if (numrows(KKs)!=n)
    error("The exit poll and count matrices must have the same number of rows");
  if (numcols(KKs)!=(numrows_pt*numcols_pt))
    error("The exit poll must have R*C columns");
  if (numrows(sr_probs)!=n || numcols(sr_probs)!=numrows_pt)
    error("sr_probs must be a matrix of dimension %u x %u",n,numrows_pt);
  if (numrows_int(sr_reps)!=n || numcols_int(sr_reps)!=numrows_pt)
    error("sr_reps must be a matrix of dimension %u x %u",n,numrows_pt);

  // Compute a conservative maxrange:
  double tmp_val, maxrange=0.0;
  for (ii=0; ii<n; ii++){
    for (jj=0; jj<nc_NNbounds; jj++){
      tmp_val = matrix_fast_get_element(NNbounds,ii,jj,n);
      if (tmp_val>maxrange)
        maxrange = tmp_val;
    }
  }
  maxrange += 1.0;
  double *ff_vec = (double *) R_alloc(maxrange,sizeof(double));
  if (ff_vec==NULL)
    error("Memory allocation failure (ff_vec)");

  const index_t maxrange_it = (index_t)maxrange;
  Matrix const * const lfactorial_vector = create_log_factorial_lookup_table(maxrange_it);

  // Unpack and initialize all matrices:
  Matrix * const mu_vec_0  = rmatrix_vector_unpack_new(getListElement(args, "mu_vec_0"));
  Matrix * const mu_vec_cu = rmatrix_vector_unpack_new(getListElement(args, "mu_vec_cu"));
  Matrix * const AUG                    = rmatrix_new(p, p+1);
  Matrix * const NNprop_vec             = rmatrix_new(1, numcells_pt);
  Matrix * const NNbounds_temp_vec      = rmatrix_new(1, numcells_pt*2);
  Matrix * const NNtots_temp_vec        = rmatrix_new(1, numrows_pt+numcols_pt);
  Matrix * const SIGMA_cu               = rmatrix_new(p, p);
  Matrix * const SIGMA_chol_cu          = rmatrix_new(p, p);
  Matrix * const SIGMA_chol_cu_temp     = rmatrix_new(p, p);
  Matrix * const THETAS_count_use_t     = rmatrix_new(1, n);
  Matrix * const THETAS_count_use_Diri  = rmatrix_new(1, n);
  Matrix * const NNs_count_use_multinom = rmatrix_new(n, numrows_pt);
  Matrix * const temp1_vec              = rmatrix_new(1, p);
  Matrix * const temp2_vec              = rmatrix_new(1, p);
  Matrix * const prop_OMEGA             = rmatrix_new(1, p);
  Matrix * const prop_THETA             = rmatrix_new(1, numcells_pt);
  Matrix * const AUG_in_matrix_inverse  = rmatrix_new(p, p);
  Matrix * const tmpMean                = rmatrix_new(p, 1);
  Matrix * const tmpOut                 = rmatrix_new(p, 1);
  Matrix * const tmpScalar              = rmatrix_new(1, 1);
  Matrix * const SIGMA_inverse_for_prop = rmatrix_new(p, p);
  Matrix * const SIGMA_dims             = rmatrix_new(p, p);
  Matrix * const tmat1                  = rmatrix_new(p, p);
  Matrix * const tmat2                  = rmatrix_new(p, p);
  Matrix * const tmat3                  = rmatrix_new(p, p);
  Matrix * const tmat4                  = rmatrix_new(p, p);
  Matrix * const tvec1_dimOMs           = rmatrix_new(1, p);
  Matrix * const tvec2_dimOMs           = rmatrix_new(1, p);
  Matrix * const tvec1_dimSI            = rmatrix_new(1, p);
  Matrix * const tvec2_dimSI            = rmatrix_new(1, p);
	Matrix * const NNs_prop               = rmatrix_new(numrows_pt,numcols_pt); 
	Matrix * const MMs_prop               = rmatrix_new(numrows_pt,numcols_pt); 
	Matrix * const NNs_curr               = rmatrix_new(numrows_pt,numcols_pt);
	Matrix * const MMs_curr               = rmatrix_new(numrows_pt,numcols_pt);
	Matrix * const tmp_KKs                = rmatrix_new(numrows_pt,numcols_pt);
	Matrix * const multinomial_parameters = rmatrix_new(1,numcols_pt);
	Matrix * const prop_row               = rmatrix_new(1,numcols_pt);
	Matrix * const curr_row               = rmatrix_new(1,numcols_pt);

  Matrix * const rho_vec = rmatrix_vector_unpack_new(getListElement(args, "rho_vec"));
  Matrix * const use_Diri_every_vec = rmatrix_vector_unpack_new(getListElement(args, "use_Diri_every_vec"));
  Matrix_int * const which_rc_int = get_which_rc(numrows_pt,numcols_pt);
  Matrix_int * const ordervec_int       = rmatrix_new_int(1, numrows_int(which_rc_int));
  Matrix_int * const zerotonumperms_int = rmatrix_new_int(1, numcols_int(ordervec_int));
  Matrix * const PSI_0 = rmatrix_diag_new(p, p, psi_0);
  
  // Allocate acceptance structures directly as R-objects to avoid wasting memory:
  PROTECT(acc_THETAS_t_vec     = allocVector(REALSXP,n));
  PROTECT(acc_THETAS_Diri_vec  = allocVector(REALSXP,n));
  PROTECT(acc_THETAS_t_matrix     = allocMatrix(REALSXP,n,num_runs));
  PROTECT(acc_THETAS_Diri_matrix  = allocMatrix(REALSXP,n,num_runs));
  nProtected += 4;

  double *acc_tv_p = REAL(acc_THETAS_t_vec);
  double *acc_dv_p = REAL(acc_THETAS_Diri_vec);
  double *acc_tm_p = REAL(acc_THETAS_t_matrix);
  double *acc_dm_p = REAL(acc_THETAS_Diri_matrix);

  memset(acc_tv_p,0,sizeof(double)*n);
  memset(acc_dv_p,0,sizeof(double)*n);
  memset(acc_tm_p,0,sizeof(double)*n*num_runs);
  memset(acc_dm_p,0,sizeof(double)*n*num_runs);

  // Store multinomial stats in a list of num_runs matrices:
  SEXP vld_NNs_multinom_list, acc_NNs_multinom_list;
  SEXP vld_NNs_multinom_mat, acc_NNs_multinom_mat;

  PROTECT(vld_NNs_multinom_list   = allocVector(VECSXP,num_runs));
  PROTECT(acc_NNs_multinom_list   = allocVector(VECSXP,num_runs));
  nProtected += 2;

	double *vld_NNs_multinomial_mat_p_holder[num_runs];
	double *acc_NNs_multinomial_mat_p_holder[num_runs];
	double *vld_mv_p;
	double *acc_mv_p;
  for (ii=0; ii<num_runs; ii++){
    PROTECT(vld_NNs_multinom_mat = allocMatrix(REALSXP,n,numrows_pt));
    PROTECT(acc_NNs_multinom_mat = allocMatrix(REALSXP,n,numrows_pt));
		SET_VECTOR_ELT(vld_NNs_multinom_list,ii, vld_NNs_multinom_mat);
		SET_VECTOR_ELT(acc_NNs_multinom_list,ii, acc_NNs_multinom_mat);
		vld_mv_p = REAL(vld_NNs_multinom_mat);
		acc_mv_p = REAL(acc_NNs_multinom_mat);
  	memset(vld_mv_p,0,sizeof(double)*n*numrows_pt);
  	memset(acc_mv_p,0,sizeof(double)*n*numrows_pt);
    vld_NNs_multinomial_mat_p_holder[ii] = vld_mv_p;
    acc_NNs_multinomial_mat_p_holder[ii] = acc_mv_p;
    nProtected += 2;
  }

  //  Initialize, set up a starting NNs matrix
  Matrix * const NNs    = rmatrix_new(n, numcells_pt);
  Matrix * const MMs = rmatrix_new(n, numcells_pt);
  Matrix * const THETAS = rmatrix_new(n, numcells_pt);
  Matrix * const OMEGAS = rmatrix_new(n, p);

  // Initialize and compute the KKtots and MMtots matrices,
  // Want to wait until the EP has been checked before doing this.
  Matrix * const KKtots = rmatrix_new(n,numrows_pt+numcols_pt);
  Matrix * const MMtots = rmatrix_new(n,numrows_pt+numcols_pt);
  initialize_KKtots_and_MMtots(KKtots,MMtots,NNtots,KKs,numrows_pt,numcols_pt);

  // Generate random starting state:
  draw_NNs_MMs_indep_start(NNs, MMs, KKs, NNbounds, MMbounds, NNtots, MMtots, 
                           NNbounds_temp_vec, NNtots_temp_vec, numrows_pt, 
                           numcols_pt, numcells_pt);

  // Check the starting states are okay:
  check_ep_validity(NNs,MMs,KKs,numcells_pt,numrows_pt,numcols_pt);
  check_bounds_all(NNs,NNbounds,numcells_pt);

  draw_THETAS_from_NNs_start(THETAS, NNs, NNtots, numrows_pt, numcols_pt);
  THETAS_to_OMEGAS(THETAS, OMEGAS, numrows_pt, numcols_pt);
  
  if (dbg){
    Rprintf("\nYou passed NNtots, the first two rows of which are as follows.\n");
    matrix_print_subset(NNtots, 0, 1, 0, numcols(NNtots)-1);
    Rprintf("\n\nYou passed NNbounds, the first two rows of which are as follows.\n");
    matrix_print_subset(NNbounds, 0, 1, 0, numcols(NNbounds)-1);
    Rprintf("You passed the following scalar parameters:\n");
    Rprintf("psi_0 = %f, nu_0 = %f, kappa_0 = %f.\n", psi_0, nu_0, kappa_0);
    Rprintf("numrows_pt = %d, numcols_pt = %d, num_iters = %d.\n",
   	    numrows_pt, numcols_pt, num_iters);
    Rprintf("nolocalmode = %f, dof = %f.\n", nolocalmode, dof);
    Rprintf("maxrange = %f, num_runs = %d, num_scans = %d.\n", maxrange, num_runs, num_scans);
    Rprintf("\n\nCreated a starting NNs matrix, the first two rows of which are as follows.\n");
    matrix_print_subset(NNs, 0, 1, 0, numcols(NNs)-1);
    Rprintf("\n\nYou passed the following mu prior parameter vector.\n");
    matrix_print_all(mu_vec_0);
    Rprintf("\n\nYou passed the following mu starting vector.\n");
    matrix_print_all(mu_vec_cu);
    Rprintf("\n\nCreated THETAS, the first two rows of which are as follows.\n");
    matrix_print_subset(THETAS, 0, 1, 0, numcols(THETAS)-1);
    Rprintf("\n\nCreated an OMEGAS matrix, the first two rows of which are as follows.\n");
    matrix_print_subset(OMEGAS, 0, 1, 0, numcols(OMEGAS)-1);
    Rprintf("\n\nUnpacked rho_vec; its first two values are as follows.\n");
    Rprintf("%f\t%f\n", matrix_get_element(rho_vec, 0, 0),
  	    matrix_get_element(rho_vec, 0, 0));
    Rprintf("\n\nUnpacked use_Diri_every_vec; its first two values are as follows.\n");
    Rprintf("%f\t%f\n", matrix_get_element(use_Diri_every_vec, 0, 0),
	    matrix_get_element(use_Diri_every_vec, 0, 0));
    Rprintf("Created the which_rc matrix.\n\n");
    Rprintf("The matrix of possible small Gibbs combinations is as follows.\n");
    matrix_print_all_int(which_rc_int);
    Rprintf("\n\nCreated the matrix PSI_0 as follows:\n");
    matrix_print_all(PSI_0);
    Rprintf("Created the temporary matrices.\n\n");
  }

  //  Error check on the starting NNs
  check_bounds_all(NNs, NNbounds, numcells_pt);

  //  RUN THE BIG LOOP
  for (kk=0; kk<num_runs; kk++){

		// Set the multinomial statistic pointers:
    vld_mv_p = vld_NNs_multinomial_mat_p_holder[kk];
    acc_mv_p = acc_NNs_multinomial_mat_p_holder[kk];
		vld_NNs_multinom_mat = VECTOR_ELT(vld_NNs_multinom_list,kk);
		acc_NNs_multinom_mat = VECTOR_ELT(acc_NNs_multinom_list,kk);

    for (ii=0; ii<num_iters; ii++){

      draw_SIGMA(SIGMA_cu, nu_0, PSI_0, mu_vec_cu, OMEGAS, 
		 tmat1, tmat2, tmat3, tmat4, AUG_in_matrix_inverse,dbg);

      matrix_cholesky(SIGMA_cu, SIGMA_chol_cu);
      m12_log_SIGMA_det_cu = -0.5*matrix_determinant(SIGMA_cu,tmat1,1);
      sample_equ_pr_wo_replace_int(zerotonumperms_int, ordervec_int);

      draw_NNs_MMs(NNs, MMs, KKs, NNprop_vec, NNbounds, ii, NNbounds_temp_vec, 
           NNtots, MMtots, KKtots, ff_vec,
	         NNtots_temp_vec, THETAS, ordervec_int, which_rc_int, nolocalmode, 
	         numrows_pt, numcols_pt, numcells_pt, num_scans,
					 NNs_prop, MMs_prop, NNs_curr, MMs_curr, tmp_KKs,
           multinomial_parameters, curr_row, prop_row, sr_probs, sr_reps,
					 vld_mv_p, acc_mv_p, NNs_count_use_multinom,
           lfactorial_vector);
      
      draw_THETAS_t_and_Dirichlet(THETAS, OMEGAS,
				  prop_THETA, prop_OMEGA,
				  SIGMA_chol_cu, temp1_vec,
				  temp2_vec, NNs,
				  mu_vec_cu, SIGMA_cu,
				  AUG, acc_tv_p,
				  rho_vec, SIGMA_chol_cu_temp,
				  acc_dv_p,
				  use_Diri_every_vec,
				  THETAS_count_use_t, THETAS_count_use_Diri,
				  dof, numrows_pt, numcols_pt,
				  numcells_pt, ii,
					SIGMA_inverse_for_prop,
					tmpMean, tmpOut, tmpScalar, SIGMA_dims);

      draw_mu(mu_vec_cu, kappa_0, mu_vec_0, OMEGAS, SIGMA_chol_cu,
	      tvec1_dimOMs, tvec2_dimOMs, tmat1, tvec1_dimSI, tvec2_dimSI);

    //  Error check
    check_bounds_all(NNs, NNbounds, numcells_pt);
    check_ep_validity(NNs,MMs,KKs,numcells_pt,numrows_pt,numcols_pt);

#ifdef _DBG_
  if (dbg){
    if ((ii%1000) == 0)
      Rprintf("Iteration %d completed\n",ii);
  }
#endif
    } // END ITERATION

    //  Turn counts into fractions
    adjust_acc_vector(acc_THETAS_t_vec, THETAS_count_use_t);
    adjust_acc_vector(acc_THETAS_Diri_vec, THETAS_count_use_Diri);
    adjust_acc_matrix(vld_NNs_multinom_mat, NNs_count_use_multinom);
    adjust_acc_matrix(acc_NNs_multinom_mat, NNs_count_use_multinom);

    //  Print some reports to the screen
    report_dry_run(acc_THETAS_t_vec, acc_THETAS_Diri_vec, kk);

    //  Save the acceptance fractions for later study (and passing back to R)
    for (ii=0; ii<n; ii++){
      acc_tm_p[ii+n*kk] = acc_tv_p[ii];
      acc_dm_p[ii+n*kk] = acc_dv_p[ii];
    }

    //  Adjust the rho_vec for the next go-around
    adjust_rho_vec(rho_vec, acc_THETAS_t_vec);

    //  Reinitialize the acc and count vectors:
    memset(acc_tv_p,0,sizeof(double)*n);
    memset(acc_dv_p,0,sizeof(double)*n);
    matrix_fastset(THETAS_count_use_t, 0);
    matrix_fastset(THETAS_count_use_Diri, 0);
    matrix_fastset(NNs_count_use_multinom, 0);
  }

  //  Package rho_vec to return to R
  rho_vec_ret = matrix_vector_repack_new(rho_vec);
  nProtected += 1;

  //  Create the final return list
  PROTECT(list_ret = allocVector(VECSXP, 5));
  ++nProtected;
  SET_VECTOR_ELT(list_ret, 0, rho_vec_ret);
  SET_VECTOR_ELT(list_ret, 1, acc_THETAS_t_matrix);
  SET_VECTOR_ELT(list_ret, 2, acc_THETAS_Diri_matrix);
  SET_VECTOR_ELT(list_ret, 3, vld_NNs_multinom_list);
  SET_VECTOR_ELT(list_ret, 4, acc_NNs_multinom_list);

  //  Assign names for the return list
  PROTECT(names_ret = allocVector(STRSXP, 5));
  ++nProtected;
  SET_STRING_ELT(names_ret, 0, mkChar("rhos"));
  SET_STRING_ELT(names_ret, 1, mkChar("acc.t"));
  SET_STRING_ELT(names_ret, 2, mkChar("acc.Diri"));
  SET_STRING_ELT(names_ret, 3, mkChar("vld.NNs"));
  SET_STRING_ELT(names_ret, 4, mkChar("acc.NNs"));

  //  Marry the names to the return list
  setAttrib(list_ret, R_NamesSymbol, names_ret);

  PutRNGstate();

  if (dbg){
    //Rprintf("\n\nReached the end of the TuneWithExitPoll function without crashing!\n\n");
    Rprintf("\n");
  }
  UNPROTECT(nProtected);
  return list_ret;
}

SEXP AnalyzeWithExitPoll(SEXP args)
{
  index_t kk;
  long nProtected = 0, ii, jj, rownum_store = 0;
  long nextplace, listlength, iter_count = 1;
  long colnum_store_NNs_internals = 0, colnum_store_THETAS= 0;
  double m12_log_SIGMA_det_cu;
  double save_every_NN_internals = R_PosInf;
  double save_every_THETAS = R_PosInf;
  SEXP THETAS_save_ret = R_NilValue;
  SEXP NNs_save_ret = R_NilValue;
  //const 
  unsigned int dbg = INTEGER(getListElement(args,"dbg"))[0];

  GetRNGstate();

  //  Print a notice
  if (dbg){
    Rprintf("**************************************************\n");
    Rprintf("Currently running the AnalyzeWithExitPoll function\n");
    Rprintf("**************************************************\n\n");
  }

  // R Constraints:
  // (1) NNs and NNbounds lists must be same length AND error checked.
  // (2) Target distribution must be first element of list, decreasing counts thereafter.
  // (3) Additional "prob_re" argument must be supplied.

  // Unpack the lists of NNtots & NNbounds matrices:
  Rprintf("Unpacking the matrix of row and column totals...\n");
  Matrix * const NNtots = rmatrix_unpack_new(getListElement(args, "NNtots"));
  Rprintf("done. Unpacking the exit poll matrix...\n");
  Matrix * const KKs = rmatrix_unpack_new(getListElement(args, "ExitPoll"));
  Rprintf("done. Unpacking the matrices of bounds...\n");
  Matrix * const NNbounds = rmatrix_unpack_new(getListElement(args, "NNbounds"));
  Matrix * const MMbounds = rmatrix_unpack_new(getListElement(args, "MMbounds"));

  Rprintf("done. Unpacking eschew option...\n");
  const int eschew = (int)INTEGER(getListElement(args,"eschew"))[0];

  // The number of precincts and (R+C):
  const index_t n = numrows(NNtots);
  const index_t nc_NNtots = numcols(NNtots); // R+C

  Rprintf("done. Setting up sr_reps and sr_probs...\n");
  
	// The frequency and number of repititions for the MH proposal:
	Matrix * const sr_probs     = rmatrix_unpack_new(getListElement(args, "sr_probs"));
	Matrix_int * const sr_reps  = rmatrix_int_unpack_new(getListElement(args, "sr_reps"));

  //  Unpack the scalar elements and report them (error-checked in R code):
  Rprintf("done. Unpacking scalar arguments...\n");
  const index_t keep_restart_info = (index_t)INTEGER(getListElement(args, "keep_restart_info"))[0];
  const index_t rstart_NNs_MMs = (index_t)INTEGER(getListElement(args, "rstart_NNs_MMs"))[0];
  const index_t rstart_THETAS=(index_t)INTEGER(getListElement(args, "rstart_THETAS"))[0];
  const index_t print_every = (index_t)INTEGER(getListElement(args, "print_every"))[0];
  const long num_iters      = (long)INTEGER(getListElement(args, "num_iters"))[0];
  const index_t numrows_pt  = (index_t)INTEGER(getListElement(args, "numrows_pt"))[0];
  const index_t numcols_pt  = (index_t)INTEGER(getListElement(args, "numcols_pt"))[0];
  const index_t num_scans   = (index_t)INTEGER(getListElement(args, "numscans"))[0];

	// Construct useful quantities involving the table dimensions:
  const index_t numcells_pt = numrows_pt*numcols_pt;
  const index_t p = (numcols_pt-1)*numrows_pt;
  Rprintf("done. Error checks...\n");

	// Error check the bounds matrix:
  const index_t nc_NNbounds = numcols(NNbounds);
  if (numrows(NNbounds)!=n || numrows(MMbounds)!=n)
    error("Bounds and totals must be matrices of with the same number of rows");
  if (numrows(KKs)!=n)
    error("The exit poll and count matrices must have the same number of rows");
  if (numcols(KKs)!=(numrows_pt*numcols_pt))
    error("The exit poll must have R*C columns");

	// Error checks on the sr_reps and sr_probs matrices:
  if (numrows(sr_probs)!=n || numcols(sr_probs)!=numrows_pt)
    error("sr_probs must be a matrix of dimension %u x %u",n,numrows_pt);
  if (numrows_int(sr_reps)!=n || numcols_int(sr_reps)!=numrows_pt)
    error("sr_reps must be a matrix of dimension %u x %u",n,numrows_pt);

	Rprintf("done. Computing maxrange...\n");
  // Compute a conservative maximum range of precinct level cell counts.
	// (this is used in the nchg draws)

  // CHANGE TO NNtots
  double tmp_val, maxrange=0.0;
  for (ii=0; ii<n; ii++){
    for (jj=0; jj<nc_NNbounds; jj++){
      tmp_val = matrix_fast_get_element(NNbounds,ii,jj,n);
      if (tmp_val>maxrange)
        maxrange = tmp_val;
    }
  }
  maxrange += 1.0;
  Rprintf("done. Extracting remaining arguments...\n");

	// Remaining arguments:
  const double how_many_NN_internals = REAL(getListElement(args, "keepNNinternals"))[0];
  const double how_many_THETAS       = REAL(getListElement(args, "keepTHETAS"))[0];

  const double psi_0       = REAL(getListElement(args, "psi_0"))[0];
  const double kappa_0     = REAL(getListElement(args, "kappa_0"))[0];
  const double nu_0        = REAL(getListElement(args, "nu_0"))[0];
  const double nolocalmode = REAL(getListElement(args, "nolocalmode"))[0];
  const double dof         = REAL(getListElement(args, "dof"))[0];
  const double save_every  = REAL(getListElement(args, "save_every"))[0];

  // Few quick error checks:
  if (psi_0<=0)
    error("psi_0 must be positive");
  if (kappa_0<=0)
    error("kappa_0 must be positive");
  if (nu_0<=0)
    error("nu_0 must be positive");
  if (dof<=0)
    error("dof must be positive");
  if (save_every<=0)
    error("save_every must be positive");

	// Allocate ff_vec for efficiency purposes (used in nchg draws):
  double * const ff_vec = (double *) R_alloc(maxrange,sizeof(double));
  if (ff_vec==NULL)
    error("Memory allocation failure (ff_vec)");
  
  //  Error check the number of iterations between saved draws:
  if (fmod(num_iters, save_every) != 0.0)
    error("num.iters not a whole multiple of save.every");

	// Compute the number of rows needed to store the draws:
  const index_t num_rows_save_matrices = (index_t)(num_iters/save_every);

  Rprintf("done. Initializing factorial lookup table...\n");

  const index_t maxrange_it = (index_t)maxrange;
  Matrix const * const lfactorial_vector = create_log_factorial_lookup_table(maxrange_it);
  if (dbg){
    Rprintf("Created a factorial lookup table from 0 to %u:\n",maxrange_it);
    Rprintf("The first few elements of the factorial lookup table are:\n");
    matrix_print_subset(lfactorial_vector,0,0,0,maxrange_it<10?maxrange_it:9);
    Rprintf("The final few elements of the factorial lookup table are:\n");
    matrix_print_subset(lfactorial_vector,0,0,maxrange_it<10?0:(maxrange_it-9),maxrange_it);
  }

  Rprintf("done. Initializing sampler...\n");

  // Unpack and initialize all matrices:
  Matrix * const mu_vec_0  = rmatrix_vector_unpack_new(getListElement(args, "mu_vec_0"));
  Matrix * const use_Diri_every_vec    = rmatrix_vector_unpack_new(getListElement(args, "use_Diri_every_vec"));
  Matrix_int * const which_rc_int      = get_which_rc(numrows_pt,numcols_pt);
  Matrix * const PSI_0                 = rmatrix_diag_new(p, p, psi_0);
  Matrix * const rho_vec               = rmatrix_vector_unpack_new(getListElement(args, "rho_vec"));
  Matrix * const THETAS_count_use_t    = rmatrix_new(1, n); // (doesn't make as much sense as a counter)
  Matrix * const THETAS_count_use_Diri = rmatrix_new(1, n); // (doesn't make as much sense as a counter)
  Matrix * const NNs_count_use_multinom= rmatrix_new(n, numrows_pt); // (doesn't make as much sense as a counter)
  Rprintf("(finished non-varying objects)\n");  

  // Objects thats are reset every iteration:
  Matrix * const AUG                   = rmatrix_new(p, p+1); 
  Matrix * const NNprop_vec            = rmatrix_new(1, numcells_pt);   
  Matrix * const NNbounds_temp_vec     = rmatrix_new(1, numcells_pt*2); 
  Matrix * const NNtots_temp_vec       = rmatrix_new(1, numrows_pt + numcols_pt); 
  Matrix * const SIGMA_chol_cu_temp    = rmatrix_new(p, p); 
  Matrix * const temp1_vec             = rmatrix_new(1, p); 
  Matrix * const temp2_vec             = rmatrix_new(1, p); 
  Matrix * const prop_OMEGA            = rmatrix_new(1, p); 
  Matrix * const prop_THETA            = rmatrix_new(1, numcells_pt); 
  Matrix * const AUG_in_matrix_inverse = rmatrix_new(p, p); 
  Matrix * const tmpMean               = rmatrix_new(p, 1); 
  Matrix * const tmpOut                = rmatrix_new(p, 1); 
  Matrix * const tmpScalar             = rmatrix_new(1, 1); 
  Matrix * const SIGMA_inverse_for_prop= rmatrix_new(p, p); 
  Matrix * const tmat1                 = rmatrix_new(p, p); 
  Matrix * const tmat2                 = rmatrix_new(p, p); 
  Matrix * const tmat3                 = rmatrix_new(p, p); 
  Matrix * const tmat4                 = rmatrix_new(p, p); 
  Matrix * const tvec1_dimOMs          = rmatrix_new(1, p); 
  Matrix * const tvec2_dimOMs          = rmatrix_new(1, p); 
  Matrix * const tvec1_dimSI           = rmatrix_new(1, p); 
  Matrix * const tvec2_dimSI           = rmatrix_new(1, p); 
  Matrix * const SIGMA_dims            = rmatrix_new(p, p); 
	Matrix * const NNs_prop              = rmatrix_new(numrows_pt,numcols_pt); 
	Matrix * const MMs_prop              = rmatrix_new(numrows_pt,numcols_pt); 
	Matrix * const NNs_curr              = rmatrix_new(numrows_pt,numcols_pt);
	Matrix * const MMs_curr              = rmatrix_new(numrows_pt,numcols_pt);
	Matrix * const tmp_KKs               = rmatrix_new(numrows_pt,numcols_pt);
	Matrix * const multinomial_parameters= rmatrix_new(1,numcols_pt); 
	Matrix * const prop_row              = rmatrix_new(1,numcols_pt); 
	Matrix * const curr_row              = rmatrix_new(1,numcols_pt); 
  Matrix_int * const ordervec_int      = rmatrix_new_int(1, numrows_int(which_rc_int)); 
  Matrix_int * const zerotonumperms_int= rmatrix_new_int(1, numcols_int(ordervec_int)); 
  Rprintf("(finished computational objects)\n");  

  SEXP checker, tptr;
  checker = getListElement(args, "mu_vec_cu");

  // Remaining objects:
  Matrix * const MMs           = rmatrix_new(n, numcells_pt);
  Matrix * const NNs           = rmatrix_new(n, numcells_pt);
  Matrix * const THETAS        = rmatrix_new(n, numcells_pt);
  Matrix * const OMEGAS        = rmatrix_new(n, p);
  Matrix * const SIGMA_cu      = rmatrix_new(p, p);
  Matrix * const SIGMA_chol_cu = rmatrix_new(p, p);
  Matrix * const mu_vec_cu     = rmatrix_vector_unpack_new(checker);

  // Initialize and compute the KKtots and MMtots matrices,
  // Want to wait until the EP has been checked before doing this.
  Matrix * const KKtots = rmatrix_new(n,nc_NNtots);
  Matrix * const MMtots = rmatrix_new(n,nc_NNtots);
  initialize_KKtots_and_MMtots(KKtots,MMtots,NNtots,KKs,numrows_pt,numcols_pt);

  Rprintf("(finished initialization of KKtots,MMtots)\n");  

  // Were the starting states of the NNs's or MMs's pre-specified?
  if (!rstart_NNs_MMs){

    // Retrieve starting states:
    double *start_NNs_p, *start_MMs_p;

      // First, extract the NNs starting state:
      checker = getListElement(args,"NNs_start");
      if ((nrows(checker)!=n) || (ncols(checker)!=numcells_pt))
        error("Invalid dimension of NNs starting state");
      start_NNs_p = REAL(checker);

      // Write it to the NNs matrix:
      for (ii=0; ii<n; ii++)
        for (jj=0; jj<numcells_pt; jj++)
          matrix_fast_set_element(NNs,ii,jj,n, start_NNs_p[ii+n*jj]);

      // Second, extract the MMs starting state:
      checker = getListElement(args,"MMs_start");
      if ((nrows(checker)!=n) || (ncols(checker)!=numcells_pt))
        error("Invalid dimension of NNs starting state");
      start_MMs_p = REAL(checker);

      // Write it to the NNs matrix:
      for (ii=0; ii<n; ii++)
        for (jj=0; jj<numcells_pt; jj++)
          matrix_fast_set_element(MMs,ii,jj,n, start_MMs_p[ii+n*jj]);

  } else {

      // Randomly start the NNs and MMs:

      // Generate random starting state:
      draw_NNs_MMs_indep_start(NNs, MMs, KKs, NNbounds, MMbounds, NNtots, MMtots, 
               NNbounds_temp_vec, NNtots_temp_vec, numrows_pt, numcols_pt, numcells_pt);

      // Check the starting states are okay:
      check_ep_validity(NNs,MMs,KKs,numcells_pt,numrows_pt,numcols_pt);
      check_bounds_all(NNs,NNbounds,numcells_pt);

  } 

  Rprintf("(finished starting state for MMs, NNs)\n");  

  // Were the starting state of the THETAS's (and OMEGAS's) pre-specified?
  if (!rstart_THETAS){

    // Retrieve starting states:
    double * start_THETAS_p;
      
    checker = getListElement(args,"THETAS_start");
    if ((nrows(checker)!=n) || (ncols(checker)!=numcells_pt))
      error("Invalid dimension of THETAS starting state");
    start_THETAS_p = REAL(checker);

    // Write to the THETAS and OMEGAS matrices:
    for (ii=0; ii<n; ii++)
      for (jj=0; jj<numcells_pt; jj++)
        matrix_fast_set_element(THETAS,ii,jj,n, start_THETAS_p[ii+n*jj]);

    // Do the deterministic transformation to OMEGAS:
    THETAS_to_OMEGAS(THETAS,OMEGAS,numrows_pt,numcols_pt);
    
  } else {

    // Randomly start the THETAS:
    draw_THETAS_from_NNs_start(THETAS, NNs, NNtots, numrows_pt, numcols_pt);
    THETAS_to_OMEGAS(THETAS, OMEGAS, numrows_pt, numcols_pt);
  
  }
  
  Rprintf("(finished starting state for THETAS, OMEGAS)\n");  

  // Acceptance vectors:
  SEXP acc_THETAS_t_vec_ret, acc_THETAS_Diri_vec_ret;

  PROTECT(acc_THETAS_t_vec_ret     = allocVector(REALSXP,n));
  PROTECT(acc_THETAS_Diri_vec_ret  = allocVector(REALSXP,n));
  nProtected+=2;

  double *acc_tv_p = REAL(acc_THETAS_t_vec_ret);
  double *acc_dv_p = REAL(acc_THETAS_Diri_vec_ret);

  memset(acc_tv_p,0,sizeof(double)*n);
  memset(acc_dv_p,0,sizeof(double)*n);
  
  // Store multinomial stats in a matrix:
	SEXP vld_NNs_multinom_mat_ret, acc_NNs_multinom_mat_ret;

  PROTECT(vld_NNs_multinom_mat_ret = allocMatrix(REALSXP,n,numrows_pt));
  PROTECT(acc_NNs_multinom_mat_ret = allocMatrix(REALSXP,n,numrows_pt));
  nProtected += 2;

	double *vld_mv_p = REAL(vld_NNs_multinom_mat_ret);
	double *acc_mv_p = REAL(acc_NNs_multinom_mat_ret);

  memset(vld_mv_p,0,sizeof(double)*n*numrows_pt);
  memset(acc_mv_p,0,sizeof(double)*n*numrows_pt);

  if (dbg)
    Rprintf("Allocating storage for draws...\n");

  int draw_list_length=0;
  SEXP mu_draws, SIGMA_draws, NNs_draws, LAMBDA_draws, TURNOUT_draws, GAMMA_draws, BETA_draws;
  PROTECT(mu_draws     =allocMatrix(REALSXP,num_rows_save_matrices,p));
  PROTECT(SIGMA_draws  =allocMatrix(REALSXP,num_rows_save_matrices,p*(p+1)/2));
  PROTECT(NNs_draws    =allocMatrix(REALSXP,num_rows_save_matrices,numrows_pt*numcols_pt));

  // Need to handle eschew argument here:
  int ncol_LAMBDA, ncol_BETA;
  if (!eschew){
    ncol_LAMBDA = (numcols_pt-1);
    ncol_BETA   = (numcols_pt-1);
  } else {
    ncol_LAMBDA = (numcols_pt-2);
    ncol_BETA   = (numcols_pt-2);
  }

  // Create the possibly eschew-modified matrices:
  PROTECT(LAMBDA_draws =allocMatrix(REALSXP,num_rows_save_matrices,numrows_pt*ncol_LAMBDA));
  PROTECT(TURNOUT_draws=allocMatrix(REALSXP,num_rows_save_matrices,numrows_pt));
  PROTECT(GAMMA_draws  =allocMatrix(REALSXP,num_rows_save_matrices,numrows_pt));
  PROTECT(BETA_draws   =allocMatrix(REALSXP,num_rows_save_matrices,numrows_pt*ncol_BETA));
  nProtected+=7;
  draw_list_length+=7;
  
  // R doesn't initialize to zero, so we need to :)
  memset(REAL(mu_draws),     0,sizeof(double)*num_rows_save_matrices*p);
  memset(REAL(SIGMA_draws),  0,sizeof(double)*num_rows_save_matrices*p*(p+1)/2);
  memset(REAL(NNs_draws),    0,sizeof(double)*num_rows_save_matrices*numrows_pt*numcols_pt);
  memset(REAL(LAMBDA_draws), 0,sizeof(double)*num_rows_save_matrices*numrows_pt*ncol_LAMBDA);
  memset(REAL(TURNOUT_draws),0,sizeof(double)*num_rows_save_matrices*numrows_pt);
  memset(REAL(GAMMA_draws),  0,sizeof(double)*num_rows_save_matrices*numrows_pt);
  memset(REAL(BETA_draws),   0,sizeof(double)*num_rows_save_matrices*numrows_pt*ncol_BETA);

  char tmpname[50];
  SEXP mu_colnames, mu_dimnames;
  PROTECT(mu_dimnames=allocVector(VECSXP,2));
  PROTECT(mu_colnames=allocVector(STRSXP,p));
  nProtected+=2;
  for (ii=0; ii<p; ii++){
    sprintf(tmpname,"mu_%u",(index_t)ii+1);
    SET_STRING_ELT(mu_colnames,ii, mkChar(tmpname));
  }
  SET_VECTOR_ELT(mu_dimnames,1,mu_colnames);
  dimnamesgets(mu_draws,mu_dimnames);

  SEXP SIGMA_colnames, SIGMA_dimnames;
  PROTECT(SIGMA_dimnames=allocVector(VECSXP,2));
  PROTECT(SIGMA_colnames=allocVector(STRSXP,p*(p+1)/2));
  nProtected+=2;
  index_t sofar=0;
  // Diagonals first:
  for (ii=0; ii<p; ii++){
    sprintf(tmpname,"sd_%u",(index_t)ii+1);
    SET_STRING_ELT(SIGMA_colnames,sofar, mkChar(tmpname));
    sofar++;
  }
  // Off-diagonals next:
  for (ii=0; ii<p; ii++){
    for (jj=(ii+1); jj<p; jj++){
      sprintf(tmpname,"corr_%u_%u",(index_t)ii+1,(index_t)jj+1);
      SET_STRING_ELT(SIGMA_colnames,sofar, mkChar(tmpname));
      sofar++;
    }
  }
  SET_VECTOR_ELT(SIGMA_dimnames,1,SIGMA_colnames);
  dimnamesgets(SIGMA_draws,SIGMA_dimnames);
  
  SEXP NNs_colnames, NNs_dimnames;
  PROTECT(NNs_dimnames=allocVector(VECSXP,2));
  PROTECT(NNs_colnames=allocVector(STRSXP,numrows_pt*numcols_pt));
  nProtected+=2;
  for (ii=0; ii<numrows_pt; ii++){
    for (jj=0; jj<numcols_pt; jj++){
      sprintf(tmpname,"NNr%uc%u",(index_t)ii+1,(index_t)jj+1);
      SET_STRING_ELT(NNs_colnames,jj+numcols_pt*ii, mkChar(tmpname));
    }
  }
  SET_VECTOR_ELT(NNs_dimnames,1,NNs_colnames);
  dimnamesgets(NNs_draws,NNs_dimnames);

  // Create and label the LAMBDAS: R*ncol_BETA of them in byrow=TRUE format:
  SEXP LAMBDA_colnames, LAMBDA_dimnames;
  PROTECT(LAMBDA_dimnames=allocVector(VECSXP,2));
  PROTECT(LAMBDA_colnames=allocVector(STRSXP,numrows_pt*ncol_LAMBDA));
  nProtected+=2;
  for (ii=0; ii<numrows_pt; ii++){
    for (jj=0; jj<ncol_LAMBDA; jj++){
      sprintf(tmpname,"LAMBDA_r%uc%u",(index_t)ii+1,(index_t)jj+1);
      SET_STRING_ELT(LAMBDA_colnames,jj+(numcols_pt-1)*ii, mkChar(tmpname));
    }
  }
  SET_VECTOR_ELT(LAMBDA_dimnames,1,LAMBDA_colnames);
  dimnamesgets(LAMBDA_draws,LAMBDA_dimnames);

  SEXP TURNOUT_colnames, TURNOUT_dimnames;
  PROTECT(TURNOUT_dimnames=allocVector(VECSXP,2));
  PROTECT(TURNOUT_colnames=allocVector(STRSXP,numrows_pt));
  nProtected+=2;
  for (ii=0; ii<numrows_pt; ii++){
    sprintf(tmpname,"TURNOUT_r%u",(index_t)ii+1);
    SET_STRING_ELT(TURNOUT_colnames,ii, mkChar(tmpname));
  }
  SET_VECTOR_ELT(TURNOUT_dimnames,1,TURNOUT_colnames);
  dimnamesgets(TURNOUT_draws,TURNOUT_dimnames);

  SEXP GAMMA_colnames, GAMMA_dimnames;
  PROTECT(GAMMA_dimnames=allocVector(VECSXP,2));
  PROTECT(GAMMA_colnames=allocVector(STRSXP,numrows_pt));
  nProtected+=2;
  for (ii=0; ii<numrows_pt; ii++){
    sprintf(tmpname,"GAMMA_r%u",(index_t)ii+1);
    SET_STRING_ELT(GAMMA_colnames,ii, mkChar(tmpname));
  }
  SET_VECTOR_ELT(GAMMA_dimnames,1,GAMMA_colnames);
  dimnamesgets(GAMMA_draws,GAMMA_dimnames);
  
	// Create and label the BETAS: R*ncol_BETA of them in byrow=TRUE format:
  SEXP BETA_colnames, BETA_dimnames;
  PROTECT(BETA_dimnames=allocVector(VECSXP,2));
  PROTECT(BETA_colnames=allocVector(STRSXP,numrows_pt*ncol_BETA));
  nProtected+=2;
  for (ii=0; ii<numrows_pt; ii++){
    for (jj=0; jj<ncol_BETA; jj++){
      sprintf(tmpname,"BETA_r%uc%u",(index_t)ii+1,(index_t)jj+1);
      SET_STRING_ELT(BETA_colnames,jj+(numcols_pt-1)*ii, mkChar(tmpname));
    }
  }
  SET_VECTOR_ELT(BETA_dimnames,1,BETA_colnames);
  dimnamesgets(BETA_draws,BETA_dimnames);

  SEXP R_draw_list;
  PROTECT(R_draw_list=allocVector(VECSXP,draw_list_length));
  nProtected++;

  SET_VECTOR_ELT(R_draw_list,0,mu_draws);
  SET_VECTOR_ELT(R_draw_list,1,SIGMA_draws);
  SET_VECTOR_ELT(R_draw_list,2,NNs_draws);
  SET_VECTOR_ELT(R_draw_list,3,LAMBDA_draws);
  SET_VECTOR_ELT(R_draw_list,4,TURNOUT_draws);
  SET_VECTOR_ELT(R_draw_list,5,GAMMA_draws);
  SET_VECTOR_ELT(R_draw_list,6,BETA_draws);
  
  SEXP R_drawnames;
  PROTECT(R_drawnames = allocVector(STRSXP, draw_list_length));
  nProtected++;
  SET_STRING_ELT(R_drawnames,0,mkChar("mu"));
  SET_STRING_ELT(R_drawnames,1,mkChar("Sigma"));
  SET_STRING_ELT(R_drawnames,2,mkChar("NNs"));
  SET_STRING_ELT(R_drawnames,3,mkChar("LAMBDA"));
  SET_STRING_ELT(R_drawnames,4,mkChar("TURNOUT"));
  SET_STRING_ELT(R_drawnames,5,mkChar("GAMMA"));
  SET_STRING_ELT(R_drawnames,6,mkChar("BETA"));
  setAttrib(R_draw_list, R_NamesSymbol, R_drawnames);

  if (dbg)
    Rprintf("done. Setting up restart object...\n");
  
  int restart_list_length = 18;
  SEXP R_restart_list = R_NilValue;
  SEXP ret_NNs = R_NilValue;
  SEXP ret_THETAS = R_NilValue;
  SEXP ret_OMEGAS = R_NilValue;
  SEXP ret_mu_vec_cu = R_NilValue;

  if (keep_restart_info){

  // Need to allocate an R object to hold the NNs (updated before returning):
  PROTECT(ret_NNs = allocMatrix(REALSXP,n,numcells_pt));
  PROTECT(ret_THETAS = allocMatrix(REALSXP,n,numcells_pt));
  PROTECT(ret_OMEGAS = allocMatrix(REALSXP,n,p));
  PROTECT(ret_mu_vec_cu = allocVector(REALSXP,numcols(mu_vec_cu)));
  nProtected+=4;

  PROTECT(R_restart_list=allocVector(VECSXP,restart_list_length));
  nProtected++;

  SET_VECTOR_ELT(R_restart_list,0, getListElement(args,"NNtots"));
  SET_VECTOR_ELT(R_restart_list,1, getListElement(args,"NNbounds"));
  SET_VECTOR_ELT(R_restart_list,2, getListElement(args,"numscans"));
  SET_VECTOR_ELT(R_restart_list,3, getListElement(args,"psi_0"));
  SET_VECTOR_ELT(R_restart_list,4, getListElement(args,"nu_0"));
  SET_VECTOR_ELT(R_restart_list,5, getListElement(args,"nolocalmode"));
  SET_VECTOR_ELT(R_restart_list,6, getListElement(args,"dof"));
  SET_VECTOR_ELT(R_restart_list,7, getListElement(args,"save_every"));
  SET_VECTOR_ELT(R_restart_list,8,getListElement(args,"keepNNinternals"));
  SET_VECTOR_ELT(R_restart_list,9,getListElement(args,"keepTHETAS"));
  SET_VECTOR_ELT(R_restart_list,10,getListElement(args,"mu_vec_0"));
  SET_VECTOR_ELT(R_restart_list,11,ret_mu_vec_cu);
  SET_VECTOR_ELT(R_restart_list,12,getListElement(args,"rho_vec"));
  SET_VECTOR_ELT(R_restart_list,13,getListElement(args,"use_Diri_every_vec"));
  SET_VECTOR_ELT(R_restart_list,14,ret_NNs);
  SET_VECTOR_ELT(R_restart_list,15,ret_THETAS);
  SET_VECTOR_ELT(R_restart_list,16,ret_OMEGAS);
  SET_VECTOR_ELT(R_restart_list,17,getListElement(args,"kappa_0"));
  
  SEXP R_restart_names;
  PROTECT(R_restart_names = allocVector(STRSXP, restart_list_length));
  nProtected++;
  SET_STRING_ELT(R_restart_names,0,mkChar("NNtots"));
  SET_STRING_ELT(R_restart_names,1,mkChar("NNbounds"));
  SET_STRING_ELT(R_restart_names,2,mkChar("numscans"));
  SET_STRING_ELT(R_restart_names,3,mkChar("psi.0"));
  SET_STRING_ELT(R_restart_names,4,mkChar("nu.0"));
  SET_STRING_ELT(R_restart_names,5,mkChar("nolocalmode"));
  SET_STRING_ELT(R_restart_names,6,mkChar("dof"));
  SET_STRING_ELT(R_restart_names,7,mkChar("save.every"));
  SET_STRING_ELT(R_restart_names,8,mkChar("keepNNinternals"));
  SET_STRING_ELT(R_restart_names,9,mkChar("keepTHETAS"));
  SET_STRING_ELT(R_restart_names,10,mkChar("mu.vec.0"));
  SET_STRING_ELT(R_restart_names,11,mkChar("mu.vec.cu"));
  SET_STRING_ELT(R_restart_names,12,mkChar("rho.vec"));
  SET_STRING_ELT(R_restart_names,13,mkChar("use.Diri.every.vec"));
  SET_STRING_ELT(R_restart_names,14,mkChar("NNs"));
  SET_STRING_ELT(R_restart_names,15,mkChar("THETAS"));
  SET_STRING_ELT(R_restart_names,16,mkChar("OMEGAS"));
  SET_STRING_ELT(R_restart_names,17,mkChar("kappa.0"));
  setAttrib(R_restart_list, R_NamesSymbol, R_restart_names);

  } // END if (keep_restart_info){...}
  
  if (dbg)
    Rprintf("done.\n");

  SEXP R_table_nrows;
  PROTECT(R_table_nrows = allocVector(INTSXP,1));
  nProtected++;
  INTEGER(R_table_nrows)[0] = INTEGER(getListElement(args,"numrows_pt"))[0];
  SEXP R_table_ncols;
  PROTECT(R_table_ncols = allocVector(INTSXP,1));
  nProtected++;
  INTEGER(R_table_ncols)[0] = INTEGER(getListElement(args,"numcols_pt"))[0];


  if (dbg){
    Rprintf("\nYou chose eschew = %d\n",eschew);
    Rprintf("\nYou passed NNtots, the first two rows of which are as follows.\n");
    matrix_print_subset(NNtots, 0, 1, 0, numcols(NNtots)-1);
    Rprintf("\n\nYou passed NNbounds, the first two rows of which are as follows.\n");
    matrix_print_subset(NNbounds, 0, 1, 0, numcols(NNbounds)-1);
    Rprintf("\nYou passed the following scalar parameters:\n");
    Rprintf("psi_0 = %g, nu_0 = %g, kappa_0 = %g.\n", psi_0, nu_0, kappa_0);
    Rprintf("numrows_pt = %d, numcols_pt = %d, num_iters = %d.\n",
  	    numrows_pt, numcols_pt, num_iters);
    Rprintf("nolocalmode = %f, dof = %f.\n", nolocalmode, dof);
    Rprintf("maxrange = %f, save_every = %f, num_scans = %d.\n", 
  	    maxrange, save_every, num_scans);
    Rprintf("how_many_NN_internals = %g, how_many_THETAS = %g\n", 
	    how_many_NN_internals, how_many_THETAS);
    Rprintf("\n\nCreated a starting NNs matrix, the first two rows of which are as follows.\n");
    matrix_print_subset(NNs, 0, 1, 0, numcols(NNs)-1);
    Rprintf("\n\nYou passed the following mu prior parameter vector.\n");
    matrix_print_all(mu_vec_0);
    Rprintf("\n\nYou passed the following mu starting vector.\n");
    matrix_print_all(mu_vec_cu);
    Rprintf("\n\nCreated THETAS, the first two rows of which are as follows.\n");
    matrix_print_subset(THETAS, 0, 1, 0, numcols(THETAS)-1);
    Rprintf("\n\nCreated an OMEGAS matrix, the first two rows of which are as follows.\n");
    matrix_print_subset(OMEGAS, 0, 1, 0, numcols(OMEGAS)-1);
    Rprintf("\n\nUnpacked rho_vec; its first two values are as follows.\n");
    Rprintf("%f\t%f\n", matrix_get_element(rho_vec, 0, 0),
  	    matrix_get_element(rho_vec, 0, 0));
    Rprintf("\n\nUnpacked use_Diri_every_vec; its first two values are as follows.\n");
    Rprintf("%f\t%f\n", matrix_get_element(use_Diri_every_vec, 0, 0),
  	    matrix_get_element(use_Diri_every_vec, 0, 0));
    Rprintf("Created the which_rc_int matrix.\n\n");
    Rprintf("The matrix of possible small Gibbs combinations is as follows.\n");
    matrix_print_all_int(which_rc_int);
    Rprintf("\n\nCreated the matrix PSI_0 as follows:\n");
    matrix_print_all(PSI_0);
    Rprintf("\n");
  }

  //  One final error check on the starting NNs
  check_ep_validity(NNs,MMs,KKs,numcells_pt,numrows_pt,numcols_pt);
  check_bounds_all(NNs,NNbounds,numcells_pt);

  //  If the user so desires, create matrices to save draws of NNs or THETAS
  Matrix *THETAS_save = NULL;
  Matrix *NNs_save = NULL;
  if (how_many_NN_internals){
    if (0.0 != fmod(num_iters, how_many_NN_internals))
      error("num.iters not divisible by how.many.NN.internals");
    
    NNs_save = rmatrix_new(n*numcells_pt, (index_t)how_many_NN_internals);
    Rprintf("\nCreated matrix to save draws of internal cell counts.\n");
    save_every_NN_internals = num_iters/how_many_NN_internals;
  }
  if (how_many_THETAS){
    if (0.0 != fmod(num_iters, how_many_THETAS))
      error("num.iters not divisible by how.many.THETAS");
    
    THETAS_save = rmatrix_new(n*numcells_pt, (index_t)how_many_THETAS);
    Rprintf("\nCreated matrix to save draws of internal cell probabilities.\n");
    save_every_THETAS = num_iters/how_many_THETAS;
  }

  Rprintf("Beginning MCMC routine...\n");

  //  RUN THE LOOP FOR ONE SET OF DRAWS
  for (ii=0; ii<num_iters; ii++){
 
        // Now do the same draw as before:
#ifdef _DBG_
        Rprintf("Drawing SIGMA...\n");
#endif
        draw_SIGMA(SIGMA_cu, nu_0, PSI_0, mu_vec_cu, OMEGAS, 
  	           tmat1, tmat2, tmat3, tmat4, AUG_in_matrix_inverse,dbg);

        //Rprintf("Cholesky and tidy-up...\n");
        matrix_cholesky(SIGMA_cu, SIGMA_chol_cu);
        m12_log_SIGMA_det_cu = -0.5*matrix_determinant(SIGMA_cu,tmat1,1);
        sample_equ_pr_wo_replace_int(zerotonumperms_int, ordervec_int);

#ifdef _DBG_
        Rprintf("Drawing MMs/NNs...\n");
#endif
        draw_NNs_MMs(NNs, MMs, KKs, NNprop_vec, NNbounds, ii, NNbounds_temp_vec, 
           NNtots, MMtots, KKtots, ff_vec,
	         NNtots_temp_vec, THETAS, ordervec_int, which_rc_int, nolocalmode, 
	         numrows_pt, numcols_pt, numcells_pt, num_scans,
					 NNs_prop, MMs_prop, NNs_curr, MMs_curr, tmp_KKs,
           multinomial_parameters, curr_row, prop_row, sr_probs, sr_reps,
					 vld_mv_p, acc_mv_p, NNs_count_use_multinom,
           lfactorial_vector);

        // Check the validity of the states:
        // May not want to check every time once code becomes stable:
#ifdef _DBG_
        Rprintf("Checking validity of NNs, MMs after draw_NNs_MMs()...\n");
        check_ep_validity(NNs,MMs,KKs,numcells_pt,numrows_pt,numcols_pt);
        Rprintf("passed.\n");
#endif

#ifdef _DBG_
        Rprintf("Drawing THETAS...\n");
#endif
        draw_THETAS_t_and_Dirichlet(THETAS, OMEGAS,	prop_THETA, prop_OMEGA,	
				    SIGMA_chol_cu, temp1_vec, temp2_vec, NNs,
				    mu_vec_cu, SIGMA_cu, AUG, acc_tv_p,
				    rho_vec, SIGMA_chol_cu_temp, acc_dv_p,
				    use_Diri_every_vec, THETAS_count_use_t, 
				    THETAS_count_use_Diri, dof, numrows_pt, 
				    numcols_pt, numcells_pt, ii,
						SIGMA_inverse_for_prop, tmpMean, tmpOut, 
				    tmpScalar, SIGMA_dims);

#ifdef _DBG_
        Rprintf("Drawing mu...\n");
#endif
        draw_mu(mu_vec_cu, kappa_0, mu_vec_0, OMEGAS, SIGMA_chol_cu, 
	        tvec1_dimOMs, tvec2_dimOMs, tmat1, tvec1_dimSI, tvec2_dimSI);

    // Should we print the iteration number to the user?
    if (dbg){
      if (((ii+1)%print_every)==0){
        Rprintf("Iter= %ld\t",ii+1);
        iter_count++;
        if(iter_count > 4){
           iter_count = 1;
           Rprintf("\n");
        }
      }
    }

    if (fmod(ii+1,save_every)==0.0){
      //  Error check the NNs:
      check_bounds_all(NNs, NNbounds, numcells_pt);
      check_ep_validity(NNs,MMs,KKs,numcells_pt,numrows_pt,numcols_pt);
      store_the_draws(mu_vec_cu,mu_draws,SIGMA_cu,SIGMA_draws,NNs,NNs_draws,
		      LAMBDA_draws,TURNOUT_draws,GAMMA_draws,BETA_draws,
		      numrows_pt,numcols_pt,&rownum_store,
          ncol_LAMBDA,ncol_BETA);
    }

    // Store the NNs state?
    if (how_many_NN_internals && (0.0 == fmod(ii+1, save_every_NN_internals))){
      check_bounds_all(NNs, NNbounds, numcells_pt);
      check_ep_validity(NNs,MMs,KKs,numcells_pt,numrows_pt,numcols_pt);
      store_internals(NNs, NNs_save, &colnum_store_NNs_internals);
    }

    // Store the THETA state?
    if (how_many_THETAS && (0.0 == fmod(ii+1, save_every_THETAS))){
      check_bounds_all(NNs, NNbounds, numcells_pt);
      check_ep_validity(NNs,MMs,KKs,numcells_pt,numrows_pt,numcols_pt);
      store_internals(THETAS, THETAS_save, &colnum_store_THETAS);
    }

  } // END LOOP FOR ONE FULL ITERATION

  Rprintf("\nFinished MCMC routine. Processing output...\n");

  //  Turn acceptance counts into fractions
  adjust_acc_vector(acc_THETAS_t_vec_ret,     THETAS_count_use_t);
  adjust_acc_vector(acc_THETAS_Diri_vec_ret,  THETAS_count_use_Diri);
  adjust_acc_matrix(vld_NNs_multinom_mat_ret, NNs_count_use_multinom);
  adjust_acc_matrix(acc_NNs_multinom_mat_ret, NNs_count_use_multinom);

  //  If desired, package draws of the internal cell probabilities
  //        for passing back to R
  if (how_many_THETAS){
    THETAS_save_ret = matrix_repack_new(THETAS_save);
    nProtected++;
  }

  //  If desired, package draws of the internal cell counts
  //        for passing back to R
  if (how_many_NN_internals){
    NNs_save_ret = matrix_repack_new(NNs_save);
    nProtected++;
  } 


  if (keep_restart_info){

  // Store the final state of the mu_vec_cu, NNs, THETAS and OMEGAS:
  const index_t nr_mu_vec_cu = numrows(mu_vec_cu);
  const index_t nc_mu_vec_cu = numcols(mu_vec_cu);
  double *ret_mu_vec_cu_p = REAL(ret_mu_vec_cu);
  double *ret_NNs_p = REAL(ret_NNs);
  double *ret_THETAS_p = REAL(ret_THETAS);
  double *ret_OMEGAS_p = REAL(ret_OMEGAS);
  for (ii=0; ii<nc_mu_vec_cu; ii++)
    ret_mu_vec_cu_p[ii] = matrix_fast_get_element(mu_vec_cu,0,ii,nr_mu_vec_cu);
  for (ii=0; ii<n; ii++){
    for (jj=0; jj<numcells_pt; jj++){
      ret_NNs_p[ii+n*jj] = matrix_fast_get_element(NNs,ii,jj,n);
      ret_THETAS_p[ii+n*jj] = matrix_fast_get_element(THETAS,ii,jj,n);
    }
  }
  for (ii=0; ii<n; ii++)
    for (jj=0; jj<p; jj++)
      ret_OMEGAS_p[ii+n*jj] = matrix_fast_get_element(OMEGAS,ii,jj,n);

  } // END if (keep_restart_info){...}

  //  Create the final return list
  listlength = 8;
  if (how_many_THETAS) 
    listlength++;
  if (how_many_NN_internals) 
    listlength++;

  SEXP list_ret;
  PROTECT(list_ret = allocVector(VECSXP, listlength));
  nProtected++;
  SET_VECTOR_ELT(list_ret, 0, R_draw_list);
  SET_VECTOR_ELT(list_ret, 1, acc_THETAS_t_vec_ret);
  SET_VECTOR_ELT(list_ret, 2, acc_THETAS_Diri_vec_ret);
  SET_VECTOR_ELT(list_ret, 3, vld_NNs_multinom_mat_ret);
  SET_VECTOR_ELT(list_ret, 4, acc_NNs_multinom_mat_ret);
  SET_VECTOR_ELT(list_ret, 5, R_table_nrows);
  SET_VECTOR_ELT(list_ret, 6, R_table_ncols);
  nextplace = 7;
  if (how_many_THETAS){
    SET_VECTOR_ELT(list_ret, nextplace, THETAS_save_ret);
    nextplace++;
  }
  if (how_many_NN_internals){
    SET_VECTOR_ELT(list_ret, nextplace, NNs_save_ret);
    nextplace++;
  }
  SET_VECTOR_ELT(list_ret, nextplace, R_restart_list);

  //  Assign names for the return list
  SEXP names_ret; 
  PROTECT(names_ret = allocVector(STRSXP, listlength));
  nProtected++;
  SET_STRING_ELT(names_ret, 0, mkChar("draws"));
  SET_STRING_ELT(names_ret, 1, mkChar("acc.t"));
  SET_STRING_ELT(names_ret, 2, mkChar("acc.Diri"));
  SET_STRING_ELT(names_ret, 3, mkChar("vld.multinom"));
  SET_STRING_ELT(names_ret, 4, mkChar("acc.multinom"));
  SET_STRING_ELT(names_ret, 5, mkChar("numrows.pt"));
  SET_STRING_ELT(names_ret, 6, mkChar("numcols.pt"));
  nextplace = 7;
  if(how_many_THETAS){
    SET_STRING_ELT(names_ret, nextplace, mkChar("THETAS.saved"));
    nextplace++;
  }
  if (how_many_NN_internals){
    SET_STRING_ELT(names_ret, nextplace, mkChar("NN.internals.saved"));
    nextplace++;
  }
  SET_STRING_ELT(names_ret, nextplace, mkChar("restart.info"));
  //  Marry the names to the return list
  setAttrib(list_ret, R_NamesSymbol, names_ret);

  PutRNGstate();

  if (dbg)
    Rprintf("\n");
  
  UNPROTECT(nProtected);
  return list_ret;
}

