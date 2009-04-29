#ifndef _JIMSRANDOM_H
#define _JIMSRANDOM_H

#include "jimsmatrix.h"
#include "jimsimatrix.h"
#include "jimsmisc.h"
#include "jimsMH.h"
#include "nchg.h"

void
draw_NNs_indep_start(Matrix * const NNprop_vec, 
	   	   Matrix * const NNbounds, 
		     Matrix * const NNbounds_temp_vec, 
		     Matrix * const NNtots, 
		     Matrix * const NNtots_temp_vec, 
		     const index_t num_p, 
		     const index_t numrows_pt, 
		     const index_t numcols_pt, 
		     const index_t numcells_pt);

void draw_NNs_MMs_indep_start(
         Matrix * const NNs,
         Matrix * const MMs,
         Matrix * const KKs,
         Matrix * const NNbounds,
         Matrix * const MMbounds,
         Matrix * const NNtots,
         Matrix * const MMtots,
         Matrix * const NNbounds_temp_vec,
         Matrix * const NNtots_temp_vec,
         const index_t numrows_pt,
         const index_t numcols_pt,
         const index_t numcells_pt);

void
draw_THETAS_from_NNs_start(Matrix * const THETAS, 
			   Matrix * const NNs, 
			   Matrix * const NNtots,
			   const index_t numrows_pt,
			   const index_t numcols_pt);

void
draw_THETAS_t_and_Dirichlet(Matrix * const THETAS, 
			    Matrix * const OMEGAS,
			    Matrix * const prop_THETA, 
			    Matrix * const prop_OMEGA,
			    Matrix * const SIGMA_chol_cu, 
			    Matrix * const temp1_vec,
			    Matrix * const temp2_vec, 
			    Matrix * const NNs,
			    Matrix * const mu_vec_cu, 
			    Matrix * const SIGMA_cu,
			    Matrix * const AUG, 
			    double * const acc_THETAS_t_vec,
			    Matrix * const rho_vec, 
			    Matrix * const SIGMA_chol_cu_temp,
			    double * const acc_THETAS_Diri_vec,
			    Matrix * const use_Diri_every_vec,
			    Matrix * const THETAS_count_use_t, 
			    Matrix * const THETAS_count_use_Diri,
			    const double dof, 
			    const index_t numrows_pt, 
			    const index_t numcols_pt,	
			    const index_t numcells_pt, 
			    const long iternum,
          Matrix * const SIGMA_inverse_for_prop,
          Matrix * const tmpMean, 
			    Matrix * const tmpOut, 
			    Matrix * const tmpScalar,
			    Matrix * const SIGMA_dims);

void
draw_THETAS_t_and_Dirichlet_with_covariates(Matrix * const THETAS, 
			    Matrix * const OMEGAS,
			    Matrix * const prop_THETA, 
			    Matrix * const prop_OMEGA,
			    Matrix * const SIGMA_chol_cu, 
			    Matrix * const temp1_vec,
			    Matrix * const temp2_vec, 
			    Matrix * const NNs,
			    Matrix * const mu_mat_cu, 
			    Matrix * const SIGMA_cu,
			    Matrix * const AUG, 
			    double * const acc_THETAS_t_vec,
			    Matrix * const rho_vec, 
			    Matrix * const SIGMA_chol_cu_temp,
			    double * const acc_THETAS_Diri_vec,
			    Matrix * const use_Diri_every_vec,
			    Matrix * const THETAS_count_use_t, 
			    Matrix * const THETAS_count_use_Diri,
			    const double dof, 
			    const index_t numrows_pt, 
			    const index_t numcols_pt,	
			    const index_t numcells_pt, 
			    const long iternum,
          Matrix * const SIGMA_inverse_for_prop,
          Matrix * const tmpMean, 
			    Matrix * const tmpOut, 
			    Matrix * const tmpScalar,
			    Matrix * const SIGMA_dims,
          Matrix * const tmp_mu_vec_cu);

void
draw_THETAS_Dirichlet_independent_one_row(Matrix * const THETAS, 
					  Matrix * const OMEGAS,
					  Matrix * const prop_THETA, 
					  Matrix * const prop_OMEGA, 
					  Matrix * const SIGMA_chol_cu, 
					  Matrix * const temp1_vec, 
					  Matrix * const temp2_vec, 
					  Matrix * const NNs, 
					  Matrix * const mu_vec_cu, 
					  Matrix * const SIGMA_cu, 
					  Matrix * const AUG, 
					  double * const acc_THETAS_Diri_vec, 
					  const index_t numrows_pt, 
					  const index_t numcols_pt, 
					  const index_t numcells_pt, 
					  const index_t which_row,
                                          Matrix * const SIGMA_inverse_for_prop,
                                          Matrix * const tmpMean, 
					  Matrix * const tmpOut, 
					  Matrix * const tmpScalar);

void
draw_THETAS_Dirichlet_independent(Matrix * const THETAS, 
				  Matrix * const OMEGAS,
				  Matrix * const prop_THETA, 
				  Matrix * const prop_OMEGA, 
				  Matrix * const SIGMA_chol_cu, 
				  Matrix * const temp1_vec,
				  Matrix * const temp2_vec, 
				  Matrix * const NNs,
				  Matrix * const mu_vec_cu, 
				  Matrix * const SIGMA_cu,
				  Matrix * const AUG, 
				  double * const acc_THETAS_Diri_vec,
				  const index_t numrows_pt, 
				  const index_t numcols_pt,
				  const index_t numcells_pt,
                                  Matrix * const SIGMA_inverse_for_prop,
                                  Matrix * const tmpMean, 
				  Matrix * const tmpOut, 
				  Matrix * const tmpScalar);

void
r_product_Dirichlet(Matrix * const prop_THETA, 
		    Matrix * const NNs, 
		    const index_t num_p,
		    const index_t numrows_pt, 
		    const index_t numcols_pt);

void
draw_THETAS_t_dependent_one_row(Matrix * const THETAS, 
				Matrix * const OMEGAS,
				Matrix * const prop_THETA, 
				Matrix * const prop_OMEGA,	
				Matrix * const SIGMA_chol_cu, 
				Matrix * const temp1_vec, 
				Matrix * const temp2_vec, 
				Matrix * const NNs,		
				Matrix * const mu_vec_cu, 
				Matrix * const SIGMA_cu,	
				Matrix * const AUG, 
				double * const acc_THETAS_t_vec,	
				Matrix * const rho_vec, 
				Matrix * const SIGMA_chol_cu_temp, 
				const double dof, 
				const index_t numrows_pt, 
				const index_t numcols_pt, 
				const index_t numcells_pt, 
				const index_t which_row,
                                Matrix * const SIGMA_inverse_for_prop,
                                Matrix * const tmpMean, 
				Matrix * const tmpOut, 
				Matrix * const tmpScalar);

void
draw_THETAS_t_dependent(Matrix * const THETAS, 
			Matrix * const OMEGAS,
			Matrix * const prop_THETA, 
			Matrix * const prop_OMEGA,
			Matrix * const SIGMA_chol_cu, 
			Matrix * const temp1_vec,
			Matrix * const temp2_vec, 
			Matrix * const NNs,
			Matrix * const mu_vec_cu, 
			Matrix * const SIGMA_cu,
			Matrix * const AUG, 
			double * const acc_THETAS_t_vec, 
			Matrix * const rho_vec, 
			Matrix * const SIGMA_chol_cu_temp,
			const double dof, 
			const index_t numrows_pt, 
			const index_t numcols_pt,
			const index_t numcells_pt,
                        Matrix * const SIGMA_inverse_for_prop,
                        Matrix * const tmpMean, 
			Matrix * const tmpOut, 
			Matrix * const tmpScalar);

void
mvrnorm_c_chol(Matrix * const xx, 
	       Matrix * const mu_vec, 
	       Matrix * const SIGMA_chol,
	       Matrix * const temp1_vec, 
	       Matrix * const temp2_vec);

void
rinvWis_c(Matrix * const xx, 
	  const double dof, 
	  Matrix * const SS,	
	  Matrix * const BB, 
	  Matrix * const CC, 
	  Matrix * const MM,
	  Matrix * const AUG_in_matrix_inverse);

void
draw_mu(Matrix * const mu_vec_cu, 
	const double kappa_0, 
	Matrix * const mu_vec_0,
	Matrix * const OMEGAS, 
	Matrix * const SIGMA_chol_cu,
	Matrix * const omegas_means_vec, 
	Matrix * const mean_param_vec,
	Matrix * const SIGMA_chol_cu_temp, 
	Matrix * const tvec1_dimSI,
	Matrix * const tvec2_dimSI);

void
draw_eta_with_covariates(Matrix * const mu_mat_cu, 
  Matrix_p * const Xmatrix_list,
  Matrix * const eta_vec_cu,
	Matrix * const eta_vec_0, // What is this?
	Matrix * const OMEGAS, 
	Matrix * const SIGMA_cu, 
	Matrix * const SIGMA_inverse,
  Matrix * const mean_vector_1_x_d, // 1 x d vector
	Matrix * const tvec1_1_x_d, // 1 x d vector
	Matrix * const tvec2_1_x_d, // 1 x d vector
  Matrix * const tmat1_d_x_d, // d x d matrix
  Matrix * const tmat2_p_x_p, // p x p matrix
  Matrix * const tmat3_p_x_d, // p x d matrix
  Matrix * const tmat4_d_x_d, // d x d matrix
  Matrix * const tmat5_d_x_p, // d x p matrix
  Matrix * const tmat6_d_x_d, // d x d matrix
  Matrix * const V_0_inverse);

void
draw_SIGMA(Matrix * const SIGMA_cu, 
	   const double nu_0, 
	   Matrix * const PSI_0, 
	   Matrix * const mu_vec_cu, 
	   Matrix * const OMEGAS, 
	   Matrix * const SS,
	   Matrix * const tmat2, 
	   Matrix * const tmat3, 
	   Matrix * const tmat4,
	   Matrix * const AUG_in_matrix_inverse,
	   const index_t dbg);

void
draw_SIGMA_with_covariates(Matrix * const SIGMA_cu, 
	   const double nu_0, 
	   Matrix * const PSI_0, 
	   Matrix * const mu_mat_cu, 
	   Matrix * const OMEGAS, 
	   Matrix * const SS,
	   Matrix * const tmat2, 
	   Matrix * const tmat3, 
	   Matrix * const tmat4,
	   Matrix * const AUG_in_matrix_inverse,
	   const index_t dbg);

double
draw_NNs_prop_anywhere(Matrix * const NNprop_vec, 
		       Matrix * const NNbounds, 
		       Matrix * const NNbounds_temp_vec,
		       Matrix * const NNtots, 
		       Matrix * const NNtots_temp_vec, 
		       const index_t num_p, 	
		       const index_t numrows_pt, 
		       const index_t numcols_pt, 
		       const index_t numcells_pt);

void
draw_NNs_anywhere(Matrix * const NNs, 
		  Matrix * const NNprop_vec, 
		  Matrix * const NNbounds, 
		  Matrix * const NNbounds_temp_vec, 
		  Matrix * const NNtots, 
		  Matrix * const NNtots_temp_vec, 
		  Matrix * const THETAS, 
		  const index_t numrows_pt, 
		  const index_t numcols_pt, 
		  const index_t numcells_pt);

void
rGibbsNNs(Matrix * const NNs, 
	  const index_t num_p, 
	  Matrix * const THETAS, 
	  Matrix_int * const which_rc_int, 
	  double * const ff_vec, 
	  const index_t whichperm, 
	  const index_t numrows_pt, 
  	  const index_t numcols_pt);

void
draw_NNs(Matrix * const NNs, 
	 Matrix * const NNprop_vec, 
	 Matrix * const NNbounds, 
	 const long iternum, 
	 Matrix * const NNbounds_temp_vec, 
	 Matrix * const NNtots, 
	 double * const ff_vec, 
	 Matrix * const NNtots_temp_vec, 
	 Matrix * const THETAS, 
	 Matrix_int * const ordervec_int,
	 Matrix_int * const which_rc_int, 
	 const double nolocalmode, 
	 const index_t numrows_pt, 
	 const index_t numcols_pt, 
	 const index_t numcells_pt, 
	 const index_t num_scans,
	 Matrix * const NNs_prop,
	 Matrix * const multinomial_parameters,
	 Matrix * const curr_row,
	 Matrix * const prop_row,
	 Matrix * const sr_probs,
	 Matrix_int * const sr_reps,
 	 double *vld_mv_p, // Counts the number of valid MH proposals for each precinct
 	 double *acc_mv_p, // Counts the number of accepted MH proposals for each precinct
	 Matrix * const NNs_count_use_multinom, // Counts the number of MH proposals for each precinct
   Matrix const * const lfactorial_vector); // Log-factorial lookup table

void // Exit poll version:
draw_NNs_MMs(Matrix * const NNs, 
   Matrix * const MMs,
   Matrix * const KKs,
	 Matrix * const NNprop_vec, 
	 Matrix * const NNbounds, 
	 const long iternum, 
	 Matrix * const NNbounds_temp_vec, 
	 Matrix * const NNtots, 
	 Matrix * const MMtots, 
	 Matrix * const KKtots, 
	 double * const ff_vec, 
	 Matrix * const NNtots_temp_vec, 
	 Matrix * const THETAS, 
	 Matrix_int * const ordervec_int,
	 Matrix_int * const which_rc_int, 
	 const double nolocalmode, 
	 const index_t numrows_pt, 
	 const index_t numcols_pt, 
	 const index_t numcells_pt, 
	 const index_t num_scans,
	 Matrix * const NNs_prop,
	 Matrix * const MMs_prop,
	 Matrix * const NNs_curr,				// Place-holder for the current state
	 Matrix * const MMs_curr,				// Place-holder for the current state
	 Matrix * const tmp_KKs,				// Place-holder for the KK state in precinct ii
	 Matrix * const multinomial_parameters,
	 Matrix * const curr_row,
	 Matrix * const prop_row,
	 Matrix * const sr_probs,
	 Matrix_int * const sr_reps,
 	 double *vld_mv_p, // Counts the number of valid MH proposals for each precinct
 	 double *acc_mv_p, // Counts the number of accepted MH proposals for each precinct
	 Matrix * const NNs_count_use_multinom, // Counts the number of MH proposals for each precinct
   Matrix const * const lfactorial_vector); // Log-factorial lookup table

void
mvrt_c_chol(Matrix * const xx, 
	    Matrix * const mu_mat, 
	    const index_t mu_row,
	    Matrix * const SIGMA_chol_cu, 
	    const double dof,
	    Matrix * const temp1_vec, 
	    Matrix * const temp2_vec);

void
mvrt_c_chol(Matrix * const xx, 
	    Matrix * const mu_mat, 
	    const index_t mu_row,
	    Matrix * const SIGMA_chol_cu, 
	    const double dof,
	    Matrix * const temp1_vec, 
	    Matrix * const temp2_vec);

void
draw_NNs_multinomial_MH(Matrix * const NNs,							// Current cell counts
												Matrix * const NNtots,					// Row and col sums
												const index_t current_precinct,	// The precinct we are updating
												const index_t special_row,			// The row to leave alone
												Matrix * const THETAS,					// The current state of THETA
												const index_t numrows_pt,				// Number of rows in each table
												const index_t numcols_pt,				// Number of cols in each table
												Matrix * const NNs_prop,				// Place-holder for the proposed state
												Matrix * const multinomial_parameters,	// Place-holder for the THETA's in each row
												Matrix * const curr_row,				// Place-holder the current cell counts in a row
												Matrix * const prop_row,				// Place-holder the new cell counts in a row
												double *vld_mv_p, // Counts the number of valid MH proposals for each precinct
												double *acc_mv_p, // Counts the number of accepted MH proposals for each precinct
												Matrix * const NNs_count_use_multinom, // Counts the MH proposals for each precinct
                        Matrix const * const lfactorial_vector); // Log-factorial lookup table

void
draw_MMs_multinomial_MH(Matrix * const NNs,							// Current cell counts
                        Matrix * const MMs,             // Non-EP-ed 'voters'
                        Matrix * const KKs,             // Exit poll numbers
												Matrix * const NNtots,					// Row and col sums
												Matrix * const MMtots,					// Row and col 'non-voter' sums
												Matrix * const KKtots,					// Row and col exit poll sums
												const index_t current_precinct,	// The precinct we are updating
												const index_t special_row,			// The row to leave alone
												Matrix * const THETAS,					// The current state of THETA
												const index_t numrows_pt,				// Number of rows in each table
												const index_t numcols_pt,				// Number of cols in each table
												Matrix * const NNs_prop,				// Place-holder for the proposed state
												Matrix * const MMs_prop,				// Place-holder for the proposed state
												Matrix * const NNs_curr,				// Place-holder for the current state
												Matrix * const MMs_curr,				// Place-holder for the current state
												Matrix * const tmp_KKs,				  // Place-holder for the KK state in precinct ii
												Matrix * const multinomial_parameters,	// Place-holder for the THETA's in each row
												Matrix * const curr_row,				// Place-holder the current cell counts in a row
												Matrix * const prop_row,				// Place-holder the new cell counts in a row
												double *vld_mv_p, // Counts the number of valid MH proposals for each precinct
												double *acc_mv_p, // Counts the number of accepted MH proposals for each precinct
												Matrix * const NNs_count_use_multinom, // Counts the MH proposals for each precinct
                        Matrix const * const lfactorial_vector); // Log-factorial lookup table

void rmultinomial(Matrix * const draw,
									Matrix * const p_vector,
									const double total_count_d);

#endif /* jimsrandom.h */
