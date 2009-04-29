
#ifndef _JIMSMH_H
#define _JIMSMH_H

#include "jimsmatrix.h"

/****************************************************************************/
/*****************  M-H CALCULATION FUNCTION DECLARATIONS  ******************/
/****************************************************************************/


double
log_THETAS_proposal_product_Dirichlet(Matrix *theta, long row,
				      Matrix *NNs,
				      long num_p);


double
log_p_target_theta_short(Matrix *in_vec_THETA, const index_t THETA_row,
			 Matrix *in_OMEGA, const index_t OMEGA_row,
                         const index_t num_p, Matrix *NNs, Matrix *mu_vec_cu,			     
			 Matrix *SIGMA_cu, Matrix *AUG, const index_t numrows_pt, 
			 const index_t numcols_pt,
                         Matrix *SIGMA_inverse_for_prop,
                         Matrix *tmpMean, Matrix *tmpOut, Matrix *tmpScalar);

double
log_THETAS_proposal_t_jacobian(Matrix *prop_OMEGA, Matrix *prop_THETA,
			       Matrix *THETAS, const index_t num_p, const index_t numrows_pt,
			       const index_t numcols_pt_m1, int is_prop);

// I don't think this function ever gets used?
double
log_p_dirichlet(double *in_vec, const index_t num_p, Matrix *NNs, 
                const index_t numrows_pt, const index_t numcols_pt);


double
log_p_target_NNs(Matrix *in, const index_t row, const index_t num_p, Matrix *THETAS, const index_t numcells_pt);


double
log_p_NNs_prop_anywhere(Matrix *NNs, Matrix *NNbounds, Matrix *NNbounds_temp_vec, 
			Matrix *NNtots, Matrix *NNtots_temp_vec, const index_t num_p, 
			const index_t numrows_pt, const index_t numcols_pt, const index_t numcells_pt);

double
log_NNs_multinomial_mh_ratio(Matrix * const curr_row,
														 Matrix * const prop_row,
														 Matrix * const multinomial_parameters,
                             Matrix const * const lfactorial_vector);

// Exit poll version is considerably more complicated since it needs all
// of the matrices. The regular one depends on the special row.

double
log_MMs_multinomial_mh_ratio(Matrix * const NNs_prop,
														 Matrix * const MMs_prop,
														 Matrix * const NNs_curr,
														 Matrix * const MMs_curr,
														 Matrix * const tmp_KKs,
														 Matrix * const THETAS,
                             const index_t current_precinct,
                             const index_t special_row,
                             const index_t numrows_pt,
                             const index_t numcols_pt,
														 Matrix * const multinomial_parameters,
                             Matrix const * const lfactorial_vector);

#endif /* jimsMH.h */
