
#include "jimsMH.h"

/****************************************************************************/
/*****************  M-H CALCULATION FUNCTION DEFINITIONS   ******************/
/****************************************************************************/


double
log_THETAS_proposal_product_Dirichlet(Matrix *theta, long row,
				      Matrix *NNs,
				      long num_p)
{
  /*  Returns the log of the kernel of the product Dirichlet proposal
      density for the THETAS.  Uses a pointer to avoid unnecessary
      copying.  Error checking done several function layers above.
      Recall that the proposal uses the value of the relevant scalar
      from the NNs matrix + .1, to assure that the density is proper.
      
      Modified: 02/22/08 -- Now takes matrix and row argument to 
      remove dependence on matrix implementation. */

  double ret_val = 0.0;
  long ii;

  for (ii=0; ii<numcols(NNs); ii++)
    ret_val += (matrix_get_element(NNs, num_p, ii)-.9)*
      log(matrix_get_element(theta,row,ii));

  return ret_val;
}

double
log_p_target_theta_short(Matrix *in_THETA, const index_t THETA_row,
			 Matrix *in_OMEGA, const index_t OMEGA_row, 
			 const index_t num_p, Matrix *NNs, Matrix *mu_vec_cu,  
			 Matrix *SIGMA_cu, Matrix *AUG, const index_t numrows_pt,
			 const index_t numcols_pt,
                         Matrix *SIGMA_inverse_for_prop,
                         Matrix *tmpMean, Matrix *tmpOut, Matrix *tmpScalar)
{
  //  Returns the log of the target density for the theta vector, conditional
  //  on all other parameters and the internal counts, under the GQ Model.
  //  The pointers *in_vec_THETA and *in_vec_OMEGA take advantage
  //  of the Greiner matrix structure, in which matrices are stored by row.
  //  Thus, one uses this function for a proposal by setting in_vec_THETA
  //  and in_vec_OMEGA to pointers to the corresponding proposal vectors.
  //  One uses this function for a current value by setting the two pointers
  //  to the location of the relevant rows in the THETAS and OMEGAS
  //  matrix.  NNs must be organized so that
  //  each row represents a precinct.

  // TODO: Modify to include covariates

  double ret_val;
  index_t ii;
  const index_t nrAUG=numrows(AUG);
  const index_t ncNNs=numcols(NNs);
#ifdef _USE_ADJUST_
  double aa_temp=0.0;
  //  First the exponentiated portion of the normal term.
  //     Fill the augmented matrix.
  matrix_get_set_block(AUG, 0,numrows(AUG)-1,0,numcols(AUG)-2, 
		       SIGMA_cu, 0,numrows(SIGMA_cu)-1,0,numcols(SIGMA_cu)-1);

  for (ii=0; ii<nrAUG; ii++){
    matrix_set_element(AUG, ii, numcols(AUG)-1, 
		       matrix_get_element(in_OMEGA,OMEGA_row,ii) - matrix_get_element(mu_vec_cu, 0, ii));
  }
  //     ADJUST the augmented matrix to find SIGMA^(-1) %*% (omega_i - mu)
  for (ii=0; ii<nrAUG; ii++)
    matrix_ADJUST(AUG, ii);
  
  //     Calculated, element by element, t(omega_i - mu) %*% vector from just above.
  for (ii=0; ii<nrAUG; ii++){
    aa_temp += (matrix_get_element(in_OMEGA,OMEGA_row,ii) - matrix_get_element(mu_vec_cu, 0, ii)) * 
		matrix_get_element(AUG, ii, numcols(AUG)-1);
  }
  ret_val = -.5 * aa_temp;
#endif
#ifndef _USE_ADJUST_
#ifdef _DBG_
  // tmpMean : nrAug-by-1 matrix holding (omega_i - mu)
  // tmpOut  : nrAug-by-1 matrix holding SIGMA^{-1} %*% tmpMean
  // tmpScalar : 1-by-1 matrix holding tmpMean %*% SIGMA^{-1} %*% tmpMean

  int dbg=0;
  if (dbg)
    Rprintf("Inside log_p_target_theta_short...\n");

#endif
  matrix_get_row(in_OMEGA,OMEGA_row,tmpMean);

#ifdef _DBG_
  if (dbg){
    Rprintf("Successfully called matrix_get_row.\n");
    Rprintf("omega_i:\n");
    matrix_print_all(tmpMean);
  }
#endif

  const index_t nr_tmpMean = numrows(tmpMean);
  const index_t nr_mvc = numrows(mu_vec_cu);
  for (ii=0; ii<nrAUG; ii++)
    matrix_fast_set_element(tmpMean,ii,0,nr_tmpMean,
      matrix_fast_get_element(tmpMean,ii,0,nr_tmpMean)-
      matrix_fast_get_element(mu_vec_cu,0,ii,nr_mvc));

#ifdef _DBG_
  if (dbg){
    Rprintf("Successfully called matrix_get_row.\n");
    Rprintf("omega_i-mu:\n");
    matrix_print_all(tmpMean);
    Rprintf("Sigma^{-1}:\n");
    matrix_print_all(SIGMA_inverse_for_prop);
  }
#endif

  matrix_multiply(SIGMA_inverse_for_prop,'N',tmpMean,'N',tmpOut);

#ifdef _DBG_
  if (dbg){
    Rprintf("Successfully called matrix_get_row.\n");
    printf("Sigma^{-1}(omega_i-mu):\n");
    matrix_print_all(tmpOut);
  }
#endif

  matrix_multiply(tmpMean,'T',tmpOut,'N',tmpScalar);

#ifdef _DBG_
  if (dbg){
    Rprintf("Successfully called matrix_get_row.\n");
    Rprintf("(omega_i-mu)^{T}Sigma^{-1}(omega_i-mu):\n");
    matrix_print_all(tmpScalar);
  }
#endif

  ret_val = -0.5*matrix_get_element(tmpScalar,0,0);

#ifdef _DBG_
  if (dbg)
    Rprintf("Finished log_p_target_theta_short.\n");
#endif
#endif // End _USE_ADJUST_
  //  Take care of the Diri-looking pieces.
  const index_t nr_NNs = numrows(NNs);
  const index_t nr_in_THETA = numrows(in_THETA);
  for (ii=0; ii<ncNNs; ii++)
    ret_val += (matrix_fast_get_element(NNs, num_p,ii,nr_NNs) - 1) * 
		log(matrix_fast_get_element(in_THETA,THETA_row,ii,nr_in_THETA));

  return ret_val;
}



double
log_THETAS_proposal_t_jacobian(Matrix *prop_OMEGA, Matrix *prop_THETA,
			       Matrix *THETAS, const index_t num_p, const index_t numrows_pt,
			       const index_t numcols_pt_m1, int is_prop)
{
  //  This function returns the sum of the negative logs of the
  //  appropriate vector of THETAS.  If is_prop == 1, then the
  //  "appropriate vector" is the M-H proposal; because such a
  //  proposal will in OMEGA form, this function first transforms
  //  the OMEGA proposal vector into a THETA proposal vector and
  //  puts it in prop_THETA.  If is_prop == 0, then the "appropriate
  //  vector" is the current THETA vector for that precinct, and the
  //  return value is taken from the THETAS matrix.
  //  Error checking done several function layers up.

  double ret_val=0.0;

  if (is_prop){  
    // Proposal calculation:
    const index_t nr_PO = numrows(prop_OMEGA);
    const index_t nr_PT = numrows(prop_THETA);
    index_t ii, jj, rmult, position = 0;
    double denom, log_denom, temp_val;
    for (ii=0; ii<numrows_pt; ii++){
      denom = 1.0;
      rmult = ii*numcols_pt_m1;
      for (jj=0; jj<numcols_pt_m1; jj++)
	denom += exp(matrix_fast_get_element(prop_OMEGA, 0,rmult+jj,nr_PO));
      log_denom = log(denom);
      for (jj=0; jj<numcols_pt_m1; jj++){
	temp_val = matrix_fast_get_element(prop_OMEGA, 0, rmult+jj,nr_PO);
	ret_val += log_denom - temp_val;// log of the relevant omega
	matrix_fast_set_element(prop_THETA, 0,position,nr_PT, exp(temp_val)/denom);
	position++;
      }
      ret_val += log_denom;
      matrix_fast_set_element(prop_THETA, 0,position,nr_PT, 1/denom);
      position++;
    }
    return ret_val;

  } else {  //current value calculation
    const index_t nc_theta = numcols(THETAS);
    const index_t nr_TH = numrows(THETAS);
    index_t ii;
    for (ii=0; ii<nc_theta; ii++)
      ret_val -= log(matrix_fast_get_element(THETAS, num_p,ii,nr_TH));
    return ret_val;
  }
}

// I don't think this function ever gets used?
double
log_p_dirichlet(double *in_vec, const index_t num_p, Matrix *NNs, 
                const index_t numrows_pt, const index_t numcols_pt)
{
  /*  Returns the log of the product of the dirichlet densities at the value
      in_vec (which has numrows_pt such densities), using the num_pth row of
      NNs as the nunrows_pt parameter vectors.  Requires the R lgamma function.  */

  index_t ii, jj;
  double ret_val = 0.0, aa_NN, aa_NN_tot = 0.0;

  for (ii=0; ii<numrows_pt; ii++){
    aa_NN_tot = 0.0;
    for (jj=0; jj<numcols_pt; jj++){
      aa_NN = matrix_get_element(NNs, num_p, (ii*numcols_pt)+jj);
      aa_NN_tot += aa_NN;
      ret_val -= lgamma(aa_NN);
      ret_val += (aa_NN-1) * log(in_vec[(ii*numcols_pt)+jj]);
    }
    ret_val += lgamma(aa_NN_tot);
  }

  return ret_val;
}


double
log_p_target_NNs(Matrix *in, const index_t row, const index_t num_p,
                 Matrix *THETAS, const index_t numcells_pt)
{
  /*  Returns the log of the value of the conditional distribution
      of the a vector of precinct NNs called in_vec (arranged as follows:
      N_11, N_12, . . . N_1C, N_21, N_22, . . . N_2C, . . ., . . ., 
      N_R1, N_R2, . . . N_RC.  Note that *in_vec is a pointer, not a
      Greiner matrix, to speed up this workhorse function by (i) avoiding
      the necessity of making a copy of the input vector, and (ii) using
      the same function for proposals and current NNs values. */

  index_t ii;
  double ret_val = 0.0;
  const index_t nr_in = numrows(in);
  const index_t nr_TH = numrows(THETAS);

  for (ii=0; ii<numcells_pt; ii++){
    ret_val += matrix_fast_get_element(in,row,ii,nr_in)*log(matrix_fast_get_element(THETAS, num_p,ii,nr_TH));
    ret_val -= lgamma(matrix_fast_get_element(in,row,ii,nr_in)+1.0);
  }
  return ret_val;
}

double
log_p_NNs_prop_anywhere(Matrix *NNs, Matrix *NNbounds, Matrix *NNbounds_temp_vec, 
			Matrix *NNtots, Matrix *NNtots_temp_vec, const index_t num_p,
			const index_t numrows_pt, const index_t numcols_pt, 
                        const index_t numcells_pt)
{
  /*  Returns the log of the probibility of drawing the CURRENT value of the
      num_pth precinct's table counts from the anywhere proposal
      distribution.  Anything ending in _vec must be a row vector.
      User allocates all memory.  Error checking to be done elsewhere.  */

  index_t rr, cc, ii, position;
  const index_t numrows_pt_m1 = numrows_pt-1;
  const index_t numcols_pt_m1 = numcols_pt-1;
  double range_NN, ret_val= 0.0, new_bound, new_tot;
  const index_t nr_NNBTV = numrows(NNbounds_temp_vec);

  //  Copy precinct's row and column totals; function will whittle these.
  matrix_get_set_block(NNtots_temp_vec, 0, 0, 0, numcols(NNtots_temp_vec)-1, 
		       NNtots, num_p, num_p, 0, numcols(NNtots)-1);

  //  Copy the precinct's bounds for each cell; function will alter these.
  matrix_get_set_block(NNbounds_temp_vec, 0, 0, 0, numcols(NNbounds_temp_vec)-1, 
		       NNbounds, num_p, num_p, 0, numcols(NNbounds)-1);

  for (rr=0; rr<numrows_pt_m1; rr++){
    for (cc=0; cc<numcols_pt_m1; cc++){
      position = rr*numcols_pt + cc; //vector location, precinct table(rr, cc)

      /*  Get the range of numbers available and the lower limit.  */
      //printf("In the loop:  rr = %d, cc = %d.\n", rr, cc);
      range_NN = matrix_fast_get_element(NNbounds_temp_vec, 0,position+numcells_pt,nr_NNBTV) - 
	matrix_fast_get_element(NNbounds_temp_vec, 0,position,nr_NNBTV);
      if (range_NN==0.0) 
	continue;
      ret_val -= log(range_NN+1.0);
    
    //  Reset the precinct row total in light of the "drawn" value
    matrix_set_element(NNtots_temp_vec, 0, rr, 
		       matrix_get_element(NNtots_temp_vec, 0, rr) - 
		       matrix_get_element(NNs, num_p, position));

    //  Calculate the new lower bound and set it
    new_bound = matrix_get_element(NNtots_temp_vec, 0, rr);
    //printf("Calculating new lower bound.\n");
    for (ii=(cc+2); ii<numcols_pt; ii++){
      //printf("In the internal loop.  ii = %d.\n", ii);
      new_bound -= matrix_get_element(NNtots_temp_vec, 0, numrows_pt+ii);
    }
    matrix_fast_set_element(NNbounds_temp_vec, 0,position+1,nr_NNBTV, max(0.0, new_bound));

    //  Calculate the new upper bound and set it
    //printf("Calculating new upper bound.\n");
    new_bound = min(matrix_get_element(NNtots_temp_vec, 0, rr),		
		    matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc+1));
    matrix_fast_set_element(NNbounds_temp_vec, 0,position+1+numcells_pt,nr_NNBTV, new_bound);
    }

    //  Recalculate column totals given the filled in values for row rr
    //printf("Recalculating column totals.\n");
    for (cc=0; cc<numcols_pt; cc++){
      //printf("In this inner loop, cc = %d.\n", cc);
      new_tot = matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc);
      new_tot -= matrix_get_element(NNs, num_p, rr*numcols_pt+cc);
      matrix_set_element(NNtots_temp_vec, 0, numrows_pt+cc, new_tot);
    }

    //  Reset the bounds for the next row
    for (cc=0; cc<numcols_pt; cc++){
      //  Reset lower bound
      new_bound = matrix_get_element(NNtots_temp_vec, 0, rr+1);  //row sum
      for (ii=0; ii<numcols_pt; ii++){
	if (ii == cc) 
	  continue;
	new_bound -= matrix_get_element(NNtots_temp_vec, 0, numrows_pt+ii);
      }
      matrix_fast_set_element(NNbounds_temp_vec, 0,
			 (rr+1)*numcols_pt+cc,nr_NNBTV, 
			 max(0.0, new_bound));
      //  Reset the upper bound
      //printf("At the very bottom, resetting the upper bounds.\n");
      new_bound = min(matrix_get_element(NNtots_temp_vec, 0, rr+1), 
		      matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc));
      matrix_fast_set_element(NNbounds_temp_vec, 0,
			(rr+1)*numcols_pt+cc+numcells_pt,nr_NNBTV,
			new_bound);
    }
  }
  return ret_val;
}

double
log_NNs_multinomial_mh_ratio(Matrix * const curr_row,
														 Matrix * const prop_row,
														 Matrix * const multinomial_parameters,
                             Matrix const * const lfactorial_vector)
{
	double ret = 0.0;
	const index_t numcols_pt = numcols(curr_row);
	double N_ri_prop, N_ri_curr;
	index_t ii;

	for (ii=0; ii<numcols_pt; ii++){

		// The (N_{r'i}^{prop}-N_{r'i}^{curr})*log(\theta_{r'i}) bit:
		N_ri_prop = matrix_fast_get_element(prop_row,0,ii,1);
		N_ri_curr = matrix_fast_get_element(curr_row,0,ii,1);
		ret += (N_ri_prop-N_ri_curr)*log(matrix_fast_get_element(multinomial_parameters,0,ii,1));

		// The log(N_{r'i}^{curr}!/N_{r'i}^{prop}!) bit:
#ifdef _DBG_
  int dbg=0;
  if (dbg){
    Rprintf("lgamma(%g)=%g, lookup_table(%u)=%g\n",
            N_ri_curr+1.0,lgamma(N_ri_curr+1.0),
            (index_t)(N_ri_curr),
            matrix_fast_get_element_const(lfactorial_vector,0,(index_t)(N_ri_curr),1));
    Rprintf("lgamma(%g)=%g, lookup_table(%u)=%g\n",
            N_ri_prop+1.0,lgamma(N_ri_prop+1.0),
            (index_t)(N_ri_prop),
            matrix_fast_get_element_const(lfactorial_vector,0,(index_t)(N_ri_prop),1));
  }
#endif
		ret += matrix_fast_get_element_const(lfactorial_vector,0,(index_t)(N_ri_curr),1); // lgamma(N_ri_curr+1.0);
		ret -= matrix_fast_get_element_const(lfactorial_vector,0,(index_t)(N_ri_prop),1); // lgamma(N_ri_prop+1.0);
	}
	return ret;
}

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
                             Matrix const * const lfactorial_vector)
{
  const index_t nr_T = numrows(THETAS);
  double tmp_kk, mm_prop, mm_curr, tmp_val;
  double logmh = 0.0;
  index_t jj;

  for (jj=0; jj<numcols_pt; jj++){
    tmp_kk  = matrix_fast_get_element(tmp_KKs,special_row,jj,numrows_pt);
    mm_prop = matrix_fast_get_element(MMs_prop,special_row,jj,numrows_pt);
    mm_curr = matrix_fast_get_element(MMs_curr,special_row,jj,numrows_pt);
#ifdef _EP_DBG_
    const index_t lfaclen = numcols(lfactorial_vector);
    if ((mm_curr<0)||(mm_curr>(lfaclen-1)))
      error("Computing MH ratio, for mm_curr wanted to look up log(%u!), the table goes from 0 to %u\n",
            mm_curr,lfaclen);
    if ((mm_prop<0)||(mm_prop>(lfaclen-1)))
      error("Computing MH ratio, for mm_prop wanted to look up log(%u!), the table goes from 0 to %u\n",
            mm_prop,lfaclen);
#endif
    logmh += (mm_prop-mm_curr)*log(matrix_fast_get_element(THETAS,current_precinct,(special_row*numcols_pt)+jj,nr_T));
    logmh += matrix_fast_get_element_const(lfactorial_vector,0,(index_t)(mm_curr),1);
    logmh -= matrix_fast_get_element_const(lfactorial_vector,0,(index_t)(mm_prop),1);
  }

/*
  for (ii=0; ii<numrows_pt; ii++){
    for (jj=0; jj<numcols_pt; jj++){

      // Pick out the relevant cell elements:
      tmp_kk  = matrix_fast_get_element(tmp_KKs,ii,jj,numrows_pt);
      mm_prop = matrix_fast_get_element(MMs_prop,ii,jj,numrows_pt);
      mm_curr = matrix_fast_get_element(MMs_curr,ii,jj,numrows_pt);

      // Compute the pieces common to every cell:
      tmp_val = tmp_kk*(log(tmp_kk+mm_prop) - log(tmp_kk+mm_curr));
      tmp_val += matrix_fast_get_element_const(lfactorial_vector,0,(index_t)(mm_curr+tmp_kk),1);
      tmp_val -= matrix_fast_get_element_const(lfactorial_vector,0,(index_t)(mm_prop+tmp_kk),1);

      // See if it is the special row or not:
      if (ii==special_row){
 
        // Add in the special row piece:
        tmp_val += (mm_prop-mm_curr)*log(matrix_fast_get_element(THETAS,current_precinct,(ii*numcols_pt)+jj,nr_T));

      } else {

        // Add in the non-special row piece:
        tmp_val += matrix_fast_get_element_const(lfactorial_vector,0,(index_t)(mm_prop),1);
        tmp_val -= matrix_fast_get_element_const(lfactorial_vector,0,(index_t)(mm_curr),1);

      }

      // Add the cell contribution to the total:
      logmh += tmp_val;
    }  
  }
*/

  // All done, pass the log-MH ratio back.
  return logmh;
}


