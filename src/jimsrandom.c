#include "jimsrandom.h"

void
draw_NNs_indep_start(Matrix * const NNprop_vec, 
		     Matrix * const NNbounds, 
		     Matrix * const NNbounds_temp_vec, 
		     Matrix * const NNtots, 
		     Matrix * const NNtots_temp_vec, 
		     const index_t num_p, 
		     const index_t numrows_pt, 
		     const index_t numcols_pt, 
		     const index_t numcells_pt)
{
  //  A function used to find starting values of NNs.  This function
  //  sets the first R-1 by C-1 cells of the precinct table to roughly
  //  the midpoint of their bounds.  It calculates the remaining value
  //  deterministically.  It respects all bounds deterministically.
  //  User allocates all memory. NNprop_vec, NNbounds_temp_vec, NNtots_temp_vec
  //  must all be row vectors.

  const index_t numrows_pt_m1 = numrows_pt-1;
  const index_t numcols_pt_m1 = numcols_pt-1;
  
  long rr, cc, ii, position;
  double range_NN, lowlim, NN_prop = 0.0, new_bound, new_tot;

#ifdef _STARTING_STATE_DBG_
  const index_t dbg=1; 
  if (dbg)
    Rprintf("Copying NNtots and NNbounds states to *_temp_vec...\n");
#endif

  // Copy precinct's row and column totals; function will whittle these
  matrix_get_set_block(NNtots_temp_vec, 0, 0, 0, numcols(NNtots_temp_vec)-1, 
		       NNtots, num_p, num_p, 0, numcols(NNtots)-1);

  //Copy the precinct's bounds for each cell; function will alter these
  matrix_get_set_block(NNbounds_temp_vec, 0, 0, 0, numcols(NNbounds_temp_vec)-1, 
		       NNbounds, num_p, num_p, 0, numcols(NNbounds)-1);

#ifdef _STARTING_STATE_DBG_
  if (dbg)
    Rprintf("done. Starting the algorithm...\n");
#endif

  for (rr=0; rr<numrows_pt_m1; rr++){
    for (cc=0; cc<numcols_pt_m1; cc++){
      position = rr*numcols_pt + cc; //vector location, precinct table(rr, cc)

      //  Get the range of numbers available and the lower limit
      lowlim = matrix_get_element(NNbounds_temp_vec, 0, position);
      range_NN = matrix_get_element(NNbounds_temp_vec, 0, position+numcells_pt) - lowlim;
      if (range_NN == 0.0)
        NN_prop = lowlim;
      else 
        NN_prop = rint(range_NN/2.0) + lowlim;
      
      //assert(check_bounds(NN_prop, NNbounds_temp_vec, position, numcells_pt));
      matrix_set_element(NNprop_vec, 0, position, NN_prop);

      //  Reset precinct row total in light of NN_prop
      matrix_set_element(NNtots_temp_vec, 0, rr,			
			 matrix_get_element(NNtots_temp_vec, 0, rr) - NN_prop);

      //  Calculate the new lower bound and set it
      new_bound = matrix_get_element(NNtots_temp_vec, 0, rr);  // new row total
      for (ii=(cc+2); ii<numcols_pt; ii++){
        new_bound -= matrix_get_element(NNtots_temp_vec, 0, numrows_pt+ii);
      }
      matrix_set_element(NNbounds_temp_vec, 0, position+1, max(0.0, new_bound));

      //  Calcuate the new upper bound and set it
      new_bound = min(matrix_get_element(NNtots_temp_vec, 0, rr),	
		      matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc+1));
      matrix_set_element(NNbounds_temp_vec, 0, position+1+numcells_pt, new_bound);
    }
    //  Set precinct table(rr, C) deterministically
    matrix_set_element(NNprop_vec, 0, (rr+1)*numcols_pt-1, 
		       matrix_get_element(NNtots_temp_vec, 0, rr));

    //  Recalculate column totals given the filled in values for row rr
    for (cc=0; cc<numcols_pt; cc++){
      new_tot = matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc);
      new_tot -= matrix_get_element(NNprop_vec, 0, rr*numcols_pt+cc);
      matrix_set_element(NNtots_temp_vec, 0, numrows_pt+cc, new_tot);
    }

    //  Reset the bounds for the next row
    for (cc=0; cc<numcols_pt; cc++){
      //  Reset lower bound
      new_bound = matrix_get_element(NNtots_temp_vec, 0, rr+1);  //row sum
      for (ii=0; ii<numcols_pt; ii++){
        if (ii == cc) 
    	    continue;
        new_bound -= matrix_get_element(NNtots_temp_vec, 0, numrows_pt + ii);
      }
      matrix_set_element(NNbounds_temp_vec, 0,
			 (rr+1)*numcols_pt+cc, max(0.0, new_bound));
      //  Rest the upper bound
      new_bound = min(matrix_get_element(NNtots_temp_vec, 0, rr+1), 
		      matrix_get_element(NNtots_temp_vec, 0, numrows_pt + cc));
      matrix_set_element(NNbounds_temp_vec, 0, 
			 (rr+1)*numcols_pt+cc+numcells_pt, new_bound);
    }
  }

#ifdef _STARTING_STATE_DBG_
  if (dbg)
    Rprintf("done. Filling in final row deterministically and exiting...\n");
#endif

  //  Set precinct row R deterministically
  for (cc=0; cc<numcols_pt; cc++){
    matrix_set_element(NNprop_vec, 0, numrows_pt_m1*numcols_pt+cc,
		       matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc));
  }
  return;
}

void draw_NNs_MMs_indep_start(
         Matrix * const NNs,
         Matrix * const MMs,
         Matrix * const KKs,
         Matrix * const NNbounds,
         Matrix * const MMbounds,
         Matrix * const NNtots,
         Matrix * const MMtots, // Is this correct?
         Matrix * const NNbounds_temp_vec,
         Matrix * const NNtots_temp_vec,
         const index_t numrows_pt,
         const index_t numcols_pt,
         const index_t numcells_pt)
{
  index_t ii, jj;
  const index_t n=numrows(NNs);
#ifdef _STARTING_STATE_DBG_
  const index_t dbg=1;
  if (dbg)
    Rprintf("Setting up temporary MM storage...\n");
#endif
  // Maybe clean this up:
  Matrix * const MMprop_vec = rmatrix_new(1, numcells_pt);

#ifdef _STARTING_STATE_DBG_
  if (dbg)
    Rprintf("Beginning loop over each precinct...\n");
#endif

  // Draw a starting state for the MMs:
  for (ii=0; ii<n; ii++){

#ifdef _STARTING_STATE_DBG_
  if (dbg){
    Rprintf("Drawing starting state for NNs in precinct %u...\n",ii);
    Rprintf("dim(MMprop_vec)        = (%u x %u)\n",numrows(MMprop_vec),numcols(MMprop_vec));
    Rprintf("dim(MMbounds)          = (%u x %u)\n",numrows(MMbounds),numcols(MMbounds));
    Rprintf("dim(NNbounds)          = (%u x %u)\n",numrows(NNbounds),numcols(NNbounds));
    Rprintf("dim(NNbounds_temp_vec) = (%u x %u)\n",numrows(NNbounds_temp_vec),numcols(NNbounds_temp_vec));
    Rprintf("dim(MMtots)            = (%u x %u)\n",numrows(MMtots),numcols(MMtots));
    Rprintf("dim(NNtots_temp_vec)   = (%u x %u)\n",numrows(NNtots_temp_vec),numcols(NNtots_temp_vec));
    Rprintf("Precinct               = (%u)\n",ii);
    Rprintf("Table Dimensions       = (%u x %u)\n",numrows_pt,numcols_pt);
    Rprintf("Cells per table        = (%u)\n",numcells_pt);
  }
#endif

    draw_NNs_indep_start(MMprop_vec, MMbounds, NNbounds_temp_vec, MMtots, 
		 	 NNtots_temp_vec, ii, numrows_pt, numcols_pt, numcells_pt);

#ifdef _STARTING_STATE_DBG_
  if (dbg)
    Rprintf("done. Recording it in the MMs matrix...\n");
#endif

    // Just set the MMs directly:
    matrix_get_set_block(MMs,ii,ii,0,numcells_pt-1, 
			 MMprop_vec,0,0,0,numcells_pt-1);

#ifdef _STARTING_STATE_DBG_
  if (dbg)
    Rprintf("done. Setting the NNs matrix as MMs+KKs...\n");
#endif

    // Set the NNs from MMs and KKs:
    for (jj=0; jj<numcells_pt; jj++)
      matrix_fast_set_element(NNs,ii,jj,n, 
			  matrix_fast_get_element(MMprop_vec,0,jj,1) +
			  matrix_fast_get_element(KKs,ii,jj,n));

#ifdef _STARTING_STATE_DBG_
  if (dbg)
    Rprintf("done. Finished precinct %u.\n",ii);
#endif

  }
  return; 
}

/*
void draw_NNs_MMs_indep_start_old(Matrix * const NNs,
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
         const index_t numcells_pt)
{
  // Function to draw a starting state for the NNs, MMs matrices:

  // For simplicity:
  index_t ii, rr, cc, position, precinct_ii, jj, kk;
  double range_of_cell_count, lower_limit_of_cell_count, KK_tmp, new_bound;
  double MM_prop, NN_prop, new_tot;
  const index_t n = numrows(NNs);
  const index_t numrows_pt_m1 = numrows_pt-1;
  const index_t numcols_pt_m1 = numcols_pt-1;
  const index_t numcols_NNtots_m1 = numcols(NNtots)-1;
  const index_t numcols_NNbounds_m1 = numcols(NNbounds)-1;
  const index_t numcols_NNtots_tv_m1 = numcols(NNtots_temp_vec)-1;
  const index_t numcols_NNbounds_tv_m1 = numcols(NNbounds_temp_vec)-1;

  // Go precinct-by-precinct:
  for (precinct_ii=0; precinct_ii<n; precinct_ii++){

    // Copy precinct's row and column totals:
    matrix_get_set_block(NNtots_temp_vec, 0, 0, 0, numcols_NNtots_tv_m1, 
		         NNtots, precinct_ii, precinct_ii, 0, numcols_NNtots_m1);

   // Copy the precinct's bounds for each cell:
    matrix_get_set_block(NNbounds_temp_vec, 0, 0, 0, numcols_NNbounds_tv_m1, 
		         NNbounds, precinct_ii, precinct_ii, 0, numcols_NNbounds_m1);

   // Now, we go cell-by-cell:
   for (rr=0; rr<numrows_pt_m1; rr++){
     for (cc=0; cc<numcols_pt_m1; cc++){

       // Time for cell (rr,cc):
       position = rr*numcols_pt + cc;

       //  Get the range of MMs+KKs available and the lower limit:
       KK_tmp = matrix_fast_get_element(KKs,precinct_ii,position,n);

       // Lower limit excluding the exit poll:
       lower_limit_of_cell_count = matrix_get_element(NNbounds_temp_vec, 0, position);

       // Lower limit after accounting for the exit poll:
       lower_limit_of_cell_count -= KK_tmp;

       // Need to stop it being negative:
       if (lower_limit_of_cell_count<0.0)
         lower_limit_of_cell_count = 0.0;

       // The range without is the upper bound minus the lower limit:
       range_of_cell_count = 
         (matrix_get_element(NNbounds_temp_vec, 0, position+numcells_pt)-KK_tmp) - // (Upper limit on NN) - KK
         lower_limit_of_cell_count; // (Lower limit on MM)
     
       // Were there more counts in the exit poll than allowed?
       if (range_of_cell_count<0.0){
         Rprintf("\n\nInvalid exit poll value drawn in draw_NNs_MMs_indep_start():\n\n");
         Rprintf("In precinct %u I hit a problem at row %u, column %u:\n\n",
                 precinct_ii,rr,cc);
         Rprintf("Exit poll for the precinct is:\n");
         for (jj=0; jj<numrows_pt; jj++){ 
           for (kk=0; kk<numcols_pt; kk++)
             Rprintf("%g\t",matrix_get_element(KKs,precinct_ii,jj*numcols_pt+kk));
           Rprintf("\n");
         }
         Rprintf("(In progress) MMs for the precinct is:\n");
         for (jj=0; jj<numrows_pt; jj++){ 
           for (kk=0; kk<numcols_pt; kk++)
             Rprintf("%g\t",matrix_get_element(MMs,precinct_ii,jj*numcols_pt+kk));
           Rprintf("\n");
         }
         Rprintf("(In progress) NNs for the precinct is:\n");
         for (jj=0; jj<numrows_pt; jj++){ 
           for (kk=0; kk<numcols_pt; kk++)
             Rprintf("%g\t",matrix_get_element(NNs,precinct_ii,jj*numcols_pt+kk));
           Rprintf("\n");
         }
         Rprintf("NNtots (row totals):\n");
         for (jj=0; jj<numrows_pt; jj++)
           Rprintf("%g\t",matrix_get_element(NNtots,precinct_ii,jj));
         Rprintf("\n");
         Rprintf("NNtots (column totals):\n");
         for (jj=0; jj<numcols_pt; jj++)
           Rprintf("%g\t",matrix_get_element(NNtots,precinct_ii,numrows_pt+jj));
         Rprintf("\n");
         Rprintf("NNtots_temp_vec:\n");
         matrix_print_all(NNtots_temp_vec);
         Rprintf("NNbounds:\n");
         matrix_print_all(NNbounds_temp_vec);
         Rprintf("Upper bound for cell [%u,%u] = %g\n",
           rr,cc,matrix_get_element(NNbounds_temp_vec, 0, position+numcells_pt));
         Rprintf("Lower bound for cell [%u,%u] = %g\n",
           rr,cc,lower_limit_of_cell_count);
         Rprintf("Exit poll for cell   [%u,%u] = %g\n",rr,cc,KK_tmp);
         Rprintf("Hence, range of MM count is %g\n",range_of_cell_count);
         error("Invalid exit poll value in draw_NNs_MMs_indep_start()");
       }
       // Okay, now lets propose an MM and NN cell count based on the bounds:
       if (range_of_cell_count == 0.0)
         MM_prop = lower_limit_of_cell_count;
       else
         MM_prop = rint(range_of_cell_count/2.0) + lower_limit_of_cell_count;
       
       // The NN value is now determined:
       NN_prop =  MM_prop + KK_tmp;

       // Put the values in the MMs and NNs matrices:
       matrix_fast_set_element(MMs,precinct_ii,position,n, MM_prop);
       matrix_fast_set_element(NNs,precinct_ii,position,n, NN_prop);

       // Adjust the temporary precinct row totals:
       matrix_set_element(NNtots_temp_vec,0,rr,			
         matrix_get_element(NNtots_temp_vec,0,rr)-NN_prop);

        //  Calculate the new lower bound on the NNs and set it:
        new_bound = matrix_get_element(NNtots_temp_vec, 0, rr);  // new row total
        for (ii=(cc+2); ii<numcols_pt; ii++)
          new_bound -= matrix_get_element(NNtots_temp_vec, 0, numrows_pt+ii);
    
        matrix_set_element(NNbounds_temp_vec, 0, position+1, max(0.0, new_bound));

        //  Calcuate the new upper bound and set it:
        new_bound = min(matrix_get_element(NNtots_temp_vec, 0, rr),	
		      matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc+1));
        matrix_set_element(NNbounds_temp_vec, 0, position+1+numcells_pt, new_bound);

     } // End column loop
    
      //  Set cell counts in location (rr,C) deterministically:
      KK_tmp  = matrix_fast_get_element(KKs,precinct_ii,rr*numcols_pt + numcols_pt_m1,n);
      NN_prop = matrix_get_element(NNtots_temp_vec, 0, rr);
      MM_prop = NN_prop-KK_tmp;

      // Check the exit poll is valid:
      if (MM_prop<0){
         Rprintf("\n\nInvalid exit poll value drawn in draw_NNs_MMs_indep_start():\n\n");
         Rprintf("In precinct %u I hit a problem at row %u, column %u:\n\n",
                 precinct_ii,rr,cc);
         Rprintf("Exit poll for the precinct is:\n");
         for (jj=0; jj<numrows_pt; jj++){ 
           for (kk=0; kk<numcols_pt; kk++)
             Rprintf("%g\t",matrix_get_element(KKs,precinct_ii,jj*numcols_pt+kk));
           Rprintf("\n");
         }
         Rprintf("(In progress) MMs for the precinct is:\n");
         for (jj=0; jj<numrows_pt; jj++){ 
           for (kk=0; kk<numcols_pt; kk++)
             Rprintf("%g\t",matrix_get_element(MMs,precinct_ii,jj*numcols_pt+kk));
           Rprintf("\n");
         }
         Rprintf("(In progress) NNs for the precinct is:\n");
         for (jj=0; jj<numrows_pt; jj++){ 
           for (kk=0; kk<numcols_pt; kk++)
             Rprintf("%g\t",matrix_get_element(NNs,precinct_ii,jj*numcols_pt+kk));
           Rprintf("\n");
         }
         Rprintf("NNtots (row totals):\n");
         for (jj=0; jj<numrows_pt; jj++)
           Rprintf("%g\t",matrix_get_element(NNtots,precinct_ii,jj));
         Rprintf("\n");
         Rprintf("NNtots (column totals):\n");
         for (jj=0; jj<numcols_pt; jj++)
           Rprintf("%g\t",matrix_get_element(NNtots,precinct_ii,numrows_pt+jj));
         Rprintf("\n");
         Rprintf("NNtots_temp_vec:\n");
         matrix_print_all(NNtots_temp_vec);
         Rprintf("NNbounds:\n");
         matrix_print_all(NNbounds_temp_vec);
         Rprintf("Upper bound for cell [%u,%u] = %g\n",
           rr,cc,matrix_get_element(NNbounds_temp_vec, 0, position+numcells_pt));
         Rprintf("Lower bound for cell [%u,%u] = %g\n",
           rr,cc,lower_limit_of_cell_count);
         Rprintf("Exit poll for cell   [%u,%u] = %g\n",rr,cc,KK_tmp);
         Rprintf("Hence, range of MM count is %g\n",range_of_cell_count);
         error("Invalid exit poll value in draw_NNs_MMs_indep_start()");
      }

      // Fill in the values: 
      matrix_fast_set_element(MMs,precinct_ii,rr*numcols_pt + numcols_pt_m1,n, MM_prop);
      matrix_fast_set_element(NNs,precinct_ii,rr*numcols_pt + numcols_pt_m1,n, NN_prop);

      //  Recalculate column totals given the filled in values for row rr
      for (cc=0; cc<numcols_pt; cc++){
        new_tot = matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc);
        new_tot -= matrix_fast_get_element(NNs,precinct_ii,rr*numcols_pt+cc,n);
        matrix_set_element(NNtots_temp_vec, 0, numrows_pt+cc, new_tot);
      }

      //  Reset the bounds for the next row
      for (cc=0; cc<numcols_pt; cc++){
        //  Reset lower bound
        new_bound = matrix_get_element(NNtots_temp_vec, 0, rr+1);  //row sum
        for (ii=0; ii<numcols_pt; ii++){
          if (ii == cc) 
      	    continue;
          new_bound -= matrix_get_element(NNtots_temp_vec, 0, numrows_pt + ii);
        }
        matrix_set_element(NNbounds_temp_vec, 0,
		  	 (rr+1)*numcols_pt+cc, max(0.0, new_bound));

        //  Reset the upper bound
        new_bound = min(matrix_get_element(NNtots_temp_vec, 0, rr+1), 
	  	      matrix_get_element(NNtots_temp_vec, 0, numrows_pt + cc));
        matrix_set_element(NNbounds_temp_vec, 0, 
			   (rr+1)*numcols_pt+cc+numcells_pt, new_bound);
      }

   } // End row loop

    //  Set precinct row R deterministically:
    for (cc=0; cc<numcols_pt; cc++){
      KK_tmp = matrix_get_element(KKs,precinct_ii,numrows_pt_m1*numcols_pt+cc);
      NN_prop = matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc);
      MM_prop = NN_prop-KK_tmp;

      // Error check the exit poll entry:
      if (MM_prop<0.0){
         Rprintf("\n\nInvalid exit poll value drawn in draw_NNs_MMs_indep_start():\n\n");
         Rprintf("In precinct %u I hit a problem at row %u, column %u:\n\n",
                 precinct_ii,rr,cc);
         Rprintf("Exit poll for the precinct is:\n");
         for (jj=0; jj<numrows_pt; jj++){ 
           for (kk=0; kk<numcols_pt; kk++)
             Rprintf("%g\t",matrix_get_element(KKs,precinct_ii,jj*numcols_pt+kk));
           Rprintf("\n");
         }
         Rprintf("(In progress) MMs for the precinct is:\n");
         for (jj=0; jj<numrows_pt; jj++){ 
           for (kk=0; kk<numcols_pt; kk++)
             Rprintf("%g\t",matrix_get_element(MMs,precinct_ii,jj*numcols_pt+kk));
           Rprintf("\n");
         }
         Rprintf("(In progress) NNs for the precinct is:\n");
         for (jj=0; jj<numrows_pt; jj++){ 
           for (kk=0; kk<numcols_pt; kk++)
             Rprintf("%g\t",matrix_get_element(NNs,precinct_ii,jj*numcols_pt+kk));
           Rprintf("\n");
         }
         Rprintf("NNtots (row totals):\n");
         for (jj=0; jj<numrows_pt; jj++)
           Rprintf("%g\t",matrix_get_element(NNtots,precinct_ii,jj));
         Rprintf("\n");
         Rprintf("NNtots (column totals):\n");
         for (jj=0; jj<numcols_pt; jj++)
           Rprintf("%g\t",matrix_get_element(NNtots,precinct_ii,numrows_pt+jj));
         Rprintf("\n");
         Rprintf("NNtots_temp_vec:\n");
         matrix_print_all(NNtots_temp_vec);
         Rprintf("NNbounds:\n");
         matrix_print_all(NNbounds_temp_vec);
         Rprintf("Upper bound for cell [%u,%u] = %g\n",
           rr,cc,matrix_get_element(NNbounds_temp_vec, 0, position+numcells_pt));
         Rprintf("Lower bound for cell [%u,%u] = %g\n",
           rr,cc,lower_limit_of_cell_count);
         Rprintf("Exit poll for cell   [%u,%u] = %g\n",rr,cc,KK_tmp);
         Rprintf("Hence, range of MM count is %g\n",range_of_cell_count);
         error("Invalid exit poll value in draw_NNs_MMs_indep_start()");
      }

      // Fill in the values: 
      matrix_fast_set_element(MMs,precinct_ii,numrows_pt_m1*numcols_pt+cc,n, MM_prop);
      matrix_fast_set_element(NNs,precinct_ii,numrows_pt_m1*numcols_pt+cc,n, NN_prop);
    }

   } // End precinct loop

  return;
}
*/

void
draw_THETAS_from_NNs_start(Matrix * const THETAS, 
			   Matrix * const NNs, 
			   Matrix * const NNtots, 
			   const index_t numrows_pt,
			   const index_t numcols_pt)
{
  //  Takes a set of NNs and deterministically draws a set
  //  of THETAS.  The draw takes each precinct row of
  //  NNs and sets the precinct row of THETAS to the MLE
  //  of a multinomial.  That is, theta_rc is set to
  //  N_rc/N_r+.  There is a small adjustment at the end
  //  to assure that no theta_rc begins at 0.0, and that
  //  each precinct row's sum-to-one constraint is observed.

  static const double smallnum = .0001;
  long ii, jj, kk, position;
  const index_t numcols_pt_m1 = numcols_pt-1;
  const double one_over_numcols_pt = 1.0/numcols_pt; 
  index_t tmp_mult;
  long nr=numrows(THETAS);
  double theta_temp = 0.0, run_tot, NN_prec_row_tot, divisor;

  for (kk=0; kk<nr; kk++){
    for (ii=0; ii<numrows_pt; ii++){
      NN_prec_row_tot = matrix_get_element(NNtots, kk, ii);
      run_tot = 0.0;
      if (NN_prec_row_tot == 0.0){
				theta_temp = one_over_numcols_pt;
				tmp_mult = ii*numcols_pt;
				for (jj=0; jj<numcols_pt_m1; jj++){
					position = tmp_mult + jj;
					run_tot += theta_temp;
					matrix_set_element(THETAS, kk, position, theta_temp);
				}
      } else {
				divisor = 1.0 + numcols_pt*smallnum;
				tmp_mult = ii*numcols_pt;
				for (jj=0; jj<numcols_pt_m1; jj++){
					position = tmp_mult + jj;
					theta_temp = smallnum + matrix_get_element(NNs, kk, position)/NN_prec_row_tot;
					theta_temp = theta_temp/divisor;
					run_tot += theta_temp;
					matrix_set_element(THETAS, kk, position, theta_temp);
				}
      }
      matrix_set_element(THETAS, kk, (ii+1)*numcols_pt - 1, 1.0 - run_tot);
    }
  }
  return;
}

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
			    Matrix * const SIGMA_dims)
{
  //  Updates a matrix of THETAS.  At present, I contemplate
  //  using the draw_THETAS_t_dependent function most of the
  //  time, using the draw_THETAS_Dirichlet_independent
  //  function only once every LLth time, where the LL can
  //  be unique to the precinct.
  //  User allocates all memory.  Lots of error-checking done here.

  index_t ii;
  const index_t nr_T=numrows(THETAS);
  const index_t nr_TcuD=numrows(THETAS_count_use_Diri);
  const index_t nr_Tcut=numrows(THETAS_count_use_t);
  const index_t nr_uDev=numrows(use_Diri_every_vec);

#ifdef _DBG_
  int dbg=0;
  if (dbg){
    Rprintf("Inside draw_THETAS_t_and_Dirichlet\n");
    Rprintf("Calling matrix_inverse()\n");
  }
#endif
  matrix_inverse(SIGMA_cu,SIGMA_inverse_for_prop,SIGMA_dims);
#ifdef _DBG_
  if (dbg)
    Rprintf("Done. Starting loop...\n");
#endif

  for (ii=0; ii<nr_T; ii++){

#ifdef _DBG_
  if (dbg){
    Rprintf("Iteration number: %ld...\n",iternum);
    Rprintf("use_Diri_every_vec:\n");
    matrix_print_all(use_Diri_every_vec);
    Rprintf("I am going to access the first %u columns in this matrix (it has %u columns).\n",
            nr_T,numcols(use_Diri_every_vec));
  }
#endif

// NOTE:
// CHANGED AS OF 02/12/09 -- WAS WANTING TO ACCESS nr_T COLUMNS OF 
// use_Diri_every_vec WHEN IT IS JUST A SCALAR. DID IT USED TO BE 
// A VECTOR?

    if (fmod(iternum, matrix_fast_get_element(use_Diri_every_vec,0,ii,nr_uDev)) == 0.0){

#ifdef _DBG_
  if (dbg)
    Rprintf("Branched to Dirichelet proposal...\n");
#endif

      matrix_fast_increment_element(THETAS_count_use_Diri,0,ii,nr_TcuD, 1.0);

      // Used less frequently:
      draw_THETAS_Dirichlet_independent_one_row(THETAS, OMEGAS,
						prop_THETA, prop_OMEGA,	
						SIGMA_chol_cu, temp1_vec, 
						temp2_vec, NNs, 
						mu_vec_cu, SIGMA_cu,
						AUG, acc_THETAS_Diri_vec, 
						numrows_pt, numcols_pt,	
						numcells_pt, ii,
                                                SIGMA_inverse_for_prop,
                                                tmpMean, tmpOut, tmpScalar);
    } else {

#ifdef _DBG_
  if (dbg)
    Rprintf("Branched to t proposal...\n");
#endif

      matrix_fast_increment_element(THETAS_count_use_t,0,ii,nr_Tcut, 1.0);

      // Used most frequently:
      draw_THETAS_t_dependent_one_row(THETAS, OMEGAS,
				      prop_THETA, prop_OMEGA,
				      SIGMA_chol_cu, temp1_vec,
				      temp2_vec, NNs,
				      mu_vec_cu, SIGMA_cu,
				      AUG, acc_THETAS_t_vec,
				      rho_vec, SIGMA_chol_cu_temp,
				      dof, numrows_pt, numcols_pt,
				      numcells_pt, ii,
                                      SIGMA_inverse_for_prop,
                                      tmpMean, tmpOut, tmpScalar);
    }
  }
#ifdef _DBG_
  if (dbg)
    Rprintf("Finished loop.\n");
#endif
  return;
}

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
          Matrix * const tmp_mu_vec_cu)
{
  //  Updates a matrix of THETAS.  At present, I contemplate
  //  using the draw_THETAS_t_dependent function most of the
  //  time, using the draw_THETAS_Dirichlet_independent
  //  function only once every LLth time, where the LL can
  //  be unique to the precinct.
  //  User allocates all memory.  Lots of error-checking done here.

  index_t ii, jj;
  const index_t nr_T=numrows(THETAS);
  const index_t nr_TcuD=numrows(THETAS_count_use_Diri);
  const index_t nr_Tcut=numrows(THETAS_count_use_t);
  const index_t nr_uDev=numrows(use_Diri_every_vec);
  const index_t p = numcols(mu_mat_cu);

#ifdef _DBG_
  int dbg=0;
  if (dbg){
    Rprintf("Inside draw_THETAS_t_and_Dirichlet_with_covariates\n");
    Rprintf("Calling matrix_inverse()\n");
  }
#endif
  matrix_inverse(SIGMA_cu,SIGMA_inverse_for_prop,SIGMA_dims);
#ifdef _DBG_
  if (dbg)
    Rprintf("Done. Starting loop...\n");
#endif

  for (ii=0; ii<nr_T; ii++){

#ifdef _DBG_
  if (dbg){
    Rprintf("Iteration number: %ld...\n",iternum);
    Rprintf("use_Diri_every_vec:\n");
    matrix_print_all(use_Diri_every_vec);
    Rprintf("I am going to access the first %u columns in this matrix (it has %u columns).\n",
            nr_T,numcols(use_Diri_every_vec));
  }
#endif

    // Set-up the precinct-specific mu_i:
    for (jj=0; jj<p; jj++){
      matrix_fast_set_element(tmp_mu_vec_cu,0,jj,1, matrix_fast_get_element(mu_mat_cu,ii,jj,nr_T));
    }

    if (fmod(iternum, matrix_fast_get_element(use_Diri_every_vec,0,ii,nr_uDev)) == 0.0){

#ifdef _DBG_
  if (dbg)
    Rprintf("Branched to Dirichelet proposal...\n");
#endif

      matrix_fast_increment_element(THETAS_count_use_Diri,0,ii,nr_TcuD, 1.0);

      // Used less frequently:
      draw_THETAS_Dirichlet_independent_one_row(THETAS, OMEGAS,
						prop_THETA, prop_OMEGA,	
						SIGMA_chol_cu, temp1_vec, 
						temp2_vec, NNs, 
						tmp_mu_vec_cu, SIGMA_cu,
						AUG, acc_THETAS_Diri_vec, 
						numrows_pt, numcols_pt,	
						numcells_pt, ii,
            SIGMA_inverse_for_prop,
            tmpMean, tmpOut, tmpScalar);

    } else {

#ifdef _DBG_
  if (dbg)
    Rprintf("Branched to t proposal...\n");
#endif

      matrix_fast_increment_element(THETAS_count_use_t,0,ii,nr_Tcut, 1.0);

      // Same as without covariates, but mu_vec_cu is replaced by the 
      // precinct-specific mu_i.

      // Used most frequently:
      draw_THETAS_t_dependent_one_row(THETAS, OMEGAS,
				      prop_THETA, prop_OMEGA,
				      SIGMA_chol_cu, temp1_vec,
				      temp2_vec, NNs,
				      tmp_mu_vec_cu, SIGMA_cu, // mu_vec_cu is replaced by precinct-specific mu_i
				      AUG, acc_THETAS_t_vec,
				      rho_vec, SIGMA_chol_cu_temp,
				      dof, numrows_pt, numcols_pt,
				      numcells_pt, ii,
              SIGMA_inverse_for_prop,
              tmpMean, tmpOut, tmpScalar);
    }
  }
#ifdef _DBG_
  if (dbg)
    Rprintf("Finished loop.\n");
#endif
  return;
}

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
					  Matrix * const tmpScalar)
{
  //  Updates matrices of THETAS and OMEGAS with a new set of draws.
  //  For each precinct, the proposals for an M-H algorithm
  //  come from numrows_pt sets of independent Dirichlets.

  const index_t numcells_pt_m1 = numcells_pt-1; 
  const index_t OM_numcells_m1 = numcols(OMEGAS)-1;
  double MHjump;

#ifdef _DBG_
  int dbg=1;
  if (dbg){
    Rprintf("Inside draw_THETAS_Dirichlet_independent_one_row\n"); 
    Rprintf("Calling r_product_dirichelet...\n");
  }
#endif

  r_product_Dirichlet(prop_THETA, NNs, which_row,
			numrows_pt, numcols_pt);

#ifdef _DBG_
  if (dbg)
    Rprintf("Done. Computing the proposed states contribution to the MH-ratio...\n");
#endif

  MHjump = log_THETAS_proposal_product_Dirichlet(THETAS,which_row,
						   NNs, which_row);
  MHjump -= log_THETAS_proposal_product_Dirichlet(prop_THETA,0, 
						    NNs, which_row);

#ifdef _DBG_
  if (dbg)
    Rprintf("Done. Calling THETAs_to_OMEGAS...\n");
#endif

  THETAS_to_OMEGAS(prop_THETA, prop_OMEGA, numrows_pt, numcols_pt);

#ifdef _DBG_
  if (dbg)
    Rprintf("Done. Computing the current states contribution to the MH-ratio...\n");
#endif

  // TODO: Will need to change to allow precinct-specific mu_i's for covariates

  MHjump += log_p_target_theta_short(prop_THETA,0, 
				       prop_OMEGA,0,
				       which_row, NNs, mu_vec_cu,
				       SIGMA_cu, AUG, numrows_pt,
				       numcols_pt,
               SIGMA_inverse_for_prop,
               tmpMean, tmpOut, tmpScalar);

  MHjump -= log_p_target_theta_short(THETAS,which_row, 
				       OMEGAS,which_row,
				       which_row, NNs, mu_vec_cu,
				       SIGMA_cu, AUG, numrows_pt,
				       numcols_pt,
               SIGMA_inverse_for_prop,
               tmpMean, tmpOut, tmpScalar);

#ifdef _DBG_
  if (dbg)
    Rprintf("Done. Performing MH-ratio test...\n");
#endif

  if (log(runif(0.0, 1.0)) < MHjump){
    matrix_get_set_block(THETAS, which_row, which_row, 0, numcells_pt_m1,
		   prop_THETA, 0, 0, 0, numcells_pt_m1);
    matrix_get_set_block(OMEGAS, which_row, which_row, 0, OM_numcells_m1,
		   prop_OMEGA, 0, 0, 0, OM_numcells_m1);
    acc_THETAS_Diri_vec[which_row] += 1.0;
  }

#ifdef _DBG_
  if (dbg)
    Rprintf("Done with draw_THETAS_Dirichlet_independent_one_row\n");
#endif

  return;
}

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
          Matrix * const tmpScalar)
{
  //  Updates matrices of THETAS and OMEGAS with a new set of draws.
  //  For each precinct, the proposals for an M-H algorithm
  //  come from numrows_pt sets of independent Dirichlets.

  index_t ii;
  const index_t numcells_pt_m1=numcells_pt-1;
  const index_t OM_numcells_m1=numcols(OMEGAS)-1;
  const index_t nr=numrows(THETAS);
  double MHjump;

  for (ii=0; ii<nr; ii++){
    r_product_Dirichlet(prop_THETA, NNs, ii, numrows_pt, numcols_pt);

    MHjump  = log_THETAS_proposal_product_Dirichlet(THETAS,ii, NNs, ii);
    MHjump -= log_THETAS_proposal_product_Dirichlet(prop_THETA,0, NNs, ii);

    THETAS_to_OMEGAS(prop_THETA, prop_OMEGA, numrows_pt, numcols_pt);

    MHjump += log_p_target_theta_short(prop_THETA,0, prop_OMEGA,0,
				       ii, NNs, mu_vec_cu,
				       SIGMA_cu, AUG, numrows_pt,
				       numcols_pt,
                                       SIGMA_inverse_for_prop,
                                       tmpMean, tmpOut, tmpScalar);
    MHjump -= log_p_target_theta_short(THETAS,ii, OMEGAS,ii,
				       ii, NNs, mu_vec_cu,
				       SIGMA_cu, AUG, numrows_pt,
				       numcols_pt,
                                       SIGMA_inverse_for_prop,
                                       tmpMean, tmpOut, tmpScalar);

    if (log(runif(0.0, 1.0)) < MHjump){
      matrix_get_set_block(THETAS, ii, ii, 0, numcells_pt_m1,
			   prop_THETA, 0, 0, 0, numcells_pt_m1);
      matrix_get_set_block(OMEGAS, ii, ii, 0, OM_numcells_m1,
			   prop_OMEGA, 0, 0, 0, OM_numcells_m1);
      acc_THETAS_Diri_vec[ii] += 1.0;
    }
  }
  return;
}


void
r_product_Dirichlet(Matrix * const prop_THETA, 
		    Matrix * const NNs, 
		    const index_t num_p,
		    const index_t numrows_pt, 
		    const index_t numcols_pt)
{
  //  Sets prop_THETA, which must be a row vector, to a draw from
  //  numrows_pt independent Dirichlet distributions with
  //  parameters from the num_pth row of NNs.  Most
  //  error checking will be done in the layers above this function.
  //  NOTE:  Because some of the NNs values may be 0, the proposal
  //  distribution used is the value of the NNs + .1, to assure that
  //  a distribution exists even when the NNs values are 0.

  index_t ii, jj, position;
  double aa_total, alpha, aa_draw;
  const index_t nr_NNs = numrows(NNs);
  const index_t nr_PT = numrows(prop_THETA);
  index_t ii_mult_nc_pt;

  for (ii=0; ii<numrows_pt; ii++){
    aa_total = 0.0;
    for (jj=0; jj<numcols_pt; jj++){
      position = ii*numcols_pt + jj;
      alpha = matrix_fast_get_element(NNs, num_p,position,nr_NNs) + 0.1;
      aa_draw = rgamma(alpha, 1.0);
      aa_total += aa_draw;
      matrix_fast_set_element(prop_THETA, 0,position,nr_PT, aa_draw);
    }
    ii_mult_nc_pt = (ii*numcols_pt);
    for (jj=0; jj<numcols_pt; jj++){
      position = ii_mult_nc_pt+jj;
      matrix_fast_set_element(prop_THETA, 0,position,nr_PT,
			 matrix_fast_get_element(prop_THETA, 0,position,nr_PT)/aa_total);
    }
  }
  return;
}


void
draw_THETAS_t_dependent_one_row(Matrix * const THETAS, 
				Matrix * const OMEGAS,
				Matrix * const prop_THETA, 
				Matrix * const prop_OMEGA,	
				Matrix * const SIGMA_chol_cu, 
				Matrix * const temp1_vec, 
				Matrix * const temp2_vec, 
				Matrix * const NNs,		
				Matrix * const mu_vec_cu, // This is either global mu (no covariates), or a precinct-specific mu (covariates) 
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
				Matrix * const tmpScalar)
{
  //  Updates one of OMEGAS and THETAS with a new set of draws.  For
  //  the particular, the proposal for an M-H algorithm
  //  come from a multivariate t distribution in the OMEGA space.
  //  The new t draw is centered at the current value in
  //  the OMEGA space.

  const index_t numcols_pt_m1 = numcols_pt-1;
  const index_t numcells_pt_m1 = numcells_pt-1;
  const index_t OM_numcells_m1 = numcols(OMEGAS)-1;
  double MHjump, sqrt_rho_temp;

  sqrt_rho_temp = sqrt(matrix_get_element(rho_vec, 0, which_row));
  matrix_scalar_multiply(SIGMA_chol_cu, sqrt_rho_temp, SIGMA_chol_cu_temp);

  mvrt_c_chol(prop_OMEGA, OMEGAS,which_row, SIGMA_chol_cu_temp, dof,
              temp1_vec, temp2_vec);

  // When the last arg is 0 the only bit that it cares about is THETAS (not prop_*)
  MHjump = log_THETAS_proposal_t_jacobian(prop_OMEGA, prop_THETA,
					    THETAS, which_row, numrows_pt,
					    numcols_pt_m1, 0);
  MHjump -= log_THETAS_proposal_t_jacobian(prop_OMEGA, prop_THETA,
					     THETAS, which_row, numrows_pt,
					     numcols_pt_m1, 1);

  // Main source of time suck: 

  // Modify log_p_target_theta_short to allow covariates:

  MHjump += log_p_target_theta_short(prop_THETA,0, 
				       prop_OMEGA,0, 
				       which_row, NNs, mu_vec_cu, 
				       SIGMA_cu, AUG, numrows_pt, 
				       numcols_pt,
               SIGMA_inverse_for_prop,
               tmpMean, tmpOut, tmpScalar);

  MHjump -= log_p_target_theta_short(THETAS,which_row, 
				       OMEGAS,which_row, 
				       which_row, NNs, mu_vec_cu,  
				       SIGMA_cu, AUG, numrows_pt,  
				       numcols_pt,
               SIGMA_inverse_for_prop,
               tmpMean, tmpOut, tmpScalar);

  if (log(runif(0.0, 1.0))<MHjump){
      acc_THETAS_t_vec[which_row] += 1.0;
      matrix_get_set_block(THETAS, which_row,which_row,0,numcells_pt_m1,
			   prop_THETA, 0,0,0,numcells_pt_m1);
      matrix_get_set_block(OMEGAS, which_row,which_row,0,OM_numcells_m1,
			   prop_OMEGA, 0,0,0,OM_numcells_m1);
  }
  return;
}


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
			Matrix * const tmpScalar)
{
  //  Updates matrices of OMEGAS and THETAS with a new set of draws.  For
  //  each precinct, the proposals for an M-H algorithm
  //  come from a multivariate t distribution in the OMEGA space.
  //  The new t draw is centered at the current value in
  //  the OMEGA space.

  index_t ii;
  const index_t numcols_pt_m1=numcols_pt-1;
  const index_t nr=numrows(THETAS);
  const index_t numcells_pt_m1 = numcells_pt-1;
  const index_t OM_numcells_m1 = numcols(OMEGAS)-1;
  double MHjump, sqrt_rho_temp;

  for (ii=0; ii<nr; ii++){
    sqrt_rho_temp = sqrt(matrix_get_element(rho_vec, 0, ii));
    matrix_scalar_multiply(SIGMA_chol_cu, sqrt_rho_temp, SIGMA_chol_cu_temp);

    mvrt_c_chol(prop_OMEGA, OMEGAS,ii, SIGMA_chol_cu_temp, dof,
		temp1_vec, temp2_vec);

    MHjump = log_THETAS_proposal_t_jacobian(prop_OMEGA, prop_THETA,
					    THETAS, ii, numrows_pt,
					    numcols_pt_m1, 0);
    MHjump -= log_THETAS_proposal_t_jacobian(prop_OMEGA, prop_THETA,
					     THETAS, ii, numrows_pt,
					     numcols_pt_m1, 1);

    MHjump += log_p_target_theta_short(prop_THETA,0, prop_OMEGA,0, 
				       ii, NNs, mu_vec_cu,
				       SIGMA_cu, AUG, numrows_pt,
				       numcols_pt,
                                       SIGMA_inverse_for_prop,
                                       tmpMean, tmpOut, tmpScalar);
    MHjump -= log_p_target_theta_short(THETAS,ii, OMEGAS,ii,
				       ii, NNs, mu_vec_cu,
				       SIGMA_cu, AUG, numrows_pt,
				       numcols_pt,
                                       SIGMA_inverse_for_prop,
                                       tmpMean, tmpOut, tmpScalar);

    if(log(runif(0.0, 1.0)) < MHjump){
      acc_THETAS_t_vec[ii] += 1.0;
      matrix_get_set_block(THETAS, ii, ii, 0, numcells_pt_m1,
			   prop_THETA, 0, 0, 0, numcells_pt_m1);
      matrix_get_set_block(OMEGAS, ii, ii, 0, OM_numcells_m1,
			   prop_OMEGA, 0, 0, 0, OM_numcells_m1);
    }
  }
  return;
}



/*************************************************************************/


void
mvrnorm_c_chol(Matrix * const xx, 
	       Matrix * const mu_vec, 
	       Matrix * const SIGMA_chol,
	       Matrix * const temp1_vec, 
	       Matrix * const temp2_vec)
{
  //  Sets xx equal to numrows(xx) draws from the multivariate normal
  //  with mean mu and Var-Covar matrix with a Cholesky decomp of
  //  SIGMA_chol.  IMPORTANT:  the function works on ROWS, not COLUMNS,
  //  so if the user desires only one draw, xx should be a ROW vector.
  //  In this function (UNLIKE the general Greiner mvrnorm_c!), mu_vec
  //  must be a ROW vector.  User allocates memory to xx.

  const index_t dim=numcols(mu_vec);
  const index_t nr_xx = numrows(xx);
  const index_t nr_tv1 = numrows(temp1_vec);
  const index_t nc_tv1 = numcols(temp1_vec);
  index_t ii, jj;

  //Initialize temp1_vec and temp2_vec to 0, for safety; remove?
  for (ii=0; ii<nc_tv1; ii++){
    matrix_fast_set_element(temp1_vec,0,ii,nr_tv1, 0.0);
    matrix_fast_set_element(temp2_vec,0,ii,nr_tv1, 0.0);
  }

  for (ii=0; ii<nr_xx; ii++){
    for (jj=0; jj<dim; jj++){
      matrix_fast_set_element(temp1_vec,0,jj,nr_tv1, rnorm(0.0, 1.0));
    }
#ifdef _MAT_MULT_DBG_
    Rprintf("Calling matrix_multiply with args:\n");
    Rprintf("Argument 1:\n");
    matrix_print_all(temp1_vec);
    Rprintf("Argument 2:\n");
    matrix_print_all(SIGMA_chol);
    Rprintf("Argument 3:\n");
    matrix_print_all(temp2_vec);
#endif
    matrix_multiply(temp1_vec,'N',SIGMA_chol,'N',temp2_vec);
#ifdef _MAT_MULT_DBG_
    Rprintf("After matrix_multiply with args:\n");
    Rprintf("Argument 1:\n");
    matrix_print_all(temp1_vec);
    Rprintf("Argument 2:\n");
    matrix_print_all(SIGMA_chol);
    Rprintf("Argument 3:\n");
    matrix_print_all(temp2_vec);
#endif
    matrix_add(mu_vec, temp2_vec, temp1_vec);
    matrix_get_set_block(xx, ii,ii,0,dim-1, temp1_vec,0,0,0,dim-1);
  }
  return;
}

/*****************************************************************************/

void
rinvWis_c(Matrix * const xx, 
	  const double dof, 
	  Matrix * const SS,	
	  Matrix * const BB, 
	  Matrix * const CC, 
	  Matrix * const MM,
	  Matrix * const AUG_in_matrix_inverse)
{
  //  Sets xx equal to a draw from the inverse Wishart distribution
  //  with degrees of freedom dof (which need not be an integer) and scale
  //  parameter SS.  User allocates all memory.
  //  The parameterization of the inv-Wish is the following:
  //  p(xx) \propto det(xx)^{-\frac{1}{2}(dof + dim(xx) + 1)}
  //                x  exp{-\frac{1}{2} tr(SS solve(xx))}
  //  For speed, this function employs a Bartlett decomposition.
  //  The idea is to draw BB as an independent set of normals and
  // square roots of chi-squares, set CC = chol(SS), then set
  //  xx = (B^T)^{-1} C.  See J.L. Schafer, Analysis of Incomplete
  //  Multivariate Data 184 (1st ed. 1997) for the details.  Note
  //    that Schafer calls my SS LAMBDA^{-1}.

  long ii, jj, nrx = numrows(xx);

  //  Initialize BB, CC, and MM to 0 matrices
  matrix_fastset(BB, 0);
  matrix_fastset(CC, 0);
  matrix_fastset(MM, 0);
  const index_t nr_B = numrows(BB);

  for (ii=0; ii<nrx; ii++){
    matrix_fast_set_element(BB,ii,ii,nr_B, sqrt(rgamma((dof - ii - 1.0)/2.0, 2.0)));
    for (jj=ii+1; jj<nrx; jj++){
      matrix_fast_set_element(BB,ii,jj,nr_B, rnorm(0.0, 1.0));
    }
  }

  matrix_transpose(BB, CC);  // Sets CC = t(BB).
  
  // matrix_inverse sets AUG_in_matrix_inverse to zero matrix:
  matrix_inverse(CC, BB, AUG_in_matrix_inverse);  // Sets BB = solve(t(BB))

  matrix_cholesky(SS, CC);     // Sets CC = chol(SS)

#ifdef _MAT_MULT_DBG_
  Rprintf("Calling matrix_multiply with args:\n");
  Rprintf("Argument 1:\n");
  matrix_print_all(BB);
  Rprintf("Argument 2:\n");
  matrix_print_all(CC);
  Rprintf("Argument 3:\n");
  matrix_print_all(MM);
#endif
  matrix_multiply(BB,'N',CC,'N',MM); // Sets MM = BB %*% CC
#ifdef _MAT_MULT_DBG_
  Rprintf("After matrix_multiply:\n");
  Rprintf("Argument 1:\n");
  matrix_print_all(BB);
  Rprintf("Argument 2:\n");
  matrix_print_all(CC);
  Rprintf("Argument 3:\n");
  matrix_print_all(MM);
#endif

  matrix_transpose(MM, CC);    // Sets CC = t(MM)

#ifdef _MAT_MULT_DBG_
  Rprintf("Calling matrix_multiply with args:\n");
  Rprintf("Argument 1:\n");
  matrix_print_all(CC);
  Rprintf("Argument 2:\n");
  matrix_print_all(MM);
  Rprintf("Argument 3:\n");
  matrix_print_all(xx);
#endif
  matrix_multiply(CC,'N',MM,'N',xx); // Sets xx = CC %*% MM 
#ifdef _MAT_MULT_DBG_
  Rprintf("After matrix_multiply:\n");
  Rprintf("Argument 1:\n");
  matrix_print_all(CC);
  Rprintf("Argument 2:\n");
  matrix_print_all(MM);
  Rprintf("Argument 3:\n");
  matrix_print_all(xx);
#endif
  
  return;
}

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
	Matrix * const tvec2_dimSI)
{
#ifdef _DBG_MODE_
  Rprintf("Executing draw_mu()...\n");
#endif
  //   Draws from N(1/(num_of_p * kappa_0 +1)*(mu_vec_0 + num_of_p*kappa_0*omegabar),
  //                      (kappa_0/(num_of_p*kappa_0 + 1))*SIGMA_cu), sets
  //   mu_vec_cu equal to the result.  OMEGAS must be a matrix with each row
  //   representing a precinct.  User allocates all memory.

  const double denom = numrows(OMEGAS)*kappa_0 + 1.0;
  double tempsum;
  index_t ii, jj;

  //  Initialize the temporary matrices to 0, for safety; delete later?
  matrix_fastset(omegas_means_vec,0);
  matrix_fastset(mean_param_vec,0);
  matrix_fastset(SIGMA_chol_cu_temp,0);

  //  Obtain the sum of the OMEGAS, which is num_of_p * the OMEGAS means.
  const index_t nr_OMEGAS = numrows(OMEGAS);
  const index_t nc_OMEGAS = numcols(OMEGAS);
  const index_t nr_OMV = numrows(omegas_means_vec);
  for (jj=0; jj<nc_OMEGAS; jj++){
    tempsum = 0.0;
    for (ii=0; ii<nr_OMEGAS; ii++)
      tempsum += matrix_fast_get_element(OMEGAS, ii,jj,nr_OMEGAS);
    matrix_fast_set_element(omegas_means_vec, 0,jj,nr_OMV, tempsum);
  }
  
  //  Set the mean parameter
  matrix_scalar_multiply(omegas_means_vec, kappa_0, omegas_means_vec);
  matrix_add(mu_vec_0, omegas_means_vec, mean_param_vec);
  matrix_scalar_multiply(mean_param_vec, 1/denom, mean_param_vec);

  //  Set the chol decomp of the var-covar parameter
  matrix_scalar_multiply(SIGMA_chol_cu, sqrt(kappa_0/denom), SIGMA_chol_cu_temp);

  //  Take the draw
#ifdef _DBG_
  int dbg=0;
  if (dbg){
  Rprintf("Calling mvrnorm_c_chol with arguments:\n");
  Rprintf("mu_vec_cu:\n");
  matrix_print_all(mu_vec_cu);
  Rprintf("mean_param_vec:\n");
  matrix_print_all(mean_param_vec);
  Rprintf("SIGMA_chol_cu_temp:\n");
  matrix_print_all(SIGMA_chol_cu_temp);
  Rprintf("tvec1_dimSI:\n");
  matrix_print_all(tvec1_dimSI);
  Rprintf("tvec2_dimSI:\n");
  matrix_print_all(tvec2_dimSI);
  }
#endif
  mvrnorm_c_chol(mu_vec_cu, mean_param_vec, SIGMA_chol_cu_temp,
		 tvec1_dimSI, tvec2_dimSI);
  return;  
}

void
draw_eta_with_covariates(Matrix * const mu_mat_cu, 
  Matrix_p * const Xmatrix_list,
  Matrix * const eta_vec_cu,
	Matrix * const eta_vec_0,
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
  Matrix * const V_0_inverse)
{
#ifdef _DBG_ETA_DRAW_
  const int dbg=1;
  if (dbg)
    Rprintf("Executing draw_eta_with_covariates()...\n\n");
#endif

  //   Draws from:
  //              N( \mu = (V^{-1} + \sum_{i=1}^{n}Z_{i}^{T}\Sigma^{-1}Z_{i})^{-1}
  //                       (V^{-1}\eta_{0} + \sum_{i=1}^{n}Z_{i}^{T}\Sigma^{-1}\omega_{i}) ,
  //                \Sigma = (V^{-1} + \sum_{i=1}^{n}Z_{i}^{T}\Sigma^{-1}Z_{i})^{-1} )
  //              
  //   and sets eta_vec_cu equal to the result. 

  index_t ii, jj, kk, ll;
  const index_t nprecincts = numrows(OMEGAS);
  const index_t p = numcols(SIGMA_cu);
  const index_t d = numcols(V_0_inverse);
  double tmp = 0.0;

  //  Initialize the temporary matrices to 0, for safety; delete later?
  matrix_fastset(mean_vector_1_x_d,0);
  matrix_fastset(tmat1_d_x_d,0);
  matrix_fastset(tmat2_p_x_p,0);
  matrix_fastset(tmat3_p_x_d,0);
  matrix_fastset(tmat4_d_x_d,0);
  matrix_fastset(tmat5_d_x_p,0);
  matrix_fastset(tmat6_d_x_d,0);
  matrix_fastset(tvec1_1_x_d,0);
  matrix_fastset(tvec2_1_x_d,0);

  // Compute \Sigma^{-1}:
  matrix_inverse(SIGMA_cu,SIGMA_inverse,tmat2_p_x_p);

  // Temporary storage for the design matrices:
  Matrix *X_matrix = NULL;

  // Let tmat1_d_x_d hold the value of (V^{-1} + \sum_{i}Z_{i}^{T}\Sigma^{-1}Z_{i}):
  matrix_copy(V_0_inverse,tmat1_d_x_d);

  // Let tvec1_1_x_d hold the value of (V^{-1}\eta_{0} + \sum_{i}Z_{i}^{T}\Sigma^{-1}\omega_{i}):
  matrix_multiply(eta_vec_0,'N',V_0_inverse,'T',tvec1_1_x_d);

#ifdef _DBG_ETA_DRAW_
  if (dbg){
    Rprintf("SIGMA:\n\n");
    matrix_print_all(SIGMA_cu);
    Rprintf("\nSIGMA_inverse:\n\n");
    matrix_print_all(SIGMA_inverse);
    Rprintf("\nV_0_inverse:\n\n");
    matrix_print_all(V_0_inverse);
    Rprintf("\nV_0_inverse \%*\% eta_0 :\n\n");
    matrix_print_all(tvec1_1_x_d);
  }
#endif

  // Compute all terms that require summing over precincts:

  for (ii=0; ii<nprecincts; ii++){

    // Extract Z_i, the design matrix for precinct i:
    X_matrix = get_mat_p_ptr(Xmatrix_list,ii);

    // Add Z_{i}^{T}\Sigma^{-1}Z_{i} to tmat1:
    matrix_multiply(SIGMA_inverse,'N',X_matrix,'N',tmat3_p_x_d); // (p x p) %*% (p x d) = (p x d)
    matrix_multiply(X_matrix,'T',  tmat3_p_x_d,'N',tmat4_d_x_d); // (d x p) %*% (p x d) = (d x d)
    matrix_add(tmat4_d_x_d,tmat1_d_x_d,tmat1_d_x_d);             // (d x d)  +  (d x d) = (d x d)

    // Add Z_{i}^{T}\Sigma^{-1}\omega_{i} to tvec1_1_x_d:
    matrix_multiply(X_matrix,'T',SIGMA_inverse,'N',tmat5_d_x_p);
    for (jj=0; jj<d; jj++){
      tmp = 0.0;
      for (kk=0; kk<p; kk++){
        tmp += matrix_fast_get_element(tmat5_d_x_p,jj,kk,d)*
               matrix_fast_get_element(OMEGAS,ii,kk,nprecincts);
      }
      matrix_fast_increment_element(tvec1_1_x_d,0,jj,1, tmp);
    }

  }

  // Now invert the tmat1_d_x_d term: (V^{-1} + \sum_{i}Z_{i}^{T}\Sigma^{-1}Z_{i}):
  matrix_inverse(tmat1_d_x_d,tmat4_d_x_d,tmat6_d_x_d);
 
  // tmat4_d_x_d now contains the inverse term, we need its Cholesky decomposition
  // for the mvrnorm_c_chol function:
  matrix_cholesky(tmat4_d_x_d,tmat1_d_x_d);

#ifdef _DBG_ETA_DRAW_
  if (dbg){
    Rprintf("(V^{-1} + sum_{i}Z_{i}^{T}Sigma^{-1}Z_{i})^{-1}:\n\n");
    matrix_print_all(tmat4_d_x_d);
    Rprintf("\n(V^{-1}eta_{0} + sum_{i}Z_{i}^{T}Sigma^{-1}omega_{i}):\n\n");
    matrix_print_all(tvec1_1_x_d);
  }
#endif

  // tmat1_d_x_d now contains the upper-triangular Cholesky decompostion of 
  // the variance matrix. Lastly, we need the mean term:
  for (jj=0; jj<d; jj++){
    tmp = 0.0;
    for (kk=0; kk<d; kk++){
      tmp += matrix_fast_get_element(tmat4_d_x_d,jj,kk,d)*
             matrix_fast_get_element(tvec1_1_x_d,0,kk,1);
    }
    matrix_fast_set_element(mean_vector_1_x_d,0,jj,1, tmp);
  }
  // mean_vector now contains the mean term, 
  // the Cholesky decomposition of the variance term is in tmat1_d_x_d.

#ifdef _DBG_ETA_DRAW_
  if (dbg){
    Rprintf("\nMean of the multivariate normal draw:\n\n");
    matrix_print_all(mean_vector_1_x_d);
    Rprintf("\nCholesky decomposition of the covariance matrix:\n\n");
    matrix_print_all(tmat1_d_x_d);
  }
#endif 

  // Make the multivariate normal draw:
  mvrnorm_c_chol(eta_vec_cu, mean_vector_1_x_d, tmat1_d_x_d, // The draw, mean and Cholesky var
                 tvec1_1_x_d, tvec2_1_x_d);                  // two dummy vectors

#ifdef _DBG_ETA_DRAW_
  if (dbg){
    Rprintf("\nNew draw for eta:\n\n");
    matrix_print_all(eta_vec_cu);
    Rprintf("\n");
  }
#endif

  // Given the new draw for eta_vec_cu, we must recompute mu_mat_cu:
  multiply_list_of_X_by_eta(mu_mat_cu,Xmatrix_list,eta_vec_cu);

  return;  
}


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
	   const index_t dbg)
{
#ifdef _DBG_SIGMA_DRAW_
  if (dbg)
    Rprintf("Executing draw_SIGMA()...\n");
#endif
  //  Draws a new SIGMA matrix from the conditional posterior.
  //  User allocates all memory.  mu_vec_cu must be a row vector.
  //  OMEGAS must be allocated so that each row is a precinct.

  //  Calculate the SS matrix, except for the prior.  
  //  Note matrix_sum_xx_m_mu begins by making SS into a 0 matrix.
  matrix_sum_xx_m_mu(SS, OMEGAS, mu_vec_cu);

  //  Add the prior.
  matrix_add(PSI_0, SS, SS);

  //  Get the draw.  
  const index_t nr = numrows(OMEGAS);

#ifdef _DBG_SIGMA_DRAW_
  if (dbg){
    Rprintf("nu_0+nr = %g\n",nu_0+nr);
    Rprintf("SS:\n");
    matrix_print_all(SS);
    Rprintf("tmat2:\n");
    matrix_print_all(tmat2);
    Rprintf("tmat3:\n");
    matrix_print_all(tmat3);
    Rprintf("tmat4:\n");
    matrix_print_all(tmat4);
    Rprintf("AUG:\n");
    matrix_print_all(AUG_in_matrix_inverse);
    Rprintf("Drawing...");
  }
#endif
  rinvWis_c(SIGMA_cu, nu_0+nr, SS, tmat2, tmat3, tmat4, AUG_in_matrix_inverse);
#ifdef _DBG_SIGMA_DRAW_
  if (dbg){
    Rprintf("SIGMA:\n");
    matrix_print_all(SIGMA_cu);
  }
#endif
  return;
}

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
	   const index_t dbg)
{
#ifdef _DBG_SIGMA_DRAW_
  if (dbg)
    Rprintf("Executing draw_SIGMA_with_covariates()...\n");
#endif
  //  Draws a new SIGMA matrix from the conditional posterior.
  //  User allocates all memory.  mu_vec_cu must be a row vector.
  //  OMEGAS must be allocated so that each row is a precinct.

  //  Calculate the SS matrix, except for the prior.  Note matrix_sum_xx_m_mu
  //       begins by making SS into a 0 matrix.

  // The matrix multiplication is slightly different with covariates:
  matrix_sum_xx_m_mu_by_precinct(SS, OMEGAS, mu_mat_cu);

  //  Add the prior.
  matrix_add(PSI_0, SS, SS);

  // Degrees of freedom are prior df + number of precincts:
  const index_t nr = numrows(OMEGAS);

  //  Get the draw.  
  rinvWis_c(SIGMA_cu, nu_0+nr, SS, tmat2, tmat3, tmat4, AUG_in_matrix_inverse);
  
  return;
}


double
draw_NNs_prop_anywhere(Matrix * const NNprop_vec, 
		       Matrix * const NNbounds, 
		       Matrix * const NNbounds_temp_vec,
		       Matrix * const NNtots, 
		       Matrix * const NNtots_temp_vec, 
		       const index_t num_p, 	
		       const index_t numrows_pt, 
		       const index_t numcols_pt, 
		       const index_t numcells_pt)
{
  /*  Sets NNprop_vec to a proposal of a new precinct table NNs drawn "at random,"
      i.e., without regard to the current state of the precinct table NNs.
      The idea is to prevent the MCMC algorithm from getting caught in a local
      mode for the NNs.  This new proposal deterministically respects the bounds.
      It fills in the first C-1 columns of each row and calculates row C
      deterministically.  It does for the first R-1 rows and calculates all
      of row R deterministically.
      The return value is the log of the probability of drawing the new 
      precint table NNs that were actually drawn.
      User allocates all memory.  Error checking to be done in function that
      calls this one.  NNprop_vec, NNbounds_temp_vec, NNtots_temp_vec
      must all be row vectors.  */

  long rr, cc, ii, numrows_pt_m1=numrows_pt-1;
  long numcols_pt_m1=numcols_pt-1, position, stop_test;
  double range_NN, lowlim, NN_prop = 0.0, ret_val = 0.0, new_bound, new_tot;

    /*  Copy precinct's row and column totals; function will whittle these. */
  matrix_get_set_block(NNtots_temp_vec, 0, 0, 0, numcols(NNtots_temp_vec)-1, 
		       NNtots, num_p, num_p, 0, numcols(NNtots)-1);

  /*  Copy the precinct's bounds for each cell; function will alter these. */
  matrix_get_set_block(NNbounds_temp_vec, 0, 0, 0, numcols(NNbounds_temp_vec)-1, 
		       NNbounds, num_p, num_p, 0, numcols(NNbounds)-1);

  for (rr=0; rr<numrows_pt_m1; rr++){
    for (cc=0; cc<numcols_pt_m1; cc++){
      position = rr*numcols_pt + cc; //vector location, precinct table(rr, cc)

      /*  Get the range of numbers available and the lower limit */
      lowlim = matrix_get_element(NNbounds_temp_vec, 0, position);
      range_NN = matrix_get_element(NNbounds_temp_vec, 0, position+numcells_pt) - lowlim;
      if (range_NN == 0.0){
	NN_prop = lowlim;
      }
      else {
	stop_test = 1;
	while(stop_test){
	  NN_prop = rint(runif(-.5, range_NN+.5));
	  if((0 <= NN_prop) && (NN_prop <= range_NN)){
	    stop_test = 0;
	    NN_prop += lowlim;
	  }
	}
	ret_val -= log(range_NN + 1.0);
      }
      matrix_set_element(NNprop_vec, 0, position, NN_prop);

      //  Reset precinct row total in light of NN_prop
      matrix_set_element(NNtots_temp_vec, 0, rr,			
			 matrix_get_element(NNtots_temp_vec, 0, rr) - NN_prop);

      //  Calculate the new lower bound and set it
      new_bound = matrix_get_element(NNtots_temp_vec, 0, rr);  // new row total
      for (ii=(cc+2); ii<numcols_pt; ii++){
	new_bound -= matrix_get_element(NNtots_temp_vec, 0, numrows_pt+ii);
      }
      matrix_set_element(NNbounds_temp_vec, 0, position+1, max(0.0, new_bound));

      //  Calcuate the new upper bound and set it
      new_bound = min(matrix_get_element(NNtots_temp_vec, 0, rr),	
		      matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc+1));
      matrix_set_element(NNbounds_temp_vec, 0, position+1+numcells_pt, new_bound);
    }
    //  Set precinct table(rr, C) deterministically
    matrix_set_element(NNprop_vec, 0, (rr+1)*numcols_pt-1, 
		       matrix_get_element(NNtots_temp_vec, 0, rr));

    //  Recalculate column totals given the filled in values for row rr
    for (cc=0; cc<numcols_pt; cc++){
      new_tot = matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc);
      new_tot -= matrix_get_element(NNprop_vec, 0, rr*numcols_pt+cc);
      matrix_set_element(NNtots_temp_vec, 0, numrows_pt+cc, new_tot);
    }

    //  Reset the bounds for the next row
    for (cc=0; cc<numcols_pt; cc++){
      //  Reset lower bound
      new_bound = matrix_get_element(NNtots_temp_vec, 0, rr+1);  //row sum
      for (ii=0; ii<numcols_pt; ii++){
	if (ii == cc) 
	  continue;
	new_bound -= matrix_get_element(NNtots_temp_vec, 0, numrows_pt + ii);
      }
      matrix_set_element(NNbounds_temp_vec, 0,
			 (rr+1)*numcols_pt+cc, max(0.0, new_bound));
      //  Rest the upper bound
      new_bound = min(matrix_get_element(NNtots_temp_vec, 0, rr+1), 
		      matrix_get_element(NNtots_temp_vec, 0, numrows_pt + cc));
      matrix_set_element(NNbounds_temp_vec, 0, 
			 (rr+1)*numcols_pt+cc+numcells_pt, new_bound);
    }
  }


  //  Set precinct row R deterministically
  for (cc=0; cc<numcols_pt; cc++){
    matrix_set_element(NNprop_vec, 0, numrows_pt_m1*numcols_pt+cc,
		       matrix_get_element(NNtots_temp_vec, 0, numrows_pt+cc));
  }

  return ret_val;
}

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
		  const index_t numcells_pt)
{
  //  Updates a matrix of NNs with new draws from the conditional posterior
  //  using a M-H algorithm with the anywhere method of drawing a proposal.
  //  User allocates all memory.  Anything ending in _vec must be a row
  //  vector.  All other matrices must be organized one precinct per row.
  double aa_temp;
  const index_t numcells_pt_m1=(numcells_pt-1);
  const index_t nr_NNs=numrows(NNs);
  index_t ii;
  for (ii=0; ii<nr_NNs; ii++){
    aa_temp = log_p_NNs_prop_anywhere(NNs, NNbounds, NNbounds_temp_vec,
				      NNtots, NNtots_temp_vec, ii,
				      numrows_pt, numcols_pt, numcells_pt);
    aa_temp -= draw_NNs_prop_anywhere(NNprop_vec, NNbounds, NNbounds_temp_vec,
				      NNtots, NNtots_temp_vec, ii,
				      numrows_pt, numcols_pt, numcells_pt);

    aa_temp += log_p_target_NNs(NNprop_vec,0, ii, THETAS, numcells_pt);
    aa_temp -= log_p_target_NNs(NNs,ii, ii, THETAS, numcells_pt);

    if (log(runif(0.0,1.0))<aa_temp)
      matrix_get_set_block(NNs,ii,ii,0,numcells_pt_m1, NNprop_vec,0,0,0,numcells_pt_m1);
  }
  return;
}

void
rGibbsNNs(Matrix * const NNs, 
	  const index_t num_p, 
	  Matrix * const THETAS, 
	  Matrix_int * const which_rc_int, 
	  double * const ff_vec, 
	  const index_t whichperm, 
	  const index_t numrows_pt, 
  	  const index_t numcols_pt)
{
  //  Updates a row of NNs by using two rows and two columns
  //  specified in the whichpermth row of which_rc_int, collapsing the resulting 4 cells
  //  into a two-by-two contingency table, then draws from
  //  a non-central hypergeometric.

  // Pick the rows and columns
  const index_t nr_WRCI = numrows_int(which_rc_int);
  const index_t r1 = matrix_fast_get_int_element(which_rc_int, whichperm, 0, nr_WRCI);
  const index_t r2 = matrix_fast_get_int_element(which_rc_int, whichperm, 1, nr_WRCI);
  const index_t c1 = matrix_fast_get_int_element(which_rc_int, whichperm, 2, nr_WRCI);
  const index_t c2 = matrix_fast_get_int_element(which_rc_int, whichperm, 3, nr_WRCI);

  // Set the corresponding positions
  const index_t pos1 = r1*numcols_pt + c1;
  const index_t pos2 = pos1+c2-c1;
  const index_t pos3 = r2*numcols_pt + c1;
  const index_t pos4 = pos3+c2-c1;

  //  Get the current NNs
  const index_t nr_NNs = numrows(NNs);
  const index_t nr_THETAS = numrows(THETAS);
  double val1 = matrix_fast_get_element(NNs, num_p, pos1, nr_NNs);
  const double val2 = matrix_fast_get_element(NNs, num_p, pos2, nr_NNs);

  //  Set the params of the noncentral hypergeometric draw.
  const double m1 = val1 + val2;
  const double n1 = val1 + matrix_fast_get_element(NNs, num_p, pos3, nr_NNs);
  const double n2 = val2 + matrix_fast_get_element(NNs, num_p, pos4, nr_NNs);
  const double psi = 
	  (matrix_fast_get_element(THETAS, num_p, pos1, nr_THETAS)* 
	   matrix_fast_get_element(THETAS, num_p, pos4, nr_THETAS))/
	  (matrix_fast_get_element(THETAS, num_p, pos2, nr_THETAS)* 
	   matrix_fast_get_element(THETAS, num_p, pos3, nr_THETAS));
  //  Decide whether to do the draw at all
  if (m1!=(n1+n2) && m1!=0.0 && n1!=0.0 && n2!=0){

    //  Draw the new val1
    val1 = rnoncenhypgeo(m1, n1, n2, psi, DBL_EPSILON, ff_vec, 0);

    //  Set the new values
    matrix_fast_set_element(NNs, num_p,pos1,nr_NNs, val1);
    matrix_fast_set_element(NNs, num_p,pos2,nr_NNs, m1-val1);
    matrix_fast_set_element(NNs, num_p,pos3,nr_NNs, n1-val1);
    matrix_fast_set_element(NNs, num_p,pos4,nr_NNs, n2-m1+val1);
  }
  return;
}

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
   Matrix const * const lfactorial_vector) // Log-factorial lookup table
{
  //  Updates a matrix of NNs with the appropriate distribution.  Once
  //  in every nolocalcalmode times it is called, this function uses
  //  the MH algorithm with an anywhere-within-the-bounds proposal
  //  distribution.  In all other calls, it executes 2 x 2 Gibbs updates
  //  num_scans times for each precinct.
  //  All matrices must be organized one-row-per precinct.

  index_t ii, jj, kk, whichpermnow;
  const index_t nr_NNs = numrows(NNs);
  const index_t nr_OVI = numrows_int(ordervec_int);
  const index_t nc_OVI = numcols_int(ordervec_int);
	const index_t numrows_pt_m1 = (numrows_pt-1);
	double tmp_u, tmp_c;
	int nreps;
	index_t special_row;

  if (fmod(iternum, nolocalmode) != 0.0){  // regular situation, use Gibbs
    // May want to randomize the ordering of the draws, i.e.,
    // Randomly pick an index from whichpermnow
    for (ii=0; ii<nr_NNs; ii++){
      for (kk=0; kk<num_scans; kk++){

				// (R-1)-product multinomial update:

				// Select the deterministic row:
				tmp_u = runif(0.0,1.0);
				tmp_c = 0.0;
				for (jj=0; jj<numrows_pt; jj++){
					tmp_c += matrix_fast_get_element(sr_probs,ii,jj,nr_NNs);
					if (tmp_u<tmp_c){
						special_row = jj;
						break;
					}
				}
				//Rprintf("For precinct %u: Selected special row %u\n",ii,special_row);
				//special_row = runif_index_t(0,numrows_pt_m1);

				// Retrieve the number of reps desired:
				nreps = matrix_fast_get_int_element(sr_reps,ii,special_row,nr_NNs);
				//Rprintf("To be repeated %u times.\n",nreps); 

				for (jj=0; jj<nreps; jj++)
					draw_NNs_multinomial_MH(NNs,NNtots,ii,special_row,THETAS,numrows_pt,
						numcols_pt,NNs_prop,multinomial_parameters,curr_row,prop_row,
						vld_mv_p,acc_mv_p,NNs_count_use_multinom,lfactorial_vector);

				// NCHG-update:
				if (nreps==0){
					//Rprintf("Reverting back to nchg update\n");
					for (jj=0; jj<nc_OVI; jj++){
						whichpermnow = matrix_fast_get_int_element(ordervec_int, 0,jj,nr_OVI);
						rGibbsNNs(NNs, ii, THETAS, which_rc_int, ff_vec, whichpermnow, numrows_pt, numcols_pt);
					}
				}

      }
      //for (kk=0; kk<num_Gibbs; kk++)
			//	rGibbsNNs_randomrc(NNs, ii, THETAS, ff_vec, numrows_pt, numcols_pt);
    }

  } else {  
		//  Use the draw_NNs_anywhere as a local mode check
    draw_NNs_anywhere(NNs, NNprop_vec, NNbounds,
		      NNbounds_temp_vec, NNtots,
		      NNtots_temp_vec, THETAS,
		      numrows_pt, numcols_pt, numcells_pt);
  }
  return;
}

// Exit poll version:
void 
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
   Matrix const * const lfactorial_vector) // Log-factorial lookup table
{
  // Some preliminaries:
  int nreps;
  double tmp_u, tmp_c, alpha, log_MH_ratio;
  const index_t n=numrows(NNs);
  index_t ii, jj, kk, special_row;

  // Go precinct by precinct:
  for (ii=0; ii<n; ii++){

      // Repeat num_scans times for each precinct:
      for (kk=0; kk<num_scans; kk++){

				// (R-1)-product multinomial update:

				// Select the deterministic row:
				tmp_u = runif(0.0,1.0);
				tmp_c = 0.0;
				for (jj=0; jj<numrows_pt; jj++){
					tmp_c += matrix_fast_get_element(sr_probs,ii,jj,n);
					if (tmp_u<tmp_c){
						special_row = jj;
						break;
					}
				}
				//Rprintf("For precinct %u: Selected special row %u\n",ii,special_row);
				//special_row = runif_index_t(0,numrows_pt_m1);

				// Retrieve the number of reps desired:
				nreps = matrix_fast_get_int_element(sr_reps,ii,special_row,n);
				//Rprintf("To be repeated %u times.\n",nreps);

				for (jj=0; jj<nreps; jj++)
					draw_MMs_multinomial_MH(NNs,MMs,KKs,NNtots,MMtots,KKtots,
            ii,special_row,THETAS,numrows_pt,
						numcols_pt,NNs_prop,MMs_prop,NNs_curr,MMs_curr,tmp_KKs,
            multinomial_parameters,curr_row,prop_row,
						vld_mv_p,acc_mv_p,NNs_count_use_multinom,lfactorial_vector);
         
      } // End kk^th scan
  } // End precinct ii
  
  return;
}

void
mvrt_c_chol(Matrix * const xx, 
	    Matrix * const mu_mat, 
	    const index_t mu_row,
	    Matrix * const SIGMA_chol_cu, 
	    const double dof,
	    Matrix * const temp1_vec, 
	    Matrix * const temp2_vec)
{
  //  Sets xx equal to numrows(xx) draws from the multivariate t
  //  with mean mu, Var-Covar matrix with a Cholesky decomp of
  //  SIGMA_chol_cu, and dof degrees of freedom.  IMPORTANT:  the
  //  function works on ROWS, not COLUMNS,
  //  so if the user desires only one draw, xx should be a ROW vector.
  //  In this function (UNLIKE the general Greiner mvrnorm_c!), mu_vec
  //  must be a ROW vector.  User allocates memory to xx.
  //  Because this function will be called so often, user must allocate
  //  memory to temp1_vec and temp2_vec, which are used in the function.
  //  In addition, to avoid a lot of copying, this function takes
  //  advantage of the Greiner matrix structure in which matrices are
  //  stored by row. 
      
  // Modified (02/22/08): Now mu is passed as a matrix, with mu_row
  // specifying the row of the matrix to use as the mean.


  double cc;
  index_t ii, jj;
  const index_t dim_m1=numcols(xx)-1;
  const index_t nrx=numrows(xx);
  const index_t ncx=numcols(xx);
  const index_t nr_t2=numrows(temp2_vec);
  const index_t nc_t2=numcols(temp2_vec);
  const index_t nr_mu=numrows(mu_mat);
  const double dof_by_2 = dof/2.0;

  for (ii=0; ii<nrx; ii++){
    cc = sqrt(dof/rgamma(dof_by_2,2.0));
    for (jj=0; jj<ncx; jj++)
      matrix_fast_set_element(temp1_vec,0,jj,nrx, rnorm(0.0, 1.0));
    
#ifdef _MAT_MULT_DBG_
    Rprintf("Calling matrix_multiply with args:\n");
    Rprintf("Argument 1:\n");
    matrix_print_all(temp1_vec);
    Rprintf("Argument 2:\n");
    matrix_print_all(SIGMA_chol_cu);
    Rprintf("Argument 3:\n");
    matrix_print_all(temp2_vec);
#endif
    matrix_multiply(temp1_vec,'N',SIGMA_chol_cu,'N',temp2_vec);
#ifdef _MAT_MULT_DBG_
    Rprintf("After matrix_multiply:\n");
    Rprintf("Argument 1:\n");
    matrix_print_all(temp1_vec);
    Rprintf("Argument 2:\n");
    matrix_print_all(SIGMA_chol_cu);
    Rprintf("Argument 3:\n");
    matrix_print_all(temp2_vec);
#endif
    matrix_scalar_multiply(temp2_vec, cc, temp2_vec);
    for (jj=0; jj<nc_t2; jj++)
      matrix_fast_increment_element(temp2_vec, 0,jj,nr_t2, 
			 matrix_fast_get_element(mu_mat,mu_row,jj,nr_mu));
    matrix_get_set_block(xx,ii,ii,0,dim_m1, temp2_vec,0,0,0,dim_m1);
  }
  return;
}

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
                        Matrix const * const lfactorial_vector) // Log-factorial lookup table
{
	// Do a Product-Multinomial MH move:
	// (1) Fix special_row, and propose multinomials with parameters
	// given by the current state of theta for precinct ii.

#ifdef _DBG_
  int dbg=0;
#endif

	index_t ii, jj, kk, tmp_mult;
	const index_t nr_mp = numrows(multinomial_parameters);
	const index_t nr_pr = numrows(prop_row);
	const index_t nr_NNs = numrows(NNs);
	const index_t nr_THETAS = numrows(THETAS);
	const index_t nr_NNtots = numrows(NNtots);
	double log_mh_ratio = 0.0;
	double NN_prec_row_tot, NN_prec_col_tot;

	// Increment the counter:
	matrix_fast_increment_element(NNs_count_use_multinom,current_precinct,special_row,nr_NNs, 1);

	// Find the parameters for the multinomial for the remaining rows:
	for (jj=0; jj<numrows_pt; jj++){

		// Skip the special row:
    if (jj==special_row)
			continue;

		// Retrieve the row total for row jj:
    NN_prec_row_tot = matrix_fast_get_element(NNtots,current_precinct,jj,nr_NNtots);

		// Save a bit of computation:
		tmp_mult = (jj*numcols_pt);

		// Set the proposed state for that row:
		for (kk=0; kk<numcols_pt; kk++)
			matrix_fast_set_element(multinomial_parameters,0,kk,nr_mp, 
											matrix_fast_get_element(THETAS,current_precinct,tmp_mult+kk,nr_THETAS));

		// Now, make the multinomial draw:
		rmultinomial(prop_row,multinomial_parameters,NN_prec_row_tot);

		// Put the draws into matrix form:
		for (kk=0; kk<numcols_pt; kk++)
			matrix_fast_set_element(NNs_prop,jj,kk,numrows_pt, matrix_fast_get_element(prop_row,0,kk,nr_pr));
		
	}

	// Determine if the proposal is valid:
	for (kk=0; kk<numcols_pt; kk++){
		NN_prec_col_tot = matrix_fast_get_element(NNtots,current_precinct,numrows_pt+kk,nr_NNtots);
		// See how many counts are left in each column:
		for (jj=0; jj<numrows_pt; jj++)
			if (jj!=special_row)
				NN_prec_col_tot -= matrix_fast_get_element(NNs_prop,jj,kk,numrows_pt);

		// Was it valid for this column?
		if (NN_prec_col_tot<0){
			// No. :(
			return;
		}
		// Yes! :)

		tmp_mult = (special_row*numcols_pt);

		// Set up the proposed values:
		matrix_fast_set_element(prop_row,0,kk,1, NN_prec_col_tot);
		matrix_fast_set_element(curr_row,0,kk,1, 
										matrix_fast_get_element(NNs,current_precinct,tmp_mult+kk,nr_NNs));
		matrix_fast_set_element(multinomial_parameters,0,kk,1,
										matrix_fast_get_element(THETAS,current_precinct,tmp_mult+kk,nr_THETAS));

		matrix_fast_set_element(NNs_prop,special_row,kk,numrows_pt, NN_prec_col_tot);
	}

	// We have a valid proposal for the whole table. :)

	// Increment the number of valid proposals for this precinct/special row:
  vld_mv_p[current_precinct+(nr_NNtots*special_row)]++;

	// Compute the log(MH-ratio):
	log_mh_ratio = log_NNs_multinomial_mh_ratio(curr_row,prop_row,multinomial_parameters,
                                              lfactorial_vector);

	if (log(runif(0.0,1.0))<log_mh_ratio){
		// ACCEPTED:
#ifdef _DBG_
		if (dbg>1){
			Rprintf("\nMultinomial proposal accepted\n");
		}
#endif

		// Increment the number of accepted proposals for this precinct/special row:
  	acc_mv_p[current_precinct+(nr_NNtots*special_row)]++;

		// Update the current state accordingly:
		for (ii=0; ii<numrows_pt; ii++){
			tmp_mult = (ii*numcols_pt);
			for (kk=0; kk<numcols_pt; kk++)
				matrix_fast_set_element(NNs,current_precinct,tmp_mult+kk,nr_NNs,
					matrix_fast_get_element(NNs_prop,ii,kk,numrows_pt));
		}

	} else {
		// REJECTED:
#ifdef _DBG_
		if (dbg>1){
			Rprintf("\nMultinomial proposal rejected\n");
		}
#endif
		
	}
	return;
}

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
                        Matrix const * const lfactorial_vector) // Log-factorial lookup table
{
	// Do a Product-Multinomial MH move on an exit-poll dataset:
	// (1) Fix special_row, and propose multinomials with parameters
	// given by the current state of theta for precinct ii.

#ifdef _DBG_
  int dbg=2;
#endif

	index_t jj, kk, tmp_mult;
	const index_t nr_mp = numrows(multinomial_parameters);
	const index_t nr_pr = numrows(prop_row);
	const index_t nr_NNs = numrows(NNs);
	const index_t nr_MMs = numrows(MMs);
	const index_t nr_KKs = numrows(KKs);
	const index_t nr_THETAS = numrows(THETAS);
	const index_t nr_NNtots = numrows(NNtots);
	const index_t nr_MMtots = numrows(MMtots);
	const index_t nr_KKtots = numrows(KKtots);
	double log_mh_ratio = 0.0;
	double NN_prec_row_tot, NN_prec_col_tot;
	double MM_prec_row_tot, MM_prec_col_tot;
  double tmp_NNs, tmp_MMs;

	// Increment the counter:
	matrix_fast_increment_element(NNs_count_use_multinom,current_precinct,special_row,nr_NNs, 1);

  // Initialize tmp_KKs, NNs_curr and MMs_curr:
  for (jj=0; jj<numrows_pt; jj++){
    for (kk=0; kk<numcols_pt; kk++){

      // Usual time-saver:
      tmp_mult = (jj*numcols_pt);

      // Initialize the KKs matrix for this precinct:
      matrix_fast_set_element(tmp_KKs,jj,kk,numrows_pt, 
        matrix_fast_get_element(KKs,current_precinct,tmp_mult+kk,nr_KKs));

      // Initialize the NNs_curr matrix for this precinct:
      matrix_fast_set_element(NNs_curr,jj,kk,numrows_pt, 
        matrix_fast_get_element(NNs,current_precinct,tmp_mult+kk,nr_NNs));

      // Initialize the MMs_curr matrix for this precinct:
      matrix_fast_set_element(MMs_curr,jj,kk,numrows_pt, 
        matrix_fast_get_element(MMs,current_precinct,tmp_mult+kk,nr_MMs));
    }
  } // End initialization

#ifdef _DBG_
  if (dbg>1){
    Rprintf("\n##### MULTINOMIAL PROPOSAL DIAGNOSTICS -- PRECINCT %u ########\n",current_precinct);
    Rprintf("KKs for precinct %u:\n",current_precinct);
    matrix_print_all(tmp_KKs);
    Rprintf("KKtots rowsums:\n");
    for (jj=0; jj<numrows_pt; jj++)
      Rprintf("%g\t",matrix_fast_get_element(KKtots,current_precinct,jj,nr_KKtots));
    Rprintf("\nKKtots colsums:\n");
    for (jj=0; jj<numcols_pt; jj++)
      Rprintf("%g\t",matrix_fast_get_element(KKtots,current_precinct,numrows_pt+jj,nr_KKtots));
    Rprintf("\nNNs_curr for precinct %u:\n",current_precinct);
    matrix_print_all(NNs_curr);
    Rprintf("NNtots rowsums:\n");
    for (jj=0; jj<numrows_pt; jj++)
      Rprintf("%g\t",matrix_fast_get_element(NNtots,current_precinct,jj,nr_NNtots));
    Rprintf("\nNNtots colsums:\n");
    for (jj=0; jj<numcols_pt; jj++)
      Rprintf("%g\t",matrix_fast_get_element(NNtots,current_precinct,numrows_pt+jj,nr_NNtots));
    Rprintf("\nMMs_curr for precinct %u:\n",current_precinct);
    matrix_print_all(MMs_curr);
    Rprintf("MMtots rowsums:\n");
    for (jj=0; jj<numrows_pt; jj++)
      Rprintf("%g\t",matrix_fast_get_element(MMtots,current_precinct,jj,nr_MMtots));
    Rprintf("\nMMtots colsums:\n");
    for (jj=0; jj<numcols_pt; jj++)
      Rprintf("%g\t",matrix_fast_get_element(MMtots,current_precinct,numrows_pt+jj,nr_MMtots));
    Rprintf("\n");
  }
#endif

	// Find the parameters for the multinomials for the remaining rows:
	for (jj=0; jj<numrows_pt; jj++){

		// Skip the special row:
    if (jj==special_row)
			continue;

		// Save a bit of computation:
		tmp_mult = (jj*numcols_pt);

		// Set the proposal parameters for that row:
		for (kk=0; kk<numcols_pt; kk++){
			matrix_fast_set_element(multinomial_parameters,0,kk,nr_mp, 
											matrix_fast_get_element(THETAS,current_precinct,tmp_mult+kk,nr_THETAS));
    }

		// Retrieve the row totals for row jj:
    MM_prec_row_tot = matrix_fast_get_element(MMtots,current_precinct,jj,nr_MMtots);
    NN_prec_row_tot = matrix_fast_get_element(NNtots,current_precinct,jj,nr_NNtots);

		// Now, make the multinomial draw on the MMs:
		rmultinomial(prop_row,multinomial_parameters,MM_prec_row_tot);

#ifdef _DBG_
  if (dbg>1){
    Rprintf("Proposed row (should sum to row %u total in MMtots):\n",jj);
    matrix_print_all(prop_row);
  }
#endif 

		// Put the draws into matrix form:
		for (kk=0; kk<numcols_pt; kk++){
      tmp_MMs = matrix_fast_get_element(prop_row,0,kk,nr_pr);
      tmp_NNs = tmp_MMs + matrix_fast_get_element(tmp_KKs,jj,kk,numrows_pt);
			matrix_fast_set_element(MMs_prop,jj,kk,numrows_pt, tmp_MMs);
			matrix_fast_set_element(NNs_prop,jj,kk,numrows_pt, tmp_NNs);
    }

	} // End loop over rows

	// Determine if the proposal is valid:
	for (kk=0; kk<numcols_pt; kk++){

    // Column kk total for MMs and NNs:
		MM_prec_col_tot = matrix_fast_get_element(MMtots,current_precinct,numrows_pt+kk,nr_MMtots);
		NN_prec_col_tot = matrix_fast_get_element(NNtots,current_precinct,numrows_pt+kk,nr_NNtots);

#ifdef _DBG_
  if (dbg>1){
    Rprintf("Starting out with %g available MM counts in col %u\n",
            MM_prec_col_tot,kk);
    Rprintf("Starting out with %g available NN counts in col %u\n",
            NN_prec_col_tot,kk);
  }
#endif

		// See how many counts are left in each column:
		for (jj=0; jj<numrows_pt; jj++){
			if (jj!=special_row){
				MM_prec_col_tot -= matrix_fast_get_element(MMs_prop,jj,kk,numrows_pt);
				NN_prec_col_tot -= matrix_fast_get_element(NNs_prop,jj,kk,numrows_pt);

#ifdef _DBG_
  if (dbg>1){
    Rprintf("Row %u used up %g available MM counts, leaving %g of them\n",
            jj,matrix_fast_get_element(MMs_prop,jj,kk,numrows_pt),MM_prec_col_tot);
    Rprintf("Row %u used up %g available NN counts, leaving %g of them\n",
            jj,matrix_fast_get_element(NNs_prop,jj,kk,numrows_pt),NN_prec_col_tot);
  }
#endif

      }
    }

    // Was it valid for this column?
    if (MM_prec_col_tot<0.0 || NN_prec_col_tot<0.0){
		 // No. :(
#ifdef _DBG_
		if (dbg>1){
			Rprintf("Exit-poll multinomial proposal invalid for precinct %u, column %u:\n",current_precinct,kk);
      Rprintf("MM_prec_col_tot = %g\n",MM_prec_col_tot);
      Rprintf("NN_prec_col_tot = %g\n",NN_prec_col_tot);
    }
#endif
      return;
    }
    // Yes! :)

		tmp_mult = (special_row*numcols_pt);

		// Set up the proposed values:
		matrix_fast_set_element(prop_row,0,kk,1, MM_prec_col_tot);
		matrix_fast_set_element(curr_row,0,kk,1, 
										matrix_fast_get_element(MMs,current_precinct,tmp_mult+kk,nr_MMs));

		matrix_fast_set_element(multinomial_parameters,0,kk,1,
										matrix_fast_get_element(THETAS,current_precinct,tmp_mult+kk,nr_THETAS));

		matrix_fast_set_element(MMs_prop,special_row,kk,numrows_pt, MM_prec_col_tot);
		matrix_fast_set_element(NNs_prop,special_row,kk,numrows_pt, NN_prec_col_tot);
	}

	// We have a valid proposal for the whole table. :)

	// Increment the number of valid proposals for this precinct/special row:
  vld_mv_p[current_precinct+(nr_NNtots*special_row)]++;

	// Compute the log(MH-ratio):
	log_mh_ratio = log_MMs_multinomial_mh_ratio(NNs_prop,MMs_prop,
                                              NNs_curr,MMs_curr,
                                              tmp_KKs,THETAS,
                                              current_precinct,
                                              special_row,
                                              numrows_pt,numcols_pt,
                                              multinomial_parameters,
                                              lfactorial_vector);

	if (log(runif(0.0,1.0))<log_mh_ratio){
		// ACCEPTED:
#ifdef _DBG_
		if (dbg>1){
			Rprintf("Multinomial proposal accepted for precinct %u\n",current_precinct);
		}
#endif

		// Increment the number of accepted proposals for this precinct/special row:
  	acc_mv_p[current_precinct+(nr_NNtots*special_row)]++;

		// Update the current state accordingly.
    // Need to update both NNs and MMs (KKs don't change).
    // We already have MMs_prop and NNs_prop in form.

    // Go cell-by-cell:
		for (jj=0; jj<numrows_pt; jj++){

      // Usual time-saver:
			tmp_mult = (jj*numcols_pt);

			for (kk=0; kk<numcols_pt; kk++){
        // First the NNs:
				matrix_fast_set_element(NNs,current_precinct,tmp_mult+kk,nr_NNs,
					matrix_fast_get_element(NNs_prop,jj,kk,numrows_pt));
        // Next, the MMs:
				matrix_fast_set_element(MMs,current_precinct,tmp_mult+kk,nr_MMs,
					matrix_fast_get_element(MMs_prop,jj,kk,numrows_pt));

      } // End column loop
		} // End row loop

	} else {
		// REJECTED:
#ifdef _DBG_
		if (dbg>1){
			Rprintf("\nMultinomial proposal rejected\n");
		}
#endif
		
	}
	return;
}

void rmultinomial(Matrix * const draw,
									Matrix * const p_vector,
									const double total_count_d)
{
  // Draw from the multinomial distribution with specified parameters.
	
	// Error check:
	const index_t k = numcols(draw);

#ifdef _MULTINOMIAL_CHK_
	const index_t nc_pars = numcols(p_vector);
	if (nc_pars!=k)
		error("Length of parameter is different to the draw holder in 'rmultinomial'");
#endif

	// Note: rbinom error checks for negative p_i's.
	int total_count = floor(total_count_d+0.5);

#ifdef _MULTINOMIAL_CHK_
	if (total_count!=total_count_d)
		error("Non-integer count supplied to rmultinomial'");
#endif

	// Short-circuit:
	if (total_count==0)
		matrix_fastset(draw,0);

	// Draw a sequence of binomials:
	double n_left = total_count_d; 	// NOTE: *binom works with doubles NOT ints.
	double tmp_draw, p_frac;
	const index_t k_m1 = (k-1);
	const index_t nr_pv = numrows(p_vector);
	const index_t nr_draw = numrows(draw);

	// p_vector should sum to 1 but we renormalize just in case:
  double p_left = 0.0;
	index_t ii;
	for (ii=0; ii<k; ii++)
			p_left += matrix_fast_get_element(p_vector,0,ii,nr_pv); 
  
	for (ii=0; ii<k_m1; ii++){
		p_frac = matrix_fast_get_element(p_vector,0,ii,nr_pv)/p_left;
		// Draw from the conditional binomial:
		tmp_draw = rbinom(n_left,p_frac);
		// Set the draw in place:
		matrix_fast_set_element(draw,0,ii,nr_draw, tmp_draw);
		// Remove p_i from the renormalization:
		p_left -= matrix_fast_get_element(p_vector,0,ii,nr_pv);
		// Remove n_i from the remaining counts:
		n_left -= tmp_draw;
	}

	// Fill in the final cell:
	matrix_fast_set_element(draw,0,k_m1,nr_draw, n_left);

	return;
}


