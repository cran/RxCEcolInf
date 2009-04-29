#ifndef _JIMSMISC_H_
#define _JIMSMISC_H_

#include <R.h>
#include <Rinternals.h>
#include <assert.h>
#include "jimsmatrix.h"
#include "jimsimatrix.h"
#include "jimsmatrix_p.h"

SEXP 
getListElement(SEXP list, char *str);

Matrix *
create_log_factorial_lookup_table(const index_t max);

void 
store_internals(Matrix * const ToBeStored, 
                Matrix * const Store_mat, 
                long * const colnum_store);

void
THETAS_to_OMEGAS(Matrix * const THETAS, 
                 Matrix * const OMEGAS, 
		             const index_t numrows_pt, 
                 const index_t numcols_pt);

index_t GQ_INLINE 
runif_index_t(const index_t lo, 
              const index_t hi)
{
  return (index_t)runif((double)lo,(double)(hi+1));
}

int GQ_INLINE
check_bounds(const double prop_val, 
             Matrix * const NNbounds,
	           const index_t num_p, 
             const index_t position, 
             const index_t numcells_pt)
{ // Returns 1 if prop_val is within the cell's bounds, 0 if outside.
  return (((prop_val < matrix_get_element(NNbounds, num_p, position)) ||
           (prop_val > matrix_get_element(NNbounds, num_p, position + numcells_pt))) ?
	0 : 1);
}

void
check_ep_validity(Matrix * const NNs,
                  Matrix * const MMs,
                  Matrix * const KKs,
                  const index_t numcells_pt,
                  const index_t numrows_pt,
                  const index_t numcols_pt);

void 
adjust_rho_vec(Matrix * const rho_vec, 
               SEXP acc_THETAS_t_vec);

void 
adjust_acc_vector(SEXP acc_vec, 
                  Matrix * const count_use_vec);

void 
adjust_acc_matrix(SEXP acc_mat, 
                  Matrix * const count_use_mat);

void
report_dry_run(SEXP acc_THETAS_t_vec, 
               SEXP acc_THETAS_Diri_vec,
               const long num_run);

void
check_bounds_all(Matrix * const NNs, 
                 Matrix * const NNbounds, 
                 const index_t numcells_pt);

void
store_draws(Matrix * const NNs, 
            Matrix * const mu_vec_cu, 
            Matrix * const SIGMA_cu,
            Matrix * const draws, 
            long * const rownum_store);

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
    const int ncol_BETA);

Matrix_int *
get_which_rc(const index_t nr_pt, const index_t nc_pt);

double
Rmatrix_get_fraction_under_c(SEXP xx, double c);

double
Rmatrix_get_fraction_over_c(SEXP xx, double c);

double
Rmatrix_max(SEXP xx);

double
Rmatrix_min(SEXP xx);

void
initialize_KKtots_and_MMtots(Matrix * const KKtots,
                             Matrix * const MMtots,
                             Matrix * const NNtots,
                             Matrix * const KKs,
                             const index_t numrows_pt,
                             const index_t numcols_pt);

void
multiply_list_of_X_by_eta(Matrix * const mu_mat_cu,
                          Matrix_p * const Xmatrix_list,
                          Matrix * const eta_vec_cu);

#endif /* jimsmisc.h */
