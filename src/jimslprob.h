#ifndef _JIMSLPROB_H_
#define _JIMSLPROB_H_

#include "jimsmatrix.h"

double 
complete_state_lprob(Matrix * THETAS,
		     Matrix * OMEGAS,
		     Matrix * mu_vec_cu,
		     Matrix * SIGMA_cu,
		     Matrix * SIGMA_chol_cu,
		     Matrix * NNs,
		     Matrix * NNtots,
		     Matrix * NNbounds);

#endif
