#ifndef _GQJRSSA_H
#define _GQJRSSA_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

#include "jimsmatrix.h"
#include "jimsimatrix.h"
#include "jimsmatrixinterface.h"
#include "jimsmatrix_p.h"
#include "jimsmisc.h"
#include "jimsMH.h"
#include "jimsrandom.h"
#include "jimslprob.h"

SEXP Tune(SEXP args);
SEXP Analyze(SEXP args);
SEXP rmultinomial_chk(SEXP pvec, SEXP n);

#endif
