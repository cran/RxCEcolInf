#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .Call calls */
extern SEXP Analyze(SEXP);
extern SEXP AnalyzeWithExitPoll(SEXP);
extern SEXP rnchg(SEXP);
extern SEXP Tune(SEXP);
extern SEXP TuneWithExitPoll(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"C_Analyze",             (DL_FUNC) &Analyze,             1},
  {"C_AnalyzeWithExitPoll", (DL_FUNC) &AnalyzeWithExitPoll, 1},
  {"C_rnchg",               (DL_FUNC) &rnchg,               1},
  {"C_Tune",                (DL_FUNC) &Tune,                1},
  {"C_TuneWithExitPoll",    (DL_FUNC) &TuneWithExitPoll,    1},
  {NULL, NULL, 0}
};

void R_init_RxCEcolInf(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
