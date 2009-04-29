#ifndef _NCHG_C
#define _NCHG_C

#include "gqjrssa.h"

#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

// Faster when not inlined:
double logsumexp_V22(double *f, int n);

// Inlined basic functions:
GQ_INLINE int gq_sign(double x){return (x<0)?-1:1;}
GQ_INLINE double gq_min(double x, double y){return (x<y)?x:y;}
GQ_INLINE double gq_max(double x, double y){return (x>y)?x:y;}

// Recursive functions:
GQ_INLINE double ri_V22(double m1,double n1,double n2,double psi,double i) 
{
	return psi*((n1-i+1.0)*(m1-i+1.0))/(i*(n2-m1+i));
}
GQ_INLINE double ru_V22(double m1,double n1,double n2,double psi,double i) 
{
	return psi*((n1-i+1.0)*(m1-i+1.0))/(i*(n2-m1+i));
}
GQ_INLINE double rl_V22(double m1,double n1,double n2,double psi,double i) 
{
	return psi*((n1-i)*(m1-i))/((i+1.0)*(n2-m1+i+1.0));
}

GQ_INLINE double pdb_ri(const double n1p1_mult_m1p1, const double n1p1_plus_m1p1, 
			const double n2mm1, const double psi, const double i)
{
	return (psi*(n1p1_mult_m1p1 - i*n1p1_plus_m1p1 + i*i)/(i*(n2mm1+i)));
}

double rnoncenhypgeo_omit(double m1, double n1, double n2, double psi, double delta, double *ff_vec, const int debug);

double jims_original_rnchg(double m1, double n1, double n2, double psi, double delta, double *ff_vec, const int debug);
double jims_from_byron_v20_rnchg(double m1, double n1, double n2, double psi, double delta, double *ff_vec, const int debug);
double byron_used_from_V20_rnchg(double n1,double n2,double m1,double psi,const int debug);
double byron_from_V22_rnchg(double m1,double n1,double n2,double psi,double delta,double *fvec,const int debug);
double rnoncenhypgeo(double m1,double n1,double n2,double psi,double delta,double *fvec,const int debug);
double old_rnoncenhypgeo(double m1,double n1,double n2,double psi,double delta, double *fvec, const int debug);

double rnoncenhypgeo(const double m1, const double n1, const double n2,
 		     const double psi, const double delta, double * const fvec,
		     const int debug);

#endif
