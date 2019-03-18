#include "nchg.h"

SEXP rnchg_test(SEXP args)
{
  index_t ii;
  const unsigned int debug = INTEGER(getListElement(args,"debug"))[0];
  const unsigned int n = INTEGER(getListElement(args,"n"))[0];  
  const char funcToUse = INTEGER(getListElement(args,"method"))[0];

  Matrix * const m1 = matrix_vector_unpack_new(getListElement(args, "m1"));
  Matrix * const n1 = matrix_vector_unpack_new(getListElement(args, "n1"));
  Matrix * const n2 = matrix_vector_unpack_new(getListElement(args, "n2"));
  Matrix * const psi= matrix_vector_unpack_new(getListElement(args, "psi"));

  // Recycle each vector as required:
  const index_t len_m1 = numcols(m1);
  const index_t len_n1 = numcols(n1);
  const index_t len_n2 = numcols(n2);
  const index_t len_psi= numcols(psi);

  static const long maxrange = 250000001;

  double * const ff_vec = (double *)R_alloc(maxrange,sizeof(double));
  if (ff_vec==NULL)
    error("Memory allocation failure (ff_vec)");
  
  if (debug){
    Rprintf("n1:\n");
    matrix_print_all(n1);
    Rprintf("n2:\n");
    matrix_print_all(n2);
    Rprintf("m1:\n");
    matrix_print_all(m1);
    Rprintf("psi:\n");
    matrix_print_all(psi);
  }
  SEXP drawsForR;
  PROTECT(drawsForR=allocVector(REALSXP,n));
  double * const rDraws = REAL(drawsForR);

  GetRNGstate();

/*
 * CHOOSE FROM THESE FUNCTIONS:
 * double       rnoncenhypgeo(double m1, double n1, double n2, double psi, double delta, double *ff_vec);
 * double jims_original_rnchg(double m1, double n1, double n2, double psi, double delta, double *ff_vec);
 * double jims_from_byron_v20_rnchg(double m1, double n1, double n2, double psi, double delta, double *ff_vec);
 * double byron_used_from_V20_rnchg(double n1, double n2, double m1, double psi);
 * double      byron_from_V22_rnchg(double m1, double n1, double n2, double psi, double delta, ,double *fvec);
 */

  // Use the version was 'c'urrent version:
  if (funcToUse == 0){
    for (ii=0; ii<n; ii++)
      rDraws[ii] = rnoncenhypgeo(matrix_fast_get_element(m1,0,ii%len_m1,1),
			 	 matrix_fast_get_element(n1,0,ii%len_n1,1),
				 matrix_fast_get_element(n2,0,ii%len_n2,1),
				 matrix_fast_get_element(psi,0,ii%len_psi,1),
				 DBL_EPSILON,ff_vec,debug);
  }

  // Use Jim's 'o'riginal function:
  if (funcToUse == 1){
    for (ii=0; ii<n; ii++)
      rDraws[ii] = jims_original_rnchg(matrix_fast_get_element(m1,0,ii%len_m1,1),
				       matrix_fast_get_element(n1,0,ii%len_n1,1),
			               matrix_fast_get_element(n2,0,ii%len_n2,1),
				       matrix_fast_get_element(psi,0,ii%len_psi,1),
				       DBL_EPSILON,ff_vec,debug);
  }

  // The 'f'inal attempt by Byron (in GQ_22):
  if (funcToUse == 2){
    for (ii=0; ii<n; ii++)
      rDraws[ii] = byron_from_V22_rnchg(matrix_fast_get_element(m1,0,ii%len_m1,1),
					matrix_fast_get_element(n1,0,ii%len_n1,1),
					matrix_fast_get_element(n2,0,ii%len_n2,1),
					matrix_fast_get_element(psi,0,ii%len_psi,1),
					DBL_EPSILON,ff_vec,debug);
  }

  // The 'i'nterim version used by Byron in V20:
  if (funcToUse == 3){
    for (ii=0; ii<n; ii++)
      rDraws[ii] = byron_used_from_V20_rnchg(matrix_fast_get_element(n1,0,ii%len_n1,1),
					     matrix_fast_get_element(n2,0,ii%len_n1,1),
				             matrix_fast_get_element(m1,0,ii%len_m1,1),
					     matrix_fast_get_element(psi,0,ii%len_psi,1),
					     debug);
  }

  matrix_free(m1);
  matrix_free(n1);
  matrix_free(n2);
  matrix_free(psi);

  PutRNGstate();

  UNPROTECT(1);
  return drawsForR;
}

SEXP rnchg(SEXP args)
{
  index_t ii;
  const unsigned int n = INTEGER(getListElement(args,"n"))[0];  
  Matrix * const m1 = matrix_vector_unpack_new(getListElement(args, "m1"));
  Matrix * const n1 = matrix_vector_unpack_new(getListElement(args, "n1"));
  Matrix * const n2 = matrix_vector_unpack_new(getListElement(args, "n2"));
  Matrix * const psi= matrix_vector_unpack_new(getListElement(args, "psi"));

  // Recycle each vector as required:
  const index_t len_m1 = numcols(m1);
  const index_t len_n1 = numcols(n1);
  const index_t len_n2 = numcols(n2);
  const index_t len_psi= numcols(psi);

  // Find required length of ff_vec:
  double tmp_val, max_range=0.0;
  for (ii=0; ii<n; ii++){
    tmp_val = matrix_fast_get_element(n1,0,ii%len_n1,1)+
	      matrix_fast_get_element(n2,0,ii%len_n2,1);
    if (tmp_val>max_range)
      max_range = tmp_val;
  }
  const long lmaxrange = (long)(max_range+1);
  double * const ff_vec = (double *)R_alloc(lmaxrange,sizeof(double));
  if (ff_vec==NULL)
    error("Memory allocation failure (ff_vec)");
  
  SEXP drawsForR;
  PROTECT(drawsForR=allocVector(REALSXP,n));
  double * const rDraws = REAL(drawsForR);

  GetRNGstate();

  for (ii=0; ii<n; ii++)
    rDraws[ii] = rnoncenhypgeo(matrix_fast_get_element(m1,0,ii%len_m1,1),
		 	       matrix_fast_get_element(n1,0,ii%len_n1,1),
			       matrix_fast_get_element(n2,0,ii%len_n2,1),
			       matrix_fast_get_element(psi,0,ii%len_psi,1),
			       DBL_EPSILON,ff_vec,0);
  matrix_free(m1);
  matrix_free(n1);
  matrix_free(n2);
  matrix_free(psi);

  PutRNGstate();

  UNPROTECT(1);
  return drawsForR;
}

double logsumexp_V22(double *f, int n)
{
	double max = f[0];
	int    i;
	for(i=1;i<n;i++) if(f[i] > max) max = f[i];
	double s = 0.0;
	for(i=0;i<n;i++) s += exp(f[i]-max);
	return max+log(s);
}

double rnoncenhypgeo_omit(double m1, double n1, double n2, double psi, double delta, double *ff_vec, const int debug)
{
        double epsilon = delta/10.0;
        double a = psi - 1.0;
        double b = -((n1+m1+2.0)*psi+n2-m1);
        double c = psi*(n1+1.0)*(m1+1.0);
        double q = -.5*(b+gq_sign(b)*sqrt(b*b-4*a*c));
        double r1 = c/q;
        double r2 = q/a;
        //double t  = m1-n2;
        int el,uu;

        if (m1>n2)        
	  el = m1-n2; 
	else 
	  el = 0;

        if (n1>m1) 
	  uu = n1; 
	else 
	  uu = m1;

        int mode = floor(r1);
        //int j, k;
	int i, top, bottom, exact = 0;

        if ( (uu<mode) || (mode<el) )
	{
                mode  = floor(r2);
                exact = 1;
        }

        ff_vec[(mode-el)] = 0.0;
        double f = 0;
        double s = 0;

	// My changes (just pre-computes things that dont change within loops):
	double n1p1    = n1+1.0;
	double m1p1    = m1+1.0;
	double n2mm1   = n2-m1;
	double n2mm1p1 = n2mm1+1.0;
	double lpsi    = log(psi);
	const double fivesixths = 5.0/6.0;
	const double sixfifths = 6.0/5.0;
	double lepsilon = log(epsilon);

        if ( (delta<=0.0) || exact ){
                for (i=(mode+1); i<=uu; i++){
                        double r = lpsi+log(n1p1-(double)i)+log(m1p1-(double)i)-log((n2mm1+(double)i))-log((double)i);
                        f += r;
                        if (isnan(f)) 
			  f = R_NegInf;
                        s  = logspace_add(s,f);
                        ff_vec[i-el] = f;
                }
                f = 0;
                for (i=(mode-1); i>=el; i--) {
                        double r = lpsi+log(n1-(double)i)+log(m1-(double)i)-log(1.0+(double)i)-log(n2mm1p1+(double)i);
                        f -= r;
                        if(isnan(f)) 
			  f = R_NegInf;                        
                        s  = logspace_add(s,f);
                        ff_vec[i-el] = f;
                }
                top    = uu;
                bottom = el;
        } else {
                i = mode+1;
                for (i=(mode+1); i<=uu; i++){
                        double r = lpsi+log(n1p1-(double)i)+log(m1p1-(double)i)-log((n2mm1+(double)i))-log((double)i);
                        f += r;
                        if (isnan(f)) 
			  f = R_NegInf;                        
                        s  = logspace_add(s,f);
                        ff_vec[i-el] = f;
                        if ( (f<=lepsilon) || (r>=fivesixths) )
			  break;
                }
                f   = 0;
                top = i;
                for (i=(mode-1); i>=el; i--){
                        double r = lpsi+log(n1-(double)i)+log(m1-(double)i)-log(1.0+(double)i)-log(n2mm1p1+(double)i);
                        f -= r;
                        if (isnan(f)) 
			  f = R_NegInf;                        
                        s  = logspace_add(s,f);
                        ff_vec[i-el] = f;                        
                        if ( (f<=lepsilon) || (r<=sixfifths) );// SHOULD THIS LINE BE LIKE THIS???
                }
                bottom = i;
        }
        if (bottom < el) 
	  bottom = el;
        if (top > uu) 
	  top = uu;
        double udraw = unif_rand();
        double cdf   = 0.0;
        for (i=bottom; i<=top; i++){
                cdf += exp(ff_vec[i-el]-s);
                if (udraw <= cdf) 
		  break;
        }
//        Rprintf("%d %f %f %d [%d:%d %d %d:%d]\n",exact,udraw,cdf,i,el,bottom,mode,top,uu);
//        if(isnan(cdf)) {
//                for(i=bottom;i<top;i++) Rprintf("%.2f\t",ff_vec[i-el]);
//                Rprintf("\n");
//                error("Dead");
//        }
        return i;
}




//  Function below was my original noncentral hypergeometric function

double
jims_original_rnchg(double m1, double n1, double n2, double psi, double delta, double *ff_vec, const int debug)
{
//      Returns a random draw from the noncentral hypergeometric distribution.
//      See J.G. LiaO & Ori Rosen, Fast And Stable Algorithms for Computing and
//      Sampling from the Noncentral Hypergeometric Distribution, 55 The Am.
//      Stat. 366 (2001), for the parameterization
//      and the approximation algorithm used.  Thanks to Kevin Quinn for
//      his C++ version, from which this function copies liberally.

  double aa, bb, cc, qq, root1, root2, el, uu, mode, ss = 1.0;
  double epsilon = delta/10.0;
  double ff = 1.0, ii, rr, udraw, psum, lower, upper, fl, fu;
  int exactcheck = 0, upcounter = 0, downcounter = 0;

  //  Calculate the mode
  aa = psi -1;
  bb = -1.0 * ((n1 + m1 + 2)*psi + n2 - m1);
  cc = psi * (n1+1) * (m1+1);
  qq = -0.5*(bb + ((gq_sign(bb))*sqrt(bb*bb - 4*aa*cc)));
  root1 = cc/qq;
  root2 = qq/aa;
  el = max(0.0, m1-n2);
  uu = min(n1, m1);
  mode = floor(root1);
  if(uu < mode || mode < el){
    mode = floor(root2);
    exactcheck = 1;
  }  // Possible inequality mistake -- mode < el+1 instead?

  // ff_vec to be arranged as follows:  ff_vec[0] corresponds to el, ff_vec[1] to
  //        el + 1, . . . ff_vec[mode-el] corresponds to mode, . . .  
  ff_vec[(int)(mode-el)] = 1.0;

  //  Seems like the following block is unnecessary
  //  compute the mass function at y
  if(delta <= 0.0 || exactcheck == 1){
    ff = 1.0;
    ss = 1.0;
    for(ii = (mode + 1.0); ii <= uu; ++ii){
      rr = psi * ((n1-ii+1.0)*(m1-ii+1.0))/(ii*(n2-m1+ii));
      ff = ff*rr;
      ss += ff;
      ff_vec[(int)(ii-el)] = ff;
    }
    //  sum from mode to el
    ff = 1.0;
    for(ii = (mode - 1.0); ii>=el; --ii){
      rr = psi*((n1-ii)*(m1-ii))/((ii+1.0)*(n2-m1+ii+1.0));
      ff = ff/rr;
      ss += ff;
      ff_vec[(int)(ii-el)] = ff;
    }
  }else{ //approximation
    //Calculate from mode+1 to upper approximation limit.
    ff = 1.0;
    ss = 1.0;
    ii = mode + 1.0;
    do{
      if(ii > uu) break;
      rr = psi * ((n1-ii+1.0)*(m1-ii+1.0))/(ii*(n2-m1+ii));
      ff = ff*rr;
      ss += ff;
      ff_vec[(int)(ii - el)] = ff;
      ++ii;
      ++upcounter;
    } while(ff >= epsilon || rr >=(5.0/6.0));

    //  Calculate from mode -1 to lower approximation limit.
    ff = 1.0;
    ii = mode - 1.0;
    do{
      if(ii < el) break;
      rr = psi * ((n1-ii)*(m1-ii))/((ii+1.0)*(n2-m1+ii+1.0));
      ff = ff/rr;
      ss += ff;
      ff_vec[(int)(ii - el)] = ff;
      --ii;
      ++downcounter;
    }while(ff >= epsilon || rr <= (6.0/5.0));
  }

  //GetRNGstate();
  udraw = runif(0.0, 1.0);
  //PutRNGstate();

  psum = ff_vec[(int)(mode - el)]/ss;
  if(udraw <= psum){
#ifdef _DBG_
    if (debug)
      Rprintf("Returning the mode %g (limits: %f,%f) in jims_rnchg...\n",mode,el,uu);
#endif
    return mode;
  }
  lower = mode - 1.0;
  upper = mode + 1.0;

  do{
    if(lower >= el) fl = ff_vec[(int)(lower - el)];
    else fl = 0.0;
    if(upper <= uu) fu = ff_vec[(int)(upper - el)];
    else fu = 0.0;

    if(fl > fu){
      psum += fl/ss;
      if(udraw <= psum){
#ifdef _DBG_
        if (debug)
          Rprintf("Returning %g (limits: %f,%f) in jims_rnchg...\n",lower,el,uu);
#endif
        return lower;
      }
      --lower;
    }
    else{
      psum += fu/ss;
      if(udraw <= psum){
#ifdef _DBG_
        if (debug)
          Rprintf("Returning %g (limits: %f,%f) in jims_rnchg...\n",upper,el,uu);
#endif
	 return upper;
      }
      ++upper;
    }
  } while(udraw > psum);

  //exit(500000);
  error("ERROR: problem with jims_original_rnchg()");
}

double
jims_from_byron_v20_rnchg(double m1, double n1, double n2, double psi, double delta, double *ff_vec, const int debug)
{
  //  Returns a random draw from the noncentral hypergeometric distribution.
  //    See J.G. LiaO & Ori Rosen, Fast And Stable Algorithms for Computing and
  //    Sampling from the Noncentral Hypergeometric Distribution, 55 The Am.
  //    Stat. 366 (2001), for the parameterization
  //    and the approximation algorithm used.  Thanks to Kevin Quinn for
  //    his C++ version, from which this function copies liberally.

  	double aa, bb, cc, qq, root1, root2, el, uu, ss = 1.0;
  	int mode;
  	double epsilon = delta/10.0, ff = 1.0, ii, rr, udraw, psum, lower, upper, fl, fu;
  	int exactcheck = 0, upcounter = 0, downcounter = 0;

  	//  Calculate the mode
  aa = psi -1;
  bb = -1.0 * ((n1 + m1 + 2)*psi + n2 - m1);
  cc = psi * (n1+1) * (m1+1);
  qq = -0.5*(bb + ((gq_sign(bb))*sqrt(bb*bb - 4*aa*cc)));
  root1 = cc/qq;
  root2 = qq/aa;
  el = max(0.0, m1-n2);
  uu = min(n1, m1);
  mode = floor(root1);
  if(uu < mode || mode < el){
    mode = floor(root2);
    exactcheck = 1;
  }  // verify this later?

  // ff_vec to be arranged as follows:  ff_vec[0] corresponds to el, ff_vec[1] to
  //        el + 1, . . . ff_vec[mode-el] corresponds to mode, . . .  
  ff_vec[(int)(mode-el)] = 1.0;

  //  Seems like the following block is unnecessary
  //  compute the mass function at y
  if(delta <= 0.0 || exactcheck == 1){
    ff = 1.0;
    ss = 1.0;
    for(ii = (mode + 1.0); ii <= uu; ++ii){
      rr = psi * (((n1-ii+1.0)*(m1-ii+1.0))/(ii*(n2-m1+ii)));
      ff = ff*rr;
      ss += ff;
      ff_vec[(int)(ii-el)] = ff;
    }
    //  sum from mode to el
    ff = 1.0;
    for(ii = (mode - 1.0); ii>=el; --ii){
      rr = psi*((n1-ii)*(m1-ii))/((ii+1.0)*(n2-m1+ii+1.0));
      ff = ff/rr;
      ss += ff;
      ff_vec[(int)(ii-el)] = ff;
    }
//	if(!R_FINITE(ss))
//		error("SS is not finite %.2f %.2f %.2f %.2f %.2f",m1,n1,n2,psi,delta);

  }else{ //approximation
    //Calculate from mode+1 to upper approximation limit.
    ff = 1.0;
    ss = 1.0;
    ii = mode + 1.0;
    do{
      if(ii > uu) break;
      rr = psi * ((n1-ii+1.0)*(m1-ii+1.0))/(ii*(n2-m1+ii));
      ff = ff*rr;
      ss += ff;
      ff_vec[(int)(ii - el)] = ff;
      ++ii;
      ++upcounter;
    } while(ff >= epsilon || rr >=(5.0/6.0));

    //  Calculate from mode -1 to lower approximation limit.
    ff = 1.0;
    ii = mode - 1.0;
    do{
      if(ii < el) break;
      rr = psi * ((n1-ii)*(m1-ii))/((ii+1.0)*(n2-m1+ii+1.0));
      ff = ff/rr;
      ss += ff;
      ff_vec[(int)(ii - el)] = ff;
      --ii;
      ++downcounter;
    }while(ff >= epsilon || rr <= (6.0/5.0));
//    if(!R_FINITE(ss))
//		error("ss is not finite %.2f %.2f %.2f %.2f %.2f",m1,n1,n2,psi,delta);
  }


  //GetRNGstate();
  udraw = unif_rand();
  //PutRNGstate();

  psum = ff_vec[(int)(mode - el)]/ss;
  if(udraw <= psum) return mode;
  lower = mode - 1.0;
  upper = mode + 1.0;

  do{
    if(lower >= el) fl = ff_vec[(int)(lower - el)];
    else fl = 0.0;
    if(upper <= uu) fu = ff_vec[(int)(upper - el)];
    else fu = 0.0;

    if(fl > fu){
      psum += fl/ss;
      if(udraw <= psum) return lower;
      --lower;
    }
    else{
      psum += fu/ss;
      if(udraw <= psum) return upper;
      ++upper;
    }
  } while(udraw > psum);
	error("Outside of CDF loop  %f [%f %f %f]",psi,lower,upper,mode);
  // Should never reach here:
  error("Functional failure");
  return R_NegInf;
}

double byron_used_from_V20_rnchg(double n1,double n2,double m1,double psi,const int debug) 
{
	double l = m1-n2 > 0 ? m1-n2 : 0;
	double u = n1 > m1 ? m1 : n1;
	if(n1 < 0 || n2 < 0 || m1 < 0) error("Invalid parameters %f %f %f",n1,n2,m1);
	if(l==u) return l; //Only one choice...
	
	double a = psi-1;
	double b = -((n1+m1+2)*psi+n2-m1);
	double c = psi*(n1+1)*(m1+1);
	double q = -0.5*(b+gq_sign(b)*sqrt(b*b-4.0*a*c));

	double mode = c/q,mode2 = q/a; 
	int   i,nu  = (int)floor(mode);
	
	//No modes. Bail.
	if(isnan(mode) && isnan(mode2)) { return R_NegInf; }
	
	//Try the other mode if we are outside of our range.	
	if((double)nu < l || (double)nu > u) {
		mode2 = floor(mode2);
		//If the other mode is ALSO outside of the range, pick the closest one 
		//to the range.
		if(mode2 < l || mode2 > u) {
			// the following 1 line were modified by Jerry Yu on 3/6/2019. 
			// WARNIG: nchg.c:555:13: warning: using integer absolute value function 'abs' when argument is of floating point type [-Wabsolute-value] 
			int d1 = abs(((double)nu > u)  ? nu - (int)u : (int) (l - nu));// Jerry Yu changed 1-nu to (int) (1-un) to make the number an integer
			int d2 = abs(mode2 > u ? (int)(mode2 - u) : (int)(l-mode2));
			//Too far away. Bail.
			if(d1 > 1000 && d2 > 1000) { return R_NegInf; }
			if(d2 < d1) nu = (int)mode2;
		} else
			nu = (int)mode2;
	}

	char *vmax = vmaxget();
	double  *f = (double*)R_alloc(sizeof(double),1+u-l);	
//	for(i=0;i<=u-l;i++) f[i] = R_NegInf;
	if((double)nu >= l && (double)nu <= u) {
		//nu is inside our range so we can do both sides
		f[nu-(int)l] = 0.0;
		for(i=nu+1;i<=(int)u;i++)
			f[i-(int)l] = f[(i-1)-(int)l]+log(ri_V22(n1,n2,m1,psi,i));
		for(i=nu-1;i>=(int)l;i--)
			f[i-(int)l] = f[(i+1)-(int)l]-log(ri_V22(n1,n2,m1,psi,i+1));
	} else {
		//nu is outside of our range
		double lf = 0.0;
		if((double)nu > u) {
			//We only need to count down from the mode
			for(i=nu-1;i>=(int)l;i--) {
				lf = lf+log(ri_V22(n1,n2,m1,psi,i+1));
				if(i<=u) f[i-(int)l] = lf;
			}
		} else {
			for(i=nu+1;i<=(int)u;i++) {
				lf = lf-log(ri_V22(n1,n2,m1,psi,i));
				if(i>=l) f[i-(int)l] = lf;
			}
		}
	}
	//Normalize
	double s = R_NegInf;
	for(i=0;i<=(int)(u-l);i++) s = logspace_add(s,f[i]);
	double U = log(unif_rand());
	double cdf = R_NegInf;
	
	for(i=0;i<u-l;i++) {
		cdf = logspace_add(cdf,f[i]-s);
		if(U < cdf) break;
	}
//	Rprintf("%f < %f [%d %f] = %d (%d)\n",U,cdf,i,l,i+(int)l,nu);
	vmaxset(vmax);
	return i+(int)l;	
	
}

//Function below was written by Byron to substitute for mine (which follows it),
//     which malfunctioned on datasets with large precinct populations,
//     probably due to underflow


//This version of rnoncenhypgeo attempts to do several things. First, it attempts to
//avoid underflow by using log-space transformations when the mode exceeds the bounds
//of the distribution. Second, when the mode doesn't, it uses non-logspace but does use
//the 5/6 optimization as well as cut-down-from-mode CDF evaluation. This should hopefully
//avoid both crashing trouble and underflow problems. Also, it allocates its own memory 
//via R_alloc, making ffvec (a source of crashing I suspect) a nonissue. 

//double byron_from_V22_rnchg(double m1,double n1,double n2,double psi,double delta, double *fvec, const int debug) 
//double rnoncenhypgeo(double m1, double n1, double n2, double psi, double delta, double *ff_vec, const int debug)
double byron_from_V22_rnchg(double m1,double n1,double n2,double psi,double delta, double *fvec, const int debug)
{
  return rnoncenhypgeo(m1,n1,n2,psi,delta,fvec,debug);
}
/*
double rnoncenhypgeo(double m1,double n1,double n2,double psi,double delta, double *fvec, const int debug)
{
	double old = old_rnoncenhypgeo(m1,n1,n2,psi,delta,fvec,debug);
	double new = new_rnoncenhypgeo(m1,n1,n2,psi,delta,fvec,debug);
	Rprintf("(old,new)=(%d,%d)\n",(int)old,(int)new);
	return old;
}
*/
//double rnoncenhypgeo(double m1,double n1,double n2,double psi,double delta, double *fvec, const int debug) 
double old_rnoncenhypgeo(double m1,double n1,double n2,double psi,double delta, double *fvec, const int debug) 
{
	//Sanity check.
	if(n1 < 0 || n2 < 0 || m1 < 0) error("Invalid parameters: %f,%f,%f",m1,n1,n2); 

	// PB: I added this check just in case:
	if(m1>(n1+n2)){ 
		Rprintf("Error: Invalid Parameters %f + %f > %f\n",n1,n2,m1);
		return R_NegInf;
	}

	//Calculate bounds
	double l = m1-n2 > 0 ? m1-n2 : 0; 
	double u = n1 > m1 ? m1 : n1;     
	
	//Easy case. :-)
	if(l==u){
		if (debug)
		  Rprintf("(l==u) Returning %f\n",l);
		return l;
	}
	
	//Try to find a mode inside our bounds
	double a = psi - 1.0;
	double b = -((n1+m1+2)*psi+n2-m1);
	double c = psi*(n1+1.0)*(m1+1.0);
	double q = -0.5*(b+gq_sign(b)*sqrt(b*b-4.0*a*c));
	double M1= floor(c/q),M2 = floor(q/a);
	if(isnan(M1) && isnan(M2)) {
		error("Both roots are NA in rnoncenhypgeo(): psi=%g,n1=%g,n2=%g,m1=%g",psi,n1,n2,m1);
		//return R_NegInf; //No modes, just bail out of this round for safety
	}
	//Find the integer mode, hopefully in the range otherwise the
	//mode closest to the actual range.
	int i,nu;

	//if(M1 < l || M1 > u) {
	//	if(M2 <l || M2 > u) {

	if (M1 < (l+1) || M1 > u) {
		if(M2 <(l+1) || M2 > u) {
			//Neither mode is in the range. Pick the closest one.
#ifdef _DBG_
			Rprintf("Neither mode is in range...\n");
#endif
			int d1 = (int)fabs(M1 > u ? M1 - u : l - M1);
			int d2 = (int)fabs(M2 > u ? M2 - u : l - M2);
			nu = (d2 < d1) ? M2 : M1;
#ifdef _DBG_
			Rprintf("Byron's mode: %d\n",nu);
#endif			
		} else
			nu = M2;
	} else
		nu = M1;
		
	//Integer bounds for ease of typing
	int     ui   = (int)u;
	int     li   = (int)l;
	
	//Density calculation vector
//	if(nu >= l && nu <= u) {
	if(nu > l && nu <= u) {
		double eps = delta/10.0;
		double S = 1.0;
		int    ustar,lstar;
#ifdef _DBG_
		if (debug)
		  Rprintf("Our mode is in the range, use the tricks...\n");
#endif
		fvec[nu-li] = 1.0;
		//Above the mode
		for(i=nu+1;i<=ui;i++) {
			double r = ri_V22(m1,n1,n2,psi,i);
			fvec[i-li] = fvec[i-1-li]*r;
			S      += fvec[i-li];
			if(fvec[i-li] <= eps && r < 5.0/6.0) break;
		}
		ustar = i > ui ? ui : i;  //This is the value we should use for our upper bound
#ifdef _DBG_
		if (debug)
	          Rprintf("Upper bound: %d\n",ustar);
#endif
		for(i=nu-1;i>=li;i--) {
			double r = ri_V22(m1,n1,n2,psi,i+1);
			fvec[i-li]  = fvec[i+1-li]/r;
			S          += fvec[i-li];
			if(fvec[i-li] <= eps && r > 6.0/5.0) break;
		}
		lstar = i < li ? li : i; //And the lower bound
#ifdef _DBG_
		if (debug)
	          Rprintf("Lower bound: %d\n",lstar);
#endif
		//Now draw using cut-down-from-mode
		double u = unif_rand();
		double cdf = fvec[nu-li]/S;
		if(u < cdf) {
#ifdef _DBG_
			if (debug)
			  Rprintf("Drawing the mode (%g<%g): Returning %d\n",u,cdf,nu);
#endif
			return nu;
		}
		int    cl= nu-1,ul=nu+1;
		double fl,fu;
		do {
			fl = (cl >= lstar) ? fvec[cl-li]/S : 0.0;
			fu = (ul <= ustar) ? fvec[ul-li]/S : 0.0;
			if(fl > fu) {
				cdf += fl;
				if(u <= cdf) {
#ifdef _DBG_
					if (debug)
					  Rprintf("Returning %d\n",cl);
#endif
					return cl;
				}
				cl--;
			} else {
				cdf += fu;
				if(u <= cdf) {
#ifdef _DBG_
					if (debug)
					  Rprintf("Returning %d\n",ul);
#endif
					return ul;
				}
				ul++;
			}
		} while( (cl >= lstar || ul <= ustar) && u > cdf);

		Rprintf("Functional error: attempt to draw a value with");
		Rprintf("zero probability. Summary:\n");
		Rprintf("(cl>=lstar,ul<=ustar,u<=cdf,(cl>=lstar||ul<=ustar)&&u<=cdf) = ");
		Rprintf("%d,%d,%d %d\n",cl >= lstar,ul <= ustar,
			u <= cdf,(cl >= lstar || ul <= ustar) && u <= cdf);
		Rprintf("U <= CDF <=> %g<=%g\n",u,cdf); 

		double q = 0.0;
		for(i=lstar;i<=ustar;i++) {
			q = fvec[i-li];	
			if(i==cl || i==ul)  
				Rprintf("<%f> ",q);
			else if(i==nu)
				Rprintf("(%f) ",q);
			else
				Rprintf("%f ",q);
		}
		Rprintf("=%f\n",q);
		error("We should never get here in rnoncenhypgeo: %f,%f,%f,%f [%d,%d,%d] %d,%d [%f < %f].",m1,n1,n2,psi,li,nu,ui,lstar,ustar,u,cdf);

	//} else if(nu < l) {
	} else if(nu <= l) {

		//Our mode is smaller than the valid range. Only need to count up.

		// I think this is a potential problem -- 
		// nu failed the condition (nu<l) so this gives positive
		// probability to i=nu+1,...,li-1 doesn't it?

#ifdef _DBG_
		Rprintf("Mode less than lower limit: %d<=%d\n",nu,li);
		Rprintf("Parameters (n1,n2,m1,psi)=(%f,%f,%f,%f)\n\n",n1,n2,m1,psi);
		Rprintf("Resetting mode to lower limit...\n");
#endif
		nu = li;

		double f = 0.0;
		for(i=nu+1;i<=ui;i++) {
			f = f+log(ru_V22(m1,n1,n2,psi,i));
			if(i>=li)
				fvec[i-li] = f;
		}
		double S = logsumexp_V22(fvec,(ui+1)-li);
		double u = unif_rand();
		double cdf = 0.0;
		for(i=0;i<(ui+1)-li;i++) {
			cdf += exp(fvec[i]-S);
			if(u <= cdf) break;
		}		
#ifdef _DBG_
		Rprintf("Returning %d\n",i+li);
#endif
		return i+li;

#ifdef _DBG_
		Rprintf("I would return %d with limits (li=%d,ui=%d)\n",i+li,li,ui);
		Rprintf("but I'm calling the original function instead though...\n");
#endif
		double retval = jims_original_rnchg(m1,n1,n2,psi,delta,fvec,1);// Call with debugging on...
#ifdef _DBG_
		Rprintf("...instead I am gonna return %d\n",(int)retval);
#endif
		return retval;

	} else {
		//Our mode is larger than the valid range. Only need to count down.
	
#ifdef _DBG_
		Rprintf("Mode greater than upper limit: %d>%d\n",nu,ui);
		Rprintf("Parameters (n1,n2,m1,psi)=(%f,%f,%f,%f)\n\n",n1,n2,m1,psi);
		Rprintf("Resetting mode to upper limit...\n");
#endif
		nu = ui;

		double f = 0.0;
		for (i=nu-1;i>=li;i--)
		{
			//f = f-log(rl_V22(m1,n2,n2,psi,i)); // BUG!!!
			f = f-log(rl_V22(m1,n1,n2,psi,i));
			if(i<=ui)
				fvec[i-li] = f;
		}
		double S = logsumexp_V22(fvec,(ui+1)-li);
		double u = unif_rand();
		double cdf = 0.0;
		for(i=ui;i>=li;i--) 
		{
			cdf += exp(fvec[i-li]-S);
			if (u <= cdf) break;
		}
#ifdef _DBG_
		Rprintf("Returning %d\n",i);
#endif
		return i;

#ifdef _DBG_
		Rprintf("I would return %d with limits (li=%d,ui=%d)\n",i,li,ui);
		Rprintf("but I'm calling the original function instead though...\n");
#endif
		double retval = jims_original_rnchg(m1,n1,n2,psi,delta,fvec,1);// Call with debugging on...
#ifdef _DBG_
		Rprintf("...instead I am gonna return %d\n",i);
#endif
		return retval;
		
	}
	error("Reached end of rnoncenhypgeo() without return value");
	return R_NegInf;
}

//double new_rnoncenhypgeo(double m1,double n1,double n2,double psi,double delta, double *fvec, const int debug) 
double rnoncenhypgeo(const double m1,
		     const double n1,
		     const double n2,
		     const double psi,
		     const double delta, 
		     double * const fvec, 
		     const int debug) 
{
	//Sanity check.
	if(n1<0 || n2<0 || m1<0 || psi<=0) 
		error("Invalid parameters: (n1,n2,m1,psi)=(%d,%d,%d,%g)",(int)n1,(int)n2,(int)m1,psi); 

	// PB: I added this check just in case:
	if(m1>(n1+n2)) 
		error("Error: Invalid Parameters %d>%d+%d\n",(int)m1,(int)n1,(int)n2);

	//Calculate bounds
	const double l = ((m1-n2)>0)?(m1-n2):0.0; 
	const double u = (n1>m1)?m1:n1;     
	
	//Easy case. :-)
	if(l==u){
#ifdef _DBG_
		if (debug)
		  Rprintf("(l==u) Returning %f\n",l);
#endif
		return l;
	}
	
	//Integer bounds for ease of typing
	const int ui = (int)u;
	const int li = (int)l;
	const int ui_p1 = ui+1;
	const int li_m1 = li-1;
	
	//Try to find a mode inside our bounds
	static const double five_sixths = (5.0/6.0);
	const double a = (psi-1.0);
	const double b = -((n1+m1+2)*psi+n2-m1);
	const double c = (psi*(n1+1.0)*(m1+1.0));
	const double q = -(b+gq_sign(b)*sqrt(b*b-4.0*a*c))/2;
	const double M1 = floor(c/q);
	const double M2 = floor(q/a);
	int i;
        int nu=ui; 
	// Use the mode-based tricks by default, unless it is not suitable:
	int use_tricks=1;
	int r_z_greater=0;
        // NOTE: nu and r_z_greater will ALWAYS be overwritten,
        // but their initializations keeps picky compilers happy.

	// For faster computation:
	const double n1p1_mult_m1p1 = (n1+1.0)*(m1+1.0);
	const double n1p1_plus_m1p1 = (n1+m1+2.0);
	const double n2mm1 = (n2-m1);	
	
	if(isnan(M1) && isnan(M2)){
		// r(z)=1 has no roots:
		// Check:
		if ((b*b-4.0*a*c)<0.0){
			// No roots -- fine, we can still proceed:
			const double r_ui = pdb_ri(n1p1_mult_m1p1,n1p1_plus_m1p1,n2mm1,psi,(double)ui);

			// NOTE: r_i_chk==1.0 should never be possible 
			// (o/w we would have a root in the range)

			if (r_ui>1.0){
				r_z_greater = 1;
				nu = ui;
			} else {
				r_z_greater = -1;
				nu = li;
			}
			// Not suitable for mode-based tricks:
//			use_tricks = 0;
		} else {
			error("Both roots are NA in rnoncenhypgeo(): psi=%g,n1=%g,n2=%g,m1=%g",psi,n1,n2,m1);
		}

	} else {

		//if ((M1<(l+1)) || (M1>u)) {
		//	if ((M2<(l+1)) || (M2>u)) {
		if ((M1<(l)) || (M1>u)) {
			if ((M2<(l)) || (M2>u)) {

				//Neither mode is in the range:
#ifdef _DBG_
			if (debug){
				Rprintf("Neither mode is in range [l+1,u]=[%d,%d]\n",li+1,ui);
				int d1 = (int)fabs(M1 > u ? M1-u : l-M1);
				int d2 = (int)fabs(M2 > u ? M2-u : l-M2);
				nu = (d2 < d1) ? M2 : M1;
				Rprintf("Byron's mode: %d\n",nu);
			}
#endif

				// Not 'peaked' enough for mode-based tricks:
//				use_tricks = 0;
	
				// Find the mode:

				// Check two end-points of the range:
				const double r_li_p1 = pdb_ri(n1p1_mult_m1p1,n1p1_plus_m1p1,n2mm1,psi,li+1.0);
				double r_ui;
				// Save a tiny bit of computation if possible:
				if (fabs(li+1.0-ui)<1e-5){
					r_ui = r_li_p1;
				} else {
					r_ui = pdb_ri(n1p1_mult_m1p1,n1p1_plus_m1p1,n2mm1,psi,(double)ui);
				}

				// NOTE: r_li_p1==1.0, and r_ui==1.0 should never be possible 
				// (o/w we would have a root in the range)

				if (r_li_p1>1.0){
					if (r_ui<1.0){
						// Note: Include equality to be cautious
						error("Impossible case! Roots: %g and %g, (r_li_p1,r_ui)=(%g,%g) [%d,%d]\n",M1,M2,r_li_p1,r_ui,li,ui);
					}
					// r(z)>=1 for all z \in [l+1,u]
#ifdef _DBG_
					if (debug)
						Rprintf("r(z)>=1 for all z in [%d,%d]\n",l+1,u);
#endif
					r_z_greater = 1;
					// Mode is the upper limit:
					nu = ui;

				} else {
					if (r_ui>1.0){
						// Note: Include equality to be cautious
						error("Impossible case! Roots: %g and %g, (r_li_p1,r_ui)=(%g,%g) [%d,%d]\n",M1,M2,r_li_p1,r_ui,li,ui);
					}
					// r(z)<=1 for all z \in [l+1,u]
#ifdef _DBG_
					Rprintf("r(z)<=1 for all z in [%d,%d]\n",li+1,ui);
#endif	
					r_z_greater = -1;
					// Mode is the lower limit:
					nu = li;
				} // END mode finding
			} else {
				nu = (int)M2;
			} // END if((M2<(li+1))...) else (...)
		} else {
			nu = (int)M1;
		} // END if ((M1<(li+1))...) else (...)
	} // END if(isnan...) else (...)

	// See if we can use the tricks... 
	// (condition can be satisfied for nu=u, so we need a separate check)

	if (use_tricks){

	//if(nu > l && nu <= u) {
	if (nu>li_m1 && nu<ui_p1) {

		//Density calculation vector
		const double eps = (delta/10.0);
		double S = 1.0;
		int ustar, lstar;
#ifdef _DBG_
		if (debug)
		  Rprintf("Our mode is in the range, use the tricks...\n");
#endif

		fvec[nu-li] = 1.0;
		double r;
		//Above the mode
		for(i=nu+1;i<ui_p1;i++) {
			r = pdb_ri(n1p1_mult_m1p1,n1p1_plus_m1p1,n2mm1,psi,i);
			fvec[i-li] = fvec[i-1-li]*r;
			S      += fvec[i-li];
			if ((r<five_sixths)&&(fvec[i-li]<eps)) 
				break;
		}
		ustar = (i>ui)?ui:i;  //This is the value we should use for our upper bound
#ifdef _DBG_
		if (debug)
	          Rprintf("Upper bound: %d\n",ustar);
#endif	
		for(i=nu-1;i>li_m1;i--) {
			r = pdb_ri(n1p1_mult_m1p1,n1p1_plus_m1p1,n2mm1,psi,i+1.0);
			fvec[i-li]  = fvec[i+1-li]/r;
			S          += fvec[i-li];
			if ((r>five_sixths)&&(fvec[i-li]<eps)) 
				break;
		}
		lstar = (i<li)?li:i; //And the lower bound
#ifdef _DBG_
		if (debug)
	          Rprintf("Lower bound: %d\n",lstar);
#endif		
		//Now draw using cut-down-from-mode
		const double u = unif_rand();
		double cdf = fvec[nu-li]/S;
		if (u<cdf) {
#ifdef _DBG_
			if (debug)
			  Rprintf("Drawing the mode (%g<%g): Returning %d\n",u,cdf,nu);
#endif
#ifdef _NCHG_BOUNDS_CHECK_
		if ((nu<li)||(nu>ui))
			error("nchg draw outside of range (mode): %d [%d,%d]",nu,li,ui);
#endif
			return nu;
		}
		int cl=(nu-1);
		int ul=(nu+1);
		double fl,fu;
		do {
			fl = (cl>=lstar)?(fvec[cl-li]/S):0.0;
			fu = (ul<=ustar)?(fvec[ul-li]/S):0.0;
			if (fl>fu) {
				cdf += fl;
				if (u<=cdf) {
#ifdef _DBG_
					if (debug)
					  Rprintf("Returning %d\n",cl);
#endif
#ifdef _NCHG_BOUNDS_CHECK_
		if ((cl<li)||(cl>ui))
			error("nchg draw outside of range (cl): %d [%d,%d]",cl,li,ui);
#endif
					return cl;
				}
				cl--;
			} else {
				cdf += fu;
				if (u<=cdf) {
#ifdef _DBG_
					if (debug)
					  Rprintf("Returning %d\n",ul);
#endif
#ifdef _NCHG_BOUNDS_CHECK_
		if ((ul<li)||(ul>ui))
			error("nchg draw outside of range (ul): %d [%d,%d]",ul,li,ui);
#endif
					return ul;
				}
				ul++;
			}
		} while( ((cl>=lstar)||(ul<=ustar)) && (u>cdf));

		Rprintf("Functional error: attempt to draw a value with");
		Rprintf("zero probability. Summary:\n");
		Rprintf("(cl>=lstar,ul<=ustar,u<=cdf,(cl>=lstar||ul<=ustar)&&u<=cdf) = ");
		Rprintf("%d,%d,%d %d\n",cl >= lstar,ul <= ustar,
			u <= cdf,(cl >= lstar || ul <= ustar) && u <= cdf);
		Rprintf("U <= CDF <=> %g<=%g\n",u,cdf); 

		double q = 0.0;
		for(i=lstar;i<=ustar;i++) {
			q = fvec[i-li];	
			if(i==cl || i==ul)  
				Rprintf("<%f> ",q);
			else if(i==nu)
				Rprintf("(%f) ",q);
			else
				Rprintf("%f ",q);
		}
		Rprintf("=%f\n",q);
		error("We should never get here in rnoncenhypgeo: %f,%f,%f,%f [%d,%d,%d] %d,%d [%f < %f].",m1,n1,n2,psi,li,nu,ui,lstar,ustar,u,cdf);

	} // END if ((nu>=l)||(nu<u))
	} // END if (use_tricks)

	// Initial value for r(i) recursion:
	fvec[nu-li] = 1.0;
	double f = 0.0;

	if (r_z_greater<0){

		// Mode on lower boundary -- count up:
		for(i=nu+1;i<ui_p1;i++) {
			f += log(pdb_ri(n1p1_mult_m1p1,n1p1_plus_m1p1,n2mm1,psi,i));
			fvec[i-li] = f;
		}
		double S = logsumexp_V22(fvec,(ui+1)-li);
		const double u = unif_rand();
		double cdf = 0.0;
		const int ui_p1_mli = ui+1-li;
		for(i=0;i<ui_p1_mli;i++) {
			cdf += exp(fvec[i]-S);
			if(u <= cdf){
#ifdef _DBG_
				Rprintf("Returning %d -- valid: [%d,%d]\n\n",i+li,li,ui);
#endif
#ifdef _NCHG_BOUNDS_CHECK_
		if (((i+li)<li)||((i+li)>ui))
			error("nchg draw outside of range (r(z)<=1): %d [%d,%d]",i+li,li,ui);
#endif
 				return i+li;
			}
		}
		if (isnan(cdf)){
			Rprintf("S = %g\n",S);
			Rprintf("fvec:\n");
			for (i=0;i<ui_p1_mli;i++)
				Rprintf("%g\t",fvec[i]);
		}
		error("No mode found (li)! U<=CDF <=> %g<=%g\n",u,cdf);

	} else {

		// Mode on upper boundary -- count down:
		for (i=nu-1;i>li_m1;i--) {
			f -= log(pdb_ri(n1p1_mult_m1p1,n1p1_plus_m1p1,n2mm1,psi,i+1.0));
			fvec[i-li] = f;
		}
		double S = logsumexp_V22(fvec,(ui+1)-li);
		const double u = unif_rand();
		double cdf = 0.0;
		for(i=ui;i>li_m1;i--) {
			cdf += exp(fvec[i-li]-S);
			if (u <= cdf){
#ifdef _DBG_
				Rprintf("Returning %d -- valid: [%d,%d]\n\n",i,li,ui);
#endif
#ifdef _NCHG_BOUNDS_CHECK_
		if ((i<li)||(i>ui))
			error("nchg draw outside of range (r(z)>=1): %d [%d,%d]",i,li,ui);
#endif
				return i;
			}
		}
		if (isnan(cdf)){
			Rprintf("S = %g\n",S);
			Rprintf("fvec:\n");
			for (i=0;i<(ui+1-li);i++)
				Rprintf("%g\t",fvec[i]);
		}
		error("No mode found (ui)! U<=CDF <=> %g<=%g [%d,%d]\n",u,cdf,li,ui);
	} // END if (r_z_greater) {...} else {...}
  // Should never reach here:
  error("Reached end of new_rnoncenhypgeo() without return value [impossible]");
  return R_NegInf;
}
