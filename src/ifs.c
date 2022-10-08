/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-2002	S. M. Iacus
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * Exports
 *	ifs_df(...)
 *  ifs_df_flex(...)
 *
 * to be called as  .C(.)  in ../R/ifs.R
 */

/*
   The relevant part of the theory of IFS related to
   this subject can be found in 
   Forte, B. and VRSCAY, E.R. (1995): "Solving the inverse 
   problem for measures using iterated function systems: a 
   new approach", Adv. Appl. Prob., 27, 800-820. Statistical
   applications of IFS to distribution function estimation
   are in Iacus, S.M. and La Torre, D. (2002).
   
   The following functional in the space of measures on [0,1]
   is used to approximate (or estimate) a measure (or a 
   distribution function):

             n
    (k)     __         -1
   T  (F) = >   p * F(w (x)) 
            --   i     i
            i=1

   where 
   { w_i, are affine maps of the form w_i = a_i * x + s_i }
   { p_i, is a ``probability distribution''}
   
   The p_i are solution of the following minimization problem
   of the quadratic form
   
   x'Qx + b'x, over sum_i x_i = 1, x_i >= 0 
   
   in which the terms of the matrix `Q' and the terms of the 
   vector `b' are built as in Forte and VRSCAY (1995).
   The choice of the family of the maps and the accuracy of the
   solution of the minimal problem are crucial to have good
   statistical estimates of the distribution function. For details
   on this see the forthcoming paper of Iacus and La Torre (2002).
   
   Thus, with good estimates of the p_i's the IFS funcional can
   be iterated starting from any distribution function on [0,1],
   for example the Uniform, and by the Banach theorem, the fixed
   point of the functional converges to F poinwisely.
   
   The estimator proposed is at least as efficient as the empirical
   distribution function and, in some norms, it is even better for 
   fixed sample size. 
   The maps proposed in Iacus and La Torre (2002) are choosen in
   a statistical perspective and thus the resulting estimator is
   generally highly efficient.
   
   IFS bases estimators are also continuous functions even tough
   not everywhere differentiable.

*/  

#include <R.h>
#include <Rmath.h>
#include <R_ext/Boolean.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Complex.h>

#define max(a, b) (a > b ? a : b) 
#define min(a, b) (a < b ? a : b) 

SEXP ifsm_setQF(SEXP X, SEXP s, SEXP a);

/* IFS estimators */      
SEXP ifs_df(SEXP x, SEXP p, SEXP s, SEXP a, SEXP k);
SEXP ifs_df_flex(SEXP x, SEXP p, SEXP s, SEXP a, SEXP k, SEXP f, SEXP rho);

SEXP ifs_ft(SEXP x, SEXP p, SEXP s, SEXP a, SEXP k);

/* Function needed to build the terms of the quadratic form
   to be minimized
*/   
SEXP ifs_setQF(SEXP mu, SEXP s, SEXP a, SEXP n);

/* 
   ATTN:   The following two are internal functions 
           NOT TO BE CALLED directly by the user !!! 
*/

double IFS( double x, int k);
double IFSflex( double x, int k, SEXP f, SEXP rho);

Rcomplex ifs_FT( Rcomplex *x, int k);


/* Two functions adapted from the code in man/doc/R-exts */
SEXP mkans(double x);
double feval(double x, SEXP f, SEXP rho);


/* Global variables 

  ps	: stores the coefficients p_i 
  cs	: stores the terms `s' in w_i = s_i * x + a_i
  ca    : stores the terms `a' in w_i = s_i * x + a_i
  mm    : stores the vector of moments of the r.v.
  nps   : the number of coefficient of the IFS

*/
  
double *ps  = NULL;  
double *cs  = NULL; 
double *ca  = NULL;
double *mm  = NULL;
int 	nps = 0;
Rboolean	firstiter = FALSE;

/* 
   IFS:
   ----
   This routine is charged to iterate the IFS. This
   is called by ifs_distfunc.
   It is assumed that the parameters `cs', `ca' and
   `ps' have been previously initialiazed 
   (ifs_df does the job).   
   Author: S. M. Iacus, Jan 9th 2002
   
   parameters:
   -----------
   x : where to estimate F
   k : number of iterations
   
   on return:
   ----------
   
   the estimated value

*/

double IFS( double x, int k)
{
  int i;
  double retval = 0;
  
  if(x<=0)
   return(0);
   
  if(x>=1) 
   return(1);

  if(k == 0)
   return(x);
  
  for(i=0; i< nps; ++i){
   if(cs[i] != 0)
    retval += (double)(ps[i] * IFS( (x - ca[i])/cs[i], k-1 ) ); 
  }
   
  return( retval ); 
}

/* 
   IFSflex:
   --------
   This routine is charged to iterate the IFS. This
   is called by ifs_df_flex.
   This function is the same as IFS but is more
   flexible in that the starting point `f' can
   be specified by the user
   It is assumed that the parameters `cs', `ca' and
   `ps' have been previously initialiazed 
   (ifs_distfunc does the job).   
   Author: S. M. Iacus, Jan 9th 2002
   
   
   parameters:
   -----------
   x : where to estimate F
   k : number of iterations
   
   on return:
   ----------
   the estimated value

*/

double IFSflex( double x, int k, SEXP f, SEXP rho )
{
  int i;
  double retval = 0;
  
  if(x<=0)
   return(0);
   
  if(x>=1) 
   return(1);

  if(k == 0)
   return(x);
  
  if(firstiter){
   for(i=0; i< nps; ++i){
    if(cs[i] != 0)
     retval += (double)(ps[i] * IFS(  feval( (x - ca[i])/cs[i] , f, rho)  , k-1 )); 
   }
   firstiter = FALSE; 
  }
  else{
   for(i=0; i< nps; ++i){
    if(cs[i] != 0)
     retval += (double)(ps[i] * IFS(  (x - ca[i])/cs[i] , k-1 ));    
   }
  }
   
  return( retval ); 
}




/* ifs_df:
   ------
   Routine that starts up paramters and calls IFS
   Author: S. M. Iacus, Jan 9th 2002
   
   parameters:
   -----------
   
   x   : where to estimate F
   p   : the n coefficients
   s   : the n coefficients of the maps w = s*x+a
   a   : the n coefficients of the maps w = s*x+a
   k   : number of iterations of the functional

   on return:
   ----------
   
   the estimated value of `F' at `x'
   
*/


SEXP ifs_df(SEXP x, SEXP p, SEXP s, SEXP a, SEXP k)
{
  SEXP ans;
  double *ics, *value;
  int *kappa, ncs;
  
  if(!isNumeric(x)) error("`x' must be numeric");
  if(!isNumeric(p)) error("`p' must be numeric");
  if(!isNumeric(s)) error("`s' must be numeric");
  if(!isNumeric(a)) error("`a' must be numeric");
  if(!isInteger(k)) error("`k' must be an integer");
 
  PROTECT(x = AS_NUMERIC(x));
  PROTECT(p = AS_NUMERIC(p));
  PROTECT(s = AS_NUMERIC(s));
  PROTECT(a = AS_NUMERIC(a));
  PROTECT(k = AS_INTEGER(k));
  
  nps = LENGTH(p); 
  ps = NUMERIC_POINTER(p);
  cs = NUMERIC_POINTER(s);
  ca = NUMERIC_POINTER(a);

  if( (ncs = LENGTH(s)) != LENGTH(a) )
   error("`a' and `s' must have same length");

  if(ncs != nps)
   error("`p', `a' and `s' must have same length");
    
  kappa   = INTEGER_POINTER(k);
  ics = NUMERIC_POINTER(x);
  
  PROTECT(ans = NEW_NUMERIC(1));
  value = NUMERIC_POINTER(ans);
  
  *value = IFS( *ics, *kappa);

  UNPROTECT(6);
  return(ans);
}


   
/* ifs_df_flex:
   ------------
   flexible version of ifs_df. This routine
   initialize parameters and calls IFSflex.
   Author: S. M. Iacus, Jan 9th 2002
   
   parameters:
   -----------
   
   x   : where to estimate F
   p   : the n coefficients
   s   : the n coefficients of the maps w = s*x+a
   a   : the n coefficients of the maps w = s*x+a
   k   : number of iterations of the functional
   f   : the initial value in the space of
         distribution functions on [0,1]
   rho : the environment where `f' is

   on return:
   ----------
   
   the estimated value of `F' at `x'
   
*/

SEXP ifs_df_flex(SEXP x, SEXP p, SEXP s, SEXP a, SEXP k, SEXP f, SEXP rho)
{
  SEXP ans;
  double *ics, *value;
  int *kappa, ncs;
  
  if(!isNumeric(x)) error("`x' must be numeric");
  if(!isNumeric(p)) error("`p' must be numeric");
  if(!isNumeric(s)) error("`s' must be numeric");
  if(!isNumeric(a)) error("`a' must be numeric");
  if(!isInteger(k)) error("`k' must be an integer");
  if(!isFunction(f)) error("`f' must be a (distribution) function");
  if(!isEnvironment(rho)) error("`rho' should be an environment");
 
  PROTECT(x = AS_NUMERIC(x));
  PROTECT(p = AS_NUMERIC(p));
  PROTECT(s = AS_NUMERIC(s));
  PROTECT(a = AS_NUMERIC(a));
  PROTECT(k = AS_INTEGER(k));
  
  nps = LENGTH(p); 
  ps = NUMERIC_POINTER(p);
  cs = NUMERIC_POINTER(s);
  ca = NUMERIC_POINTER(a);

  if( (ncs = LENGTH(s)) != LENGTH(a) )
   error("`a' and `s' must have same length");

  if(ncs != nps)
   error("`p', `a' and `s' must have same length");
    
  kappa   = INTEGER_POINTER(k);
  ics = NUMERIC_POINTER(x);
  
  PROTECT(ans = NEW_NUMERIC(1));
  value = NUMERIC_POINTER(ans);
  
  firstiter = TRUE;
  
  *value = IFSflex( *ics, *kappa, f, rho );

  UNPROTECT(6);
  return(ans);
}


/* 
   FT:
   ----
   This routine is charged to iterate the IFS in the
   space of Fourier transforms. This
   is called by ifs_ft.
   It is assumed that the parameters `cs', `ca' and
   `ps' have been previously initialiazed 
   (ifs_ft does the job).   
   Author: S. M. Iacus, Jan 31st 2002
   Fixed on Feb 22nd 2002. FT return an Rcomplex
   instead of *Romcplex. Thanks to Brian Ripley.
   
   parameters:
   -----------
   x : where to estimate FT
   k : number of iterations
   
   on return:
   ----------
   
   the estimated complex value

*/


Rcomplex FT(Rcomplex *x, int k)
{
  int i;
  Rcomplex tempval;
  Rcomplex retval;
  Rcomplex y;
  double u, a, b, c, d;
  
  retval.r = 0;
  retval.i = 0;
  
  if(k == 0){
   for(i=0; i< nps; ++i){
    if(ps[i] != 0){
     u = x->r * cs[i];
     y.r = sin(u) / u ;
     y.i = ( cos(u) - 1) / u ;
     a = ps[i] * cos( x->r * ca[i] );
     b = -ps[i] * sin( x->r * ca[i] );
     c = y.r;
     d = y.i;
     retval.r +=  a*c-b*d;
     retval.i += b*c+a*d;
     }
   }
  }
  else{
   for(i=0; i< nps; ++i){
    if(ps[i] != 0){
     y.r  = x->r * cs[i];
     y.i = 0;
     tempval = FT( &y , k-1);
     a = ps[i] * cos( x->r * ca[i] );
     b = -ps[i] * sin( x->r * ca[i] );
     c = tempval.r;
     d = tempval.i;
     retval.r +=  a*c-b*d;
     retval.i += b*c+a*d;
    }
   }
  } 
  return( retval ); 
}




/* ifs_ft:
   ------
   Routine that starts up paramters and calls FT
   Author: S. M. Iacus, Jan 31st 2002
   
   parameters:
   -----------
   
   x   : where to estimate the FT of F
   p   : the n coefficients
   s   : the n coefficients of the maps w = s*x+a
   a   : the n coefficients of the maps w = s*x+a
   k   : number of iterations of the functional

   on return:
   ----------
   
   the estimated complex value of `FT_F(x)' at `x'
   
   where FT_F(x) = int_0^1 e^{-iux} dF(u) 
   
*/

SEXP ifs_ft(SEXP x, SEXP p, SEXP s, SEXP a, SEXP k)
{
  SEXP ans;
  double *ics;
  Rcomplex *value, temp;
  int *kappa, ncs;
  
  if(!isNumeric(x)) error("`x' must be numeric");
  if(!isNumeric(p)) error("`p' must be numeric");
  if(!isNumeric(s)) error("`s' must be numeric");
  if(!isNumeric(a)) error("`a' must be numeric");
  if(!isInteger(k)) error("`k' must be an integer");
 
  PROTECT(x = AS_NUMERIC(x));
  PROTECT(p = AS_NUMERIC(p));
  PROTECT(s = AS_NUMERIC(s));
  PROTECT(a = AS_NUMERIC(a));
  PROTECT(k = AS_INTEGER(k));
  
  nps = LENGTH(p); 
  ps = NUMERIC_POINTER(p);
  cs = NUMERIC_POINTER(s);
  ca = NUMERIC_POINTER(a);

  if( (ncs = LENGTH(s)) != LENGTH(a) )
   error("`a' and `s' must have same length");

  if(ncs != nps)
   error("`p', `a' and `s' must have same length");
    
  kappa   = INTEGER_POINTER(k);
  ics = NUMERIC_POINTER(x);
  
  PROTECT(ans = NEW_COMPLEX(1));
  value = COMPLEX_POINTER(ans);
  
  firstiter = TRUE;
   
  temp.r = *ics;
  temp.i = 0;
  value[0] = FT( &temp, *kappa);

  UNPROTECT(6);
  return(ans);
}


/*
   This function is charged to build the Quadratic Form
   to be minimized but any optim algorithm.
   
   Reference:  Forte, B. and VRSCAY, E.R. (1995): 
   "Solving the inverse problem for measures using iterated 
   function systems: a new approach", Adv. Appl. Prob., 27, 
   800-820.
   
   Ingredients are:
   
   mu : the moments' vector of the target measure of the 
        r.v. X or the estimated moments, s.t. mu[0] = 1, 
        m[1] = mean(X), and so forth. At most M terms.
   s  : coefficientfs of the affine maps w_i = s_i x + a_i
        here s[0] = s_1, etc. At least N terms, N < M. At
        most M terms.        
   a  : coefficientfs of the affine maps w_i = s_i x + a_i
        here a[0] = a_1, etc. At least N terms, N < M. At
        most M terms.        
   N  : the dimension of the output QF, i.e. the number
        of maps used in the IFS iterator.
        
   On exit, it returns the matrix (N by N) Q and the (N by 1)
   vector b of x'Q'x + b'x and the matrix A (M by N) of the
   linear operator h = A*m[1:M]
       
*/

SEXP ifs_setQF(SEXP mu, SEXP s, SEXP a, SEXP n)
{
    SEXP ans, names, Q, b, A;
    int m, *nc, ndim, mdim, i, j, ncs, nu, k;
  
    if(!isNumeric(mu)) error("`mu' must be numeric");
    if(!isNumeric(s)) error("`s' must be numeric");
    if(!isNumeric(a)) error("`a' must be numeric");
    if(!isInteger(n)) error("`n' must be an integer");
 
    PROTECT(mu = AS_NUMERIC(mu));
    PROTECT(s = AS_NUMERIC(s));
    PROTECT(a = AS_NUMERIC(a));
    PROTECT(n = AS_INTEGER(n));
  
    m = LENGTH(mu);
    nc = INTEGER_POINTER(n);
    ndim = *nc;
    mdim = m-1;
  
    if(mdim < ndim)
        error("`n' length is too high with respect to `mu' one");

    cs = NUMERIC_POINTER(s);
    ca = NUMERIC_POINTER(a);
    mm = NUMERIC_POINTER(mu);

    if( (ncs = LENGTH(s)) != LENGTH(a) )
        error("`a' and `s' must have same length");
  
    PROTECT(A = allocMatrix(REALSXP, mdim, ndim));
  
    for (i = 0; i < mdim; i++)
        for (j = 0; j < ndim; j++)
            REAL(A)[i + j * mdim] = 0.0;

    for(nu=0; nu < mdim; nu++)
        for(i=0; i< ndim ;i++)
            for(k=0; k <= nu+1; k++)
                REAL(A)[nu + i * mdim] += choose((double)(nu+1.),(double)k) * R_pow(cs[i],(double)k) * R_pow(ca[i],(double)(nu+1-k)) * mm[k];
    
    PROTECT(ans = allocVector(VECSXP,3));
    PROTECT(names = allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("Q"));
    SET_STRING_ELT(names, 1, mkChar("b"));
    SET_STRING_ELT(names, 2, mkChar("A"));

    SET_VECTOR_ELT(ans, 0, PROTECT(Q = allocMatrix(REALSXP, ndim, ndim)));

    for (i = 0; i < ndim; i++)
        for (j = 0; j < ndim; j++)
            REAL(Q)[i + j * ndim] = 0.0;


    for(i=0; i < ndim ; i++)
        for(j=0; j < ndim ; j++)
            for(nu=0; nu < mdim ; nu++)
                REAL(Q)[i + j * ndim] += REAL(A)[nu + i * mdim] *  REAL(A)[nu + j * mdim] / R_pow((double)(nu+1),2.);

    SET_VECTOR_ELT(ans, 1, PROTECT(b = allocVector(REALSXP, ndim)));
    
    for(i=0; i < ndim ; i++)
        REAL(b)[i] = 0.0;
   
    for(i=0; i < ndim ; i++){
        for(nu=0; nu < mdim ; nu++)
            REAL(b)[i] += mm[nu+1] * REAL(A)[nu + i * mdim] / R_pow((double)(nu+1),2.);
        REAL(b)[i] = -2.0 * REAL(b)[i];
    }
    
    SET_VECTOR_ELT(ans, 2, A);
 
    setAttrib(ans, R_NamesSymbol, names);
  
    UNPROTECT(9);
    return(ans);
}


double proc(double t, double *X, int nx);


/*
   This function is charged to build the Quadratic Form
   to be minimized but any optim algorithm.
   
   Reference:  Forte, B., Vrscay, E.R. (1994) Solving the 
   inverse problem for function/image approximation using iterated 
	function systems, I. Theoretical basis, Fractal, 2, 3 
	(1994), 325-334
   
   Ingredients are:
   
   X  : the vector of simulated trajectory
   s  : coefficientfs of the affine maps w_i = s_i x + a_i
        here s[0] = s_1, etc. At least N terms, N < M. At
        most M terms.        
   a  : coefficients of the affine maps w_i = s_i x + a_i
        here a[0] = a_1, etc. At least N terms, N < M. At
        most M terms.        
   N  : the dimension of the output QF, i.e. the number
        of maps used in the IFS iterator.
        
   On exit, it returns the matrix (N by N) Q and the (N by 1)
   vector b of x'Q'x + b'x 
       
*/

SEXP ifsm_setQF(SEXP X, SEXP s, SEXP a)
{
    SEXP ans, names, Q, b, mids, L1, L2, M1;
    int na, ns, n, N, i, j, ncs, nx;
 //   double l2=0, l1=0, m1=0, *xx;
    double l2=0, l1=0, m1=0;
    double tmin, tmax, Dt, t, S, tmp;

  
    if(!isNumeric(X)) error("`X' must be numeric");
    if(!isNumeric(s)) error("`s' must be numeric");
    if(!isNumeric(a)) error("`a' must be numeric");
 
    PROTECT(X = AS_NUMERIC(X));
    PROTECT(s = AS_NUMERIC(s));
    PROTECT(a = AS_NUMERIC(a));
  
    nx = LENGTH(X);
    na = LENGTH(a);
    ns = LENGTH(s);

    if( na != ns )
        error("`a' and `s' must have same length");

    n = na;
    N = 2*n;
  
    cs = NUMERIC_POINTER(s);
    ca = NUMERIC_POINTER(a);
    //xx = NUMERIC_POINTER(X);

    if( (ncs = LENGTH(s)) != LENGTH(a) )
        error("`a' and `s' must have same length");
  
    PROTECT(mids = allocVector(REALSXP,nx));
    PROTECT(ans = allocVector(VECSXP,5));
    PROTECT(names = allocVector(STRSXP, 5));
    SET_STRING_ELT(names, 0, mkChar("Q"));
    SET_STRING_ELT(names, 1, mkChar("b"));
    SET_STRING_ELT(names, 2, mkChar("L1"));
    SET_STRING_ELT(names, 3, mkChar("L2"));
    SET_STRING_ELT(names, 4, mkChar("M1"));

    SET_VECTOR_ELT(ans, 0, PROTECT(Q = allocMatrix(REALSXP, N, N)));

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            REAL(Q)[i + j * N] = 0.0;

    for(i=0; i<nx-1; i++){
        tmin = 0.0;
        tmax = 1.0;
        Dt = (tmax-tmin)/100.0;
        t = tmin+Dt/2.0;
        while(t < tmax){
            tmp = proc(t, REAL(X), nx);
            l2 += Dt*tmp*tmp;
            m1 += Dt*tmp;
            l1 += Dt*fabs(tmp);
            t += Dt;
        }
    }

    for(i=0; i < n ; i++){
        for(j=i; j < n ; j++){
            tmin = max(REAL(a)[i],REAL(a)[j]);
            tmax = min(REAL(a)[i]+REAL(s)[i],REAL(a)[j]+REAL(s)[j]);
            Dt = (tmax-tmin)/100.0;
            t = tmin+Dt/2.0;
            S = 0.0;
            while(t < tmax){
                S += Dt*proc((t-REAL(a)[i])/REAL(s)[i], REAL(X), nx)*proc((t-REAL(a)[j])/REAL(s)[j], REAL(X), nx);
                t += Dt;
            }
            REAL(Q)[i + j * N] = S;
            REAL(Q)[j + i * N] = S;
        }
    }
   

    for(i=0; i < n ; i++){
        for(j=i; j < n ; j++){
            tmin = REAL(a)[j];
            tmax = REAL(a)[j]+REAL(s)[j];
            Dt = (tmax-tmin)/100.0;
            t = tmin+Dt/2.0;
            S = 0.0;
            while(t < tmax){
                S += Dt*proc((t-REAL(a)[i])/REAL(s)[i], REAL(X), nx);
                t += Dt;
            }
            REAL(Q)[i + (j+n) * N] = S;
            REAL(Q)[(n+i) + j * N] = S;
        }
    }

    for(i=0; i < n ; i++){
        for(j=i; j < n ; j++){
            if((REAL(s)[i] + REAL(a)[i] > REAL(a)[j]) && (REAL(a)[i]< REAL(a)[j]+REAL(s)[j])){
                tmin = max(REAL(a)[i],REAL(a)[j]);
                tmax = min(REAL(a)[i]+REAL(s)[i],REAL(a)[j]+REAL(s)[j]);
                REAL(Q)[(n+i) + (j+n) * N] = tmax-tmin;
            }
        }
    }


    SET_VECTOR_ELT(ans, 1, PROTECT(b = allocVector(REALSXP, N)));

    for(i=0; i < N ; i++)
        REAL(b)[i] = 0.0;
       
    for(i=0; i < n ; i++){
        tmin = 0.0;
		tmax = 1.0;
		Dt = (tmax-tmin)/100.0;
		t = tmin+Dt/2.0;
		S = 0.0;
		while(t < tmax){
			S += Dt*proc(t, REAL(X), nx)*proc((t-REAL(a)[i])/REAL(s)[i], REAL(X), nx);
			t += Dt;
		}
		REAL(b)[i] = -2*S;
    }
   
	
	for(i=0; i < n ; i++){
		tmin = REAL(a)[i];
		tmax = REAL(a)[i] + REAL(s)[i];
		Dt = (tmax-tmin)/100.0;
		t = tmin+Dt/2.0;
		S = 0.0;
		while(t < tmax){
			S += Dt*proc(t, REAL(X), nx);
			t += Dt;
		}
		REAL(b)[i+n] = -2*S;
    }
   

    SET_VECTOR_ELT(ans, 2, PROTECT(L1 = allocVector(REALSXP, 1)));
    REAL(L1)[0] = l1;

    SET_VECTOR_ELT(ans, 3, PROTECT(L2 = allocVector(REALSXP, 1)));
    REAL(L2)[0] = l2;
  
    SET_VECTOR_ELT(ans, 4, PROTECT(M1 = allocVector(REALSXP, 1)));
    REAL(M1)[0] = m1;
    setAttrib(ans, R_NamesSymbol, names);
  
    UNPROTECT(11);
    return(ans);
}

double proc(double t, double *X, int nx){
    if((t<0.0) | (t>1.0))
	 return(0.0);
	return( X[ (int)(round(t*(nx-1))) ] ); 
}


static R_CMethodDef R_CDef[] = {
   {"ifs_df", (DL_FUNC)&ifs_df, 5},
   {"ifs_df_flex", (DL_FUNC)&ifs_df_flex, 7},
   {"ifs_setQF", (DL_FUNC)&ifs_setQF, 4},
   {"ifs_ft", (DL_FUNC)&ifs_ft, 5},
   {"ifsm_setQF", (DL_FUNC)&ifsm_setQF, 4},
   {NULL, NULL, 0},
};

void
R_init_ifs(DllInfo *info)
{
    R_registerRoutines(info, R_CDef, NULL, NULL, NULL);
    R_useDynamicSymbols(info, TRUE);
}




SEXP mkans(double x)
{
    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = x;
    UNPROTECT(1);
    return ans;
}

double feval(double x, SEXP f, SEXP rho)
{
    double val;
    SEXP R_fcall;    
    defineVar(install("x"), mkans(x), rho);
    PROTECT(R_fcall = lang2(f, mkans(x)));
    val = *NUMERIC_POINTER(eval(R_fcall, rho));
    UNPROTECT(1);
    return(val);
}


 
