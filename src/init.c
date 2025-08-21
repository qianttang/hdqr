#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/
#define FDEF(name){#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}

/* .Fortran calls */
void F77_NAME(lqr_hd)(double *alpha, double *lam2, double *hval, 
  int *nobs, int *nvars, double *x, double *y, double *tau, int *jd, int *pfncol, 
  double *pf, double *pf2, int *dfmax, int *pmax, int *nlam, double *flmin, 
  double *ulam, double *eps, int *isd, int *maxit, int *nalam, 
  double *b0, double *beta, int *ibeta, int *nbeta, double *alam, 
  int *npass, int *jerr, double *sigma, int *is_exact);

static R_NativePrimitiveArgType lqr_hd_t[] = {REALSXP, REALSXP, REALSXP, 
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, 
  REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, 
  INTSXP, INTSXP, REALSXP, INTSXP};

static R_FortranMethodDef FortranEntries[] = {
    FDEF(lqr_hd) ,
    {NULL, NULL, 0}
};


void R_init_hdqr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
