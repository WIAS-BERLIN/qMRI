#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#include <Rinternals.h> // for SEXP
#include <R_ext/RS.h>
void F77_NAME(hg1f1)(double* a, double* b, double* z, int* n, double* fz);
void F77_NAME(estatics1)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(estatics1fixedr2)(double* th, double* r2star, double* des, int* n,
  double* fval, double* grad);
void F77_NAME(estatics2)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(estatics2fixedr2)(double* th, double* r2star, double* des, int* n,
  double* fval, double* grad);
void F77_NAME(estatics3)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(estatics3fixedr2)(double* th, double* r2star, double* des, int* n,
  double* fval, double* grad);

static R_NativePrimitiveArgType hg1f1_t[]={REALSXP, REALSXP, REALSXP, INTSXP,
  REALSXP};
static R_NativePrimitiveArgType estatics1_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType estatics1fixedr2_t[]={REALSXP, REALSXP, REALSXP, INTSXP,
  REALSXP, REALSXP};
static R_NativePrimitiveArgType estatics2_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType estatics2fixedr2_t[]={REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP};
static R_NativePrimitiveArgType estatics3_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType estatics3fixedr2_t[]={REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP};

static const R_FortranMethodDef fmethods[] = {
            {"hg1f1", (DL_FUNC) &hg1f1_ , 5, hg1f1_t},
            {"estatics1", (DL_FUNC) &estatics1_ , 5, estatics1_t},
            {"estatics1fixedr2", (DL_FUNC) &estatics1fixedr2_ , 6, estatics1fixedr2_t},
            {"estatics2", (DL_FUNC) &estatics2_ , 5, estatics2_t},
            {"estatics2fixedr2", (DL_FUNC) &estatics2fixedr2_ , 6, estatics2fixedr2_t},
            {"estatics3", (DL_FUNC) &estatics3_ , 5, estatics3_t},
            {"estatics3fixedr2", (DL_FUNC) &estatics3fixedr2_ , 6, estatics3fixedr2_t},
            {NULL, NULL, 0,NULL}
};

void R_init_qMRI(DllInfo *dll)
         {
             R_registerRoutines(dll, NULL, NULL, fmethods , NULL);
             R_useDynamicSymbols(dll,FALSE);
             R_forceSymbols(dll,TRUE);
         }
