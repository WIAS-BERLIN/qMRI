#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#include <Rinternals.h> // for SEXP
#include <R_ext/RS.h>
void F77_NAME(gethani)(double* x, double* y, int* kern, double* value,
  double* wght, double* eps, double* bw);
void F77_NAME(hg1f1)(double* a, double* b, double* z, int* n, double* fz);
void F77_NAME(mediansm)(double* y, int* mask, int* n1, int* n2, int* n3,
  int* ind, int* nind, double* work, int* ncores, double* yout);
void F77_NAME(paramw3)(double* h, double* vext, int* indn, double* w, int* n);
void F77_NAME(pvawsme)(double* y, double* yd, int* pos, int* nv, int* nvd, int* nd, int* n1, int* n2,
  int* n3, double* hakt, double* lambda, double* theta, double* bi, double* bin,
  double* thnew, double* ydnew, double* invcov, int* ncores, double* spmin, double* lwght,
  double* wght, double* swjy, double* swjd, int* np1, int* np2, int* np3);
void F77_NAME(pvawsm2)(double* y, int* pos, int* nv, int* nvd, int* n1, int* n2,
  int* n3, double* hakt, double* lambda, double* theta, double* bi, double* bin,
  double* thnew, double* invcov, int* ncores, double* spmin, double* lwght,
  double* wght, double* swjy, int* np1, int* np2, int* np3);
void F77_NAME(qflashm0)(double* th, double* des, int* n, double* fval);
void F77_NAME(qflashm1)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(estatics3)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(qflashp0)(double* th, double* r2star, double* des, int* n,
  double* fval, double* grad);
void F77_NAME(estatics2)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(qflashp20)(double* th, double* r2star, double* des, int* n,
  double* fval, double* grad);
void F77_NAME(estatics1)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(vaws2)(double* y, int* mask, int* nv, int* n1, int* n2, int* n3,
  double* hakt, double* lambda, double* theta, double* s2, double* bi,
  double* thnew, double* s2new, int* ncores, double* lwght, double* wght,
  double* swjy);

static R_NativePrimitiveArgType gethani_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType hg1f1_t[]={REALSXP, REALSXP, REALSXP, INTSXP,
  REALSXP};
static R_NativePrimitiveArgType mediansm_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType paramw3_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
  INTSXP};
static R_NativePrimitiveArgType pvawsme_t[]={REALSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType pvawsm2_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType qflashm0_t[]={REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType qflashm1_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType estatics3_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType qflashp0_t[]={REALSXP, REALSXP, REALSXP, INTSXP,
  REALSXP, REALSXP};
static R_NativePrimitiveArgType estatics2_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP};
  static R_NativePrimitiveArgType qflashp20_t[]={REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP};
static R_NativePrimitiveArgType estatics1_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType vaws2_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};

static const R_FortranMethodDef fmethods[] = {
            {"gethani", (DL_FUNC) &gethani_ , 7, gethani_t},
            {"hg1f1", (DL_FUNC) &hg1f1_ , 5, hg1f1_t},
            {"mediansm", (DL_FUNC) &mediansm_ , 10, mediansm_t},
            {"paramw3", (DL_FUNC) &paramw3_ , 5, paramw3_t},
            {"pvawsme", (DL_FUNC) &pvawsme_ , 26, pvawsme_t},
            {"pvawsm2", (DL_FUNC) &pvawsm2_ , 22, pvawsm2_t},
            {"qflashm0", (DL_FUNC) &qflashm0_ , 4, qflashm0_t},
            {"qflashm1", (DL_FUNC) &qflashm1_ , 5, qflashm1_t},
            {"estatics3", (DL_FUNC) &estatics3_ , 5, estatics3_t},
            {"qflashp0", (DL_FUNC) &qflashp0_ , 6, qflashp0_t},
            {"qflashp20", (DL_FUNC) &qflashp20_ , 6, qflashp20_t},
            {"estatics2", (DL_FUNC) &estatics2_ , 5, estatics2_t},
            {"estatics1", (DL_FUNC) &estatics1_ , 5, estatics1_t},
            {"vaws2", (DL_FUNC) &vaws2_ ,17, vaws2_t},
            {NULL, NULL, 0,NULL}
};

void R_init_qMRI(DllInfo *dll)
         {
             R_registerRoutines(dll, NULL, NULL, fmethods , NULL);
             R_useDynamicSymbols(dll,FALSE);
             R_forceSymbols(dll,TRUE);
         }
