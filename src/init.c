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
void F77_NAME(pvaws2)(double* y, int* mask, int* nv, int* nvd, int* n1, int* n2,
  int* n3, double* hakt, double* lambda, double* theta, double* bi, double* bin,
  double* thnew, double* invcov, int* ncores, double* spmin, double* lwght,
  double* wght, double* swjy, int* np1, int* np2, int* np3);
void F77_NAME(qflashm0)(double* th, double* des, int* n, double* fval);
void F77_NAME(qflashm1)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(qflashpl)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(qflashp0)(double* th, double* r2star, double* des, int* n,
  double* fval, double* grad);
void F77_NAME(qflashpl2)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(qflashpl3)(double* th, double* des, int* n, double* fval, double* grad);
void F77_NAME(vaws)(double* y, int* mask, int* nv, int* n1, int* n2, int* n3,
  double* hakt, double* lambda, double* theta, double* si2, double* bi,
  double* thnew, int* ncores, double* lwght, double* wght, double* swjy);
void F77_NAME(vawsext)(double* y, int* mask, int* nv, int* n1, int* n2, int* n3,
  double* yext, int* nve, double* hakt, double* lambda, double* theta, double* si2,
  double* bi, double* thnew, double* thext, int* ncores, double* lwght, double* wght,
  double* swjy, double* swjye);
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
  static R_NativePrimitiveArgType pvaws2_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
  static R_NativePrimitiveArgType qflashm0_t[]={REALSXP, REALSXP, INTSXP, REALSXP};
  static R_NativePrimitiveArgType qflashm1_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP};
  static R_NativePrimitiveArgType qflashpl_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP};
  static R_NativePrimitiveArgType qflashp0_t[]={REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP};
  static R_NativePrimitiveArgType qflashpl2_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP};
  static R_NativePrimitiveArgType qflashpl3_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP};
  static R_NativePrimitiveArgType vaws_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP, REALSXP};
  static R_NativePrimitiveArgType vawsext_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP};
  static R_NativePrimitiveArgType vaws2_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};

static const R_FortranMethodDef fmethods[] = {
            {"gethani", (DL_FUNC) &gethani_ , 7, gethani_t},
            {"hg1f1", (DL_FUNC) &hg1f1_ , 5, hg1f1_t},
            {"mediansm", (DL_FUNC) &mediansm_ , 10, mediansm_t},
            {"paramw3", (DL_FUNC) &paramw3_ , 5, paramw3_t},
            {"pvaws2", (DL_FUNC) &pvaws2_ , 22, pvaws2_t},
            {"qflashm0", (DL_FUNC) &qflashm0_ , 4, qflashm0_t},
            {"qflashm1", (DL_FUNC) &qflashm1_ , 5, qflashm1_t},
            {"qflashpl", (DL_FUNC) &qflashpl_ , 5, qflashpl_t},
            {"qflashp0", (DL_FUNC) &qflashp0_ , 6, qflashp0_t},
            {"qflashpl2", (DL_FUNC) &qflashpl2_ , 5, qflashpl2_t},
            {"qflashpl3", (DL_FUNC) &qflashpl3_ , 5, qflashpl3_t},
            {"vaws", (DL_FUNC) &vaws_ , 16, vaws_t},
            {"vawsext", (DL_FUNC) &vawsext_ , 20, vawsext_t},
            {"vaws2", (DL_FUNC) &vaws2_ ,17, vaws2_t},
            {NULL, NULL, 0,NULL}
};

void R_init_aws(DllInfo *dll)
         {
             R_registerRoutines(dll, NULL, NULL, fmethods , NULL);
             R_useDynamicSymbols(dll,FALSE);
             R_forceSymbols(dll,TRUE);
         }
