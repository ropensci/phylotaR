#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP cFindPrids(SEXP, SEXP, SEXP);
extern SEXP cGetNdmtrx(SEXP, SEXP, SEXP);
extern SEXP cGetNdPrids(SEXP, SEXP);
extern SEXP cGetNdPtids(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cFindPrids",  (DL_FUNC) &cFindPrids,  3},
    {"cGetNdmtrx",  (DL_FUNC) &cGetNdmtrx,  3},
    {"cGetNdPrids", (DL_FUNC) &cGetNdPrids, 2},
    {"cGetNdPtids", (DL_FUNC) &cGetNdPtids, 2},
    {NULL, NULL, 0}
};

void R_init_phylotaR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}