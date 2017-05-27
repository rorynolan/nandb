#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP nandb_ExpSmooth(SEXP, SEXP, SEXP);
extern SEXP nandb_ExpSmoothNaive(SEXP, SEXP);
extern SEXP nandb_ExpSmoothPillars(SEXP, SEXP);
extern SEXP nandb_ExpSmoothRows(SEXP, SEXP, SEXP);
extern SEXP nandb_MeanPillars(SEXP);
extern SEXP nandb_MedianFilterB(SEXP, SEXP, SEXP, SEXP);
extern SEXP nandb_MedianPillars(SEXP);
extern SEXP nandb_MedReflectExtend(SEXP, SEXP, SEXP);
extern SEXP nandb_MedReflectExtendRows(SEXP, SEXP, SEXP);
extern SEXP nandb_ReflectIndexMed(SEXP, SEXP, SEXP);
extern SEXP nandb_Smooth(SEXP);
extern SEXP nandb_SmoothFilterB(SEXP, SEXP, SEXP, SEXP);
extern SEXP nandb_SpreadSpecificHelper(SEXP, SEXP, SEXP);
extern SEXP nandb_VarPillars(SEXP);
extern SEXP nandb_WhichIntervalC(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"nandb_ExpSmooth",            (DL_FUNC) &nandb_ExpSmooth,            3},
  {"nandb_ExpSmoothNaive",       (DL_FUNC) &nandb_ExpSmoothNaive,       2},
  {"nandb_ExpSmoothPillars",     (DL_FUNC) &nandb_ExpSmoothPillars,     2},
  {"nandb_ExpSmoothRows",        (DL_FUNC) &nandb_ExpSmoothRows,        3},
  {"nandb_MeanPillars",          (DL_FUNC) &nandb_MeanPillars,          1},
  {"nandb_MedianFilterB",        (DL_FUNC) &nandb_MedianFilterB,        4},
  {"nandb_MedianPillars",        (DL_FUNC) &nandb_MedianPillars,        1},
  {"nandb_MedReflectExtend",     (DL_FUNC) &nandb_MedReflectExtend,     3},
  {"nandb_MedReflectExtendRows", (DL_FUNC) &nandb_MedReflectExtendRows, 3},
  {"nandb_ReflectIndexMed",      (DL_FUNC) &nandb_ReflectIndexMed,      3},
  {"nandb_Smooth",               (DL_FUNC) &nandb_Smooth,               1},
  {"nandb_SmoothFilterB",        (DL_FUNC) &nandb_SmoothFilterB,        4},
  {"nandb_SpreadSpecificHelper", (DL_FUNC) &nandb_SpreadSpecificHelper, 3},
  {"nandb_VarPillars",           (DL_FUNC) &nandb_VarPillars,           1},
  {"nandb_WhichIntervalC",       (DL_FUNC) &nandb_WhichIntervalC,       2},
  {NULL, NULL, 0}
};

void R_init_nandb(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
