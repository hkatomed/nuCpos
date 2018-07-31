#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(hba_3)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(localhba_3)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nucpos_1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nucpos_2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nucpos2_1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nucpos2_2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"hba_3",      (DL_FUNC) &F77_NAME(hba_3),       9},
    {"localhba_3", (DL_FUNC) &F77_NAME(localhba_3), 33},
    {"nucpos_1",   (DL_FUNC) &F77_NAME(nucpos_1),   13},
    {"nucpos_2",   (DL_FUNC) &F77_NAME(nucpos_2),   13},
    {"nucpos2_1",  (DL_FUNC) &F77_NAME(nucpos2_1),  16},
    {"nucpos2_2",  (DL_FUNC) &F77_NAME(nucpos2_2),  16},
    {NULL, NULL, 0}
};

void R_init_nuCpos(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
