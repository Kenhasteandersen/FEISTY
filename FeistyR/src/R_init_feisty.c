#ifndef R_R_H
# include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif


#include <Rinternals.h>
#include <stdlib.h> // for NULL

/* register native routines ------------------------------------------------ */


/* NO .C calls */
/* NO .Call calls */

/* .Fortran calls */
void F77_NAME(initfeisty)(void (* steadyparms)(int *, double *));     /*customize*/
void F77_NAME(initfeistysetupbasic)(void (* steadyparms)(int *, double *));        /* setupbasic (Petrik et al., 2019) */
void F77_NAME(initfeistysetupbasic2)(void (* steadyparms)(int *, double *));       /* setupbasic2 */
void F77_NAME(initfeistysetupvertical)(void (* steadyparms)(int *, double *));     /* setupvertical (van Denderen et al., 2021)  */
//void F77_NAME(initfeistysetupvertical2)(void (* steadyparms)(int *, double *));    /* setupvertical2                  */
void F77_NAME(runfeisty) (int *, double *, double *, double *, double *, int *);

R_FortranMethodDef FEntries[] = {
    {"initfeisty",    (DL_FUNC) &F77_SUB(initfeisty),   1},
    {"initfeistysetupbasic",    (DL_FUNC) &F77_SUB(initfeistysetupbasic),   1},
    {"initfeistysetupbasic2",    (DL_FUNC) &F77_SUB(initfeistysetupbasic2),   1},
    {"initfeistysetupvertical",    (DL_FUNC) &F77_SUB(initfeistysetupvertical),   1},
   // {"initfeistysetupvertical2",    (DL_FUNC) &F77_SUB(initfeistysetupvertical2),   1},
    {"runfeisty",     (DL_FUNC) &F77_SUB(runfeisty),    6},
    {NULL, NULL, 0}
};

/* Initialization ---------------------------------------------------------- */
void R_init_feisty(DllInfo *dll) {

  R_registerRoutines(dll, NULL, NULL, FEntries, NULL);

  // the following line protects against accidentially finding entry points

  R_useDynamicSymbols(dll, FALSE); // disable dynamic searching
}
