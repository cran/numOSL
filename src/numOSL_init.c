#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(caled_fort)(void *, void *, void *, void *, void *, void *, void *, 
                                 void *, void *, void *, void *, void *, void *, 
                                 void *, void *, void *, void *, 
                                 void *, void *);

extern void F77_NAME(caled1)(void *, void *, void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *);

extern void F77_NAME(calsgced_fort)(void *, void *, void *, void *, void *, void *, void *,
                                    void *, void *, void *, void *, void *);

extern void F77_NAME(decomp_fort)(void *, void *, void *, void *, void *, void *, void *, void *, void *, 
                                  void *, void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(fitgok)(void *, void *, void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *);

extern void F77_NAME(fitgrowth_fort)(void *, void *, void *, void *, void *, void *, void *, 
                                     void *, void *, void *, void *, void *);

extern void F77_NAME(calcsf)(void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *, void *);

extern void F77_NAME(mccam)(void *, void *, void *, void *, void *, 
                            void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(mcfmm2)(void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(mcfmm3)(void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(mcfmm4)(void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(mcmam3)(void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(mcmam4)(void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(apmamstd)(void *, void *, void *, void *, 
                               void *, void *, void *);

extern void F77_NAME(goodcomp)(void *, void *, void *, void *, 
                               void *, void *, void *);

extern void F77_NAME(comped)(void *, void *, void *, void *, void *, 
                             void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"caled_fort",     (DL_FUNC) &F77_NAME(caled_fort),     19},
    {"caled1",         (DL_FUNC) &F77_NAME(caled1),         18},
    {"calsgced_fort",  (DL_FUNC) &F77_NAME(calsgced_fort),  12},
    {"decomp_fort",    (DL_FUNC) &F77_NAME(decomp_fort),    17},
    {"fitgok",         (DL_FUNC) &F77_NAME(fitgok),         11},
    {"fitgrowth_fort", (DL_FUNC) &F77_NAME(fitgrowth_fort), 12},
    {"calcsf",         (DL_FUNC) &F77_NAME(calcsf),         10},
    {"mccam",          (DL_FUNC) &F77_NAME(mccam),          12},
    {"mcfmm2",         (DL_FUNC) &F77_NAME(mcfmm2),         12},
    {"mcfmm3",         (DL_FUNC) &F77_NAME(mcfmm3),         12},
    {"mcfmm4",         (DL_FUNC) &F77_NAME(mcfmm4),         12},
    {"mcmam3",         (DL_FUNC) &F77_NAME(mcmam3),         13},
    {"mcmam4",         (DL_FUNC) &F77_NAME(mcmam4),         13},
    {"apmamstd",       (DL_FUNC) &F77_NAME(apmamstd),        7},
    {"goodcomp",       (DL_FUNC) &F77_NAME(goodcomp),        7},
    {"comped",         (DL_FUNC) &F77_NAME(comped),         10},
    {NULL, NULL, 0}
};

void R_init_numOSL(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
