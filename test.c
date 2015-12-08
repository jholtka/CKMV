/*--------------------- revision history -----------------------------------
 * 21-11-95 (JAD): Created
 *--------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "oxexport.h"

#if defined (__WATCOMC__)
  /* and OpenWatcom Fortran */
  #define FORTRAN __fortran
#else
  /* e.g. gfortran */
  #define FORTRAN
  #define TEST test_
#endif

extern void FORTRAN TEST(double (FORTRAN *func)(double *), int *, double *, int *);

double FORTRAN CALL2C(double *pd)
{
    return 2 * *pd;
}

void OXCALL FnFortranTest(OxVALUE *rtn, OxVALUE *pv, int cArg)
{
    int i, fl;  double d;

    OxLibCheckType(OX_ARRAY, pv, 0, 2);

    fl = i = 0; d = 0;
    TEST(CALL2C, &fl, &d, &i); /* TEST changes these values to 1, -5, 1000 */

    OxSetInt( OxArray(pv,0), 0, fl);
    OxSetDbl( OxArray(pv,1), 0, d);
    OxSetInt( OxArray(pv,2), 0, i);
}

