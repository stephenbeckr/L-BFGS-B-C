#include "lbfgsb.h"
#include <time.h>

int timer(double *ttime)
{
    clock_t temp;

/*     this routine computes cpu time in double precision; it makes use of */
/*     the intrinsic f90 cpu_time therefore a conversion type is */
/*     needed. */

/*           j.l morales  departamento de Matematicas, */
/*                        instituto Tecnologico Autonomo de Mexico */
/*                        mexico D.F. */

/*           j.l nocedal  department of Electrical Engineering and */
/*                        computer Science. */
/*                        northwestern University. Evanston, IL. USA */

/*                        january 21, 2011 */

/*     temp = (real) (*ttime); */
/*     *ttime = (double) temp; */
    temp    = clock();
    *ttime  = ((double) temp)/CLOCKS_PER_SEC;
    return 0;
} /* timer */

