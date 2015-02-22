/* Downloaded from http://www.netlib.org/clapack/cblas/
 * This is a merger of the 4 following cblas files:
 *   daxpy.c  dcopy.c  ddot.c  dscal.c
 * These files were long ago taken from the fortran source
 * using f2c.
 * Use this file for the simple, unoptimized implementation
 * of these basic BLAS functions. These are "level-1" BLAS
 * so optimization does not make a huge difference, so there
 * is not a major reason not to use this. If for some reason
 * you think these are a bottleneck of computation,
 * then link with an optimized BLAS (gotoBLAS, ATLAS, Intel MKL).
 * But, be warned that there some BLAS expect 32-bit/4 byte integers,
 * and some expect 64-bit/8-byte integers, and they are passed by
 * reference, so this can cause crashes if you get it wrong!
 *
 * Feb 28 2014, Stephen Becker
  
   According to wikipedia, most of netlib work is in the public domain,
   so we redistribute it here.
   See also http://www.netlib.org/clapack/faq.html#1.2
  
   LAPACK, CLAPACK and CLBAS are the work of many authors, most of them
   involved in the LAPACK user guide 3rd edition:
   @BOOK{laug,
      AUTHOR = {Anderson, E. and Bai, Z. and Bischof, C. and
                Blackford, S. and Demmel, J. and Dongarra, J. and
                Du Croz, J. and Greenbaum, A. and Hammarling, S. and
                McKenney, A. and Sorensen, D.},
      TITLE = {{LAPACK} Users' Guide},
      EDITION = {Third},
      PUBLISHER = {Society for Industrial and Applied Mathematics},
      YEAR = {1999},
      ADDRESS = {Philadelphia, PA},
      ISBN = {0-89871-447-8 (paperback)} }

*
* */

/* 
 * In lbfgsb.h, we can choose to redefine the standard
 * names, e.g, "daxpy", to "daxpyRef"
 * ("Ref" stands for "Reference implementation", since
 * these are the reference implementation from Netlib)
 * */
#include "lbfgsb.h"

/* Subroutine */ int daxpyRef(integer *n, double *da, double *dx, 
	integer *incx, double *dy, integer *incy)
{

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     constant times a vector plus a vector.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   



       Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
        return 0;
    }
    if (*da == 0.) {
        return 0;
    }
    if (*incx == 1 && *incy == 1) {
        goto L20;
    }

    /*        code for unequal increments or equal increments   
              not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
        DY(iy) += *da * DX(ix);
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    return 0;

    /*        code for both increments equal to 1   


              clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
        DY(i) += *da * DX(i);
        /* L30: */
    }
    if (*n < 4) {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 4) {
        DY(i) += *da * DX(i);
        DY(i + 1) += *da * DX(i + 1);
        DY(i + 2) += *da * DX(i + 2);
        DY(i + 3) += *da * DX(i + 3);
        /* L50: */
    }
    return 0;
} /* daxpyRef */


/* Subroutine */ int dcopyRef(integer *n, double *dx, integer *incx, 
        double *dy, integer *incy)
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


    /*     copies a vector, x, to a vector, y.   
           uses unrolled loops for increments equal to one.   
           jack dongarra, linpack, 3/11/78.   
           modified 12/3/93, array(1) declarations changed to array(*)   



           Parameter adjustments   
           Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
        return 0;
    }
    if (*incx == 1 && *incy == 1) {
        goto L20;
    }

    /*        code for unequal increments or equal increments   
              not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
        DY(iy) = DX(ix);
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    return 0;

    /*        code for both increments equal to 1   


              clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
        DY(i) = DX(i);
        /* L30: */
    }
    if (*n < 7) {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 7) {
        DY(i) = DX(i);
        DY(i + 1) = DX(i + 1);
        DY(i + 2) = DX(i + 2);
        DY(i + 3) = DX(i + 3);
        DY(i + 4) = DX(i + 4);
        DY(i + 5) = DX(i + 5);
        DY(i + 6) = DX(i + 6);
        /* L50: */
    }
    return 0;
} /* dcopyRef */

double ddotRef(integer *n, double *dx, integer *incx, double *dy, 
        integer *incy)
{


    /* System generated locals */
    integer i__1;
    double ret_val;

    /* Local variables */
    static integer i, m;
    static double dtemp;
    static integer ix, iy, mp1;


    /*     forms the dot product of two vectors.   
           uses unrolled loops for increments equal to one.   
           jack dongarra, linpack, 3/11/78.   
           modified 12/3/93, array(1) declarations changed to array(*)   



           Parameter adjustments   
           Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
        return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
        goto L20;
    }

    /*        code for unequal increments or equal increments   
              not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
        dtemp += DX(ix) * DY(iy);
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    ret_val = dtemp;
    return ret_val;

    /*        code for both increments equal to 1   


              clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
        dtemp += DX(i) * DY(i);
        /* L30: */
    }
    if (*n < 5) {
        goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 5) {
        dtemp = dtemp + DX(i) * DY(i) + DX(i + 1) * DY(i + 1) + DX(i + 2) * 
            DY(i + 2) + DX(i + 3) * DY(i + 3) + DX(i + 4) * DY(i + 4);
        /* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddotRef */

/* Subroutine */ int dscalRef(integer *n, double *da, double *dx, 
        integer *incx)
{


    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


    /*     scales a vector by a constant.   
           uses unrolled loops for increment equal to one.   
           jack dongarra, linpack, 3/11/78.   
           modified 3/93 to return if incx .le. 0.   
           modified 12/3/93, array(1) declarations changed to array(*)   



           Parameter adjustments   
           Function Body */
#define DX(I) dx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
        return 0;
    }
    if (*incx == 1) {
        goto L20;
    }

    /*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
        DX(i) = *da * DX(i);
        /* L10: */
    }
    return 0;

    /*        code for increment equal to 1   


              clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
        goto L40;
    }
    i__2 = m;
    for (i = 1; i <= m; ++i) {
        DX(i) = *da * DX(i);
        /* L30: */
    }
    if (*n < 5) {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= *n; i += 5) {
        DX(i) = *da * DX(i);
        DX(i + 1) = *da * DX(i + 1);
        DX(i + 2) = *da * DX(i + 2);
        DX(i + 3) = *da * DX(i + 3);
        DX(i + 4) = *da * DX(i + 4);
        /* L50: */
    }
    return 0;
} /* dscalRef */

