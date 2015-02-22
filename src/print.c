#include "lbfgsb.h"

int prn1lb(integer *n, integer *m, double *l, 
	double *u, double *x, integer *iprint, fileType itfile, 
	double *epsmch)
{
    /*
    ************ 

    Subroutine prn1lb 

    This subroutine prints the input data, initial point, upper and 
      lower bounds of each variable, machine precision, as well as 
      the headings of the output. 


                          *  *  * 

    NEOS, November 1994. (Latest revision June 1996.) 
    Optimization Technology Center. 
    Argonne National Laboratory and Northwestern University. 
    Written by 
                       Ciyou Zhu 
    in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. 


        ************ 
    */
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --x;
    --u;
    --l;

    /* Function Body */
    if (*iprint >= 0) {
        printf("           * * *\n");
        printf("        RUNNING THE L-BFGS-B CODE\n");
        printf("           * * *\n");
        printf("Machine precision = %.2e\n", *epsmch );
        printf(" N = %10ld\n M = %10ld\n", *n, *m );
        if (*iprint >= 1) {
            if (*iprint > 100) {
                printf("L  =");
                i__1 = *n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    printf("%.2e ", l[i__] );
                }
                printf("\n");
                printf("X0 =");
                i__1 = *n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    printf("%.2e ", x[i__] );
                }
                printf("\n");
                printf("U  =");
                i__1 = *n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    printf("%.2e ", u[i__] );
                }
                printf("\n");
            }
        }
    }
    return 0;
} /* prn1lb */


/* Subroutine */ int prn2lb(integer *n, double *x, double *f, 
	double *g, integer *iprint, fileType itfile, integer *iter, 
	integer *nfgv, integer *nact, double *sbgnrm, integer *nseg, integer  
	*word, integer *iword, integer *iback, double *stp, double *
	xstep, ftnlen word_len)
{
    /*
    ************ 

    Subroutine prn2lb 

    This subroutine prints out new information after a successful 
      line search. 


                          *  *  * 

    NEOS, November 1994. (Latest revision June 1996.) 
    Optimization Technology Center. 
    Argonne National Laboratory and Northwestern University. 
    Written by 
                       Ciyou Zhu 
    in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. 


    ************ 
    */
/*           'word' records the status of subspace solutions. */
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, imod;

    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    if (*iword == 0) {
        *word = WORD_CON;
        /*  the subspace minimization converged. */
/*         s_copy(word, "con", (ftnlen)3, (ftnlen)3); */
    } else if (*iword == 1) {
        *word = WORD_BND;
        /*  the subspace minimization stopped at a bound. */
/*         s_copy(word, "bnd", (ftnlen)3, (ftnlen)3); */
    } else if (*iword == 5) {
        *word = WORD_TNT;
        /*  the truncated Newton step has been used. */
/*         s_copy(word, "TNT", (ftnlen)3, (ftnlen)3); */
    } else {
        *word = WORD_DEFAULT;
/*         s_copy(word, "---", (ftnlen)3, (ftnlen)3); */
    }
    if (*iprint >= 99) {
        printf("LINE SEARCH %ld times; norm of step = %.2e\n", *iback, *xstep );
        printf("At iterate %5ld, f(x)= %5.2e, ||proj grad||_infty = %.2e\n",*iter,*f,*sbgnrm );
        if (*iprint > 100) {
            printf("X =");
            i__1 = *n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                printf("%.2e ", x[i__] );
            }
            printf("\n");
            printf("G =");
            i__1 = *n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                printf("%.2e ", g[i__] );
            }
            printf("\n");
        }
    } else if (*iprint > 0) {
        imod = *iter % *iprint;
        if (imod == 0) {
            printf("At iterate %5ld, f(x)= %5.2e, ||proj grad||_infty = %.2e\n",*iter,*f,*sbgnrm );
        }
    }
    return 0;
} /* prn2lb */

int prn3lb(integer *n, double *x, double *f, integer *
	task, integer *iprint, integer *info, fileType itfile, integer *iter, 
	integer *nfgv, integer *nintol, integer *nskip, integer *nact, 
	double *sbgnrm, double *time, integer *nseg, integer *word, 
	integer *iback, double *stp, double *xstep, integer *k, 
	double *cachyt, double *sbtime, double *lnscht, ftnlen 
	task_len, ftnlen word_len)
{
    /*
    ************ 

    Subroutine prn3lb 

    This subroutine prints out information when either a built-in 
      convergence test is satisfied or when an error message is 
      generated. 


                          *  *  * 

    NEOS, November 1994. (Latest revision June 1996.) 
    Optimization Technology Center. 
    Argonne National Laboratory and Northwestern University. 
    Written by 
                       Ciyou Zhu 
    in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. 


    ************ 
    */
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    /* Parameter adjustments */
    --x;

    /* Function Body */
/*     if (s_cmp(task, "ERROR", (ftnlen)5, (ftnlen)5) == 0) { */
    if ( IS_ERROR(*task) ) {
        goto L999;
    }
    if (*iprint >= 0) {


        printf("           * * * \n");
        printf("Tit   = total number of iterations\n");
        printf("Tnf   = total number of function evaluations\n");
        printf("Tnint = total number of segments explored during Cauchy searches\n");
        printf("Skip  = number of BFGS updates skipped\n");
        printf("Nact  = number of active bounds at final generalized Cauchy point\n");
        printf("Projg = norm of the final projected gradient\n");
        printf("F     = final function value\n");
        printf("           * * * \n");

        printf("   N    Tit   Tnf  Tnint  Skip  Nact      Projg        F\n");
        printf("%5ld %5ld %5ld %5ld %5ld %5ld\t%6.2e %9.5e\n", 
                *n, *iter, *nfgv, *nintol, *nskip, *nact, *sbgnrm, *f );
        if (*iprint >= 100) {
            printf("X = ");
            i__1 = *n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                printf(" %.2e", x[i__] );
            }
            printf("\n");
        }
        if (*iprint >= 1) {
            printf("F(x) = %.9e\n", *f );
        }
    }
L999:
    if (*iprint >= 0) {
        printf("%ld\n",*task);
        if (*info != 0) {
            if (*info == -1) {
                printf(" Matrix in 1st Cholesky factorization in formk is not Pos. Def.\n");
            }
            if (*info == -2) {
                printf(" Matrix in 2nd Cholesky factorization in formk is not Pos. Def.\n");
            }
            if (*info == -3) {
                printf(" Matrix in the Cholesky factorization in formt is not Pos. Def.\n");
            }
            if (*info == -4) {
                printf(" Derivative >= 0, backtracking line search impossible.\n");
                printf("  Previous x, f and g restored.\n");
                printf(" Possible causes: 1 error in function or gradient evaluation;\n");
                printf("                  2 rounding errors dominate computation.\n");
            }
            if (*info == -5) {
                printf(" Warning:  more than 10 function and gradient\n");
                printf("   evaluations in the last line search.  Termination\n");
                printf("   may possibly be caused by a bad search direction.\n");
            }
            if (*info == -6) {
                printf(" Input nbd(%ld) is invalid\n", *k );
            }
            if (*info == -7) {
                printf(" l(%ld) > u(%ld). No feasible solution.\n", *k, *k );
            }
            if (*info == -8) {
                printf(" The triangular system is singular.\n");
            }
            if (*info == -9) {
                printf(" Line search cannot locate an adequate point after 20 function\n");
                printf("  and gradient evaluations.  Previous x, f and g restored.\n");
                printf(" Possible causes: 1 error in function or gradient evaluation;\n");
                printf("                  2 rounding error dominate computation.\n");
            }
        }

        if (*iprint >= 1) {
            printf("Cauchy                time %.3e seconds.\n", *cachyt );
            printf("Subspace minimization time %.3e seconds.\n", *sbtime );
            printf("Line search           time %.3e seconds.\n", *lnscht );
        }
        printf(" Total User time %.3e seconds.\n", *time );
    }
    return 0;
} /* prn3lb */


int errclb(integer *n, integer *m, double *factr, 
	double *l, double *u, integer *nbd, integer *task, integer *info,
	 integer *k, ftnlen task_len)
{
    /*
    ************ 

    Subroutine errclb 

    This subroutine checks the validity of the input data. 


                          *  *  * 

    NEOS, November 1994. (Latest revision June 1996.) 
    Optimization Technology Center. 
    Argonne National Laboratory and Northwestern University. 
    Written by 
                       Ciyou Zhu 
    in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. 


    ************ 
    */
/*     Check the input arguments for errors. */
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --nbd;
    --u;
    --l;

    /* Function Body */
    if (*n <= 0) *task = ERROR_N0;
    if (*m <= 0) *task = ERROR_M0;
    if (*factr < 0.)  *task = ERROR_FACTR;
    /*     Check the validity of the arrays nbd(i), u(i), and l(i). */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (nbd[i__] < 0 || nbd[i__] > 3) {
            /*   return */
            *task = ERROR_NBD;
            *info = -6;
            *k = i__;
        }
        if (nbd[i__] == 2) {
            if (l[i__] > u[i__]) {
                /* return */
                *task = ERROR_FEAS;
                *info = -7;
                *k = i__;
            }
        }
        /* L10: */
    }
    return 0;
} /* errclb */

