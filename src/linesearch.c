#include "lbfgsb.h"
static integer c__1 = 1;


 int lnsrlb(integer *n, double *l, double *u, 
	integer *nbd, double *x, double *f, double *fold, 
	double *gd, double *gdold, double *g, double *d__, 
	double *r__, double *t, double *z__, double *stp, 
	double *dnorm, double *dtd, double *xstep, double *
	stpmx, integer *iter, integer *ifun, integer *iback, integer *nfgv, 
	integer *info, integer *task, logical *boxed, logical *cnstnd, integer *
	csave, integer *isave, double *dsave) /* ftnlen task_len, ftnlen 
	csave_len) */
{
    /*
    ********** 

    Subroutine lnsrlb 

    This subroutine calls subroutine dcsrch from the Minpack2 library 
      to perform the line search.  Subroutine dscrch is safeguarded so 
      that all trial points lie within the feasible region. 

    Subprograms called: 

      Minpack2 Library ... dcsrch. 

      Linpack ... dtrsl, ddot. 


                          *  *  * 

    NEOS, November 1994. (Latest revision June 1996.) 
    Optimization Technology Center. 
    Argonne National Laboratory and Northwestern University. 
    Written by 
                       Ciyou Zhu 
    in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. 


    ********** 
    */


    /* Table of constant values */
    static double c_b14 = FTOL;
    static double c_b15 = GTOL;
    static double c_b16 = XTOL;
    static double c_b17 = STEPMIN;
    /* System generated locals */
    integer i__1;
    double d__1;


    /* Local variables */
    static integer i__;
    static double a1, a2;

    /* Parameter adjustments */
    --z__;
    --t;
    --r__;
    --d__;
    --g;
    --x;
    --nbd;
    --u;
    --l;
    --isave;
    --dsave;

    /* Function Body */
    if ( *task == FG_LN ) { 
        goto L556;
    }
    *dtd = ddot(n, &d__[1], &c__1, &d__[1], &c__1);
    *dnorm = sqrt(*dtd);
    /*     Determine the maximum step length. */
    *stpmx = 1e10;
    if (*cnstnd) {
        if (*iter == 0) {
            *stpmx = 1.;
        } else {
            i__1 = *n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                a1 = d__[i__];
                if (nbd[i__] != 0) {
                    if (a1 < 0. && nbd[i__] <= 2) {
                        a2 = l[i__] - x[i__];
                        if (a2 >= 0.) {
                            *stpmx = 0.;
                        } else if (a1 * *stpmx < a2) {
                            *stpmx = a2 / a1;
                        }
                    } else if (a1 > 0. && nbd[i__] >= 2) {
                        a2 = u[i__] - x[i__];
                        if (a2 <= 0.) {
                            *stpmx = 0.;
                        } else if (a1 * *stpmx > a2) {
                            *stpmx = a2 / a1;
                        }
                    }
                }
                /* L43: */
            }
        }
    }
    if (*iter == 0 && ! (*boxed)) {
        /* Computing MIN */
        d__1 = 1. / *dnorm;
        *stp = fmin(d__1,*stpmx);
    } else {
        *stp = 1.;
    }
    dcopy(n, &x[1], &c__1, &t[1], &c__1);
    dcopy(n, &g[1], &c__1, &r__[1], &c__1);
    *fold = *f;
    *ifun = 0;
    *iback = 0;
    *csave = START;
L556:
    *gd = ddot(n, &g[1], &c__1, &d__[1], &c__1);
    if (*ifun == 0) {
        *gdold = *gd;
        if (*gd >= 0.) {
            /*  the directional derivative >=0. */
            /*  Line search is impossible. */
            printf("ascend direction in projection gd = %.2e\n", *gd ); 
            *info = -4;
            return 0;
        }
    }
    dcsrch(f, gd, stp, &c_b14, &c_b15, &c_b16, &c_b17, stpmx, csave, &isave[
            1], &dsave[1]);/* (ftnlen)60);*/
    *xstep = *stp * *dnorm;
    if (  !(IS_WARNING(*csave)) && !(IS_CONVERGED(*csave)) )  {
/*     if (     !( (csave>=WARNING)&&(csave<=WARNING_END) )   &&   */
/*              !( (csave>=CONVERGENCE)&&(csave<=CONVERGENCE_END) )   ) { */
/*     if (s_cmp(csave, "CONV", (ftnlen)4, (ftnlen)4) != 0 && s_cmp(csave, "WARN" */
/*                 , (ftnlen)4, (ftnlen)4) != 0) { */
        *task = FG_LNSRCH;
        ++(*ifun);
        ++(*nfgv);
        *iback = *ifun - 1;
        if (*stp == 1.) {
            dcopy(n, &z__[1], &c__1, &x[1], &c__1);
        } else {
            i__1 = *n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                x[i__] = *stp * d__[i__] + t[i__];
                /* L41: */
            }
        }
    } else {
        *task = NEW_X;
    }
    return 0;
} /* lnsrlb */

/* ======================= The end of lnsrlb ============================= */
int dcsrch(double *f, double *g, double *stp, 
        double *ftol, double *gtol, double *xtol, double *
        stpmin, double *stpmax, integer *task, integer *isave, double *
        dsave) /* ftnlen task_len) */
{
    /* System generated locals */
    double d__1;


    /* Local variables */
    static double fm, gm, fx, fy, gx, gy, fxm, fym, gxm, gym, stx, sty;
    static integer stage;
    static double finit, ginit, width, ftest, gtest, stmin, stmax, width1;
    static logical brackt;

    /*
     ********** 

     Subroutine dcsrch 

     This subroutine finds a step that satisfies a sufficient 
     decrease condition and a curvature condition. 

     Each call of the subroutine updates an interval with 
     endpoints stx and sty. The interval is initially chosen 
     so that it contains a minimizer of the modified function 

     psi(stp) = f(stp) - f(0) - ftol*stp*f'(0). 

     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the 
     interval is chosen so that it contains a minimizer of f. 

     The algorithm is designed to find a step that satisfies 
     the sufficient decrease condition 

     f(stp) <= f(0) + ftol*stp*f'(0), 

     and the curvature condition 

     abs(f'(stp)) <= gtol*abs(f'(0)). 

     If ftol is less than gtol and if, for example, the function 
     is bounded below, then there is always a step which satisfies 
     both conditions. 

     If no step can be found that satisfies both conditions, then 
     the algorithm stops with a warning. In this case stp only 
     satisfies the sufficient decrease condition. 

     A typical invocation of dcsrch has the following outline: 

     task = 'START' 
     10 continue 
     call dcsrch( ... ) 
     if (task .eq. 'FG') then 
     Evaluate the function and the gradient at stp 
     goto 10 
     end if 

NOTE: The user must no alter work arrays between calls. 

The subroutine statement is 

subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax, 
task,isave,dsave) 
where 

f is a double precision variable. 
On initial entry f is the value of the function at 0. 
On subsequent entries f is the value of the 
function at stp. 
On exit f is the value of the function at stp. 

g is a double precision variable. 
On initial entry g is the derivative of the function at 0. 
On subsequent entries g is the derivative of the 
function at stp. 
On exit g is the derivative of the function at stp. 

stp is a double precision variable. 
On entry stp is the current estimate of a satisfactory 
step. On initial entry, a positive initial estimate 
must be provided. 
On exit stp is the current estimate of a satisfactory step 
if task = 'FG'. If task = 'CONV' then stp satisfies 
the sufficient decrease and curvature condition. 

    ftol is a double precision variable. 
        On entry ftol specifies a nonnegative tolerance for the 
        sufficient decrease condition. 
        On exit ftol is unchanged. 

        gtol is a double precision variable. 
        On entry gtol specifies a nonnegative tolerance for the 
        curvature condition. 
        On exit gtol is unchanged. 

        xtol is a double precision variable. 
        On entry xtol specifies a nonnegative relative tolerance 
        for an acceptable step. The subroutine exits with a 
            warning if the relative difference between sty and stx 
                is less than xtol. 
                On exit xtol is unchanged. 

                stpmin is a double precision variable. 
                On entry stpmin is a nonnegative lower bound for the step. 
                On exit stpmin is unchanged. 

                stpmax is a double precision variable. 
                On entry stpmax is a nonnegative upper bound for the step. 
                On exit stpmax is unchanged. 

                task is a character variable of length at least 60. 
                On initial entry task must be set to 'START'. 
                On exit task indicates the required action: 

                If task(1:2) = 'FG' then evaluate the function and 
                derivative at stp and call dcsrch again. 

                If task(1:4) = 'CONV' then the search is successful. 

                If task(1:4) = 'WARN' then the subroutine is not able 
                to satisfy the convergence conditions. The exit value of 
                stp contains the best point found during the search. 

                If task(1:5) = 'ERROR' then there is an error in the 
                input arguments. 

                On exit with convergence, a warning or an error, the 
                variable task contains additional information. 

                isave is an integer work array of dimension 2. 

                dsave is a double precision work array of dimension 13. 

                Subprograms called 

                MINPACK-2 ... dcstep 

                MINPACK-1 Project. June 1983. 
                Argonne National Laboratory. 
                Jorge J. More' and David J. Thuente. 

                MINPACK-2 Project. October 1993. 
                Argonne National Laboratory and University of Minnesota. 
                Brett M. Averick, Richard G. Carter, and Jorge J. More'. 

                ********** 
                */
                /*     Initialization block. */
                /* Parameter adjustments */
                --dsave;
    --isave;

    /* Function Body */
    if ( *task == START ) {
        /*        Check the input arguments for errors.  See lbfgsb.h for messages */
        if (*stp < *stpmin)  *task=ERROR_SMALLSTP;
        if (*stp > *stpmax)  *task=ERROR_LARGESTP;
        if (*g >= 0.)        *task=ERROR_INITIAL;
        if (*ftol < 0.)      *task=ERROR_FTOL;
        if (*gtol < 0.)      *task=ERROR_GTOL;
        if (*xtol < 0.)      *task=ERROR_XTOL;
        if (*stpmin < 0.)    *task=ERROR_STP0;
        if (*stpmax < *stpmin) *task=ERROR_STP1;
        /*        Exit if there are errors on input. */
        if ( IS_ERROR(*task) ) {
            return 0;
        }
        /*        Initialize local variables. */
        brackt = FALSE_;
        stage = 1;
        finit = *f;
        ginit = *g;
        gtest = *ftol * ginit;
        width = *stpmax - *stpmin;
        width1 = width / .5;
        /*        The variables stx, fx, gx contain the values of the step, */
        /*        function, and derivative at the best step. */
        /*        The variables sty, fy, gy contain the value of the step, */
        /*        function, and derivative at sty. */
        /*        The variables stp, f, g contain the values of the step, */
        /*        function, and derivative at stp. */
        stx = 0.;
        fx = finit;
        gx = ginit;
        sty = 0.;
        fy = finit;
        gy = ginit;
        stmin = 0.;
        stmax = *stp + *stp * 4.;
/*         s_copy(task, "FG", task_len, (ftnlen)2); */
        *task = FG;
        goto L1000;
    } else {
        /*        Restore local variables. */
        if (isave[1] == 1) {
            brackt = TRUE_;
        } else {
            brackt = FALSE_;
        }
        stage = isave[2];
        ginit = dsave[1];
        gtest = dsave[2];
        gx = dsave[3];
        gy = dsave[4];
        finit = dsave[5];
        fx = dsave[6];
        fy = dsave[7];
        stx = dsave[8];
        sty = dsave[9];
        stmin = dsave[10];
        stmax = dsave[11];
        width = dsave[12];
        width1 = dsave[13];
    }
    /*     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
    /*     algorithm enters the second stage. */
    ftest = finit + *stp * gtest;
    if (stage == 1 && *f <= ftest && *g >= 0.) {
        stage = 2;
    }
    /*     Test for warnings. */
    if (brackt && (*stp <= stmin || *stp >= stmax)) {
        *task = WARNING_ROUND;
    }
    if (brackt && stmax - stmin <= *xtol * stmax) {
        *task = WARNING_XTOL;
    }
    if (*stp == *stpmax && *f <= ftest && *g <= gtest) {
        *task = WARNING_STPMAX;
    }
    if (*stp == *stpmin && (*f > ftest || *g >= gtest)) {
        *task = WARNING_STPMIN;
    }
    /*     Test for convergence. */
    if (*f <= ftest && fabs(*g) <= *gtol * (-ginit)) {
        *task = CONVERGENCE;
    }
    /*     Test for termination. */
/*     if (     ( (task>=WARNING)&&(task<=WARNING_END) )    || */
/*              ( (task>=CONVERGENCE)&&(task<=CONVERGENCE_END) )   ) { */
    if ( (IS_WARNING(*task)) || (IS_CONVERGED(*task) ) ) {
        goto L1000;
    }
    /*     A modified function is used to predict the step during the */
    /*     first stage if a lower function value has been obtained but */
    /*     the decrease is not sufficient. */
    if (stage == 1 && *f <= fx && *f > ftest) {
        /*        Define the modified function and derivative values. */
        fm = *f - *stp * gtest;
        fxm = fx - stx * gtest;
        fym = fy - sty * gtest;
        gm = *g - gtest;
        gxm = gx - gtest;
        gym = gy - gtest;
        /*        Call dcstep to update stx, sty, and to compute the new step. */
        dcstep(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &
                stmin, &stmax);
        /*        Reset the function and derivative values for f. */
        fx = fxm + stx * gtest;
        fy = fym + sty * gtest;
        gx = gxm + gtest;
        gy = gym + gtest;
    } else {
        /*       Call dcstep to update stx, sty, and to compute the new step. */
        dcstep(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &
                stmax);
    }
    /*     Decide if a bisection step is needed. */
    if (brackt) {
        if ((d__1 = sty - stx, fabs(d__1)) >= width1 * .66) {
            *stp = stx + (sty - stx) * .5;
        }
        width1 = width;
        width = (d__1 = sty - stx, fabs(d__1));
    }
    /*     Set the minimum and maximum steps allowed for stp. */
    if (brackt) {
        stmin = fmin(stx,sty);
        stmax = fmax(stx,sty);
    } else {
        stmin = *stp + (*stp - stx) * 1.1;
        stmax = *stp + (*stp - stx) * 4.;
    }
    /*     Force the step to be within the bounds stpmax and stpmin. */
    *stp = fmax(*stp,*stpmin);
    *stp = fmin(*stp,*stpmax);
    /*     If further progress is not possible, let stp be the best */
    /*     point obtained during the search. */
    /*     if (brackt && (*stp <= stmin || *stp >= stmax) || brackt && stmax - stmin  */
    /* SRB: guess the precedence. && precedence over || */
    if ( (brackt && (*stp <= stmin || *stp >= stmax) ) || (brackt && stmax - stmin 
                <= *xtol * stmax)) {
        *stp = stx;
    }
    /*     Obtain another function and derivative. */
/*     s_copy(task, "FG", task_len, (ftnlen)2); */
    *task = FG;
L1000:
    /*     Save local variables. */
    if (brackt) {
        isave[1] = 1;
    } else {
        isave[1] = 0;
    }
    isave[2] = stage;
    dsave[1] = ginit;
    dsave[2] = gtest;
    dsave[3] = gx;
    dsave[4] = gy;
    dsave[5] = finit;
    dsave[6] = fx;
    dsave[7] = fy;
    dsave[8] = stx;
    dsave[9] = sty;
    dsave[10] = stmin;
    dsave[11] = stmax;
    dsave[12] = width;
    dsave[13] = width1;
    return 0;
} /* dcsrch */

/* ====================== The end of dcsrch ============================== */
/* Subroutine */ int dcstep(double *stx, double *fx, double *dx, 
        double *sty, double *fy, double *dy, double *stp, 
        double *fp, double *dp, logical *brackt, double *stpmin, 
        double *stpmax)
{
    /* System generated locals */
    double d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static double p, q, r__, s, sgnd, stpc, stpf, stpq, gamma, theta;

    /*
     ********** 

     Subroutine dcstep 

     This subroutine computes a safeguarded step for a search 
     procedure and updates an interval that contains a step that 
     satisfies a sufficient decrease and a curvature condition. 

     The parameter stx contains the step with the least function 
     value. If brackt is set to .true. then a minimizer has 
     been bracketed in an interval with endpoints stx and sty. 
     The parameter stp contains the current step. 
     The subroutine assumes that if brackt is set to .true. then 

     min(stx,sty) < stp < fmax(stx,sty), 

     and that the derivative at stx is negative in the direction 
     of the step. 

     The subroutine statement is 

     subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt, 
     stpmin,stpmax) 

     where 

     stx is a double precision variable. 
     On entry stx is the best step obtained so far and is an 
     endpoint of the interval that contains the minimizer. 
     On exit stx is the updated best step. 

     fx is a double precision variable. 
     On entry fx is the function at stx. 
     On exit fx is the function at stx. 

     dx is a double precision variable. 
     On entry dx is the derivative of the function at 
     stx. The derivative must be negative in the direction of 
     the step, that is, dx and stp - stx must have opposite 
     signs. 
     On exit dx is the derivative of the function at stx. 

     sty is a double precision variable. 
     On entry sty is the second endpoint of the interval that 
     contains the minimizer. 
     On exit sty is the updated endpoint of the interval that 
     contains the minimizer. 

     fy is a double precision variable. 
     On entry fy is the function at sty. 
     On exit fy is the function at sty. 

     dy is a double precision variable. 
     On entry dy is the derivative of the function at sty. 
     On exit dy is the derivative of the function at the exit sty. 

     stp is a double precision variable. 
     On entry stp is the current step. If brackt is set to .true. 
     then on input stp must be between stx and sty. 
     On exit stp is a new trial step. 

     fp is a double precision variable. 
     On entry fp is the function at stp 
     On exit fp is unchanged. 

     dp is a double precision variable. 
     On entry dp is the the derivative of the function at stp. 
     On exit dp is unchanged. 

     brackt is an logical variable. 
    On entry brackt specifies if a minimizer has been bracketed. 
        Initially brackt must be set to .false. 
        On exit brackt specifies if a minimizer has been bracketed. 
        When a minimizer is bracketed brackt is set to .true. 

        stpmin is a double precision variable. 
        On entry stpmin is a lower bound for the step. 
        On exit stpmin is unchanged. 

        stpmax is a double precision variable. 
        On entry stpmax is an upper bound for the step. 
        On exit stpmax is unchanged. 

        MINPACK-1 Project. June 1983 
        Argonne National Laboratory. 
        Jorge J. More' and David J. Thuente. 

        MINPACK-2 Project. October 1993. 
        Argonne National Laboratory and University of Minnesota. 
        Brett M. Averick and Jorge J. More'. 

        ********** 
        */
        sgnd = *dp * (*dx / fabs(*dx));
    /*     First case: A higher function value. The minimum is bracketed. */
    /*     If the cubic step is closer to stx than the quadratic step, the */
    /*     cubic step is taken, otherwise the average of the cubic and */
    /*     quadratic steps is taken. */
    if (*fp > *fx) {
        theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
        /* Computing MAX */
        d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = fmax(d__1,d__2), 
          d__2 = fabs(*dp);
        s = fmax(d__1,d__2);
        /* Computing 2nd power */
        d__1 = theta / s;
        gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
        if (*stp < *stx) {
            gamma = -gamma;
        }
        p = gamma - *dx + theta;
        q = gamma - *dx + gamma + *dp;
        r__ = p / q;
        stpc = *stx + r__ * (*stp - *stx);
        stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2. * (*stp 
                - *stx);
        if ((d__1 = stpc - *stx, fabs(d__1)) < (d__2 = stpq - *stx, fabs(d__2)))
        {
            stpf = stpc;
        } else {
            stpf = stpc + (stpq - stpc) / 2.;
        }
        *brackt = TRUE_;
        /*     Second case: A lower function value and derivatives of opposite */
        /*     sign. The minimum is bracketed. If the cubic step is farther from */
        /*     stp than the secant step, the cubic step is taken, otherwise the */
        /*     secant step is taken. */
    } else if (sgnd < 0.) {
        theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
        /* Computing MAX */
        d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = fmax(d__1,d__2), 
          d__2 = fabs(*dp);
        s = fmax(d__1,d__2);
        /* Computing 2nd power */
        d__1 = theta / s;
        gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
        if (*stp > *stx) {
            gamma = -gamma;
        }
        p = gamma - *dp + theta;
        q = gamma - *dp + gamma + *dx;
        r__ = p / q;
        stpc = *stp + r__ * (*stx - *stp);
        stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
        if ((d__1 = stpc - *stp, fabs(d__1)) > (d__2 = stpq - *stp, fabs(d__2)))
        {
            stpf = stpc;
        } else {
            stpf = stpq;
        }
        *brackt = TRUE_;
        /*     Third case: A lower function value, derivatives of the same sign, */
        /*     and the magnitude of the derivative decreases. */
    } else if (fabs(*dp) < fabs(*dx)) {
        /*        The cubic step is computed only if the cubic tends to infinity */
        /*        in the direction of the step or if the minimum of the cubic */
        /*        is beyond stp. Otherwise the cubic step is defined to be the */
        /*        secant step. */
        theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
        /* Computing MAX */
        d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = fmax(d__1,d__2), 
          d__2 = fabs(*dp);
        s = fmax(d__1,d__2);
        /*        The case gamma = 0 only arises if the cubic does not tend */
        /*        to infinity in the direction of the step. */
        /* Computing MAX */
        /* Computing 2nd power */
        d__3 = theta / s;
        d__1 = 0., d__2 = d__3 * d__3 - *dx / s * (*dp / s);
        gamma = s * sqrt((fmax(d__1,d__2)));
        if (*stp > *stx) {
            gamma = -gamma;
        }
        p = gamma - *dp + theta;
        q = gamma + (*dx - *dp) + gamma;
        r__ = p / q;
        if (r__ < 0. && gamma != 0.) {
            stpc = *stp + r__ * (*stx - *stp);
        } else if (*stp > *stx) {
            stpc = *stpmax;
        } else {
            stpc = *stpmin;
        }
        stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
        if (*brackt) {
            /*           A minimizer has been bracketed. If the cubic step is */
            /*           closer to stp than the secant step, the cubic step is */
            /*           taken, otherwise the secant step is taken. */
            if ((d__1 = stpc - *stp, fabs(d__1)) < 
                (d__2 = stpq - *stp, fabs(d__2))) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            if (*stp > *stx) {
                /* Computing MIN */
                d__1 = *stp + (*sty - *stp) * .66;
                stpf = fmin(d__1,stpf);
            } else {
                /* Computing MAX */
                d__1 = *stp + (*sty - *stp) * .66;
                stpf = fmax(d__1,stpf);
            }
        } else {
            /*           A minimizer has not been bracketed. If the cubic step is */
            /*           farther from stp than the secant step, the cubic step is */
            /*           taken, otherwise the secant step is taken. */
            if ((d__1 = stpc - *stp, fabs(d__1)) > 
                (d__2 = stpq - *stp, fabs(d__2))) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            stpf = fmin(*stpmax,stpf);
            stpf = fmax(*stpmin,stpf);
        }
        /*     Fourth case: A lower function value, derivatives of the same sign, */
        /*     and the magnitude of the derivative does not decrease. If the */
        /*     minimum is not bracketed, the step is either stpmin or stpmax, */
        /*     otherwise the cubic step is taken. */
    } else {
        if (*brackt) {
            theta = (*fp - *fy) * 3. / (*sty - *stp) + *dy + *dp;
            /* Computing MAX */
            d__1 = fabs(theta), d__2 = fabs(*dy), d__1 = fmax(d__1,d__2), 
              d__2 = fabs(*dp);
            s = fmax(d__1,d__2);
            /* Computing 2nd power */
            d__1 = theta / s;
            gamma = s * sqrt(d__1 * d__1 - *dy / s * (*dp / s));
            if (*stp > *sty) {
                gamma = -gamma;
            }
            p = gamma - *dp + theta;
            q = gamma - *dp + gamma + *dy;
            r__ = p / q;
            stpc = *stp + r__ * (*sty - *stp);
            stpf = stpc;
        } else if (*stp > *stx) {
            stpf = *stpmax;
        } else {
            stpf = *stpmin;
        }
    }
    /*     Update the interval which contains a minimizer. */
    if (*fp > *fx) {
        *sty = *stp;
        *fy = *fp;
        *dy = *dp;
    } else {
        if (sgnd < 0.) {
            *sty = *stx;
            *fy = *fx;
            *dy = *dx;
        }
        *stx = *stp;
        *fx = *fp;
        *dx = *dp;
    }
    /*     Compute the new step. */
    *stp = stpf;
    return 0;
} /* dcstep */

