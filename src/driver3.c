#include "lbfgsb.h"


/*
                            DRIVER 3 in Fortran 77 
    -------------------------------------------------------------- 
           TIME-CONTROLLED DRIVER FOR L-BFGS-B (version 3.0) 
    -------------------------------------------------------------- 

       L-BFGS-B is a code for solving large nonlinear optimization 
            problems with simple bounds on the variables. 

       The code can also be used for unconstrained problems and is 
       as efficient for these problems as the earlier limited memory 
                         code L-BFGS. 

       This driver shows how to terminate a run after some prescribed 
       CPU time has elapsed, and how to print the desired information 
       before exiting. 

    References: 

       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited 
       memory algorithm for bound constrained optimization'', 
       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. 

       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN 
       Subroutines for Large Scale Bound Constrained Optimization'' 
       Tech. Report, NAM-11, EECS Department, Northwestern University, 
       1994. 


         (Postscript files of these papers are available via anonymous 
          ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) 

                             *  *  * 

        February 2011   (latest revision) 
        Optimization Center at Northwestern University 
        Instituto Tecnologico Autonomo de Mexico 

        Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778: 
        L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained 
        Optimization"  (2011). To appear in  ACM Transactions on 
        Mathematical Software, 
        */


/*     ************** */

int main(void)
{

    /* System generated locals */
    integer i__1;
    double d__1, d__2;

    /* Local variables */
    static double f, g[1024];
    static integer i__, j;
    static double l[1024];
    static integer m, n;
    static double u[1024], x[1024], t1, t2, wa[43251];
    static integer nbd[1024], iwa[3072];
/*     static char task[60]; */
    static integer taskValue;
    static integer *task=&taskValue; /* must initialize !! */
/*      http://stackoverflow.com/a/11278093/269192 */
    static double time1, time2, factr;
/*     static char csave[60]; */
    static integer csaveValue;
    static integer *csave=&csaveValue;
    static double dsave[29];
    static integer isave[44];
    static logical lsave[4];
    static double pgtol;
    static double tlimit;
    static integer iprint;


    /*
    This time-controlled driver shows that it is possible to terminate 
    a run by elapsed CPU time, and yet be able to print all desired 
    information. This driver also illustrates the use of two 
    stopping criteria that may be used in conjunction with a limit 
    on execution time. The sample problem used here is the same as in 
    driver1 and driver2 (the extended Rosenbrock function with bounds 
    on the variables). 
    nmax is the dimension of the largest problem to be solved. 
    mmax is the maximum number of limited memory corrections. 
    Declare the variables needed by the code. 
    A description of all these variables is given at the end of 
    driver1. 
    Declare a few additional variables for the sample problem 
    and for keeping track of time. 
    */

    /*     We specify a limite on the CPU time (in seconds). */
    tlimit = .2f;
    /*     We suppress the default output.  (The user could also elect to */
    /*       use the default output by choosing iprint >= 0.) */
    iprint = -1;
    /*     We suppress the code-supplied stopping tests because we will */
    /*       provide our own termination conditions */
    factr = 0.;
    pgtol = 0.;
    /*     We specify the dimension n of the sample problem and the number */
    /*        m of limited memory corrections stored.  (n and m should not */
    /*        exceed the limits nmax and mmax respectively.) */
    n = 1000;
    m = 10;
    /*     We now specify nbd which defines the bounds on the variables: */
    /*                    l   specifies the lower bounds, */
    /*                    u   specifies the upper bounds. */
    /*     First set bounds on the odd-numbered variables. */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
        nbd[i__ - 1] = 2;
        l[i__ - 1] = 1.;
        u[i__ - 1] = 100.;
        /* L10: */
    }
    /*     Next set bounds on the even-numbered variables. */
    i__1 = n;
    for (i__ = 2; i__ <= i__1; i__ += 2) {
        nbd[i__ - 1] = 2;
        l[i__ - 1] = -100.;
        u[i__ - 1] = 100.;
        /* L12: */
    }
    /*     We now define the starting point. */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        x[i__ - 1] = 3.;
        /* L14: */
    }
    /*     We now write the heading of the output. */
    printf("     Solving sample problem.\n");
    printf("      (f = 0.0 at the optimal solution.)\n");
    /*     We start the iteration by initializing task. */

    *task = START;
    /*        ------- the beginning of the loop ---------- */
    /*     We begin counting the CPU time. */
    timer(&time1);
L111:
    /*     This is the call to the L-BFGS-B code. */
    setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
            iprint, csave, lsave, isave, dsave);
    if (IS_FG(*task) ) {
        /*        the minimization routine has returned to request the */
        /*        function f and gradient g values at the current x. */
        /*        Before evaluating f and g we check the CPU time spent. */
        timer(&time2);
        if (time2 - time1 > tlimit) {
            *task = STOP_CPU;
/*             s_copy(task, "STOP: CPU EXCEEDING THE TIME LIMIT.", (ftnlen)60, ( */
/*                         ftnlen)35); */
            /*          Note: Assigning task(1:4)='STOP' will terminate the run; */
            /*          setting task(7:9)='CPU' will restore the information at */
            /*          the latest iterate generated by the code so that it can */
            /*          be correctly printed by the driver. */
            /*          In this driver we have chosen to disable the */
            /*          printing options of the code (we set iprint=-1); */
            /*          instead we are using customized output: we print the */
            /*          latest value of x, the corresponding function value f and */
            /*          the norm of the projected gradient |proj g|. */
            /*          We print out the information contained in task. */

            /*          We print the latest iterate contained in wa(j+1:j+n), where */

            /* Computing 2nd power */
            j = n * 3 + (m << 1) * n + m * m * 11;
            printf("Latest iterate X=");
            for (i__ = j + 1; i__ <= m ; ++i__) {
                printf("%.2e ", wa[i__ -1] );
            }
            printf("\n");
            /*          We print the function value f and the norm of the projected */
            /*          gradient |proj g| at the last iterate; they are stored in */
            /*          dsave(2) and dsave(13) respectively. */
            printf("At latest iterate   f = %.2e, |proj g| = %.2e\n", dsave[1], dsave[12] );
        } else {
            /*          The time limit has not been reached and we compute */
            /*          the function value f for the sample problem. */
            /* Computing 2nd power */
            d__1 = x[0] - 1.;
            f = d__1 * d__1 * .25;
            i__1 = n;
            for (i__ = 2; i__ <= i__1; ++i__) {
                /* Computing 2nd power */
                d__2 = x[i__ - 2];
                /* Computing 2nd power */
                d__1 = x[i__ - 1] - d__2 * d__2;
                f += d__1 * d__1;
                /* L20: */
            }
            f *= 4.;
            /*          Compute gradient g for the sample problem. */
            /* Computing 2nd power */
            d__1 = x[0];
            t1 = x[1] - d__1 * d__1;
            g[0] = (x[0] - 1.) * 2. - x[0] * 16. * t1;
            i__1 = n - 1;
            for (i__ = 2; i__ <= i__1; ++i__) {
                t2 = t1;
                /* Computing 2nd power */
                d__1 = x[i__ - 1];
                t1 = x[i__] - d__1 * d__1;
                g[i__ - 1] = t2 * 8. - x[i__ - 1] * 16. * t1;
                /* L22: */
            }
            g[n - 1] = t1 * 8.;
        }
        /*          go back to the minimization routine. */
        goto L111;
    }

    if ( *task == NEW_X ) {
        /*        the minimization routine has returned with a new iterate. */
        /*        The time limit has not been reached, and we test whether */
        /*        the following two stopping tests are satisfied: */
        /*        1) Terminate if the total number of f and g evaluations */
        /*             exceeds 900. */
        if (isave[33] >= 900) {
            *task = STOP_ITER;
/*             s_copy(task, "STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIM" */
/*                     "IT", (ftnlen)60, (ftnlen)52); */
        }
        /*        2) Terminate if  |proj g|/(1+|f|) < 1.0d-10. */
        if (dsave[12] <= (fabs(f) + 1.) * 1e-10) {
            *task = STOP_GRAD;
/*             s_copy(task, "STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL", */
/*                     (ftnlen)60, (ftnlen)50); */
        }
        /*        We wish to print the following information at each iteration: */
        /*          1) the current iteration number, isave(30), */
        /*          2) the total number of f and g evaluations, isave(34), */
        /*          3) the value of the objective function f, */
        /*          4) the norm of the projected gradient,  dsve(13) */

        /*        See the comments at the end of driver1 for a description */
        /*        of the variables isave and dsave. */
        printf("Iterate %ld, nfg = %ld, f = %.2e, |proj g| = %.2e\n", isave[29], isave[33], f, dsave[12] );
        /*        If the run is to be terminated, we print also the information */
        /*        contained in task as well as the final value of x. */
/*         if (s_cmp(task, "STOP", (ftnlen)4, (ftnlen)4) == 0) { */
        if (IS_STOP(*task)) {
            i__1 = n;
            printf("Final X = ");
            if ( n < 30 ) {
                for (i__ = 1; i__ <= i__1; ++i__) {
                    printf("%.2e ", x[i__ -1] );
                }
                printf("\n");
            } else {
                /* Print out just some */
                for (i__ = 1; i__ <= 10 ; ++i__) {
                    printf("%.2e ", x[i__ -1] );
                }
                printf("... \n ...(suppressing some output)...\n ... ");
                for (i__ = n-9; i__ <= n ; ++i__) {
                    printf("%.2e ", x[i__ -1] );
                }
                printf("\n");
            }
        }
        /*          go back to the minimization routine. */
        goto L111;
    }
    /*           ---------- the end of the loop ------------- */
    /*     If task is neither FG nor NEW_X we terminate execution. */
    return 0;
} /* MAIN__ */

