/* Regression test for concurrent, independent setulb calls. */
#define _XOPEN_SOURCE 700

#include "lbfgsb.h"

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>

#define N 64
#define M 5
#define THREADS 8
#define ROUNDS 16
#define WA_SIZE (2 * M * N + 11 * M * M + 5 * N + 8 * M)

typedef struct {
    double x[N];
    double f;
    integer task;
} result;

typedef struct {
    pthread_barrier_t *start;
    result *results;
    int index;
} worker_args;

static void rosenbrock(const double *x, double *f, double *g)
{
    int i;
    double t1;

    *f = .25 * (x[0] - 1.) * (x[0] - 1.);
    for (i = 1; i < N; ++i) {
        double residual = x[i] - x[i - 1] * x[i - 1];
        *f += residual * residual;
    }
    *f *= 4.;

    t1 = x[1] - x[0] * x[0];
    g[0] = 2. * (x[0] - 1.) - 16. * x[0] * t1;
    for (i = 1; i < N - 1; ++i) {
        double t2 = t1;
        t1 = x[i + 1] - x[i] * x[i];
        g[i] = 8. * t2 - 16. * x[i] * t1;
    }
    g[N - 1] = 8. * t1;
}

static result optimize(void)
{
    integer n = N, m = M, nbd[N], iwa[3 * N];
    integer task = START, iprint = -1, csave = 0, isave[44];
    logical lsave[4];
    double x[N], l[N], u[N], g[N], wa[WA_SIZE], dsave[29];
    double f = 0., factr = 1e7, pgtol = 1e-8;
    int i;
    result output;

    for (i = 0; i < N; ++i) {
        x[i] = 3.;
        l[i] = (i % 2 == 0) ? 1. : -100.;
        u[i] = 100.;
        nbd[i] = 2;
    }
    memset(iwa, 0, sizeof(iwa));
    memset(isave, 0, sizeof(isave));
    memset(lsave, 0, sizeof(lsave));
    memset(wa, 0, sizeof(wa));
    memset(dsave, 0, sizeof(dsave));

    for (i = 0; i < 10000; ++i) {
        setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa,
               &task, &iprint, &csave, lsave, isave, dsave);
        if (IS_FG(task)) {
            rosenbrock(x, &f, g);
        } else if (task != NEW_X) {
            break;
        }
    }

    memcpy(output.x, x, sizeof(x));
    output.f = f;
    output.task = task;
    return output;
}

static void *worker(void *opaque)
{
    worker_args *args = (worker_args *)opaque;
    int round;

    for (round = 0; round < ROUNDS; ++round) {
        pthread_barrier_wait(args->start);
        args->results[round * THREADS + args->index] = optimize();
    }
    return NULL;
}

static int matches(const result *actual, const result *expected)
{
    int i;

    if (!IS_CONVERGED(actual->task) || actual->f > 1e-8) {
        return 0;
    }
    if (fabs(actual->f - expected->f) > 1e-12) {
        return 0;
    }
    for (i = 0; i < N; ++i) {
        if (fabs(actual->x[i] - expected->x[i]) > 1e-10) {
            return 0;
        }
    }
    return 1;
}

int main(void)
{
    pthread_t threads[THREADS];
    worker_args args[THREADS];
    result results[THREADS * ROUNDS];
    result expected = optimize();
    pthread_barrier_t start;
    int i;

    if (!matches(&expected, &expected)) {
        fprintf(stderr, "serial reference optimization did not converge "
                "(task=%ld, f=%.17g)\n", (long)expected.task, expected.f);
        return 1;
    }
    if (pthread_barrier_init(&start, NULL, THREADS) != 0) {
        perror("pthread_barrier_init");
        return 1;
    }
    for (i = 0; i < THREADS; ++i) {
        args[i].start = &start;
        args[i].results = results;
        args[i].index = i;
        if (pthread_create(&threads[i], NULL, worker, &args[i]) != 0) {
            perror("pthread_create");
            return 1;
        }
    }
    for (i = 0; i < THREADS; ++i) {
        pthread_join(threads[i], NULL);
    }
    pthread_barrier_destroy(&start);

    for (i = 0; i < THREADS * ROUNDS; ++i) {
        if (!matches(&results[i], &expected)) {
            fprintf(stderr, "parallel optimization %d disagreed with serial result "
                    "(task=%ld, f=%.17g)\n", i, (long)results[i].task,
                    results[i].f);
            return 1;
        }
    }
    puts("parallel setulb regression test passed");
    return 0;
}
