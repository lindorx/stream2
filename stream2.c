#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

#define NMIN 30
#define NMAX 25000000
#define NTIMES 10
#define NUMSIZES 32
#define NPAD 5

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

double mysecond()
{

    struct timeval tp;
    struct timezone tzp;
    int i;

    i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

double sum, start, finish;
int inner;
static double a[NMAX + NPAD], b[NMAX + NPAD];
static double time[4][NTIMES], scalar;
static double rate[4], besttime[4], bytes[4] = {8, 16, 24, 16};

void s_fill(int k, int i, int M);
void s_copy(int k, int i, int M);
void s_daxpy(int k, int i, int M);
void s_dot(int k, int i, int M);

int main()
{
    double exp1, tdelta;
    char ALLTIMES = 0;
    register int i, j, k, l, M;
#pragma omp parallel for
    for (i = 0; i < MIN(10000, NMAX); ++i)
    {
        a[i] = 0.0e0;
    }
    for (i = 0; i < MIN(10000, NMAX); ++i)
    {
        a[i] = mysecond();
    }
    tdelta = 1e36;
    for (i = 0; i < MIN(10000, NMAX) - 1; ++i)
    {
        if (a[i + 1] != a[i])
        {
            tdelta = MIN(tdelta, fabs(a[i + 1] - a[i]));
        }
    }
    printf("Smallest time delta is %g\n", tdelta);
    printf("%8s%8s%10s%10s%10s%10s\n", "Size", "Iter", "FILL", "COPY", "DAXPY", "DOT");
    for (j = 0; j < NUMSIZES; ++j)
    {
        exp1 = log10((double)NMIN) + (double)(j) / (double)(NUMSIZES - 1) *
                                         (log10((double)(NMAX)) - log10((double)(NMIN)));
        M = (unsigned int)(pow((double)10, exp1));
#pragma omp parallel for
        for (i = 0; i < M; ++i)
        {
            a[i] = 0.0e0;
            b[i] = 0.0e0;
        }
        for (k = 0; k < NTIMES; ++k)
        {
            inner = NMAX / M;
            s_fill(k, i, M);
            s_copy(k, i, M);
            s_daxpy(k, i, M);
            s_dot(k, i, M);
        }
        for (i = 0; i < 4; ++i)
        {
            besttime[i] = 1e+36;
            for (k = 0; k < NTIMES; ++k)
            {
                besttime[i] = MIN(besttime[i], time[i][k]);
                if (ALLTIMES)
                    printf("%f\t%f\t%f\t\n", i, k, time[i][k]);
            }
            rate[i] = (double)M * bytes[i] / besttime[i] / 1e6;
        }
        printf("%8d%8d%10.1lf%10.1lf%10.1lf%10.1lf%10.1lf\n", M, NTIMES, rate[0], rate[1], rate[2], rate[3], tdelta / besttime[0]);
    }

    return 0;
}

void s_fill(int k, int i, int M)
{
    int l;
    start = mysecond();
#pragma omp parallel for
    for (l = 0; l < inner; ++l)
    {
        scalar = (double)(k + l + 1);
        for (i = 0; i < M; ++i)
        {
            a[i] = scalar;
        }
    }
    finish = mysecond();
    time[0][k] = (finish - start) / (double)inner;
}

void s_copy(int k, int i, int M)
{
    int l;
    start = mysecond();
#pragma omp parallel for
    for (l = 0; l < inner; ++l)
    {
        a[l] = 1.0e0;
        for (i = 0; i < M; ++i)
        {
            b[i] = a[i];
        }
    }
    finish = mysecond();
    time[1][k] = (finish - start) / (double)inner;
}

void s_daxpy(int k, int i, int M)
{
    int l;
    start = mysecond();
#pragma omp parallel for
    for (l = 0; l < inner; ++l)
    {
        a[l] = 1.0e0;
        for (i = 0; i < M; ++i)
        {
            b[i] = b[i] + scalar * a[i];
        }
    }
    finish = mysecond();
    time[2][k] = (finish - start) / (double)inner;
}

void s_dot(int k, int i, int M)
{
    int l;
    double sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7;
    start = mysecond();
#pragma omp parallel for
    for (l = 0; l < inner; ++l)
    {
        b[l] = 1.0e0;
        sum0 = 0.0e0;
        sum1 = 0.0e0;
        sum2 = 0.0e0;
        sum3 = 0.0e0;
        sum4 = 0.0e0;
        sum5 = 0.0e0;
        sum6 = 0.0e0;
        sum7 = 0.0e0;
        for (i = 0; i < M; i += 8)
        {
            sum0 = sum0 + a[i + 0] * b[i + 0];
            sum1 = sum1 + a[i + 1] * b[i + 1];
            sum2 = sum2 + a[i + 2] * b[i + 2];
            sum3 = sum3 + a[i + 3] * b[i + 3];
            sum4 = sum4 + a[i + 4] * b[i + 4];
            sum5 = sum5 + a[i + 5] * b[i + 5];
            sum6 = sum6 + a[i + 6] * b[i + 6];
            sum7 = sum7 + a[i + 7] * b[i + 7];
        }
    }
    sum = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7;
    finish = mysecond();
    time[3][k] = (finish - start) / (double)inner;
}