#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Referenced from
Chi Leung Lau,
Algorithm AS 147:
A Simple Series for the Incomplete Gamma Integral,
Applied Statistics,
Volume 29, Number 1, 1980, pages 113-114.
double l_gamma_norm ( double x, double p){
C++ Implementation referenced from John Burkardt

A Functions that returns the normalized lower gamma function.

Parameters:

Input, double X, P, the arguments of the incomplete
Gamma integral.  X and P must be greater than 0.

Output, double l_gamma_norm, the value of the incomplete
Gamma integral.
*/
double l_gamma_norm(double s, double x)
{
    double a;
    double arg;
    double c;
    double e = 1.0E-09;
    double f;
    double uflo = 1.0E-37;
    double value;

    if (x <= 0)
    {
        value = 0;
        return value;
    }

    if (s <= 0)
    {
        value = 0;
        return value;
    }

    arg = s * log(x) - lgamma(s + 1.0) - x;

    if ( arg < log(uflo))
    {
        value = 0;
        return value;
    }

    f = exp(arg);

    if (f == 0)
    {
        value = 0;
        return value;
    }

    c = 1.0;
    value = 1.0;
    a = s;

    while(1)
    {
        a = a + 1.0;
        c = c * x / a;
        value = value + c;

        if ( c <= e * value ) break;
    }

    value = value * f;
    return value;
}

double u_gamma(double s, double x){
    double rv;

    rv = tgamma(s) * (1 - l_gamma_norm(s, x));
    return rv;
}

double u_gamma_log(double s, double x){
    double rv;

    rv = tgamma(s) * (1 - l_gamma_norm(s, x));
    return log(rv);
}

double logGeneralizedIntegral(double p, double z){
    double rv = (p - 1) * log(z) + u_gamma_log(1 - p, z);

    return rv;
}

int main(int argc, char ** argv){
    int a, b;
    int i;

    for(i = 0; i < 1000000; i++){
    a = rand() % 30;
    b = rand() % 40;

    /*printf("a = %i, b = %i\n",a,b );
    printf("IGamma = %f\n", u_gamma(a, b));
    printf("IGamma = %f\n", u_gamma_log(a, b));
    printf("\n");*/
    }
}
