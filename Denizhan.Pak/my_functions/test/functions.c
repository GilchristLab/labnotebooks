#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//Continued fractions helper function for upper incomplete gamma
double u_gamma_helper(double s, double x)
{
    double rv;
    int i;

    rv = 10000.0 / x;
    for(i = 10000; i > 0; i--){
        if(i % 2 == 0) rv = (double) ((int) (i / 2)) / (x + rv);
        if(i % 2 != 0) rv = ((double) (((int) (i / 2) + 1) - s)) / (1 + rv);
    }

    return x + rv;
}

//Upper incomplete gamma function
double u_gamma(double s, double x)
{
    //int i;
    double d, rv;

    rv = pow(x, s) * exp(0 - x);

    d = u_gamma_helper(s, x);

    return rv/d;

}

//log of upper incomplete gamma function
double u_gamma_log(double s, double x)
{
    double rv, d;

    rv = s * log(x) - x;
    d = log(u_gamma_helper(s, x));

    return rv - d;

}

//Generalized integral function
double generalized_integral(double p, double z){
    return pow(z, p - 1.0) * u_gamma(1.0 - p, z);
}

//Log of generalized integral function
double generalized_integral_log(double p, double z){
    double co;

    co = p - 1.0;
    return co * log(z) + u_gamma_log(1.0 - p, z);
}

//Calculation of the probability of elongation at current codon
double prob_elongation(double curralpha, double currlambda, double currv){
    return exp(currlambda * currv) * generalized_integral(curralpha, currlambda * currv);
}

//Log probability of elongation at current codon
double prob_elongation_log(double curralpha, double currlambda, double currv){
    double val1, val2;

    val1 = currlambda * currv;
    val2 = generalized_integral_log(curralpha, currlambda * currv);

    return val1 + val2;
}

double delta_g(int i, int g, double *lambda, double *v_g, double *alpha){
    int j;
    double sum = 0;
    double product = 1;

    for(j = 0; j < i; j++){
        sum += lambda[j] * v_g[j];
    }

    for(j = 0; j < i; j++){
        product *= generalized_integral(lambda[j], v_g[j]);
    }

    return exp(sum) * product;
}

double delta_g_log(int i, int g, double *lambda, double *v_g, double *alpha){
    int j;
    double sum = 0;
    double product = 0;

    for(j = 0; j < i; j++){
        sum += lambda[j] * v_g[j];
    }

    for(j = 0; j < i; j++){
        product += generalized_integral_log(lambda[j], v_g[j]);
    }

    return sum + product;
}

double prob_Y_g(double curralpha, int sample_size, double lambda_prime, double psi, double prevdelta){
    double term1, term2, term3;

    term1 = tgamma(curralpha + sample_size) / tgamma(curralpha);
    term2 = psi * prevdelta / (lambda_prime + (psi * prevdelta));
    term3 = lambda_prime / (lambda_prime + (psi * prevdelta));
    
    term2 = pow(term2, sample_size);
    term3 = pow(term3, curralpha);

    return term2 * term3 * term1;
}

double prob_Y_g_log(double curralpha, int sample_size, double lambda_prime, double psi, double prevdelta){
    double term1, term2, term3;

    term1 = lgamma(curralpha + sample_size) - lgamma(curralpha);
    term2 = log(psi) + log(prevdelta) - log(lambda_prime + (psi * prevdelta));
    term3 = log(lambda_prime) - log(lambda_prime + (psi * prevdelta));

    term2 *= sample_size;
    term3 *= curralpha;

    return term1 + term2 + term3;
}
