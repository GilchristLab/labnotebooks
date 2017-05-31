#include <math.h>
#include <stdio.h>

double cont_frac(double v, double x, int flag){
    double c, rv;

    if(flag == 10000) return x;

    if(flag % 2 == 0){
        c = x;
        rv = c + ((double)((flag / 2) + 1 - v) / cont_frac(v, x, flag + 1));
    }

    else{
        c = 1;
        rv = c + ((double)((flag / 2) + 1) / cont_frac(v, x, flag + 1));
    }
    return rv;
}

double gamma_inc(double v, double x){
    int i;
    double d, rv;

    rv = pow(x, v) * exp(0 - x);

    d = cont_frac(v, x, 0);

    return rv/d;
}

