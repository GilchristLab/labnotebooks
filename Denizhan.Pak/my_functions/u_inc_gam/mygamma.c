#include <math.h>
#include <stdio.h>

double cont_frac(double v, double x, int flag){
    double c, rv;

    if(flag == 15000) return x;

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

double cont_frac_rev(double v, double x, int flag, double rv){

    if (flag == 0) return x + rv;

    if(flag % 2 == 0){
        rv = (double) ((int) (flag / 2)) / (x + rv);
    }
    else{
        rv = ((double) (((int) (flag / 2) + 1) - v)) / (1 + rv);
    }

    return cont_frac_rev(v, x, flag - 1, rv);
}

double cont_frac_it(double v, double x){
    double rv;
    int i;

    rv = 10000.0 / x;
    for(i = 10000; i > 0; i--){
        if(i % 2 == 0) rv = (double) ((int) (i / 2)) / (x + rv);
        if(i % 2 != 0) rv = ((double) (((int) (i / 2) + 1) - v)) / (1 + rv);
    }

    return x + rv;
}

double gamma_inc(double v, double x){
    int i;
    double d, rv;

    rv = pow(x, v) * exp(0 - x);

    //d = cont_frac_rev(v, x, 10000, 10000.0 / x);
    //d = cont_frac(v, x, 0);
    d = cont_frac_it(v, x);

    return rv/d;
}

