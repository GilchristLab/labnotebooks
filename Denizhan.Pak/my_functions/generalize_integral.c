#include <math.h>
double logGeneralizedIntegral(double p, double z){
    double rv = (p - 1) * log(z) + u_gamma_log(1 - p, z);

    return rv;
}
