#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>

int main(int argc, char **argv){
    double *s, *x;
    int sample, seed;

    if(argc > 3){
        printf("usage: gamma_test [-size] [-seed]\n");
        return -1;
    }

    if(argc > 2){
        //User set flag to determine sample size
        if(strcmp(argv[1], "-size") == 0 || strcmp(argv[2], "-size") == 0){
            printf("Please enter sample size: ");
            scanf("%i", &sample);
        }
        else{
            //Generic sample size is 1000
            sample = 1000;
        }

        //User set flag to determine seed
        if(strcmp(argv[1], "-seed") == 0 || strcmp(argv[2], "-seed") == 0){
            printf("Please enter seed: ");
            scanf("%i", &seed);
        }
        else{
            //Generic seed is 0
            seed = 0;
        }

    }
        printf("The sample was %i, The seed was %i\n", sample, seed);
}
