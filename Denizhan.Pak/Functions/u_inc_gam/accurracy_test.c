#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <gsl/gsl_sf_gamma.h>

#include "mygamma.h"

/*Description: Creates a log file with cuurrent time and date as name*/
FILE *makeLog(){
    time_t now;
    struct tm *date;
    char *fname, buf[1000];

    date = (struct tm *) malloc(sizeof(struct tm));
    now = time(NULL);
    date = localtime(&now);
    sprintf(buf, "Log.%i.%i.%i-%i:%i:%i", date->tm_mday,
                                       date->tm_mon + 1,
                                       date->tm_year + 1900,
                                       date->tm_hour,
                                       date->tm_min,
                                       date->tm_sec);
    fname = buf;
    return fopen(buf, "w");
}

double compare(double theirs, double mine, FILE *log, int iteration, double s, double x){
    double accurracy;

    accurracy = (mine - theirs) / theirs;
    if(log != NULL){
        fprintf(log, "Iteration %i: \nInputs: s = %.20f x = %.20f\nTheir value: %.20f\nMy value: %.20f\nDifference: %.20f\nError: %.20f\n\n",
                iteration, s, x, theirs, mine, mine - theirs, accurracy);
    }
        
    return accurracy;
}

int main(int argc, char **argv){
    FILE *log;
    int size, seed, i, j;
    double s, x, sum, acc;

    if(argc > 5 || (argc != 1 && argc % 2 == 0)){
        printf("usage: gamma_test [-size i] [-seed r]\n");
        printf("    -size i+: A positive integer for iteration count\n");
        printf("    -seed ui: An unsigned inteeger seed for srand() (-1 will result in TIME being used)\n");
        return -1;
    }

    size = seed = -1;

    for(i = 1; i < argc; i++){
        if(strcmp(argv[i], "-size") == 0){
            size = atoi(argv[(i + 1) % argc]);
        }
        
        if(strcmp(argv[i], "-seed") == 0){
            seed = atoi(argv[(i + 1) % argc]);
        }
    }

    if(size < 1) size = 100;
    if(seed == -1) seed = time(NULL);

    sum = 0;
    srand(seed);
    log = makeLog();
    fprintf(log, "The sample size was %i, The seed was %i\n", size, seed);

    for(i = j = 0; i < size; i++){
        do{
            x = (double)rand()/RAND_MAX;
            s = (double)rand()/RAND_MAX*2.0-1.0;
        } while(x <= s);
        acc = compare(gsl_sf_gamma_inc(s, x), gamma_inc(s, x), log, i, s, x);
        if(acc < 0) acc = 0.0 - acc;
        sum += acc;
    }

    acc = sum / size;
    fprintf(log, "\nThe average error is %.20f\n", acc); 
    printf("\nThe average error is %.20f\n", acc); 
}
