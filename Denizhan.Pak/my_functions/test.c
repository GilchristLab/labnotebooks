#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mygamma.h"
#include <gsl/gsl_sf_gamma.h>

/*Description: Creates a log file with cuurrent time and date as name*/
FILE *makeLog(){
    time_t now;
    struct tm *date;
    char *fname, buf[1000];

    date = (struct tm *) malloc(sizeof(struct tm));
    now = time(NULL);
    date = localtime(&now);
    sprintf(buf, "Log.%i.%i.%i-%i:%i", date->tm_mday,
                                       date->tm_mon + 1,
                                       date->tm_year + 1900,
                                       date->tm_hour,
                                       date->tm_min);
    fname = buf;
    return fopen(buf, "w");
}

double compare(double theirs, double mine, FILE *log, int iteration, double s, double x){
    double accurracy;

    accurracy = (mine - theirs) / theirs;
    if(log != NULL){
        fprintf(log, "Iteration %i: \nInputs: s = %f x = %f\nTheir value: %f\nMy value: %f\nDifference: %f\nError: %f\n\n",
                iteration, s, x, theirs, mine, mine - theirs, accurracy);
    }
        
    return accurracy;
}

int main(int argc, char **argv){
    FILE *log;
    int size, seed, i;
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

    if(size < 1) size = 200;
    if(seed == -1) seed = time(NULL);

    sum = 0;
    srand(seed);
    log = makeLog();
    fprintf(log, "The sample size was %i, The seed was %i\n", size, seed);

    for(i = 0; i < size; i++){
        do{
            x = (double)rand()/RAND_MAX + rand() % 10;
            s = (double)rand()/RAND_MAX*2.0-1.0;
        } while(x == 0);
        acc = compare(gsl_sf_gamma_inc(s, x), gamma_inc(s, x), log, i, s, x);
        sum += acc;
    }

    acc = sum / size;
    fprintf(log, "The average error is %f\n", acc); 
}
