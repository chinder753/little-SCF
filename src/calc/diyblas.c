#include "diyblas.h"

#include <stdio.h>

void spprintf(const int n, const char str[], const double * const a){
    printf("%s\n", str);
    for(int i=0; i<n; i++){
        for(int j=0; j<i+1; j++) {
            printf("%15e", a[i*(i+1)/2+j]);
        }
        printf("\n");
    }
}

void geprintf(const int n, const char str[], const double * const a){
    printf("\n%s\n", str);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++) {
            printf("%15e", a[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}

int sp2sy(const int n
          , const double * const x
          , double * const y
          ){
    for(int i=0; i<n; i++){
        y[i*n+i] = x[i*(i+1)/2+i];
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<i; j++){
            y[j*n+i] = x[i*(i+1)/2+j];
        }
    }
}



int sy2sp(const int n
          , const double * const x
          , double * const y
          ){
    for(int i=0; i<n; i++){
        for(int j=i; j<n; j++){
            y[j*(j+1)/2+i] = x[i*n+j];
        }
    }
}


int tran(const int n
        , const double * const x
        , double * const y
        ){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            y[j*n+i] = x[i*n+j];
        }
    }
}