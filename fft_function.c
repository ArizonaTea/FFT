#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <iostream>

#include "complex.h"

#define Pi 3.14159265358979

complex* FFT(complex *,int,int);
complex* fftshift(complex *,int);

int main(void)
{
    //Testing.....
}




complex* FFT(complex *a,int N,int sign)//sign=-1: FFT; sign=1: Inverse FFT
{
    int m,mr,l,i;
    complex T;
    mr=0;
    for(m=1;m<=N-1;m++){
        l=N;
        do{
            l /=2;
        }while(mr+l>=N);
        mr=(int)fmod(mr,l)+l;
        if(mr>m){
            T=*(a+m);
            *(a+m)=*(a+mr);
            *(a+mr)=T;
        }
    }
    l=1;
    while(l<N){
        for(m=0;m<=l-1;m++){
            for(i=m;i<=N-1;i +=2*l){
                T=*(a+i+l)*exp(complex(0.,m*Pi*sign/(double)l));
                *(a+i+l)=*(a+i)-T;
                *(a+i) +=T;
            }
        }
        l=2*l;
    }
    if(sign==1) for(i=0;i<=N-1;i++) *(a+i) /=N;
    return a;
}


double* fftshift(double *A,int N)
{
    int i;
    double temp;
    for(i=0;i<=N/2-1;i++){
        temp=*(A+i);
        *(A+i)=*(A+N/2+i);
        *(A+N/2+i)=temp;
    }
    return A;
}

complex* fftshift(complex *A,int N)
{
    int i;
    complex temp;
    for(i=0;i<=N/2-1;i++){
        temp=*(A+i);
        *(A+i)=*(A+N/2+i);
        *(A+N/2+i)=temp;
    }
    return A;
}
