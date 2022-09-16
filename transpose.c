#include <stdio.h>
#include <stdlib.h>
#include "mathutils.h"
void transpose_square(int n,double (*a)[n]){
    double tmp;
    for(int ii=0;ii<n;ii++)
        for(int jj=(ii+1);jj<n;jj++){
            tmp=a[ii][jj];
            a[ii][jj]=a[jj][ii];
            a[jj][ii]=tmp;
        }
}

double* transpose_general(int nlines,int ncols,double (*a)[ncols]){
    double* b=(double*)malloc(nlines*ncols*sizeof(double));
    for(int jj=0;jj<ncols;jj++)
        for(int ii=0;ii<nlines;ii++)    
            b[nlines*jj+ ii]=a[ii][jj];
    return b;
}
