#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Function Sgn()
 * 
 *	Parameters:
 *		x (input): Double number
 *
 * 	Return: abs(x)
 *
 * 	Description: Function sgn takes a double number and returns
 * 				 its absolute values.
 *
 * */

int Sgn(double x){
	int signal;
	int xx=(x>0)?ceil(x):floor(x);
	signal=(xx!=0)?abs(xx)/xx:0;
	return signal;
}

/* 	Function TranposeSquare()
 *
 *	Parameters:
 *		n (input): Square matrix dimension
 *		a (input|output): Pointer to matrix that will be transposed
 *
 * 	Return: 
 * 		None
 *		
 *	Description: Given input pointer of a square matrix it returns the
 *				 transposed matrix
 * */



void TransposeSquare(int n,double (*a)[n]){
    double tmp;
    for(int ii=0;ii<n;ii++)
        for(int jj=(ii+1);jj<n;jj++){
            tmp=a[ii][jj];
            a[ii][jj]=a[jj][ii];
            a[jj][ii]=tmp;
        }
}

/*	Function TranposeGeneral()
 *
 *	Parameters:
 *		nlines (input): Number of lines of the matrix
 *		ncols  (input): Number of columns of the matrix
 *		a 	   (input): Pointer to matrix
 * 	Return
 *		b: Tranposed matrix allocated in the Heap
 *
 *	Description: Compute the tranpose of a general matrix A.
 *
 * */



double* TransposeGeneral(int nlines,int ncols,double (*a)[ncols]){
    double* b=(double*)malloc(nlines*ncols*sizeof(double));
    for(int jj=0;jj<ncols;jj++)
        for(int ii=0;ii<nlines;ii++)    
            b[nlines*jj+ ii]=a[ii][jj];
    return b;
}

/*	Function Linspace()
 *		
 *	 Parameters:
 *		s(output):	Data to be filled in 
 *		s0(input):	Initial value
 *		ds(input):  Sampling of values 
 *		 n(input):	Number of samples
 *
 *	Return:
 *		None
 *
 *	Description: Arrange values linearly in an array.
 *
 * */


void Linspace(double* s,double s0,double ds,unsigned n){
	for(int ii=0;ii<n;ii++)
		s[ii]=s0 + ii*ds;
}



