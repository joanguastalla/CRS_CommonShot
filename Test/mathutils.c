#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Function Sgn()
 * 
 *	Parameters:
 *		x (input): float number
 *
 * 	Return: 
 * 		signal: abs(x)
 *
 * 	Description: Function sgn takes a float number and returns
 * 				 its absolute values.
 *
 * */

int Sgn(float x){
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



void TransposeSquare(unsigned n,double (*a)[n]){
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



float* TransposeGeneral(int nlines,int ncols,float (*a)[ncols]){
    float* b=(float*)malloc(nlines*ncols*sizeof(float));
    for(int jj=0;jj<ncols;jj++)
        for(int ii=0;ii<nlines;ii++)    
            b[nlines*jj+ ii]=a[ii][jj];
    return b;
}

/*	Function Linspace()
 *		
 *	 Parameters:
 *		s(input|output):	Data to be filled in 
 *		s0(input):			Initial value
 *		ds(input):  		Sampling of values 
 *		 n(input):			Number of samples
 *
 *	Return:
 *		None
 *
 *	Description: Arrange values linearly in an array.
 *
 * */


void Linspace(float* s,float s0,float ds,unsigned n){
	for(int ii=0;ii<n;ii++)
		s[ii]=s0 + ii*ds;
}

void INTLinspace(int* s,int s0,int ds,unsigned n){
	for(int ii=0;ii<n;ii++)
		s[ii]=s0 + ii*ds;
}
/* Function Arraymerge() 
 *	
 *	Parameters:
 *		m(input|output): 			Data to be filled with two arrays
 *		arrayOne(input):	First array to be merged
 *		arrayTwo(input):	Second array to be merged
 *		n1(input):			Length of arrayOne
 *		n2(input):			Length of arrayTwo
 *
 * Return:
 * 		None
 *
 * Description: Merges two arrays in a single array in the order which they 
 * 				are past.
 *
 * Obs: Array m must be allocated to be of size n1+n2 before doing function 
 *		call.
 *
 * */



void Arraymerge(float* m,float* arrayOne,float* arrayTwo,unsigned n1,unsigned n2){
	for(int ii=0;ii < n1;ii++)
		m[ii]=arrayOne[ii];
	for(int jj=n1;jj<(n1+n2);jj++)
		m[jj]=arrayTwo[jj-n1];
}

void INTArraymerge(int* m,int* arrayOne,int* arrayTwo,unsigned n1,unsigned n2){
	for(int ii=0;ii < n1;ii++)
		m[ii]=arrayOne[ii];
	for(int jj=n1;jj<(n1+n2);jj++)
		m[jj]=arrayTwo[jj-n1];
}

/*	Function Vscale()
 *
 *	 Parameters:
 *	 	ss(output):			Scaled array version of s
 *		s(input): 		    Array to be scaled
 *		ns(input):			Number of samples
 *		factor(input):		Factor to be multiplied in each entry of the array
 *
 *	 Return:
 *	 	None
 *
 *	 Description: Scales an array by a factor
 * */


void Vscale(float* ss,float* s,unsigned ns,float factor){
	for(int ii=0;ii<ns;ii++)
		ss[ii]=s[ii]*factor;
}


/*	Function VsumScalar()
 *
 *	 Parameters:
 *	 	sums(output): 		Scaled array s summed with scalar
 *		s(input):  	 	 	Array to be scaled
 *		ns(input):			Number of samples
 *		scalar(input):		Scalar number to be summed in each entry of the array
 *
 *	 Return:
 *	 	None
 *
 *	 Description: Sum a scalar in each entry of a given array
 *
 * */

void VsumScalar(float* sums,float* s,unsigned ns,float scalar){
	for(int ii=0;ii<ns;ii++)
		sums[ii]=s[ii]+scalar;
}

/*	Function Vpower()
 *
 *	 Parameters:
 *	 	ps(output):			Raised powers of array s
 *		s(input):   		Array to be scaled
 *		ns(input):			Number of samples
 *		power(input):		Power for each exponentiation of the array elements
 *
 *	 Return:
 *	 	None
 *
 *	 Description: Raise to the given power each element of an array
 *
 * */

void Vpower(float* ps,float* s,unsigned ns,float power){
	for(int ii=0;ii<ns;ii++)
		ps[ii]=pow(s[ii],power);
}

/*	Function Vsqrt()
 *
 *	 Parameters:
 *	 	sqrts(output):		Array of square roots of s
 *		s(input):    		Array to be scaled
 *		ns(input):			Number of samples
 *
 *	 Return:
 *	 	None
 *
 *	 Description: Calculate square root of each element in the array
 *
 * */

void Vsqrt(float* sqrts,float* s,unsigned ns){
	for(int ii=0;ii<ns;ii++)
		sqrts[ii]=sqrt(s[ii]);
}

/*	Function Vsum()
 *
 *	 Parameters:
 *	 	result(output): 	Array resulting of summing two arrays
 *		arrayOne(input):	First array
 *		arrayTwo(input):	Second array
 *		ns(input):			Length of both arrays
 *
 *	Return:
 *		None
 *
 *	Description: Given two arrays of the same length, this function
 *				 gives the entrywise sum array. 
 *
 * */


void Vsum(float* result,float* arrayOne,float* arrayTwo,unsigned ns){
	for(int ii=0;ii<ns;ii++)
		result[ii]=arrayOne[ii] + arrayTwo[ii];

}

/*	Function Vabs()
 *
 *	 Parameters:
 *	 	abs_s(output): 		Array of entrywise absolute value
 *		s(input):			Array input
 *		ns(input):			Length of the array
 *
 *	Return:
 *		None
 *
 *	Description: Return array of entrywise absolute value 
 *
 * */

void Vabs(float* abs_s,float* s,unsigned ns){
	for(int ii=0;ii<ns;ii++)
		abs_s[ii]=fabs(s[ii]);
}













