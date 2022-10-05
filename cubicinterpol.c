#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lapack.h>
#include "mathutils.h"
#ifndef max
#define max(a,b) (((a)>(b))? (a):(b))
#endif
#ifndef min
#define min(a,b) (((a)>(b))? (b):(a))
#endif
#define N 4
#define DERIVATECOND 2

/* Function Cubicinterpol()
 *
 *  Parameters:
 *      len (input): Length of data used for interpolation
 *      dt  (input): Data sampling 
 *      s   (input): Original data
 *      t   (input): Time for interpolating
 *      
 *  Return:
 *      float: Cubic interpolation of data in time t
 *
 *
 *  Description: With the ind=floor(t/dt) this function compute 
 *               cubic interpolation of data such that at ind
 *               and (ind+1) the values and centered finite dife-
 *               rence is the same.
 *
 *  
 *   s[start],s[end],s'[start],s'[end] 	
 *
 *	where,
 *	
 *	start=floor(t/dt), end=ceil(t/dt)
 *	
 *  and s'[start],s'[end] are centered finite diffrence derivatives, i.e., 
 *	
 *	s'[start]=(s[end]-s[start-1])/(2*dt), s'[end]=(s[end+1]-s[start])/(2*dt) 
 * 
 * */


double Cubicinterpol(int len,double dt,float* s,double t){
    int ind,n = N,info,nrhs=1;
	double alpha[2];
    double A[N][N];
    double b[N];
	double tmp,piso;
    int ipiv[N];
	tmp=t/dt;
	piso=floor(tmp);
	ind=(int)tmp;
	if(t==piso*dt){
		return s[ind];
    } 
    else{
		ind=min(ind,len-2);
		alpha[0]=(s[ind+1]-s[ind-1])/(2*dt);
		alpha[1]=(s[ind+2]-s[ind])/(2*dt);
       	for(int ii=0;ii<DERIVATECOND;ii++){
            for(int jj=(N-1);jj>=0;jj-=1){
                A[ii][(N-1)-jj]=pow((ind+ii)*dt,jj);
            }
        }  
       	for(int ii=DERIVATECOND;ii<N;ii++){
            for(int jj=0;jj<N;jj+=1){
                A[ii][jj]=((N-1)-jj)*pow((ind +(ii-DERIVATECOND))*dt,abs(DERIVATECOND-jj));
            }   
        }
        TransposeSquare(N,A);
		b[0]=s[ind];
        b[1]=s[ind+1];
        b[2]=alpha[0];
        b[3]=alpha[1];
        
		LAPACK_dgesv(&n,&nrhs,&A[0][0],&n,ipiv,b,&n,&info);
		return (b[0]*pow(t,3) + b[1]*pow(t,2)+ b[2]*t + b[3]);
	}
		
}

