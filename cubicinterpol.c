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

void printer(int size1,int size2,double* a){
    if(size1 > 1){
        for(int ii=0;ii<size1;ii++){
            for(int jj=0;jj<size2;jj++)
                printf("%f\t",a[size1*ii + jj]);
            printf("\n");
        }
    }
    else{
       for (int jj=0;jj<size2;jj++)
            printf("%f\n",a[jj]);
    }
}


double cubicinterpol(int len,double dt,double* s,double t){
    int ind,n = N,info,nrhs=1;
    double rest=0;
	double alpha[2];
    double A[N][N];
    double b[N],c[N];
	static double r[N];
    int ipiv[N];
	//t=fabs((sgn(t)+1)*t/2);
	ind=t/dt;
    if((double) ind==t/dt){
		return s[ind];
    } 
    else{
		ind=min(ind,len-2);
		alpha[0]=(ind>0)?(s[ind+1]-s[ind-1])/(2*dt):(s[ind+1]-s[ind])/dt;
		alpha[1]=(ind<(len-2))?(s[ind+2]-s[ind])/(2*dt):(s[ind+1]-s[ind])/dt; 
		 

       	for(int ii=0;ii<2;ii++){
            for(int jj=3;jj>=0;jj--){
                A[ii][3-jj]=pow((ind+ii)*dt,jj);
            }
        }  
       	for(int ii=2;ii<4;ii++){
            for(int jj=0;jj<3;jj++){
                A[ii][jj]=(3-jj)*pow((ind +(ii-2))*dt,2-jj);
            }   
        }
        
        transpose_square(N,A);
        b[0]=s[ind];
        b[1]=s[ind+1];
        b[2]=alpha[0];
        b[3]=alpha[1];
		c=b;
		LAPACK_dgesv(&n,&nrhs,&A[0][0],&n,ipiv,b,&n,&info);
		for(int ii=0;ii<N;ii++){
			for(int jj=0;jj<N;jj++){
				r[ii]+=A[ii][jj]*b[jj];
			}
			rest+=pow((r[ii] - c[ii]),2);
		}
		printf("%f\n",sqrt(rest));
    	return (b[0]*pow(t,3) + b[1]*pow(t,2)+b[2]*t + b[3]);
	}
		
}

void wfunc(char* filename,int len,double* x,double* y){
	char filepath[100]="./Fig/";
	FILE* outfile=fopen(strcat(filepath,filename),"w");
	for(int ii=0;ii<len;ii++){
		fprintf(outfile,"%f %f\n",x[ii],y[ii]);
	}
}

int main(){
    int len=5;
	const int len_t=401;
	double inc_t=0.01,dt=1.0;
    double t[len_t];
	double cubic[len_t];
    double s[]={2.8,3,1,3,4};
	double* p=s;
	char* cubic_data="cubic_data.txt";
	for(int ii=0;ii<len_t;ii++){
		t[ii]=ii*inc_t;
		cubic[ii]=cubicinterpol(len,dt,p,t[ii]);
	}
	wfunc(cubic_data,len_t,t,cubic);
    return 0;
}
