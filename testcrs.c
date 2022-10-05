#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mathutils.h"
#include "seismicutils.h"
#define PI 3.14159265359


int main(){
	clock_t start,end;
	float* u;
	double* Alpha;
	double* Vel;
	double feet2km,ds,h0,dh,dt,m0,f,w,alphamin,alphamax,dalpha,v0,tdown,tup;
	unsigned hyperbolaJump,nh,ns,nt,sizeBytes,tt;
	FILE* seismicData;
	FILE* CRS;
	int* traces;
	CRS=fopen("./Fig/CRS.txt","w");
	crs_t crsdata;
	seis_t* s;
	seis_t* h;
	feet2km=0.0003048;
	ds=150*feet2km;
	ns=694;
	nh=314;
	nt=1126;
	dh=75*feet2km;
	h0=36*dh;
	dt=0.008;
	sizeBytes=4;
	m0=300*ds;
	f=15; 					//Frequency in Hz
	w=sqrt(6)/(PI*f);
	alphamin=-10*PI/180;
	alphamax=20*PI/180;
	dalpha=0.5*PI/180;
	v0=1.5;
	tdown=0;
	tup=2.5;
	hyperbolaJump=2;
	u=malloc(sizeBytes*nt*ns*nh);
	seismicData=fopen("pluto.bin","rb");
	fread(u,sizeBytes,nt*nh*ns,seismicData);
	s=malloc(sizeof(crs_t*));
	h=malloc(sizeof(crs_t*));
	s->n=ns;
	s->delta=ds;
	s->data=malloc(sizeof(double)*ns);
	h->n=nh;
	h->delta=dh;
	h->data=malloc(sizeof(double)*nh);
	Alpha=malloc(sizeof(double)*nt);
	Vel=malloc(sizeof(double)*nt);
	traces=malloc(sizeof(int)*nh);
	Linspace(s->data,0,ds,ns);
	Linspace(h->data,h0,dh,nh);
	if(!seismicData){
		printf("Error reading file!!\n");
		exit(EXIT_FAILURE);
	}else{
		printf("Sucess opening file\n");
		fclose(seismicData);
	}
	start=clock();

//	# pragma omp parallel for num_threads(2) \
		default(none) private(tt,crsdata,getsemblance) shared(m0,nt,u,h,s,w,dt,alphamin,alphamax,dalpha,v0,tdown,tup,hyperbolaJump,Semblance,Alpha,Vel)
	

		for(tt=133;tt<134;tt++){
			crsdata=CRSCommonShot(m0,tt*dt,nt,u,h,s,w,dt,alphamin,alphamax,dalpha,v0,tdown,tup,hyperbolaJump);
			printf("semblance: %f\n",crsdata.semblance);
			printf("alpha: %f\n",crsdata.alpha);
			printf("vel: %f\n",crsdata.vel);
			fprintf(CRS,"%f\t\t",crsdata.semblance);
			fprintf(CRS,"%f\t\t",crsdata.alpha);
			fprintf(CRS,"%f\n",crsdata.vel);
		}	
	end=clock();
	printf("Time used:%f seconds\n",(double) (end-start)/CLOCKS_PER_SEC );

	free(s);
	free(h);
	free(Alpha);
	free(Vel);
	free(traces);
	return 0;



}


