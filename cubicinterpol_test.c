#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double Cubicinterpol(int len,double dt,float* s,double t);
int main(){
	FILE* seismicData;
	FILE* Interpolation;
	float* u;
	int nt,nh,ns,sizeBytes;
	double dt,amp;
	nt=1126;
	nh=314;
	ns=694;
	dt=0.008;
	u=malloc(sizeof(float)*nt*nh*ns);
	seismicData=fopen("pluto.bin","rb");
	Interpolation=fopen("./Fig/Interpolated.txt","w");
	sizeBytes=4;
	fread(u,sizeBytes,nt*nh*ns,seismicData);
	for(int tt=1;tt<4*(nt-1);tt++){
		amp=Cubicinterpol(nt,dt,(u+314*1126),tt*dt/4);
		fprintf(Interpolation,"%f\n",amp);
	}
	return 0;
}
