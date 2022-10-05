#include <stdio.h>
#include <math.h>
#define PI 3.14159265359
void HyperbolaPlot(float t0,float v0,float* h,const int hlen,float alpha,float tf){
	printf("tf: %f\n",tf);
	FILE* HyperbolaPlot;
	HyperbolaPlot=fopen("./Fig/HyperbolasPlot.txt","a");
	float hyperbolaCurve[hlen];
	float vel;
	for(int ii=0;ii<hlen;ii++){
		hyperbolaCurve[ii]=pow((t0 + 2*h[ii]*sin(alpha)/v0),2);
	}
	vel=(pow(tf,2) - hyperbolaCurve[hlen-1])/(pow(h[hlen-1],2));
	for(int ii=0;ii<hlen;ii++){
		hyperbolaCurve[ii]+= vel*pow(h[ii],2);
		hyperbolaCurve[ii]=sqrt(hyperbolaCurve[ii]);
		fprintf(HyperbolaPlot,"%f\n",hyperbolaCurve[ii]);
	}
	

}


int main(){
	int hlen,hyperbolaJump,nsamples;
	float h0,t0,dt,dh,v0,alpha,timeup,feet2km;
	dt=0.008;
	t0=134*dt;
	timeup=3;
	hyperbolaJump=16;
	v0=1.5;
	alpha=20*PI/180;
	hlen=314;
	float h[hlen];
	feet2km=0.0003048;
	h0=18*150*feet2km;
	dh=75*feet2km;
	for(int hh=0;hh<hlen;hh++){
		h[hh]=(h0 + hh*dh)/2; 
	}
	printf("t0: %f\n",t0);
	nsamples=(int) (timeup/dt)/hyperbolaJump;
	for(int cc=1;cc<nsamples;cc+=1){
		HyperbolaPlot(t0,v0,h,hlen,alpha,t0+cc*dt*hyperbolaJump);
	}
	return 0;
}
