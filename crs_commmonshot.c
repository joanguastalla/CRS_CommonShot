#include <stdio.h>
#include <math.h>
#include "mathutils.h"
#ifndef max
#define max(a,b) (((a)>(b))? (a):(b))
#endif
#ifndef min
#define min(a,b) (((a)>(b))? (b):(a))
#endif
typedef struct{
	unsigned n;
	float  delta;
	float* data;
}seis_t;
void CRSCommonShot(float m0,float t0,float* u,seis_t* h,seis_t* s,float w,float dt,float alphamin,float alphamax,float dalpha,float v0,float tdown,float tup,unsigned hyperbola_jump){
	unsigned shotind;
	float ndown,nup,nsamples;
	float geos[(int)(h->n/(s->delta/h->delta))];
	ndown=min(t0/dt,(int)(tdown/dt));
	nup=(int)(tup/dt);
	nsamples=nup+ndown;
	shotind=(int)(m0/s->delta)*h->n;
	lengeos=CommonReceiver(s->data,geos,h->data,s->delta,h->delta,h->n,m0);
	


}


