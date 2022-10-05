#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathutils.h"
#include "seismicutils.h"
#ifndef max
#define max(a,b) (((a)>(b))? (a):(b))
#endif
#ifndef min
#define min(a,b) (((a)>(b))? (b):(a))
#endif
#define indmax(a,b) (((a)>(b))? (0) : (1))
#define PI 3.142857


/*	Function CRSCommonShot()
 *
 *	 Parameters:
 *		m0(input): 				Central Zero-Offset(ZO) midpoint
 *		t0(input): 				Traveltime of the ZO ray
 *		nt(input): 				Number of time samples in the seismic section recorded 
 *		u(input):				Common-Shot recorded data
 *		h(input):				Offset structure 
 *		s(input): 				Shots structure
 *		w(input):				Wavelet duration[s]
 *		dt(input):				Time sampling in the recorded data
 *		alphamin(input):		Minimum angle of emergence of the central ZO ray(degrees)	
 *		alphamax(input):   		Maximum angle of emergence of the central ZO ray(degrees)
 *		dalpha(input):			Sampling of emergence angle(degrees)
 *		v0(input): 				Velocity at emergence ZO position at surface (x=m0,z=0)
 *		tdown(input):			Relates to time given by t=t0-tdown in which hyperbola fits at the maximum data offset 
 *		tup(input):				Relates to time given by t=t0+tup in which hyperbola fits at the maximum data offset 
 *		hyperbolaJump(input):	Jump in time samples between (t0-tdown,t0+tup) which hyperbola fits at maximum data offset 
 *
 *
 *	Return:
 *		paramCRS: Structure that contains the greatest found semblance, and the two hyperbola parameters.
 *				  
 *
 *  Description: Computation of DCRS parameters, associated to a NIP diffraction point, from a given zero-offset
 *  			 trace and traveltime (m0,t0). Configuration of the traces is supossed to be Common Shot,linearly 
 *  			 distributed with traces sorted as:
 *
 *	
 *    			 U[:,index]= trace parametrized by the functions s(index) and h(index), calculated by
 *
 *  		     s(index)=(index//len(h))*ds
 *    			 h(index)=h[0] + mod(index,len(h))*dh
 *
 *    			 where,
 *
 *    			 ds=s[1]-s[0]
 *    			 dh=h[1]-h[0]
 *
 *	  			 Hyperbola fitting mathematical expression is given by:
 *
 *	  		     T^2(m,h;alfa,v)=(t0 + 2*sin(alfa)/v0*h)^2 + v*h^2
 *			
 *				 Function work is to estimate, using Semblance as coherency measure, the best alfa and v in a given grid of values.
 *
 * */



crs_t CRSCommonShot(double m0,double t0,const unsigned nt,float *u,seis_t* h,seis_t* s,double w,double dt,double alphamin,double alphamax,double dalpha,
				   double v0,double tdown,double tup,unsigned hyperbolaJump){

	crs_t paramCRS;
	FILE* Plots;
	Plots=fopen("./Fig/Curves.txt","w");
	int ww,nsw,ii,jump,contador,sz;
	unsigned ind,shotind,lengeos,nh;
	int	traces[h->n];
	int geos[(int)(h->n/(s->delta/h->delta))];
	double ndown,nsamples,factor,sumaux,squareAmps,innerAmps,Amp,semblanceMax,alpha;
	double v[2],alfaMax[2],velMax[2],semblance[2];
	sz=140;
	semblanceMax=0;
	lengeos=CommonReceiver(s->data,geos,h->data,s->delta,h->delta,h->n,m0);
	nh=lengeos + h->n;
	double hminus[lengeos];
	// Assignment of negative offsets for a common-receiver section at m0
	Linspace(hminus,-h->data[0]-(lengeos-1)*s->delta,s->delta,lengeos);
	shotind=(int)(m0/(s->delta))*h->n;
	INTLinspace(traces,shotind,1,sz);
	ndown=min(t0/dt,(int)(tdown/dt));
	nsamples=(int)(tup/dt)+ ndown; 
	// semblance windowing size
	nsw=w/dt;
	// Manual allocation of memory for vectorized functions
	double offsets[nh];
	double squareHolder[sz];
	double factorTime[sz];
	int   stackTraces[nh];
	double alphaCurve[sz];
	double dCurve[sz];
	Arraymerge(offsets,hminus,h->data,lengeos,h->n);
	INTArraymerge(stackTraces,geos,traces,lengeos,h->n);
	Vpower(squareHolder,h->data,sz,2);
	//	
	//		Loops in alpha and jump are called parameter loops, since that is all we need 
	//	to build hyperbolic traveltime of Common-Shot sections in a 2D Medium.
	//
	//
	//	!! Index "0" is always reserved to store CRS parameters that gives the greater
	//	   Semblance inside the parameter loops.
	//
	//	!! Index "1" is used to store the actual looping value inside parameter loops.   
	//
	// **Loop jump**
	// semblance[2] - Used to compare two semblance values inside loop "jump" for fixed alpha
	// v[2] - 		  [v(semblance[0]), actual v given jump]
	// 
	// **Loop alpha**
	// semblanceMax - Always store maximum Semblance value until next alpha change in the loop. 
	// velMax[2] - 	  Stores velocity of previous greatest Semblance(semblanceMax) and actual v[0] 	
	// alfaMax[2] -   Stores angles of previous greatest Semblance and the actual looping alpha.
	contador=0;
	for(alpha=alphamin;alpha<alphamax;alpha+=dalpha){
		printf("%d\n",contador);
		semblance[0]=0;
		factor=sin(alpha)/v0;
		Vscale(alphaCurve,h->data,sz,factor);
		VsumScalar(alphaCurve,alphaCurve,sz,t0);
		Vpower(alphaCurve,alphaCurve,sz,2);
		for(jump=0;jump<nsamples;jump+=hyperbolaJump){
			squareAmps=0;
			innerAmps=0;
			v[1]=(pow((t0 +(jump-ndown)*dt),2) - alphaCurve[sz-1])/squareHolder[sz-1];
			Vscale(factorTime,squareHolder,sz,v[1]);
			Vsum(dCurve,alphaCurve,factorTime,sz);
			Vabs(dCurve,dCurve,sz);
			Vsqrt(dCurve,dCurve,sz);
			for(ww=-nsw;ww<=nsw;ww++){
				sumaux=0;
				for(ii=0;ii<sz;ii++){
					Amp=Cubicinterpol(nt,dt,u+nt*traces[ii],dCurve[ii] + ww*dt);	
			//		Amp=(*(u + nt*traces[ii] + (int)floor(dCurve[ii]/dt) + ww) + *(u + nt*traces[ii] + (int)ceil(dCurve[ii]/dt) + ww))/2; 
					squareAmps+=pow(Amp,2);
					sumaux+=Amp;
				}
				innerAmps+=pow(sumaux,2);
		   	}
			semblance[1]=innerAmps/(sz*squareAmps + 1e-40);
			ind=indmax(semblance[0],semblance[1]);
			semblance[0]=semblance[ind];
			v[0]=v[ind];
		}
		velMax[1]=v[0];
		alfaMax[1]=alpha;
		semblanceMax=max(semblanceMax,semblance[0]);
		ind=indmax(semblanceMax,semblance[0]);
		velMax[0]=velMax[ind];
		alfaMax[0]=alfaMax[ind];
		contador+=1;
	}
	paramCRS.semblance=semblanceMax;
	paramCRS.alpha=alfaMax[0];
	paramCRS.vel=velMax[0];
	return paramCRS;

}
