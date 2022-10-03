#include <stdio.h>
#include <stdlib.h>
#ifndef max
#define max(a,b) (((a)>(b))? (a) : (b))
#endif
/* Function CommonReceiver()
 *
 * 	Parameters:
 * 		s(input): 		Source coordinates 
 *		geos(output):	Index of shots for a common-receiver at position g
 *		h(input):		Offsets
 *		nh(input):		Number of different offsets
 *		g(input):		Coordinate of the common-receiver
 *
 *	Return:
 *		lengeos: Number of offsets in the common-receiver section
 *
 *	Description: Given sets of possible offsets and source positions in a seismic experiment 
 *				 this function computes possible source indexes, in such a way that we have a 
 *				 common-receiver section at g.
 * */


unsigned CommonReceiver(double* s,int* geos,double* h,double ds,double dh,unsigned nh,double g){
	double ratio;
	unsigned geoind,bigoffset,lengeos,smalloffset,ii;
	ratio=ds/dh;
	geoind=(g-s[0])/ds;
	smalloffset=geoind - (int)(h[0]/ds);
	bigoffset=max(0,smalloffset - (int)((nh-1)/ratio));
	lengeos=smalloffset - bigoffset+1;
	for(ii=0;ii<lengeos;ii++)
		geos[ii]=(bigoffset+ii)*nh + (lengeos-ii-1)*ratio;
	
	return lengeos;
}


