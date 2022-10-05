#include <stdio.h>
#include "../mathutils.h"


int main(){
	float x=-0.47898;
	float v[10];
	float w[10];
	float merge[10];
	unsigned len=10;
	float dv=0.1;
	float dw=0.5;
	Linspace(v,0,dv,len);
	Linspace(w,0,dw,len);
	Arraymerge(merge,v,w,len,len);
	for(int ii=0;ii<2*len;ii++)
		printf("%f\n",merge[ii]);
	//printf("O valor absoluto de %f é:%d\n",x,Sgn(x));


	return 0;
}
