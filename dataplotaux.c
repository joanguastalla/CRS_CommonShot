#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void WriteGraph(char* filename,int len,double* x,double* y){
	char filepath[100]="./Fig/";
	FILE* outfile=fopen(strcat(filepath,filename),"w");
	for(int ii=0;ii<len;ii++){
		fprintf(outfile,"%6f %6f\n",x[ii],y[ii]);
	}
}
