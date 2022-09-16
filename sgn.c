#include <stdlib.h>
#include <stdio.h>
#include <math.h>
int sgn(double x){
	int signal;
	int xx=(x>0)?ceil(x):floor(x);
	signal=(xx!=0)?abs(xx)/xx:0;
	return signal;
}
