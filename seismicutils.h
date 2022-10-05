typedef struct{
	unsigned n;
	double  delta;
	double* data;
}seis_t;

typedef struct{
	double semblance;
	double alpha;
	double vel;
}crs_t;

unsigned CommonReceiver(double* s,int* geos,double* h,double ds,double dh,unsigned nh,double g);
crs_t CRSCommonShot(double m0,double t0,const unsigned nt,float *u,seis_t* h,seis_t* s,double w,double dt,double alphamin,double alphamax,double dalpha,
				   double v0,double tdown,double tup,unsigned hyperbolaJump);
