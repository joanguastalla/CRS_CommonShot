void TransposeSquare(unsigned n,double (*a)[n]);
int Sgn(double x);
double* TransposeGeneral(int nlines,int ncols,double (*a)[ncols]);
void Linspace(double* s,double s0,double ds,unsigned n);
void INTLinspace(int* s,int s0,int ds,unsigned n);
void Arraymerge(double* m,double* arrayOne,double* arrayTwo,unsigned n1,unsigned n2);
void INTArraymerge(int* m,int* arrayOne,int* arrayTwo,unsigned n1,unsigned n2);
void Vscale(double* ss,double* s,unsigned ns,double factor);
void VsumScalar(double* sums,double* s,unsigned ns,double scalar);
void Vpower(double* ps,double* s,unsigned ns,double power);
void Vsqrt(double* sqrts,double* s,unsigned ns);
void Vsum(double* result,double* arrayOne,double* arrayTwo,unsigned ns);
void Vabs(double* abs_s,double* s,unsigned ns);
double Cubicinterpol(int len,double dt,float* s,double t);

