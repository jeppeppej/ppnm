#include <stdlib.h>
#include <assert.h>
#include <math.h>
typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;

cspline* cspline_alloc(int n, double *x, double *y){
	cspline *s = (cspline*)malloc(sizeof(cspline));
    // coefficients
	s->b = (double*)malloc((n)*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->d = (double*)malloc((n-1)*sizeof(double));
    // data
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->n = n;
	// into spline
	for(int i=0;i<n;i++){
		s->x[i] = x[i];s->y[i] = y[i];
	}

	double p[n-1],h[n-1];
	for(int i=0;i<n-1;i++){
		h[i]=x[i+1]-x[i];
		assert(h[i]>0);
		p[i]=(y[i+1]-y[i])/h[i];
	}

    // gauss elimination and backsubstitution
	double D[n], Q[n-1], B[n];
	D[0] = 2;
	for(int i=0;i<n-2;i++){
		D[i+1]=2*h[i]/h[i+1]+2;
	}
	D[n-1]=2;
	Q[0]=1;
	for(int i=0;i<n-2;i++){
		Q[i+1]=h[i]/h[i+1];
	}
	for(int i=0;i<n-2;i++){
		B[i+1]=3*(p[i]+p[i+1]*h[i]/h[i+1]);
	}
	B[0]=3*p[0]; B[n-1]=3*p[n-2];

	//gauss
	for(int i=1;i<n;i++){
		D[i]-=Q[i-1]/D[i-1];
		B[i]-=B[i-1]/D[i-1];
	}

	s->b[n-1] = B[n-1]/D[n-1];
	for(int i=n-2;i>=0;i--){
		s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];
	}

	for(int i=0;i<n-1;i++){
		s->c[i] = (-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
		s->d[i] = (s->b[i]+s->b[i+1]-2*p[i])/h[i]/h[i];
	}

	return s;
}

void cspline_free(cspline *s){
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s->d);
	free(s);
}

double cspline_eval(cspline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	// Binary search
	int i = 0;
    int j = s->n-1;
	while(j - i > 1) {
		int m = (i + j)/2;
		if(z>s->x[m]) i=m;
		else j = m;
	}
	double h=z-s->x[i];
	return s->y[i] + h*s->b[i] + h*h*s->c[i] + h*h*h*s->d[i];
}

double cspline_integ(cspline* s, double z) {
	assert(z>= s->x[0] && z<= s->x[s->n -1]);
	// Binary search
	int i = 0;
    int j = s->n-1;
	while(j - i > 1) {
		int m = (i + j)/2;
		if(z>s->x[m]) i=m;
		else j = m;
	}
	double integral = 0;
	for(j = 0; j < i; j++) {
		double dx = s->x[j+1] - s->x[j];
		integral += dx*s->y[j] + 1.0/2*dx*dx*s->b[j] + 1.0/3*pow(dx,3)*s->c[j] + 1.0/4*pow(dx,4)*s->d[j];
	}
	double h = z - s->x[i];
	return integral +=h*s->y[i] + 1.0/2*h*h*s->b[i] + 1.0/3*pow(h,3)*s->c[i] + 1.0/4*pow(h,4)*s->d[i];
}


double cspline_deriv(cspline* s, double z) {
	assert(z>= s->x[0] && z<= s->x[s->n -1]);
		// Binary search
	int i = 0;
    int j = s->n-1;
	while(j - i > 1) {
		int m = (i + j)/2;
		if(z>s->x[m]) i=m;
		else j = m;
	}
	double h = z - s->x[i];
	return s->b[i] + 2.0*h*s->c[i] + 3.0*h*h*s->d[i];
}
