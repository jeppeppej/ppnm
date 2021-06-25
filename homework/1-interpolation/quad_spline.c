#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>

typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline* qspline_alloc(int n, double* x, double* y){
	qspline* s = malloc(sizeof(qspline));

	s->b = malloc((n-1)*sizeof(double));
	s->c = malloc((n-1)*sizeof(double));
	// data
	s->x = malloc(n*sizeof(double));
	s->y = malloc(n*sizeof(double));
	s->n = n;
	// into spline
	for(int i=0; i<n; i++) {
		s->x[i] = x[i];
		s->y[i] = y[i];
	}
	int i;
	double p[n-1];
	double h[n-1];
	//coefficients

	for(i=0; i<n-1; i++) {
		h[i]=x[i+1]-x[i];
		p[i]=(y[i+1]-y[i])/h[i];
	}

	s->c[0]=0;
	for(i=0; i<n-2; i++) {
		s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
    }

	s->c[n-2]/=2;
	for(i=n-3; i>=0; i--) {
		s->c[i]=(p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
    }


	for(i=0; i<n-1; i++) {
		s->b[i]=p[i]-s->c[i]*h[i];
    }
	return s;
}

double qspline_eval(qspline *s, double z) {
	assert(z>=s->x[0] && z<=s->x[s->n-1]);

    // Binary search
	int i = 0;
    int j = s->n-1;
	while(j - i > 1) {
		int m = (i + j)/2;
		if(z>s->x[m]) i=m;
		else j = m;
	}
	double h = z - s->x[i];
	return s->y[i] + h*(s->b[i] + h*s->c[i]);
}

void qspline_free(qspline *s) {
	free(s->x);
	free(s->y);
    	free(s->b);
    	free(s->c);
    	free(s);
}

// now the integral
double qspline_integ(qspline* s, double z){
	// Binary search
	int i = 0;
    int j = s->n-1;
	while(j - i > 1) {
		int m = (i + j)/2;
		if(z>s->x[m]) i=m;
		else j = m;
	}
	double integral = 0;
	for(j=0;j<i;++j){
		double h = s->x[j+1]-s->x[j];
		integral += s->y[j]*(h) + 0.5*s->b[j]*(h*h) + 1.0/3*s->c[j]*pow(h,3);
	}
	double h = z - s->x[i];
	integral += s->y[i]*h + 1.0/2*s->b[i]*h*h + 1.0/3*s->c[i]*pow(h,3);
	return integral;
}

// and derivative
double qspline_deriv(qspline* s, double z){
	// Binary search
	int i = 0;
    int j = s->n-1;
	while(j - i > 1) {
		int m = (i + j)/2;
		if(z>s->x[m]) i=m;
		else j = m;
	}
	double h = z-s->x[i];
	return s->b[i]+2*h*s->c[i];
}
