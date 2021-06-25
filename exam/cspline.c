#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include"cspline.h"

int binary_search(int n, double * x, double z){

        int i = 0;
        int j = n-1;

        while ( j-i > 1 ){
                int m = (i+j)/2;
                if(z >= x[m]) i = m;
                else j = m;
        }
        return i;
}

/*the cspline struct is {in n; double *x,*y,*x,*b,*c,*d;}*/

cspline* cspline_alloc(int n, double *x, double *y){
	cspline* s = (cspline *) malloc(sizeof(cspline));
	s->x = (double *)malloc(n * sizeof(double));
	s->y = (double *)malloc(n * sizeof(double));
	s->b = (double *)malloc(n * sizeof(double));
	s->c = (double *)malloc((n-1) * sizeof(double));
	s->d = (double *)malloc((n-1) * sizeof(double));
	s->n = n;


	for(int i= 0; i<n; ++i){
		s->x[i]=x[i];
		s->y[i]=y[i];
	}

	double h[n-1], p[n-1];

	for(int i=0; i<n-1;++i){
		h[i] = x[i+1]-x[i];
		assert(h[i]>0);
	}
	for (int i = 0; i<n-1; ++i){
		p[i] = (y[i+1] - y[i]) / h[i];
	}
	double D[n], Q[n-1], B[n];
	D[0]=2;
	for(int i=0; i<n-2; i++){
		D[i+1] = 2* h[i] / h[i+1] +2;
	}
	D[n-1] = 2;
	Q[0] = 1;
	for(int i = 0; i < n-2; i++){
		Q[i+1] = h[i]/h[i+1];
	}
	for(int i = 0; i < n-2; i++){
		B[i+1] = 3*(p[i] + p[i+1] * h[i] / h[i+1]);
	}
	B[0] = 3*p[0];
	B[n-1] = 3 * p[n-2];
	for( int i =1; i<n; i++){
		D[i] -= Q[i-1]/D[i-1];
		B[i] -= B[i-1]/D[i-1];
	}


	s->b[n-1] = B[n-1] / D[n-1];


	for(int i = n-2; i >= 0; i--){
		s->b[i] = (B[i] - Q[i] * s->b[i+1]) / D[i];
	}


	for(int i=0; i < n-1; i++){
		s->c[i] = (-2*s->b[i] - s->b[i+1] + 3* p[i]) / h[i];
		s->d[i] = (s->b[i] + s->b[i+1] - 2*p[i]) / h[i] / h[i];
	}
	return s;
}

double cspline_eval(cspline* s, double z){
	int i = binary_search(s->n, s->x, z);
	double h = z - s->x[i];
	return s->y[i] + h*s->b[i] + h*h*s->c[i] + h*h*h*s->d[i];
}

double cspline_integrate(cspline* s, double z){
	int i = binary_search(s->n, s->x, z);
	double value = 0;
	for(int j = 0; j < i; j++){
		double dx = s->x[j+1]-s->x[j];
		value += dx*s->y[j] + 1.0/2*dx*dx*s->b[j] + 1.0/3*dx*dx*dx*s->c[j] + 1.0/4*dx*dx*dx*dx*s->d[j];
	}
	double h = z - s->x[i];
	return value +=h*s->y[i] + 1.0/2*h*h*s->b[i] + 1.0/3*h*h*h*s->c[i] + 1.0/4*h*h*h*h*s->d[i];
}

double cspline_derive(cspline* s, double z){
	int i = binary_search(s->n, s->x, z);
	double h = z - s->x[i];
	return s->b[i] + 2.0*h*s->c[i] + 3.0*h*h*s->d[i];
}

void cspline_free(cspline* s){
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s->d);
	free(s);
	return;
}
