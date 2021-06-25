#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"aspline.h"

aspline* aspline_alloc(int n, double* x, double* y, double* y_prime){
	if(n<3){fprintf(stderr,"Too few points passed to cspline_w_derives_alloc to construct spline.\n");return NULL;}
	double p[n-1], dx[n-1];
	for(int i=0; i<n-1; ++i){dx[i] = x[i+1]-x[i];if(dx<=0){fprintf(stderr,"Negative or 0 dx stepsize passed to cspline_w_derives_alloc.\n");return NULL;}}
	for(int i=0; i<n-1; ++i){p[i] = (y[i+1]-y[i])/dx[i];}
	aspline* s = (aspline*)malloc(sizeof(aspline));
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->b = (double*)malloc(n*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->d = (double*)malloc((n-1)*sizeof(double));
	s->n = n; for(int i=0; i<n; ++i){s->x[i]=x[i];s->y[i]=y[i];s->b[i]=y_prime[i];}
	for(int i=0; i<n-1; ++i){
		s->c[i] =(3*p[i]-2*s->b[i]-s->b[i+1])/dx[i];
                s->d[i]=(s->b[i+1]+s->b[i]-2*p[i])/dx[i]/dx[i];

	}
	return s;
}

double aspline_eval(aspline* s, double z){
	if(z<s->x[0] || z>s->x[s->n-1]){fprintf(stderr,"Point passed to cspline_w_derives_eval is outside the x-range of the data and spline interpolation.\n");return 1;}
	int i = 0, j=s->n-1;
	while(j-i>1){int m = (i+j)/2; if(z>s->x[m]) i=m; else j=m;}
	double h= z - s->x[i];
	return s->y[i]+h*s->b[i]+h*h*s->c[i]+h*h*h*s->d[i];
}

double aspline_derive(aspline* s, double z){
	if(z<s->x[0] || z>s->x[s->n-1]){fprintf(stderr,"Point passed to cspline_w_derives_derive is outside the x-range of the data and spline interpolation.\n");return 1;}
	int i = 0, j=s->n-1;
        while(j-i>1){int m = (i+j)/2; if(z>s->x[m]) i=m; else j=m;}
	double h = z-s->x[i];
	return s->b[i] + 2.0*h*s->c[i] + 3.0*h*h*s->d[i];
}

double aspline_integrate(aspline* s, double z){
	if(z<s->x[0] || z>s->x[s->n-1]){fprintf(stderr,"Point passed to cspline_w_derives_integrate is outside the x-range of the data and spline interpolation.\n");return 1;}
        int i = 0, j=s->n-1;
        while(j-i>1){int m = (i+j)/2; if(z>s->x[m]) i=m; else j=m;}
	double value = 0.0;
	double dx = 0.0;
	for(int j = 0; j<i; ++j){
		dx = s->x[j+1]-s->x[j];
		value += dx*(s->y[j] + dx*(1.0/2*s->b[j] + dx*(1.0/3*s->c[j] + dx*1.0/4*s->d[j])));
	}
        double h = z - s->x[i];
	value += h*(s->y[i] + h*(1.0/2*s->b[i] + h*(1.0/3*s->c[i] + h*1.0/4*s->d[i])));
	return	value;
}


void aspline_free(aspline* s){
	free(s->x);free(s->y);free(s->b);free(s->c);free(s->d);free(s);}
