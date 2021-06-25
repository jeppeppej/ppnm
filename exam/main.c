#include<stdio.h>
#include<math.h>
#include"aspline.h"
#include"cspline.h"

//My chosen test function
double test_function(double x){
	return 2.0/(1.0+exp(-8*(x)))-1;
}


//and its derivative
double test_function_derivative(double x){
	return 5*exp(-5*(x-5))/pow(1.0+exp(-5*(x-5)),2);
}

int main(){
	int points = 6;
	double start = -4.0, stop = 4.0;
	int plotpoints = 1025;
	double step = (stop-start)/(points-1.0);
	double x[points], y[points], y_prime[points];
	for(int i=0; i<points; ++i){
		x[i] = start + i*step;
		y[i] = test_function(x[i]);
		y_prime[i] = test_function_derivative(x[i]);
	}
	aspline* aspline = aspline_alloc(points,x,y,y_prime);
	cspline* cspline = cspline_alloc(points,x,y);
	FILE* data = fopen("data.txt","w");
	for(int i=0; i<points; ++i){
		fprintf(data,"%g %g\n",x[i],y[i]);
	}
	fprintf(data,"\n\n\n");
	double x_plot;
	double plotStep = (stop-start)/(plotpoints-1.0);
	for(int i=0; i<plotpoints;++i){
		x_plot = start+i*plotStep;
		fprintf(data,"%g %g %g %g\n",x_plot, aspline_eval(aspline,x_plot), test_function(x_plot),cspline_eval(cspline,x_plot));
	}
	fclose(data);

	aspline_free(aspline);
	cspline_free(cspline);

	printf("Multiple plots comparing this implementation of both a cubic sub-spline and a quartic sub-spline have been produced and included in the exam_report.pdf.\n");
	return 0;
}
