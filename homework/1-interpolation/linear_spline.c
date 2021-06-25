#include<stdio.h>
#include<math.h>
#include<assert.h>

double linterp(int n, double* x, double* y, double z) {

	assert(n > 1 && z >= x[0] && z <= x[n-1]);
	int i = 0;
    int j = n - 1;
	while(j - i > 1) {
        int mid = (i + j) / 2;
        if(z > x[mid]) { i = mid; }
        else{ j = mid; }
    }
	assert(x[i+1]>x[i]);

	double p = (y[i+1] - y[i]) / (x[i+1] - x[i]);

	return y[i] + p * ( z - x[i] );
}

double linterp_integ(int n, double x[], double y[], double z) {
	assert(n > 1 && z >= x[0] && z <= x[n-1]);

	int i = 0;
    int j = n - 1;
	while(j - i > 1) {
        int mid = (i + j) / 2;
        if(z > x[mid]) { i = mid; }
        else{ j = mid; }
    }

	double integral=0;
	for(int i=0; i < j; i++){
		double p = (y[i+1] - y[i]) / (x[i+1] - x[i]);
		double lhs = y[i] * x[i] + p * (x[i] * x[i] / 2 - x[i] * x[i]);
		if(z >= x[i+1]) {
		    double rhs = y[i] * x[i+1] + p * (x[i+1] * x[i+1] / 2 - x[i] * x[i+1]);
		    integral += rhs - lhs;
		}
		else {
		    double rhs = y[i] * z + p * (z * z / 2 - x[i] * z);
		    integral += rhs - lhs;
		}
	}
	return integral;
}
