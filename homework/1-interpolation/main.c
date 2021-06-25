#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_interp.h>

// linear header
double linterp(int n, double* x, double* y, double z);
double linterp_integ(int n, double x[], double y[], double z);

//quadratic header
typedef struct {int n; double *x,*y,*b,*c;} qspline;
qspline* qspline_alloc(int n, double* x, double* y);
double qspline_eval(qspline *s, double z);
void qspline_free(qspline *s);
double qspline_integ(qspline *s, double z);
double qspline_deriv(qspline *s, double z);

//cubic header
typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;
cspline* cspline_alloc(int n, double *x, double *y);
void cspline_free(cspline *s);
double cspline_eval(cspline *s, double z);
double cspline_integ(cspline *s, double z);
double cspline_deriv(cspline *s, double z);


// some points
int n = 10;
double x[] = {0, 1, 2, 3,  4, 5,  6, 7, 8, 9, 10};
double y[] = {0, 1, 2, 4, 10, 6, 10, 4, 2, 1, 0};


void linear_spline() {
	// Write points to a file
	FILE* points = fopen("points.txt", "w");
	for(int i = 0; i < n; i++) {
		fprintf(points, "%10.4f %10.4f\n", x[i], y[i]);
	}

	gsl_interp* gsl_linterp = gsl_interp_alloc(gsl_interp_linear, n);
	gsl_interp_init(gsl_linterp, x, y, n);

	int N = 100;
	double stepsize = 9.0 / N;
	double z = 0;
	for(int i = 0; i <= N; i++) {
		// my linear interpolant
		double yi = linterp(n, x, y, z);
		double yi_integ = linterp_integ(n, x, y, z);

		// GSL's linear interpolant
		double yi_gsl = gsl_interp_eval(gsl_linterp, x, y, z, NULL);
		double yi_gsl_integ = gsl_interp_eval_integ(gsl_linterp, x, y, x[0], z, NULL);

		printf("%10.4f %10.4f %10.4f %10.4f %10.4f\n", z, yi, yi_integ, yi_gsl, yi_gsl_integ);
		z = z+stepsize;
	}

	gsl_interp_free(gsl_linterp);

}

int main() {
	linear_spline();

	printf("\n\n\n");
	qspline *qsp = qspline_alloc(n, x, y);
	int N = 100;
	double zmin = x[0];
	double zmax = x[n-1];
	double dz = (zmax - zmin)/N;
	double z = zmin;
	for(int i = 0; i<N; i++) {
		double qs = qspline_eval(qsp, z);
		double qsi = qspline_integ(qsp, z);
		double qsd = qspline_deriv(qsp, z);
		printf("%10.4f %10.4f %10.4f %10.4f\n", z, qs, qsi, qsd);
		z += dz;
	}
	qspline_free(qsp);


	printf("\n\n\n");
	cspline* csp = cspline_alloc(n, x, y);

	z = zmin;
	for(int i = 0; i<N; i++) {
		double cs = cspline_eval(csp, z);
		double csi = cspline_integ(csp, z);
		double csd = cspline_deriv(csp, z);
		printf("%10.4f %10.4f %10.4f %10.4f\n", z, cs, csi, csd);
		z += dz;
	}
	cspline_free(csp);

	return 0;
}
