#include <stdio.h>
#include <math.h>
#include <complex.h>

int main(){
	double a = tgamma(5);
	printf("Gamma(5)=%g\n", a);

	double b = j1(0.5);
	printf("J_1(0.5)=%g\n", b);

	double complex c = csqrt(-2);
	printf("sqrt(-2)=%g+%gi\n", creal(c), cimag(c));

	double complex d = cexpl(I*M_PI);
	printf("e^(i*pi)=%g+%gi\n", creal(d), cimag(d));

	double complex e = cexpl(I);
	printf("e^(i)=%g+%gi\n", creal(e), cimag(e));

	double complex f = cpow(I, M_E);
	printf("i^(e)=%g+%gi\n", creal(f), cimag(f));

	double complex g = cpow(I, I);
	printf("i^(i)=%g+%gi\n", creal(g), cimag(g));


	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;
	printf("\nHow many significant digits can the types hold?\n%.25g\n%.25lg\n%.25Lg\n", x_float, x_double, x_long_double);
return 0;
}
