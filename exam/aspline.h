#ifndef HAVE_ASPLINES_H
#define HAVE_ASPLINES_H
typedef struct {int n; double *x, *y, *b, *c, *d;} aspline;
aspline* aspline_alloc(int n, double* x, double* y, double* y_prime);
double aspline_eval(aspline* s, double z);
void aspline_free(aspline* s);
double aspline_derive(aspline* s, double z);
double aspline_integrate(aspline* s, double z);
#endif
