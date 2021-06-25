#ifndef HAVE_CSPLINE_H
#define HAVE_CSPLINE_H
typedef struct {int n; double *x, *y, *b, *c, *d;} cspline;
cspline* cspline_alloc(int n, double *x, double *y); 
double cspline_eval(cspline* s, double z);
double cspline_integrate(cspline* s, double z);
double cspline_derive(cspline* s, double z);
void cspline_free(cspline* s);
int binary_search(int n, double * x, double z);
#endif
