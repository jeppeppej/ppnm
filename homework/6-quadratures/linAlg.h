#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX


double norm(double x);
double cdot(gsl_vector* A, gsl_vector* B);
void show_matrix(gsl_matrix* A);
void show_vector(gsl_vector* V);
void GS_bak(gsl_matrix* R, gsl_vector* x);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void GS_decomp(gsl_matrix* A,gsl_matrix* R);
void GS_inv(gsl_matrix* A,gsl_matrix* Inv);
void lsfit(int m, double f(int i, double x), gsl_vector* x, gsl_vector* y,gsl_vector* dy,gsl_vector* c, gsl_matrix* S);
void timesJ(gsl_matrix* A, int p, int q, double theta);
void Jtimes(gsl_matrix* A, int p, int q, double theta);
void jacobi_diag(gsl_matrix* A, gsl_matrix* V);
