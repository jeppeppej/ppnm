#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "jacobi_diag.h"


void print_matrix(gsl_matrix* M)
{
    for(int i = 0; i < M->size1; i++)
    {
        for(int j = 0; j < M->size2; j++)
        {
            printf("%10g\t", gsl_matrix_get(M, i, j));
        }
        printf("\n");
    }
}


void print_vector(gsl_vector* v)
{
    for(int i = 0; i < v->size; i++) printf("%g\n", gsl_vector_get(v, i));
}

gsl_matrix* random_matrix(int m, int n)
{
    gsl_matrix* M = gsl_matrix_alloc(m, n);
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            gsl_matrix_set(M, i, j, (double)rand()*10.0/RAND_MAX);
        }
    }
    return M;
}

gsl_matrix* random_symm_matrix(int n)
{
    gsl_matrix* M = gsl_matrix_alloc(n, n);
    for(int i = 0; i < n; i++)
    {
        for(int j = i; j < n; j++)
        {
            double Mij = (double)rand()*10.0/RAND_MAX;
            gsl_matrix_set(M, i, j, Mij);
            gsl_matrix_set(M, j, i, Mij);
        }
    }
    return M;
}

void timesJ(gsl_matrix* A, int p, int q, double theta)
{
	double c = cos(theta), s = sin(theta);
	for(int i = 0; i < A->size1; i++)
    {
		double new_aip = c*gsl_matrix_get(A, i, p) - s*gsl_matrix_get(A, i, q);
		double new_aiq = s*gsl_matrix_get(A, i, p) + c*gsl_matrix_get(A, i, q);
		gsl_matrix_set(A, i, p, new_aip);
		gsl_matrix_set(A, i, q, new_aiq);
	}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta)
{
	double c = cos(theta), s = sin(theta);
	for(int j = 0; j < A->size2; j++)
    {
		double new_apj = c*gsl_matrix_get(A, p, j) + s*gsl_matrix_get(A, q, j);
		double new_aqj = -s*gsl_matrix_get(A, p, j) + c*gsl_matrix_get(A, q, j);
		gsl_matrix_set(A, p, j, new_apj);
		gsl_matrix_set(A, q, j, new_aqj);
	}
}

void jacobi_diag(gsl_matrix* A, gsl_matrix* V)
{
    int changed,  n = A->size1;
    gsl_matrix_set_identity(V);
    do
    {
        changed = 0;
	    for(int p = 0; p < n-1; p++)
        {
            for(int q = p+1; q < n; q++)
            {
                double apq = gsl_matrix_get(A, p, q);
                double app = gsl_matrix_get(A, p, p);
                double aqq = gsl_matrix_get(A, q, q);
                double theta = 0.5*atan2(2*apq, aqq-app);
                double c = cos(theta), s = sin(theta);
                double new_app = c*c*app - 2*s*c*apq + s*s*aqq;
                double new_aqq = s*s*app + 2*s*c*apq + c*c*aqq;
                if(new_app != app || new_aqq != aqq)
                {
                    changed = 1;
                    timesJ(A, p, q, theta);
                    Jtimes(A, p, q, -theta);
                    timesJ(V, p, q, theta);
                }
            }

	    }
    }
    while(changed != 0);
}
