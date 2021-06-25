#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include "jacobi_diag.h"

#define M_PI 3.14159265358979323846

int main()
{
    //PART A
    gsl_matrix* A = gsl_matrix_alloc(2, 2);
    gsl_matrix* A_copy = gsl_matrix_alloc(2, 2);
    gsl_matrix* V = gsl_matrix_alloc(2, 2);

    gsl_matrix_set(A, 0, 0, 0.0);
    gsl_matrix_set(A, 0, 1, 2.0);
    gsl_matrix_set(A, 1, 0, 2.0);
    gsl_matrix_set(A, 1, 1, 3.0);
    gsl_matrix_memcpy(A_copy, A);

    //Making matrices
    gsl_matrix* VTV = gsl_matrix_alloc(2, 2);
    gsl_matrix* VVT = gsl_matrix_alloc(2, 2);
    gsl_matrix* VTAV = gsl_matrix_alloc(2, 2);
    gsl_matrix* AV = gsl_matrix_alloc(2, 2);
    gsl_matrix* VDVT = gsl_matrix_alloc(2, 2);
    gsl_matrix* DVT = gsl_matrix_alloc(2, 2);

    printf("PART A\n");
    printf("A matrix: \n");
    print_matrix(A);

    //Diagonalize
    jacobi_diag(A, V);

    printf("Diagonal matrix: \n");
    print_matrix(A);
    printf("The orthogonalized eigenvectors: \n");
    print_matrix(V);

    //And a final check
    printf("Now some checks.\nV^T*V = \n");
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, V, 0, VTV);
    print_matrix(VTV);

    printf("V*V^T = \n");
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, V, V, 0, VVT);
    print_matrix(VVT);

    printf("V^T*A*V = \n");
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A_copy, V, 0, AV);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, AV, 0, VTAV);
    print_matrix(VTAV);

    printf("V*D*V^T = \n");
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, A, V, 0, DVT);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, V, DVT, 0, VDVT);
    print_matrix(VDVT);

    //PART B
    printf("PART B\n");
    //Hamiltonian
    int n = 100;
    double s = 1.0/(n + 1);
    gsl_matrix* H = gsl_matrix_alloc(n, n);
    for(int i = 0; i < n-1; i++)
    {
        gsl_matrix_set(H, i, i, -2);
        gsl_matrix_set(H, i, i + 1, 1);
        gsl_matrix_set(H, i+1, i, 1);
    }
    gsl_matrix_set(H, n-1, n-1, -2);
    gsl_matrix_scale(H, -1/s/s);

    //Diagonalize it
    gsl_matrix* eigenvectors = gsl_matrix_alloc(n, n);
    jacobi_diag(H, eigenvectors);

    printf("Reulting energies: \n");
    printf("Nr. \t Calculated \t Exact\n");
    for (int k=0; k < n/3; k++)
    {
        double exact = M_PI*M_PI*(k + 1)*(k + 1);
        double calculated = gsl_matrix_get(H, k, k);
        printf("%i \t %g \t %g\n", k, calculated, exact);
    }

    printf("A plot has also been produced, showing the  wavefunctions\n");
    //Data to plot
    FILE* data = fopen("data.txt", "w");
    for(int i = 0; i < n; i++)
    {
        double Vi1 = gsl_matrix_get(eigenvectors, i, 1);
        double Vi2 = gsl_matrix_get(eigenvectors, i, 2);
        double Vi3 = gsl_matrix_get(eigenvectors, i, 3);

        fprintf(data, "%10g %10g %10g %10g\n", (i + 1.0)/(n + 1), Vi1, Vi2, Vi3);
    }

    gsl_matrix_free(H);
    gsl_matrix_free(eigenvectors);

    return 0;
}
