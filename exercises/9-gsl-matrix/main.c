#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

int main(){
	printf("A: solve the following system of equations\n");
	printf("[ 6.13   -2.90   5.86 ] [x0]   [6.23]\n");
	printf("[ 8.08   -6.31  -3.89 ] [x1] = [5.37]\n");
	printf("[-4.36    1.00   0.19 ] [x2]   [2.29]\n");
	double a_data[] = {
		6.13, -2.90, 5.86,
		8.08, -6.31, -3.89,
		-4.36, 1.00, 0.19
	};
	double b_data[] = {6.23, 5.37, 2.29};
	gsl_matrix_view A = gsl_matrix_view_array(a_data, 3, 3);
	gsl_vector_view b = gsl_vector_view_array(b_data, 3);
	gsl_vector *x= gsl_vector_alloc(3);
	gsl_matrix *A_copy = gsl_matrix_alloc(3, 3);
	gsl_matrix_memcpy(A_copy, &A.matrix);

	gsl_linalg_HH_solve(&A.matrix, &b.vector, x);

	gsl_vector *c = gsl_vector_alloc(3);

	gsl_blas_dgemv(CblasNoTrans /*no transformation*/, 1 /*no scaling*/, 
A_copy/*the matrix*/, x/*the vector*/, 0 /*we do not have an extra term*/, c
/*matrix-vector product*/);//blas bruger mærkelige navne, double genereal matrix vector

	printf("The vector x:\n");
	gsl_vector_fprintf(stdout,x,"%.5g");

	printf("The matrix-vector product:\n");
	gsl_vector_fprintf(stdout,c,"%.5g");

	printf("It is the same as the intial vector.\n");

	return 0;
}


