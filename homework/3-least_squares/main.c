#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <assert.h>

void print_vector(gsl_vector *v);
void print_matrix(gsl_matrix *A);
void gen_rand_vector(gsl_vector *v);
void gen_rand_matrix(gsl_matrix *A);
int equals(double a, double b, double tau, double epsilon);
int vector_equals(gsl_vector *v, gsl_vector *u, double tau, double eps);
int matrix_equals(gsl_matrix *A, gsl_matrix *B, double tau, double eps);
int check_identity(gsl_matrix *A, double tau, double eps);
void backsub(gsl_matrix *U, gsl_vector *c);
void GS_decomp(gsl_matrix *A, gsl_matrix *R);
void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x);
void GS_inverse(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B);

void print_vector(gsl_vector *v){
    int n = v->size;
    for(int i = 0; i<n; i++){
        printf("%10.7f ", gsl_vector_get(v, i));
    }
    printf("\n");
}

void print_matrix(gsl_matrix *A){
    int n = A->size1;
    for(int i = 0; i<n; i++){
        gsl_vector_view v  = gsl_matrix_row(A, i);
        print_vector(&v.vector);
    }
}

void GS_decomp(gsl_matrix *A, gsl_matrix *R){
        int n = A->size1;
        int m = A->size2;
        assert(R->size1==R->size2 && R->size1 == m);
        assert(n>=m);

        for(int i = 0; i<m; i++){
                gsl_vector_view a_i = gsl_matrix_column(A, i);

                double a_i_norm = gsl_blas_dnrm2(&a_i.vector);
                gsl_matrix_set(R, i, i, a_i_norm);

                gsl_vector_scale(&a_i.vector, 1./a_i_norm);

                for(int j = i+1; j<m; j++){
                        gsl_vector_view a_j = gsl_matrix_column(A, j);

                        double qa = 0;
                        gsl_blas_ddot(&a_i.vector, &a_j.vector, &qa);
                        gsl_matrix_set(R, i, j, qa);

                        gsl_blas_daxpy(-qa, &a_i.vector, &a_j.vector);
                }
        }
}

void backsub(gsl_matrix *U, gsl_vector *c){
        assert(U->size1 == U->size2);
        int n = U->size1;
        for(int i = n-1; i>=0; i--){
                double s = gsl_vector_get(c, i);
                for(int k = i+1; k<n; k++){
                        s -= gsl_matrix_get(U, i, k)*gsl_vector_get(c, k);
                }
                gsl_vector_set(c, i, s/gsl_matrix_get(U,i,i));
        }
}

void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x){
        //Ax=(QR)x=b, QTQ=I
        //Rx=QTb
        //R is upper triangular -> use backsub
        gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
        backsub(R, x);
}

void GS_inverse(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B){
        assert(Q->size1 == Q->size2 && R->size1 == R->size2 && Q->size1 == R->size1);
        int n = Q->size1;
        assert(B->size1 == n && B->size2 == n);
        gsl_vector *e_i = gsl_vector_alloc(n);
        gsl_vector *x = gsl_vector_alloc(n);
        for(int i = 0; i<n; i++){ // Solve Ax_i = e_i where A⁻¹ = {x_i}, A = QR
                gsl_vector_set_basis(e_i, i);
                gsl_vector_set_zero(x);

                gsl_vector_set(e_i, i, 1);
                GS_solve(Q, R, e_i, x);
                gsl_matrix_set_col(B, i, x);
        }
        gsl_vector_free(e_i);
        gsl_vector_free(x);
}


double f0(double x){
	return 1.;
}
double f1(double x){
	return x;
}

void lsfit(double *x, double *y, double *yerr, double (**fs)(double), gsl_vector *v, gsl_matrix *M){
	gsl_matrix *A = gsl_matrix_alloc(8, 2);
	gsl_matrix *B = gsl_matrix_alloc(2, 2);
	gsl_vector *u = gsl_vector_alloc(8);

	//Translating the Python implementation

	for(int i = 0; i<8; i++){
		gsl_vector_set(u, i, y[i]/yerr[i]);
		for(int j = 0; j<2; j++){
			gsl_matrix_set(A, i, j, fs[j](x[i])/yerr[i]);
		}
	}

	//Now the fitting
	GS_decomp(A, B);
	GS_solve(A, B, u, v);

	gsl_matrix *B_inv = gsl_matrix_alloc(2, 2);
	gsl_matrix *I = gsl_matrix_alloc(2, 2);
	gsl_matrix_set_identity(I);
	GS_inverse(I, B, B_inv);

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, B_inv, B_inv, 0, M);

	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_matrix_free(I);
	gsl_vector_free(u);
}

int main(){
	double t[9] = {1,2,3,4,6,9,10,13,15};// the radioactivity data
	double A[9] = {117,100,88,72,53,29.5,25.2,15.2,11.1};

	gsl_vector *v = gsl_vector_alloc(2);
	gsl_matrix *M = gsl_matrix_alloc(2, 2);

	double (*fs[])(double) = {f0, f1};

	double logA[9];
	double logyerr[9];
	for(int i = 0; i<8; i++){
		logA[i] = log(A[i]);
		logyerr[i] = 1./20;
	}

	lsfit(t, logA, logyerr, fs, v, M);

	double a = exp(gsl_vector_get(v, 0));
	double a_err = a*sqrt(gsl_matrix_get(M, 0, 0));
	double g = -gsl_vector_get(v, 1);
	double g_err = sqrt(gsl_matrix_get(M, 1, 1));

	printf("\nThis gives the function: f(x) = %g*exp(-%g*x)\n", a, g);

	printf("\nCovariance matrix S\n");
	print_matrix(M);
	printf("\nThis gives a relative activity @ time 0: a = %g +- %g\n", a, a_err);
	printf("And gamma = %g +- %g\n", g, g_err);
	printf("Which gives a half-life of (%g +- %g) d\n", log(2)/g, log(2)/(g*g)*g_err);
	printf("Half-life for 224Ra = 3.63 d\nWhich is not terribly far away from our fitted value.\nHowever it still outside the estimated uncertainty.\n");
	gsl_vector_free(v);
	gsl_matrix_free(M);
	//free(t);
	//free(A);

return 0;
}
