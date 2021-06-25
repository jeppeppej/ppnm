#include<assert.h>
#include<gsl/gsl_blas.h>
#include<math.h>

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

void random_matrix(gsl_matrix *A){
	int n = A->size1;
	int m = A->size2;
	for(int i = 0; i<n; i++){
		for(int j = 0; j<m; j++){
			gsl_matrix_set(A, i, j, 1.0*rand()/RAND_MAX);
		}
	}
}

void random_vector(gsl_vector *v){
	for(int i = 0; i<v->size; i++){
		gsl_vector_set(v, i, 1.0*rand()/RAND_MAX);
	}
}

int main(){
//	srand(time(NULL));
	double n = 8;
	double m = 4;

	gsl_matrix *A = gsl_matrix_alloc(n,m);
	random_matrix(A);

	gsl_matrix *Q = gsl_matrix_alloc(n,m);
	gsl_matrix_memcpy(Q, A);

	gsl_matrix *R = gsl_matrix_alloc(m,m);
	GS_decomp(Q,R);

	//Now check results
	gsl_matrix *QR = gsl_matrix_alloc(n,m);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Q, R, 0, QR);
	gsl_matrix *QTQ = gsl_matrix_alloc(m,m);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, Q, Q, 0, QTQ);

	printf("Random A: \n");
	print_matrix(A);
	printf("\n Random Q: \n");
	print_matrix(Q);
	printf("\n Upper triangular R: \n");
	print_matrix(R);
	printf("\n QR=A: \n");
	print_matrix(QR);
	printf("\n Identity=QTQ: \n");
	print_matrix(QTQ);

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(Q);
	gsl_matrix_free(QR);
	gsl_matrix_free(QTQ);

	A = gsl_matrix_alloc(n,n);
	random_matrix(A);
	gsl_vector *b = gsl_vector_alloc(n);
	random_vector(b);
	Q = gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(Q,A);
	R =  gsl_matrix_alloc(n,n);
	GS_decomp(Q,R);
	gsl_vector *x = gsl_vector_alloc(n);
	GS_solve(Q,R,b,x);
	gsl_vector *Ax = gsl_vector_alloc(n);
	gsl_blas_dgemv(CblasNoTrans, 1, A, x, 0, Ax);

	printf("b: \n");
	print_vector(b);
	printf("\n Ax=b: \n");
	print_vector(Ax);

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(Q);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(Ax);

	return 0;
}
