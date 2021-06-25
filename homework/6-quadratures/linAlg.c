#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX

double norm(double x){
	double result = sqrt(pow(x,2));
	return result;
}

double cdot(gsl_vector* A, gsl_vector* B){
	assert(A->size==B->size);
	double result=0;
	for(int i = 0;i<A->size;i++){
		double b = gsl_vector_get(A,i)*gsl_vector_get(B,i);
		result += b;
	}
	return result;
}

void show_matrix(gsl_matrix* A){
	int n = A->size1;
	int m = A->size2;
	fprintf(stderr,"\n");
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			double x = gsl_matrix_get(A,i,j);
				if(norm(x)<10e-10){
					fprintf(stderr,"%9i ",0);
				}
				else{
					fprintf(stderr,"%9.3g ",x);
				}
		}
		printf("\n");
	}
	printf("\n");
}

void show_vector(gsl_vector* V){
	int n = V->size;
	printf("\n");
	for(int i=0;i<n;i++){
		fprintf(stderr,"%g\n",gsl_vector_get(V,i));
	}
	printf("\n");
}
void GS_bak(gsl_matrix* R, gsl_vector* x){
	int m=R->size1;
	for(int i=m-1;i>=0;i--){
		double xi = gsl_vector_get(x,i);
		for(int j=i+1;j<m;j++)
			xi-=gsl_matrix_get(R,i,j)*gsl_vector_get(x,j);
		gsl_vector_set(x,i,xi/gsl_matrix_get(R,i,i));
	}
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);
	GS_bak(R,x);
}

void GS_decomp(gsl_matrix* A,gsl_matrix* R){
	assert(A->size2 == R->size1);
	int N = A->size1;
	int M = A->size2;
	gsl_vector* ai = gsl_vector_alloc(N);
	gsl_vector* aj = gsl_vector_alloc(N);
	for(int i=0;i<M;i++){
                gsl_matrix_get_col(ai,A,i);
		double Rii = gsl_blas_dnrm2(ai);
		gsl_matrix_set(R,i,i,Rii);
		gsl_vector_scale(ai,1/Rii);
		gsl_matrix_set_col(A,i,ai);
		for(int j=i+1;j<M;j++){
	  		gsl_matrix_get_col(aj,A,j);
			double Rij = cdot(ai,aj);
			gsl_blas_daxpy(-Rij,ai,aj);
			gsl_matrix_set(R,i,j,Rij);
			gsl_matrix_set_col(A,j,aj);
		}
	}
	gsl_vector_free(ai);
	gsl_vector_free(aj);
}

void GS_inv(gsl_matrix* A,gsl_matrix* Inv){
	int n = A->size2;
	assert(n==A->size1);
	gsl_matrix* I = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_matrix* Q = gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(Q,A);
	gsl_matrix_set_identity(I);
	gsl_vector* ei = gsl_vector_alloc(n);
	gsl_vector* x = gsl_vector_alloc(n);
	GS_decomp(Q,R);
	for(int i = 0;i<n;i++){
		gsl_matrix_get_col(ei,I,i);
		GS_solve(Q,R,ei,x);
		gsl_matrix_set_col(Inv,i,x);
		}
	gsl_matrix_free(I);
	gsl_matrix_free(R);
	gsl_matrix_free(Q);
	gsl_vector_free(ei);
	gsl_vector_free(x);
}

void lsfit(int m, double f(int i, double x), gsl_vector* x, gsl_vector* y,gsl_vector* dy,gsl_vector* c, gsl_matrix* S){
	int n = y->size;
	assert(y->size==x->size);
   	assert(c->size==m);
	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_matrix* R = gsl_matrix_alloc(m,m);
	gsl_vector* b = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		double xi = gsl_vector_get(x ,i);
		gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
		for(int k=0;k<m;k++){
			gsl_matrix_set(A,i,k,f(k,xi)/gsl_vector_get(dy,i));
		}
	}
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,S);
	GS_decomp(A,R);
	GS_solve(A,R,b,c);
	GS_inv(S,S);
	gsl_matrix_free(R);
	gsl_matrix_free(A);
	gsl_vector_free(b);
}
void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c= cos(theta), s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip = -s*gsl_matrix_get(A,i,q) + c*gsl_matrix_get(A,i,p);
		double new_aiq = c*gsl_matrix_get(A,i,q) + s*gsl_matrix_get(A,i,p);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta), s=sin(theta);
	for(int j=0;j<A->size2;j++){
		double new_apj=c*gsl_matrix_get(A,p,j) + s*gsl_matrix_get(A,q,j);
		double new_aqj=-s*gsl_matrix_get(A,p,j) + c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_apj);
		gsl_matrix_set(A,q,j,new_aqj);
	}
}

void jacobi_diag(gsl_matrix* A, gsl_matrix* V){
	int n=A->size1;
	int changed;
	do{
		changed = 0;
		for(int p=0;p<n-1;p++){
			for(int q=p+1;q<n;q++){
				double app=gsl_matrix_get(A,p,p);
				double apq=gsl_matrix_get(A,p,q);
				double aqq=gsl_matrix_get(A,q,q);
				double theta=0.5*atan2(2*apq,(aqq-app));
				double c=cos(theta),s=sin(theta);
				double new_aqq=s*s*app-2*s*c*apq+c*c*aqq;
				double new_app=c*c*app-2*s*c*apq+s*s*aqq;
				if(new_app!=app || new_aqq!=aqq){
					changed=1;
					timesJ(A,p,q,theta);
					Jtimes(A,p,q,-theta);
					timesJ(V,p,q,theta);
				}
			}
		}
	}while(changed!=0);
}
