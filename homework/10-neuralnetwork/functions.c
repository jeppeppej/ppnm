#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<stdio.h>
#include<float.h>
#include<stdlib.h>
#define delx  sqrt(DBL_EPSILON)

void vector_print(const char* s,gsl_vector * v){
	printf("%s ( ",s);
	int n=v->size;
	for(int i=0;i<n;i++){
		printf("%.3g ",gsl_vector_get(v,i));
	}
	printf(")\n");
}
void matrix_print(const char* s,gsl_matrix * M){
	printf("%s\n",s);
	int n=M->size1,m=M->size2;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
		printf("%.3g ",gsl_matrix_get(M,i,j));
		}
		printf("\n");
	}
}

void ngradient(double F(gsl_vector * x),
		gsl_vector * x,
		gsl_vector * grad)
		{
	int dim = x->size;
	for(int i = 0;i<dim;i++){
		double f_of_x = F(x);
		double dx,xi = gsl_vector_get(x,i);
		if(fabs(xi)<delx) dx = delx;
		else dx = fabs(xi)*delx;
		gsl_vector_set(x,i,xi+dx);
		double f_of_xpdx = F(x);
		double dFdxi = (f_of_xpdx - f_of_x) / dx;
		gsl_vector_set(grad,i,dFdxi);
		gsl_vector_set(x,i,xi);
	}
}

int qnewton(double F(gsl_vector * x),
		gsl_vector * x,
		double tol)
		{
int dim = x->size,nsteps=0;
gsl_matrix * B = gsl_matrix_alloc(dim,dim);
gsl_vector * deltax   = gsl_vector_alloc(dim);
gsl_vector * gradient = gsl_vector_alloc(dim);
gsl_vector * z        = gsl_vector_alloc(dim);
gsl_vector * u        = gsl_vector_alloc(dim);

gsl_matrix_set_identity(B);

while(nsteps<1000){
	nsteps++;
ngradient(F,x,gradient);
gsl_blas_dgemv(CblasNoTrans,-1,B,gradient,0,deltax);

double lambda = 1,lambda_min = DBL_EPSILON;
gsl_vector_memcpy(z,x);
gsl_vector_add(z,deltax);
double fz = F(z),fx = F(x);
double sTg;
gsl_blas_ddot(deltax,gradient,&sTg);

while(fx+0.01*sTg<=fz){
	if(lambda<lambda_min) {
		gsl_matrix_set_identity(B); break;
		}
	lambda *= 0.5;
	gsl_vector_scale(deltax,0.5);
	gsl_vector_memcpy(z,x);
	gsl_vector_add(z,deltax);
	fz = F(z);
	gsl_blas_ddot(deltax,gradient,&sTg);
}

gsl_vector_memcpy(x,z);
ngradient(F,x,z);
if(gsl_blas_dnrm2(z)<tol){
	break;
}

gsl_vector_sub(z,gradient);
gsl_vector_memcpy(u,deltax);
gsl_blas_dgemv(CblasNoTrans,-1,B,z,1,u);
double sTy,uTy;
gsl_blas_ddot(deltax,z,&sTy);
if(fabs(sTy)>1e-6){
	gsl_blas_ddot(u,z,&uTy);
	double gamma=uTy/2/sTy;
	gsl_blas_daxpy(-gamma,deltax,u);
	gsl_blas_dger(1./sTy,u,deltax,B);
	gsl_blas_dger(1./sTy,deltax,u,B);
}
}


gsl_matrix_free(B);
gsl_vector_free(deltax);
gsl_vector_free(gradient);
gsl_vector_free(z);
gsl_vector_free(u);
return nsteps;
}






void vector_print(const char* s,gsl_vector * v);

void reflection(gsl_vector* hi, gsl_vector* ce){
	gsl_vector_scale(hi,-1);
	gsl_blas_daxpy(2,ce,hi);
}

void expansion(gsl_vector * hi, gsl_vector* ce){
	gsl_vector_scale(hi,-2);
	gsl_blas_daxpy(2,ce,hi);
}

void contraction(gsl_vector * hi,gsl_vector * ce){
	gsl_vector_scale(hi,0.5);
	gsl_blas_daxpy(0.5,ce,hi);
}

void reduction(gsl_matrix * simplex,
	       gsl_vector * lo,
	       gsl_vector * x,
	       int kmin,
	       int n){

	for(int i=0;i<n+1;i++){
		if(i!=kmin){
		gsl_matrix_get_col(x,simplex,i);
		gsl_vector_add(x,lo);
		gsl_vector_scale(x,0.5);
		gsl_matrix_set_col(simplex,i,x);
		}
	}
}

double size_of_simplex(gsl_matrix * S){
	int n=S->size1;

	double s = 0,dist;
for(int k=0;k<n;k++){
	for(int i=k+1;i<n+1;i++){
		dist = 0;
		for(int l=0;l<n;l++){
		double pk = gsl_matrix_get(S,l,k);
		double pi = gsl_matrix_get(S,l,i);
		dist+=(pk-pi)*(pk-pi);
		}
		if(dist>s*s) s=sqrt(dist);
	}
}
return s;
}

void downhill_simplex(double f(gsl_vector * x),
		gsl_matrix * simplex,
		double tolerance,
		gsl_vector* minimum){
	int n = simplex->size1;

	gsl_vector* x  = gsl_vector_alloc(n);
	gsl_vector* hi = gsl_vector_alloc(n);
	gsl_vector* lo = gsl_vector_alloc(n);
	gsl_vector* ct = gsl_vector_alloc(n);
	gsl_vector* reflected  = gsl_vector_alloc(n);
	gsl_vector* expanded   = gsl_vector_alloc(n);
	gsl_vector* contracted = gsl_vector_alloc(n);
	gsl_vector* reduced    = gsl_vector_alloc(n);
	gsl_matrix_get_col(hi,simplex,0);
	       	gsl_vector_memcpy(lo,hi);
	double fmax = f(hi),newval,fmin = f(lo);
	int kmax = 0,kmin=0;
for(int k=1;k<n+1;k++){
	gsl_matrix_get_col(x,simplex,k);
	newval = f(x);
       	if(newval>fmax){
	       	gsl_vector_memcpy(hi,x);
		fmax = newval;kmax=k;
	}
	if(newval<fmin){
		gsl_vector_memcpy(lo,x);
		fmin = newval;kmin = k;
	}
}
for(int k=0;k<n+1;k++){
	if(k!=kmax){
	gsl_matrix_get_col(x,simplex,k);
	gsl_vector_add(ct,x);
	}
}
gsl_vector_scale(ct,1./n);

double size = size_of_simplex(simplex);

while(size>tolerance){
gsl_vector_memcpy(reflected,hi);
reflection(reflected,ct);
if(f(reflected)<f(lo)){
	gsl_vector_memcpy(expanded,hi);
	expansion(expanded,ct);
	if(f(expanded)<f(reflected)){
		gsl_matrix_set_col(simplex,kmax,expanded);
	} else gsl_matrix_set_col(simplex,kmax,reflected);
}
else{
	if(f(reflected)<f(hi)){
		gsl_matrix_set_col(simplex,kmax,reflected);
	}
	else{
		gsl_vector_memcpy(contracted,hi);
		contraction(contracted,ct);
		if(f(contracted)<f(hi)){
			gsl_matrix_set_col(simplex,kmax,contracted);
		}
		else {reduction(simplex,lo,x,kmin,n);
		}
	}
}
gsl_matrix_get_col(hi,simplex,0);
	gsl_vector_memcpy(lo,hi);
fmax = f(hi);
fmin = f(lo);

kmax = 0,kmin=0;
for(int k=1;k<n+1;k++){
	gsl_matrix_get_col(x,simplex,k);
	newval = f(x);
       	if(newval>fmax){
	       	gsl_vector_memcpy(hi,x);
		fmax = newval;kmax=k;
	}
	if(newval<fmin){
		gsl_vector_memcpy(lo,x);
		fmin = newval;kmin = k;
	}
}
gsl_vector_set_zero(ct);
for(int k=0;k<n+1;k++){
	if(k!=kmax){
	gsl_matrix_get_col(x,simplex,k);
	gsl_vector_add(ct,x);
	}
}
gsl_vector_scale(ct,1./n);

size = size_of_simplex(simplex);
}
gsl_vector_memcpy(minimum,lo);
	gsl_vector_free( x);
	gsl_vector_free(hi);
	gsl_vector_free(lo);
	gsl_vector_free(ct);
	gsl_vector_free( reflected);
	gsl_vector_free(  expanded);
	gsl_vector_free(contracted);
	gsl_vector_free(   reduced);
}
