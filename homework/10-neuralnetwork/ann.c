#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"ann.h"

void ngradient(double F(gsl_vector * x),
		gsl_vector * x,
		gsl_vector * grad);

int qnewton(double F(gsl_vector * x),
		gsl_vector * x,
		double eps);

void downhill_simplex(double f(gsl_vector* x),
		gsl_matrix* simplex,
		double tolerance,
		gsl_vector* minimum);

ANN* ANN_alloc(int n,double (*f)(double)){
	ANN* network = malloc(sizeof(ANN));
	network->n=n;
	network->f=f;
	network->params = gsl_vector_alloc(3*n);
	return network;
}

void ANN_free(ANN* network){
	gsl_vector_free(network->params);
	free(network);
}

double ANN_response(double x,ANN* network){
	double y=0;
	for(int i=0;i<network->n;i++){
		double a = gsl_vector_get(network->params,3*i+0);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		y+=network->f((x-a)/b)*w;
	}
	return y;
}

void ANN_learn(ANN* network,
		gsl_vector* xs,
		gsl_vector* ys){
	gsl_vector* p = gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);

	double cost(gsl_vector * p){
		gsl_vector_memcpy(network->params,p);
		double s = 0;
		for(int k=0;k<xs->size;k++){
		double Fk = ANN_response(gsl_vector_get(xs,k),network);
		double yk = gsl_vector_get(ys,k);
		s+=(Fk-yk)*(Fk-yk);
		}
		return s;
	}

	double tolerance = 0.001;
	qnewton(cost,p,tolerance);

	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}

double ANN_eval_int(ANN* network,double F(double x),double x){
	double sum = 0;
	for(int i=0;i<network->n;i++){
		double a = gsl_vector_get(network->params,3*i+0);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		sum+=w*b*F((x-a)/b);
	}
	return sum;
}

double ANN_eval_deriv(ANN* network,double fp(double x),double x){
	double sum = 0;
	for(int i=0;i<network->n;i++){
		double a = gsl_vector_get(network->params,3*i+0);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		sum+=w*fp((x-a)/b)/b;
	}
	return sum;
}
