#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

typedef struct {int n;double (*f)(double);gsl_vector * params;} ANN;

ANN* ANN_alloc(int n,double (*f)(double));

void ANN_free(ANN* network);

double ANN_response(double x,ANN* network);

void ANN_learn(ANN* network,
		gsl_vector* xs,
		gsl_vector* ys);

double ANN_eval_int(ANN* network,double F(double x),double x);

double ANN_eval_deriv(ANN* network,double fp(double x),double x);
