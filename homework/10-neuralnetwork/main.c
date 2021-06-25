#include"ann.h"

double activation(double x){
	return x*exp(-x*x);
}
double activation_int(double x){
	return -0.5*exp(-x*x);
}
double activation_deriv(double x){
	return exp(-x*x)-2*x*x*exp(-x*x);
}

double some_function(double x){
	return sin(5*x-1)*exp(-x*x);
}
double some_function_deriv(double x){
	return -exp(-x*x)*(5*sin(5*x-1)+2*x*cos(5*x-1));
}
double some_function_int(double x){
	return cos(5*x-1)*exp(-x*x);
}

int main(){
	int n = 4;
	ANN* my_network = ANN_alloc(n,activation);

	int N = 30;
	gsl_vector* xs = gsl_vector_alloc(N);
	gsl_vector* ys = gsl_vector_alloc(N);
	double a=-1,b=1,delx=(b-a)/(N-1),x;
	x = a;
	for(int i=0;i<N;i++){
		gsl_vector_set(xs,i,x);
		gsl_vector_set(ys,i,some_function(x));
		x+=delx;
	}

	for(int i=0;i<my_network->n;i++){
	gsl_vector_set(my_network->params,3*i+0,a+i*(b-a)/(my_network->n-1));
	gsl_vector_set(my_network->params,3*i+1,0.5);
	gsl_vector_set(my_network->params,3*i+2,0.5);
	}

	ANN_learn(my_network,xs,ys);

	for(int i=0;i<N;i++){
		double xval = gsl_vector_get(xs,i);
		double yval = gsl_vector_get(ys,i);
		printf("%g %g\n",xval,yval);
	}

	printf("\n\n");

	double xval = a;N=100;
	for(int i=0;i<N;i++){
		double Fval  = ANN_response(xval,my_network);
		double integ = ANN_eval_int(my_network,activation_int,xval);
		double deriv = ANN_eval_deriv(my_network,activation_deriv,xval);
		double trued = some_function_deriv(xval);
		printf("%g %g %g %g %g\n",xval,Fval,integ,deriv,trued);
		xval+=(b-a)/(N-1);
	}

	ANN_free(my_network);
}
