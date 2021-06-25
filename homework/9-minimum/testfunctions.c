#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<assert.h>

double Rosenbrock(gsl_vector* vec){
	double x=gsl_vector_get(vec,0);
	double y=gsl_vector_get(vec,1);
	double value=(1-x)*(1-x)+100*(y-x*x)*(y-x*x);
return value;
}

double Himmelblau(gsl_vector* vec){
	double x=gsl_vector_get(vec,0);
	double y=gsl_vector_get(vec,1);
	double value=(x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
return value;
}


double Breit_Wigner(double E,gsl_vector* vec){
	double m=gsl_vector_get(vec,0);
	double L=gsl_vector_get(vec,1);
	double A=gsl_vector_get(vec,2);
	double value=(double) A/((E-m)*(E-m)+L*L/4);
return value;
}


double deviation_function(gsl_vector* vec){
	double E[30], C[30], dC[30];

	FILE *my_file;
	my_file = fopen("data.txt","r");
	if(my_file == NULL)printf("Can't open file\n");
	else{
		for(int i=0;i<30;i++){
			assert(fscanf(my_file,"%lf",&E[i]) == 1);
			assert(fscanf(my_file,"%lf",&C[i]) == 1);
			assert(fscanf(my_file,"%lf",&dC[i]) == 1);
			}
	}
	fclose(my_file);



	double value=0;
	for(int i=0;i<30;i++){
		double s=(Breit_Wigner(E[i],vec)-C[i])*(Breit_Wigner(E[i],vec)-C[i])/(dC[i]*dC[i]);
		value += s;}

return value;
}
