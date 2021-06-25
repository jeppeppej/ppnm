#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));
	printf("\n");
}

int qnewton(double F(gsl_vector*), gsl_vector* x, double eps);

double fun(gsl_vector* x_vec){
	double x=gsl_vector_get(x_vec,0);
	double value=(x*x+2*x-1);
return value;
}

double Rosenbrock(gsl_vector*);
double Himmelblau(gsl_vector*);
double Breit_Wigner(double,gsl_vector*);
double deviation_function(gsl_vector*);


int main(){
int steps,stepsR,stepsH,stepsB;
gsl_vector* x=gsl_vector_calloc(2);

steps=qnewton(fun,x,0.0005);

printf("A and B of minimization");
vector_print("A chosen test function f(x,y)=x^2+2x-1 has the minimum",x);
printf("which fits the expected. Steps used = %i\n\n",steps);


//Rosenbrock
gsl_vector* R=gsl_vector_calloc(2);
gsl_vector_set(R,0,1);
gsl_vector_set(R,1,0);

stepsR=qnewton(Rosenbrock,R,0);
vector_print("Moving on to Rosenbrock's valley funktionen.\n Minimum is:",R);
printf("The expected value is (1,1) for a=1, b=100, so this is good enough. Initial guess was (x,y)=(1,0) and number of steps=%i\n\n",stepsR);

//Himmelblau
gsl_vector* H=gsl_vector_calloc(2);
gsl_vector_set(H,0,3);gsl_vector_set(H,0,0);
stepsH=qnewton(Himmelblau,H,0.0001);
vector_print("Himmelblau's minimum:",H);
printf("Which exactly is one of them. Initial guess: (x,y)=(3,0) - steps: %i \n\n",stepsH);
printf("On to task B.\n\n");
//Berit-Wigner
gsl_vector* B=gsl_vector_calloc(3);
gsl_vector_set(B,0,122);gsl_vector_set(B,1,0.02);gsl_vector_set(B,2,0.001);
vector_print("Fitting to Berit Wigner, with initial guess:",B);
stepsB=qnewton(deviation_function,B,0.0000001);
vector_print("The found  minimum on the form (mass,width,A) is",B);
printf("Found in %i steps\n",stepsB);
printf("CERN's value = 125.3(6), so not far off.\n");
printf("A fit and CERN's data is plotted in the png-file.\n\n");

FILE *plot;
plot = fopen("plot.txt","w");
for(double e=95;e<170;e++)fprintf(plot,"%6g %6g\n",e,Breit_Wigner(e,B));
fclose(plot);

return 0;
}
