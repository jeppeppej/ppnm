#include <stdio.h>
#include <limits.h>
#include <float.h>

int equal(double a, double b, double tau, double epsilon); //declaring function prototype

int main(){
	
	printf("INT_MAX = %i\n",INT_MAX);

	int i=1;
	while(i+1>i) {i++;}
	printf("in a while-loop, my max int = %i\n",i);

	int j=1;
	for(j; j+1>j; j++) {;}
	printf("in a for-loop, my max int = %i\n",j);

	int k=1;
	do{k++;}
	while(k+1>k);
	printf("in a do-while-loop, my max int = %i\n",k);

	printf("\nINT_MIN = %i\n",INT_MIN);

	int l=1;
	while(l-1<l) {l--;}
	printf("in a while-loop, my min int = %i\n",l);

	int m=1;
	for(m; m-1<m; m--) {;}
	printf("in a for-loop, my min int = %i\n",m);

	int n=1;
	do{n--;}
	while(n-1<n);
	printf("in a do-while-loop, my min int = %i\n",n);

	printf("\nFLT_EPSILON = %g\n",FLT_EPSILON);
	printf("DBL_EPSILON = %g\n",DBL_EPSILON);
	printf("LDBL_EPSILON = %Lg\n",LDBL_EPSILON);

	float x=1;
	while(1+x!=1){x/=2;}
	x*=2;
	printf("\nin a while-loop, my machine epsilon for floats = %g\n",x);

	double y=1;
	for(y; 1+y!=1; y/=2){;}
	y*=2;
	printf("in a for-loop, my machine epsilon for doubles = %g\n",y);

	long double z=1;
	do{z/=2;}
	while(1+z!=1);
	z*=2;
	printf("in a do-while-loop, my machine epsilon for long doubles = %Lg\n",z);



	int max=INT_MAX/2;
	float f=10;

	for(int i=1;i<=max;i++){
		f += 1.0f/i;
	}
	printf("sum_up_float = %g\n",f);
	f=10;
	for(int i=0;i<max;i++){
		f += 1.0f/(max-i);
	}
	printf("sum_down_float = %g\n",f);

	double d=10;

	for(int i=1;i<=max;i++){
		d += 1.0d/i;
	}
	printf("\nsum_up_double = %g\n",d);
	d=10;
	for(int i=0;i<max;i++){
		d += 1.0d/(max-i);
	}
	printf("sum_down_float = %g\n",d);

	double a = 12;
	double b = 60;
	double tau = 10;
	double epsilon = 10;

	printf("\nAre the numbers %g and %g equal with absolute precision tau=%g or relative precision epsilon=%g?\n", a, b, tau, epsilon);

	if(equal(a, b, tau, epsilon)==1)
		printf("Yes.\n");
	else
		printf("No.\n");
}
