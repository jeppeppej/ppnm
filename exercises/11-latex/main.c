#include<stdio.h>
#include<math.h>
double myexp(double);
int main(){
	for(double x=1./16-100;x<=100;x+=1./16)
		printf("%g %g %g\n",x,myexp(x),exp(x));
return 0;
}
