#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_sf_gamma.h>

double myerf(double x){
/// single precision error function (Abramowitz and Stegun, from Wikipedia)
	if(x<0) return -myerf(-x);
	double a[]={0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429};
	double t=1/(1+0.3275911*x);
	double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));/* the right thing */
	return 1-sum*exp(-x*x);
}

double mygamma(double x){
///single precision gamma function (Gergo Nemes, from Wikipedia)
	if(x<0)return M_PI/sin(M_PI*x)/mygamma(1-x);
	if(x<9)return mygamma(x+1)/x;
	double lnmygamma=x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
	return exp(lnmygamma);
}


int main(){
	FILE *errorStream = fopen("error.txt","w");
	double xmin=-2,xmax=2;
	for(double x=xmin;x<=xmax;x+=1.0/8){
		fprintf(errorStream,"%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),myerf(x));}
	fclose(errorStream);
	FILE *gammaStream = fopen("gamma.txt","w");
	xmin=0.1, xmax=6;
	for(double x=xmin;x<=xmax;x+=1.0/8){
		fprintf(gammaStream,"%10g %10g %10g %10g\n",x,tgamma(x),gsl_sf_gamma(x),mygamma(x));}
	fclose(gammaStream);
	return 0;
}
