#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<assert.h>
#define RND ((double)rand()/RAND_MAX)

double corput(int n, int base){
double q=0, bk=(double)1/base ;
while (n>0){q+=(n%base)*bk; n/=base; bk/=base ; }
return q;
}

void halton(int n,gsl_vector *x){
	int base[]={2, 3 , 5 , 7 , 11 , 13 , 17 , 19 , 23 , 29 , 31 , 37 , 41 , 43 , 47 , 53 , 59 , 61 , 67} ;
	int maxd=sizeof(base)/sizeof(int);assert(x->size <= maxd);
	for(int i=0;i<x->size;i++) gsl_vector_set(x,i,corput(n,base[i]));
}

void plainmcA(double f(gsl_vector* x),gsl_vector* a,gsl_vector* b,int N,gsl_vector * results){
        double V=1; for(int i=0;i<a->size;i++)V*=gsl_vector_get(b,i)-gsl_vector_get(a,i);

		double sum=0,sum2=0;

		gsl_vector * x=gsl_vector_alloc(a->size);

		for(int i=0;i<N;i++){
                for(int i=0;i<x->size;i++)gsl_vector_set(x,i,gsl_vector_get(a,i)+RND*(gsl_vector_get(b,i)-gsl_vector_get(a,i)));
                double fx=f(x); sum+=fx; sum2+=fx*fx;
                }

		double mean=sum/N, sigma=sqrt(sum2/N-mean*mean);
        gsl_vector_set(results,0,mean*V);
        gsl_vector_set(results,1,sigma*V/sqrt(N));
}

void plainmcB(double f(gsl_vector* x),gsl_vector* a,gsl_vector* b,int N,gsl_vector * results){
        double V=1; for(int i=0;i<a->size;i++)V*=gsl_vector_get(b,i)-gsl_vector_get(a,i);

		double sum=0;
		double sum2=0;

		gsl_vector * x=gsl_vector_alloc(a->size);

		for(int i=0;i<N/2;i++){
                halton(i,x);   // gsl_vector_get(a,i)+halton(i,x)*(gsl_vector_get(b,i)-gsl_vector_get(a,i)))
				for(int i=0;i<x->size;i++) gsl_vector_set(x,i,gsl_vector_get(a,i)+gsl_vector_get(x,i)*(gsl_vector_get(b,i)-gsl_vector_get(a,i)));
                double fx=f(x); sum+=fx;
                }
		for(int i=N/2;i<N;i++){
            halton(i,x);   // gsl_vector_get(a,i)+halton(i,x)*(gsl_vector_get(b,i)-gsl_vector_get(a,i)))
			for(int i=0;i<x->size;i++) gsl_vector_set(x,i,gsl_vector_get(a,i)+gsl_vector_get(x,i)*(gsl_vector_get(b,i)-gsl_vector_get(a,i)));
            double fx=f(x); sum2+=fx;
            }

		double mean=(sum/(N/2)+sum2/(N/2))/2;
		double err=sum/(N/2)-sum2/(N/2);
        gsl_vector_set(results,0,mean*V);
        gsl_vector_set(results,1,sqrt(err*err));
}
