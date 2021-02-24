#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

//struct vector {int n, double a[]}
//struct vector my_vector;
typedef struct {int n, double v[]} vector;
vector* vector_allocate(int n){
	vector* v=malloc(sizeof(vector));
	(*v).n=n;
	('v).data=malloc(n*sizeof(double));
	return v;
	}
void vector_free(vector* v){
	free((*v).data);
	free(v);
	}
void vector_set(vector* v,int i, double value){
	if(i<0)printf("i<0\n");
	if(i>n-1)printf("i>n-1\n");
	(*v).
	}
double vector_get(vector* v,int i){
	
	}

void set0(double x){x=0; }

void set0p(double *x){ (*x)=0; }

void print_array(int n, double a[]){
	for(int i=0;i<n;i++)printf("print_array: a[%d]=%g\n",i,a[i]);
	}

int main(){
	double y=1;
	set0(y);
	printf("y after set0  = %g\n",y);
	set0p(&y);
	printf("y after set0p = %g\n",y);
	int n=5;
	double v[5]; //en array med 5 indgange af typen double
	for(int i=0;i<n;i++){
		v[i]=i;
		}
	int i=0; while(i<n) {
		printf("v[%d]=%g\n",i,v[i]);
		i++;
	}
	v[100]=1
	int N=7;
	double *a = malloc(N*sizeof(double));
	for(int i=0;i<N;i++) a[i]=i+100;

free(a);
return 0;
}
