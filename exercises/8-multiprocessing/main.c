#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <pthread.h>

typedef struct {unsigned int seed; int N_in;} args;


void *estimate_pi(void *temp){
	int N = 1e4;
	args *arg = (args*)temp;
	unsigned int seed = (*arg).seed;
	int N_in = 0;
	for(int i = 0; i<N; i++){
		double x = (double)rand_r(&seed)/RAND_MAX;
		seed++;
		double y = (double)rand_r(&seed)/RAND_MAX;
		seed++;
		if(x*x+y*y < 1.0){
		N_in++;
		}
	}
	(*arg).N_in=N_in;
	return NULL;
}

int main(){
	//int CPUs = 4;
	int N = 1e4;
	pthread_t threadx;//, thready;
	void* attributes = NULL;
	unsigned int seed0 = time(NULL);
	printf("hello");
	args targs;
	targs.seed = seed0;
	args margs;
	margs.seed = seed0+2*N;
	pthread_create(&threadx, attributes, estimate_pi, (void*)&targs);
	/*
	void* retvaly = NULL
	pthread_create(&thready, attributes, estimate_pi, (void*)&seed);
	seed+=N;
	*/
	estimate_pi((void*)&margs);
	pthread_join(threadx,attributes);
	printf("Pi: %.10g\n", 4.0*targs.N_in/N);
	printf("Pi: %.10g\n", 4.0*margs.N_in/N);
	return 0;
}

