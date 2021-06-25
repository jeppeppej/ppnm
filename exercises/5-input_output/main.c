#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char** argv){//the note "command line" says to use this
	if(argc<2) fprintf(stderr,"%s: there were no arguments\n",argv[0]);
	else{
		printf("\tx\t|\tsin(x)\t|\tcos(x)\t\n");
		printf("\t\t|\t\t|\t\t\n");
		for(int i=0;i<argc;i++){
			double x = atof(argv[i]);
			printf("\t%.4g\t|\t%.4g\t|\t%.4g\t\n",x,cos(x),sin(x));
		}
}
}
