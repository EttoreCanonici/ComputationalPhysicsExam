#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


int main(int argc, char ** argv)  {
	
	double pi;
	double a,l;
	double x;
	int i,t,T;
	double nc,nd;
	double delta;
	double theta;	
	double scartoperc;
    FILE * gp;
    
    
    gp = popen("gnuplot","w");
    fprintf(gp, "set title 'pi'\n");
  	fprintf(gp,"plot '-' u 1:2 w l,pi\n");
	
	 	
	srand(time(NULL));
	
		
	//INIZIALIZZAZIONE
	nc=0;
	l=1.5;
	a=2;
	T = 1000000;
	
	//T = atoi(argv[1]);
	//a = atof(argv[2]);
	//l = atof(argv[3]);

	
	
	if( l>a){
		printf("l deve essere minore di a.\n");
		exit(1);
	}
	
	for(t=1; t<=T; t++){
		
		x = -a +2*a*drand48();
		theta = -M_PI + drand48()*2*M_PI;
		
		if( fabs(l*cos(theta)) > fabs(x) ){
			nc++;
		}
		
		if(t%1000 == 0){	
			pi = ((double)t/nc)*(2*l/a); 
			scartoperc = fabs((M_PI - pi)/M_PI);
			fprintf(gp,"%d	%f	%f\n", t, pi, scartoperc);
		}
	}
	
	 printf("\npi: %f	nc/n: %f	scarto: %f\n",pi, (double) nc/T, scartoperc);
	
	 fprintf(gp, "e\n");
	 fflush(gp);
	 getchar();
     fclose(gp);

	
	
	
	
}
	
	
