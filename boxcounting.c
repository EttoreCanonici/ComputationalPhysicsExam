#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


double tent_map(double r, double x);

int main(int argc, char ** argv) {
  double x,x1,r,N,eps;
  int i,l,n;
  int * count;
  int K;
  int conteggi = 0;
  double boxcounting = 0;
  //FILE *gp;
  
  
  //gp = popen("gnuplot","w");
  srand48(time(NULL));

  //INIZIALIZZAZIONE
	n=10000;
	r = 2.5;
	eps = 0.01;
	//x = drand48();
	//x = M_PI/4;
	x = 0.7;
	conteggi = 0;
  
  
  	//gp = popen("gnuplot", "w");
	//fprintf(gp, "set xrange [0:1]\n");
    //fprintf(gp,"set yrange [%lf:%lf]\n",0.,1.);
    //fprintf(gp, "set size square\n");   
	//fprintf(gp, "set title 'r = %lf'\n", r);
	//fprintf(gp, "plot '-' u 1:2 w p\n"); 
  
  
  
  
	N = 1/eps;
	K = (int)N;
	
	count = calloc(K, sizeof(int));
  
	for(i=0; i<K; i++){
		count[i] = 0;
	}
  
  
	for(i=0; i<n; i++){
		
		x1 = tent_map(r,x);		
		//fprintf(gp, "%lf %lf\n", x, x1);
		x = x1; 
		
		for(l=0; l<K; l++){
			if( (x >= l*K) && (x < (l+1)*K) ){
				count[l] = count[l] + 1;
			}
		}
		
	}
	
	
	
	for(i=0; i<K; i++){
		if( count[i] > 0){
			conteggi += 1;
		}
	}
	
	
	boxcounting = - log(conteggi
	)/log(eps);
	printf("DIMENSIONE FRATTALE: %lf, epsilon: %lf\n", boxcounting, eps);
	//fprintf(gp,"%lf	%lf\n", eps, boxcounting);

  //fprintf(gp, "end\n");
  //fflush(gp);
  //getchar();
  //fclose(gp);
 
  
  
    
 }


double tent_map(double r, double x){ 
	
	double x1;

	if(x <0.5 && x >= 0){
			x1 =  r*x;
	}else{
			x1 = r*(1-x);
	}

return x1;
}
