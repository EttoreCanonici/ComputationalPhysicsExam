#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


int main(int argc, char ** argv) {
  double r;
  int i,n;
  double x, x0, x1;
  FILE *gp;


	srand48(time(NULL));

	//INIZIALIZZAZIONE:
	n=1000000;
	r = 3;
	//x0 = drand48();
	x0 = 0.7;
	//x0 = M_PI/4;
	
	
	
	gp = popen("gnuplot", "w");
	fprintf(gp, "set xrange [0:1]\n");
   // fprintf(gp,"set yrange [-%lf:%lf]\n",r,r);
    fprintf(gp,"set yrange [%lf:%lf]\n",0.,1.);
    fprintf(gp, "set size square\n");
   
	//fprintf(gp, "set title 'r = %lf, x0 = %lf'\n", r, x0);
	fprintf(gp, "set title 'r = %lf, [(r-r^2)/2, r/2] = [%lf,%lf], x0 = %lf'\n", r, (r-r*r)/2, r/2, x0);
	//fprintf(gp, "set title 'r = %lf, r/(r+1) = %lf, x0 = %lf'\n", r, r/(r+1), x0);
	 fprintf(gp, "set xlabel 'x_n'\n");
    fprintf(gp, "set ylabel 'f(x_n)'\n");
	fprintf(gp, "plot '-' u 1:2 w p\n"); 
  
  
	x = x0;
	x1 = 0;
	
	fprintf(gp, "%lf %lf\n", 0., x);
	
	for(i=0; i<n; i++){
	  
		if(x <0.5 && x >= 0){
			x1 =  r*x;
		}else{
			x1 = r*(1-x);
		}
	
		
		fprintf(gp, "%lf %lf\n", x, x1);
		//printf( "%lf %lf\n", x, x1);
		
		x = x1; 
	}
	

  fprintf(gp, "end\n");
  fflush(gp);
  getchar();
  fclose(gp);


}
