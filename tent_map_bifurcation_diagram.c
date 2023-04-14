#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(int argc, char ** argv) {
  double r, Deltar;
  double rmin, rmax;
  int i,j,k,n;
  double x, x0, x1;
  FILE *gp;
  
  gp = popen("gnuplot","w");
  
   
	//INIZIALIZZAZIONE:
	n=10000;
	rmin = 1;
	rmax = 2;
	Deltar = 0.01;
	//x0 = drand48();
	x0 = 0.5;
	printf("%lf\n",x0);


  	fprintf(gp, "set xrange [%lf:%lf]\n",rmin,rmax);
    fprintf(gp,"set yrange [0:1]\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "set title 'Bifurcation diagram at x = %lf'\n", x0);
    fprintf(gp, "set xlabel 'r'\n");
    fprintf(gp, "set ylabel 'x'\n");
    fprintf(gp, "plot '-' u 2:1 w p\n"); 

		
	
	for(r=1; r<=2; r += Deltar){

		x = x0;
		fprintf(gp, "%lf %lf\n", x, r);
		
		for(i=0; i<n; i++){
		  
			if(x <0.5 && x >= 0){
				x1 =  r*x;
			}else{
				x1 = r*(1-x);
			}
		
			
			fprintf(gp, "%lf %lf\n", x1, r);
			x = x1; 
		}
	}

  fprintf(gp, "end\n");
  fflush(gp);
  getchar();
  fclose(gp);


}
