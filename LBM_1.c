 /*  Poiseuille flow
 *
 *                             lx
 *      --------------------------------------------------
 *      |                                                |
 *      |                                                |
 *      |                                                |
 *      |                                                |
 *      |                                                |    ly
 *      |                                                |
 *      |                                                |
 *      |                                                |
 *      |                                                |
 *      |                                                |
 *      --------------------------------------------------
 *
 * To compile and run
 *  $ gcc lbe-00.c -lm -o lbe-00.e
 *  $ ./lbe-00.e 500 51 10000  1.0 1.0 2000
 *
 *  Initial condition:
 *  rho_0 throughout the cavity
 *  u = 0 throughout the cavity
 *
 *  Boundary conditions:
 *  half way bounce back on top and bottom horizontal walls. 
 *  Periodic boundary conditions on the vertical walls,
 *  The corners
 *  have directions going in the cavity, directions with halfway bounce 
 *  back, and directions with periodic boundary conditions.  
 *  
 *  Evaluate:
 *  u(lx/2,ly/2,t) as a function of time. This is in 
 *  lbe-00_ut_rho_x500_y51_t1p00_n0p17_Re2p00e+03.dat
 *
 *  horizontal velocity in a vertical line in the middle of the cavity
 *  at t_max. This is in 
 *  lbe-00_uy_rho_x500_y51_t1p00_n0p1667_Re2p00e+03.dat
 *
 *  These files have names that make sense and change if the parameters are
 *  changed.
 *  
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include "lbe-ale.c"

int main(int argc,char *argv[]) {
  char s_file[80];
  char s_file1[80];
  int dlta,j,lx,ly,n_dim,n_dir,n_ptos;
  int *opp_of,**c;
  int t,t_max,t_res;
  double cs,gx,nu,Re,rho_0,tau,umax;
  double *f_loc,*rho,*rho_u,*u,*w;
  double ***fv,***fn,***tmp;
  FILE *fout;
  FILE *gp, *gp1;

  // input
  if (argc <7) {
    fprintf(stderr, "usage: %s LX LY TMAX TAU RHO_0 Re\n", argv[0]);\
    exit(1);
  }
  lx = atoi(argv[1]); // 500
  ly = atoi(argv[2]); // 51
  t_max = atoi(argv[3]);  // 20000
  tau = atof(argv[4]);  // 0.5 < tau < 10
  rho_0 = atof(argv[5]); // 1
  Re = atof(argv[6]); // 0< Re < 1e6
  //initial values
  cs = 1./sqrt(3.);
  dlta = lx/100;
  n_dim = 2;
  n_dir = 9;
  n_ptos = 100;
  nu = (tau-0.5)/3.;
  gx = 8.*nu*nu*Re/ly/ly/ly;
  umax = gx*ly*ly*ly/nu/nu/8;
  printf("nu=%f  gx=%f\n",nu,gx);  
  if (t_max<50)
    t_res = 1;
  else t_res = t_max/n_ptos;
  //start arrays
  c = start_c(n_dir,n_dim);
  opp_of = start_opp_of(n_dir);
  f_loc = one_d_double_array(n_dir);
  rho_u = one_d_double_array(3);
  rho = one_d_double_array(lx);
  u = one_d_double_array(ly);
  w = start_w(n_dir);
  fv = three_d_double_array(lx,ly,n_dir);
  fn = three_d_double_array(lx,ly,n_dir);
  //files
  sprintf(s_file,"lbe-00_ut_rho_x%d_y%d_t%1.2f_n%1.2f_Re%1.2e.dat",lx,ly,tau,nu,Re);
  printf("evoluz. temporale: %s\n", s_file);
  fout = fopen(s_file,"w");
  // initial condition
  initial_condition_const(lx,ly,n_dir,rho_0,w,fn,fv);
  // time evolution
  for (t=0;t<=t_max;++t) { 
    evolve_pbc_hlfway(lx,ly,n_dir,opp_of,c,3.*gx,tau,f_loc,rho_u,w,fn,fv);
    if ((t%t_res)==0) {
      get_velocity_middle(lx/2,ly,n_dir,c,f_loc,rho_u,u,fn);
      fprintf(fout,"%f\t%f\n",(double) t*nu/ly/ly,u[ly/2]*ly/nu);
    }
    tmp = fv;
    fv = fn;
    fn = tmp;
  }
  fclose(fout);
  // u as function of y at t=t_max
  sprintf(s_file1,"lbe-00_uy_rho_x%d_y%d_t%1.2f_n%1.4f_Re%1.2e.dat",lx,ly,tau,nu,Re);
  printf("profile vel: %s\n",s_file1);
  fout = fopen(s_file1,"w");
  for (j=0;j<ly;++j) {
    fprintf(fout,"%f\t%f\n",(double) j/ly,u[j]*ly/nu);
  }
  fclose(fout);
  gp = popen("gnuplot", "w");
  fprintf(gp, "plot '%s' w l\n", s_file);
  fflush(gp);
  gp1 = popen("gnuplot", "w");
  fprintf(gp1, "plot '%s' w l\n", s_file1);
  printf("premi un tasto per uscire\n");
  fflush(gp1);
  getchar();
  free(f_loc);
  free(rho_u);
  free(u);
  free(w);
  free_two_d_int_array(n_dir,c);
  free_three_d_double_array(lx,ly,fv);
  free_three_d_double_array(lx,ly,fn);
  return 0;
}


