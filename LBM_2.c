 /*  lbe.c     */
/*  12-01-06  */
/*  25-01-07  */

/*
 *                
 *                  6   2   5
 *                   \  |  /
 *                    \ | /
 *                     \|/
 *                 3----0----1
 *                     /|\
 *                    / | \
 *                   /  |  \
 *                  7   4   8
 *
 */

# define LBOUND 1.e-15

////////////////  ARREGLOS 

int *one_d_int_array(int n)
{
  int *ptr;

  ptr = (int *) calloc(n,sizeof(int));
  return ptr;
}

double *one_d_double_array(int n)
{
  double *ptr;

  ptr = (double *) calloc(n,sizeof(double));
  return ptr;
}

/*	s[i][j], 0<=i<lx, 0<=j<ly	*/
int **two_d_int_array(int lx, int ly) 
{
  int **ptr, i;
  
  ptr = (int **) calloc(lx,sizeof(int*));
  for(i=0;i<lx;++i) 
    ptr[i] = (int *) calloc(ly,sizeof(int));
  return ptr;
}

/*	s[i][j], 0<=i<lx, 0<=j<ly	*/
double **two_d_double_array(int lx, int ly) 
{
  int i;
  double  **ptr;

  ptr = (double **) calloc(lx,sizeof(double*));
  for(i=0;i<lx;++i) 
    ptr[i] = (double*) calloc(ly,sizeof(double));
  return ptr;
}

void free_two_d_int_array(int n,int **s)
{
   int i;
   
   for (i=0;i<n;++i)
      free(s[i]);
}

void free_two_d_double_array(int n,double **s)
{
   int i;
   
   for (i=0;i<n;++i)
      free(s[i]);
}

double ***three_d_double_array(int lx, int ly, int lz) 
{
  int i,j;
  double ***ptr;
  
  ptr = (double ***) calloc(lx,sizeof(double**));
  for(i=0;i<lx;++i) {
    ptr[i] = (double **) calloc(ly,sizeof(double*));
    for(j=0;j<ly;++j) {
      ptr[i][j] = (double *) calloc(lz,sizeof(double));
    }
  }
  return ptr;
}

void free_three_d_double_array(int lx, int ly, double ***s)
{
  int i,j;
   
  for (i=0;i<lx;++i) {
    for (j=0;j<ly;++j) {
      free(s[i][j]);
    }
    free(s[i]);
  }
}

/////////////////////////////////////////////////////////////////////////////////

// OPEN ARRAYS


double *start_w(int n_dir)
{
  double *w;

  w = one_d_double_array(n_dir);
  w[0] = 4.0/9.0;
  w[1] = w[2] = w[3] = w[4] = 1./9.;
  w[5] = w[6] = w[7] = w[8] = 1./36.;
  return w;
}

int **start_c(int n_dir,int n_dim)
{
  int **c;

  c = two_d_int_array(n_dir,n_dim);
  c[0][0] = +0;    c[0][1] = +0;
  c[1][0] = +1;    c[1][1] = +0;
  c[2][0] = +0;    c[2][1] = +1;   
  c[3][0] = -1;    c[3][1] = +0;
  c[4][0] = +0;    c[4][1] = -1;
  c[5][0] = +1;    c[5][1] = +1;
  c[6][0] = -1;    c[6][1] = +1;
  c[7][0] = -1;    c[7][1] = -1;
  c[8][0] = +1;    c[8][1] = -1;
  return c;
}

int *start_opp_of(int n_dir)
{
  int *opp_of;

  opp_of = one_d_int_array(n_dir);
  opp_of[0] = 0;
  opp_of[1] = 3;
  opp_of[2] = 4;
  opp_of[3] = 1;
  opp_of[4] = 2;
  opp_of[5] = 7;
  opp_of[6] = 8;
  opp_of[7] = 5;
  opp_of[8] = 6;
  return opp_of;
}

///////////////////////////////////////////////////////////////////

void initial_condition_const(int lx,int ly,int n_dir,double rho_0,
		  double *w,double ***fn,double ***fv)
{
  int i,j,k;

  for (i=0;i<lx;++i) {
    for (j=0;j<ly;++j) {
      for (k=0;k<n_dir;++k) {
	fn[i][j][k] = 0.;
	fv[i][j][k] = w[k]*rho_0;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
//rho_u[0] densidad, rho_u[1] velocidad en x, rho_u[2] velocidad en y
void get_density_velocity(int n_dir,int **c,double *f_loc,double *rho_u)
{
  int k;

  for (k=0;k<3;++k)
    rho_u[k] = 0;
  for (k=0;k<n_dir;++k) {
    rho_u[0] += f_loc[k];
    rho_u[1] += f_loc[k]*c[k][0];
    rho_u[2] += f_loc[k]*c[k][1];
  }
  if (rho_u[0]>LBOUND) {
    rho_u[1] /= rho_u[0];
    rho_u[2] /= rho_u[0];
  } else {
    rho_u[0] = 0.;
    rho_u[1] = 0.;
    rho_u[2] = 0.;
  }
}

//external force in the x direction to simulate Poiseuille flow
//periodic boundary conditions on vertical walls, half way bounce back on
//horizontal walls
void evolve_pbc_hlfway(int lx,int ly,int n_dir,int *opp_of,int **c,double gx,
                       double tau,double *f_loc,double *rho_u,double *w,
                       double ***fn,double ***fv)
{
  int i,j,k;
  double feq,p;

  // bulk and pbc on vertical walls
  for (i=0;i<lx;++i){
    for (j=1;j<ly-1;++j) {
      f_loc = fv[i][j];
      get_density_velocity(n_dir,c,f_loc,rho_u);
      for (k=0;k<n_dir;++k) {
        p = c[k][0]*rho_u[1]+c[k][1]*rho_u[2];
        feq = w[k]*rho_u[0]*(1.+3.*p+4.5*p*p-1.5*(rho_u[1]*rho_u[1]+
                                                  rho_u[2]*rho_u[2]));
        fn[(lx+i+c[k][0])%lx][j+c[k][1]][k] = f_loc[k]+
          (feq-f_loc[k])/tau+gx*w[k]*c[k][0];
     }
    }
  }
  //halfway bounce back
  //top wall, left corner
  i = 0;
  j = ly-1;
  f_loc = fv[i][j];
  get_density_velocity(n_dir,c,f_loc,rho_u);
  for (k=0;k<n_dir;++k) {
    p = c[k][0]*rho_u[1]+c[k][1]*rho_u[2];
    feq = w[k]*rho_u[0]*(1.+3.*p+4.5*p*p-1.5*(rho_u[1]*rho_u[1]+
                                              rho_u[2]*rho_u[2]));
    switch (k) {
    case 0: case 1: case 4: case 8: //in
      fn[i+c[k][0]][j+c[k][1]][k] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    case 2: case 5: case 6: //half way bounce back
      fn[i][j][opp_of[k]] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    case 3: case 7: //pbc
      fn[(lx+i+c[k][0])%lx][j+c[k][1]][k] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    }
  }
  //  printf("Top left OK\n");
  //top wall, center
  j = ly-1;
  for (i=1;i<lx-1;++i) {
    //    printf("top i=%d\n",i);
    f_loc = fv[i][j];
    get_density_velocity(n_dir,c,f_loc,rho_u);
    for (k=0;k<n_dir;++k) {
      p = c[k][0]*rho_u[1]+c[k][1]*rho_u[2];
      feq = w[k]*rho_u[0]*(1.+3.*p+4.5*p*p-1.5*(rho_u[1]*rho_u[1]+
                                                rho_u[2]*rho_u[2]));
      switch (k) {
      case 0: case 1: case 3: case 4: case 7: case 8: //in
        fn[i+c[k][0]][j+c[k][1]][k] = f_loc[k]+(feq-f_loc[k])/tau+
          gx*w[k]*c[k][0];
        break;
     case 2: case 5: case 6: //out
        fn[i][j][opp_of[k]] = f_loc[k]+(feq-f_loc[k])/tau+
          gx*w[k]*c[k][0];
        break;
      }
    }
  }
  //  printf("Top center OK\n");
  //top wall, right corner
  i = lx-1;
  j = ly-1;
  f_loc = fv[i][j];
  get_density_velocity(n_dir,c,f_loc,rho_u);
  for (k=0;k<n_dir;++k) {
    p = c[k][0]*rho_u[1]+c[k][1]*rho_u[2];
    feq = w[k]*rho_u[0]*(1.+3.*p+4.5*p*p-1.5*(rho_u[1]*rho_u[1]+
                                              rho_u[2]*rho_u[2]));
    switch (k) {
    case 0: case 3: case 4: case 7: //in
      fn[i+c[k][0]][j+c[k][1]][k] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    case 2: case 5: case 6: //half way bounce back
      fn[i][j][opp_of[k]] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    case 1: case 8: // pbc
      fn[(i+c[k][0])%lx][j+c[k][1]][k] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    }
  }
  //bottom wall, left corner
  i = 0;
  j = 0;
  f_loc = fv[i][j];
  get_density_velocity(n_dir,c,f_loc,rho_u);
  for (k=0;k<n_dir;++k) {
    p = c[k][0]*rho_u[1]+c[k][1]*rho_u[2];
    feq = w[k]*rho_u[0]*(1.+3.*p+4.5*p*p-1.5*(rho_u[1]*rho_u[1]+
                                              rho_u[2]*rho_u[2]));
    switch (k) {
    case 0: case 1: case 2: case 5://in
      fn[i+c[k][0]][j+c[k][1]][k] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    case 4: case 7: case 8: //half way bounce back
      fn[i][j][opp_of[k]] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    case 3: case 6: // pbc
      fn[(lx+i+c[k][0])%lx][j+c[k][1]][k] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    }
  }
  //  printf("Bottom left OK\n");
  //  getchar();
  //bottom wall, center
  j = 0;
  for (i=1;i<lx-1;++i) {
    f_loc = fv[i][j];
    get_density_velocity(n_dir,c,f_loc,rho_u);
    for (k=0;k<n_dir;++k) {
      p = c[k][0]*rho_u[1]+c[k][1]*rho_u[2];
      feq = w[k]*rho_u[0]*(1.+3.*p+4.5*p*p-1.5*(rho_u[1]*rho_u[1]+
                                                rho_u[2]*rho_u[2]));
      switch (k) {
      case 0: case 1: case 2: case 3: case 5: case 6: //in
        fn[i+c[k][0]][j+c[k][1]][k] = f_loc[k]+(feq-f_loc[k])/tau+
          gx*w[k]*c[k][0];
        break;
      case 4: case 7: case 8: //out
        fn[i][j][opp_of[k]] = f_loc[k]+(feq-f_loc[k])/tau+
          gx*w[k]*c[k][0];
        break;
      }
    }
  }
  //  printf("Bottom center OK\n");
  //  getchar();
  // bottom wall, right corner
  i = lx-1;
  j = 0;
  f_loc = fv[i][j];
  get_density_velocity(n_dir,c,f_loc,rho_u);
  for (k=0;k<n_dir;++k) {
    p = c[k][0]*rho_u[1]+c[k][1]*rho_u[2];
    feq = w[k]*rho_u[0]*(1.+3.*p+4.5*p*p-1.5*(rho_u[1]*rho_u[1]+
                                              rho_u[2]*rho_u[2]));
    switch (k) {
    case 0: case 2: case 3: case 6: //in
      fn[i+c[k][0]][j+c[k][1]][k] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    case 4: case 7: case 8: //half way bounce back
      fn[i][j][opp_of[k]] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    case 1: case 5: // pbc
      fn[(i+c[k][0])%lx][j+c[k][1]][k] = f_loc[k]+(feq-f_loc[k])/tau+
        gx*w[k]*c[k][0];
      break;
    }
  }
}
  
/////////////////  get

//velocidad en la direcci'on horizontal en lx/2
void get_velocity_middle(int lx2,int ly,int n_dir,int **c,double *f_loc,
			 double *rho_u,double *u,double ***fv)
{
  int j;

  for (j=0;j<ly;++j) {
    f_loc = fv[lx2][j];
    get_density_velocity(n_dir,c,f_loc,rho_u);
    u[j] = rho_u[1];
  }
}

