/* program to calculate weighting potentials of coaxial Ge detectors
   by relaxation

   loosely based on David's ppc8.c

   to run: ./fields3d geometry_setup.dat field_calc_setup.dat ev.dat [wp.dat]
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
float fminf(float x, float y);
#include <unistd.h>
#include <ctype.h>

#include "field_init.h"
#include "signal_calc_util.h"
#include "point.h"
#ifdef GRETINA
#include "assym_detector.h"
#endif

#define e_over_E (1.413*8/3)
#define NUMITER 10000
#define ITER_REPEAT_N 1
#define MAX_FNAME_LEN 512
#define DEFAULT_FIELD_SETUP "field_setup.dat"


static int grid_init(char *field_fname, char *geometry_fname);
static int ev_calc(char *ev_fname);
static int wp_calc(char *wp_fname);
static int malloc_arrays(void);
static int fname_insert(char fname[MAX_FNAME_LEN], char dest[MAX_FNAME_LEN], 
		 char *str);
static int write_ev(char *fname);
static int write_wp(char *fname);
int do_relaxation(int numiter, int efld_calc);
static int free_arrays(void);

static int numiter = NUMITER;
static int iter_repeat_n = ITER_REPEAT_N;

float xmin, xmax, xstep, ymin, ymax, ystep, zmin, zmax, zstep;
int numx, numy, numz;
static float Vb;
static float rho0, drho, lambda, rho_b, rho_c;
static float ***v[2];
char ***point_type;
float ***fraction_in_det;
static int ncontacts;

int main(int argc, char *argv[]){
  char opt;
#ifdef GRETINA
  char optstring[] = "n:r:g:f:v:w:ht:q:";
  char ctype;
  int qnum;
  char *endp;
#else
  char optstring[] = "n:r:g:f:v:w:h";
#endif
  extern char *optarg;
  extern int optind, opterr, optopt;
  char geometry_fname[MAX_FNAME_LEN], field_fname[MAX_FNAME_LEN],
    ev_fname[MAX_FNAME_LEN], wp_fname[MAX_FNAME_LEN];

  *geometry_fname = *field_fname = *ev_fname = *wp_fname = '\0';
  while ((opt = getopt(argc, argv, optstring)) >= 0){
    if (opt == 'h' || opt == '?'){
#ifdef GRETINA
      printf("usage: %s [-n numiter] [-r iter_repeat_n] [-g geometry.dat]\n"
	     "        [-f field_calc.dat] [-v ve.out] [-w wp_out] -t[A|B] -q[1|n]\n", argv[0]);
#else
      printf("usage: %s [-n numiter] [-r iter_repeat_n] [-g geometry.dat]\n"
	     "        [-f field_calc.dat] [-v ve.out] [-w wp_out]\n", argv[0]);
#endif
      return 0;
    }
    if (opt == 'n'){
      numiter = atoi(optarg);
      if (numiter < 0){
	printf("invalid number of iterations: %s\n", optarg);
	return 1;
      }	
    }else if (opt == 'r'){
      iter_repeat_n = atoi(optarg);
      if (iter_repeat_n < 0){
	printf("invalid number of repeats: %s\n", optarg);
	return 1;
      }
    }else if (opt == 'g'){
      strcpy(geometry_fname, optarg);
    }else if (opt == 'f'){
      strcpy(field_fname, optarg);
    }else if (opt == 'v'){
      strcpy(ev_fname, optarg);
    }else if (opt == 'w'){
      strcpy(wp_fname, optarg);
#ifdef GRETINA
    }else if (opt == 't'){
      ctype = toupper(*optarg) - 'A';
      set_crystal_geometry(ctype);
    }else if (opt == 'q'){
      qnum = strtol(optarg, &endp, 10);
      if (endp ==optarg){
	printf("failed to parse quad number: %s (should be integer > 0)\n",optarg);
	return 1;
      }       
      set_quad_number(qnum);
#endif
    }else{
      printf("unknown option: %c\n", opt);
    }
  }
  if (!*ev_fname && !*wp_fname){
    printf("error! no output files supplied\n");
    return 1;
  }
  printf("number of iterations: %d X %d\n", numiter, iter_repeat_n);
  if (*geometry_fname)
    printf("geometry file name: %s\n", geometry_fname);
  else
    printf("will use default geometry conf. file\n");
  if (*field_fname)
    printf("field calc. file name: %s\n", field_fname);
  else
    printf("will use default field calc. conf. file\n");
  if (*ev_fname)
    printf("field output file name: %s\n", ev_fname);
  else
    printf("will not calculate E,V\n");
  if (*wp_fname)
    printf("wp output file name: %s\n", wp_fname);
  else
    printf("will not calculate wp\n");

  if (grid_init(field_fname, geometry_fname) != 0){
    error("failed to init field calculations\n");
    return 1;
  }
  if (*ev_fname) ev_calc(ev_fname);
  if (*wp_fname) wp_calc(wp_fname);
     
  return 0;
}


static int grid_init(char *field_fname, char *geometry_fname){
  FILE *fp;
  char line[MAX_LINE];

  if (*field_fname){
    if ((fp = fopen(field_fname,"r")) == NULL){
      error("Failed to open field calc. setup file: %s\n", field_fname);
      return -1;
    }
  }else{
    if ((fp = fopen(DEFAULT_FIELD_SETUP,"r")) == NULL){
      error("Failed to open field calc. setup file: %s\n", DEFAULT_FIELD_SETUP);
      return -1;
    }
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f %f %f", &xmin, &xmax, &xstep) != 3
      || read_setup_line(fp, line) != 0
      || sscanf(line, "%f %f %f", &ymin, &ymax, &ystep) != 3
      || read_setup_line(fp, line) != 0
      || sscanf(line, "%f %f %f", &zmin, &zmax, &zstep) != 3){
    error("Failed to read geometry boundary data from file: %s\n", field_fname);
    fclose(fp);
    return -1;
  }
  
  numx = (int)rint((xmax - xmin)/xstep) + 1;
  numy = (int)rint((ymax - ymin)/ystep) + 1;
  numz = (int)rint((zmax - zmin)/zstep) + 1;
  
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &Vb) != 1){
    error("failed to read bias voltage from file: %s\n", field_fname);
    fclose(fp);
    return -1;
  }

  lambda = rho_b = rho_c = 0.0;
  if (read_setup_line(fp, line) != 0
      || ((sscanf(line, "%f %f %f %f",  
		  &rho0, &drho, &lambda, &rho_b) != 4)
	  &&  sscanf(line, "%f %f", &rho0, &drho) != 2)){
    error("failed to read impurity & gradient from file: %s\n", field_fname);
    fclose(fp);
    return -1;
  }
  if (malloc_arrays() < 0) return -1;

  if ((ncontacts = geometry_init(geometry_fname)) <0){
    error("failed to init geometry\n");
    return 1;
  }

  if (lambda != 0 || rho_b != 0){
    rho_c = -drho/10 -rho_b/zlen_detector()*(1-exp(-zlen_detector()/lambda));
    printf("rho_c = %f\n", rho_c);
  }

  fclose(fp);
  
  printf("\n\n"
	 "xmin, xmax, xstep: %f %f %f\n"
	 "ymin, ymax, ystep: %f %f %f\n"
	 "zmin, zmax, zstep: %f %f %f\n"
	 "Vb, rho0, drho: %f %f %f\n"
	 "lambda, rho_b, rho_c: %f %f %f\n", 
	 xmin, xmax, xstep, ymin, ymax, ystep, zmin, zmax, zstep,
	 Vb, rho0, drho, lambda, rho_b, rho_c);
  printf("\n");



   return 0;
}


int do_relaxation(int numiter, int efld_calc){
  int i,j,k;
  int old,new;
  int iter;
  float sum_dif, max_dif, dif;
  float mean, denom;
  float z, x, y, r;
  char undep;
  char vchanged;

  old = 1; new = 0;
  for (i = 0; i < numx; i++){
    for (j = 0; j < numy; j++){
      for (k = 0; k < numz; k++){
	if (point_type[i][j][k] != INSIDE
#ifdef CALC_VFIELD
	    && point_type[i][j][k] != OUTSIDE
	    
#endif
	    ) continue;//contacts are @ fixed potential
	if (i - 1 < 0 || i+1 >= numx || j -1 < 0 || j+1 >= numy
	    || k - 1 < 0 || k +1 >= numz){
	  printf("error for point %d %d %d\n", i,j,k);
	  exit(1);
	}
#ifndef CALC_VFIELD
	if ((point_type[i-1] == OUTSIDE
	     && point_type[i+1] == OUTSIDE)
	    ||(point_type[j-1] == OUTSIDE
	       && point_type[j+1] == OUTSIDE)
	    ||(point_type[k-1] == OUTSIDE
	       && point_type[k+1] == OUTSIDE)){
	  printf("error2 for point %d %d %d\n", i,j,k);
	}
#endif
      }
    }
  }


  printf("starting relaxation...\n");
  for (iter = 0; iter < numiter; iter++){
    old = old == 0;
    new = new == 0;
    sum_dif = 0.0;
    max_dif = 0.0;
    vchanged = 0;

    for (i = 0; i < numx; i++){
      x = xmin + i*xstep;
      for (j = 0; j < numy; j++){
	y = ymin + j*ystep;
	r = sqrt(x*x + y*y);
	for (k = 0; k < numz; k++){
	  z = zmin + k*zstep;
	  if (point_type[i][j][k] != INSIDE
#ifdef CALC_VFIELD
	      && point_type[i][j][k] != OUTSIDE

#endif
	      ) continue;//contacts are @ fixed potential
	  mean = 0;
	  denom = 0;
#ifdef CALC_VFIELD
	  /*assume that boundary between detector and vaccuum 
	    is always on voxel centers, perp. to z axis,
	  and that the boundary is on an "inside" voxel*/
	  if (point_type[i][j][k] == INSIDE
	      && point_type[i][j][k+1] == OUTSIDE
	      && point_type[i][j][k-1] == INSIDE){
	    mean = v[old][i-1][j][k]*16
	      + v[old][i+1][j][k]*16
	      + v[old][i][j-1][k]*16
	      + v[old][i][j+1][k]*16
	      + v[old][i][j][k-1]*16
	      + v[old][i][j][k+1]*1;
	    denom = 16*5+1;
	  }else	if (point_type[i][j][k] == INSIDE
		    && point_type[i][j][k+1] == INSIDE
		    && point_type[i][j][k-1] == OUTSIDE){
	    mean = v[old][i-1][j][k]*16
	      + v[old][i+1][j][k]*16
	      + v[old][i][j-1][k]*16
	      + v[old][i][j+1][k]*16
	      + v[old][i][j][k-1]*1
	      + v[old][i][j][k+1]*16;
	    denom = 16*5+1;
	  }else{
	    mean = v[old][i-1][j][k]
	      + v[old][i+1][j][k]
	      + v[old][i][j-1][k]
	      + v[old][i][j+1][k]
	      + v[old][i][j][k-1]
	      + v[old][i][j][k+1];
	    denom = 6;
	  }
#else
	  denom = 6;
	  if (point_type[i-1][j][k] == OUTSIDE){
	    mean += v[old][i+1][j][k];
	  }else{
	    mean += v[old][i-1][j][k];
	  }
	  if (point_type[i+1][j][k] == OUTSIDE){
	    mean += v[old][i-1][j][k];
	  }else{
	    mean += v[old][i+1][j][k];
	  }
	  if (point_type[i][j-1][k] == OUTSIDE){
	    mean += v[old][i][j+1][k];
	  }else{
	    mean += v[old][i][j-1][k];
	  }
	  if (point_type[i][j+1][k] == OUTSIDE){
	    mean += v[old][i][j-1][k];
	  }else{
	    mean += v[old][i][j+1][k];
	  }
	  if (point_type[i][j][k-1] == OUTSIDE){
	    mean += v[old][i][j][k+1];
	  }else{
	    mean += v[old][i][j][k-1];
	  }
	  if (point_type[i][j][k+1] == OUTSIDE){
	    mean += v[old][i][j][k-1];
	  }else{
	    mean += v[old][i][j][k+1];
	  }
#endif	
#ifdef CALC_VFIELD
	  if (fraction_in_det[i][j][k] == 0.0){
#else	    
	  if (point_type[i][j][k] == OUTSIDE){
#endif
	    v[new][i][j][k]  = mean/denom;
	  }else{
	    if (efld_calc){
	      if (lambda == 0.0 || rho_b == 0.0)
		v[new][i][j][k] = mean/denom 
		  + (rho0+drho*z/10)
#ifdef CALC_VFIELD
		  *fraction_in_det[i][j][k]
#endif
		  *(xstep*xstep + ystep*ystep + zstep*zstep)
		  *e_over_E/denom;
	      else
		v[new][i][j][k] = mean/denom 
		  + (rho0-rho_b*(1-exp(-z/lambda))-rho_c*z) 
#ifdef CALC_VFIELD
		  *fraction_in_det[i][j][k]
#endif
		  *(xstep*xstep + ystep*ystep + zstep*zstep)
		  *e_over_E/denom;
	    }else{
	      v[new][i][j][k] = mean/denom; 
	    }
	  }
	  if (v[new][i][j][k] != v[old][i][j][k]) vchanged=1;
	  dif = v[old][i][j][k] - v[new][i][j][k];
	  if (dif < 0.0) dif = -dif;
	  sum_dif += dif;
	  if (max_dif < dif) max_dif = dif;
	}
      }
    }
    if ((iter < 1000 && iter %100 == 0) || iter%1000 == 0){
      printf("%5d %.3e %.3e\n", 
	     iter, max_dif, sum_dif/(numx*numy*numz));
      fflush(stdout);
    }
    if (max_dif == 0.0f) break;
    if (!vchanged) break;
  }
  if (!vchanged && max_dif > 0) 
    printf("No change in last iteration!\n");
  printf("%5d %.3e %.3e\n", 
	 iter, max_dif, sum_dif/(numx*numy*numz));

  undep = 0;
  for (i = 0; i < numx && !undep; i++){
    for (j = 0; j < numy && !undep; j++){
      for (k = 0; k < numz && !undep; k++){
	if (v[new][i][j][k] > Vb){
	  printf("detector is undepleted!\n");
	  undep = 1;
	}
      }
    }
  }
  printf("...  relaxation done\n");

  return 0;
}


static int ev_calc(char *ev_fname){
  int i;
  char fname[MAX_FNAME_LEN], str[MAX_FNAME_LEN];
  
  printf("\n\n ---- starting EV calculation --- \n");
  fflush(stdout);
  init_ev_calc(v, Vb);
  snprintf(str, MAX_FNAME_LEN, "-%d", 0);
  fname_insert(ev_fname, fname, str);
  printf("saving initial guess for e, v in %s\n", fname);
  write_ev(fname);

  for (i = 0; i < iter_repeat_n; i++){
    printf("starting ev calculations for iteration number: %d\n", i*numiter+1);
    do_relaxation(numiter, 1);
    if (iter_repeat_n > 1){
      snprintf(str, MAX_FNAME_LEN, "-%d", (i+1)*numiter);
      fname_insert(ev_fname, fname, str);
      printf("saving current e, v  in %s\n", fname);
      write_ev(fname);
    }
  }
  if (iter_repeat_n > 0){
    printf("saving final e, v  in %s\n", ev_fname);
    write_ev(ev_fname);
  }

  return 0;
}

static int wp_calc(char *wp_fname){
  int i, nc;
  char fname[MAX_FNAME_LEN], str[MAX_FNAME_LEN];

  printf("\n\n ---- starting WP calculation --- \n");
  for (nc = 0; nc < ncontacts; nc++){
    //for (nc = 33; nc < ncontacts; nc++){
    init_wp_calc(v,nc);
    snprintf(str, MAX_FNAME_LEN, "-seg%d_%d", nc, 0);
    fname_insert(wp_fname, fname, str);
    printf("saving initial guess for wp in %s\n", fname);
    write_wp(fname);
    
    for (i = 0; i < iter_repeat_n; i++){
      printf("starting wp calculations for iteration number: %d\n", i*numiter+1);
      do_relaxation(numiter, 0);
      if (iter_repeat_n > 1){
	snprintf(str, MAX_FNAME_LEN, "-seg%d_%d", nc, (i+1)*numiter);
	fname_insert(wp_fname, fname, str);
	printf("saving current wp  in %s\n", fname);
	write_wp(fname);
      }
    }
    if (iter_repeat_n > 0){
      snprintf(str, MAX_FNAME_LEN, "-seg%d", nc);
      fname_insert(wp_fname, fname, str);
      printf("saving final wp in %s\n", fname);
      write_wp(fname);
    }
  }

  geometry_finalize();
  free_arrays();

  return 0;
}

int write_ev(char *fname){
  int i,j,k;
  FILE *fp = NULL;
  struct efld{float x; float y; float z;} ***e;


  if (fname != NULL){
    if ((fp = fopen(fname, "w")) == NULL){
      error("Failed to open output file: %s\n", fname);
    }
  }else{
    fp = stdout;
  }

  if ((e = malloc(numx*sizeof(*e))) == NULL){
    error("malloc failed\n");
    fclose(fp);
    return -1;
  }
  for (i = 0; i < numx; i++){
    if ((e[i] = malloc(numy*sizeof(*e[i]))) == NULL){
      error("malloc failed\n");
      fclose(fp);
      return -1;
    }
    for (j = 0; j < numy; j++){
      if ((e[i][j] = calloc(numz, sizeof(*e[i][j]))) == NULL){
	error("malloc failed\n");
	fclose(fp);
	return -1;
      }
    }
  }
  for (i = 0; i < numx; i++){
    for (j = 0; j < numy; j++){
      for (k = 0; k < numz; k++){
	if (point_type[i][j][k] == OUTSIDE) continue;
	//if (point_type[i][j][k] == CONTACT_0) continue;
	//if (point_type[i][j][k] == CONTACT_VB) continue;
	if (i == 0 || (point_type[i-1][j][k] == OUTSIDE))
	  e[i][j][k].x = -(v[0][i+1][j][k] - v[0][i][j][k])/xstep;
	else if (i == numx -1 || (point_type[i+1][j][k] == OUTSIDE))
	  e[i][j][k].x = -(v[0][i][j][k] - v[0][i-1][j][k])/xstep;
	else
	  e[i][j][k].x = -(v[0][i+1][j][k] - v[0][i-1][j][k])/xstep/2;
	if (j == 0 || (point_type[i][j-1][k] == OUTSIDE))
	  e[i][j][k].y = -(v[0][i][j+1][k] - v[0][i][j][k])/ystep;
	else if (j == numy -1 || point_type[i][j+1][k] == OUTSIDE)
	  e[i][j][k].y = -(v[0][i][j][k] - v[0][i][j-1][k])/ystep;
	else
	e[i][j][k].y = -(v[0][i][j+1][k] - v[0][i][j-1][k])/ystep/2;
	if (k == 0 || (point_type[i][j][k-1] == OUTSIDE))
	  e[i][j][k].z = -(v[0][i][j][k+1] - v[0][i][j][k])/zstep;
	else if (k == numz -1 || (point_type[i][j][k+1] == OUTSIDE))
	  e[i][j][k].z = -(v[0][i][j][k] - v[0][i][j][k-1])/zstep;
	else
	  e[i][j][k].z = -(v[0][i][j][k+1] - v[0][i][j][k-1])/zstep/2;	
      }
    }
  }
  //fprintf(fp, "#x y z V |E| Ex Ey Ez\n");
  fprintf(fp, "#x y z Ex Ey Ez |E| V\n");
#define SQ(x) ((x)*(x))
  for (i = 0; i < numx; i++){
    for (j = 0; j < numy; j++){
      for (k = 0; k < numz; k++){
	fprintf(fp, "%.5f %.5f %.5f\t %10e %10e %10e %10e %10e\n",
		(xmin+i*xstep)/1000, 
		(ymin+j*ystep)/1000, (zmin+k*zstep)/1000, 		
		e[i][j][k].x*1000, e[i][j][k].y*1000, e[i][j][k].z*1000,
		sqrt(SQ(e[i][j][k].x)+ SQ(e[i][j][k].y) 
		     + SQ(e[i][j][k].z))*1000,
		v[0][i][j][k]*1000);
      }
    }
  }
  for (i = 0; i < numx; i++){
    for (j = 0; j < numy; j++){ 
      free(e[i][j]);
    }
    free(e[i]);
  }
  free(e);
  
  if (fp != stdout) fclose(fp);
  return 0;
}

int write_wp(char *fname){
  int i,j,k;
  FILE *fp = NULL;

  if (fname != NULL){
    if ((fp = fopen(fname, "w")) == NULL){
      error("Failed to open output file: %s\n", fname);
    }
  }else{
    fp = stdout;
  }

  fprintf(fp, "#x y z WP\n");
  for (i = 0; i < numx; i++){
    for (j = 0; j < numy; j++){
      for (k = 0; k < numz; k++){
	fprintf(fp, "%.5f %.5f %.5f\t %10e\n",
		(xmin+i*xstep)/1000, (ymin+j*ystep)/1000, (zmin+k*zstep)/1000, 
		v[0][i][j][k]);
      }
    }
  }
  
  if (fp != stdout) fclose(fp);
  return 0;

}


static int malloc_arrays(void){
  int i,j;

  if ((v[0] = malloc(numx*sizeof(*v[0]))) == NULL
      ||(v[1] = malloc(numx*sizeof(*v[1]))) == NULL
      || (point_type = malloc(numx*sizeof(*point_type))) == NULL
      || (fraction_in_det = malloc(numx*sizeof(*fraction_in_det))) == NULL){
    error("malloc failed\n");
    return -1;
  }
  for (i = 0; i < numx; i++){
    if ((v[0][i] = malloc(numy*sizeof(*v[0][i]))) == NULL
	|| (v[1][i] = malloc(numy*sizeof(*v[1][i]))) == NULL
	|| (point_type[i] = malloc(numy*sizeof(point_type[i]))) == NULL
	|| (fraction_in_det[i] = malloc(numy*sizeof(fraction_in_det[i]))) == NULL){
      error("malloc failed\n");
      return -1;
    }
    for (j = 0; j < numy; j++){
      if ((v[0][i][j] = malloc(numz*sizeof(*v[0][i][j]))) == NULL
	  || (v[1][i][j] = malloc(numz*sizeof(*v[1][i][j]))) == NULL
	  || (point_type[i][j] = calloc(numz,sizeof(point_type[i][j]))) 
	  == NULL
	  || (fraction_in_det[i][j] = calloc(numz,sizeof(fraction_in_det[i][j]))) 
	  == NULL){
	error("malloc failed\n");
	return -1;
      }     
    }
  }

  return 0;
}

static int free_arrays(void){
  int i, j;
  for (i = 0; i < numx; i++){
    for (j = 0; j < numy; j++){
      free(v[0][i][j]);
      free(v[1][i][j]);
      free(point_type[i][j]);
      free(fraction_in_det[i][j]);
    }
    free(v[0][i]);
    free(v[1][i]);
    free(point_type[i]);
    free(fraction_in_det[i]);
  }
  free(v[0]);
  free(v[1]);
  free(point_type);
  free(fraction_in_det);

  return 0;
}

static int fname_insert(char fname[MAX_FNAME_LEN], char dest[MAX_FNAME_LEN], 
		 char *str){
  char *cp1, *cp2;

  strncpy(dest, fname, MAX_FNAME_LEN);
  cp1 = rindex(dest,'.');
  cp2 = rindex(fname,'.');
  if (cp1 == NULL) cp1 = dest + strlen(dest);
  strncpy(cp1, str, MAX_FNAME_LEN - (cp1 - dest));
  cp1 = dest + strlen(dest);
  if (cp2 != NULL)
    strncpy(cp1, cp2, MAX_FNAME_LEN - (cp1 - dest));

  return 0;
}
