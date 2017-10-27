#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "field_init.h"
#include "signal_calc_util.h"
#include "point.h"

static int init_point_fractions(void);
static float get_fraction(int i, int j, int k, int nsteps);
static float in_detector(point pt);

static float zlen;
static float radius;
static float hole_r;
static float hole_start_z;
extern float xmin, ymin, zmin, xstep, ystep, zstep, xmax, ymax, zmax;
extern int numx, numy, numz; 
extern char ***point_type;

int geometry_init(char *geometry_fn){

  FILE *fp;
  char line[MAX_LINE];
  

  if ((fp = fopen(geometry_fn, "r")) == NULL){
    error("failed to open geometry configuration file: %s\n", geometry_fn);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &zlen) != 1){
    error("Failed to read z length from file: %s\n");
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &radius) != 1){
    error("Failed to read radius from file: %s\n");
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &hole_r) != 1){
    error("Failed to read hole radius from file: %s\n");
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &hole_start_z) != 1){
    error("Failed to read hole start z from file: %s\n");
    fclose (fp);
    return -1;
  }
  fclose(fp);

  init_point_fractions();

  return 2;
}

int init_ev_calc(float ***v[2], float Vb){
  int i, j, k;
  point pt;
  float r;

  for (i = 0; i < numx; i++){
    pt.x = xmin + i*xstep;
    for (j = 0; j < numy; j++){
      pt.y = ymin + j*ystep;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < numz; k++){
	pt.z = zmin + k*zstep;
	if (point_type[i][j][k] == CONTACT_0){
	  v[0][i][j][k] = v[1][i][j][k] = 0.0;
	}else if (point_type[i][j][k] == CONTACT_VB){
	  v[0][i][j][k] = v[1][i][j][k] = Vb;
	}else if (point_type[i][j][k]  == OUTSIDE){
	  v[0][i][j][k] = v[1][i][j][k] = 0.0;	  
#ifdef CALC_VFIELD
	}else if (point_type[i][j][k]  == VACUUM){
	  v[0][i][j][k] = v[1][i][j][k] = 0.0;	  
#endif
	}else{
	  if (pt.z > hole_start_z)
	    v[0][i][j][k] = v[1][i][j][k] = 
	      Vb*(r - hole_r)/(radius - hole_r); 
	  else
	    v[0][i][j][k] = v[1][i][j][k] = 
	      Vb*(r+1)/radius;
	  if (v[0][i][j][k] >= Vb*0.95) v[0][i][j][k] = Vb*0.95;
	  if (v[0][i][j][k] <= Vb*0.05) v[0][i][j][k] = Vb*0.05;
	}
      }
    }
  }
  printf("ev ... init done\n"); fflush(stdout);
  return 0;
}

int init_wp_calc(float ***v[2], int cno){
  int i, j, k;
  point pt;
  float r;

  for (i = 0; i < numx; i++){
    pt.x = xmin + i*xstep;
    for (j = 0; j < numy; j++){
      pt.y = ymin + j*ystep;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < numz; k++){
	pt.z = zmin + k*zstep;
	if (point_type[i][j][k] == CONTACT_0){
	  if (cno == 1)
	    v[0][i][j][k] = v[1][i][j][k] = 1.0;
	  else
	    v[0][i][j][k] = v[1][i][j][k] = 0.0;
	}else if (point_type[i][j][k] == CONTACT_VB){
	  if (cno == 1)
	    v[0][i][j][k] = v[1][i][j][k] = 0.0;
	  else
	    v[0][i][j][k] = v[1][i][j][k] = 1.0;
	}else if (point_type[i][j][k]  == OUTSIDE){
	  v[0][i][j][k] = v[1][i][j][k] = 0.0;	  
#ifdef CALC_VFIELD
	}else if (point_type[i][j][k]  == VACUUM){
	  v[0][i][j][k] = v[1][i][j][k] = 0.0;	  
#endif
	}else{
	  if (pt.z > hole_start_z)
	    if (cno == 1)
	      v[0][i][j][k] = v[1][i][j][k] = 
		(radius -r)/(radius - hole_r); 
	    else
	      v[0][i][j][k] = v[1][i][j][k] = 
		(r - hole_r)/(radius - hole_r); 
	  else
	    if (cno == 1)
	      v[0][i][j][k] = v[1][i][j][k]  = (radius-r)/radius;
	    else
	      v[0][i][j][k] = v[1][i][j][k] = r/radius;
	  if (v[0][i][j][k] >= 0.95) v[0][i][j][k] = 0.95;
	  if (v[0][i][j][k] <= 0.05) v[0][i][j][k] = 0.05;
	}
      }
    }
  }
  printf("wp ... init done\n"); fflush(stdout);
  return 0;
}

int geometry_finalize(void){
  return 0;
}


static int init_point_fractions(void){
  int i, j, k;
  point pt;;
  float r, f;

  for (i = 0; i < numx; i++){
    printf("\r %d/%d", i, numx-1);fflush(stdout);
    pt.x = xmin + i*xstep;
    for (j = 0; j < numy; j++){
      pt.y = ymin + j*ystep;
      for (k = 0; k < numz; k++){
	pt.z = zmin + k*zstep;
#ifdef CALC_VFIELD
	if (i == 0 || j == 0 || k == 0 
	    || i == numx -1 || j == numy -1 || k == numz -1){
	  point_type[i][j][k] = CONTACT_0;
	  continue;
	}
#endif
	f = get_fraction(i,j,k,1); 
	if (f == 1.0) {
	  point_type[i][j][k] = INSIDE;
	}else if (f == 0.0){
#ifdef CALC_VFIELD
	  point_type[i][j][k] = VACUUM;
#else
	  point_type[i][j][k] = OUTSIDE;
#endif
	}else{
	  r = sqrt(pt.x*pt.x + pt.y*pt.y);
	  if (pt.z < zstep){
	    point_type[i][j][k] = CONTACT_VB;	    
	  }else if (radius - r < xstep){
	    point_type[i][j][k] = CONTACT_VB;	    	    
	  }else if (r - hole_r < xstep){
	    point_type[i][j][k] = CONTACT_0;
	  }else if (zlen - pt.z < zstep && r > hole_r){
#ifdef CALC_VFIELD
	    point_type[i][j][k] = INSIDE;	    	    
#else
	    if (point_type[i][j][k-1] == Z_SURF)
	      point_type[i][j][k] = OUTSIDE;
	    else
	      point_type[i][j][k] = Z_SURF;
#endif
	  }else{
	    error("point %d %d %d not accounted for!\n", i,j ,k);
	  }
	}
      }
    }
  }
  printf("\n");
  return 0;
}

float get_fraction(int i, int j, int k, int nsteps){
  int ii, jj, kk;
  point pt;
  int n, m; 
  float fraction_in_det;

  n = m = 0;
  fraction_in_det = 0.0;
  for (ii = -nsteps; ii <= nsteps; ii+=1){
    pt.x = xmin + (i + ii/2.0/nsteps)*xstep;
    if (pt.x < xmin || pt.x > xmax){
      n++;
      continue;
    }
    for (jj = -nsteps; jj <=nsteps; jj+=1){
      pt.y = ymin + (j + jj/2.0/nsteps)*ystep;
      if (pt.y < ymin || pt.y > ymax){
	n++;
	continue;
      }
      for (kk = -nsteps; kk <= nsteps; kk+=1){
	pt.z = zmin + (k + kk/2.0/nsteps)*zstep;
	if (pt.z < zmin || pt.z > zmax){
	  n++;
	  continue;
	}
	if (in_detector(pt)) {
	  m++;
	  fraction_in_det += 1;
	}
	n++;
      }
    }
  }
  if (m == n) return 1.0;
  return fraction_in_det / n;
}

static float in_detector(point pt){
  float r, z;
  
  r = sqrt(pt.x*pt.x + pt.y*pt.y);
  z = pt.z;

  if (z < 0 || z > zlen || r > radius)
    return 0;
  if (r < hole_r && z > hole_start_z)
    return 0;
  
  return 1;
}
