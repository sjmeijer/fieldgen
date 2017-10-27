#include <stdio.h>
#include <stdlib.h>
float strtof(const char *nptr, char **endptr);
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
static float hole_tap_r; /*addition to r at z=zlen by tapering*/
static float contact_r;
static float contact_l;
static int nseg_phi;
static float seg_start_r;
static float ditch_depth;
static float ditch_width;
static int nseg_z;
static float *seg_no_z, *seg_no_phi;
static int ncontacts;
static float rear_seg_r;

extern float xmin, ymin, zmin, xstep, ystep, zstep, xmax, ymax, zmax;
extern int numx, numy, numz; 
extern char ***point_type;

int geometry_init(char *geometry_fn){
  int i,j;
  FILE *fp;
  char line[MAX_LINE], *ep, *cp;
  float *zmax_segment;

  if ((fp = fopen(geometry_fn, "r")) == NULL){
    error("failed to open geometry configuration file: %s\n", geometry_fn);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &zlen) != 1){
    error("Failed to read z length from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &radius) != 1){
    error("Failed to read radius from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &hole_r) != 1){
    error("Failed to read hole radius from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &hole_start_z) != 1){
    error("Failed to read hole start z from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &hole_tap_r) != 1){
    error("Failed to read hole taper r from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &contact_r) != 1){
    error("Failed to read central contact r from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &contact_l) != 1){
    error("Failed to read contact depth from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%d", &nseg_phi) != 1){
    error("Failed to read # azimutal segments from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if ((seg_no_phi = malloc(360*sizeof(*seg_no_phi))) == NULL){
    error("malloc failed\n");
    exit(1);
  }
  for (i = 0; i < 360; i++){
    seg_no_phi[i] = (int)(nseg_phi*i/360.0);
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &seg_start_r) != 1){
    error("Failed to read azimutal segmentation inner r from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &ditch_depth) != 1){
    error("Failed to read front ditch depth from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &ditch_width) != 1){
    error("Failed to read front ditch width from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }
  if ((zmax_segment = malloc(zlen*sizeof(*zmax_segment))) == NULL){
    error("Malloc failed\n", geometry_fn);
    exit(1);
  }
  if (read_setup_line(fp, line) != 0){
    error("Failed to read z segmentation from file: %s\n", geometry_fn);
    fclose(fp);
    return -1;
  }
  cp = line;
  for (nseg_z = 0; nseg_z < zlen; nseg_z++){
    zmax_segment[nseg_z] = strtof(cp, &ep);
    if (ep == cp) break;
    if (zmax_segment[nseg_z] <= 0 || zmax_segment[nseg_z] > zlen){
      error("invalid segment width: %f\n",zmax_segment[nseg_z]);
      fclose(fp);
      return -1;
    }
    cp = ep;
  }
  if (nseg_z == zlen){
    error("too many segments in z direction\n");
    fclose(fp);
    return -1;
  }
  for (i = 1; i < nseg_z; i++) zmax_segment[i] += zmax_segment[i-1];
  if (zmax_segment[nseg_z-1] > zlen){
    error("sum of segment thicknesses is larger than detector size\n");
    fclose(fp);
    return -1;
  }
  if ((seg_no_z = malloc((int)(zlen+1)*sizeof(*seg_no_z))) == NULL){
    error("Realloc failed\n");
    exit(1);
  }
  j = 0;
  for (i = 0; i < nseg_z; i++){
    for (j; j < zmax_segment[i] && j < zlen+1; j++){
      seg_no_z[j] = i;
    }
  }
  for (j; j < zlen+1; j++){
    seg_no_z[j] = nseg_z-1;
  }

  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &rear_seg_r) != 1){
    error("Failed to read rear segmentation r from file: %s\n",geometry_fn);
    fclose (fp);
    return -1;
  }

  ncontacts = nseg_phi + nseg_z + 3;
  if ( rear_seg_r > 0 )
    ncontacts += 1;
  printf("contacts: z %d phi %d\n", nseg_z, nseg_phi);
  printf("ncontacts: %d\n", ncontacts);
  free(zmax_segment);


  fclose(fp);

  init_point_fractions();

  return ncontacts;
}

int init_ev_calc(float ***v[2], float Vb){
  int i, j, k;
  point pt;
  float r;

#ifdef FRONT_WP_HACK
  for (i = 0; i < numx; i++){
    pt.x = xmin + i*xstep;
    for (j = 0; j < numy; j++){
      pt.y = ymin + j*ystep;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < numz; k++){
	pt.z = zmin + k*zstep;
	if (pt.z < zstep && r <= seg_start_r 
	    && point_type[i][j][k] == CONTACT_0){
#ifdef CALC_VFIELD
	    point_type[i][j][k] = INSIDE;
#else
	    point_type[i][j][k] = Z_SURF;
#endif
	}
      }
    }
  }
#endif

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
	  v[0][i][j][k] = v[1][i][j][k] = 
	    Vb*(1- fabsf(hole_r + (radius- hole_r)/2 -r)/radius*2);
	  if (Vb > 0){
	    if (v[0][i][j][k] >= Vb*0.95) v[0][i][j][k] = Vb*0.95;
	    if (v[0][i][j][k] <= Vb*0.05) v[0][i][j][k] = Vb*0.05;
	  }else{
	    if (v[0][i][j][k] <= Vb*0.95) v[0][i][j][k] = Vb*0.95;
	    if (v[0][i][j][k] >= Vb*0.05) v[0][i][j][k] = Vb*0.05;	  
	  }
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
  float r, phi;


#ifdef FRONT_WP_HACK
  for (i = 0; i < numx; i++){
    pt.x = xmin + i*xstep;
    for (j = 0; j < numy; j++){
      pt.y = ymin + j*ystep;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < numz; k++){
	pt.z = zmin + k*zstep;
	if (pt.z < zstep && r <= seg_start_r 
	    && point_type[i][j][k] == CONTACT_0){
#ifdef CALC_VFIELD
	    point_type[i][j][k] = INSIDE;
#else
	    point_type[i][j][k] = Z_SURF;
#endif
	}
      }
    }
  }
#endif

  for (i = 0; i < numx; i++){
    pt.x = xmin + i*xstep;
    for (j = 0; j < numy; j++){
      pt.y = ymin + j*ystep;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < numz; k++){
	pt.z = zmin + k*zstep;
	if (point_type[i][j][k] == CONTACT_VB){
	  if (cno == ncontacts -1)
	    v[0][i][j][k] = v[1][i][j][k] = 1.0;
	  else
	    v[0][i][j][k] = v[1][i][j][k] = 0.0;
	}else if (point_type[i][j][k] == CONTACT_0){	  
	  if (cno == ncontacts -1){
	    v[0][i][j][k] = v[1][i][j][k] = 0.0;
	  }else if (cno < nseg_phi){// phi segment
	    if (pt.x != 0.0){
	      phi = atan(pt.y/pt.x); /*-PI/2 -- PI/2*/
	      if (pt.x < 0) phi += M_PI;
	      if (pt.x > 0 && pt.y < 0) phi += 2*M_PI;
	    }else{
	      if (pt.y >= 0) phi = M_PI/2;
	      else phi = 3*M_PI/2;
	    }
	    phi = phi/M_PI*180.0;	    
	    if (pt.z < zstep && r > seg_start_r && seg_no_phi[(int)phi] == cno){
		v[0][i][j][k] = v[1][i][j][k] = 1.0;
	    }else{
	      v[0][i][j][k] = v[1][i][j][k] = 0.0;
	    }
	  }else if (cno < nseg_phi + nseg_z){// z segment
	    if (fabsf(r - radius) < xstep 
		&& pt.z >= zstep
		&& pt.z <= zlen - zstep
		&& seg_no_z[(int)pt.z] +nseg_phi == cno){
	      v[0][i][j][k] = v[1][i][j][k] = 1.0;	      
	    }else{
	      v[0][i][j][k] = v[1][i][j][k] = 0.0;	      
	    }
	  }else if (cno == nseg_phi + nseg_z){// back of det.
	    if (pt.z > zlen - zstep && 
		(rear_seg_r == 0 || (r > rear_seg_r && r <= radius))){
	      v[0][i][j][k] = v[1][i][j][k] = 1.0;	      
	    }else{
	      v[0][i][j][k] = v[1][i][j][k] = 0.0;	      
	    }
	  }else if (rear_seg_r > 0 && cno == nseg_phi + nseg_z + 1){ 
	    if (pt.z > zlen - zstep && r <= rear_seg_r && r > hole_r){
	      v[0][i][j][k] = v[1][i][j][k] = 1.0;	      
	    }else{
	      v[0][i][j][k] = v[1][i][j][k] = 0.0;	      
	    }
	  }else if (r - hole_r < xstep && //central contact
		    (pt.z <= zlen - zstep)){
	    v[0][i][j][k] = v[1][i][j][k] = 1.0;	      
	  }else{
	    v[0][i][j][k] = v[1][i][j][k] = 0.0;	      	  
	  }
	}else if (point_type[i][j][k]  == OUTSIDE){
	  v[0][i][j][k] = v[1][i][j][k] = 0.0;	  
#ifdef CALC_VFIELD
	}else if (point_type[i][j][k]  == VACUUM){
	  v[0][i][j][k] = v[1][i][j][k] = 0.0;	  
#endif
	}else{
	  v[0][i][j][k] = v[1][i][j][k] = 0.5;	  	  
	}
      }
    }
  }

#ifdef FRONT_WP_HACK
  for (i = 0; i < numx; i++){
    pt.x = xmin + i*xstep;
    for (j = 0; j < numy; j++){
      pt.y = ymin + j*ystep;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < numz; k++){
	pt.z = zmin + k*zstep;
	if (point_type[i][j][k] != INSIDE
	    && point_type[i][j][k] != Z_SURF)
	  continue;
	  if (pt.x != 0.0){
	    phi = atan(pt.y/pt.x); /*-PI/2 -- PI/2*/
	    if (pt.x < 0) phi += M_PI;
	    if (pt.x > 0 && pt.y < 0) phi += 2*M_PI;
	  }else{
	    if (pt.y >= 0) phi = M_PI/2;
	    else phi = 3*M_PI/2;
	  }
	  phi = phi/M_PI*180.0;	    
	  if (pt.z < zstep && seg_no_phi[(int)phi] == cno
	      && cno < nseg_phi){	    
	    if ((i > numx/2 && point_type[i+2][j][k] == CONTACT_0)
		|| (i < numx/2 && point_type[i-2][j][k] == CONTACT_0)
		||(j > numy/2 && point_type[i][j+2][k] == CONTACT_0)
		|| (j < numy/2 && point_type[i][j-2][k] == CONTACT_0)){
	      v[0][i][j][k] = v[1][i][j][k] = 1.0;	  
	      point_type[i][j][k] = CONTACT_0;
	    }
	  }
      }
    }
  }
#endif
  return 0;
}

int geometry_finalize(void){
  return 0;
}


static int init_point_fractions(void){
  int i, j, k;
  point pt;;
  float r, f;
  float hr;

  hr = hole_r;

  for (i = 0; i < numx; i++){
    printf("\r %d/%d", i, numx-1);fflush(stdout);
    pt.x = xmin + i*xstep;
    for (j = 0; j < numy; j++){
      pt.y = ymin + j*ystep;
      for (k = 0; k < numz; k++){
	pt.z = zmin + k*zstep;
	if (pt.z > hole_start_z)
	  hr = hole_r + hole_tap_r*(pt.z - hole_start_z)/(zlen - hole_start_z);
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
	  if(pt.z >= hole_start_z && pt.z <= zlen && r <= hr)
	    point_type[i][j][k] = CONTACT_0;
	  else
#ifdef CALC_VFIELD
	    point_type[i][j][k] = VACUUM;
#else
	    point_type[i][j][k] = OUTSIDE;
#endif
	}else{
	  f = get_fraction(i,j,k, 2);
	  r = sqrt(pt.x*pt.x + pt.y*pt.y);
	  if (pt.z <= contact_l && r -contact_r< xstep){
	    point_type[i][j][k] = CONTACT_VB;
	  }else if (ditch_depth > 0 
		    && r >= seg_start_r - ditch_width 
		    && r < seg_start_r 
		    && pt.z - zstep <= ditch_depth){
	    point_type[i][j][k] = VACUUM;
	  }else if (pt.z < zstep && r < seg_start_r){
#ifdef CALC_VFIELD
	    point_type[i][j][k] = INSIDE;
#else
	    if (k > 0 && point_type[i][j][k-1] == Z_SURF)
	      point_type[i][j][k] = INSIDE;
	    else
	      point_type[i][j][k] = Z_SURF;
#endif
	  }else if (pt.z < zstep){
	    point_type[i][j][k] = CONTACT_0;
	  }else if (radius - r < xstep){
	    point_type[i][j][k] = CONTACT_0;	    	    
	  }else if (zlen - pt.z < zstep && r > hr){	    
	    point_type[i][j][k] = CONTACT_0;	    	    
	  }else if (r - hr < xstep){
	    if (f >= 0.5)
	      point_type[i][j][k] = INSIDE;
	    else
	      point_type[i][j][k] = CONTACT_0;
	  }else if (ditch_depth > 0 
		    && r + xstep >= seg_start_r - ditch_width 
		    && r - xstep < seg_start_r 
		    && pt.z - zstep <= ditch_depth){
	    point_type[i][j][k] = INSIDE;
	  }else{
	    error("point %d %d %d (%.2f %.2f %.2f) not accounted for!\n", 
		  i,j ,k, pt.x, pt.y, pt.z);
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
  float hr;

  r = sqrt(pt.x*pt.x + pt.y*pt.y);
  z = pt.z;
  if (z > hole_start_z)
    hr = hole_r + hole_tap_r*(z - hole_start_z)/(zlen - hole_start_z);
  else
    hr = hole_r;

  if (z < 0 || z > zlen || r > radius)
    return 0;
  if (z < contact_l && r < contact_r)
    return 0;
  if (r < hr && z > hole_start_z)
    return 0;
  if (ditch_width > 0 && ditch_depth > 0)
    if (r < seg_start_r && r >= seg_start_r - ditch_width 
	&& z <= ditch_depth)
      return 0;

  return 1;
}
