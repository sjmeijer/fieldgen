#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
float fminf(float x, float y);

#include "field_init.h"
#include "point.h"
#include "signal_calc_util.h"
#include "assym_detector.h"

static float zlen;          /*z = 0 => front of detector, zlen => length */
static float max_r;         /* maximum radius */
static float hole_r;        /* radius for central hole*/
static float hole_start_z;  /* z coordinate where central hole starts*/
static float hole_end_curv_r; /*needs better name! hole end curvature radius 
                               at r=rmax*/
static float hole_max_width; /*width of central hole @ back of crystal*/
static float hole_taper_ang; /*angle (deg) for widening of central hole @back*/
static int *seg_no_z;       /*segment number as a function of (integer) z*/
static int *seg_no_phi;     /*segment number as a function of phi in degrees*/
static int nseg_z;          /*number of segments in z direction*/
static int nseg_phi;        /*number of segments in phi direction*/
static float *zmax_segment; /*max z for each segment in z dir*/
static int crystal_type = DEFAULT_CRYSTAL_TYPE;
static int quad_number = 2;/* 2 and onwards have the same orientation*/
static struct point virtual_crystal_corners[N_CRYSTAL_TYPES][2][NCORNERS/2];

extern float xmin, ymin, zmin, xstep, ystep, zstep, xmax, ymax, zmax;
extern int numx, numy, numz; 
extern char ***point_type;
extern float ***fraction_in_det;
static int in_crystal(point pt);
static int init_point_fractions(void);
static float get_fraction(int i, int j, int k, int nsteps);
static int project_to_edges(struct point pt, struct point projections[NCORNERS/2]);
static int segment_number(point pt);

float zlen_detector(void){
  return zlen;
}

int geometry_init(char *geometry_fname){
  FILE *fp;
  char line[MAX_LINE], *cp, *next;
  int i, j, k;
  struct point pt;

  if ((fp = fopen(geometry_fname, "r")) == NULL){
    error("failed to open geometry configuration file: %s\n", geometry_fname);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%d %d", &nseg_z, &nseg_phi) != 2){
    error("Failed to read number of segments from %s\n", geometry_fname);
    return -1;
  }
  printf( "Segments in z-direction: %d, phi: %d\n", nseg_z, nseg_phi);
  if (read_setup_line(fp, line) != 0){
    error("Failed to read segment thickness info from file: %s\n", 
	  geometry_fname);
    return -1;
  }
  if (( zmax_segment = malloc(nseg_z*sizeof(*zmax_segment))) == NULL){
    error("malloc failed in geometry_init\n");
    return -1; /*FIXME -- or exit?*/
  }
  cp = line;
  for (i = 0; i < nseg_z; i++){
    zmax_segment[i] = strtod(cp, &next);
    if (next == cp){
      error("Failed to parse segment thickness info from file: %s\n",
	    geometry_fname);
      return -1;
    }
    cp = next;
  }
  for (i = 1; i < nseg_z; i++){
    zmax_segment[i] += zmax_segment[i-1];
  } 
  printf( "Segment max z, as a fn of segment in z direction: ");
  for (i = 0; i < nseg_z; i++){
    printf( "%.1f ", zmax_segment[i]);
  }
  printf("\n");

  if (read_setup_line(fp, line) != 0
	|| sscanf(line, "%f", &zlen) != 1){
    error("Failed to read crystal length from file: %s\n", geometry_fname);
    return -1;
  }
  printf( "crystal length: %.1f\n", zlen);
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &max_r) != 1){
    error("failed to read crystal radius from file: %s\n", geometry_fname);
    return -1;
  }
  printf( "crystal radius: %.1f\n", max_r);
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &hole_r) != 1){
    error("failed to read hole radius from file: %s\n", geometry_fname);
    return -1;
  }
  printf( "hole radius: %.1f\n", hole_r);
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &hole_start_z) != 1){
    error("failed to read z coordinate for start of central hole from file: %s\n",
	  geometry_fname);
    return -1;
  }
  printf( "central hole starts at z=%.1f\n", hole_start_z);
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f", &hole_end_curv_r) != 1){
    error("Failed to read hole end curvature from file: %s\n", 
	  geometry_fname);
    return -1;
  }
  printf( "The radius of curvature at the bottom "
       "of the central hole is %.1f\n", hole_end_curv_r);
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f %f", &hole_max_width, &hole_taper_ang) != 2){
    error("failed to read hole opening taper parameters from file: %s\n", 
	  geometry_fname);
    return -1;
  }
  printf ("Hole opening taper parameters: %.1f %.1f\n", 
	hole_max_width, hole_taper_ang);
  for (i = 0; i < N_CRYSTAL_TYPES; i++){
    for (j = 0; j < 2; j++){
      for (k = 0; k < NCORNERS/2; k++){
	if (read_setup_line(fp,line) != 0
	    || sscanf(line, "%f %f %f", &pt.x, &pt.y, &pt.z) != 3){
	  error("failed to read crystal corner (%d,%d,%d) from file: %s\n",
		i,j,k,geometry_fname);
	  return -1;
	}
	virtual_crystal_corners[i][j][k] = pt;
      }
    }
  }
  printf( "Succesfully read %d corner positions\n", 
       N_CRYSTAL_TYPES*2*NCORNERS/2);
  fclose(fp);

  if (seg_no_z) free(seg_no_z);
  if (seg_no_phi) free(seg_no_phi);

  init_point_fractions();

  return nseg_phi*nseg_z + 1;
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
	}else{
	  v[0][i][j][k] = v[1][i][j][k] = Vb*(max_r - r)/(max_r - hole_r); 
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
    for (j = 0; j < numy; j++){
      for (k = 0; k < numz; k++){
	v[0][i][j][k] = v[1][i][j][k] = 0.0;
      }
    }
  }

  for (i = 0; i < numx; i++){
    pt.x = xmin + i*xstep;
    for (j = 0; j < numy; j++){
      pt.y = ymin + j*ystep;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < numz; k++){
	pt.z = zmin + k*zstep;
	if (point_type[i][j][k] == CONTACT_VB){
	  if (cno == nseg_z*nseg_phi)
	    v[0][i][j][k] = v[1][i][j][k] = 1.0;
	  else
	    v[0][i][j][k] = v[1][i][j][k] = 0.0;
	}else if (point_type[i][j][k] == CONTACT_0){
	  if (cno < nseg_z*nseg_phi){
	    if (segment_number(pt) == cno)
	      v[0][i][j][k] = v[1][i][j][k] = 1.0;	      
	  }else{
	    v[0][i][j][k] = v[1][i][j][k] = 0.0;
	  }
	}else if (point_type[i][j][k]  == OUTSIDE){
	  v[0][i][j][k] = v[1][i][j][k] = 0.0;	  
	}else{
	  v[0][i][j][k] = v[1][i][j][k] = 0.5;
	}
      }
    }
  }
  printf("wp ... init done\n"); fflush(stdout);
  return 0;


}


/* segment_number
   returns the (geometrical) segment number at point pt, or -1 
   if outside crystal
*/
static int segment_number(point pt){
  int seg_phi, seg_z;
  float thp, thn, th, costh;
  struct point xyvect, edge_p[NCORNERS/2], side_pxy, cp, cpold;
  int i, j, seg_phi_n, seg_phi_p;

  //if (!in_crystal(pt)) return -1;

  xyvect = pt;
  xyvect.z = 0.0;
  project_to_edges(pt, edge_p);
  seg_z = 0;
  /*find segment number in z direction*/
  for (i = 0; i < nseg_z; i++){
    if (pt.z <= zmax_segment[i]){
      seg_z = i;
      break;
    }
  }
  /*find segment number in phi direction*/
  thp = thn = M_PI;
  side_pxy.z = 0;
  seg_phi = seg_phi_p = seg_phi_n = 0;
  for (i = 0; i < NCORNERS/2; i++){
    j = (i+1)%6;
    side_pxy.x = (edge_p[i].x + edge_p[j].x)/2;
    side_pxy.y = (edge_p[i].y + edge_p[j].y)/2;
    costh = dot_prod(side_pxy, xyvect)
      /sqrt(dot_prod(side_pxy, side_pxy)*dot_prod(xyvect, xyvect));
    th = acos(costh);
    cp = cross_prod(side_pxy, xyvect);
    if (i == 0){
      cpold = cp;
    }      
    if (cpold.z*cp.z >= 0){
      if (th <= thp){
	thp = th;
	seg_phi_p = i;
      }
    }else{
      if (th <= thn){
	thn  = th;
	seg_phi_n = i;
      }
    }
  }
  if (seg_phi_p - seg_phi_n == 1) seg_phi = seg_phi_p;
  if (seg_phi_n - seg_phi_p == 1) seg_phi = seg_phi_n;
  if (abs(seg_phi_n - seg_phi_p) == 5) seg_phi = 0;
  //if (sign == 0.0) seg_phi = 0; //what?
  if (crystal_type == CRYSTAL_A)
    seg_phi = (seg_phi + 4) % 6;
  if (quad_number > 1){/*need to rotate B crystal for quad 2 onwards*/
    if (crystal_type == CRYSTAL_B)
      seg_phi = (seg_phi + 5) %6;
  }
  if (seg_phi < 0) seg_phi += 6;
  
  return nseg_phi*seg_z + seg_phi;
}


int geometry_finalize(void){
  if (seg_no_z)   free(seg_no_z);
  if (seg_no_phi) free(seg_no_phi);
  seg_no_z = seg_no_phi = NULL;

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
	  fraction_in_det[i][j][k] = 0.0;
	  continue;
	}
#endif
	f = fraction_in_det[i][j][k] = get_fraction(i,j,k,3); 
	if (f == 1.0) {
	  point_type[i][j][k] = INSIDE;
	}else if (f == 0.0){
	  point_type[i][j][k] = OUTSIDE;
	}else{
	  r = sqrt(pt.x*pt.x + pt.y*pt.y);
	  if (pt.z < zstep){/*front face*/
	    point_type[i][j][k] = CONTACT_0;	    
	  }else if(r < hole_max_width){/*central contact*/
	    point_type[i][j][k] = CONTACT_VB;	
	  }else if (pt.z >= zlen - zstep && r <= max_r - xstep){
	      point_type[i][j][k] = INSIDE;	    
	  }else{/*outside contact*/
	    point_type[i][j][k] = CONTACT_0;	    
	  }
	}
      }
    }
  }
  for (i = 1; i < numx-1; i++){
    for (j = 1; j < numy-1; j++){
      for (k = 1; k < numz-1; k++){
	if (point_type[i][j][k] != INSIDE) continue;
	if (i < numx/2 && point_type[i-1][j][k] == OUTSIDE){
	  point_type[i-1][j][k] = CONTACT_0;
	}
	if (i > numx/2 && point_type[i+1][j][k] == OUTSIDE){
	  point_type[i+1][j][k] = CONTACT_0;
	}
	if (j < numy/2 && point_type[i][j-1][k] == OUTSIDE){
	  point_type[i][j-1][k] = CONTACT_0;
	}
	if (j > numy/2 && point_type[i][j+1][k] == OUTSIDE){
	  point_type[i][j+1][k] = CONTACT_0;
	}
      }
    }
  }

  printf("\n");
  return 0;
}

static float get_fraction(int i, int j, int k, int nsteps){
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
	if (in_crystal(pt)) {
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



#define SQ(x) ((x)*(x))
/*returns 0 (false) or 1 (true) depending on whether pt is inside the crystal*/
static int in_crystal(point pt){
  float r, rmin;
  struct point ep[NCORNERS/2], vij, vpj, cp, cpold;
  int i,j;
  float hz, dotp;
  float rmax, zmin, zmax, hr, htapd, htapl;

  zmin = 0;
  zmax = zlen;
  rmax = max_r;
  hr = hole_r;
  htapd = hole_max_width;

  /*outside crystal in z direction*/
  if (pt.z < zmin || pt.z >= zmax) return 0;

  /*outside maximum radius*/
  r = sqrt(SQ(pt.x) + SQ(pt.y));
  if (r > rmax) return 0;

  /*inside central hole*/
  hz = hole_start_z + hole_end_curv_r;
  if (pt.z >= hz && r < hr) 
    return 0;

  /*rounded end of hole*/
  rmin = 0;
  if (pt.z < hz && pt.z >= hole_start_z){
    rmin = sqrt(SQ(hole_end_curv_r) - SQ(hz - pt.z)) + hr - hole_end_curv_r;
  }
  if (r < rmin) return 0;

  /*widening at hole opening -- only for nonsymm crystal*/
  htapl = htapd/tan(hole_taper_ang*M_PI/180);
  if (pt.z > zlen - htapl){
    rmin = htapd - (zlen - pt.z)*tan(hole_taper_ang*M_PI/180);
    if (r < rmin) return 0;
  }

  /*find points corresponding to same "z" between pairs 
    of virtual corners @ front & back of detector => points at 
    "edges" of crystal*/
  project_to_edges(pt, ep);

  /* is pt inside the hexagon formed by ep?*/
  for (i = 0; i < NCORNERS/2; i++){
    j = (i + 1)%6;
    vij = vector_sub(ep[j], ep[i]);
    vpj = vector_sub(pt, ep[j]);
    cp = cross_prod(vij,vpj);
    if (i == 0){
      cpold = cp;
    }else{
      dotp = dot_prod(cp,cpold);
      if (dotp < 0) return 0; /*outside crystal*/
    }
  }

  return 1;
}
#undef SQ


static int project_to_edges(struct point pt, struct point projections[NCORNERS/2]){
  int i;
  struct point vc1, vc2, dvc;
  float f;

  for (i = 0; i < NCORNERS/2; i++){
    vc1 = virtual_crystal_corners[crystal_type][0][i];
    vc2 = virtual_crystal_corners[crystal_type][1][i];
    f = (pt.z - vc1.z)/(vc2.z - vc1.z);
    dvc = vector_sub(vc2, vc1);
    projections[i] = vector_add(vc1, vector_scale(dvc,f));
  }
  return 0;
}

int set_crystal_geometry(int type){
  if (type < 0 || type >= N_CRYSTAL_TYPES){
    printf("invalid crystal type: %c\n", type + 'A');
    return -1;
  }
  crystal_type = type;
  printf("crystal type set to %c\n", crystal_type + 'A');
  return 0;
}

int set_quad_number(int n){
  if (n < 0){
    printf("quad number %d is invalid\n", n);
    return -1;
  }
  quad_number = n;
  printf("quad number set to %d\n", quad_number);
  return 0;
}
