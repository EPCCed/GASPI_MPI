/*****************************************************************************
 *
 *  bath.c
 *
 *  Compile and link with -lm for the maths library. 
 *
 *  This utility produces an output file suitable for initialising
 *  a capillary structure within a bath in Ludwig. 
 *
 *  See also capillary.c
 *
 *  $Id: bath.c,v 1.2 2010-08-11 08:21:39 kevin Exp $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2008 The University of Edinburgh
 *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../src/site_map.h"

/* SYSTEM SIZE */
/* Set the system size as desired. Clearly, this must match the system
 * set in the main input file for Ludwig. */

const int xmax = 20;
const int ymax = 20;
const int zmax = 20;

/* CROSS SECTION and DIAMETER */
/* You can choose a square or circular cross section */

enum {CIRCLE, SQUARE};
const int xsection = CIRCLE;
const int diameter = 14;

/* FREE ENERGY PARAMETERS */
/* Set the fluid and solid free energy parameters. The fluid parameters
 * must match those used in the main calculation. See Desplat et al.
 * Comp. Phys. Comm. (2001) for details. */

const double kappa = 0.04;
const double B = 0.0625;
const double H = 0.016;

/* CAPILLARY IMMERSED IN BATH. */
/* Sets the extent of the capillary section. There is a gap at
 * each end. All the tube wets in the same way, inside and outside. */

const int z1 = 10;    /* both inclusive */
const int z2 = 12;

/* OUTPUT */
/* You can generate a file with solid/fluid status information only,
 * or one which includes the wetting parameter H. */

enum {STATUS_ONLY, STATUS_WITH_H};
const int output_type = STATUS_WITH_H;

/* OUTPUT FILENAME */

const char * filename = "bath.001-001";


static void result_height(const char *);

/*****************************************************************************
 *
 *  main program
 *
 *****************************************************************************/

int main(int argc, char ** argv) {

  char * map_in;
  FILE * fp_orig;
  int i, j, k, n;
  int i0, j0;
  int nsolid = 0;

  double rc = 0.5*diameter;
  double x0 = 0.5*xmax + 0.5;
  double y0 = 0.5*ymax + 0.5;

  double map_h;
  double x, y, r;
  double h, h1, theta;

  assert(z1 < zmax);
  assert(z2 < zmax);

  if (argc == 2) result_height(argv[1]);

  printf("Free energy parameters:\n");
  printf("free energy parameter kappa = %f\n", kappa);
  printf("free energy parameter B     = %f\n", B);
  printf("surface free energy   H     = %f\n", H);
  h = H*sqrt(1.0/(kappa*B));
  printf("dimensionless parameter h   = %f\n", h);
  h1 = 0.5*(-pow(1.0 - h, 1.5) + pow(1.0 + h, 1.5));
  theta = acos(h1);
  printf("contact angle theta         = %f radians\n", theta);
  theta = theta*180.0/(4.0*atan(1.0));
  printf("                            = %f degrees\n", theta);


  map_in = (char *) malloc(xmax*ymax*zmax*sizeof(char));
  if (map_in == NULL) exit(-1);

  /* Initialise all points fluid */

  for (i = 0; i < xmax; i++) {
    for (j = 0; j < ymax; j++) {
      for (k = 0; k < zmax; k++) {
	n = ymax*zmax*i + zmax*j + k;
	map_in[n] = FLUID;
      }
    }
  }

  if (xsection == CIRCLE) {

    for (i = 0; i < xmax; i++) {
      x = 1.0 + i - x0;
      for (j = 0; j < ymax; j++) {
        y = 1.0 + j - y0;
        for (k = z1; k <= z2; k++) {
          n = ymax*zmax*i + zmax*j + k;

          map_in[n] = BOUNDARY;
          r = sqrt(x*x + y*y);
          if (r <= rc || r > (rc + sqrt(2.0))) {
            map_in[n] = FLUID;
          }
	}
      }
    }

  }
  else {
    /* Square */

    i0 = (xmax - diameter - 2) / 2;
    j0 = (ymax - diameter - 2) / 2;

    /* i = i0 and i = i0 + diameter + 1 */
 
    for (j = j0; j <= j0 + diameter + 1; j++) {
      for (k = z1; k <= z2; k++) {
	n = ymax*zmax*i0 + zmax*j + k;
	map_in[n] = BOUNDARY;
	n = ymax*zmax*(i0 + diameter + 1) + zmax*j + k;
	map_in[n] = BOUNDARY;
      }
    }

    /* j = j0 and j = j0 + diameter + 1 */

    for (i = i0; i <= i0 + diameter + 1; i++) {
      for (k = z1; k <= z2; k++) {
	n = ymax*zmax*i + zmax*j0 + k;
	map_in[n] = BOUNDARY;
	n = ymax*zmax*i + zmax*(j0 + diameter + 1) + k;
	map_in[n] = BOUNDARY;
      }
    }

  }

  /* picture */

  printf("\nCross section (%d = fluid, %d = solid)\n", FLUID, BOUNDARY);

  k = z1;
  for (i = 0; i < xmax; i++) {
    for (j = 0; j < ymax; j++) {
	n = ymax*zmax*i + zmax*j + k;
      
	if (map_in[n] == BOUNDARY) printf(" %d", BOUNDARY);
	if (map_in[n] == FLUID)    printf(" %d", FLUID);
    }
    printf("\n");
  }

  printf("Upways:\n");

  i0 = xmax/2;
  for (k = zmax-1; k >= 0; k--) {
    for (j = 0; j < ymax; j++) {
	n = ymax*zmax*i0 + zmax*j + k;
      
	if (map_in[n] == BOUNDARY) printf(" %d", BOUNDARY);
	if (map_in[n] == FLUID)    printf(" %d", FLUID);
    }
    printf("\n");
  }

  /* Write new data as char */

  fp_orig = fopen(filename, "w");
  if (fp_orig == NULL) {
    printf("Cant open output\n");
    exit(-1);
  }

  for (i = 0; i < xmax; i++) {
    for (j = 0; j < ymax; j++) {
      for (k = 0; k < zmax; k++) {
	n = ymax*zmax*i + zmax*j + k;

	map_h = 0.0;
	if (map_in[n] == BOUNDARY) {
	  map_h = H;
	  nsolid++;
	}

	/* Put in non-wetting bottom wall */
	if (k == 0) {
	  map_in[n] = BOUNDARY;
	  nsolid++;
	}
	fputc(map_in[n], fp_orig);

	if (output_type == STATUS_WITH_H) {
	  fwrite(&map_h, sizeof(double), 1, fp_orig);
	}
      }
    }
  }

  fclose(fp_orig);

  printf("n = %d nsolid = %d nfluid = %d\n", xmax*ymax*zmax, nsolid,
	 xmax*ymax*zmax - nsolid);

  free(map_in);

  return 0;
}

/*****************************************************************************
 *
 *  result_height
 *
 *  This is to read a file containing order parameter results
 *  and compute the heights inside and outside the tube.
 *
 *****************************************************************************/

void result_height(const char * filename) {

  int * n_in;
  int * n_out;
  int i, j, k;
  double * phi_in;
  double * phi_out;
  double r, x, y;
  double phi;
  double hplus, dh;

  double rc = 0.5*diameter;
  double x0 = 0.5*xmax + 0.5;
  double y0 = 0.5*ymax + 0.5;
  FILE * fp_phi;

  /* Allocate z-wise arrays and initialise */

  n_in = (int *) calloc(zmax, sizeof(int));
  assert(n_in);
  n_out = (int *) calloc(zmax, sizeof(int));
  assert(n_out);
  phi_in = (double *) malloc(zmax*sizeof(double));
  assert(phi_in);
  phi_out = (double *) malloc(zmax*sizeof(double));
  assert(phi_out);

  for (k = 0; k < zmax; k++) {
    n_in[k] = 0;
    n_out[k] = 0;
    phi_in[k] = 0.0;
    phi_out[k] = 0.0;
  }

  /* Read the file */

  fp_phi = fopen(filename, "r");
  if (fp_phi == NULL) {
    printf("fopen(%s) failed\n", filename);
    exit(-1);
  }

  if (xsection == CIRCLE) {
    for (i = 0; i < xmax; i++) {
      x = 1.0 + i - x0;
      for (j = 0; j < ymax; j++) {
        y = 1.0 + j - y0;
        for (k = 0; k < zmax; k++) {

	  fread(&phi, sizeof(double), 1, fp_phi);

          r = sqrt(x*x + y*y);
          if (r <= rc) {
            /* inside */
	    phi_in[k] += phi;
	    n_in[k] += 1;
          }
          if (r > (rc + sqrt(2.0))) {
            /* outside */
	    phi_out[k] += phi;
	    n_out[k] += 1;
          }

	}
      }
    }
  }

  fclose(fp_phi);

  /* Find the mean order parameter profile, and locate the position
   * of the zero to get the heights inside and out */

  for (k = z1; k <= z2; k++) {
    phi_in[k] = phi_in[k] / n_in[k];
    phi_out[k] = phi_out[k] / n_out[k];
  }

  /* Now, with phi = -1 at the bottom of the tube. We need the height
   * from the bottom of the tube to the outside surface (h_+) and the
   * height from the outside surface to the inside surface */

  for (k = z1; k <= z2; k++) {
    if (phi_out[k] < 0.0 && phi_out[k+1] > 0.0) {
      /* Linear interpolation to get surface position */
      dh = phi_out[k] / (phi_out[k] - phi_out[k+1]);
      hplus = 1.0 + (k - z1) + dh;
      printf("Outside (tube base -> surface, h+) = %f\n", hplus);
    }
  }

  for (k = z1; k <= z2; k++) {
    if (phi_in[k] < 0.0 && phi_in[k+1] > 0.0) {
      dh = phi_in[k] / (phi_in[k] - phi_in[k+1]);
      dh = 1.0 + (k - z1) - hplus + dh;
      printf("Inside (from outside surface, h) = %f\n", dh);
    }
  }

  free(phi_in);
  free(phi_out);
  free(n_in);
  free(n_out);
  exit(0);

  return;
}
