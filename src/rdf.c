/* 2020-03-16
 * by dianmo
 * Radial distribution function g(r)
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define gbox 8
#define nbox 500

#define vshell 4 * M_PI * (pow((i + 1) * dbox, 3) - pow(i * dbox, 3)) / 3.

int main(int argc, char *argv[])
{
  FILE *fp;
  double *xyz;
  int nstep = 1, natom = 43, ipbc = 0, ibox, istep, i, j;
  double gr[nbox] = {0};
  double dbox, dx, dy, dz, dr, box2, pbc_box, scale = 1, rho = 1;

  for (i = 1; i < argc; i += 2)
  {
    if (strcmp(argv[i], "nstep") == 0) {
      nstep = atoi(argv[i+1]);
    } else if (strcmp(argv[i], "natom") == 0) {
      natom = atoi(argv[i+1]);
    } else if (strcmp(argv[i], "pbc") == 0) {
      ipbc = 1;
      pbc_box = atof(argv[i+1]);
      scale = 3 * pbc_box * pbc_box;  
      rho = natom / pow(pbc_box, 3);
    }
  }

  if ((fp = fopen("md.out","rb")) == 0)
  {
    printf("Error: file 'md.out' does not exist!\n");
    exit(1);
  }
  xyz = (double *)malloc(sizeof(double) * 3 * natom);
  dbox = (float)gbox / nbox;
  box2 = gbox * gbox / scale;

  for (istep = 0; istep < nstep; ++istep)
  {
    /* Read coordinate of atoms */
    fread(xyz, sizeof(double), 3 * natom, fp);

    /* Calculate distance */
    for (i = 0; i < natom - 1; ++i)
    {
      for (j = i + 1; j < natom; ++j)
      {
        dx = *(xyz + 3 * i) - *(xyz + 3 * j);
        dy = *(xyz + 3 * i + 1) - *(xyz + 3 * j + 1); 
        dz = *(xyz + 3 * i + 2) - *(xyz + 3 * j + 2);

        /* Periodic boundary condition */
        if (ipbc) {
          if (dx > 0.5) dx -= 1;
          else if (dx < -0.5) dx += 1; 
          if (dy > 0.5) dy -= 1;
          else if (dy < -0.5) dy += 1;
          if (dz > 0.5) dz -= 1;
          else if (dz < -0.5 * pbc_box) dz += 1; 
        }

        dr = dx * dx + dy * dy + dz * dz;
      }

      if (dr > box2) continue;
      dr = sqrt(dr * scale);
      ibox = (int)(dr / dbox);
      ++gr[ibox];
    }
  }
  free(xyz);
  fclose(fp);

  /* Correlation and output */
  fp = fopen("rdf.dat","w");
  for (i = 0; i < nbox; ++i)
  {
    gr[i] /= vshell * nstep * rho;
    fprintf(fp,"%-.3f   %-.6f\n", (i + 0.5) * dbox, gr[i]);
  }
  fclose(fp);

  return 0;
}
