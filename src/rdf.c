/* 2020-03-16
 * by dianmo
 * Radial distribution function g(r)
 */

#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define rcut 10
#define nbin 500

#define vshell 4 * M_PI * (pow((i + 1) * dbox, 3) - pow(i * dbox, 3)) / 3.

int main(int argc, char *argv[])
{
  FILE *fp;
  double *xyz, gr[nbin] = {0};
  int nthreads = 4, nstep = 1, natom = 43, ipbc = 0, ibox, istep, i, j;
  double dbox, pbc_x = 1, pbc_y = 1, pbc_z = 1;
  double rcut2 = rcut * rcut, dx, dy, dz, dr, rho = 1;

  for (i = 1; i < argc; i += 2)
  {
    if (strcmp(argv[i], "nstep") == 0) {
      nstep = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "natom") == 0) {
      natom = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "pbc") == 0) {
      ipbc = 1;
      pbc_x = atof(argv[i + 1]);
      pbc_y = pbc_x;
      pbc_z = pbc_x;
      rho = natom / (pbc_x * pbc_y * pbc_z); 
    } else if (strcmp(argv[i], "nthreads") == 0) {
      nthreads = atoi(argv[i + 1]);
    }
  }

  if ((fp = fopen("md.out","rb")) == 0)
  {
    printf("Error: file 'md.out' does not exist!\n");
    exit(1);
  }

  xyz = (double *)malloc(sizeof(double) * 3 * natom);
  dbox = (float)rcut / nbin;

  for (istep = 0; istep < nstep; ++istep)
  {
    /* Read coordinate of atoms */
    fread(xyz, sizeof(double), 3 * natom, fp);

    /* Calculate distance */
    omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(guided) \
        private(j, dx, dy, dz, dr, ibox) firstprivate(dbox, ipbc, pbc_x, pbc_y, pbc_z)
    for (i = 0; i < natom - 1; ++i)
    {
      for (j = i + 1; j < natom; ++j)
      {
        dx = *(xyz + 3 * i) - *(xyz + 3 * j);
        dy = *(xyz + 3 * i + 1) - *(xyz + 3 * j + 1);
        dz = *(xyz + 3 * i + 2) - *(xyz + 3 * j + 2);

        /* Periodic boundary condition */
        if (ipbc) {
          if (dx > 0.5) --dx;
          else if (dx < -0.5) ++dx;
          if (dy > 0.5) --dy;
          else if (dy < -0.5) ++dy;
          if (dz > 0.5) --dz;
          else if (dz < -0.5) ++dz;
          dx *= pbc_x;
          dy *= pbc_y;
          dz *= pbc_z;
        }
        dr = dx * dx + dy * dy + dz * dz;
      }

      if (dr > rcut2) continue;
      dr = sqrt(dr);
      ibox = (int)(dr / dbox);
#pragma omp critcal
      ++ *(gr + ibox);
    }
  }
  free(xyz);
  fclose(fp);

  /* Correlation and output */
  fp = fopen("rdf.dat","w");
  for (i = 0; i < nbin; ++i)
  {
    gr[i] /= vshell * nstep * rho;
    fprintf(fp,"%-.3f   %-.6f\n", (i + 0.5) * dbox, gr[i]);
  }
  fclose(fp);

  return 0;
}
