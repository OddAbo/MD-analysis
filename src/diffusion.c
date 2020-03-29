/* 2020-03-24
 * by dianmo
 * Mean-squared displacement
 * Velocity auto-correlation function
 */

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define msd_time 20
#define omega_max 20
#define WrapPBC             \
{ if (dx >= 0.5) --dx;      \
  else if (dx < -0.5) ++dx; \
  if (dy >= 0.5) --dy;      \
  else if (dy < -0.5) ++dy; \
  if (dz >= 0.5) --dz;      \
  else if (dz < -0.5) ++dz; \
  dx *= pbc_x, dy *= pbc_y, dz *= pbc_z; }

int nstep, natom, ipbc, nthread, msd_step, omega_step, ithread, i, j;
double dx, dy, dz, dr, dt, pbc_x, pbc_y, pbc_z;
double *xyz, *xyz_new, *msd, *msd_tmp, *omega;

int main(int argc, char* argv[])
{
  FILE *fp;

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
    } else if (strcmp(argv[i], "nthread") == 0) {
      nthread = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "dt") == 0) {
      dt = atof(argv[i + 1]);
    }
  }

  msd_step = msd_time / dt;

  if ((fp = fopen("md.out","rb")) == NULL)
  {
    printf("Error: file 'md.out' does not exist!\n");
    exit(1);
  }

  if (msd_step > nstep)
  {
    printf("Error: value of 'nstep' is too small!\n");
    exit(1);
  }

  xyz = (double *)malloc(sizeof(double) * 3 * natom);
  xyz_new = (double *)malloc(sizeof(double) * 3 * natom);
  msd = (double *)malloc(sizeof(double) * msd_step);
  msd_tmp = (double *)malloc(sizeof(double) * nthread);
  omega = (double *)malloc(sizeof(double) * omega_step);
  memset(msd, 0, sizeof(double) * msd_step);

  for (i = 0; i < (int)(nstep / msd_step); ++i)
  {
    fread(xyz, sizeof(double), 3 * natom, fp);
    for(j = 1; j < msd_step; ++j)
    {
      fread(xyz_new, sizeof(double), 3 * natom, fp);
      memset(msd_tmp, 0, sizeof(double) * nthread);

#pragma omp parallel for schedule(dynamic) \
        private(dx, dy, dz, dr, ithread) \
        firstprivate(ipbc, pbc_x, pbc_y, pbc_z)
      for (int iatom = 0; iatom < natom; ++iatom)
      {
        ithread = omp_get_thread_num();
        dx = *(xyz + 3 * iatom) - *(xyz_new + 3 * iatom);
        dy = *(xyz + 3 * iatom + 1) - *(xyz_new + 3 * iatom + 1);
        dz = *(xyz + 3 * iatom + 2) - *(xyz_new + 3 * iatom + 2);
        if (ipbc) WrapPBC;
        dr = sqrt(dx * dx + dy * dy + dz * dz);
        *(msd_tmp + ithread) += dr;
      }

      for (ithread = 0; ithread < nthread; ++ithread)
      {
        *(msd + j) += *(msd_tmp + ithread);
      }
    }
  }

  for (i = 1; i < msd_step; ++i)
  {
    *()
  }













  return 0;
}
