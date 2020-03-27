/* 2020-03-16
 * by dianmo
 * Radial distribution function g(r)
 */

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define rcut 8
#define nbin 500
#define list_max 500
#define vshell 4 * M_PI * (pow((i + 1) * dbox, 3) - pow(i * dbox, 3)) / 3.
#define rlistshell 0.3 * rcut

#define WrapPBC             \
{ if (dx >= 0.5) --dx;      \
  else if (dx < -0.5) ++dx; \
  if (dy >= 0.5) --dy;      \
  else if (dy < -0.5) ++dy; \
  if (dz >= 0.5) --dz;      \
  else if (dz < -0.5) ++dz; \
  dx *= pbc_x, dy *= pbc_y, dz *= pbc_z; }

int natom, ilist, ipbc = 0, nthread = 4, \
    nstep, ibin, istep, ifresh, ithread, i, j;
double dr, sum_dr = 0, dbox = (double)rcut / nbin, \
       rlistshell2 = (rcut + rlistshell) * (rcut + rlistshell), \
       pbc_x = 1, pbc_y = 1, pbc_z = 1, dx, dy, dz, \
       rcut2 = rcut * rcut, rho;
int *list;
double *xyz, *xyz_old, *gr, *dr_max;

int main(int argc, char *argv[])
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
      rho = natom / (pbc_x * pbc_y * pbc_z); 
    } else if (strcmp(argv[i], "nthread") == 0) {
      nthread = atoi(argv[i + 1]);
    }
  }

  /* Check file status */
  if ((fp = fopen("md.out","rb")) == 0)
  {
    printf("Error: file 'md.out' does not exist!\n");
    exit(1);
  }

  /* Check thread setting */
  if (nthread > omp_get_max_threads())
  {
    printf("Error: require too many threads!\n");
    exit(1);
  }

  /* Set how many threads will be used (4 by default) */
  omp_set_num_threads(nthread);

  /* Initialize array */
  dr_max = (double *)malloc(sizeof(double) * nthread);
  gr = (double *)malloc(sizeof(double) * nbin * nthread);
  list = (int *)malloc(sizeof(double) * list_max * natom);
  xyz = (double *)malloc(sizeof(double) * 3 * natom);
  xyz_old = (double *)malloc(sizeof(double) * 3 * natom);
  memset(xyz, 0, sizeof(double) * 3 * natom);
  memset(gr, 0, sizeof(double) * nbin * nthread);

  for (istep = 0; istep < nstep; ++istep)
  {
    /* Read coordinate of atoms */
    memcpy(xyz_old, xyz, sizeof(double) * 3 * natom);
    memset(dr_max, 0, sizeof(double) * (nthread + 1));
    fread(xyz, sizeof(double), 3 * natom, fp);

    /* Calculate max displacement
     * To increase the speed, every thread has its own 'max displacement'
     * "The" max displacement will be determined after parallel part
     */
#pragma omp parallel for schedule(dynamic) \
        private(ithread, dx, dy, dz, dr) \
        firstprivate(pbc_x, pbc_y, pbc_z, natom)
    for (i = 0; i < natom; ++i)
    {
      ithread = omp_get_thread_num();
      dx = *(xyz + 3 * i) - *(xyz_old + 3 * i);
      dy = *(xyz + 3 * i + 1) - *(xyz_old + 3 * i + 1);
      dz = *(xyz + 3 * i + 2) - *(xyz_old + 3 * i + 2);
      if (ipbc) WrapPBC;
      dr = sqrt(dx * dx + dy * dy + dz * dz);
      *(dr_max + ithread) = (dr > *(dr_max + ithread)) ? \
                             dr : *(dr_max + ithread);
    }

    for (ithread = 1; ithread < nthread; ++ithread)
    {
      *dr_max = (*(dr_max + ithread) > *dr_max) ? \
                 *(dr_max + ithread) : *dr_max;
    }
    sum_dr += *dr_max;

    /* Neighbor-list algorithm is used in this code */
    if (sum_dr > 0.5 * rlistshell) 
    {
      /* Reset sum of max diplacement and neighbor-list */
      sum_dr = 0;
      memset(list, 0, sizeof(int) * list_max * natom);

#pragma omp parallel for schedule(dynamic) \
        private(j, dx, dy, dz, dr, ilist) \
        firstprivate(natom, ipbc, pbc_x, pbc_y, pbc_z, rlistshell2)
      for (i = 0; i < natom - 1; ++i)
      {
        ilist = 1;
        for (j = i + 1; j < natom; ++j)
        {
          dx = *(xyz + 3 * i) - *(xyz + 3 * j);
          dy = *(xyz + 3 * i + 1) - *(xyz + 3 * j + 1);
          dz = *(xyz + 3 * i + 2) - *(xyz + 3 * j + 2);
          if (ipbc) WrapPBC;
          dr = dx * dx + dy * dy + dz * dz;
          if (dr > rlistshell2) continue;
          *(list + i * list_max + ilist) = j;
          ++ilist;
        }
        /* Total number of list
         * To be aware: index 'ilist' starts from list[i][1]
         */
        *(list + i * list_max) = ilist;
      }
    }

    /* Calculate distance */
#pragma omp parallel for schedule(dynamic) \
        private(ithread, ilist, j, dx, dy, dz, dr, ibin) \
        firstprivate(dbox, ipbc, pbc_x, pbc_y, pbc_z)
    for (i = 0; i < natom - 1; ++i)
    {
      ithread = omp_get_thread_num();
      for (ilist = 1; ilist < *(list + i * list_max); ++ilist)
      {
        j = *(list + i * list_max + ilist);
        dx = *(xyz + 3 * i) - *(xyz + 3 * j);
        dy = *(xyz + 3 * i + 1) - *(xyz + 3 * j + 1);
        dz = *(xyz + 3 * i + 2) - *(xyz + 3 * j + 2);
        if (ipbc) WrapPBC;
        dr = dx * dx + dy * dy + dz * dz;
        if (dr > rcut2) continue;
        dr = sqrt(dr);
        ibin = (int)(dr / dbox);
        /* To increase the speed, every thread has its own g(r)
         * After parallel part, g(r) will be summed
         */
        ++ *(gr + ithread * nbin + ibin);
      }
    }
  }

  /* Sum g(r) */
#pragma omp parallel for schedule(dynamic) \
        private(ithread)
  for (ibin = 0; ibin < nbin; ++ibin)
  {
    for (ithread = 1; ithread < nthread; ++ithread)
    {
      *(gr + ibin) += *(gr + ithread * nbin + ibin);
    }
  }

  free(xyz);
  free(xyz_old);
  free(list);
  fclose(fp);

  /* Normalize and output */
  fp = fopen("rdf.dat","w");
  for (i = 0; i < nbin; ++i)
  {
    *(gr + i) /= 0.5 * natom * vshell * nstep * rho;
    fprintf(fp,"%-.3f   %-.6f\n", (i + 0.5) * dbox, *(gr + i));
  }
  fclose(fp);

  return 0;
}
