/* 2020-03-22
 * by dianmo
 * Self-intermediate scattering funtion Fs(q, t)
 * Sets several start points, calculates Fs until it reduces to zero
 * Run with command:
 * ./SISF  natom 1024  nstep 500000  dt 0.001  nstart 1000  stepskip 5000 \
 *         pbc 16.20  max 2.327
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
  FILE *fp;
  double *xyz, *xyznew, *Fs;
  double dt, dx, dy, dz, dr, qmax, Fs_tmp, pbc_x, pbc_y, pbc_z;
  int natom, nstep, nstart, stepskip;
  int ipbc, istart, iatom, nFs, iFs;

  /* Checks input arguments */
  for (int i = 1; i < argc; i += 2)
  {
    if (strcmp(argv[i], "natom") == 0) {
      natom = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "nstep") == 0) {
      nstep = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "nstart") == 0) {
      nstart = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "stepskip") == 0) {
      stepskip = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "qmax") == 0) {
      qmax = atof(argv[i + 1]);
    } else if (strcmp(argv[i], "pbc") == 0) {
      pbc_x = atof(argv[i + 1]);
      pbc_y = pbc_x;
      pbc_z = pbc_x;
      ipbc = 1;
    } else if (strcmp(argv[i], "dt") == 0) {
      dt = atof(argv[i + 1]);
      nFs = 1000 / dt;
    }
  }

  if (nstart * nFs > nstep)
  {
    printf("Error: too many start points!\n");
    exit(1);
  }

  if ((fp = fopen("md.out","rb")) == NULL)
  {
    printf("Error: file does not exist!\n");
    exit(1);
  }

  xyz = (double *)malloc(3 * natom * sizeof(double));
  xyznew = (double *)malloc(3 * natom * sizeof(double));
  Fs = (double *)malloc(nFs * sizeof(double));

  for (istart = 0; istart < nstart; ++istart)
  {
    fread(xyz, sizeof(xyz), 1, fp);
    for (iFs = 0; iFs < nFs; ++iFs)
    {
      fread(xyznew, sizeof(xyz), 1, fp);
      Fs_tmp = 0;

#pragma omp paralell for schedule(guided) reduction(Fs_tmp) private(dx, dy, dz, dr)
      for (iatom = 0; iatom < natom; ++iatom)
      {
        dx = *(xyz + 3 * iatom) - *(xyznew + 3 * iatom);
        dy = *(xyz + 3 * iatom + 1) - *(xyznew + 3 * iatom + 1);
        dz = *(xyz + 3 * iatom + 2) - *(xyznew + 3 * iatom + 2);
        if (ipbc)
        {
          if (dx > 0.5) dx -= 1;
          else if (dx < -0.5) dx += 1;
          if (dy > 0.5) dy -= 1;
          else if (dy < -0.5) dy += 1;
          if (dz > 0.5) dz -= 1;
          else if (dz < -0.5) dz += 1;
          dx *= pbc_x;
          dy *= pbc_y;
          dz *= pbc_z;
        }
        dr = dx * dx + dy * dy + dz * dz;
        Fs_tmp += cos(qmax * dr);
      }

      *(Fs + iFs) += Fs_tmp;
      if (Fs_tmp < 0.001) break;
    }
  }

  fclose(fp);
  free(xyz);
  free(xyznew);

#pragma omp paralell for
  for (iFs = 0; iFs < nFs; ++iFs)
  {
    *(Fs + iFs) /= nstart;
  }

  fp = fopen("sisf.dat","w");
  for (iFs = 0; iFs < nFs; ++iFs)
  {
    fprintf(fp, "%-.3f    %-.6f\n", iFs * dt, *(Fs + iFs));
  }
  fclose(fp);
  free(Fs);

  return 0;
}
