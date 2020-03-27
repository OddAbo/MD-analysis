/* 2020-03-22
 * by dianmo
 * Self-intermediate scattering funtion Fs(q, t)
 * Sets several start points, calculates Fs until it reduces to zero
 * Run with command:
 * ./SISF  natom 1024  nstep 5000000  dt 0.001  nstart 1000 \
 *         pbc 16.20  qmax 2.327
 */

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
  fpos_t *fpos;
  FILE *fp;
  double *xyz, *xyznew, *Fs;
  double dt, dx, dy, dz, dr, qmax, Fs_tmp, pbc_x, pbc_y, pbc_z;
  int natom, nstep, nstart, ipbc, istart, iatom, nFs, nFs_max, iFs;

  /* Checks input arguments */
  for (int i = 1; i < argc; i += 2)
  {
    if (strcmp(argv[i], "natom") == 0) {
      natom = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "nstep") == 0) {
      nstep = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "nstart") == 0) {
      nstart = atoi(argv[i + 1]);
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
    printf("Error: file 'md.out' does not exist!\n");
    exit(1);
  }
  
  /* Initialize array */
  Fs = (double *)malloc(sizeof(double) * nFs);
  fpos = (fpos_t *)malloc(sizeof(fpos_t) * nstart);
  xyz = (double *)malloc(sizeof(double) * 3 * natom);
  xyznew = (double *)malloc(sizeof(double) * 3 * natom);
  memset(Fs, 0, sizeof(double) * nFs);

  /* Initialization of file pointer postion:
   * Sets start postion, initializes corresponding pointer position
   */
  for (istart = 0; istart < nstart; ++istart)
  {
    fgetpos(fp, (fpos + istart));
    fseek(fp, (long)sizeof(xyz) * nstep / nstart, SEEK_CUR);
  }

  /* Sets position and starts calculation */
  for (istart = 0; istart < nstart; ++istart)
  {
    fsetpos(fp, (fpos + istart));
    fread(xyz, sizeof(double), 3 * natom, fp);
    for (iFs = 1; iFs < nFs; ++iFs)
    {
      fread(xyznew, sizeof(double), 3 * natom, fp);

#pragma omp parallel for schedule(dynamic) \
        reduction(+:Fs_tmp) private(dx, dy, dz, dr)
      for (iatom = 0; iatom < natom; ++iatom)
      {
        Fs_tmp = 0;
        dx = *(xyz + 3 * iatom) - *(xyznew + 3 * iatom);
        dy = *(xyz + 3 * iatom + 1) - *(xyznew + 3 * iatom + 1);
        dz = *(xyz + 3 * iatom + 2) - *(xyznew + 3 * iatom + 2);
        if (ipbc)
        {
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
        Fs_tmp += cos(qmax * dr);
      }

      *(Fs + iFs) += Fs_tmp;
      if (Fs_tmp < 0.001) {
        nFs_max = (nFs_max > nFs) ? nFs_max : nFs;
        break;
      }
    }
  }

  fclose(fp);
  free(xyz);
  free(xyznew);

  /* Normalizes and outputs */
#pragma omp paralell for schedule(guided)
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
