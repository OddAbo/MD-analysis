/* 2020-03-14
 * by dianmo
 * This is a simple example of C language reads binary file formed by Fortran.
 * In Fortran, the file should be opened with access = 'stream'.
 */

#include <stdio.h>
#include <stdlib.h>

int main()
{
  /* Size of array is settled when define for convinience, but I'd recommend:
   * define pointers to arrays first,
   * then allocate memory for arrays depending on input arguments.
   */
  FILE* fp = fopen("md.out","rb");
  const int natom = 2048;
  double xyz[natom][3];

  fread(xyz, sizeof(double), 3 * natom, fp);
  fclose(fp);

  fp = fopen("xyz","w");
  for (int i = 0; i < natom; ++i)
  {
    fprintf(fp,"%-.6f  %-.6f  %-.6f\n", xyz[i][0], xyz[i][1], xyz[i][2]);
  }

  /* Always free memory when no longer used */
  free(xyz);

  return 0;
}
