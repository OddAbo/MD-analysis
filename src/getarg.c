/* 2020-03-19
 * by dianmo
 * This is a simple example of allocating memory for pointer to array
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
  double *xyz;
  int natom, nstep;
  double dt;
  
  for (int i = 1; i < argc; i += 2)
  {
    if (strcmp(argv[i], "natom") == 0) {
      natom = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "nstep") == 0) {
      nstep = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "dt") == 0) {
      dt = atof(argv[i + 1]);
    }
  }

  /* Allocate memory */
  xyz = (double *)malloc(3 * natom * sizeof(double));

  /* Then you can use 'xyz' as an array */
  xyz[0] = 3.1415926;
  xyz[1] = 2.7182818;
  xyz[2] = 42;

  free(xyz);

  return 0;
}
