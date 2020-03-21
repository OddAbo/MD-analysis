/* 2020-03-19
 * by dianmo
 * This is a simple example of allocating memory for pointer to array.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
  double* xyz;
  int natom, nstep;
  double dt;

  for (int i = 1; i < argc; i += 2)
  {
    if (strcmp(argv[i], "natom") == 0) {
      natom = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "nstep") == 0) {
      nstep = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "pbc") == 0) {
      pbc_box = atof(argv[i + 1]);
      ipbc = 1;
    }
  }
  return;
}
