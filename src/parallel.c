/* 2020-03-23
 * by dianmo
 * This is a simple test of omp paralell coding
 */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
  int i;
  double *xyz;

  xyz = (double *)malloc(50 * sizeof(double));
  
  for (i = 0; i < 50; ++i)
  {
    *(xyz + i) = i;
  }

#pragma omp parallel for schedule(dynamic)
  for (i = 0; i < 50; ++i)
  {
    *(xyz + i) += i;
  }
  
  for (i = 0; i < 50; ++i)
  {
    printf("%f\n",*(xyz + i));
  }
  return 0;
}
