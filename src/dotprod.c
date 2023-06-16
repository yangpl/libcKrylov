
#include "solver.h"

//compute dot product s=<x,y>
double dotprod(int n, double *a, double *b)
{
  int i;
  double s;
  
  for(i=0, s=0; i<n; i++) s += a[i]*b[i];
  return s;  
}
