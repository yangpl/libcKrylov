#include <complex.h>

//compute dot product s=<x,y>=x^H*y
complex cdotprod(int n, complex *a, complex *b)
{
  int i;
  complex s;
  
  for(i=0, s=0; i<n; i++) s += conj(a[i])*b[i];
  return s;  
}
