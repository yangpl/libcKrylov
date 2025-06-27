#include <complex.h>

//compute dot product s=<x,y>=x^T*conj(y), this definiton is required according to Saad book
complex cdotprod(int n, complex *a, complex *b)
{
  int i;
  complex s;
  
  for(i=0, s=0; i<n; i++) s += a[i]*conj(b[i]);
  return s;  
}
