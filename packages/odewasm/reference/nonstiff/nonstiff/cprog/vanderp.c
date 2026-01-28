#include <math.h>
#include <stdio.h>
#include "dop853.h"


#define  ndgl      2
#define  nrdens    2

char format99[] = "x=%f  y=%12.10f %12.10f  nstep=%li\r\n";


void fvpol (unsigned n, double x, double *y, double *f)
{
  const double eps = 1.0E-3;;

  f[0] = y[1];
  f[1] = ((1.0 - y[0]*y[0]) * y[1] - y[0]) / eps;

} /* fvpol */


void solout (long nr, double xold, double x, double* y, unsigned n, int* irtrn)
{
  static double xout; 

  if (nr == 1)
  { 
    printf ( "x=%f  y=%12.10f %12.10f  nstep=%li\r\n", x, y[0], y[1], nr-1);
    xout = x + 0.1;
  }
  else 
    while (x >= xout)
    {
      printf (format99, xout, contd8(0,xout), contd8(1,xout), nr-1);
      xout += 0.1;
    }
    
} /* solout */


int main (void)
{
  double   y[ndgl];
  int      res, iout, itoler;
  double   x, xend, atoler, rtoler;

  iout = 2;
  x = 0.0;
  y[0] = 2.0;
  y[1] = 0.0;
  xend = 2.0;
  itoler = 0;
  rtoler = 1.0E-6;
  atoler = rtoler;
  
  res = dop853 (ndgl, fvpol, x, y, xend, &rtoler, &atoler, itoler, solout, iout,
		stdout, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 1, ndgl, NULL, 0);

  printf ("x=xend  y=%12.10f %12.10f\r\n", y[0], y[1]);
  printf ("rtol=%12.10f   fcn=%li   step=%li   accpt=%li   rejct=%li\r\n",
	  rtoler, nfcnRead(), nstepRead(), naccptRead(), nrejctRead());

  return 0;

} /* main */


