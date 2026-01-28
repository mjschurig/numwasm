#include <math.h>
#include <stdio.h>
#include "dopri5.h"


#define  ndgl      4
#define  nrdens    2
#define  licont    nrdens

char format99[] = "x=%f  y=%12.10f %12.10f  nstep=%li\r\n";


void faren (unsigned n, double x, double *y, double *f)
{
  double amu, amup, r1, r2, sqr;

  amu = 0.012277471;
  amup = 1.0 - amu;
  f[0] = y[2];
  f[1] = y[3];
  sqr = y[0] + amu;
  r1 = sqr*sqr + y[1]*y[1];
  r1 = r1 * sqrt(r1);
  sqr = y[0] - amup;
  r2 = sqr*sqr + y[1]*y[1];
  r2 = r2 * sqrt(r2);
  f[2] = y[0] + 2.0 * y[3] - amup * (y[0]+amu) / r1 - amu * (y[0]-amup) / r2;
  f[3] = y[1] - 2.0 * y[2] - amup * y[1] / r1 - amu * y[1] / r2;

} /* faren */


void solout (long nr, double xold, double x, double* y, unsigned n, int* irtrn)
{
  static double xout;

  if (nr == 1)
  {
    printf (format99, x, y[0], y[1], nr-1);
    xout = x + 2.0;
  }
  else
    while (x >= xout)
    {
      printf (format99, xout, contd5(0,xout), contd5(1,xout), nr-1);
      xout += 2.0;
    }

} /* solout */


int main (void)
{
  double   y[ndgl];
  unsigned icont[licont], i;
  int      res, iout, itoler;
  double   x, xend, atoler, rtoler;

  iout = 2;
  x = 0.0;
  y[0] = 0.994;
  y[1] = 0.0;
  y[2] = 0.0;
  y[3] = -2.00158510637908252240537862224;
  xend = 17.0652165601579625588917206249;
  itoler = 0;
  rtoler = 1.0E-7;
  atoler = rtoler;
  icont[0] = 0;
  icont[1] = 1;

  res = dopri5 (ndgl, faren, x, y, xend, &rtoler, &atoler, itoler, solout, iout,
		stdout, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, ndgl, NULL, licont);

  printf ("x=xend  y=%12.10f %12.10f\r\n", y[0], y[1]);
  printf ("rtol=%12.10f   fcn=%li   step=%li   accpt=%li   rejct=%li\r\n",
	  rtoler, nfcnRead(), nstepRead(), naccptRead(), nrejctRead());

  return 0;

} /* main */



