/* driver for retard.c */

#include <math.h>
#include <stdio.h>
#include "retard.h"

#define ndgl  3
#define ngrid 11
#define nrdens 1

char format99[] = "x=%f  y=%12.10f   nstep=%li\r\n";


void solout (long nr, double xold, double x, double* y, unsigned n, int* irtrn)
{
  static double xout;

  if (nr == 1)
  {
    printf (format99, x, y[0], nr-1);
    xout = x + 5.0;
  }
  else
    while (x >= xout)
    {
      printf (format99, x, y[0], nr-1);
      xout += 5.0;
    }

} /* solout */


double phi (unsigned i, double x)
{
  if (i == 1)
    return 0.1;
  else
    return 0.0;

} /* phi */


void fcn (unsigned n, double x, double* y, double* f)
{
  double y1l1, y1l10;

  y1l1 = ylag (1, x-1.0, phi);
  y1l10 = ylag (1, x-10.0, phi);
  f[0] = -y[0] * y1l1 + y1l10;
  f[1] = y[0] * y1l1 - y[1];
  f[2] = y[1] - y1l10;

} /* fcn */


int main(void)
{
  double y[ndgl];
  double grid[ngrid+1];
  unsigned icont[1];
  int i, res;
  int iout = 1;
  int itoler = 0;
  unsigned licont = 1;
  double x = 0.0;
  double xend = 40.0;
  double rtoler = 1.0E-5;
  double atoler = rtoler;

  y[0] = 5.0;
  y[1] = 0.1;
  y[2] = 1.0;

  icont[0] = 1;

  for (i = 0; i < ngrid-1; i++)
    grid[i] = i+1;
  grid[ngrid-1] = 20.0;

  res = retard (ndgl, fcn, x, y, xend, &rtoler, &atoler, itoler, solout, iout,
		stdout, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 100, nrdens,
		icont, licont, ngrid, grid);

  printf ("x=xend  y=%12.10f %12.10f %12.10f\r\n", y[0], y[1], y[2]);
  printf ("rtol=%12.10f   fcn=%li   step=%li   accpt=%li   rejct=%li\r\n",
	  rtoler, nfcnRead(), nstepRead(), naccptRead(), nrejctRead());

  return 0;

} /* main */



