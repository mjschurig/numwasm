#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <limits.h>
#include <memory.h>
#include "retard.h"


static long      nfcn, nstep, naccpt, nrejct;
static double    hout, xold, xout, x0, uround, hmax;
static int       last, ipos, idif, iact, iout;
static unsigned  nrds, *indir, mxst;
static int       irtrn;
static FILE*     fileout;
static double    *yy1, *k1, *k2, *k3, *k4, *k5, *k6, *ysti, *rcont;


long nfcnRead (void)
{
  return nfcn;

} /* nfcnRead */


long nstepRead (void)
{
  return nstep;

} /* stepRead */


long naccptRead (void)
{
  return naccpt;

} /* naccptRead */


long nrejctRead (void)
{
  return nrejct;

} /* nrejct */


double hRead (void)
{
  return hout;

} /* hRead */


double xRead (void)
{
  return xout;

} /* xRead */


static double sign (double a, double b)
{
  return (b > 0.0) ? fabs(a) : -fabs(a);

} /* sign */


static double min_d (double a, double b)
{
  return (a < b)?a:b;

} /* min_d */


static double max_d (double a, double b)
{
  return (a > b)?a:b;

} /* max_d */


static double hinit (unsigned n, FcnEqDiff fcn, double x, double* y,
	      double posneg, double* f0, double* f1, double* yy1, int iord,
	      double* atoler, double* rtoler, int itoler)
{
  double   dnf, dny, atoli, rtoli, sk, h, h1, der2, der12, sqr;
  unsigned i;

  dnf = 0.0;
  dny = 0.0;
  atoli = atoler[0];
  rtoli = rtoler[0];

  if (!itoler)
    for (i = 0; i < n; i++)
    {
      sk = atoli + rtoli * fabs(y[i]);
      sqr = f0[i] / sk;
      dnf += sqr*sqr;
      sqr = y[i] / sk;
      dny += sqr*sqr;
    }
  else
    for (i = 0; i < n; i++)
    {
      sk = atoler[i] + rtoler[i] * fabs(y[i]);
      sqr = f0[i] / sk;
      dnf += sqr*sqr;
      sqr = y[i] / sk;
      dny += sqr*sqr;
    }

  if ((dnf <= 1.0E-10) || (dny <= 1.0E-10))
    h = 1.0E-6;
  else
    h = sqrt (dny/dnf) * 0.01;

  h = min_d (h, hmax);
  h = sign (h, posneg);

  /* perform an explicit Euler step */
  for (i = 0; i < n; i++)
    yy1[i] = y[i] + h * f0[i];
  fcn (n, x+h, yy1, f1);

  /* estimate the second derivative of the solution */
  der2 = 0.0;
  if (!itoler)
    for (i = 0; i < n; i++)
    {
      sk = atoli + rtoli * fabs(y[i]);
      sqr = (f1[i] - f0[i]) / sk;
      der2 += sqr*sqr;
    }
  else
    for (i = 0; i < n; i++)
    {
      sk = atoler[i] + rtoler[i] * fabs(y[i]);
      sqr = (f1[i] - f0[i]) / sk;
      der2 += sqr*sqr;
    }
  der2 = sqrt (der2) / h;

  /* step size is computed such that h**iord * max_d(norm(f0),norm(der2)) = 0.01 */
  der12 = max_d (fabs(der2), sqrt(dnf));
  if (der12 <= 1.0E-15)
    h1 = max_d (1.0E-6, fabs(h)*1.0E-3);
  else
    h1 = pow (0.01/der12, 1.0/(double)iord);
  h = min_d (100.0 * h, min_d (h1, hmax));

  return sign (h, posneg);

} /* hinit */


/* core integrator */
static int retcor (unsigned n, FcnEqDiff fcn, double x, double* y,
		   double xend, double h, double* rtoler, double* atoler,
		   int itoler, SolTrait solout,
		   long nmax, int meth, long nstiff, double safe,
		   double beta, double fac1, double fac2, unsigned* icont,
		   unsigned ngrid, double* grid)
{
  double   facold, expo1, fac, facc1, facc2, fac11, posneg, xph;
  double   atoli, rtoli, hlamb, err, sk, hnew, yd0, ydiff, bspl;
  double   stnum, stden, sqr;
  int      iasti, iord, reject, nonsti, igrid;
  unsigned i, j, nrdl;
  double   c2, c3, c4, c5, e1, e3, e4, e5, e6, e7, d1, d3, d4, d5, d6, d7;
  double   a21, a31, a32, a41, a42, a43, a51, a52, a53, a54;
  double   a61, a62, a63, a64, a65, a71, a73, a74, a75, a76;

  /* initialisations */
  switch (meth)
  {
    case 1:

      c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0;
      a21=0.2, a31=3.0/40.0, a32=9.0/40.0;
      a41=44.0/45.0, a42=-56.0/15.0; a43=32.0/9.0;
      a51=19372.0/6561.0, a52=-25360.0/2187.0;
      a53=64448.0/6561.0, a54=-212.0/729.0;
      a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0;
      a64=49.0/176.0, a65=-5103.0/18656.0;
      a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0;
      a75=-2187.0/6784.0, a76=11.0/84.0;
      e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0;
      e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
      d1=-12715105075.0/11282082432.0, d3=87487479700.0/32700410799.0;
      d4=-10690763975.0/1880347072.0, d5=701980252875.0/199316789632.0;
      d6=-1453857185.0/822651844.0, d7=69997945.0/29380423.0;

      break;
  }

  facold = 1.0E-4;
  expo1 = 0.2 - beta * 0.75;
  facc1 = 1.0 / fac1;
  facc2 = 1.0 / fac2;
  posneg = sign (1.0, xend-x);

  /* initial preparations */
  iact = 1;
  ipos = 1;
  x0 = x;
  xend = grid[0];
  igrid = 0;
  uround *= 10.0;
  for (i = 0; i < mxst; i++)
    rcont[idif*i+1] = x;
  atoli = atoler[0];
  rtoli = rtoler[0];
  last  = 0;
  hlamb = 0.0;
  iasti = 0;
  hmax = fabs (hmax);
  irtrn = 2;
  fcn (n, x, y, k1);
  irtrn = 1;
  iord = 5;
  if (h == 0.0)
    h = hinit (n, fcn, x, y, posneg, k1, k2, k3, iord, atoler, rtoler, itoler);
  nfcn += 2;
  reject = 0;
  xold = x;
  if (iout)
  {
    hout = h;
    xout = x;
    solout (naccpt+1, xold, x, y, n, &irtrn);
    if (irtrn < 0)
    {
      if (fileout)
	fprintf (fileout, "Exit of retard at x = %.16e\r\n", x);
      return 2;
    }
  }

  /* basic integration step */
  while (1)
  {
    if (nstep > nmax)
    {
      if (fileout)
	fprintf (fileout, "Exit of retard at x = %.16e, more than nmax = %li are needed\r\n", x, nmax);
      xout = x;
      hout = h;
      return -2;
    }

    if (fabs(h) <= fabs(x) * uround)
    {
      if (fileout)
	fprintf (fileout, "Exit of retard at x = %.16e, step size too small h = %.16e\r\n", x, h);
      xout = x;
      hout = h;
      return -3;
    }

    if ((x + 1.01*h - xend) * posneg > 0.0)
    {
      h = xend - x;
      last = 1;
    }
    else if ((x+1.8*h-xend)*posneg > 0.0)
      h = (xend-x)*0.55;

    nstep++;

    /* the first 6 stages */
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * a21 * k1[i];
    fcn (n, x+c2*h, yy1, k2);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a31*k1[i] + a32*k2[i]);
    fcn (n, x+c3*h, yy1, k3);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a41*k1[i] + a42*k2[i] + a43*k3[i]);
    fcn (n, x+c4*h, yy1, k4);
    for (i = 0; i <n; i++)
      yy1[i] = y[i] + h * (a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
    fcn (n, x+c5*h, yy1, k5);
    for (i = 0; i < n; i++)
      ysti[i] = y[i] + h * (a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
    xph = x + h;
    fcn (n, xph, ysti, k6);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a71*k1[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
    irtrn = 1;
    fcn (n, xph, yy1, k2);

    /* prepare dense output */
    nrdl = 4 * nrds + iact;
    if (nrds == n)
      for (i = 0; i < n; i++)
      {
	rcont[nrdl+i+1] = h * (d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*k2[i]);
      }
    else
      for (j = 0; j < nrds; j++)
      {
	i = icont[j];
	rcont[nrdl+j+1] = h * (d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*k2[i]);
      }

    for (i = 0; i < n; i++)
      k4[i] = h * (e1*k1[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*k2[i]);
    nfcn += 6;

    /* error estimation */
    err = 0.0;
    if (!itoler)
      for (i = 0; i < n; i++)
      {
	sk = atoli + rtoli * max_d (fabs(y[i]), fabs(yy1[i]));
	sqr = k4[i] / sk;
	err += sqr*sqr;
      }
    else
      for (i = 0; i < n; i++)
      {
	sk = atoler[i] + rtoler[i] * max_d (fabs(y[i]), fabs(yy1[i]));
	sqr = k4[i] / sk;
	err += sqr*sqr;
      }
    err = sqrt (err / (double)n);

    /* computation of hnew */
    fac11 = pow (err, expo1);
    /* Lund-stabilization */
    fac = fac11 / pow(facold,beta);
    /* we require fac1 <= hnew/h <= fac2 */
    fac = max_d (facc2, min_d (facc1, fac/safe));
    hnew = h / fac;

    if (err <= 1.0)
    {
      /* step accepted */

      facold = max_d (err, 1.0E-4);
      naccpt++;

      /* stiffness detection */
      if (!(naccpt % nstiff) || (iasti > 0))
      {
	stnum = 0.0;
	stden = 0.0;
	for (i = 0; i < n; i++)
	{
	  sqr = k2[i] - k6[i];
	  stnum += sqr*sqr;
	  sqr = yy1[i] - ysti[i];
	  stden += sqr*sqr;
	}
	if (stden > 0.0)
	  hlamb = h * sqrt (stnum / stden);
	if (hlamb > 3.25)
	{
	  nonsti = 0;
	  iasti++;
	  if (iasti == 15)
	    if (fileout)
	      fprintf (fileout, "The problem seems to become stiff at x = %.16e\r\n", x);
	    else
	    {
	      xout = x;
	      hout = h;
	      return -4;
	    }
	}
	else
	{
	  nonsti++;
	  if (nonsti == 6)
	    iasti = 0;
	}
      }

      /* compute dense output */
      if (nrds == n)
	for (i = 0; i < n; i++)
	{
	  ydiff = yy1[i] - y[i];
	  bspl = h * k1[i] - ydiff;
	  rcont[iact+i+1] = y[i];
	  rcont[iact+nrds+i+1] = ydiff;
	  rcont[iact+2*nrds+i+1] = bspl;
	  rcont[iact+3*nrds+i+1] = -h * k2[i] + ydiff - bspl;
	}
      else
	for (j = 0; j < nrds; j++)
	{
	  i = icont[j];
	  ydiff = yy1[i] - y[i];
	  bspl = h * k1[i] - ydiff;
	  rcont[iact+j+1] = y[i];
	  rcont[iact+nrds+j+1] = ydiff;
	  rcont[iact+2*nrds+j+1] = bspl;
	  rcont[iact+3*nrds+j+1] = -h * k2[i] + ydiff - bspl;
	}

      rcont[iact] = x;
      iact += idif;
      rcont[iact-1] = h;
      if (iact+idif-1 > mxst*idif)
	iact = 1;

      memcpy (k1, k2, n * sizeof(double));
      memcpy (y, yy1, n * sizeof(double));
      xold = x;
      x = xph;

      if (irtrn == 3)
      {
	irtrn = 4;
	fcn (n, x, y, k1);
	nfcn++;
	irtrn = 1;
      }

      if (iout)
      {
	hout = h;
	xout = x;
	solout (naccpt+1, xold, x, y, n, &irtrn);
	if (irtrn < 0)
	{
	  if (fileout)
	    fprintf (fileout, "Exit of retard at x = %.16e\r\n", x);
	  return 2;
	}
      }

      /* normal exit */
      if (last)
	if (igrid == ngrid-1)
	{
	  hout=hnew;
	  xout = x;
	  return 1;
	}
	else
	{
	  igrid++;
	  last = 0;
	  xend = grid[igrid];
	  hnew = 0.9 * hnew;
	}

      if (fabs(hnew) > hmax)
	hnew = posneg * hmax;
      if (reject)
	hnew = posneg * min_d (fabs(hnew), fabs(h));

      reject = 0;
    }
    else
    {
      /* step rejected */
      if (irtrn < 0)
      {
	xout = x;
	hout = h;
	if (fileout)
	  fprintf (fileout, "Exit of retard at x = %.16e\r\n", x);
	return 2;
      }
      hnew = h / min_d (facc1, fac11/safe);
      reject = 1;
      if (naccpt >= 1)
	nrejct++;
      last = 0;
    }

    h = hnew;
    if (irtrn < 0)
    {
      xout = x;
      hout = h;
      return -5;
    }
  }

} /* retcor */


/* front-end */
int retard
 (unsigned n, FcnEqDiff fcn, double x, double* y, double xend, double* rtoler,
  double* atoler, int itoler, SolTrait solout, int iout_i, FILE* fileout_i, double uround_i,
  double safe, double fac1, double fac2, double beta, double hmax_i, double h,
  long nmax, int meth, long nstiff, unsigned maxbst, unsigned nrdens, unsigned* icont,
  unsigned licont, unsigned ngrid, double* grid)
{
  int       arret, idid;
  unsigned  i, lrcont;
  double    xuro;

  /* initialisations */
  nfcn = nstep = naccpt = nrejct = arret = 0;
  rcont = NULL;
  indir = NULL;
  fileout = fileout_i;

  /* n, the dimension of the system */
  if (n == UINT_MAX)
  {
    if (fileout)
      fprintf (fileout, "System too big, max. n = %u\r\n", UINT_MAX-1);
    arret = 1;
  }

  /* nmax, the maximal number of steps */
  if (!nmax)
    nmax = 100000;
  else if (nmax <= 0)
  {
    if (fileout)
      fprintf (fileout, "Wrong input, nmax = %li\r\n", nmax);
    arret = 1;
  }

  /* meth, coefficients of the method */
  if (!meth)
    meth = 1;
  else if ((meth <= 0) || (meth >= 2))
  {
    if (fileout)
      fprintf (fileout, "Curious input, meth = %i\r\n", meth);
    arret = 1;
  }

  /* nstiff, parameter for stiffness detection */
  if (!nstiff)
    nstiff = 1000;
  else if (nstiff < 0)
    nstiff = nmax + 10;

  /* iout, switch for calling solout */
  iout = iout_i;
  if ((iout < 0) || (iout > 1))
  {
    if (fileout)
      fprintf (fileout, "Wrong input, iout_i = %i\r\n", iout);
    arret = 1;
  }

  /* nrdens, number of dense output components */
  if (nrdens > n)
  {
    if (fileout)
      fprintf (fileout, "Curious input, nrdens = %u\r\n", nrdens);
    arret = 1;
  }
  else if (nrdens)
  {
    idif = 5*nrdens + 2;
    lrcont = (maxbst+1) * idif; /* +1 to keep the fortran indexes */
    mxst = maxbst;

    /* is there enough memory to allocate rcont&indir ? */
    rcont = (double*) malloc (lrcont*sizeof(double));
    if (nrdens < n)
      indir = (unsigned*) malloc (n*sizeof(unsigned));

    if (!rcont || (!indir && (nrdens < n)))
    {
      if (fileout)
	fprintf (fileout, "Not enough free memory for rcont&indir\r\n");
      arret = 1;
    }

    /* control of length of icont */
    if (nrdens == n)
    {
      if (icont && fileout)
	fprintf (fileout, "Warning : when nrdens = n there is no need allocating memory for icont\r\n");
      nrds = n;
    }
    else if (licont < nrdens)
    {
      if (fileout)
	fprintf (fileout, "Insufficient storage for icont, min. licont = %u\r\n", nrdens);
      arret = 1;
    }
    else
    {
      nrds = nrdens;
      for (i = 0; i < n; i++)
	indir[i] = UINT_MAX;
      for (i = 0; i < nrdens; i++)
	indir[icont[i]] = i;
    }
  }

  /* uround, smallest number satisfying 1.0+uround > 1.0 */
  uround = uround_i;
  if (uround == 0.0)
    uround = 2.3E-16;
  else if ((uround <= 1.0E-35) || (uround >= 1.0))
  {
    if (fileout)
      fprintf (fileout, "Which machine do you have ? Your uround_i was : %.16e\r\n", uround);
    arret = 1;
  }

  /* safety factor */
  if (safe == 0.0)
    safe = 0.9;
  else if ((safe >= 1.0) || (safe <= 1.0E-4))
  {
    if (fileout)
      fprintf (fileout, "Curious input for safety factor, safe = %.16e\r\n", safe);
    arret = 1;
  }

  /* fac1, fac2, parameters for step size selection */
  if (fac1 == 0.0)
    fac1 = 0.2;
  if (fac2 == 0.0)
    fac2 = 10.0;

  /* beta for step control stabilization */
  if (beta == 0.0)
    beta = 0.04;
  else if (beta < 0.0)
    beta = 0.0;
  else if (beta > 0.2)
  {
    if (fileout)
      fprintf (fileout, "Curious input for beta : beta = %.16e\r\n", beta);
    arret = 1;
  }

  /* maximal step size */
  hmax = hmax_i;
  if (hmax == 0.0)
    hmax = xend - x;

  /* grid with discontinuities */
  xuro = 100.0 * uround * fabs(xend);
  if (grid[ngrid-1] - xend >= xuro)
  {
    if (fileout)
      fprintf (fileout, "grid[ngrid-1] has to be <= xend\r\n");
    arret = 1;
  }
  if (fabs(grid[ngrid-1] - xend) >= xuro)
    ngrid++;
  grid[ngrid-1] = xend;

  /* is there enough free memory for the method ? */
  yy1 = (double*) malloc (n*sizeof(double));
  k1 = (double*) malloc (n*sizeof(double));
  k2 = (double*) malloc (n*sizeof(double));
  k3 = (double*) malloc (n*sizeof(double));
  k4 = (double*) malloc (n*sizeof(double));
  k5 = (double*) malloc (n*sizeof(double));
  k6 = (double*) malloc (n*sizeof(double));
  ysti = (double*) malloc (n*sizeof(double));

  if (!yy1 || !k1 || !k2 || !k3 || !k4 || !k5 || !k6 || !ysti)
  {
    if (fileout)
      fprintf (fileout, "Not enough free memory for the method\r\n");
    arret = 1;
  }

  /* when a failure has occured, we return -1 */
  if (arret)
  {
    if (ysti)
      free (ysti);
    if (k6)
      free (k6);
    if (k5)
      free (k5);
    if (k4)
      free (k4);
    if (k3)
      free (k3);
    if (k2)
      free (k2);
    if (k1)
      free (k1);
    if (yy1)
      free (yy1);
    if (indir)
      free (indir);
    if (rcont)
      free (rcont);

    return -1;
  }
  else
  {
    idid = retcor (n, fcn, x, y, xend, h, rtoler, atoler, itoler,
		   solout, nmax, meth, nstiff, safe, beta,
		   fac1, fac2, icont, ngrid, grid);
    free (ysti);
    free (k6);
    free (k5);    /* reverse order freeing too increase chances */
    free (k4);    /* of efficient dynamic memory managing       */
    free (k3);
    free (k2);
    free (k1);
    free (yy1);
    if (indir)
      free (indir);
    if (rcont)
      free (rcont);

    return idid;
  }

} /* retard */


/* dense output function */
double ylag (unsigned ii, double x, InitFuncPhi phi)
{
  unsigned i, j;
  int      inext;
  double   res, theta, theta1, compar, xright;

  /* initial phase */
  compar = uround * max_d (fabs(x), fabs(x0));
  if (x-x0 <= compar)
    if (irtrn <= 3)
    {
      res = phi (ii, x);
      if (irtrn == 2)
	hmax = min_d (hmax, x0-x);
      if (x0-x <= compar)
	irtrn = 3;
      return res;
    }
    else if (x0-x > compar)
      return phi (ii, x);

  /* compute the ii-th component place */
  i = UINT_MAX;

  if (!indir)
    i = ii;
  else
    i = indir[ii];

  if (i == UINT_MAX)
  {
    printf ("No dense output available for %uth component\r\n", ii);
    return 0.0;
  }

  /* compute the position of x */
  if (x-rcont[iact] < -compar)
  {
    if (fileout)
      fprintf (fileout, "Memory full, maxbst = %u\r\n", mxst);
    irtrn = -1;
    return 0.0;
  }

  inext = iact - idif;
  if (inext < 1)
    inext = (mxst-1) * idif + 1;
  xright = rcont[inext] + rcont[inext+idif-1];
  if (x-xright > uround*max_d(fabs(x),fabs(xright)))
  {
    if (fileout)
      fprintf (fileout, "Dont use advanced arguments\r\n");
    irtrn = -1;
    return 0.0;
  }

  while (x-rcont[ipos] < -compar)
  {
    ipos -= idif;
    if (ipos < 1)
      ipos = (mxst-1) * idif + 1;
  }

  inext = ipos + idif;
  if (inext > (mxst-1)*idif+1)
    inext = 1;

  while ((x > rcont[inext]) && (inext != iact))
  {
    ipos = inext;
    inext = ipos + idif;
    if (inext > (mxst-1)*idif+1)
      inext = 1;
  }

  /* compute the desired approximation */
  theta = (x - rcont[ipos]) / rcont[ipos+idif-1];
  theta1 = 1.0 - theta;
  i += ipos + 1;

  return rcont[i] + theta*(rcont[nrds+i] + theta1*(rcont[2*nrds+i] + theta*(rcont[3*nrds+i] + theta1*rcont[4*nrds+i])));

} /* ylag */

