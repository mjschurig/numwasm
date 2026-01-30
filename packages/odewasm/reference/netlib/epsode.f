      subroutine epsode (diffun, pederv, n, t0, h0, y0, tout, eps,
     1                   ierror, mf, index)
c this is the june 24, 1975  version of
c episode..  experimental package for integration of
c systems of ordinary differential equations,
c    dy/dt = f(y,t),  y = (y(1),y(2),...,y(n)) transpose,
c given the initial value of y.
c this code is for the ibm 370/195 at argonne national laboratory
c and is a modification of earlier versions by g.d.byrne
c and a.c.hindmarsh.
c                 references
c 1.  g. d. byrne and a. c. hindmarsh, a polyalgorithm for the
c       numerical solution of ordinary differential equations,
c       ucrl-75652, lawrence livermore laboratory, p. o. box 808,
c       livermore, ca 94550, april 1974. also in acm transactions
c       on mathematical software, 1 (1975), pp. 71-96.
c
c 2.  a. c. hindmarsh and g. d. byrne, episode.. an experimental
c       package for the integration of systems of ordinary
c       differential equations, ucid-30112, l.l.l., may, 1975.
c
c 3.  a. c. hindmarsh, gear.. ordinary differential equation
c       system solver, ucid-30001, rev. 3, l.l.l., december, 1974.
c
c-----------------------------------------------------------------------
c epsode is a driver subroutine for the episode package.
c epsode is to be called once for each output value of t.
c it then makes repeated calls to the core integrator
c subroutine, tstep.
c
c the input parameters are as follows.
c   diffun=  name of double precision subroutine diffun(n,t,y,ydot)
c            which the user must supply and declare external. see
c            note below.
c   pederv=  name of double precision subroutine pederv(n,t,y,pd,n0)
c            which the user must supply and declare external. see
c            note below.
c   n     =  the number of differential equations (used only on
c              first call, unless index = -1).  n must never be
c              increased during a given problem.
c   t0    =  the initial value of t, the independent variable
c              (used for input only on first call).
c   h0    =  the step size in t (used for input only on the
c              first call, unless index = 3 on input).  when
c              index = 3, h0 is the maximum absolute value of
c              the step size to be used.
c   y0    =  a vector of length n containing the initial values of
c              y (used for input only on first call).
c   tout  =  the value of t at which output is desired next.
c              integration will normally go beyond tout and
c              interpolate to t = tout.  (used only for input.)
c   eps   =  the relative error bound (used only on first call,
c              unless index = -1).  this bound is used as follows.
c              let r(i) denote the estimated relative local error
c              in y(i), i.e. the error relative to ymax(i), as
c              measured per step (of size h) or per ss units of t.
c              then eps is a bound on the root-mean-square norm
c              of the vector r, i.e.
c                      n
c              sqrt ( sum ( r(i)**2 )/n ) .lt. eps.
c                     i=1
c              the vector ymax is computed in epsode as described
c              under ierror below.
c              if error control per ss units of t is desired, set ss
c              to a positive number after statement 10 (where it is
c              now set to zero) and update it after statement 60.
c              see also the comments on ss and ymax below.
c  ierror =  the error flag with values and meanings as follow.
c         1    absolute error is controlled.  ymax(i) = 1.0.
c         2    error relative to abs(y) is controlled.  if y(i) = 0.0
c              a divide error will occur.  ymax(i) = abs(y(i)).
c         3    error relative to the largest value of abs(y(i)) seen
c              so far is controlled.  if the initial value of y(i) is
c              0.0, then ymax(i) is set to 1.0 initially and remains
c              at least 1.0.
c   mf    =  the method flag (used only on first call, unless
c              index = -1).  allowed values are 10, 11, 12, 13,
c              20, 21, 22, 23.  mf is an integer with two decimal
c              digits, meth and miter (mf = 10*meth + miter).  (mf
c              can be thought of as the ordered pair (meth,miter).)
c              meth is the basic method indicator.
c                meth = 1 indicates variable-step size, variable-
c                         order adams method, suitable for non-
c                         stiff problems.
c                meth = 2 indicates variable-step size, variable-
c                         order backward differentiation method,
c                         suitable for stiff problems.
c              miter indicates the method of iterative correction
c                (nonlinear system solution).
c                miter = 0 indicates functional iteration (no
c                          partial derivatives needed).
c                miter = 1 indicates a chord or semi-stationary
c                          newton method with closed form (exact)
c                          jacobian, which is computed in the
c                          user supplied subroutine
c                          pederv(n,t,y,pd,n0) described below.
c                miter = 2 indicates a chord  or semi-stationary
c                          newton method with an internally
c                          computed finite difference approximation
c                          to the jacobian.
c                miter = 3 indicates a chord or semi-stationary
c                          newton method with an internally
c                          computed diagonal matrix approximation
c                          to the jacobian, based on a directional
c                          derivative.
c   index =  integer used on input to indicate type of call,
c              with the following values and meanings..
c         1    this is the first call for this problem.
c         0    this is not the first call for this problem,
c              and integration is to continue.
c        -1    this is not the first call for the problem,
c              and the user has reset n, eps, and/or mf.
c         2    same as 0 except that tout is to be hit
c              exactly (no interpolation is done).
c              assumes tout .ge. the current t.
c         3    same as 0 except control returns to calling
c              program after one step.  tout is ignored.
c            since the normal output value of index is 0,
c            it need not be reset for normal continuation.
c
c after the initial call, if a normal return occurred and a normal
c continuation is desired, simply reset tout and call again.
c all other parameters will be ready for the next call.
c a change of parameters with index = -1 can be made after
c either a successful or an unsuccessful return.
c
c the output parameters are as follows.
c   t0    =  the output value of t.  if integration was successful,
c            t0 = tout.  otherwise, t0 is the last value of t
c            reached successfully.
c   h0    =  the step size h used last, whether successfully or not.
c   y0    =  the computed values of y at t = t0.
c   index =  integer used on output to indicate results,
c              with the following values and meanings..
c         0    integration was completed to tout or beyond.
c        -1    the integration was halted after failing to pass the
c              error test even after reducing h by a factor of
c              1.e10 from its initial value.
c        -2    after some initial success, the integration was
c              halted either by repeated error test failures or
c              by a test on eps.  possibly too much accuracy has
c              been requested, or a bad choice of mf was made.
c        -3    the integration was halted after failing to achieve
c              corrector convergence even after reducing h by a
c              factor of 1.e10 from its initial value.
c        -4    immediate halt because of illegal values of input
c              parameters.  see printed message.
c        -5    index was -1 on input, but the desired changes of
c              parameters were not implemented because tout
c              was not beyond t.  interpolation to t = tout was
c              performed as on a normal return.  to continue,
c              simply call again with index = -1 and a new tout.
c        -6    index was 2 on input, but tout was not beyond t.
c              no action was taken.
c
c in addition to epsode, the following subroutines are used by and
c provided in this package:
c   interp(tout,y,n0,y0)  interpolates to give output values at
c                         t = tout  by using data in the y array.
c   tstep(diffun,pederv,y,n0)  is the core integration subroutine, which
c                integrates over a single step and does associated error
c                control.
c   coset  sets coefficients for use in tstep.
c   adjust(y,n0)  adjusts the history array y on reduction of order.
c   pset(diffun,pederv,y,n0,con,miter,ier)  computes and processes the
c                             jacobian matrix, j = df/dy.
c   dec(n,n0,a,ip,ier)  performs the lu decomposition of a matrix.
c   sol(n,n0,a,b,ip)  solves a linear system a*x = b, after dec
c                       has been called for the matrix a.
c note:  pset, dec, and sol are called if and only if miter = 1
c        or miter = 2.
c
c the user must furnish the following double precision subroutines
c      and declare them external in his calling program:
c   diffun(n,t,y,ydot)  computes the function  ydot = f(y,t),
c                       the right hand side of the ordinary
c                       differential equation system, where y
c                       and ydot are vectors of length n.
c   pederv(n,t,y,pd,n0)  computes the n by n jacobian matrix of
c                        partial derivatives and stores it in pd as
c                        an n0 by n0 array.   pd(i,j) is to be set
c                        to the partial derivative of ydot(i) with
c                        respect to y(j).  pederv is called if and
c                        only if miter = 1. for other values of
c                        miter, pederv can be a dummy subroutine.
c
c caution:  at the present time the maximum number of differential
c           equations, which can be solved by episode, is 20.  to
c           change this number to a new value, say nmax, change
c           y(20,13) to y(nmax,13), ymax(20) to ymax(nmax),
c           error(20) to error(nmax), save1(20) to save1(nmax),
c           save2(20) to save2(nmax), pw(400) to pw(nmax*nmax),
c           and ipiv(20) to ipiv(nmax) in the common and dimension
c           statements below.  also change the argument in the
c           if...go to 440 statement (after the common statements)
c           from 20 to nmax.  no other changes need to be made to
c           any other subroutine in this package when the maximum
c           number of equations is changed.  elsewhere, the column
c           length of the y array is n0 instead of 20.  the row
c           length of y can be reduced from 13 to 6 if meth = 2.
c           the array ipiv is used if and only if miter = 1 or
c           miter = 2.  the size of the pw array can be reduced
c           to 1 if miter = 0 or to n if miter = 3.
c
c the common block epcom2 can be accessed externally by the user,
c if he desires to modify the error control behaviour of the package.
c the block contains ymax(20) used to control the error. see above
c and reference (2).
c
c the common block epcom9 can be accessed externally by the user,
c if he desires.  it contains the step size last used successfully
c (hused), the order last used successfully (nqused), the
c number of steps taken so far (nstep), the number of function
c evaluations (diffun calls) so far (nfe), and the number of
c jacobian evaluations so far (nje).
c
c in a data statement below, lout is set to the logical unit number
c for the output of messages during integration.  currently, lout
c = 6.
c-----------------------------------------------------------------------
c
      external diffun, pederv
      integer ierror, index, mf, n
      integer ipiv, jstart, kflag, mfc, nc, nfe, nje,
     1        nqused, nsq, nstep
      integer i, kgo, nhcut, n0
      integer lout
      double precision eps, h0, tout, t0, y0
      double precision epsc, epsj, error, hmax, h, hmin, hused,
     1                 pw, save1, save2, ss, t, uround, ymax
      double precision ayi, d, epsj, h0, t0p, y
      double precision hcut
      double precision four, hundrd, one, ten, zero
      dimension y(20,13)
      dimension y0(n)
c
      common /epcom1/t,h,hmin,hmax,epsc,ss,uround,nc,mfc,kflag,jstart
      common /epcom2/ ymax(20)
      common /epcom3/ error(20)
      common /epcom4/ save1(20)
      common /epcom5/ save2(20)
      common /epcom6/ pw(400)
      common /epcom7/ ipiv(20)
      common /epcom8/ epsj,nsq
      common /epcom9/ hused,nqused,nstep,nfe,nje
c
      data lout /6/
      data hcut /0.1d0/
      data four /4.0d0/, hundrd /1.0d2/, one  /1.0d0/,
     1     ten  /1.0d1/, zero   /0.0d0/
      if (index .eq. 0) go to 20
      if (index .eq. 2) go to 25
      if (index .eq. -1) go to 30
      if (index .eq. 3) go to 40
      if (index .ne. 1) go to 430
      if (eps .le. zero) go to 400
      if (n .le. 0) go to 410
      if (n .gt. 20) go to 440
      if ((t0-tout)*h0 .ge. zero) go to 420
c-----------------------------------------------------------------------
c if initial values for ymax other than those below are desired,
c they should be set here.  all ymax(i) must be positive.  if
c values for hmin or hmax, the bounds on the absolute value of h,
c other than those below, are desired, they also should be set here.
c if error per ss units of t is to be controlled, ss should be set
c to a positive value below.  error per unit step is controlled
c when ss = 1.  the default value for ss is 0 and yields control
c of error per step.
c-----------------------------------------------------------------------
c set uround, the machine roundoff constant, here.
c use statement below for short precision on ibm 360 or 370.
c     uround = 9.53674e-7
c use statement below for single precision on cdc 7600 or 6600.
c     uround = 7.105427406e-15
c use statement below for long precision on ibm 360 or 370.
c     uround = 2.220446047d-16
      uround = d1mach(4)
      do 10 i = 1,n
        go to (5, 6, 7), ierror
c ierror   =   1, 2, 3  ------------------------------------------------
  5     ymax(i) = one
        go to 10
  6     ymax(i) = dabs(y0(i))
        go to 10
  7     ymax(i) = dabs(y0(i))
        if (ymax(i) .eq. zero) ymax(i) = one
 10     y(i,1) = y0(i)
      nc = n
      t = t0
      h = h0
      if ((t+h) .eq. t) write(lout,15) t
 15   format(/47h---  message from subroutine epsode in episode,,
     1       24h the o.d.e. solver.  ---/22h warning.. t + h = t =,
     2       e18.8,18h in the next step./)
      hmin = dabs(h0)
      hmax = dabs(t0 - tout)*ten
      epsc = eps
      mfc = mf
      jstart = 0
      ss = zero
      n0 = n
      nsq = n0*n0
      epsj = dsqrt(uround)
      nhcut = 0
      go to 50
c t0p is the previous output value of t0 for use in hmax. --------------
  20  hmax = dabs(tout - t0p)*ten
      go to 80
  25  hmax = dabs(tout - t0p)*ten
c
       if ((t-tout)*h .ge. zero) go to 460
       go to 85
c
 30   if ((t-tout)*h .ge. zero) go to 450
      if (mf .ne. mfc) jstart = -1
      nc = n
      epsc = eps
      mfc = mf
      go to 45
c
 40   hmax = h0
c
 45   if ((t+h) .eq. t) write(lout,15) t
c
 50   call tstep (diffun, pederv, y, n0)
c
      kgo = 1 - kflag
      go to (60, 100, 200, 300), kgo
c kflag  =   0,  -1,  -2,  -3  -----------------------------------------
c
 60   continue
c-----------------------------------------------------------------------
c normal return from tstep.
c
c the weights ymax(i) are updated.  if different values are desired,
c they should be set here.  if ss is to be updated for control of
c error per ss units of t, it should also be done here.  a test is
c made to determine if eps is too small for machine precision.
c
c any other tests or calculations that are required after each step
c should be inserted here.
c
c if index = 3, y0 is set to the current y values on return.
c if index = 2, h is controlled to hit tout (within roundoff
c error), and then the current y values are put in y0 on
c return.  for any other value of index, control returns to
c the integrator unless tout has been reached.  then
c interpolated values of y are computed and stored in y0 on
c return.
c if interpolation is not desired, the call to interp should
c be deleted and control transferred to statement 500 instead
c of 520.
c-----------------------------------------------------------------------
      d = zero
      do 70 i = 1,n
        ayi = dabs(y(i,1))
        go to (70, 66, 67), ierror
c ierror  =     1,  2,  3  ---------------------------------------------
 66     ymax(i) = ayi
        go to 70
  67    ymax(i) = dmax1(ymax(i), ayi)
  70    d = d + (ayi/ymax(i))**2
      d = d*(uround/eps)**2
      if (d .gt. dfloat(n)) go to 250
      if (index .eq. 3) go to 500
      if (index .eq. 2) go to 85
 80   if ((t-tout)*h .lt. zero) go to 45
      call interp (tout, y, n0, y0)
      t0 = tout
      go to 520
 85   if (((t+h)-tout)*h .le. zero) go to 45
      if (dabs(t-tout) .le. hundrd*uround*hmax) go to 500
      if ((t-tout)*h .ge. zero) go to 500
      h = (tout - t)*(one - four*uround)
      jstart = -1
      go to 45
c-----------------------------------------------------------------------
c on an error return from tstep, an immediate return occurs if
c kflag = -2, and recovery attempts are made otherwise.
c to recover, h and hmin are reduced by a factor of .1 up to 10
c times before giving up.
c-----------------------------------------------------------------------
 100  write (lout,101)
 101  format(/47h---  message from subroutine epsode in episode,,
     1        24h the o.d.e. solver.  ---/)
      write(lout,105) t,hmin
 105  format(//35h kflag = -1 from integrator at t = ,e18.8/
     1       40h  error test failed with abs(h) = hmin =,e18.8/)
 110  if (nhcut .eq. 10) go to 150
      nhcut = nhcut + 1
      hmin = hcut*hmin
      h = hcut*h
      write (lout,115) h
 115  format(24h  h has been reduced to ,e18.8,
     1       26h  and step will be retried//)
      jstart = -1
      go to 45
c
 150  write (lout,155)
 155  format(//44h problem appears unsolvable with given input//)
      go to 500
c
 200  write (lout,101)
      write (lout,205) t,h,eps
 205  format(//16h kflag = -2  t =,e18.8,4h h =,e18.8,6h eps =,e18.8/
     1       50h  the requested error is too small for integrator.//)
      go to 500
c
 250  write (lout,101)
      write (lout,255) t,eps
 255  format(//47h integration halted by subroutine epsode at t =,
     1       e18.8/43h eps is too small for machine precision and/
     2       29h problem being solved.  eps =,e18.8//)
      kflag = -2
      go to 500
c
 300  write (lout,101)
      write (lout,305) t
 305  format(//34h kflag = -3 from integrator at t =,e18.8/
     1       45h  corrector convergence could not be achieved/)
      go to 110
c
 400  write (lout,101)
      write (lout,405) eps
 405  format(//35h illegal input.. eps .le. 0. eps = ,e18.8//)
      index = -4
      return
c
 410  write (lout,101)
      write (lout,415) n
 415  format(//31h illegal input.. n .le. 0. n = ,i8//)
      index = -4
      return
c
 420  write (lout,101)
      write (lout,425) t0,tout,h0
 425  format(//39h illegal input.. (t0 - tout)*h0 .ge. 0./
     1       5h t0 =,e18.8,7h tout =,e18.8,5h h0 =,e18.8//)
      index = -4
      return
c
 430  write (lout,101)
      write (lout,435) index
 435  format(//24h illegal input.. index =,i8//)
      index = -4
      return
c
 440  write (lout,101)
      write (lout,445) n
 445  format (//39h illegal input.  the number of ordinary/
     1        43h differential equations being solved is n =, i6/
     2        43h storage allocation in subroutine epsode is/
     3        47h too small.  see writeups#epsode public......../)
      index = -4
      return
c
 450  write (lout,101)
      write (lout,455) t,tout,h
 455  format(//46h index = -1 on input with (t - tout)*h .ge. 0./
     1       44h interpolation was done as on normal return./
     2       41h desired parameter changes were not made./
     3       4h t =,e18.8,7h tout =,e18.8,4h h =,e18.8//)
      call interp (tout, y, n0, y0)
      t0 = tout
      index = -5
      return
c
 460  write (lout,101)
      write (lout,465) t,tout,h
 465  format(//45h index = 2 on input with (t - tout)*h .ge. 0./
     1       4h t =,e18.8,7h tout =,e18.8,4h h =,e18.8//)
      index = -6
      return
c
 500  t0 = t
      do 510 i = 1,n
 510    y0(i) = y(i,1)
 520  index = kflag
      t0p = t0
      h0 = hused
      if (kflag .ne. 0) h0 = h
      return
c-----------------------  end of subroutine epsode ---------------------
      end
      subroutine interp (tout, y, n0, y0)
c-----------------------------------------------------------------------
c subroutine interp computes interpolated values of the dependent
c variable y and stores them in y0.  the interpolation is to the
c point t = tout and uses the nordsieck history array y as follows..
c                             nq
c                  y0(i)  =  sum  y(i,j+1)*s**j ,
c                            j=0
c where s = -(t-tout)/h.
c-----------------------------------------------------------------------
c caution:  not all members of epcom1 are used in this subroutine.
c-----------------------------------------------------------------------
      integer n0
      integer jstart, kflag, mf, n
      integer i, j, l
      double precision tout, y, y0
      double precision eps, h, hmax, hmin, ss, t, uround
      double precision s, s1
      double precision one
      dimension y0(n0),y(n0,13)
c
      common /epcom1/ t,h,hmin,hmax,eps,ss,uround,n,mf,kflag,jstart
      data one /1.0d0/
      do 10 i = 1,n
 10     y0(i) = y(i,1)
      l = jstart + 1
      s = (tout - t)/h
      s1 = one
      do 30 j = 2,l
        s1 = s1*s
        do 20 i = 1,n
 20       y0(i) = y0(i) + s1*y(i,j)
 30     continue
      return
c-----------------------  end of subroutine interp  --------------------
      end
      subroutine tstep (diffun, pederv, y, n0)
c-----------------------------------------------------------------------
c tstep performs one step of the integration of an initial value
c problem for a system of ordinary differential equations.
c communication with tstep is via the following variables..
c
c   y       an n0 by lmax array containing the dependent variables
c             and their scaled derivatives.  lmax is currently 6 for
c             the variable step backward differentiation formulas,
c             and 13 for the variable step adams formulas.
c             (lmax -1) = maxder, the maximum order used.
c             see subroutine coset.  y(i,j+1) contains the
c             j-th derivative of y(i), scaled by h**j/factorial(j)
c             for j = 0,1,...,nq, where nq is the current order.
c   n0      a constant integer .ge. n, used for dimensioning
c             purposes.
c   t       the independent variable, updated on each step taken.
c   h       the step size to be attempted on the next step.
c             h is altered by the error control algorithm during
c             the solution of the problem. h can be either positive
c             or negative, but its sign must remain constant
c             throughout the problem run.
c   hmin,   the minimum and maximum absolute values of the step
c    hmax     size to be used for the step. these may be changed at
c             any time, but the change will not take effect until the
c             next change in h is made.
c   eps     the relative error bound.  see description in
c             subroutine epsode.
c   ss      the size of the time interval to be used for error
c             control.  a default value of 0 is used to produce
c             control of error per step.  see subroutine epsode.
c   uround  the unit of roundoff for the computer being used.
c   n       the number of first order ordinary differential
c             equations being solved.
c   mf      the method flag.  see description in subroutine epsode.
c   kflag   a completion code with the following meanings..
c                     0  the step was succesful.
c                    -1  the requested error could not be achieved
c                          with abs(h) = hmin.
c                    -2  the requested error is smaller than can
c                          be handled for this problem.
c                    -3  corrector convergence could not be
c                          achieved for abs(h) = hmin.
c             on a return with kflag negative, the values of t and
c             the y array are as of the beginning of the last
c             step and h is the last step size attempted.
c   jstart  an integer used on input and output.
c             on input, it has the following values and meanings..
c                     0  perform the first step.
c                 .gt.0  take a new step continuing from the last.
c                 .lt.0  take the next step with a new value of
c                          h and/or mf.
c             on exit, jstart is set to nq, the current order of the
c             method.
c   ymax      an array of n elements with which the estimated local
c               errors in y are compared.
c   error     an array of n elements.  error(i)/tq(2) is the
c               estimated local error in y(i) per ss units of
c               t or per step (of size h).
c   save1,    two arrays for working storage,
c     save2     each of length n.
c   pw        a block of locations used for the partial derivatives
c               of f with respect to y, if miter is not o.  see
c               description in subroutine epsode.
c   ipiv      an integer array of length n, which is used for pivot
c               information for the linear algebraic system in the
c               correction process, when miter = 1 or 2.
c
c the common block epcm10, declared below, is primarily intended
c for internal use, but it can be accessed externally.
c-----------------------------------------------------------------------
      external pederv
      integer n0
      integer ipiv, jstart, kflag, l, lmax, meth, mf, n, nfe, nje,
     1        nq, nqindx, nqused, nstep
      integer i, iback, ier, iredo, j, j1, j2, m, mfold, mio,
     1        miter, miter1, newj, nstepj
      integer istepj, kfc, kfh, maxcor
      double precision y
      double precision el, eps, error, h, hmax, hmin, hused, pw,
     1                 save1, save2, ss, t, tau, tq, uround, ymax
      double precision bnd, cnquot, con, conp, crate, d, drc,
     1                 d1, e, edn, eta, etamax, etamin, etaq, etaqm1,
     2                 etaqp1, eup, flotl, flotn, hold, hrl1, phrl1,
     3                 prl1, r, rc, rl1, r0, r1, told
      double precision addon, bias1, bias2, bias3, crdown, delrc,
     1                 etacf, etamin, etamxf, etamx1, etamx2,
     2                 etamx3, onepsm, short, thresh
      double precision one, pt5, zero
      dimension y(n0,13)
c
      common /epcom1/ t,h,hmin,hmax,eps,ss,uround,n,mf,kflag,jstart
      common /epcom2/ ymax(1)
      common /epcom3/ error(1)
      common /epcom4/ save1(1)
      common /epcom5/ save2(1)
      common /epcom6/ pw(1)
      common /epcom7/ ipiv(1)
      common /epcom9/ hused,nqused,nstep,nfe,nje
      common /epcm10/ tau(13),el(13),tq(5),lmax,meth,nq,l,nqindx
c
      data istepj /20/, kfc /-3/, kfh /-7/, maxcor /3/
      data addon  /1.0d-6/,    bias1  /2.5d1/, bias2  /2.5d1/,
     1     bias3  /1.0d2/,     crdown /0.1d0/, delrc  /0.3d0/,
     2     etacf  /0.25d0/,    etamin /0.1d0/, etamxf /0.2d0/,
     3     etamx1 /1.0d4/,     etamx2 /1.0d1/, etamx3 /1.5d0/,
     4     onepsm /1.00001d0/, short  /0.1d0/, thresh /1.3d0/
      data one /1.0d0/, pt5 /0.5d0/, zero /0.0d0/
      kflag = 0
      told = t
      flotn = dfloat(n)
      if (jstart .gt. 0) go to 200
      if (jstart .ne. 0) go to 150
c-----------------------------------------------------------------------
c on the first call, the order is set to 1 and the initial
c derivatives are calculated.  etamax is the maximum ratio by
c which h can be increased in a single step.  it is 1.e04 for the
c first step to compensate for the small initial h, then 10 for
c the next 10 steps, and then 1.5 thereafter.  if a failure
c occurs (in corrector convergence or error test), etamax is set at 1
c for the next increase.  etamin = .1 is the minimum ratio by which
c h can be reduced on any retry of a step.
c-----------------------------------------------------------------------
      call diffun (n, t, y, save1)
      do 110 i = 1,n
 110    y(i,2) = h*save1(i)
      meth = mf/10
      miter = mf - 10*meth
      miter1 = miter + 1
      mfold = mf
      nq = 1
      l = 2
      tau(1) = h
      prl1 = one
      rc = zero
      etamax = etamx1
      nqindx = 2
      nstep = 0
      nstepj = 0
      nfe = 1
      nje = 0
      go to 200
c-----------------------------------------------------------------------
c if the user has changed h, then y must be rescaled.  if the
c user has changed miter, then newj is set to miter to force
c the partial derivativees to be updated, if they are being used.
c-----------------------------------------------------------------------
 150  if (mf .eq. mfold) go to 170
      mio = miter
      meth = mf/10
      miter = mf - 10*meth
      mfold = mf
      if (miter .eq. mio) go to 170
      newj = miter
      miter1 = miter + 1
 170  if (h .eq. hold) go to 200
      eta = h/hold
      h = hold
      iredo = 3
      go to 185
  180 eta = dmax1(eta,hmin/dabs(h),etamin)
  185 eta = dmin1(eta,hmax/dabs(h),etamax)
      r1 = one
      do 190 j = 2,l
        r1 = r1*eta
        do 190 i = 1,n
 190      y(i,j) = y(i,j)*r1
      h = h*eta
      rc = rc*eta
      if (iredo .eq. 0) go to 690
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the y array by the pascal triangle matrix.  then
c coset is called to obtain el, the vector of coefficients of
c length nq + 1.  rc is the ratio of new to old values of the
c coefficient h/el(2).  when rc differs from 1 by more than
c delrc, newj is set to miter to force the partial derivatives
c to be updated, if used.  delrc is 0.3.  in any case, the partial
c derivatives are updated at least every 20-th step.
c-----------------------------------------------------------------------
 200  t = t + h
      do 210 j1 = 1,nq
        do 210 j2 = j1,nq
          j = (nq + j1) - j2
          do 210 i = 1,n
 210        y(i,j) = y(i,j) + y(i,j+1)
      call coset
      bnd = flotn*(tq(4)*eps)**2
      rl1 = one/el(2)
      rc = rc*(rl1/prl1)
      prl1 = rl1
      if (nstep .ge. nstepj+istepj) newj = miter
      drc = dabs(rc-one)
      if (drc .le. delrc) go to 215
      newj = miter
      crate = one
      rc = one
      go to 220
  215 if ((miter .ne. 0) .and. (drc .ne. zero)) crate = one
c-----------------------------------------------------------------------
c up to 3 corrector iterations are taken.  a convergence test is made
c on the root mean square norm of each correction, using bnd, which
c is dependent on eps.  the sum of the corrections is accumulated in
c the vector error.  the y array is not altered in the corrector
c loop.  the updated y vector is stored temporarily in save1.
c-----------------------------------------------------------------------
 220  do 230 i = 1,n
 230    error(i) = zero
      m = 0
      call diffun (n, t, y, save2)
      nfe = nfe + 1
      if (newj .le. 0) go to 290
c-----------------------------------------------------------------------
c if indicated, the matrix p = i - h*rl1*j is reevaluated before
c starting the corrector iteration.  newj is set to 0 as an
c indicator that this has been done.  if miter = 1 or 2, p is
c computed and processed in pset.  if miter = 3, the matrix  is
c p = i - h*rl1*d, where d is a diagonal matrix.  rl1 is 1/el(2).
c-----------------------------------------------------------------------
      newj = 0
      rc = one
      nje = nje + 1
      nstepj = nstep
      go to (250, 240, 260), miter
 240  nfe = nfe + n
 250  con = -h*rl1
      call pset(diffun, pederv, y, n0, con, miter, ier)
      if (ier .ne. 0) go to 420
      go to 350
 260  r = rl1*short
      do 270 i = 1,n
 270    pw(i) = y(i,1) + r*(h*save2(i) - y(i,2))
      call diffun(n, t, pw, save1)
      nfe = nfe + 1
      hrl1 = h*rl1
      do 280 i = 1,n
        r0 = h*save2(i) - y(i,2)
        pw(i) = one
        d = short*r0 - h*(save1(i) - save2(i))
        save1(i) = zero
        if (dabs(r0) .lt. uround*ymax(i)) go to 280
        if (dabs(d) .eq. zero) go to 420
        pw(i) = short*r0/d
        save1(i) = pw(i)*rl1*r0
 280    continue
      go to 370
 290  go to (295, 350, 350, 310), miter1
c-----------------------------------------------------------------------
c in the case of functional iteration, y is updated directly from
c the result of the last diffun call.
c-----------------------------------------------------------------------
 295  d = zero
      do 300 i = 1,n
        r = rl1*(h*save2(i) - y(i,2))
        d = d + ((r - error(i))/ymax(i))**2
        save1(i) = y(i,1) + r
 300    error(i) = r
      go to 400
c-----------------------------------------------------------------------
c in the case of a chord method, the residual -g(y sub n(m))
c is computed and the linear system with that as right-hand side
c and p as coefficient matrix is solved, using the lu decomposition
c of p if miter = 1 or 2.  if miter = 3 the scalar h*rl1 is updated.
c-----------------------------------------------------------------------
 310  phrl1 = hrl1
      hrl1 = h*rl1
      if (hrl1 .eq. phrl1) go to 330
      r = hrl1/phrl1
      do 320 i = 1,n
        d = one - r*(one - one/pw(i))
        if (dabs(d) .eq. zero) go to 440
 320    pw(i) = one/d
 330  do 340 i = 1,n
 340    save1(i) = pw(i)*(rl1*h*save2(i) - (rl1*y(i,2) + error(i)))
      go to 370
 350  do 360 i = 1,n
 360    save1(i) = rl1*h*save2(i) - (rl1*y(i,2) + error(i))
      call sol (n, n0, pw, save1, ipiv)
 370  d = zero
      do 380 i = 1,n
        error(i) = error(i) + save1(i)
        d = d + (save1(i)/ymax(i))**2
 380    save1(i) = y(i,1) + error(i)
c-----------------------------------------------------------------------
c test for convergence.  if m .gt. 0, an estimate of the square of
c the convergence rate constant is stored in crate, and this is used
c in the test.
c-----------------------------------------------------------------------
  400 if (m .ne. 0) crate = dmax1(crdown*crate,d/d1)
      if (d*dmin1(one,crate) .le. bnd) go to 450
      d1 = d
      m = m + 1
      if (m .eq. maxcor) go to 410
      call diffun (n, t, save1, save2)
      go to (295, 350, 350, 310), miter1
c-----------------------------------------------------------------------
c the corrector iteration failed to converge in 3 tries. if partial
c derivatives are involved but are not up to date, they are
c reevaluated for the next try.  otherwise the y array is restored
c to its values before prediction, and h is reduced,
c if possible.  if not, a no-convergence exit is taken.
c-----------------------------------------------------------------------
 410  nfe = nfe + maxcor - 1
      if (newj .eq. -1) go to 440
 420  t = told
      etamax = one
      do 430 j1 = 1,nq
        do 430 j2 = j1,nq
          j = (nq + j1) - j2
          do 430 i = 1,n
 430        y(i,j) = y(i,j) - y(i,j+1)
      if (dabs(h) .le. hmin*onepsm) go to 680
      eta = etacf
      iredo = 1
      go to 180
 440  newj = miter
      go to 220
c-----------------------------------------------------------------------
c the corrector has converged.  newj is set to -1 if partial
c derivatives were used, to signal that they may need updating on
c subsequent steps.  the error test is made and control passes to
c statement 500 if it fails.
c-----------------------------------------------------------------------
 450  if (miter .ne. 0) newj = -1
      nfe = nfe + m
      d = zero
      do 460 i = 1,n
 460    d = d + (error(i)/ymax(i))**2
      e = flotn*(tq(2)*eps)**2
      if (d .gt. e) go to 500
c-----------------------------------------------------------------------
c after a successful step, the y array, tau, nstep, and nqindx are
c updated, and a new value of h at order nq is computed.
c the vector tau contains the nq + 1 most recent values of h.
c a change in nq up or down by 1 is considered if nqindx = 0.
c if nqindx = 1 and nq .lt. maxder, then error is saved
c for use in a possible order increase on the next step.
c a change in h or nq is made only of the increase in h
c is by a factor of at least 1.3.
c if not, nqindx is set to 2 to prevent testing for that many
c steps.  if nq is changed, nqindx is set to nq + 1 (new value).
c-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nstep = nstep + 1
      hused = h
      nqused = nq
      do 470 iback = 1,nq
        i = l - iback
 470    tau(i+1) = tau(i)
      tau(1) = h
      do 480 j = 1,l
        do 480 i = 1,n
 480      y(i,j) = y(i,j) + error(i)*el(j)
      nqindx = nqindx - 1
      if ((l .eq. lmax) .or. (nqindx .ne. 1)) go to 495
      do 490 i = 1,n
 490    y(i,lmax) = error(i)
      conp = tq(5)
 495  if (etamax .ne. one) go to 520
      if (nqindx .lt. 2) nqindx = 2
      go to 690
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c t and the y array are restored to their previous values.  a new
c h for a retry of the step is computed.  the order is kept fixed.
c-----------------------------------------------------------------------
 500  kflag = kflag - 1
      t = told
      do 510 j1 = 1,nq
        do 510 j2 = j1,nq
          j = (nq + j1) - j2
          do 510 i = 1,n
 510        y(i,j) = y(i,j) - y(i,j+1)
      newj = miter
      etamax = one
      if (dabs(h) .le. hmin*onepsm) go to 660
      if (kflag .le. kfc) go to 630
      iredo = 2
c compute ratio of new h to current h at the current order. ------------
  520 flotl = dfloat(l)
      etaq = one/((bias2*d/e)**(pt5/flotl) + addon)
      if ((nqindx .ne. 0) .or. (kflag .ne. 0)) go to 580
      etaqm1 = zero
      if (nq .eq. 1) go to 540
c compute ratio of new h to current h at the current order less one. ---
      d = zero
      do 530 i = 1,n
 530    d = d + (y(i,l)/ymax(i))**2
      edn = flotn*(tq(1)*eps)**2
      etaqm1 = one/((bias1*d/edn)**(pt5/(flotl - one)) + addon)
 540  etaqp1 = zero
      if (l .eq. lmax) go to 560
c compute ratio of new h to current h at current order plus one. -------
      cnquot = (tq(5)/conp)*(h/tau(2))**l
      d = zero
      do 550 i = 1,n
 550    d = d + ((error(i) - cnquot*y(i,lmax))/ymax(i))**2
      eup = flotn*(tq(3)*eps)**2
      etaqp1 = one/((bias3*d/eup)**(pt5/(flotl + one)) + addon)
 560  nqindx = 2
      if (etaq .ge. etaqp1) go to 570
      if (etaqp1 .gt. etaqm1) go to 600
      go to 590
 570  if (etaq .lt. etaqm1) go to 590
 580  if ((etaq .lt. thresh) .and. (kflag .eq. 0)) go to 690
      eta = etaq
      if ((kflag .le. -2) .and. (eta .gt. etamxf)) eta = etamxf
      go to 180
 590  if (etaqm1 .lt. thresh) go to 690
      call adjust (y, n0)
      l = nq
      nq = nq - 1
      eta = etaqm1
      nqindx = l
      go to 180
 600  if (etaqp1 .lt. thresh) go to 690
      nq = l
      eta = etaqp1
      l = l + 1
      do 610 i = 1,n
 610    y(i,l) = zero
      nqindx = l
      go to 180
c-----------------------------------------------------------------------
c control reaches this section if 3 or more consecutive failures
c have occurred.  it is assumed that the elements of the y array
c have accumulated errors of the wrong order.  the order is reduced
c by one, if possible.  then h is reduced by a factor of 0.1 and
c the step is retried.  after a total of 7 consecutive failures,
c an exit is taken with kflag = -2.
c-----------------------------------------------------------------------
 630  if (kflag .eq. kfh) go to 670
      if (nq .eq. 1) go to 640
      eta = etamin
      call adjust (y, n0)
      l = nq
      nq = nq - 1
      nqindx = l
      go to 180
  640 eta = dmax1(etamin,hmin/dabs(h))
      h = h*eta
      call diffun (n, t, y, save1)
      nfe = nfe + 1
      do 650 i = 1,n
 650    y(i,2) = h*save1(i)
      nqindx = 10
      go to 200
c-----------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
 660  kflag = -1
      go to 700
 670  kflag = -2
      go to 700
 680  kflag = -3
      go to 700
 690  etamax = etamx3
      if (nstep .le. 10) etamax = etamx2
 700  hold = h
      jstart = nq
      return
c-----------------------  end of subroutine tstep  ---------------------
      end
      subroutine coset
c-----------------------------------------------------------------------
c coset is called by tstep and sets coefficients for use there.
c
c for each order nq, the coefficients in el are calculated by use of
c  the generating polynomial lambda(x), with coefficients el(i):
c      lambda(x) = el(1) + el(2)*x + ... + el(nq+1)*(x**nq).
c for the backward differentiation formulas,
c                    nq
c      lambda(x) = product (1 + x/xi(i) ) .
c                   i = 1
c for the adams formulas,
c                              nq-1
c      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
c                              i = 1
c      lambda(-1) = 0,    lambda(0) = 1,
c where c is a normalization constant.
c in both cases, xi(i) is defined by
c      h*xi(i) = t sub n  -  t sub (n-i)
c              = h + tau(1) + tau(2) + ... tau(i-1).
c
c coset also sets maxder, the maximum order of the formulas
c available. currently this is 5 for the backward differentiation
c formulas, and 12 for the adams formulas.  to use different
c values (.le. 13),  change the numbers in statements 1 and 2 below.
c
c in addition to variables described previously, communication
c with coset uses the following..
c   tau    = a vector of length 13 containing the past nq values
c            of h.
c   el     = a vector of length 13 in which coset stores the
c            coefficients for the corrector formula.
c   tq     = a vector of length 5 in which coset stores constants
c            used for the convergence test, the error test, and
c            selection of h at a new order.
c   lmax   = maxder + 1, where maxder is the maximum order
c            available.  lmax is the maximum number of columns
c            of the y array to be used.
c   meth   = the basic method indicator.
c   nq     = the current order.
c   l      = nq + 1, the length of the vector stored in el, and
c            the number of columns of the y array being used.
c   nqindx = a counter controlling the frequency of order changes.
c            an order change is about to be considered if
c            nqindx = 1.
c-----------------------------------------------------------------------
c caution:  not all members of epcom1 are used in this subroutine.
c-----------------------------------------------------------------------
      integer jstart, kflag, l, lmax, meth, mf, n, nq, nqindx
      integer i, iback, j, jp1, jstart, kflag, l, lmax, maxder,
     1        meth, mf,n, nq, nqindx, nqm1
      double precision el, eps, h, hmax, hmin, ss, t, tau, tq,
     1                 uround
      double precision ahdss, cnqm1, csum, elp, em, em0, floti,
     1                 flotl, flotnq, hsum, hsum1, prod, rxi, s, xi
      double precision cortes
      double precision one, six, two, zero
      dimension em(13)
c
      common /epcom1/ t,h,hmin,hmax,eps,ss,uround,n,mf,kflag,jstart
      common /epcm10/ tau(13),el(13),tq(5),lmax,meth,nq,l,nqindx
      data cortes /0.1d0/
      data one  /1.0d0/, six /6.0d0/, two /2.0d0/, zero /0.0d0/
      ahdss = one
      if (ss .ne. zero) ahdss = dabs(h)/ss
      flotl = dfloat(l)
      nqm1 = nq - 1
      go to (1, 2), meth
 1    maxder = 12
      go to 100
c
 2    maxder = 5
      go to 200
c
 100  if (nq .ne. 1) go to 110
      el(1) = one
      el(2) = one
      tq(1) = one
      tq(2) = two*ahdss
      tq(3) = six*tq(2)
      tq(5) = one
      go to 300
 110  hsum = h
      em(1) = one
      flotnq = flotl - one
      do 115 i = 2,l
 115    em(i) = zero
      do 150 j = 1,nqm1
        if ((j .ne. nqm1) .or. (nqindx .ne. 1)) go to 130
        s = one
        csum = zero
        do 120 i = 1,nqm1
          csum = csum + s*em(i)/dfloat(i+1)
 120      s = -s
        tq(1) = ahdss*em(nqm1)/(flotnq*csum)
 130    rxi = h/hsum
        do 140 iback = 1,j
          i = (j + 2) - iback
 140      em(i) = em(i) + em(i-1)*rxi
 150    hsum = hsum + tau(j)
c compute integral from -1 to 0 of polynomial and of x times it. -------
      s = one
      em0 = zero
      csum = zero
      do 160 i = 1,nq
        floti = dfloat(i)
        em0 = em0 + s*em(i)/floti
        csum = csum + s*em(i)/(floti+1)
 160    s = -s
c in el, form coefficients of normalized integrated polynomial. --------
      s = one/em0
      el(1) = one
      do 170 i = 1,nq
  170   el(i+1) = s*em(i)/dfloat(i)
      xi = hsum/h
      tq(2) = ahdss*xi*em0/csum
      tq(5) = xi/el(l)
      if (nqindx .ne. 1) go to 300
c for higher order control constant, multiply polynomial by 1+x/xi(q). -
      rxi = one/xi
      do 180 iback = 1,nq
        i = (l + 1) - iback
 180    em(i) = em(i) + em(i-1)*rxi
c compute integral of polynomial. --------------------------------------
      s = one
      csum = zero
      do 190 i = 1,l
      csum = csum + s*em(i)/dfloat(i+1)
 190    s = -s
      tq(3) = ahdss*flotl*em0/csum
      go to 300
c
 200  do 210 i = 3,l
 210    el(i) = zero
      el(1) = one
      el(2) = one
      hsum = h
      hsum1 = zero
      prod = one
      rxi = one
      if (nq .eq. 1) go to 240
      do 230 j = 1,nqm1
c in el, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
        hsum = hsum + tau(j)
        hsum1 = hsum1 + tau(j)
        prod = prod*(hsum/hsum1)
        rxi = h/hsum
        jp1 = j + 1
        do 220 iback = 1,jp1
          i = (j + 3) - iback
 220      el(i) = el(i) + el(i-1)*rxi
 230    continue
 240  tq(2) = ahdss*el(2)*(one + prod)
      tq(5) = (one + prod)/el(l)
      if (nqindx .ne. 1) go to 300
      cnqm1 = rxi/el(l)
      elp = el(2) - rxi
      tq(1) = ahdss*elp/cnqm1
      hsum = hsum + tau(nq)
      rxi = h/hsum
      elp = el(2) + rxi
      tq(3) = ahdss*elp*rxi*(one + prod)*(flotl + one)
 300  tq(4) = cortes*tq(2)
      lmax = maxder + 1
      return
c-----------------------  end of subroutine coset  ---------------------
      end
      subroutine adjust (y, n0)
c-----------------------------------------------------------------------
c this subroutine adjusts the y array on reduction of order.
c see reference 1 for details.
c-----------------------------------------------------------------------
c caution:  not all members of epcom1 are used in this subroutine.
c-----------------------------------------------------------------------
      integer n0
      integer jstart, kflag, l, lmax, meth, mf, n, nq, nqindx
      integer i, iback, j, jp1, nqm1, nqm2
      double precision y
      double precision el, eps, h, hmax, hmin, ss, t, tau, tq, uround
      double precision hsum, xi
      double precision one, zero
      dimension y(n0,13)
c
      common /epcom1/ t,h,hmin,hmax,eps,ss,uround,n,mf,kflag,jstart
      common /epcm10/ tau(13),el(13),tq(5),lmax,meth,nq,l,nqindx
      data one /1.0d0/, zero /0.0d0/
      if (nq .eq. 2) return
      nqm1 = nq - 1
      nqm2 = nq - 2
      go to (100, 200), meth
c
 100  do 110 j = 1,lmax
 110    el(j) = zero
      el(2) = one
      hsum = zero
      do 130 j = 1,nqm2
c construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
        hsum = hsum + tau(j)
        xi = hsum/h
        jp1 = j + 1
        do 120 iback = 1,jp1
          i = (j + 3) - iback
 120      el(i) = el(i)*xi + el(i-1)
 130    continue
c construct coefficients of integrated polynomial. ---------------------
      do 140 j = 2,nqm1
  140   el(j+1) = dfloat(nq)*el(j)/dfloat(j)
      go to 300
c
 200  do 210 j = 1,lmax
 210    el(j) = zero
      el(3) = one
      hsum = zero
      do 230 j = 1,nqm2
c construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
        hsum = hsum + tau(j)
        xi = hsum/h
        jp1 = j + 1
        do 220 iback = 1,jp1
          i = (j + 4) - iback
 220      el(i) = el(i)*xi + el(i-1)
 230    continue
c
c subtract correction terms from y array. ------------------------------
 300  do 320 j = 3,nq
        do 310 i = 1,n
 310      y(i,j) = y(i,j) - y(i,l)*el(j)
 320    continue
      return
c-----------------------  end of subroutine adjust  --------------------
      end
      subroutine pset (diffun, pederv, y, n0, con, miter, ier)
c-----------------------------------------------------------------------
c pset is called by tstep to compute and to process the matrix
c p = i - (h/el(2))*j, where j is an approximation to the
c jacobian.  j is computed by either the user supplied
c subroutine pederv, when miter = 1, or by finite differences,
c when miter = 2.  j is stored in pw and replaced by p, using
c con = -h/el(2).  then p is subjected to an lu decomposition
c for later solution of linear algebraic systems with p as the
c coefficient matrix.
c
c in addition to variables described previously, communication
c with pset uses the following..
c   epsj = sqrt(uround), used in the numerical jacobian increments.
c   nsq  = n0**2.
c-----------------------------------------------------------------------
c caution:  not all epcom1 variables are used inthis subroutine.
c-----------------------------------------------------------------------
      integer ier, miter, n0
      integer ipiv, jstart, kflag, mf, n, nsq
      integer i, ier, j, j1, n
      double precision con, y
      double precision eps, epsj, h, hmax, hmin, pw, save1, save2,
     1                 ss, t, uround, ymax
      double precision d, r, r0, t, yj
      double precision one, rep, zero
      dimension y(n0,1)
c
      common /epcom1/ t,h,hmin,hmax,eps,ss,uround,n,mf,kflag,jstart
      common /epcom2/ ymax(1)
      common /epcom4/ save1(1)
      common /epcom5/ save2(1)
      common /epcom6/ pw(1)
      common /epcom7/ ipiv(1)
      common /epcom8/ epsj,nsq
      data one /1.0d0/, rep /1.0d-3/, zero /0.0d0/
      if (miter .eq. 2) go to 20
c if miter = 1, call pederv and multiply by a scalar.  -----------------
      call pederv (n, t, y, pw, n0)
      do 10 i = 1,nsq
  10    pw(i) = pw(i)*con
      go to 60
c if miter = 2, make n calls to diffun to approximate j. ---------------
  20  d = zero
      do 30 i = 1,n
  30    d = d + save2(i)**2
      r0 = dabs(h)*dsqrt(d)*uround/rep
      j1 = 0
      do 50 j = 1,n
        yj = y(j,1)
        r = epsj*ymax(j)
        r = dmax1(r,r0)
        y(j,1) = y(j,1) + r
        d = con/r
        call diffun (n, t, y, save1)
        do 40 i = 1,n
  40      pw(i+j1) = (save1(i) - save2(i))*d
        y(j,1) = yj
        j1 = j1 + n0
  50    continue
c add on the identity matrix.  -----------------------------------------
  60  j = 1
      do 70 i = 1,n
        pw(j) = pw(j) + one
  70    j = j + (n0 + 1)
c get lu decomposition of p. -------------------------------------------
      call dec (n, n0, pw, ipiv, ier)
      return
c-----------------------  end of subroutine pset -----------------------
      end
      subroutine dec (n, ndim, a, ip, ier)
c-----------------------------------------------------------------------
c matrix triangularization by gaussian elimination.
c input..
c    n = order of matrix.
c    ndim = declared dimension of array  a .
c    a = matrix to be triangularized.
c output..
c    a(i,j), i.le.j = upper triangular factor, u .
c    a(i,j), i.gt.j = multipliers = lower triangular factor, i - l.
c    ip(k), k.lt.n = index of k-th pivot row.
c    ip(n) = (-1)**(number of interchanges) or o .
c    ier = 0 if a nonsingular, or k if a found to be
c          singular at stage k.
c use  sol  to obtain solution of linear system.
c determ(a) = ip(n)*a(1,1)*a(2,2)*...*a(n,n).
c if ip(n)=0, a is singular, sol will divide by zero.
c interchanges finished in u , only partly in l .
c
c reference..
c c. b. moler, algorithm 423, linear equation solver,
c comm. assoc. comput. mach., 15 (1972), p. 274.
c-----------------------------------------------------------------------
      integer ier, ip, n, ndim
      integer i, j, k, kp1, m, nm1
      double precision a
      double precision t
      double precision one, zero
      dimension a(ndim,n),ip(n)
      data one /1.0d0/, zero /0.0d0/
      ier = 0
      ip(n) = 1
      if (n .eq. 1) go to 70
      nm1 = n - 1
      do 60 k = 1,nm1
        kp1 = k + 1
        m = k
        do 10 i = kp1,n
  10      if (dabs(a(i,k)) .gt. dabs(a(m,k))) m = i
        ip(k) = m
        t = a(m,k)
        if (m .eq. k) go to 20
        ip(n) = -ip(n)
        a(m,k) = a(k,k)
        a(k,k) = t
 20     if (t .eq. zero) go to 80
        t = one/t
        do 30 i = kp1,n
 30       a(i,k) = -a(i,k)*t
        do 50 j = kp1,n
          t = a(m,j)
          a(m,j) = a(k,j)
          a(k,j) = t
          if (t .eq. zero) go to 50
          do 40 i = kp1,n
 40         a(i,j) = a(i,j) + a(i,k)*t
 50       continue
 60     continue
 70   k = n
      if (a(n,n) .eq. zero) go to 80
      return
 80   ier = k
      ip(n) = 0
      return
c-----------------------  end of subroutine dec  -----------------------
      end
      subroutine sol (n, ndim, a, b, ip)
c-----------------------------------------------------------------------
c solution of linear system, a*x = b .
c input..
c   n = order of matrix.
c   ndim = declared dimension of array  a .
c   a = triangularized matrix obtained from dec.
c   b = right hand side vector.
c   ip = pivot vector obtained from dec.
c do not use if dec has set ier .ne. 0.
c output..
c   b = solution vector, x .
c-----------------------------------------------------------------------
      integer ip, n, ndim
      integer i, k, kb, km1, kp1, m, nm1
      double precision a, b
      double precision t
      dimension a(ndim, n), b(n), ip(n)
c
      if (n .eq. 1) go to 50
      nm1 = n - 1
      do 20 k = 1,nm1
        kp1 = k + 1
        m = ip(k)
        t = b(m)
        b(m) = b(k)
        b(k) = t
        do 10 i = kp1,n
 10       b(i) = b(i) + a(i,k)*t
 20     continue
      do 40 kb = 1,nm1
        km1 = n - kb
        k = km1 + 1
        b(k) = b(k)/a(k,k)
        t = -b(k)
        do 30 i = 1,km1
 30       b(i) = b(i) + a(i,k)*t
 40     continue
 50   b(1) = b(1)/a(1,1)
      return
c-----------------------  end of subroutine sol  -----------------------
      end
