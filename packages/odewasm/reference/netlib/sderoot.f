      subroutine deroot(f,neqn,y,t,tout,relerr,abserr,iflag,
     *                  g,reroot,aeroot)
      external f
      integer neqn, iflag
      real y(neqn), t, tout, relerr, abserr, g,
     *                 reroot, aeroot
c
c   subroutine  deroot  integrates a system of up to 20 first order
c   ordinary differential equations of the form
c             dy(i)/dt = f(t,y(1),...,y(neqn))
c             y(i) given at t.
c   the subroutine integrates from  t  in the direction of  tout  until
c   it locates the first root of the nonlinear equation
c         g(t,y(1),...,y(neqn),yp(1),...,yp(neqn)) = 0.
c   upon finding the root, the code returns with all parameters in the
c   call list set for continuing the integration to the next root or
c   the first root of a new function  g  .  if no root is found, the
c   integration proceeds to  tout .  again all parameters are set to
c   continue.
c
c   the differential equations are actually solved by a suite of codes,
c   deroot , step , and  intrp .  deroot  is a supervisor which directs
c   the integration.  it calls on  step  to advance the solution and
c   intrp  to interpolate the solution and its derivative.  step  uses
c   a modified divided difference form of the adams pece formulas and
c   local extrapolation.  it adjusts the order and step size to control
c   the local error per unit step in a generalized sense.  normally each
c   call to  step  advances the solution one step in the direction of
c   tout .  for reasons of efficiency  deroot  integrates beyond  tout
c   internally, though never beyond t+10*(tout-t), and calls  intrp  to
c   interpolate the solution and derivative at  tout .  an option is
c   provided to stop the integration at  tout  but it should be used
c   only if it is impossible to continue the integration beyond  tout .
c
c   after each internal step,  deroot  evaluates the function  g  and
c   checks for a change in sign in the function value from the
c   preceding step.  such a change indicates a root lies in the
c   interval of the step just completed.  deroot  then calls subroutine
c   root  to reduce the bracketing interval until the root is
c   determined to the desired accuracy.  subroutine  root  uses a
c   combination of the secant rule and bisection to do this.  the
c   solution and derivative values required are obtained by
c   interpolation with  intrp .  the code locates only those roots
c   for which  g  changes sign in  (t,tout)  and for which a
c   bracketing interval exists.  in particular, it does not locate
c   a root at the initial point  t .
c
c   the codes  step  and  intrp  and that portion of  deroot  which
c   directs the integration are completely explained and documented in
c   the text, computer solution of ordinary differential equations,
c   the initial value problem by l. f. shampine and m. k. gordon.
c   subroutine  root  is a slightly modified version of the root-solver
c   discussed in the text, numerical computing, an introduction by
c   l. f. shampine and r. c. allen.
c
c   the parameters for deroot are
c      f -- subroutine f(t,y,yp) to evaluate derivatives yp(i)=dy(i)/dt
c      neqn -- number of equations to be integrated
c      y(*) -- solution vector at  t
c      t -- independent variable
c      tout -- arbitrary point beyond the root desired
c      relerr,abserr -- relative and absolute error tolerances for local
c           error test.  at each step the code requires
c             abs(local error) .le. abs(y)*relerr + abserr
c           for each component of the local error and solution vectors
c      iflag -- indicates status of integration
c      g - function of t, y(*), yp(*) whose root is desired.
c      reroot, aeroot -- relative and absolute error tolerances for
c           accepting the root.  the interval containing the root is
c           reduced until it satisfies
c            0.5*abs(length of interval) .le. reroot*abs(root)+aeroot
c           where root is that endpoint yielding the smaller value of
c           g  in magnitude.  pure relative error is not recommended
c           if the root might be zero.
c
c   first call to  deroot  --
c
c   the user must provide storage in his calling program for the
c   array in the call list,
c              y(neqn)
c   and declare  f , g  in an external statement.  he must supply the
c   subroutine  f(t,y,yp)  to evaluate
c           dy(i)/dt = yp(i) = f(t,y(1),...,y(neqn))
c   and the function  g(t,y,yp)  to evaluate
c           g = g(t,y(1),...,y(neqn),yp(1),...,yp(neqn)).
c   note that the array yp is an input argument and should not be
c   computed in the function subprogram.  finally the user must
c   initialize the parameters
c      neqn -- number of equations to be integrated
c      y(*) -- vector of initial conditions
c      t -- starting point of integration
c      tout -- arbitrary point beyond the root desired
c      relerr,abserr -- relative and absolute local error tolerances
c                       for integrating the equations
c      iflag -- +1,-1.  indicator to initialize the code.  normal input
c           is +1.  the user should set iflag=-1 only if it is
c           impossible to continue the integration beyond  tout .
c      reroot,aeroot -- relative and absolute error tolerances for
c                       computing the root of  g
c
c   all parameters except f, g, neqn, tout, reroot and aeroot may be
c   altered by the code on output so must be variables in the calling
c   program.
c
c   output from  deroot  --
c
c      neqn -- unchanged
c      y(*) -- solution at  t
c      t -- last point reached in integration.  normal return has
c           t = tout or t = root
c      tout -- unchanged
c      relerr,abserr -- normal return has tolerances unchanged.  iflag=3
c           signals tolerances increased
c      iflag = 2 -- normal return.  integration reached  tout
c            = 3 -- integration did not reach  tout  because error
c                   tolerances too small.  relerr ,  abserr  increased
c                   appropriately for continuing
c            = 4 -- integration did not reach  tout  because more than
c                   maxnum  steps needed
c            = 5 -- integration did not reach  tout  because equations
c                   appear to be stiff
c            = 6 -- invalid input parameters (fatal error)
c            = 7 -- normal return.  a root was found which satisfied
c                   the error criterion or had a zero residual
c            = 8 -- abnormal return.  an odd order pole of  g  was
c                   found.
c            = 9 -- abnormal return.  too many evaluations of  g  were
c                   required (as programmed 500 are allowed.)
c           the value of  iflag  is returned negative when the input
c           value is negative and the integration does not reach
c           tout , i.e., -3,-4,-5,-7,-8,-9.
c      reroot,aeroot -- unchanged
c
c   subsequent calls to  deroot  --
c
c   subroutine  deroot  returns with all information needed to continue
c   the integration.  if the integration did not reach  tout  and the
c   user wants to continue, he just calls again.  if the integration
c   reached  tout , the user need only define a new  tout  and call
c   again.  the output value of  iflag  is the appropriate input value
c   for subsequent calls.  the only situation in which it should be
c   altered is to stop the integration internally at the new  tout ,
c   i.e., change output  iflag=2  to input  iflag=-2 .  only the error
c   tolerances and the function  g  may be changed by the user before
c   continuing.  all other parameters must remain unchanged.
c
      logical stiff,crash,start
      real fouru, eps, gxold, gx, x, delsgn,
     *                 b, c, gc, del, tend, releps, abseps, troot,
     *                 absdel, h, hold, told, gt
      real yy(20),wt(20),phi(20,16),p(20),yp(20),ypout(20),
     *                 psi(12)
      real abs, amax1, amin1, r1mach, sign                              ***
c
c***********************************************************************
c*  the only machine dependent constant is based on the machine unit   *
c*  roundoff error  u  which is the smallest positive number such that *
c*  1.0+u .gt. 1.0 .  u  must be calculated and  fouru=4.0*u  inserted *
c*  in the following statement before using  deroot .  the subroutine  *
c*  machin  calculates  u .  fouru  and  twou=2.0*u  must also be      *
c*  inserted in subroutine  step  before calling  deroot .             *
c     data fouru/8.8d-16/                                               ***
c***********************************************************************
c
c   the constant  maxnum  is the maximum number of steps allowed in one
c   call to  deroot .  the user may change this limit by altering the
c   following statement
      data maxnum/500/
c
c   this version of  deroot  allows a maximum of 20 equations.  to
c   increase this number, only the number 20 in the dimension statement
c   and in the following statement need be changed
c            ***            ***            ***
c   test for improper parameters
c
c-----------------------------------------------------------------
      fouru = 4.0*r1mach(4)                                             ***
      if(neqn .lt. 1  .or.  neqn .gt. 20) go to 10
      if(t .eq. tout) go to 10
      if(relerr .lt. 0.0  .or.  abserr .lt. 0.0) go to 10
      eps = amax1(relerr,abserr)
      if(eps .le. 0.0) go to 10
      if(reroot .lt. 0.0  .or.  aeroot .lt. 0.0) go to 10
      if(reroot+aeroot .le. 0.0) go to 10
      if(iflag .eq. 0) go to 10
      isn = isign(1,iflag)
      iflag = iabs(iflag)
      if(iflag .eq. 1) go to 20
      if(t .ne. told) go to 10
      if(iflag .ge. 2  .and.  iflag .le. 5) go to 15
      if(iflag .ge. 7  .and.  iflag .le. 9) go to 15
   10 iflag = 6
      return
c
c   for a new function g, check for a root in interval of step
c   just completed
c
   15 gxold = gx
      gx = g(x,yy,yp)
      if(gx .eq. gxold) go to 20
      if(iflag .gt. 2  .and.  iflag .le. 5) go to 20
      if(isnold .lt. 0  .or.  delsgn*(tout-t) .lt. 0.0) go to 20
      jflag = 1
      b = x
      c = t
      call intrp(x,yy,c,y,ypout,neqn,kold,phi,psi)
      gc = g(c,y,ypout)
      if(sign(1.0,gc)*sign(1.0,gx) .lt. 0.0) go to 140
      if(gc .eq. 0.0  .or.  gx .eq. 0.0) go to 140
c
c   on each call set interval of integration and counter for number of
c   steps.  adjust input error tolerances to define weight vector for
c   subroutine  step
c
   20 del = tout - t
      absdel = abs(del)
      tend = t + 10.0*del
      if(isn .lt. 0) tend = tout
      nostep = 0
      kle4 = 0
      stiff = .false.
      releps = relerr/eps
      abseps = abserr/eps
      if(iflag .eq. 1) go to 30
      if(isnold .lt. 0) go to 30
      if(delsgn*del .gt. 0.0) go to 50
c
c   on start and restart also set work variables x and yy(*), store the
c   direction of integration, initialize the step size, and evaluate  g
c
   30 start = .true.
      x = t
      troot = t
      do 40 l = 1,neqn
   40   yy(l) = y(l)
      delsgn = sign(1.0,del)
      h = sign(amax1(abs(tout-x),fouru*abs(x)),tout-x)
      call f(x,yy,yp)
      gx = g(x,yy,yp)
c
c   if already past output point, interpolate and return
c
   50 gxold = gx
      if(abs(x-t) .lt. absdel) go to 60
      call intrp(x,yy,tout,y,ypout,neqn,kold,phi,psi)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
c
c   if cannot go past output point and sufficiently close,
c   extrapolate and return
c
   60 if(isn .gt. 0  .or.  abs(tout-x) .ge. fouru*abs(x)) go to 80
      h = tout - x
      call f(x,yy,yp)
      do 70 l = 1,neqn
   70   y(l) = yy(l) + h*yp(l)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
c
c   test for too much work
c
   80 if(nostep .lt. maxnum) go to 100
      iflag = isn*4
      if(stiff) iflag = isn*5
      do 90 l = 1,neqn
   90   y(l) = yy(l)
      t = x
      told = t
      isnold = 1
      return
c
c   limit step size, set weight vector and take a step
c
  100 h = sign(amin1(abs(h),abs(tend-x)),h)
      do 110 l = 1,neqn
  110   wt(l) = releps*abs(yy(l)) + abseps
      call step(x,yy,f,neqn,h,eps,wt,start,
     *  hold,k,kold,crash,phi,p,yp,psi)
c
c   test for tolerances too small.  if so, set the derivative at x
c   before returning
c
      if(.not.crash) go to 130
      iflag = isn*3
      relerr = eps*releps
      abserr = eps*abseps
      do 120 l = 1,neqn
        yp(l) = phi(l,1)
  120   y(l) = yy(l)
      t = x
      told = t
      isnold = 1
      return
c
c   augment counter on work and test for stiffness.  also test for a
c   root in the step just completed
c
  130 nostep = nostep + 1
      kle4 = kle4 + 1
      if(kold .gt. 4) kle4 = 0
      if(kle4 .ge. 50) stiff = .true.
      gx = g(x,yy,yp)
      if(sign(1.0,gxold)*sign(1.0,gx) .lt. 0.0) go to 135
      if(gx .eq. 0.0  .or.  gxold .eq. 0.0) go to 135
      go to 50
c
c   locate root of g.  interpolate with  intrp  for solution and
c   derivative values
c
  135 b = x
      c = x - hold
      jflag = 1
  140 call root(t,gt,b,c,reroot,aeroot,jflag)
      if(jflag .gt. 0) go to 150
      call intrp(x,yy,t,y,ypout,neqn,kold,phi,psi)
      gt = g(t,y,ypout)
      go to 140
  150 continue
      iflag = jflag + 6
      if(jflag .eq. 2  .or.  jflag .eq. 4) iflag = 7
      if(jflag .eq. 3) iflag = 8
      if(jflag .eq. 5) iflag = 9
      iflag = iflag*isn
      call intrp(x,yy,b,y,ypout,neqn,kold,phi,psi)
      t = b
      if(abs(t-troot) .le. reroot*abs(t) + aeroot) go to 50
      troot = t
      told = t
      isnold = 1
      return
      end
      subroutine intrp(x,y,xout,yout,ypout,neqn,kold,phi,psi)
      integer neqn, kold
      real x, y(neqn), xout, yout(neqn), ypout(neqn),
     *                 phi(20,16), psi(12)
c
c   the methods in subroutine  step  approximate the solution near  x
c   by a polynomial.  subroutine  intrp  approximates the solution at
c   xout  by evaluating the polynomial there.  information defining this
c   polynomial is passed from  step  so  intrp  cannot be used alone.
c
c   this code is completely explained and documented in the text,
c   computer solution of ordinary differential equations   the initial
c   value problem  by l. f. shampine and m. k. gordon.
c
c   input to intrp --
c
c   the user provides storage in the calling program for the arrays in
c   the call list
c   and defines
c      xout -- point at which solution is desired.
c   the remaining parameters are defined in  step  and passed to  intrp
c   from that subroutine
c
c   output from  intrp --
c
c      yout(*) -- solution at  xout
c      ypout(*) -- derivative of solution at  xout
c   the remaining parameters are returned unaltered from their input
c   values.  integration with  step  may be continued.
c
      real hi, temp1, term, psijm1, gamma, eta, temp2,
     *                 temp3
      real g(13),w(13),rho(13)
      data g(1)/1.0/,rho(1)/1.0/
c
      hi = xout - x
      ki = kold + 1
      kip1 = ki + 1
c
c   initialize w(*) for computing g(*)
c
      do 5 i = 1,ki
        temp1 = i
    5   w(i) = 1.0/temp1
      term = 0.0
c
c   compute g(*)
c
      do 15 j = 2,ki
        jm1 = j - 1
        psijm1 = psi(jm1)
        gamma = (hi + term)/psijm1
        eta = hi/psijm1
        limit1 = kip1 - j
        do 10 i = 1,limit1
   10     w(i) = gamma*w(i) - eta*w(i+1)
        g(j) = w(1)
        rho(j) = gamma*rho(jm1)
   15   term = psijm1
c
c   interpolate
c
      do 20 l = 1,neqn
        ypout(l) = 0.0
   20   yout(l) = 0.0
      do 30 j = 1,ki
        i = kip1 - j
        temp2 = g(i)
        temp3 = rho(i)
        do 25 l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
   25     ypout(l) = ypout(l) + temp3*phi(l,i)
   30   continue
      do 35 l = 1,neqn
   35   yout(l) = y(l) + hi*yout(l)
      return
      end
      subroutine root(t,ft,b,c,relerr,abserr,iflag)
      real t, ft, b, c, relerr, abserr
      integer iflag
c
c  root computes a root of the nonlinear equation f(x)=0
c  where f(x) is a continuous real function of a single real
c  variable x.  the method used is a combination of bisection
c  and the secant rule.
c
c  normal input consists of a continuous function f and an
c  interval (b,c) such that f(b)*f(c).le.0.0.  each iteration
c  finds new values of b and c such that the interval (b,c) is
c  shrunk and f(b)*f(c).le.0.0.  the stopping criterion is
c
c         abs(b-c).le.2.0*(relerr*abs(b)+abserr)
c
c  where relerr=relative error and abserr=absolute error are
c  input quantities.  set the flag, iflag, positive to initialize
c  the computation.  as b,c and iflag are used for both input and
c  output, they must be variables in the calling program.
c
c  if 0 is a possible root, one should not choose abserr=0.0.
c
c  the output value of b is the better approximation to a root
c  as b and c are always redefined so that abs(f(b)).le.abs(f(c)).
c
c  to solve the equation, root must evaluate f(x) repeatedly. this
c  is done in the calling program.  when an evaluation of f is
c  needed at t, root returns with iflag negative.  evaluate ft=f(t)
c  and call root again.  do not alter iflag.
c
c  when the computation is complete, root returns to the calling
c  program with iflag positive:
c
c     iflag=1  if f(b)*f(c).lt.0 and the stopping criterion is met.
c
c          =2  if a value b is found such that the computed value
c              f(b) is exactly zero.  the interval (b,c) may not
c              satisfy the stopping criterion.
c
c          =3  if abs(f(b)) exceeds the input values abs(f(b)),
c              abs(f(c)).   in this case it is likely that b is close
c              to a pole of f.
c
c          =4  if no odd order root was found in the interval.  a
c              local minimum may have been obtained.
c
c          =5  if too many function evaluations were made.
c              (as programmed, 500 are allowed.)
c
c  this code is a modification of the code  zeroin  which is completely
c  explained and documented in the text,  numerical computing:  an
c  introduction  by l. f. shampine and r. c. allen.
c
      real u, re, ae, acbs, a, fa, fb, fc, fx, cmb,
     *                 acmb, tol, p, q
      real abs, amax1, r1mach, sign
c***********************************************************************
c*  the only machine dependent constant is based on the machine unit   *
c*  roundoff error  u  which is the smallest positive number such that *
c*  1.0+u .gt. 1.0 .  u  must be calculated and inserted in the        *
c*  following data statement before using  root .  the routine  machin *
c*  calculates  u .                                                    *
c     data u/2.2d-16/
c***********************************************************************
c
      u = r1mach(4)
      if(iflag.ge.0) go to 100
      iflag=iabs(iflag)
      go to (200,300,400), iflag
  100 re=amax1(relerr,u)
      ae=amax1(abserr,0.0)
      ic=0
      acbs=abs(b-c)
      a=c
      t=a
      iflag=-1
      return
  200 fa=ft
      t=b
      iflag=-2
      return
  300 fb=ft
      fc=fa
      kount=2
      fx=amax1(abs(fb),abs(fc))
    1 if(abs(fc).ge.abs(fb))go to 2
c
c  interchange b and c so that abs(f(b)).le.abs(f(c)).
c
      a=b
      fa=fb
      b=c
      fb=fc
      c=a
      fc=fa
    2 cmb=0.5*(c-b)
      acmb=abs(cmb)
      tol=re*abs(b)+ae
c
c  test stopping criterion and function count.
c
      if(acmb.le.tol)go to 8
      if(kount.ge.500)go to 12
c
c  calculate new iterate implicitly as b+p/q
c  where we arrange p.ge.0.  the implicit
c  form is used to prevent overflow.
c
      p=(b-a)*fb
      q=fa-fb
      if(p.ge.0.0)go to 3
      p=-p
      q=-q
c
c  update a, check if reduction in the size of bracketing
c  interval is satisfactory.  if not, bisect until it is.
c
    3 a=b
      fa=fb
      ic=ic+1
      if(ic.lt.4)go to 4
      if(8.0*acmb.ge.acbs)go to 6
      ic=0
      acbs=acmb
c
c  test for too small a change.
c
    4 if(p.gt.abs(q)*tol)go to 5
c
c  increment by tolerance.
c
      b=b+sign(tol,cmb)
      go to 7
c
c  root ought to be between b and (c+b)/2.
c
    5 if(p.ge.cmb*q)go to 6
c
c  use secant rule.
c
      b=b+p/q
      go to 7
c
c  use bisection.
c
    6 b=0.5*(c+b)
c
c  have completed computation for new iterate b.
c
    7 t=b
      iflag=-3
      return
  400 fb=ft
      if(fb.eq.0.0)go to 9
      kount=kount+1
      if(sign(1.0,fb).ne.sign(1.0,fc))go to 1
      c=a
      fc=fa
      go to 1
c
c  finished.  set iflag.
c
    8 if(sign(1.0,fb).eq.sign(1.0,fc))go to 11
      if(abs(fb).gt.fx)go to 10
      iflag=1
      return
    9 iflag=2
      return
   10 iflag=3
      return
   11 iflag=4
      return
   12 iflag=5
      return
      end
      subroutine step(x,y,f,neqn,h,eps,wt,start,
     *                hold,k,kold,crash,phi,p,yp,psi)
      external f
      integer neqn, k, kold
      logical start, crash
      real x, y(neqn), h, eps, wt(20), hold,
     *                 phi(20,16), p(20), yp(20), psi(12)
c
c   subroutine  step  integrates a system of first order ordinary
c   differential equations one step, normally from x to x+h, using a
c   modified divided difference form of the adams pece formulas.  local
c   extrapolation is used to improve absolute stability and accuracy.
c   the code adjusts its order and step size to control the local error
c   per unit step in a generalized sense.  special devices are included
c   to control roundoff error and to detect when the user is requesting
c   too much accuracy.
c
c   this code is completely explained and documented in the text,
c   computer solution of ordinary differential equations   the initial
c   value problem  by l. f. shampine and m. k. gordon.
c
c
c   the parameters represent
c      x -- independent variable
c      y(*) -- solution vector at x
c      yp(*) -- derivative of solution vector at  x  after successful
c           step
c      neqn -- number of equations to be integrated
c      h -- appropriate step size for next step.  normally determined by
c           code
c      eps -- local error tolerance.  must be variable
c      wt(*) -- vector of weights for error criterion
c      start -- logical variable set .true. for first step,  .false.
c           otherwise
c      hold -- step size used for last successful step
c      k -- appropriate order for next step (determined by code)
c      kold -- order used for last successful step
c      crash -- logical variable set .true. when no step can be taken,
c           .false. otherwise.
c   the arrays  phi, psi  are required for the interpolation subroutine
c   intrp .  the array  p  is internal to the code.
c
c   input to  step
c
c      first call --
c
c   the user must provide storage in his driver program for all arrays
c   in the call list, namely
c
c     dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
c
c   the user must also declare  start  and  crash  logical variables
c   and  f  an external subroutine, supply the subroutine  f(x,y,yp)
c   to evaluate
c      dy(i)/dx = yp(i) = f(x,y(1),y(2),...,y(neqn))
c   and initialize only the following parameters
c      x -- initial value of the independent variable
c      y(*) -- vector of initial values of dependent variables
c      neqn -- number of equations to be integrated
c      h -- nominal step size indicating direction of integration
c           and maximum size of step.  must be variable
c      eps -- local error tolerance per step.  must be variable
c      wt(*) -- vector of non-zero weights for error criterion
c      start -- .true.
c
c   step  requires the l2 norm of the vector with components
c   local error(l)/wt(l)  be less than  eps  for a successful step.  the
c   array  wt  allows the user to specify an error test appropriate
c   for his problem.  for example,
c      wt(l) = 1.0  specifies absolute error,
c            = abs(y(l)) error relative to the most recent value of the
c                 l-th component of the solution,
c            = abs(yp(l))  error relative to the most recent value of
c                 the l-th component of the derivative,
c            = amax1(wt(l),abs(y(l)))  error relative to the largest
c                 magnitude of l-th component obtained so far,
c            = abs(y(l))*relerr/eps + abserr/eps  specifies a mixed
c                 relative-absolute test where  relerr  is relative
c                 error,  abserr  is absolute error and  eps =
c                 amax1(relerr,abserr) .
c
c      subsequent calls --
c
c   subroutine  step  is designed so that all information needed to
c   continue the integration, including the step size  h  and the order
c   k , is returned with each step.  with the exception of the step
c   size, the error tolerance, and the weights, none of the parameters
c   should be altered.  the array  wt  must be updated after each step
c   to maintain relative error tests like those above.  normally the
c   integration is continued just beyond the desired endpoint and the
c   solution interpolated there with subroutine  intrp .  if it is
c   impossible to integrate beyond the endpoint, the step size may be
c   reduced to hit the endpoint since the code will not take a step
c   larger than the  h  input.  changing the direction of integration,
c   i.e., the sign of  h , requires the user set  start = .true. before
c   calling  step  again.  this is the only situation in which  start
c   should be altered.
c
c   output from  step
c
c      successful step --
c
c   the subroutine returns after each successful step with  start  and
c   crash  set .false. .  x  represents the independent variable
c   advanced one step of length  hold  from its value on input and  y
c   the solution vector at the new value of  x .  all other parameters
c   represent information corresponding to the new  x  needed to
c   continue the integration.
c
c      unsuccessful step --
c
c   when the error tolerance is too small for the machine precision,
c   the subroutine returns without taking a step and  crash = .true. .
c   an appropriate step size and error tolerance for continuing are
c   estimated and all other information is restored as upon input
c   before returning.  to continue with the larger tolerance, the user
c   just calls the code again.  a restart is neither required nor
c   desirable.
c
      logical phase1,nornd
      real twou, fouru, p5eps, round, sum, absh, realns,
     *                 temp1, temp2, reali, temp3, temp4, temp5, temp6,
     *                 tau, xold, erkm1, erkm2, erk, err, rho, erkp1,
     *                 r, hnew
      real alpha(12),beta(12),sig(13),w(12),v(12),g(13),
     *                 gstr(13),two(13)
      real abs, amax1, amin1, r1mach, sign, sqrt
c***********************************************************************
c*  the only machine dependent constants are based on the machine unit *
c*  roundoff error  u  which is the smallest positive number such that *
c*  1.0+u .gt. 1.0  .  the user must calculate  u  and insert          *
c*  twou=2.0*u  and  fouru=4.0*u  in the data statement before calling *
c*  the code.  the routine  machin  calculates  u .                    *
c     data twou,fouru/4.4d-16, 8.8d-16/                                 ***
c***********************************************************************
      data two/2.0,4.0,8.0,16.0,32.0,64.0,128.0,256.0,
     *  512.0,1024.0,2048.0,4096.0,8192.0/
      data gstr/0.500,0.0833,0.0417,0.0264,0.0188,0.0143,
     *  0.0144,0.00936,0.00789,0.00679,0.00592,
     *  0.00524,0.00468/
      data g(1),g(2)/1.0,0.5/,sig(1)/1.0/
c
c
c       ***     begin block 0     ***
c   check if step size or error tolerance is too small for machine
c   precision.  if first step, initialize phi array and estimate a
c   starting step size.
c                   ***
c
c   if step size is too small, determine an acceptable one
c
      twou = 2.0*r1mach(4)                                              ***
      fouru = 2.0*twou                                                  ***
      crash = .true.
      if(abs(h) .ge. fouru*abs(x)) go to 5
      h = sign(fouru*abs(x),h)
      return
    5 p5eps = 0.5*eps
c
c   if error tolerance is too small, increase it to an acceptable value
c
      round = 0.0
      do 10 l = 1,neqn
   10   round = round + (y(l)/wt(l))**2
      round = twou*sqrt(round)
      if(p5eps .ge. round) go to 15
      eps = 2.0*round*(1.0 + fouru)
      return
   15 crash = .false.
      if(.not.start) go to 99
c
c   initialize.  compute appropriate step size for first step
c
      call f(x,y,yp)
      sum = 0.0
      do 20 l = 1,neqn
        phi(l,1) = yp(l)
        phi(l,2) = 0.0
   20   sum = sum + (yp(l)/wt(l))**2
      sum = sqrt(sum)
      absh = abs(h)
      if(eps .lt. 16.0*sum*h*h) absh = 0.25*sqrt(eps/sum)
      h = sign(amax1(absh,fouru*abs(x)),h)
      hold = 0.0
      k = 1
      kold = 0
      start = .false.
      phase1 = .true.
      nornd = .true.
      if(p5eps .gt. 100.0*round) go to 99
      nornd = .false.
      do 25 l = 1,neqn
   25   phi(l,15) = 0.0
   99 ifail = 0
c       ***     end block 0     ***
c
c       ***     begin block 1     ***
c   compute coefficients of formulas for this step.  avoid computing
c   those quantities not changed when step size is not changed.
c                   ***
c
  100 kp1 = k+1
      kp2 = k+2
      km1 = k-1
      km2 = k-2
c
c   ns is the number of steps taken with size h, including the current
c   one.  when k.lt.ns, no coefficients change
c
      if(h .ne. hold) ns = 0
      if (ns .le. kold) ns=ns+1
      nsp1 = ns+1
      if (k .lt. ns) go to 199
c
c   compute those components of alpha(*),beta(*),psi(*),sig(*) which
c   are changed
c
      beta(ns) = 1.0
      realns = ns
      alpha(ns) = 1.0/realns
      temp1 = h*realns
      sig(nsp1) = 1.0
      if(k .lt. nsp1) go to 110
      do 105 i = nsp1,k
        im1 = i-1
        temp2 = psi(im1)
        psi(im1) = temp1
        beta(i) = beta(im1)*psi(im1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        reali = i
  105   sig(i+1) = reali*alpha(i)*sig(i)
  110 psi(k) = temp1
c
c   compute coefficients g(*)
c
c   initialize v(*) and set w(*).  g(2) is set in data statement
c
      if(ns .gt. 1) go to 120
      do 115 iq = 1,k
        temp3 = iq*(iq+1)
        v(iq) = 1.0/temp3
  115   w(iq) = v(iq)
      go to 140
c
c   if order was raised, update diagonal part of v(*)
c
  120 if(k .le. kold) go to 130
      temp4 = k*kp1
      v(k) = 1.0/temp4
      nsm2 = ns-2
      if(nsm2 .lt. 1) go to 130
      do 125 j = 1,nsm2
        i = k-j
  125   v(i) = v(i) - alpha(j+1)*v(i+1)
c
c   update v(*) and set w(*)
c
  130 limit1 = kp1 - ns
      temp5 = alpha(ns)
      do 135 iq = 1,limit1
        v(iq) = v(iq) - temp5*v(iq+1)
  135   w(iq) = v(iq)
      g(nsp1) = w(1)
c
c   compute the g(*) in the work vector w(*)
c
  140 nsp2 = ns + 2
      if(kp1 .lt. nsp2) go to 199
      do 150 i = nsp2,kp1
        limit2 = kp2 - i
        temp6 = alpha(i-1)
        do 145 iq = 1,limit2
  145     w(iq) = w(iq) - temp6*w(iq+1)
  150   g(i) = w(1)
  199   continue
c       ***     end block 1     ***
c
c       ***     begin block 2     ***
c   predict a solution p(*), evaluate derivatives using predicted
c   solution, estimate local error at order k and errors at orders k,
c   k-1, k-2 as if constant step size were used.
c                   ***
c
c   change phi to phi star
c
      if(k .lt. nsp1) go to 215
      do 210 i = nsp1,k
        temp1 = beta(i)
        do 205 l = 1,neqn
  205     phi(l,i) = temp1*phi(l,i)
  210   continue
c
c   predict solution and differences
c
  215 do 220 l = 1,neqn
        phi(l,kp2) = phi(l,kp1)
        phi(l,kp1) = 0.0
  220   p(l) = 0.0
      do 230 j = 1,k
        i = kp1 - j
        ip1 = i+1
        temp2 = g(i)
        do 225 l = 1,neqn
          p(l) = p(l) + temp2*phi(l,i)
  225     phi(l,i) = phi(l,i) + phi(l,ip1)
  230   continue
      if(nornd) go to 240
      do 235 l = 1,neqn
        tau = h*p(l) - phi(l,15)
        p(l) = y(l) + tau
  235   phi(l,16) = (p(l) - y(l)) - tau
      go to 250
  240 do 245 l = 1,neqn
  245   p(l) = y(l) + h*p(l)
  250 xold = x
      x = x + h
      absh = abs(h)
      call f(x,p,yp)
c
c   estimate errors at orders k,k-1,k-2
c
      erkm2 = 0.0
      erkm1 = 0.0
      erk = 0.0
      do 265 l = 1,neqn
        temp3 = 1.0/wt(l)
        temp4 = yp(l) - phi(l,1)
        if(km2)265,260,255
  255   erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
  260   erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
  265   erk = erk + (temp4*temp3)**2
      if(km2)280,275,270
  270 erkm2 = absh*sig(km1)*gstr(km2)*sqrt(erkm2)
  275 erkm1 = absh*sig(k)*gstr(km1)*sqrt(erkm1)
  280 temp5 = absh*sqrt(erk)
      err = temp5*(g(k)-g(kp1))
      erk = temp5*sig(kp1)*gstr(k)
      knew = k
c
c   test if order should be lowered
c
      if(km2)299,290,285
  285 if(amax1(erkm1,erkm2) .le. erk) knew = km1
      go to 299
  290 if(erkm1 .le. 0.5*erk) knew = km1
c
c   test if step successful
c
  299 if(err .le. eps) go to 400
c       ***     end block 2     ***
c
c       ***     begin block 3     ***
c   the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
c   if third consecutive failure, set order to one.  if step fails more
c   than three times, consider an optimal step size.  double error
c   tolerance and return if estimated step size is too small for machine
c   precision.
c                   ***
c
c   restore x, phi(*,*) and psi(*)
c
      phase1 = .false.
      x = xold
      do 310 i = 1,k
        temp1 = 1.0/beta(i)
        ip1 = i+1
        do 305 l = 1,neqn
  305     phi(l,i) = temp1*(phi(l,i) - phi(l,ip1))
  310   continue
      if(k .lt. 2) go to 320
      do 315 i = 2,k
  315   psi(i-1) = psi(i) - h
c
c   on third failure, set order to one.  thereafter, use optimal step
c   size
c
  320 ifail = ifail + 1
      temp2 = 0.5
      if(ifail - 3) 335,330,325
  325 if(p5eps .lt. 0.25*erk) temp2 = sqrt(p5eps/erk)
  330 knew = 1
  335 h = temp2*h
      k = knew
      if(abs(h) .ge. fouru*abs(x)) go to 340
      crash = .true.
      h = sign(fouru*abs(x),h)
      eps = eps + eps
      return
  340 go to 100
c       ***     end block 3     ***
c
c       ***     begin block 4     ***
c   the step is successful.  correct the predicted solution, evaluate
c   the derivatives using the corrected solution and update the
c   differences.  determine best order and step size for next step.
c                   ***
  400 kold = k
      hold = h
c
c   correct and evaluate
c
      temp1 = h*g(kp1)
      if(nornd) go to 410
      do 405 l = 1,neqn
        rho = temp1*(yp(l) - phi(l,1)) - phi(l,16)
        y(l) = p(l) + rho
  405   phi(l,15) = (y(l) - p(l)) - rho
      go to 420
  410 do 415 l = 1,neqn
  415   y(l) = p(l) + temp1*(yp(l) - phi(l,1))
  420 call f(x,y,yp)
c
c   update differences for next step
c
      do 425 l = 1,neqn
        phi(l,kp1) = yp(l) - phi(l,1)
  425   phi(l,kp2) = phi(l,kp1) - phi(l,kp2)
      do 435 i = 1,k
        do 430 l = 1,neqn
  430     phi(l,i) = phi(l,i) + phi(l,kp1)
  435   continue
c
c   estimate error at order k+1 unless
c     in first phase when always raise order,
c     already decided to lower order,
c     step size not constant so estimate unreliable
c
      erkp1 = 0.0
      if(knew .eq. km1  .or.  k .eq. 12) phase1 = .false.
      if(phase1) go to 450
      if(knew .eq. km1) go to 455
      if(kp1 .gt. ns) go to 460
      do 440 l = 1,neqn
  440   erkp1 = erkp1 + (phi(l,kp2)/wt(l))**2
      erkp1 = absh*gstr(kp1)*sqrt(erkp1)
c
c   using estimated error at order k+1, determine appropriate order
c   for next step
c
      if(k .gt. 1) go to 445
      if(erkp1 .ge. 0.5*erk) go to 460
      go to 450
  445 if(erkm1 .le. amin1(erk,erkp1)) go to 455
      if(erkp1 .ge. erk  .or.  k .eq. 12) go to 460
c
c   here erkp1 .lt. erk .lt. amax1(erkm1,erkm2) else order would have
c   been lowered in block 2.  thus order is to be raised
c
c   raise order
c
  450 k = kp1
      erk = erkp1
      go to 460
c
c   lower order
c
  455 k = km1
      erk = erkm1
c
c   with new order determine appropriate step size for next step
c
  460 hnew = h + h
      if(phase1) go to 465
      if(p5eps .ge. erk*two(k+1)) go to 465
      hnew = h
      if(p5eps .ge. erk) go to 465
      temp2 = k+1
      r = (p5eps/erk)**(1.0/temp2)
      hnew = absh*amax1(0.5,amin1(0.9,r))
      hnew = sign(amax1(hnew,fouru*abs(x)),h)
  465 h = hnew
      return
c       ***     end block 4     ***
      end
      subroutine intrp(x,y,xout,yout,ypout,neqn,kold,phi,psi)
c
c   the methods in subroutine  step  approximate the solution near  x
c   by a polynomial.  subroutine  intrp  approximates the solution at
c   xout  by evaluating the polynomial there.  information defining this
c   polynomial is passed from  step  so  intrp  cannot be used alone.
c
c   this code is completely explained and documented in the text,
c   computer solution of ordinary differential equations:  the initial
c   value problem  by l. f. shampine and m. k. gordon.
c
c   input to intrp --
c
c   all floating point variables are real
c   the user provides storage in the calling program for the arrays in
c   the call list
       dimension y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),psi(12)
c   and defines
c      xout -- point at which solution is desired.
c   the remaining parameters are defined in  step  and passed to  intrp
c   from that subroutine
c
c   output from  intrp --
c
c      yout(*) -- solution at  xout
c      ypout(*) -- derivative of solution at  xout
c   the remaining parameters are returned unaltered from their input
c   values.  integration with  step  may be continued.
c
      dimension g(13),w(13),rho(13)
      data g(1)/1.0/,rho(1)/1.0/
c
      hi = xout - x
      ki = kold + 1
      kip1 = ki + 1
c
c   initialize w(*) for computing g(*)
c
      do 5 i = 1,ki
        temp1 = i
    5   w(i) = 1.0/temp1
      term = 0.0
c
c   compute g(*)
c
      do 15 j = 2,ki
        jm1 = j - 1
        psijm1 = psi(jm1)
        gamma = (hi + term)/psijm1
        eta = hi/psijm1
        limit1 = kip1 - j
        do 10 i = 1,limit1
   10     w(i) = gamma*w(i) - eta*w(i+1)
        g(j) = w(1)
        rho(j) = gamma*rho(jm1)
   15   term = psijm1
c
c   interpolate
c
      do 20 l = 1,neqn
        ypout(l) = 0.0
   20   yout(l) = 0.0
      do 30 j = 1,ki
        i = kip1 - j
        temp2 = g(i)
        temp3 = rho(i)
        do 25 l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
   25     ypout(l) = ypout(l) + temp3*phi(l,i)
   30   continue
      do 35 l = 1,neqn
   35   yout(l) = y(l) + hi*yout(l)
      return
      end
