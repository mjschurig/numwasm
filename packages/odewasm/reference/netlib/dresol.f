c
c latest revision date: august 30, 1993
c
c-----------------------------------------------------------------------
c    *** author ***
c
c author and contact 
c
c                      luca dieci (*)
c                      school of mathematics
c                      georgia institute of technology
c                      atlanta, ga 30332-0160
c                      phone: (404) 853-9209
c                      e-mail: dieci@math.gatech.edu
c
c (*) please communicate promptly
c     any mistake you should find.
c
c-----------------------------------------------------------------------
c
      subroutine dresol(neq,neqlen,x,nx,y,ny,mf,t,tout,rtol,atol,
     *       itask,istate,iopt,rwork,lrex,iwork,liw,probl,rarr,iarr)
c
c Permission to use, copy, modify, and distribute this software for any
c purpose without fee is hereby granted, provided that this entire notice
c is included in all copies of any software which is or includes a copy
c or modification of this software and in all copies of the supporting
c documentation for such software.
c THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
c WARRANTY.  IN PARTICULAR, NEITHER THE AUTHOR NOR AT&T MAKE ANY
c REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
c OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
c
      external probl
      integer neq, neqlen, nx, ny, mf, itask, istate, iopt, lrex
      integer iwork, liw, iarr
      real x, y, t, tout, rtol, atol, rwork, rarr
      dimension neq(neqlen), iwork(liw), iarr(*)
      dimension x(nx,*), y(ny), rwork(lrex), rarr(*)
c----------------------------------------------------------------------
c     *** name ****
c we will refer to the full body of the programs here as the
c
c                   ----------------------
c                   |   DRESOL PACKAGE   |
c                   ----------------------
c
c-----------------------------------------------------------------------
c    *** purpose ***
c
c dresol solves the initial value problem for stiff or nonstiff
c systems of Matrix Differential Riccati Equations
c
c     dX/dt = A21(t) + A22(t)*X -X*A11(t) - X*A12(t)*X
c     + appropriate ICs.
c
c there are a number of options provided with this solver, for
c symmetric and unsymmetric riccati equations (see documentation
c below).  We will refer to the following partition of the underlying
c matrix A(t) in what follows.
c
c               nnq     nnp     
c             |             |
c        nnq{ |A(t)  | A(t) |
c             | 11   |  12  |  
c  A(t):=     | -----|----  |  
c             |      |      |  
c        nnp{ |A(t)  | A(t) |  
c             | 21   |  22  |  
c
c as usual, the purpose of the riccati equation is that one of
c decoupling the above matrix A(t).  In particular, to bring it
c to block upper triangular structure.  That is, if X(t) solves the
c riccati equation, then A(t) is transformed in the following
c
c             |~       ~    |
c             |A(t)  | A(t) |
c  ~          | 11   |  12  |          ~     
c  A(t):=     | -----|----  |  , where A(t) = A(t)+A(t)*X(t)
c             |      | ~    |           11     11     12    
c             |  0   | A(t) |      ~                      ~ 
c             |      |  22  |      A(t) = A(t)-X(t)*A(t), A(t) = A(t)
c                                   22     22        12    12     12
c 
c in what follows, we will refer to the above partitionings.
c--------------------------------------------------------------------------
c    *** history ***
c
c the DRESOL package is a collection of subroutines for the
c numerical integration of matrix Riccati differential equations.
c the basic ODE solver on which this solver is based is the well-known
c stiff/nonstiff solver LSODE of Hindmarsh, written at
c LLL in the early '80s.  The modifications to the original
c code of Hindmarsh have been quite extensive, but 
c nonetheless much of the original philosophy and heuristics for
c order selection, local error estimatation 
c (hence stepsize change), and so forth, remain unchanged.
c 
c
c subroutines dresol calls the subroutine drequs and they both 
c constitute merely an interface between the user's driver and
c the differential riccati equation solver, subroutine drequ.
c This interface has been created to relieve the user of
c a tedious calling sequence and to keep the calling
c sequence to dresol conforming to the style of the lsode
c package of hindmarsh.
c
c
c-----------------------------------------------------------------------
c
c author: luca dieci, school of mathematics, Ga Tech
c         atlanta, ga 30332-0160
c         phone: (404) 853-9209
c         e-mail: dieci@math.gatech.edu
c
c     please address all correspondance related to the riccati
c     equation to this author, who is uniquely responsible
c     for the possible mulfunctioning of the
c     riccati solver.  please communicate promptly
c     any errors you should find.
cc
c-----------------------------------------------------------------------
c    *** references ***
c
c           Alan C. Hindmarsh, "Odepack, a systematized collection of 
c           ode solvers", in Scientific Computing, 
c           R. S. Stepleman et al. (eds.),
c           north-holland, Amsterdam, 1983, pp. 55-64.
c
c           Luca Dieci, "DRESOL. A differential riccati equations' solver:
c                        user's guide.  Version 1.0,"
c           Technical Report # 071090-52 
c           School of Mathematics, Georgia Institute of Technology
c           Atlanta, GA 30332-0160
c
c           Luca Dieci, "Numerical integration of the differential
c           Riccati equation and some related issues,"
c           SIAM J. Numer. Analysis 29-3 (1992), pp. 781-815.
c
c------------------------------------------------------------------------
c    *** summary of usage ***
c
c communication between the user and the dresol package, for normal
c situations, is handled trough calls to the subroutines dresol.
c this routine perform some of the array splittings and arrays
c checkings. (other splittings and/or checkings performed in drequ).
c if the user desires, he might call drequ directly.
c we discourage this. In any case, here below we describe basically 
c the full set of options available.  (see the description
c under the subroutine drequ for an extra avilable option on the type 
c of error tolerances allowed, if needed).
c
c this summary indicates the typical sequence of operations required
c from the user in order to properly call dresol, hence drequs and drequ
c
c a. first provide a subroutine of the form..
c               subroutine probl (t,a,rarr,iarr)
c supplying the problem matrix in a 
c
c
c b. next determine (or guess) whether or not the problem is stiff.
c stiffness occurs when the frechet derivative (jacobian) of the riccati
c equation has one of its eigenvalues with real part negative and
c large in magnitude, compared to the reciprocal of the t span of interest.
c if the problem is nonstiff, use a method flag  mf=10 or 20 (or mf=11,21)
c if it is stiff, we suggest the basic option provided by dresol
c which consists in using the backward differentiation formulas
c for carrying out one step of the integration, and perform a newton
c like iteration for the nolinear system.  when coupled with the chord
c method for the nonlinear iteration, this is the standard choice.
c this is obtained by requiring mf = 21 (and neq(5) = 0).  more choices
c are however available depending on the required strategy for solving the
c nonlinear system (the algebraic riccati equation) arising at each
c discretization step.  see the parameter neq(5) in the description below.
c
c
c c. write a main program which calls subroutine dresol once for
c each point at which answers are desired.  this should also provide
c for possible use of logical unit 6 for output of error messages.
c On the first call to drequ, supply arguments as explained in the 
c discussion which follows.  
c
c
c d. on return, the user needs to monitor the solution x (typically,
c computed at the requested tout point) and the status flag istate
c if the run has not been succesfull (as signalled by a negative value
c for istate, see below) the user might need to remove the causes for the
c insucces before calling dresol again.
c
c
c e. in case of a succesful run (positive istate), the user needs only to 
c choose a new output point (tout) and call dresol again.
c
c
c N.B.: The user should set the value of the data ncpw in the routine 
c       xerrwv accordingly to his own arithmetic environment,
c       and set the format of label 10 in subroutine xerrwv accordingly.
c       This ensures proper writing of error messages.
c       The default value at this moment is for computers with 4 bytes per
c       integer (ibm-360, ibm-370, ibm-303, ibm-43xx, ibm-pc, sun, sparc,..).
c --------------------------------------------------------------------------
c  *** interface ***
c 
c the interface to dresol has been kept basically the same as
c the original interface of lsode for conformity and safety reasons.
c so doing, we have internally generated the necessary calling sequence
c for drequ.  the user of the dresol package should not attempt
c to call drequ directly, but call subroutine dresol.
c
c the interface to dresol consists of the following parts.
c
c i.   the call sequence to subroutine dresol, which is a driver
c      routine for the solver.  this includes descriptions of the call
c      sequence arguments and of internally supplied routines.
c      following these descriptions is a description of
c      optional inputs available through the call sequence, and then
c      a description of optional outputs (in the work arrays).
c
c ii.  descriptions of other routines in the dresol package that may be
c      (optionally) called by the user.  these provide the ability to
c      alter error message handling, save and restore the internal
c      common, and obtain specified derivatives of the solution
c
c iii. descriptions of common blocks to be declared in overlay
c      or similar environments, or to be saved when doing an interrupt
c      of the problem and continued solution later.
c
c iv.  description of two routines in the dresol package, either of
c      which the user may replace with his own version, if desired.
c      these relate to the measurement of errors.
c      we recommend not to alter these routines unless absolutely
c      necessary as it might affect the performance of the method
c
c-----------------------------------------------------------------------
c **** part i.  call sequence. ***
c
c the call sequence parameters used for input only are
c     neq, neqlen, nx, ny, mf, tout, rtol, atol, itask, iopt, 
c     lrex, liw, probl,
c and those used for both input and output are
c     x, (y), t, istate.
c the work arrays rwork and iwork are also used for conditional and
c optional inputs and optional outputs.  (the term output here refers
c to the return from subroutine dresol to the user-s calling program.)
c and the arrays rarr and iarr can also be used for communication
c by the user.
c
c the legality of input parameters will be thoroughly checked on the
c initial call for the problem, but not checked thereafter unless a
c change in input parameters is flagged by istate = 3 on input.
c 
c -----------------------------------------------------------------------
c  *** explanation of parameters and basic options ***
c
c here below we provide a datailed explanation of all of the parameters in the
c calling sequence to dresol and of most of the available options.  a number
c of extra options are also available, but they should not be needed,
c and are summarized at the end of this block.
c 
c probl  = name of subroutine giving the problem matrix a.
c          it must be declared external in the calling program,
c          and it must have the following form
c
c               subroutine probl (t,a,rarr,iarr)
c               real t, rarr, a
c               integer iarr
c               dimension a(nn,nn), rarr(*), iarr(*)
c	   where nn=nnp+nnq (see neq(2), neq(3))
c
c	   in this routine we must supply al of the entries of the
c          matrix a, regardless of its structure and of icontr (see neq(4))
c          we set a by loading a(i,j) with the (i,j)-th entry of
c          the (time dependent) matrix a.  the matrix a must be
c          be supplied entirely, even if some storage savings (due to 
c          possible symmetries) had been possible.  
c          The real and integer arrays rarr and iarr are provided 
c          for possible communication between the user and the solver.
c          If the user decides to use them, then he needs to dimension 
c          them oppurtunely in the calling program.  Otherwise, just 
c          treat them as dummy arguments.  
c         (N.B.: They could be possibly used for the reimbedding strategy.)
c
c neq    = integer array of dimension at least neqlen, used for
c          communication with the solver, whose entries are as follows:
c          neq(1)=neqn, this is defined as neqn=nnp*nnq, where
c          neq(2)=nnq, this is the dimension of the square matrix a11
c          neq(3)=nnp, this is the dimension of the square matrix a22
c          neq(4)=icontr, this must be set to 0 if the riccati equation is
c                 of the type arising in optimal control from a hamiltonian
c                 matrix with guaranteed existence and stability properties.
c                 in this case, explicit exploitation of the structure
c                 is performed
c          neq(5)=nonlin, this parameter specifies how the nonlinear
c                 system arising at each step should be solved in case the
c                 user has chosen mf=11 or 21.  allowed values are
c                 nonlin=0,1 where 0 corresponds to the chord method, and
c                 1 corresponds to the quasi-newton method. a value of 1 is
c                 recommended when more accurate eigenvalue information
c                 is desired to possibly monitor stability or the like.
c                 usually, nonlin=0 is a cheaper choice (the default).
c          neq(6)=nsuper, this parameter specifies whether we want (nsuper=1)
c                 or not (nsuper=0) to perform the superstability check. 
c                 it makes sense only for mf=11,21, because the test relies
c                 on eigenvalues of the jacobian, and we further suggest 
c                 to use it only with mf=21. 
c                 This test tries to avoid that with mf=21 (or 11)
c                 the solver steps on a wrong solution branch, missing a 
c                 transition region.  The user must posses a complete 
c                 understanding of the implications of using nsuper=1, if
c                 he is going to use this value (see the ref.s).  Otherwise,
c                 we recommend to let nsuper=0 (the default) and 
c                 monitor iwork(2) (see optional output).
c                 With nsuper=1, the code attempts to avoid
c                 superstability occurrences both for neq(5)=0,1 but 
c                 it works much better for neq(5)=1, because the code has
c                 current eigenvalues information.  We suggest 
c                 to use nsuper=1 only with neq(5)=1.
c                 The superstability test performs satisfactorily
c                 when (i) there is a stable trajectory for the solution 
c                 and the solution is (initially) along this trajectory 
c                 (ii) the solution is along an unstable trajectory and we 
c                 want to be cautious.  Again, we refer to the papers given
c                 in the introduction for a more detailed explanation.
c           N.B.: the solver tries to detect/avoid the superstability 
c                 occurrences (see the documentation in the subroutine stode)
c                 on a given step, up to 6 step reductions, and unless the
c                 stepsize has become so small to lead to a loss of meaning.
c                 if the code fails to recover (due to these occurrences)
c                 from the possible superstability, then it does not try
c                 any more to recover for the given run and it signals its
c                 failing by setting nsuper=0, so that the user should 
c                 monitor this quantity if he originally set nsuper=1.
c                 (the code prints a message saying that this occurred).
c
c neqlen = dimension of neq.  it must have neqlen.ge.6.  
c
c x      = riccati matrix of dimension (nx,nnq) containing the initial
c          values of the riccati solution.  user must supply all entries of
c          x even if icontr=0, and hence there could have been symmetries in x
c          (which will be exploited internally by the code, in this case)
c          on output, x contains the matrix of computed values x(t).
c
c nx     = row dimension of the matrix x, must have nx.ge.nnp
c
c y      = work array of length ny used internally by the code, whose 
c          content is usually of no relevance to the user.  the first neqn 
c          components of y contain the riccati matrix x unrolled
c          according to y(k)=x(i,j), k=nnp*(j-1)+i, i=1,nnp, j=1,nnq.
c          the array y is passed as the y argument in all calls.
c          y(neqn+1),... are used to store the problem matrix 
c          and its partitioned form as needed.
c
c ny	 = dimension of y of length at least neqn+2*nn*nn
c
c t      = on input, the initial value of the independent variable.
c          on output, value of independent variable (normally tout),
c          where the solution has been evaluated.  on an error return, 
c          t is simply the farthest point reached.
c
c tout   = first point where output is desired, used only for input.
c          notice that we allow for either t.gt.tout or t.lt.tout
c          and so both forward and backward integration are allowed.
c          when starting the problem (istate = 1), tout may be equal
c          to t for one call, then should .ne. t for the next call.
c          for the initial t, an input value of tout .ne. t is used
c          in order to determine the direction of the integration
c          (i.e. the algebraic sign of the step sizes) and the rough
c          scale of the problem. 
c          if itask = 2 or 5 (one-step modes), tout is ignored after
c          the first call (i.e. the first call with tout .ne. t).
c          otherwise, tout is required on every call.
c          if itask = 1, 3, or 4, the values of tout need not be
c          monotone, but a value of tout which backs up is limited
c          to the current internal t interval, whose endpoints are
c          tcur - hu and tcur (see optional outputs, below, for
c          tcur and hu).
c
c rtol   = relative tolerance parameter (scalar).
c
c atol   = absolute tolerance parameter (scalar).
c          the estimated local error in y(i) will be controlled so as
c          to be roughly less (in magnitude) than
c             ewt(i) = rtol*abs(y(i)) + atol
c          thus the local error test passes if, in each component,
c          either the absolute error is less than atol,
c          or the relative error is less than rtol.
c          use rtol = 0.0 for pure absolute error control, and
c          use atol = 0.0 for pure relative error control (not to use if
c          some solution component is known to vanish).
c          the values of rtol and atol should all be non-negative.
c          if the above choice (with atol, rtol scalar) is not suitable, 
c          more general choices can be obtained.  we refer to the
c          supplementary documentation at the beginning of the
c          subroutine drequ below.
c          caution.. actual (global) errors may exceed these
c          local tolerances, so choose them conservatively.
c          if global errors are to be estimated by making a repeated
c          run on the same problem with smaller tolerances, then all
c          components of rtol and atol (i.e. of ewt) should be scaled
c          down uniformly.
c 
c itask  = an index specifying the task to be performed.
c          input only.  itask has the following values and meanings.
c          1  means normal computation of output values of at
c             t = tout (by overshooting and interpolating).
c          2  means take one step only and return.
c          3  means stop at the first internal mesh point at or
c             beyond t = tout and return.
c          4  means normal computation of output values at
c             t = tout but without overshooting t = tcrit.
c             tcrit must be input as rwork(1).  tcrit may be equal to
c             or beyond tout, but not behind it in the direction of
c             integration.  this option is useful if the problem
c             has a singularity at or beyond t = tcrit.
c          5  means take one step, without passing tcrit, and return.
c             tcrit must be input as rwork(1).
c          note..  if itask = 4 or 5 and the solver reaches tcrit
c          (within roundoff), it will return t = tcrit (exactly) to
c          indicate this (unless itask = 4 and tout comes before tcrit,
c          in which case answers at t = tout are returned first).
c
c istate = an index used for input and output to specify the
c          state of the calculation.
c          on input, the values of istate are as follows.
c          1  means this is the first call for the problem
c             (initializations will be done).  see note below.
c          2  means this is not the first call, and the calculation
c             is to continue normally, with no change in any input
c             parameters except possibly tout and itask.
c             (if itol, rtol, and/or atol are changed between calls
c             with istate = 2, the new values will be used but not
c             tested for legality.)
c          3  means this is not the first call, and the
c             calculation is to continue normally, but with
c             a change in input parameters other than
c             tout and itask.  changes are allowed in
c             neq(4), neq(5), neq(6), 
c             rtol, atol, iopt, lrex, liw, mf, 
c             and any of the optional inputs except h0.
c             (of course, above changes can be made subject to mutual
c             compatibilities). 
c          note..  a preliminary call with tout = t is not counted
c          as a first call here, as no initialization or checking of
c          input is done.  (such a call is sometimes useful for the
c          purpose of outputting the initial conditions.)
c          thus the first call for which tout .ne. t requires
c          istate = 1 on input.
c          on output, istate has the following values and meanings.
c           1  means nothing was done, as tout was equal to t with
c              istate = 1 on input.  (however, an internal counter was
c              set to detect and prevent repeated calls of this type.)
c           2  means the integration was performed successfully.
c          -1  means an excessive amount of work (more than mxstep
c              steps) was done on this call, before completing the
c              requested task, but the integration was otherwise
c              successful as far as t.  (mxstep is an optional input
c              and is normally 500.)  to continue, the user may
c              simply reset istate to a value .gt. 1 and call again
c              (the excess work step counter will be reset to 0).
c              in addition, the user may increase mxstep to avoid
c              this error return (see below on optional inputs).
c          -2  means too much accuracy was requested for the precision
c              of the machine being used.  this was detected before
c              completing the requested task, but the integration
c              was successful as far as t.  to continue, the tolerance
c              parameters must be reset, and istate must be set
c              to 3.  the optional output tolsf may be used for this
c              purpose.  (note.. if this condition is detected before
c              taking any steps, then an illegal input return
c              (istate = -3) occurs instead.)
c          -3  means illegal input was detected, before taking any
c              integration steps.  see written message for details.
c              note..  if the solver detects an infinite loop of calls
c              to the solver with illegal input, it will cause
c              the run to stop.
c          -4  means there were repeated error test failures on
c              one attempted step, before completing the requested
c              task, but the integration was successful as far as t.
c              the problem may have a singularity, or the input
c              may be inappropriate.
c          -5  means there were repeated convergence test failures on
c              one attempted step, before completing the requested
c              task, but the integration was successful as far as t.
c              this may be caused by an inaccurate jacobian matrix,
c              if one is being used.
c          -6  means ewt(i) became zero for some i during the
c              integration.  pure relative error control (atol=0.0)
c              was requested, but a variable has now vanished.
c              the integration was successful as far as t.
c          note..  since the normal output value of istate is 2,
c          it does not need to be reset for normal continuation.
c          also, since a negative input value of istate will be
c          regarded as illegal, a negative output value requires the
c          user to change it, and possibly other inputs, before
c          calling the solver again.
c
c iopt   = an integer flag to specify whether or not any optional
c          inputs are being used on this call.  input only.
c          the optional inputs are listed separately below.
c          iopt = 0 (the default) means no optional inputs are being used.
c                   default values will be used in all cases.
c          iopt = 1 means one or more optional inputs are being used.
c
c rwork  = a real work array of length at least lrex
c          the content of rwork can be very relevant to the user,
c          and some of it is explained in the optional output later on.
c          the first 20 words of rwork are reserved for conditional
c          and optional inputs and optional outputs.
c
c lrex   = declared length of rwork; we must have lrex at least
c             20 + neqn*(maxord + 4) + lns    where
c          neqn = neq(1) 
c          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
c                   smaller value is given as an optional input),
c          lns = 0                           if miter = 0,
c          lns = 2+2*(nnp**2+nnq**2)+3*nn    if miter = 1
c          regardless of nonlin=0,1.
c          (see the mf description for meth and miter.)
c          this paramter lrex will be checked by the solver.
c          thus if maxord has its default value this length is..
c             20 + 16*neqn                  for mf = 10,
c             20 + 16*neqn + lns            for mf = 11 
c             20 +  9*neqn                  for mf = 20,
c             20 +  9*neqn + lns            for mf = 21 
c note..  the work arrays must not be altered between calls to drequ
c for the same problem.
c
c iwork  = integer work array of length at least liw
c          the first few words of iwork are used for conditional
c          and optional input and optional output
c
c liw    = declared length of iwork; we must have liw at least
c             20        for mf = 10, 20
c             20 + nn   for mf = 11, 21
c          this parameter liw will be checked by the solver
c
c mf     = the method flag.  used only for input.  the legal values of
c          mf are 10, 11, 20, 21.
c          mf has decimal digits meth and miter.. mf = 10*meth + miter.
c          meth indicates the basic linear multistep method..
c            meth = 1 means the implicit adams method.
c            meth = 2 means the method based on backward
c                     differentiation formulas (bdf-s).
c          miter indicates the corrector iteration method..
c            miter = 0 means functional iteration (no jacobian matrix
c                      is involved).
c            miter = 1 means chord iteration with an internally generated 
c                      jacobian, kept in compact matrix form.
c          if miter = 1, the user must also specify the particular 
c          choice of nonlinear solver he desires, as explained under
c          the parameter neq(5).
c          To summarize, the following are available and suggested options
c          10 for nonstiff (adams) method, no jacobian used, simple
c             functional fixed point iteration for nonlinear system
c          11 for adams method, but with jacobian used in a newton type
c             iteration (see neq(5))
c          20 for nonstiff (bdf) method, functional fixed point iteration
c          21 for stiff (bdf) method, but with jacobian used (recommended
c             value for bdf method).
c          We recommend to use mf=10 or 20 for nonstiff Riccati 
c          equations and mf=21 for stiif equations.
c
c rarr,  = real and integer arrays which can be used by the user for
c   iarr     communication between the calling program and routines 
c          probl (hence all of those called by probl).  dresol does not
c          touch these arrays, but it passes them.  If the user does
c          not need them ,just treat them as dummy arguments.  if user
c          decides to use them, they must be declared as arrays of
c          appropriate length in the calling program and in probl.
c          they can be put to use if user wants to perform a reimbedding
c          strategy, in the usual ways for riccati equations.
c
c note that the main program must declare arrays neq, x, y, rwork, iwork, and
c rarr and iarr if they are used.
c
c-----------------------------------------------------------------------
c  *** optional inputs ***
c
c the following is a list of the optional inputs provided for in the
c call sequence.  (see also part ii.)  for each such input variable,
c this table lists its name as used in this documentation, its
c location in the call sequence, its meaning, and the default value.
c the use of any of these inputs requires iopt = 1, and in that
c case all of these inputs are examined.  a value of zero for any
c of these optional inputs will cause the default value to be used.
c thus to use a subset of the optional inputs, simply preload
c locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
c then set those of interest to nonzero values.
c
c name    location      meaning and default value
c
c tcrit   rwork(1)  critical value of t which the solver
c                   is not to overshoot.  required if itask is
c                   4 or 5, and ignored otherwise.  (see itask.)
c
c h0      rwork(5)  the step size to be attempted on the first step.
c                   the default value is determined by the solver.
c
c hmax    rwork(6)  the maximum absolute step size allowed.
c                   the default value is infinite.
c
c hmin    rwork(7)  the minimum absolute step size allowed.
c               note: the default value is 0.0, so be careful.     
c                   (this lower bound is not enforced on the final
c                   step before reaching tcrit when itask = 4 or 5.)
c
c maxord  iwork(5)  the maximum order to be allowed.  the default
c                   value is 12 if meth = 1, and 5 if meth = 2.
c                   if maxord exceeds the default value, it will
c                   be reduced to the default value.
c                   if maxord is changed during the problem, it may
c                   cause the current order to be reduced.
c
c mxstep  iwork(6)  maximum number of (internally defined) steps
c                   allowed during one call to the solver.
c                   the default value is 500.
c
c mxhnil  iwork(7)  maximum number of messages printed (per problem)
c                   warning that t + h = t on a step (h = step size).
c                   this must be positive to result in a non-default
c                   value.  the default value is 10.
c-----------------------------------------------------------------------
c *** optional outputs ***
c
c as optional additional output from dresol, the variables listed
c below are quantities related to the performance of dresol
c which are available to the user.  these are communicated by way of
c the work arrays, but also have internal mnemonic names as shown.
c except where stated otherwise, all of these outputs are defined
c on any successful return from drequ, and on any return with
c istate = -1, -2, -4, -5, or -6.  on an illegal input return
c (istate = -3), they will be unchanged from their existing values
c (if any), except possibly for tolsf, lenrw, and leniw.
c on any error return, outputs relevant to the error will be defined,
c as noted below.
c
c name    location      meaning
c
c hu      rwork(11) the step size in t last used (successfully).
c
c hcur    rwork(12) the step size to be attempted on the next step.
c
c tcur    rwork(13) the current value of the independent variable
c                   which the solver has actually reached, i.e. the
c                   current internal mesh point in t.  on output, tcur
c                   will always be at least as far as the argument
c                   t, but may be farther (if interpolation was done).
c
c tolsf   rwork(14) a tolerance scale factor, greater than 1.0,
c                   computed when a request for too much accuracy was
c                   detected (istate = -3 if detected at the start of
c                   the problem, istate = -2 otherwise).  
c                   if rtol and atol are uniformly
c                   scaled up by a factor of tolsf for the next call,
c                   then the solver is deemed likely to succeed.
c                   (the user may also ignore tolsf and alter the
c                   tolerance parameters in any other way appropriate.)
c
c iwm(1)   iwork(1)  if miter=1 (mf=11 or 21), this entry is set to
c                    1 if during integration the solver found purely 
c                    imaginary eigenvalues.  it is set to 0 otherwise.
c                    of course, this is not an error condition, but it 
c                    might indicate that the user might want to use 
c                    the adams formulas if he is using the backward-
c                    differentiation formulas.
c
c iwm(2)   iwork(2)  if miter=1 (mf=11 or 21), this entry os set to
c                    1 if during integration the solver found unstable
c                    eigenvalues with respect to the direction of 
c                    integration.  it is set to 0 otherwise.
c                    this is not an error condition, but the user should    
c                    monitor the eigenvalues to be certain that the
c                    integrator is performing adequately.
c
c nst     iwork(11) the number of steps taken for the problem so far.
c
c nfe     iwork(12) the number of f evaluations for the problem so far.
c
c nje     iwork(13) the number of jacobian evaluations so far.
c                   in this code, every jacobian evaluation corresponds 
c                   to the solution of a sylvester equation, hence
c                   to a factorization of the blocks 
c                      ~    ~    
c                      A  , A   via the bartels-stewart algorithm
c                       11   22
c                   (see the documentation for prepj).
c
c nqu     iwork(14) the method order last used (successfully).
c
c nqcur   iwork(15) the order to be attempted on the next step.
c
c imxer   iwork(16) the index of the component of largest magnitude in
c                   the weighted local error vector ( e(i)/ewt(i) ),
c                   on an error return with istate = -4 or -5.
c
c lenrw   iwork(17) the length of rwork actually required.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c leniw   iwork(18) the length of iwork actually required.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c
c the following arrays are segments of the rwork array which
c can be of interest to the user as optional outputs.
c for each array, the table below gives its internal name,
c its base address in rwork, and its description.
c 
c
c N.B.: Important eigenvalue information is in these arrays in the entries
c       era, erb, eia, aib, itypea, itypeb, as described below.
c
c name    base address      description
c
c in the rwork array
c
c yh      21             the nordsieck history array, of size neqn by
c                        (nqcur + 1). for j = 0,1,...,nqcur, column j+1
c                        of yh contains hcur**j/factorial(j) times
c                        the j-th derivative of the interpolating
c                        polynomial currently representing the solution,
c                        evaluated at t = tcur.
c
c era      lwm+nttp      only if miter=1 (mf=11,21), and we have set
c                        lwm=20+(maxord+1)*neqn, nttp=3+2*nnp**2+2*nnq**2+nnp.
c                        array of dimension nnp, of real parts of ev.s
c                    --> if icontr=0, it contains the real parts of the
c                        eigenvalues of the closed loop spectrum matrix
c                        ~  ^            ^  
c                        A (t   ), where t    is the farthest point 
c                         11 n+1          n+1
c                        reached in the integration (not necessarily
c                        the same as tout), when nonlin=1, and otherwise  
c                        is the last point where the code performed a 
c                        jacobian evaluation/factorization if nonlin=0.
c                        These eigenvalues are ordered from smallest
c                        to largest real parts
c                    --> if icontr=1, it contains the real parts of the 
c                        eigenvalues of the (nnp*nnp) lower matrix
c                        ~  ^            ^  
c                        A (t   ), where t    is as above. 
c                         22 n+1          n+1
c                        The eigenvalues are ordered from smallest to 
c                        largest real part.
c
c eia      lwm+ntt2p     only when miter=1, and we set ntt2p=nttp+nnp.
c                        array of dimension nnp, containing the imaginary
c                        parts of the eigenvalues whose real parta are era
c                        (see aabove entry), and they are stored as in the
c                        way explained in era, that is eia(i) is imaginary
c                        part corresponding to era(i).
c
c erb      lwm+ntt1q     only when miter=1 and icontr=1, and we have
c                        set ntt1q=lwm+nttp+2*nnp+nnq.
c                        array of dimension nnq, containing the real
c                        parts of the eigenvalues of the upper matrix
c                        ~  ^            ^  
c                        A (t   ), where t    is as above. 
c                         11 n+1          n+1
c                        The eigenvalues are ordered from smallest to 
c                        largest real part.
c
c eib       lwm+ntt2q    only when miter=1 and icontr=1, where
c                        we set ntt2q=ntt1q+nnq.
c                        array of dimension nnq, containing the imaginary
c                        parts eib(i), of the eigenvalues whose real parts 
c                        are given by erb(i).
c
c acor     lenrw-neq+1   array of size neqn used for the accumulated
c                        corrections on each step, scaled on output
c                        to represent the estimated local error in x (i.e.
c                        in y) on the last step.  this is the vector e in
c                        the description of the error control.  it is
c                        defined only on a successful return from dresol.
c
c in the iwork array
c
c itypea   21            array of size nnp, whose entries specify whether
c                        the i-th eigenvalue corresponding to (era(i),eia(i))
c                        is real, purely imaginary, or complex with nonzero
c                        real part.  see the documnetation to the
c                        subroutines axpxb, atxpxa, hqr3.
c
c itypeb   21+nnp        array of size nnq, whose entries are as for
c                        itypea, but now they are relative to the
c                        eigenvalues (erb(i),eib(i)).
c
c N.B.: further information which might be of some use is 
c       given by the orthogonal matrices reducing the sylvester equation
c       blocks, when miter=1.  if interested, see description in prepj
c       and connected routines.
c
c ---------------------------------------------------------------------
c *** part ii.  other routines callable ***
c
c the following are optional calls which the user may make to
c gain additional capabilities in conjunction with dresol.
c (the routines xsetun and xsetf are designed to conform to the
c slatec error handling package.)
c
c     form of call                  function
c   call xsetun(lun)          set the logical unit number, lun, for
c                             output of messages from dresol, if
c                             the default is not desired.
c                             the default value of lun is 6.
c
c   call xsetf(mflag)         set a flag to control the printing of
c                             messages by dresol.
c                             mflag = 0 means do not print. (danger..
c                             this risks losing valuable information.)
c                             mflag = 1 means print (the default).
c
c                             either of the above calls may be made at
c                             any time and will take effect immediately.
c
c   call srcom(rsav,isav,job) saves and restores the contents of
c                             the internal common blocks used by
c                             dresol (see part iii below).
c                             rsav must be a real array of length 222
c                             or more, and isav must be an integer
c                             array of length 66 or more.
c                             job=1 means save common into rsav/isav.
c                             job=2 means restore common from rsav/isav.
c                                srcom is useful if one is
c                             interrupting a run and restarting
c                             later, or alternating between two or
c                             more problems solved with dresol.
c
c   call intdy(,,,,,)         provide derivatives of x (i.e. of y), of 
c        (see below)          various orders, at a specified point t, if
c                             desired.  it may be called only after
c                             a successful return from dresol.
c
c the detailed instructions for using intdy are as follows.
c the form of the call is..
c
c   call intdy (t, k, rwork(21), neqn, dky, iflag)
c
c the input parameters are..
c
c t         = value of independent variable where answers are desired
c             (normally the same as the t last returned by drequ).
c             for valid results, t must lie between tcur - hu and tcur.
c             (see optional outputs for tcur and hu.)
c k         = integer order of the derivative desired.  k must satisfy
c             0 .le. k .le. nqcur, where nqcur is the current order
c             (see optional outputs).  the capability corresponding
c             to k = 0, i.e. computing y(t), is already provided
c             by drequ directly.  since nqcur .ge. 1, the first
c             derivative dx/dt is always available with intdy.
c rwork(21) = the base address of the history array yh.
c neqn      = column length of yh, equal to the value of neq(1).
c
c the output parameters are..
c
c dky       = a real array of length neqn containing the computed value
c             of the k-th derivative of the riccati solution x(t).
c iflag     = integer flag, returned as 0 if k and t were legal,
c             -1 if k was illegal, and -2 if t was illegal.
c             on an error return, a message is also written.
c-----------------------------------------------------------------------
c *** part iii.  common blocks ***
c
c if dresol is to be used in an overlay situation, the user
c must declare, in the primary overlay, the variables in..
c   (1) the call sequence to dresol,
c   (2) the three internal common blocks
c         /ls0001/  of length  257  (218 single precision words
c                         followed by 39 integer words),
c         /eh0001/  of length  2 (integer words),
c         /intdre/ of length 29 (4 single precision words
c                         followed by 25 integer words).
c
c if dresol is used on a system in which the contents of internal
c common blocks are not preserved between calls, the user should
c declare the above three common blocks in his main program to insure
c that their contents are preserved.
c
c if the solution of a given problem by dresol is to be interrupted
c and then later continued, such as when restarting an interrupted run
c or alternating between two or more problems, the user should save,
c following the return from the last dresol call prior to the
c interruption, the contents of the call sequence variables and the
c internal common blocks, and later restore these values before the
c next dresol call for that problem.  to save and restore the common
c blocks, use subroutine srcom (see part ii above).
c
c-----------------------------------------------------------------------
c *** part iv.  optionally replaceable solver routines ***
c
c below are descriptions of two routines in the dresol package which
c relate to the measurement of errors.  either routine can be
c replaced by a user-supplied version, if desired.  however, since such
c a replacement may have a major impact on performance, it should be
c done only when absolutely necessary, and only with great caution.
c (note.. the means by which the package version of a routine is
c superseded by the user-s version may be system-dependent.)
c
c (a) ewset.
c the following subroutine is called just before each internal
c integration step, and sets the array of error weights, ewt, as
c described under rtol atol above.. (see documentation of drequ if
c need further atol/rtol possibilities, such as vector tolerances)
c     subroutine ewset (neqn, itol, rtol, atol, ycur, ewt)
c where neqn, rtol, and atol are as in the dresol call sequence,
c itol = 1 in this version (unless user decides to call drequ directly,
c and set itol to a different value, but this is discouraged),
c ycur contains the current dependent variable vector, and
c ewt is the array of weights set by ewset.
c
c if the user supplies this subroutine, it must return in ewt(i)
c (i = 1,...,neqn) a positive quantity suitable for comparing errors
c in y(i) to (recall that the riccati matrix x is stored internally
c in the vector y, columnwise).  
c the ewt array returned by ewset is passed to the vnorm 
c routine (see below), and also used by dresol in the computation
c of the optional output imxer.
c
c in the user-supplied version of ewset, it may be desirable to use
c the current values of derivatives of the riccati solution x.
c Derivatives up to order nq are available from the history array yh, 
c described above under optional outputs.  
c In ewset, yh is identical to the ycur array,
c extended to nq + 1 columns with a column length of neqn and scale
c factors of h**j/factorial(j).  on the first call for the problem,
c given by nst = 0, nq is 1 and h is temporarily set to 1.0.
c the quantities nq, neqn, h, and nst can be obtained by including
c in ewset the statements..
c     common /ls0001/ rls(218),ils(39)
c     nq = ils(35)
c     neqn = ils(14)
c     nst = ils(36)
c     h = rls(212)
c thus, for example, the current value of dy/dt can be obtained as
c ycur(neqn+i)/h  (i=1,...,neq)  (and the division by h is
c unnecessary when nst = 0).
c
c (b) vnorm.
c the following is a real function routine which computes the weighted
c max-norm norm of a vector v..
c     d = vnorm (n, v, w)
c where..
c   n = the length of the vector,
c   v = real array of length n containing the vector,
c   w = real array of length n containing weights,
c   d = max ( abs(v(i)*w(i)), i=1, n )
c vnorm is called with n = neqn and with w(i) = 1.0/ewt(i), where
c ewt is as set by subroutine ewset.
c
c if the user supplies this function, it should return a non-negative
c value of vnorm suitable for use in the error control in drequ.
c none of the arguments should be altered by vnorm.
c for example, a user-supplied vnorm routine might..
c   -substitute the root-mean-square norm for the max-norm, or
c   -ignore some components of v in the norm, with the effect of
c    suppressing the error control on those components of y.
c-----------------------------------------------------------------------
c *** other routines in the dresol package ***
c
c in addition to subroutine dresol, drequs and drequ, which effectively 
c perform initializations, decide upon returns-errors, and so forth,
c in such a way to keep the code as close to lsode as possible,
c the dresol package also includes the following subroutines and 
c function routines.. (some of these are straight out of the lsode 
c package, other are modified from lsode, and other are original to dresol).
c  intdy    computes an interpolated value of the y vector (which is
c           the unrolled riccati matrix x) at t = tout.
c  stode    is the core integrator, which does one step of the
c           integration and the associated error control
c           and eigenvalues checking.  if nsuper=1 also the superstability
c           testing is performed here.
c  cfode    sets all method coefficients and test constants.
c  f, doblks, and foftex compute the right hand side of the riccati
c           equation by maximizing efficiency (also depending on icontr)
c  prepj, jac  manage computation and preprocess of the jacobian matrix and  
c           the newton iteration matrix p = i - h*l0*j.  there are several
c           important technical details here.  these routines call 
c           the routines axpxbd, atxxad, orthes, ortran, hqr3, exchng,
c           qrstep, split, to perform an ordered eigenvalue decomposition.
c  solsy    manages solution of linear system in chord iteration and in the
c           quasi-newton iteration (nonlin=0,1).  it calls the routines
c           axpxbs, shrslv, atxxas, symslv, sysslv, to solve the decomposed
c           sylvester equations.
c  ewset    sets the error weight vector ewt before each step.
c  vnorm    computes the weighted max norm of a vector.
c  srcom    is a user-callable routine to save and restore
c           the contents of the internal common blocks.
c  r1mach   computes the unit roundoff in a machine-independent manner.
c  xerrwv, xsetun, and xsetf   handle the printing of all error
c           messages and warnings.  xerrwv is machine-dependent.
c note..  vnorm and r1mach are function routines.
c all the others are subroutines.
c
c there are also two block data in dresol, atlan1 and atlan2 which set
c the constant illin, ntrep (in atlan1) and msflg, lunit (in atlan2)
c (see part ii. other routines callable, for replacing these entries).
c
c the intrinsic and external routines used by dresol are..
c abs, amax1, amin1, float, max0, min0, mod, sign, sqrt, and write.
c *** end of prologue ***
      integer neqn, ney
      integer nnq, nnp, icontr, nonlin, nsuper, nn, nns, nqs, nps,
     1   neqn1, neqn11, neqn12, neqn21, neqn22, np3, np23, np23q,
     2   ntt, nttp, ntt2p, ntt1, ntt1q, ntt2q, nflevc, nflevp
      real reuvct, reuvpt, tcj, tpj
c-----------------------------------------------------------------------
c the following internal common block contains
c (a) variables which are communicated between subroutines.  these are
c variables needed for arraies spitting and other checkings.
c the structure of the block is as follows..  all real variables are
c listed first, followed by all integers.  the block is declared in 
c subroutines dresol, drequs, drequ, stode, f, prepj, and solsy.  
c most of the variables of this block are of no interest to the
c user, but he might find interesting the quantity tcj and tpj:
c they are the last (and the previous to the last) time at which
c a jacobian was decomposed/updated.
c-----------------------------------------------------------------------
      common /intdre/ reuvct, reuvpt, tcj, tpj, nnq, nnp, icontr, 
     1   nonlin, nsuper, nn, nns, nqs, nps, neqn1, neqn11, neqn12, 
     2   neqn21, neqn22, np3, np23, np23q, ntt, nttp, ntt2p, ntt1, 
     3   ntt1q, ntt2q, nflevc, nflevp
c here we perform an array splitting and check for the first time
c some of the arrays dimensions
c (further checking will be done in drequ)
c the checking is done only for first call or if input has been changed
      if (istate.ne.1.and.istate.ne.3) goto 10
      neqn=neq(1)
      nnq=neq(2)
      nnp=neq(3)
      icontr=neq(4)
      if (icontr .ne. 0) icontr = 1
      if (icontr.eq.0.and.nnq.ne.nnp) then
         istate=-3
         call xerrwv(45hdresol--  wrong icontr-nnp-nnq: run aborted  ,
     1      45,1,2,0,0,0,0,0.e0,0.e0)
         return
      endif
      nonlin=neq(5)
      if (nonlin.ne.1) nonlin=0
      nsuper=neq(6)
      if (nsuper.ne.1) nsuper=0
      nn=nnp+nnq
      if (neqlen.lt.6) then
         istate=-3
         call xerrwv(40hdresol-- neqlen too small: run aborted  ,
     1      40,1,2,0,0,0,0,0.e0,0.e0)
         return
      endif
      nns=nn*nn
      ney=neqn+2*nns
      if (ny.lt.ney) then
         istate=-3
         call xerrwv(40hdresol-- ny too small: run aborted      ,
     1      40,1,2,0,0,0,0,0.e0,0.e0)
         return
      endif
      nqs=nnq*nnq
      nps=nnp*nnp
      iwork(1)=0
      iwork(2)=0
      nflevc=1
      nflevp=1
      neqn1=neqn+1
      neqn11=neqn1+nns
      neqn12=neqn11+nqs
      neqn21=neqn12+neqn
      neqn22=neqn21+neqn
      np3=nps+3
      np23=np3+nps
      np23q=np23+nqs
      ntt=np23q+nqs
      nttp=ntt+nnp
      ntt2p=nttp+nnp
      ntt1=ntt2p+nnp
      ntt1q=ntt1+nnq
      ntt2q=ntt1q+nnq
c splitting of the various arrays completed
c ready to call the solver
 10   call drequs(neq,x,nx,y,mf,t,tout,rtol,atol,itask,istate,iopt,
     *            rwork,lrex,iwork,liw,probl,rarr,iarr)
      return
c ----------------- end of subroutine dresol -----------------------------
      end
      subroutine drequs(neq,x,nx,y,mf,t,tout,
     *                  rtol,atol,itask,istate,iopt,
     *                  rwork,lrw,iwork,liw,probl,rarr,iarr)
c
c revision date: July 10 1990
c author: Luca Dieci, School of Mathematics, GaTech
c
      external probl,f,jac
      integer mf, itask, istate, iopt, lrw, iwork, liw, iarr
      real x, y, t, tout, rtol, atol, rwork, rarr
      dimension neq(*),x(nx,*),y(*),rwork(lrw), 
     1   iwork(liw), rarr(*), iarr(*)
c subroutine drequs completes the interface to the solver drequ
c it does further array splitting and sets the input to conform
c to the LSODE solver of hindmarsh.
      integer i, j, inpj
      integer nnq, nnp, icontr, nonlin, nsuper, nn, nns, nqs, nps,
     1   neqn1, neqn11, neqn12, neqn21, neqn22, np3, np23, np23q,
     2   ntt, nttp, ntt2p, ntt1, ntt1q, ntt2q, nflevc, nflevp
      real reuvct, reuvpt, tcj, tpj
      common /intdre/ reuvct, reuvpt, tcj, tpj, nnq, nnp, icontr, 
     1   nonlin, nsuper, nn, nns, nqs, nps, neqn1, neqn11, neqn12, 
     2   neqn21, neqn22, np3, np23, np23q, ntt, nttp, ntt2p, ntt1, 
     3   ntt1q, ntt2q, nflevc, nflevp
      if (istate.ne.1.and.istate.ne.3) goto 10
      do 5 j=1,nnq
         do 5 i=1,nnp
            inpj=i+nnp*(j-1)
            y(inpj)=x(i,j)
  5   continue
c force scalar tolerances to drequ (the revised lsode)
 10   itol=1
      call drequ(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,
     *           rwork,lrw,iwork,liw,jac,mf,probl,rarr,iarr)
      do 14 j=1,nnq
         do 14 i=1,nnp
            inpj=i+nnp*(j-1)
 14         x(i,j)=y(inpj)
      return
c ---------- end of subroutine drequs -------------------------------
      end
      subroutine drequ (f, neq, y, t, tout, itol, rtol, atol, itask,
     1                  istate, iopt, rwork, lrw, iwork, liw, jac, mf,
     2                  probl, rarr, iarr)
c
c revision date: July 10 1990 
c author of modifications: Luca Dieci, School of Mathematics, GaTech
c
c history: this is essentially "subroutine lsode" written by Alan
c          Hindmarsh at LLL, with the necessary modifications for
c          differential Riccati equations
c
      external f, jac, probl
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf, iarr
      real y, t, tout, rtol, atol, rwork
      dimension neq(*), y(*), rtol(*), atol(*), rwork(lrw), iwork(liw)
      dimension rarr(*), iarr(*)
c ***********************************************************
c       This program drequ constitutes a modification of the
c       lsode package made by Luca Dieci at georgia tech 
c       for solving Differential Riccati Equations.
c       The modified code mantains many of the features
c       of the lsode package, but has a totally redesigned
c       nonlinear system solver and error control
c       (being based also upon the eigenvalues of the Jacobian)
c ************************************************************
c
c drequ solves the initial value problem for stiff or nonstiff
c systems of first order Riccati Differential Equations
c
c     dX/dt = A21(t) + A22(t)*X -X*A11(t) - X*A12(t)*X
c     + appropriate ICs.
c
c      The program actually rewrites (because of the
c      column orientation of fortran, but matrix arithmetic
c      is performed throughout) the riccati equation as
c
c     dy/dt = f(t,y) ,  or, in component form,
c     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neqn)) (i = 1,...,neqn).
c
c the user, however, should never call drequ directly, but
c rather use the calling sequence as described in drequ
c the next information should be read only to use the higher capabilities
c of the dresol package, but the extra available options should not
c be necessary.
c
c-----------------------------------------------------------------------
c author and contact..
c                      (original author of the LSODE package)
c                      alan c. hindmarsh,
c                      computing and mathematics research div., l-316
c                      lawrence livermore national laboratory
c                      livermore, ca 94550.
c
c                      (author of the solver for riccati equations)
c                      luca dieci (*)
c                      school of mathematics
c                      georgia institute of technology
c                      atlanta, ga 30332
c
c (*) author to whom correspondance related to the riccati
c     equation should be addressed, and unique responsible
c     for the possible mulfunctioning of this
c     riccati solver.  please communicate promptly
c     any mistake you should find.
c
c ----------------------------------------------------------
c the descriptions of the arguments which are not examined in
c the introduction to dresol follows. we refer to the introduction
c to dresol for a full explanation of the other parameters.
c
c f      = the name of the subroutine defining the riccati ode system.
c          the system has been internally put in the first-order
c          form dy/dt = f(t,y), where f is a vector-valued function
c          of the scalar t and the vector y through the related
c          subroutines doblks and foftex (user might want to
c          look at these routines if he desires to call drequ 
c          directly).  in any case, the purpose is to provide dresol
c          with a way to compute the function f.
c	   the subroutine f has the form
c		subroutine f(t, y, ydot, probl, rarr, iarr)
c               dimension y(*), ydot(*), neq(*), rarr(*), iarr(*)
c	   where t, probl and y are input, and ydot = f(t,y) is output.
c	   for the dimensions of y needed, see later (ydot does not
c	   need to be explicitely dimensioned by the calling program as
c	   it is always made up by pieces of rwork)
c          subroutine f does not alter y(1),...,y(neq(1)).
c	   but it alters the other entries of y.
c	   f is declared external in the calling program (i.e. drequ)
c          notice that subroutine f accesses internally defined
c	   quantities in y(neq(1)+1),...
c          if the derivative dy/dt is needed, use intdy instead
c          of calling f for this purpose
c
c itol   = an indicator for the type of error control.  see
c          description below under atol.  used only for input.
c          we have forced itol=1 in the routine drequs.
c          this is reasonable for riccati equations.
c
c rtol   = a relative error tolerance parameter, either a scalar or
c          an array of length neq(1).  see description below under atol.
c          input only.
c
c atol   = an absolute error tolerance parameter, either a scalar or
c          an array of length neq(1).  input only.
c          the input parameters itol, rtol, and atol determine
c          the error control performed by the solver.  the solver will
c          control the vector e = (e(i)) of estimated local errors
c          in y, according to an inequality of the form
c                      max-norm of ( e(i)/ewt(i) )   .le.   1,
c          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),
c          and the max-norm is max-norm(v) = max(abs(v(i), i=1,neqn)  
c          aslo, ewt = (ewt(i))
c          is a vector of weights which must always be positive, and
c          the values of rtol and atol should all be non-negative.
c          the following table gives the types (scalar/array) of
c          rtol and atol, and the corresponding form of ewt(i).
c             itol    rtol       atol          ewt(i)
c              1     scalar     scalar     rtol*abs(y(i)) + atol
c              2     scalar     array      rtol*abs(y(i)) + atol(i)
c              3     array      scalar     rtol(i)*abs(y(i)) + atol
c              4     array      array      rtol(i)*abs(y(i)) + atol(i)
c          when either of these parameters is a scalar, it need not
c          be dimensioned in the user-s calling program.
c          if none of the above choices (with itol, rtol, and atol
c          fixed throughout the problem) is suitable, more general
c          error controls can be obtained by substituting
c          user-supplied routines for the setting of ewt and/or for
c          the norm calculation, as already explained in the      
c          introduction of dresol.
c
c jac    = the name of the internally generated routine (miter=1) to
c          compute the jacobian matrix. 
c          this routine sets up the block of the sylvester 
c          equation to solve at each newton iteration and it has the form
c               subroutine jac (nnq, nnp, x, a011, a012, a022, icontr)
c          jac must be declared external in the calling program.
c          the routine jac is called by the subroutine prepj, and 
c          we refer to the documentation in prepj for further explanation
c          of the available strategy for the newton iteration
c
c-----------------------------------------------------------------------
      external prepj, solsy
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer nnq, nnp, icontr, nonlin, nsuper, nn, nns, nqs, nps,
     1   neqn1, neqn11, neqn12, neqn21, neqn22, np3, np23, np23q,
     2   ntt, nttp, ntt2p, ntt1, ntt1q, ntt2q, nflevc, nflevp
      integer i, i1, i2, iflag, imxer, kgo, lf0,
     1   leniw, lenrw, lenwm,  mord,  mxhnl0, mxstp0, nstext
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real reuvct, reuvpt, tcj, tpj
      real atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0,
     2   r1mach, vnorm
      dimension mord(2)
      logical ihit
c-----------------------------------------------------------------------
c the following internal common block contains
c (a) variables which are local to any subroutine but whose values must
c     be preserved between calls to the routine (own variables), and
c (b) variables which are communicated between subroutines.
c the structure of the block is as follows..  all real variables are
c listed first, followed by all integers.  within each type, the
c variables are grouped with those local to subroutine drequ first,
c then those local to subroutine stode, and finally those used
c for communication.  the block is declared in subroutines
c drequ, intdy, stode, prepj, and solsy.  groups of variables are
c replaced by dummy arrays in the common declarations in routines
c where those variables are not used.
c-----------------------------------------------------------------------
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c
      common /intdre/ reuvct, reuvpt, tcj, tpj, nnq, nnp, icontr, 
     1   nonlin, nsuper, nn, nns, nqs, nps, neqn1, neqn11, neqn12, 
     2   neqn21, neqn22, np3, np23, np23q, ntt, nttp, ntt2p, ntt1, 
     3   ntt1q, ntt2q, nflevc, nflevp
c
      data  mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
c-----------------------------------------------------------------------
c block a.
c this code block is executed on every call.
c it tests istate and itask for legality and branches appropriately.
c if istate .gt. 1 but the flag init shows that initialization has
c not yet been done, an error return occurs.
c if istate = 1 and tout = t, jump to block g and return immediately.
c-----------------------------------------------------------------------
      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .eq. 1) go to 10
      if (init .eq. 0) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   init = 0
      if (tout .eq. t) go to 430
 20   ntrep = 0
c-----------------------------------------------------------------------
c block b.
c the next code block is executed for the initial call (istate = 1),
c or for a continuation call with parameter changes (istate = 3).
c it contains checking of all inputs and various initializations.
c
c first check legality of the non-optional inputs neq, itol, iopt,
c mf.
c-----------------------------------------------------------------------
      if (neq(1) .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .ne. 1) iopt = 0
      meth = mf/10
      miter = mf - 10*meth
      if (meth .lt. 1 .or. meth .gt. 2) go to 608
      if (miter .lt. 0 .or. miter .gt. 1) go to 608
c next process and check the optional inputs. --------------------------
      if (iopt .eq. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      if (istate .eq. 1) h0 = 0.0e0
      hmxi = 0.0e0
      hmin = 0.0e0
      go to 60
 40   maxord = iwork(5)
      if (maxord .lt. 0) go to 611
      if (maxord .eq. 0) maxord = 100
      maxord = min0(maxord,mord(meth))
      mxstep = iwork(6)
      if (mxstep .lt. 0) go to 612
      if (mxstep .eq. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      if (mxhnil .lt. 0) go to 613
      if (mxhnil .eq. 0) mxhnil = mxhnl0
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      if ((tout - t)*h0 .lt. 0.0e0) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. 0.0e0) go to 615
      hmxi = 0.0e0
      if (hmax .gt. 0.0e0) hmxi = 1.0e0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0e0) go to 616
c-----------------------------------------------------------------------
c set work array pointers and check lengths lrw and liw.
c pointers to segments of rwork and iwork are named by prefixing l to
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).
c segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
c-----------------------------------------------------------------------
 60   lyh = 21
      if (istate .eq. 1) nyh = n
      lwm = lyh + (maxord + 1)*nyh
      nstext=3*(nnq+nnp)+2*(nqs+nps)
      if (miter .eq. 0) lenwm = 0
      if (miter .eq. 1) lenwm = nstext + 2
      lewt = lwm + lenwm
      lsavf = lewt + n
      lacor = lsavf + n
      lenrw = lacor + n - 1
      iwork(17) = lenrw
      liwm = 1
      leniw = 20 + nnq +nnp
      if (miter .eq. 0) leniw = 20
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618
c check rtol and atol for legality. ------------------------------------
      rtoli = rtol(1)
      atoli = atol(1)
      do 70 i = 1,n
        if (itol .ge. 3) rtoli = rtol(i)
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        if (rtoli .lt. 0.0e0) go to 619
        if (atoli .lt. 0.0e0) go to 620
 70   continue
      if (istate .eq. 1) go to 100
c if istate = 3, set flag to signal parameter changes to stode. --------
      jstart = -1
      if (nq .le. maxord) go to 90
c maxord was reduced below nq.  copy yh(*,maxord+2) into savf. ---------
      do 80 i = 1,n
 80     rwork(i+lsavf-1) = rwork(i+lwm-1)
c reload wm(1) = rwork(lwm), since lwm may have changed. ---------------
 90   if (miter .gt. 0) rwork(lwm) = sqrt(uround)
      if (n .eq. nyh) go to 200
c neq was reduced.  zero part of yh to avoid undefined references. -----
      i1 = lyh + l*nyh
      i2 = lyh + (maxord + 1)*nyh - 1
      if (i1 .gt. i2) go to 200
      do 95 i = i1,i2
 95     rwork(i) = 0.0e0
      go to 200
c-----------------------------------------------------------------------
c block c.
c the next block is for the initial call only (istate = 1).
c it contains all remaining initializations, the initial call to f,
c and the calculation of the initial step size.
c the error weights in ewt are inverted after being loaded.
c-----------------------------------------------------------------------
 100  uround = r1mach(4)
      tn = t
      if (itask .ne. 4 .and. itask .ne. 5) go to 110
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. 0.0e0) go to 625
      if (h0 .ne. 0.0e0 .and. (t + h0 - tcrit)*h0 .gt. 0.0e0)
     1   h0 = tcrit - t
 110  jstart = 0
      if (miter .gt. 0) rwork(lwm) = sqrt(uround)
      nhnil = 0
      nst = 0
      nje = 0
      nslast = 0
      hu = 0.0e0
      nqu = 0
      ccmax = 0.3e0
      maxcor = 3
      msbp = 20
      mxncf = 10
c initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      call f (t, y, rwork(lf0), probl, rarr, iarr)
      nfe = 1
c load the initial value vector in yh. ---------------------------------
      do 115 i = 1,n
 115    rwork(i+lyh-1) = y(i)
c load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      nq = 1
      h = 1.0e0
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 120 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 621
 120    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
c-----------------------------------------------------------------------
c the coding below computes the step size, h0, to be attempted on the
c first step, unless the user has supplied a value for this.
c first check that tout - t differs significantly from zero.
c a scalar tolerance quantity tol is computed, as max(rtol(i))
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
c so as to be between 100*uround and 1.0e-3.
c then the computed value h0 is given by..
c                                      neq
c   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
c                                       1
c where   w0     = max ( abs(t), abs(tout) ),
c         f(i)   = i-th component of initial value of f,
c         ywt(i) = ewt(i)/tol  (a weight for y(i)).
c the sign of h0 is inferred from the initial values of tout and t.
c-----------------------------------------------------------------------
      if (h0 .ne. 0.0e0) go to 180
      tdist = abs(tout - t)
      w0 = amax1(abs(t),abs(tout))
      if (tdist .lt. 2.0e0*uround*w0) go to 622
      tol = rtol(1)
      if (itol .le. 2) go to 140
      do 130 i = 1,n
 130    tol = amax1(tol,rtol(i))
 140  if (tol .gt. 0.0e0) go to 160
      atoli = atol(1)
      do 150 i = 1,n
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        ayi = abs(y(i))
        if (ayi .ne. 0.0e0) tol = amax1(tol,atoli/ayi)
 150    continue
 160  tol = amax1(tol,100.0e0*uround)
      tol = amin1(tol,0.001e0)
      sum = vnorm (n, rwork(lf0), rwork(lewt))
      sum = 1.0e0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0e0/sqrt(sum)
      h0 = amin1(h0,tdist)
      h0 = sign(h0,tout-t)
c adjust h0 if necessary to meet hmax bound. ---------------------------
 180  rh = abs(h0)*hmxi
      if (rh .gt. 1.0e0) h0 = h0/rh
c load h with h0 and scale yh(*,2) by h0. ------------------------------
      h = h0
      do 190 i = 1,n
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)
      go to 270
c-----------------------------------------------------------------------
c block d.
c the next code block is for continuation calls only (istate = 2 or 3)
c and is to check stop conditions before taking a step.
c-----------------------------------------------------------------------
 200  nslast = nst
      go to (210, 250, 220, 230, 240), itask
 210  if ((tn - tout)*h .lt. 0.0e0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 220  tp = tn - hu*(1.0e0 + 100.0e0*uround)
      if ((tp - tout)*h .gt. 0.0e0) go to 623
      if ((tn - tout)*h .lt. 0.0e0) go to 250
      go to 400
 230  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0e0) go to 624
      if ((tcrit - tout)*h .lt. 0.0e0) go to 625
      if ((tn - tout)*h .lt. 0.0e0) go to 245
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 240  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0e0) go to 624
 245  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0e0 + 4.0e0*uround)
      if ((tnext - tcrit)*h .le. 0.0e0) go to 250
      h = (tcrit - tn)*(1.0e0 - 4.0e0*uround)
      if (istate .eq. 2) jstart = -2
c-----------------------------------------------------------------------
c block e.
c the next block is normally executed for all calls and contains
c the call to the one-step core integrator stode.
c
c this is a looping point for the integration steps.
c
c first check for too many steps being taken, update ewt (if not at
c start of problem), check for too much accuracy being requested, and
c check for h below the roundoff level in t.
c-----------------------------------------------------------------------
 250  continue
      if ((nst-nslast) .ge. mxstep) go to 500
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 260 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 510
 260    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
 270  tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))
      if (tolsf .le. 1.0e0) go to 280
      tolsf = tolsf*2.0e0
      if (nst .eq. 0) go to 626
      go to 520
 280  if ((tn + h) .ne. tn) go to 290
      nhnil = nhnil + 1
      if (nhnil .gt. mxhnil) go to 290
      call xerrwv(50hdresol-- warning..internal t (=r1) and h (=r2) are,
     1   50, 101, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(
     1  60h      such that in the machine, t + h = t on the next step  ,
     1   60, 101, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      (h = step size). solver will continue anyway,
     1   50, 101, 0, 0, 0, 0, 2, tn, h)
      if (nhnil .lt. mxhnil) go to 290
      call xerrwv(50hdresol--  above warning has been issued i1 times  ,
     1   50, 102, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      it will not be issued again for this problem,
     1   50, 102, 0, 1, mxhnil, 0, 0, 0.0e0, 0.0e0)
 290  continue
c-----------------------------------------------------------------------
c     call stode(neq,y,yh,nyh,yh,ewt,savf,acor,wm,iwm,f,
c    1   jac,prepj,solsy,probl,rarr,iarr)
c-----------------------------------------------------------------------
      call stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),
     1   rwork(lsavf), rwork(lacor), rwork(lwm), iwork(liwm),
     2   f, jac, prepj, solsy, probl, rarr, iarr)
      kgo = 1 - kflag
      go to (300, 530, 540), kgo
c-----------------------------------------------------------------------
c block f.
c the following block handles the case of a successful return from the
c core integrator (kflag = 0).  test for stop conditions.
c-----------------------------------------------------------------------
 300  init = 1
      go to (310, 400, 330, 340, 350), itask
c itask = 1.  if tout has been reached, interpolate. -------------------
 310  if ((tn - tout)*h .lt. 0.0e0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
c itask = 3.  jump to exit if tout was reached. ------------------------
 330  if ((tn - tout)*h .ge. 0.0e0) go to 400
      go to 250
c itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
 340  if ((tn - tout)*h .lt. 0.0e0) go to 345
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
 345  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0e0 + 4.0e0*uround)
      if ((tnext - tcrit)*h .le. 0.0e0) go to 250
      h = (tcrit - tn)*(1.0e0 - 4.0e0*uround)
      jstart = -2
      go to 250
c itask = 5.  see if tcrit was reached and jump to exit. ---------------
 350  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
c-----------------------------------------------------------------------
c block g.
c the following block handles all successful returns from drequ.
c if itask .ne. 1, y is loaded from yh and t is set accordingly.
c istate is set to 2, the illegal input counter is zeroed, and the
c optional outputs are loaded into the work arrays before returning.
c if istate = 1 and tout = t, there is a return with no action taken,
c except that if this has happened repeatedly, the run is terminated.
c-----------------------------------------------------------------------
 400  do 410 i = 1,n
 410    y(i) = rwork(i+lyh-1)
      t = tn
      if (itask .ne. 4 .and. itask .ne. 5) go to 420
      if (ihit) t = tcrit
 420  istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c
 430  ntrep = ntrep + 1
      if (ntrep .lt. 5) return
      call xerrwv(
     1  60hdresol--  repeated calls with istate = 1 and tout = t (=r1) ,
     1   60, 301, 0, 0, 0, 0, 1, t, 0.0e0)
      go to 800
c-----------------------------------------------------------------------
c block h.
c the following block handles all unsuccessful returns other than
c those for illegal input.  first the error message routine is called.
c if there was an error test or convergence test failure, imxer is set.
c then y is loaded from yh, t is set to tn, and the illegal input
c counter illin is set to 0.  the optional outputs are loaded into
c the work arrays before returning.
c-----------------------------------------------------------------------
c the maximum number of steps was taken before reaching tout. ----------
 500  call xerrwv(50hdresol--  at current t (=r1), mxstep (=i1) steps  ,
     1   50, 201, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      taken on this call before reaching tout     ,
     1   50, 201, 0, 1, mxstep, 0, 1, tn, 0.0e0)
      istate = -1
      go to 580
c ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  ewti = rwork(lewt+i-1)
      call xerrwv(50hdresol-- at t (=r1), ewt(i1) has become r2 .le. 0.,
     1   50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 580
c too much accuracy requested for machine precision. -------------------
 520  call xerrwv(50hdresol--  at t (=r1), too much accuracy requested ,
     1   50, 203, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      for precision of machine..  see tolsf (=r2) ,
     1   50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580
c kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
 530  call xerrwv(50hdresol-- at t(=r1) and step size h(=r2), the error,
     1   50, 204, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      test failed repeatedly or with abs(h) = hmin,
     1   50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 560
c kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
 540  call xerrwv(50hdresol--  at t (=r1) and step size h (=r2), the   ,
     1   50, 205, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      corrector convergence failed repeatedly     ,
     1   50, 205, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(30h      or with abs(h) = hmin   ,
     1   30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
c compute imxer if relevant. -------------------------------------------
 560  big = 0.0e0
      imxer = 1
      do 570 i = 1,n
        size = abs(rwork(i+lacor-1)*rwork(i+lewt-1))
        if (big .ge. size) go to 570
        big = size
        imxer = i
 570    continue
      iwork(16) = imxer
c set y vector, t, illin, and optional outputs. ------------------------
 580  do 590 i = 1,n
 590    y(i) = rwork(i+lyh-1)
      t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c-----------------------------------------------------------------------
c block i.
c the following block handles all error returns due to illegal input
c (istate = -3), as detected before calling the core integrator.
c first the error message routine is called.  then if there have been
c 5 consecutive such returns just before this call to the solver,
c the run is halted.
c-----------------------------------------------------------------------
 601  call xerrwv(30hdresol--  istate (=i1) illegal,
     1   30, 1, 0, 1, istate, 0, 0, 0.0e0, 0.0e0)
      go to 700
 602  call xerrwv(30hdresol--  itask (=i1) illegal ,
     1   30, 2, 0, 1, itask, 0, 0, 0.0e0, 0.0e0)
      go to 700
 603  call xerrwv(50hdresol--  istate .gt. 1 but dresol not initialized,
     1   50, 3, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      go to 700
 604  call xerrwv(30hdresol--  neq (=i1) .lt. 1    ,
     1   30, 4, 0, 1, neq(1), 0, 0, 0.0e0, 0.0e0)
      go to 700
 605  call xerrwv(50hdresol--  istate = 3 and neq increased (i1 to i2) ,
     1   50, 5, 0, 2, n, neq(1), 0, 0.0e0, 0.0e0)
      go to 700
 606  call xerrwv(30hdresol--  itol (=i1) illegal  ,
     1   30, 6, 0, 1, itol, 0, 0, 0.0e0, 0.0e0)
      go to 700
 608  call xerrwv(30hdresol--  mf (=i1) illegal    ,
     1   30, 8, 0, 1, mf, 0, 0, 0.0e0, 0.0e0)
      go to 700
 611  call xerrwv(30hdresol--  maxord (=i1) .lt. 0 ,
     1   30, 11, 0, 1, maxord, 0, 0, 0.0e0, 0.0e0)
      go to 700
 612  call xerrwv(30hdresol--  mxstep (=i1) .lt. 0 ,
     1   30, 12, 0, 1, mxstep, 0, 0, 0.0e0, 0.0e0)
      go to 700
 613  call xerrwv(30hdresol--  mxhnil (=i1) .lt. 0 ,
     1   30, 13, 0, 1, mxhnil, 0, 0, 0.0e0, 0.0e0)
      go to 700
 614  call xerrwv(40hdresol--  tout (=r1) behind t (=r2)     ,
     1   40, 14, 0, 0, 0, 0, 2, tout, t)
      call xerrwv(50h      integration direction is given by h0 (=r1)  ,
     1   50, 14, 0, 0, 0, 0, 1, h0, 0.0e0)
      go to 700
 615  call xerrwv(30hdresol--  hmax (=r1) .lt. 0.0 ,
     1   30, 15, 0, 0, 0, 0, 1, hmax, 0.0e0)
      go to 700
 616  call xerrwv(30hdresol--  hmin (=r1) .lt. 0.0 ,
     1   30, 16, 0, 0, 0, 0, 1, hmin, 0.0e0)
      go to 700
 617  call xerrwv(
     1  60hdrequ--  rwork length needed, lenrw (=i1), exceeds lrw (=i2),
     1   60, 17, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 618  call xerrwv(
     1  60hdrequ--  iwork length needed, leniw (=i1), exceeds liw (=i2),
     1   60, 18, 0, 2, leniw, liw, 0, 0.0e0, 0.0e0)
      go to 700
 619  call xerrwv(40hdrequ--  rtol(i1) is r1 .lt. 0.0        ,
     1   40, 19, 0, 1, i, 0, 1, rtoli, 0.0e0)
      go to 700
 620  call xerrwv(40hdrequ--  atol(i1) is r1 .lt. 0.0        ,
     1   40, 20, 0, 1, i, 0, 1, atoli, 0.0e0)
      go to 700
 621  ewti = rwork(lewt+i-1)
      call xerrwv(40hdresol--  ewt(i1) is r1 .le. 0.0        ,
     1   40, 21, 0, 1, i, 0, 1, ewti, 0.0e0)
      go to 700
 622  call xerrwv(
     1  60hdresol-- tout (=r1) too close to t(=r2) to start integration,
     1   60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  call xerrwv(
     1  60hdresol--  itask = i1 and tout (=r1) behind tcur - hu (= r2) ,
     1   60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624  call xerrwv(
     1  60hdresol--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)  ,
     1   60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  call xerrwv(
     1  60hdresol--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)  ,
     1   60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  call xerrwv(50hdresol--  at start of problem, too much accuracy  ,
     1   50, 26, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(
     1  60h      requested for precision of machine..  see tolsf (=r1) ,
     1   60, 26, 0, 0, 0, 0, 1, tolsf, 0.0e0)
      rwork(14) = tolsf
      go to 700
 627  call xerrwv(50hdresol-- trouble from intdy. itask = i1, tout = r1,
     1   50, 27, 0, 1, itask, 0, 1, tout, 0.0e0)
c
 700  if (illin .eq. 5) go to 710
      illin = illin + 1
      istate = -3
      return
 710  call xerrwv(50hdresol--  repeated occurrences of illegal input   ,
     1   50, 302, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
c
 800  call xerrwv(50hdresol--  run aborted.. apparent infinite loop    ,
     1   50, 303, 2, 0, 0, 0, 0, 0.0e0, 0.0e0)
      return
c----------------------- end of subroutine drequ ----------------------
      end
      block data atlan1
c
c revision date: July 10 1990 
c author: Luca Dieci, School of Mathematics, GaTech
c
c this block data is needed for initialization of the 
c quantities illin and ntrep.  In the original LSODE
c these quantities were initialized via a DATA statement, but this
c is not allowed on the SUN nor on the IBM (it is fine on the CDC)
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      data illin /0/, ntrep /0/
c----------------------- end of block-data atlan1 ----------------------
      end
      subroutine intdy (t, k, yh, nyh, dky, iflag)
c 
c author: Alan Hindmarsh, LLL, in the original LSODE package
c
      integer k, nyh, iflag
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      real t, yh, dky
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real c, r, s, tp
      dimension yh(nyh,*), dky(*)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c intdy computes interpolated values of the k-th derivative of the
c dependent variable vector y, and stores it in dky.  this routine
c is called within the package with k = 0 and t = tout, but may
c also be called by the user for any k up to the current order.
c (see detailed instructions in the usage documentation.)
c-----------------------------------------------------------------------
c the computed values in dky are gotten by interpolation using the
c nordsieck history array yh.  this array corresponds uniquely to a
c vector-valued polynomial of degree nqcur or less, and dky is set
c to the k-th derivative of this polynomial at t.
c the formula for dky is..
c              q
c  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
c             j=k
c where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
c the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
c communicated by common.  the above sum is done in reverse order.
c iflag is returned negative if either k or t is out of bounds.
c-----------------------------------------------------------------------
      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu -  100.0e0*uround*(tn + hu)
      if ((t-tp)*(t-tn) .gt. 0.0e0) go to 90
c
      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = float(ic)
      do 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1,j
 30       ic = ic*jj
 35     c = float(ic)
        do 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
 60     dky(i) = r*dky(i)
      return
c
 80   call xerrwv(30hintdy--  k (=i1) illegal      ,
     1   30, 51, 0, 1, k, 0, 0, 0.0e0, 0.0e0)
      iflag = -1
      return
 90   call xerrwv(30hintdy--  t (=r1) illegal      ,
     1   30, 52, 0, 0, 0, 0, 1, t, 0.0e0)
      call xerrwv(
     1  60h      t not in interval tcur - hu (= r1) to tcur (=r2)      ,
     1   60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      return
c----------------------- end of subroutine intdy -----------------------
      end
      subroutine stode (neq, y, yh, nyh, yh1, ewt, savf, acor,
     1   wm, iwm, f, jac, pjac, slvs, probl, rarr, iarr)
c
c revision date: July 10 1990
c modifications' author: Luca Dieci, School of Mathematics, GaTech
c
c history: based on the routine by the same name, which is in the 
c          lsode package of A.Hindmarsh.  This version has many 
c          changes 
c 
      external f, jac, pjac, slvs, probl
      integer neq, nyh, iwm, iarr
      integer iownd, ialth, ipup, lmax, meo, nqnyh, nslp,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer nnq, nnp, icontr, nonlin, nsuper, nn, nns, nqs, nps,
     1   neqn1, neqn11, neqn12, neqn21, neqn22, np3, np23, np23q,
     2   ntt, nttp, ntt2p, ntt1, ntt1q, ntt2q, nflevc, nflevp
      integer i, i1, iredo, iret, j, jb, m, ncf, newq
      integer ksup
      real y, yh, yh1, ewt, savf, acor, wm, rarr
      real conit, crate, el, elco, hold, rmax, tesco,
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real reuvct, reuvpt, tcj, tpj
      real dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
     1   r, rh, rhdn, rhsm, rhup, told, vnorm
      real tpjph, signh
      dimension neq(*), y(*), yh(nyh,*), yh1(*), ewt(*), savf(*),
     1   acor(*), wm(*), rarr(*), iwm(*), iarr(*)
      common /ls0001/ conit, crate, el(13), elco(13,12),
     1   hold, rmax, tesco(3,12),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14),
     3   ialth, ipup, lmax, meo, nqnyh, nslp,
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /intdre/ reuvct, reuvpt, tcj, tpj, nnq, nnp, icontr, 
     1   nonlin, nsuper, nn, nns, nqs, nps, neqn1, neqn11, neqn12, 
     2   neqn21, neqn22, np3, np23, np23q, ntt, nttp, ntt2p, ntt1, 
     3   ntt1q, ntt2q, nflevc, nflevp
c-----------------------------------------------------------------------
c stode performs one step of the integration of an initial value
c problem for a system of ordinary differential equations.
c note.. stode is independent of the value of the iteration method
c indicator miter, when this is .ne. 0, and hence is independent
c of the type of chord method used.
c communication with stode is done with the following variables..
c
c neq    = integer array containing problem size in neq(1), and
c          passed as the neq argument in all calls to f and jac.
c y      = an array of length .ge. n used as the y argument in
c          all calls to f and jac.
c yh     = an nyh by lmax array containing the dependent variables
c          and their approximate scaled derivatives, where
c          lmax = maxord + 1.  yh(i,j+1) contains the approximate
c          j-th derivative of y(i), scaled by h**j/factorial(j)
c          (j = 0,1,...,nq).  on entry for the first step, the first
c          two columns of yh must be set from the initial values.
c nyh    = a constant integer .ge. n, the first dimension of yh.
c yh1    = a one-dimensional array occupying the same space as yh.
c ewt    = an array of length n containing multiplicative weights
c          for local error measurements.  local errors in y(i) are
c          compared to 1.0/ewt(i) in various error tests.
c savf   = an array of working storage, of length n.
c          also used for input of yh(*,maxord+2) when jstart = -1
c          and maxord .lt. the current order nq.
c acor   = a work array of length n, used for the accumulated
c          corrections.  on a successful return, acor(i) contains
c          the estimated one-step local error in y(i).
c wm,iwm = real and integer work arrays associated with matrix
c          operations in chord iteration (miter .ne. 0).
c pjac   = name of routine to evaluate and preprocess jacobian matrix
c          and p = i - h*el0*jac, if a chord method is being used.
c slvs   = name of routine to solve linear system in chord iteration.
c    N.B.: the routines pjac and slvs actually consist in simple
c          interfaces calling various other routines.  See their
c          own documentation 
c ccmax  = maximum relative change in h*el0 before pjac is called.
c h      = the step size to be attempted on the next step.
c          h is altered by the error control algorithm during the
c          problem.  h can be either positive or negative, but its
c          sign must remain constant throughout the problem.
c hmin   = the minimum absolute value of the step size h to be used.
c hmxi   = inverse of the maximum absolute value of h to be used.
c          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
c          hmin and hmxi may be changed at any time, but will not
c          take effect until the next change of h is considered.
c tn     = the independent variable. tn is updated on each step taken.
c jstart = an integer used for input only, with the following
c          values and meanings..
c               0  perform the first step.
c           .gt.0  take a new step continuing from the last.
c              -1  take the next step with a new value of h, maxord,
c                    n, meth, miter, and/or matrix parameters.
c              -2  take the next step with a new value of h,
c                    but with other inputs unchanged.
c          on return, jstart is set to 1 to facilitate continuation.
c kflag  = a completion code with the following meanings..
c               0  the step was succesful.
c              -1  the requested error could not be achieved.
c              -2  corrector convergence could not be achieved.
c              -3  fatal error in pjac or slvs.
c          a return with kflag = -1 or -2 means either
c          abs(h) = hmin or 10 consecutive failures occurred.
c          on a return with kflag negative, the values of tn and
c          the yh array are as of the beginning of the last
c          step, and h is the last step size attempted.
c maxord = the maximum order of integration method to be allowed.
c maxcor = the maximum number of corrector iterations allowed.
c msbp   = maximum number of steps between pjac calls (miter .gt. 0).
c mxncf  = maximum number of convergence failures allowed.
c meth/miter = the method flags.  see description in driver.
c n      = the number of first-order differential equations.
c
c --- a number of additions have been made in this routine with respect       
c  the original one of LSODE.  In particular, we handle here the option
c  of how to solve the nolinear system and the superstability option/check
c  based upon eigenvalue information.  Communication is provided by the 
c  following variables in the block /intdre/
c
c reuvct = most unstable real part of the eigenvalues at current time t
c reuvpt = most unstable real part of the eigenvalues at previous time t
c tcj = time at which the most current evaluation/factorization of the
c       Jacobian matrix has been made
c tpj = next to the last time at which an evaluation/factorization of the
c       Jacobian martix has been made
c nsuper = integer flag signalling whether the code should (nsuper = 1) 
c          perform a superstability check or not (nsuper = 0).  For the 
c          criterion adopted we refer to the documentation of DRESOL
c nflevc = integer flaf monitoring whether at the current step we had 
c          (nflevc = 1) or not (nflevc = 0) unstable eigenvalues with
c          respect to the direction of integration and the machine precision
c nflevp = integer flag as nflevc, but recording unstable/stable eigenvalues
c          at the previous step
c nonlin = the way we should solve the nonlinear system of equations in
c          case miter = 1 (i.e. mf=11 or 21).  We set nonlin = 1 if we
c          want a real quasi-Newton iteration (that is, a factorizatio
c          at each step and then freeze the jacobian for that one step)
c          and we set nonlin = 0 if we want a chord iteration with a 
c          factored jacobian being kept across many steps (depending on msbp)
c
c all of the other parameters are as explained in dresol
c-----------------------------------------------------------------------
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0e0
      if (jstart .eq. 0) tcj = tn
      tpj = tcj
      nflevp = nflevc
      reuvpt = reuvct
      ksup = 0
      if (jstart .gt. 0) go to 200
      if (jstart .eq. -1) go to 100
      if (jstart .eq. -2) go to 160
c-----------------------------------------------------------------------
c on the first call, the order is set to 1, and other variables are
c initialized.  rmax is the maximum ratio by which h can be increased
c in a single step.  it is initially 1.e4 to compensate for the small
c initial h, but then is normally equal to 10.  if a failure
c occurs (in corrector convergence or error test), rmax is set at 2
c for the next increase.
c-----------------------------------------------------------------------
      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0e0
      rc = 0.0e0
      el0 = 1.0e0
      crate = 0.7e0
      hold = h
      meo = meth
      nslp = 0
      ipup = miter
      iret = 3
      go to 140
c-----------------------------------------------------------------------
c the following block handles preliminaries needed when jstart = -1.
c ipup is set to miter to force a matrix update.
c if an order increase is about to be considered (ialth = 1),
c ialth is reset to 2 to postpone consideration one more step.
c if the caller has changed meth, cfode is called to reset
c the coefficients of the method.
c if the caller has changed maxord to a value less than the current
c order nq, nq is reduced to maxord, and a new h chosen accordingly.
c if h is to be changed, yh must be rescaled.
c if h or meth is being changed, ialth is reset to l = nq + 1
c to prevent further changes in h for that many steps.
c-----------------------------------------------------------------------
 100  ipup = miter
      lmax = maxord + 1
      if (ialth .eq. 1) ialth = 2
      if (meth .eq. meo) go to 110
      call cfode (meth, elco, tesco)
      meo = meth
      if (nq .gt. maxord) go to 120
      ialth = l
      iret = 1
      go to 150
 110  if (nq .le. maxord) go to 160
 120  nq = maxord
      l = lmax
      do 125 i = 1,l
 125    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)
      ddn = vnorm (n, savf, ewt)/tesco(1,l)
      exdn = 1.0e0/float(l)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
      rh = amin1(rhdn,1.0e0)
      iredo = 3
      if (h .eq. hold) go to 170
      rh = amin1(rh,abs(h/hold))
      h = hold
      go to 175
c-----------------------------------------------------------------------
c cfode is called to get all the integration coefficients for the
c current meth.  then the el vector and related constants are reset
c whenever the order nq is changed, or at the start of the problem.
c-----------------------------------------------------------------------
 140  call cfode (meth, elco, tesco)
 150  do 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)
      go to (160, 170, 200), iret
c-----------------------------------------------------------------------
c if h is being changed, the h ratio rh is checked against
c rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
c l = nq + 1 to prevent a change of h for that many steps, unless
c forced by a convergence or error test failure.
c-----------------------------------------------------------------------
 160  if (h .eq. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = amax1(rh,hmin/abs(h))
 175  rh = amin1(rh,rmax)
      rh = rh/amax1(1.0e0,abs(h)*hmxi*rh)
      r = 1.0e0
      do 180 j = 2,l
        r = r*rh
        do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 690
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the yh array by the pascal triangle matrix.
c rc is the ratio of new to old values of the coefficient  h*el(1).
c when rc differs from 1 by more than ccmax, ipup is set to miter
c to force pjac to be called, if a jacobian is involved.
c in any case, pjac is called at least every msbp steps.
c if nonlin = 1 pjac is called every step.
c-----------------------------------------------------------------------
 200  if (abs(rc-1.0e0) .gt. ccmax) ipup = miter
      if (nst .ge. nslp+msbp) ipup = miter
c next line forces a genuine quasi-newton iteration, with 
c the jacobian updated at each new step, if nonlin=1  
      if (nonlin .eq. 1) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do 215 jb = 1,nq
        i1 = i1 - nyh
        do 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    continue
c-----------------------------------------------------------------------
c up to maxcor corrector iterations are taken.  a convergence test is
c made on the max norm of each correction, weighted by the error
c weight vector ewt.  the sum of the corrections is accumulated in the
c vector acor(i).  the yh array is not altered in the corrector loop.
c-----------------------------------------------------------------------
 220  m = 0
      do 230 i = 1,n
 230    y(i) = yh(i,1)
      call f (tn, y, savf, probl, rarr, iarr)
      nfe = nfe + 1
      if (ipup .le. 0) go to 250
c-----------------------------------------------------------------------
c if indicated, the matrix p = i - h*el(1)*j is reevaluated and
c preprocessed before starting the corrector iteration.  ipup is set
c to 0 as an indicator that this has been done.  see pjac for how this is done.
c-----------------------------------------------------------------------
      call pjac (y,wm,iwm,probl,jac,rarr,iarr)
      ipup = 0
      rc = 1.0e0
      nslp = nst
      crate = 0.7e0
      if (ierpj .ne. 0) go to 430
 250  do 260 i = 1,n
 260    acor(i) = 0.0e0
 270  if (miter .ne. 0) go to 350
c-----------------------------------------------------------------------
c in the case of functional iteration, update y directly from
c the result of the last function evaluation.
c-----------------------------------------------------------------------
      do 290 i = 1,n
        savf(i) = h*savf(i) - yh(i,2)
 290    y(i) = savf(i) - acor(i)
      del = vnorm (n, y, ewt)
      do 300 i = 1,n
        y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400
c-----------------------------------------------------------------------
c in the case of the chord or quasi-Newton method, compute the corrector 
c error, and solve the linear system with that as right-hand side and
c p as coefficient matrix.
c-----------------------------------------------------------------------
 350  do 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))
      call slvs (wm, y)
      if (iersl .lt. 0) go to 430
      if (iersl .gt. 0) go to 410
      del = vnorm (n, y, ewt)
      do 380 i = 1,n
        acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)
c-----------------------------------------------------------------------
c test for convergence.  if m.gt.0, an estimate of the convergence
c rate constant is stored in crate, and this is used in the test.
c-----------------------------------------------------------------------
 400  continue
      if (m .ne. 0) crate = amax1(0.2e0*crate,del/delp)
      dcon = del*amin1(1.0e0,1.5e0*crate)/(tesco(2,nq)*conit)
      if (dcon .le. 1.0e0) go to 450
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. 2.0e0*delp) go to 410
      delp = del
      call f (tn, y, savf, probl, rarr, iarr)
      nfe = nfe + 1
      go to 270
c-----------------------------------------------------------------------
c the corrector iteration failed to converge.
c if miter .ne. 0 and the jacobian is out of date, pjac is called for
c the next try.  otherwise the yh array is retracted to its values
c before prediction, and h is reduced, if possible.  if h cannot be
c reduced or mxncf failures have occurred, exit with kflag = -2.
c-----------------------------------------------------------------------
 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0e0
      tn = told
      i1 = nqnyh + 1
      do 445 jb = 1,nq
        i1 = i1 - nyh
        do 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    continue
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
      if (abs(h) .le. hmin*1.00001e0) go to 670
      if (ncf .eq. mxncf) go to 670
      rh = 0.25e0
      ipup = miter
      iredo = 1
      go to 170
c-----------------------------------------------------------------------
c the corrector has converged.  jcur is set to 0
c to signal that the jacobian involved may need updating later.
c the local error test is made and control passes to statement 500
c if it fails.
c-----------------------------------------------------------------------
 450  jcur = 0
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vnorm (n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0e0) go to 500
c-----------------------------------------------------------------------
c after a successful step, update the yh array.
c consider changing h if ialth = 1.  otherwise decrease ialth by 1.
c if ialth is then 1 and nq .lt. maxord, then acor is saved for
c use in a possible order increase on the next step.
c if a change in h is considered, an increase or decrease in order
c by one is considered also.  a change in h is made only if it is by a
c factor of at least 1.1.  if not, ialth is set to 3 to prevent
c testing for that many steps.
c-----------------------------------------------------------------------
c	       ----- superstability/eigenvalues check -----
c the last decomposition of the jacobian did not occur at current t,
c or we do not have the eigenvalues altogether
      if (tn .ne. tcj .or. miter .eq. 0) goto 466
c ready to check the sign of eigenvalues of jacobian
c do this regardless of superstability check, to signal
c possible instabilities in integration to the user
c in particular, iwm(1) = 1 means purely imaginary e.vs.  it is not
c error but user might want to use adams formulas if he is using
c the bdf.  iwm(2) = 1 means unstable real parts of eigenvalues
c with respect to direction of integration: sign(h) * real-part(e.vs) > 0
c initially, have iwm(1) = iwm(2) = 0.
      if (iwm(1) .eq. 1) goto 457
      iwm(1) = 20
      if (icontr .eq. 1) goto 454
      do 452 i = 1,nnp
            if ((iwm(20+i) .ne. 0) .and. (wm(nttp+i-1) .eq. 0))
     1           iwm(1) = 10
 452  continue
      goto 457
 454  do 455 i = 1,nnp
         do 455 j = 1,nnq
            if ((iwm(20+i) + iwm(20+nnp+j) .ne. 0) .and.
     1   (wm(nttp+i-1) - wm(ntt1q+j-1) .eq. 0)) iwm(1) = 10
 455  continue
c above, we were looking at itypea(i), itypeb(j) of the routine hqr3 (see pjac)
 457  tpjph = tpj + h
      if ((iwm(2) .eq. 1) .and. (nsuper .eq. 0 .or.
     1   tcj .ne. tpjph)) goto 466
      if (iwm(2) .ne. 1) iwm(2) = 20
      signh = -1.e0
      if (h .ge. 0.e0) signh = 1.e0
      nflevc = 0
      if (signh .lt. 0.e0) goto 460
c we check whether the most unstable eigenvalue is unstable with
c respect to the direction of integration
c in case of icontr=0, the eigenvalues are ordered in incresing
c order of real parts, and they have the reversed sign.
      if (icontr .eq. 0) reuvct = - wm(nttp)
      if (icontr .eq. 1) reuvct = wm(nttp+nnp-1)-wm(ntt1q) 
      goto 464
 460  if (icontr .eq. 0) reuvct = - wm(nttp+nnp-1) 
      if (icontr .eq. 1) reuvct = wm(nttp)-wm(ntt1q+nnq-1) 
 464  reuvct = signh * reuvct
      if (reuvct .gt. 1.0e-3) nflevc = 1
      if (nflevc .eq. 1 .and. iwm(2) .ne. 1) iwm(2) = 10
      if (nsuper .eq. 0 .or. tcj .ne. tpjph) goto 466
c check for possible superstable behavior if using the BDF or Adams
c with the Newton-like iteration (no e.vs are available, when
c using the functional fixed point iteration).
c the logic behind this check relies on (i) the assumed existence of
c a stable trajectory, or (ii) that, if we are along an unstable trajectory,
c then we have to proceed cautiously; in particular, we are not allowed
c to increase the stepsize if the degree of instability (as measured
c by the most unstable real part of the eigenvalues) is increased 
c by 100% or more. If neither of these two things is true, the code might 
c end up working hard for no reason at all.  We have judged as unstable
c an eigenvalue if it is positive w.r.t. h, within heuristic accuracy
c bound (1.0e-3, above).  
      if (nflevc .eq. 1 .and. (abs(h) .ge. abs(hu) .or. 
     1   ksup .gt. 0) .and. ((nflevp .eq. 0 .and. reuvct .gt. 1.e-2)
     2   .or. reuvct .gt. (2.0e0*reuvpt))) then
         ksup = ksup+1
         goto 465
      endif
c in the process of avoiding superstability, we have reduced the stepsize
c within machine precision without recovering from the possible occurrence
c of the phenomenon. To avoid further testing, we set nsuper to 0.  The
c user should monitor this, if wants to change it again to snuper=1 later
      if (tpjph .eq. tpj) then
         call xerrwv(47hdresol-- in subroutine stode reset nsuper to 0 ,
     1      47,11,1,0,0,0,0,0.e0,0.e0)
         call xerrwv(47h    as after i1 reductions h (=r1) is too small,
     1      47,12,1,1,ksup,0,1,h,0.e0)
         nsuper = 0
         neq(6) = 0
      endif
      ksup = 0
c if cannot recover from superstable behavior, after 6 steps reductions
c signal it to the user by letting nsuper = 0
 465  if (ksup .gt. 6) then
         call xerrwv(47hdresol-- in subroutine stode reset nsuper to 0 ,
     1      47,13,1,0,0,0,0,0.e0,0.e0)
         call xerrwv(47h    after 6 tries still superstability possible,
     1      47,14,1,0,0,0,0,0.e0,0.e0)
         nsuper = 0
         neq(6) = 0
         ksup = 0
      endif
      if (ksup .eq. 0) goto 466
c above, either the superstability check failed, or we made it more
c than 6 times without recovering.  in these cases, accept the step
c because the error tols were fine.  otherwise, force a step-size
c reduction.  after 3 failures (because of ksup or err-tols) we force
c a reduction by at least a factor of 10 at order 1.
c if we detect superstability, the step-size reduction is enforced  by
c setting dsm=10.  this should give -at first- reduction at same order
      dsm = 10
      goto 500
c ----- end of superstability/eigenvalues' check -----
 466  kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      do 470 j = 1,l
        do 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)
      ialth = ialth - 1
      if (ialth .eq. 0) go to 520
      if (ialth .gt. 1) go to 700
      if (l .eq. lmax) go to 700
      do 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 700
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c restore tn and the yh array to their previous values, and prepare
c to try the step again.  compute the optimum step size for this or
c one lower order.  after 2 or more failures, h is forced to decrease
c by a factor of 0.2 or less.
c-----------------------------------------------------------------------
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
        i1 = i1 - nyh
        do 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    continue
      rmax = 2.0e0
      if (abs(h) .le. hmin*1.00001e0) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0e0
      go to 540
c-----------------------------------------------------------------------
c regardless of the success or failure of the step, factors
c rhdn, rhsm, and rhup are computed, by which h could be multiplied
c at order nq - 1, order nq, or order nq + 1, respectively.
c in the case of failure, rhup = 0.0 to avoid an order increase.
c the largest of these is determined and the new order chosen
c accordingly.  if the order is to be increased, we compute one
c additional scaled derivative.
c-----------------------------------------------------------------------
 520  rhup = 0.0e0
      if (l .eq. lmax) go to 540
      do 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vnorm (n, savf, ewt)/tesco(3,nq)
      exup = 1.0e0/float(l+1)
      rhup = 1.0e0/(1.4e0*dup**exup + 0.0000014e0)
 540  exsm = 1.0e0/float(l)
      rhsm = 1.0e0/(1.2e0*dsm**exsm + 0.0000012e0)
      rhdn = 0.0e0
      if (nq .eq. 1) go to 560
      ddn = vnorm (n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0e0/float(nq)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
 560  if (rhsm .ge. rhup) go to 570
      if (rhup .gt. rhdn) go to 590
      go to 580
 570  if (rhsm .lt. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      if (kflag .lt. 0 .and. rh .gt. 1.0e0) rh = 1.0e0
      go to 620
 590  newq = l
      rh = rhup
      if (rh .lt. 1.1e0) go to 610
      r = el(l)/float(l)
      do 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 700
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1e0)) go to 610
      if (kflag .le. -2) rh = amin1(rh,0.2e0)
c-----------------------------------------------------------------------
c if there is a change of order, reset nq, l, and the coefficients.
c in any case h is reset according to rh and the yh array is rescaled.
c then exit from 690 if the step was ok, or redo the step otherwise.
c-----------------------------------------------------------------------
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
c-----------------------------------------------------------------------
c control reaches this section if 3 or more failures have occured.
c if 10 failures have occurred, exit with kflag = -1.
c it is assumed that the derivatives that have accumulated in the
c yh array have errors of the wrong order.  hence the first
c derivative is recomputed, and the order is set to 1.  then
c h is reduced by a factor of 10, and the step is retried,
c until it succeeds or h reaches hmin.
c-----------------------------------------------------------------------
 640  if (kflag .eq. -10) go to 660
      rh = 0.1e0
      rh = amax1(hmin/abs(h),rh)
      h = h*rh
      do 645 i = 1,n
 645    y(i) = yh(i,1)
      call f (tn, y, savf, probl, rarr, iarr)
      nfe = nfe + 1
      do 650 i = 1,n
 650    yh(i,2) = h*savf(i)
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
c-----------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
 660  kflag = -1
      go to 720
 670  kflag = -2
      go to 720
 680  kflag = -3
      go to 720
 690  rmax = 10.0e0
 700  r = 1.0e0/tesco(2,nqu)
      do 710 i = 1,n
 710    acor(i) = acor(i)*r
 720  hold = h
      jstart = 1
      if (iwm(1) .eq. 1 .and. iwm(2) .eq. 1) return
      if (iwm(1) .eq. 10) iwm(1) = 1
      if (iwm(1) .eq. 20) iwm(1) = 0
      if (iwm(2) .eq. 10) iwm(2) = 1
      if (iwm(2) .eq. 20) iwm(2) = 0
      return
c----------------------- end of subroutine stode -----------------------
      end
      subroutine cfode (meth, elco, tesco)
c
c author: Alan Hindmarsh, LLL
c
      integer meth
      integer i, ib, nq, nqm1, nqp1
      real elco, tesco
      real agamq, fnq, fnqm1, pc, pint, ragq,
     1   rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)
c-----------------------------------------------------------------------
c cfode is called by the integrator routine to set coefficients
c needed there.  the coefficients for the current method, as
c given by the value of meth, are set for all orders and saved.
c the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
c (a smaller value of the maximum order is also allowed.)
c cfode is called once at the beginning of the problem,
c and is not called again unless and until meth is changed.
c
c the elco array contains the basic method coefficients.
c the coefficients el(i), 1 .le. i .le. nq+1, for the method of
c order nq are stored in elco(i,nq).  they are given by a genetrating
c polynomial, i.e.,
c     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
c for the implicit adams methods, l(x) is given by
c     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
c for the bdf methods, l(x) is given by
c     l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
c where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
c
c the tesco array contains test constants used for the
c local error test and the selection of step size and/or order.
c at order nq, tesco(k,nq) is used for the selection of step
c size at order nq - 1 if k = 1, at order nq if k = 2, and at order
c nq + 1 if k = 3.
c-----------------------------------------------------------------------
      dimension pc(12)
c
      go to (100, 200), meth
c
 100  elco(1,1) = 1.0e0
      elco(2,1) = 1.0e0
      tesco(1,1) = 0.0e0
      tesco(2,1) = 2.0e0
      tesco(1,2) = 1.0e0
      tesco(3,12) = 0.0e0
      pc(1) = 1.0e0
      rqfac = 1.0e0
      do 140 nq = 2,12
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq-1).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        rq1fac = rqfac
        rqfac = rqfac/float(nq)
        nqm1 = nq - 1
        fnqm1 = float(nqm1)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq-1). ----------------------------------
        pc(nq) = 0.0e0
        do 110 ib = 1,nqm1
          i = nqp1 - ib
 110      pc(i) = pc(i-1) + fnqm1*pc(i)
        pc(1) = fnqm1*pc(1)
c compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        pint = pc(1)
        xpin = pc(1)/2.0e0
        tsign = 1.0e0
        do 120 i = 2,nq
          tsign = -tsign
          pint = pint + tsign*pc(i)/float(i)
 120      xpin = xpin + tsign*pc(i)/float(i+1)
c store coefficients in elco and tesco. --------------------------------
        elco(1,nq) = pint*rq1fac
        elco(2,nq) = 1.0e0
        do 130 i = 2,nq
 130      elco(i+1,nq) = rq1fac*pc(i)/float(i)
        agamq = rqfac*xpin
        ragq = 1.0e0/agamq
        tesco(2,nq) = ragq
        if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/float(nqp1)
        tesco(3,nqm1) = ragq
 140    continue
      return
c
 200  pc(1) = 1.0e0
      rq1fac = 1.0e0
      do 230 nq = 1,5
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        fnq = float(nq)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq). ------------------------------------
        pc(nqp1) = 0.0e0
        do 210 ib = 1,nq
          i = nq + 2 - ib
 210      pc(i) = pc(i-1) + fnq*pc(i)
        pc(1) = fnq*pc(1)
c store coefficients in elco and tesco. --------------------------------
        do 220 i = 1,nqp1
 220      elco(i,nq) = pc(i)/pc(2)
        elco(2,nq) = 1.0e0
        tesco(1,nq) = rq1fac
        tesco(2,nq) = float(nqp1)/elco(1,nq)
        tesco(3,nq) = float(nq+2)/elco(1,nq)
        rq1fac = rq1fac/fnq
 230    continue
      return
c----------------------- end of subroutine cfode -----------------------
      end
c
      subroutine f(t,y,ydot,probl,rarr,iarr)
c
c revision date: July 10 1990
c author: Luca Dieci, School of Mathematics, GaTech
c
      external probl
      integer iarr
      real t, y, ydot, rarr
      dimension y(*),ydot(*),rarr(*),iarr(*)
c
c we give the function here
c
      integer nnq, nnp, icontr, nonlin, nsuper, nn, nns, nqs, nps,
     1   neqn1, neqn11, neqn12, neqn21, neqn22, np3, np23, np23q,
     2   ntt, nttp, ntt2p, ntt1, ntt1q, ntt2q, nflevc, nflevp
      real reuvct, reuvpt, tcj, tpj
      common /intdre/ reuvct, reuvpt, tcj, tpj, nnq, nnp, icontr, 
     1   nonlin, nsuper, nn, nns, nqs, nps, neqn1, neqn11, neqn12, 
     2   neqn21, neqn22, np3, np23, np23q, ntt, nttp, ntt2p, ntt1, 
     3   ntt1q, ntt2q, nflevc, nflevp
      call probl(t,y(neqn1),rarr,iarr)
      call doblks(nn,nnq,nnp,y(neqn1),y(neqn11),y(neqn12),
     *            y(neqn21),y(neqn22))
      call foftex(nn,nnq,nnp,y(neqn11),y(neqn12),y(neqn21),y(neqn22),
     *            y(1),ydot,y(neqn1),icontr)
      return
      end
c
      subroutine doblks(nn,nnq,nnp,a,a11,a12,a21,a22)
c
c revision date: July 10 1990
c author: Luca Dieci, School of Mathematics, GaTech
c
      integer nn, nnq, nnp
      real a, a11, a12, a21, a22
      dimension a(nn,nn),a11(nnq,nnq),a12(nnq,nnp),
     *          a21(nnp,nnq),a22(nnp,nnp)
c
c form blocks of matrix here
c
      integer i, j
      do 5 j=1,nnq
         do 5 i=1,nnq
 5          a11(i,j)=a(i,j)
      do 10 j=1,nnp
         do 10 i=1,nnq
 10         a12(i,j)=a(i,nnq+j)
      do 15 j=1,nnq
         do 15 i=1,nnp
 15         a21(i,j)=a(nnq+i,j)
      do 20 j=1,nnp
         do 20 i=1,nnp
 20         a22(i,j)=a(nnq+i,nnq+j)
      return
      end
c
      subroutine foftex(n,nnq,nnp,a11,a12,a21,a22,
     1       xsave,effe,como,icontr)
c
c revision date: July 10 1990
c author: Luca Dieci, School of Mathematics, GaTech
c
      integer n,nnq,nnp,icontr
      real a11, a12, a21, a22, xsave, effe, como
      dimension a11(nnq,nnq),a12(nnq,nnp),a21(nnp,nnq),a22(nnp,nnp),
     *   xsave(nnp,nnq),como(n,n),effe(nnp,nnq)
c here we form F(t,X)=A21+A22*X-X*A11-X*A12*X, where all quantities have
c been evaluated at t beforehand.
c we attempt to economize on forming F(t,X) according to whether
c NP.GT.NQ or not, and if ICONTR=0,1
      integer i, j, k
      real prec
c
      if (nnp.lt.nnq) goto 50
      do 10 j=1,nnq
         do 10 i=1,nnq
            prec=0.e0
            do 5 k=1,nnp
               if (icontr.eq.1) prec=prec+a12(i,k)*xsave(k,j)
               if (icontr.eq.0) prec=prec+a12(i,k)*xsave(k,j)*.5
 5          continue
            prec=prec+a11(i,j)
            como(i,j)=prec
 10   continue
      do 20 j=1,nnq
         do 20 i=1,nnp
            prec=0.e0
            do 15 k=1,nnq
               prec=prec-xsave(i,k)*como(k,j)
 15         continue
            effe(i,j)=prec
 20   continue
      if (icontr.eq.0) goto 35
      do 30 j=1,nnq
         do 30 i=1,nnp
            prec=0.e0
            do 25 k=1,nnp
               prec=prec+a22(i,k)*xsave(k,j)
 25         continue
            effe(i,j)=effe(i,j)+prec+a21(i,j)
 30   continue
      return
 35   continue
      do 40 j=1,nnq
         do 40 i=1,nnp
            como(i,j)=effe(j,i)
 40   continue
      do 45 j=1,nnq
         do 45 i=1,nnp
            effe(i,j)=effe(i,j)+como(i,j)+a21(i,j)
 45   continue
      return
 50   continue
      do 60 j=1,nnp
         do 60 i=1,nnp
            prec=0.e0
            do 55 k=1,nnq
               prec=prec+xsave(i,k)*a12(k,j)
 55         continue
            prec=a22(i,j)-prec
            como(i,j)=prec
 60   continue
      do 70 j=1,nnq
         do 70 i=1,nnp
            prec=0.e0
            prec=a21(i,j)
            do 65 k=1,nnp
               prec=prec+como(i,k)*xsave(k,j)
 65         continue
            effe(i,j)=prec
 70   continue
      do 80 j=1,nnq
         do 80 i=1,nnp
            prec=0.e0
            do 75 k=1,nnq
               prec=prec+xsave(i,k)*a11(k,j)
 75         continue
            effe(i,j)=effe(i,j)-prec
 80   continue
      return
      end
c
      subroutine prepj (y, wm, iwm, probl, jac, rarr, iarr)
c
c revision date: July 10 1990
c author: Luca Dieci, School of Mathematics, GaTech
c
c history: a routine by the same name is in the original lsode of 
c          Alan Hindmarsh (LLL), but this is basically all new 
c
      external probl, jac
      integer iwm, iarr
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer nnq, nnp, icontr, nonlin, nsuper, nn, nns, nqs, nps,
     1   neqn1, neqn11, neqn12, neqn21, neqn22, np3, np23, np23q,
     2   ntt, nttp, ntt2p, ntt1, ntt1q, ntt2q, nflevc, nflevp
      integer i, ier,  j,  lenp, nnqp1, npp1
      real y, wm, rarr
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real reuvct, reuvpt, tcj, tpj
      real con, hl0  
      dimension y(*), wm(*), rarr(*), iwm(*), iarr(*)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /intdre/ reuvct, reuvpt, tcj, tpj, nnq, nnp, icontr, 
     1   nonlin, nsuper, nn, nns, nqs, nps, neqn1, neqn11, neqn12, 
     2   neqn21, neqn22, np3, np23, np23q, ntt, nttp, ntt2p, ntt1, 
     3   ntt1q, ntt2q, nflevc, nflevp
c-----------------------------------------------------------------------
c prepj is called by stode to compute and process the matrix
c p = i - h*el(1)*j , where j refers to the jacobian.
c the jacobian j does not need to be supplied, as it
c is (exactly) computed from the Riccati equation itself, 
c via the routine jac
c
c the computation of the jacobian is separated in two distinct 
c phases: (i) the computation of the blocks of the Sylvester equation
c ~         ~                                                         
c A22*X - X*A11 = C , to solve, and their decomposition , and 
c (ii) the update of these blocks by premultiplication by h*el(1).  
c (i) this phase is done by first forming the original matrix blocks
c (with calls to probl and doblks) and then forming the blocks of
c the transformed matrix as A22-X*A12 and A11+A12*X respectively,
c where X is the approximate solution as computed by the predictor
c formulas. The decomposition is done by our version of the 
c Bartels-Stewart algorithm (with calls to axpxbd or atxxad), 
c so that the transformed blocks are eigendecomposed.  
c (ii) This update is postponed after the decomposition 
c because of numerical accuracy and efficiency.  
c
c the solution of the linear system is administered by subroutine solsy
c with calls to the Sylvester-Lyapunov solver we wrote.
c
c
c in addition to variables described previously, communication
c with prepj uses the following..
c y     = array containing predicted values on entry.
c wm    = real work space for matrices.  on output it contains the
c         decomposed blocks of the jacobian (kept in compact matrix    
c         notation) and the accumulated orthogonal transformations
c         which we have used for the solution of the sylvester 
c         equation and also eigenvalue information which has emerged
c         during the decomposition.
c         storage of matrix elements starts at wm(3).
c iwm   = integer work space containing eigenvalue information 
c         starting at iwm(21). iwm also contains some stability 
c         parameters ml = iwm(1) and mu = iwm(2).
c el0   = el(1) (input).
c ierpj = output error flag,  = 0 if no trouble, .gt. 0 if
c         p matrix found to be singular.
c jcur  = output flag = 1 to indicate that the jacobian matrix
c         is now current.
c this routine also uses the common variables el0, h, tn, uround,
c miter, n, nfe, and nje, and several from the block /intdre/.
c-----------------------------------------------------------------------
c next 2 lines added for the superstability check 
      tcj = tn
      if (jstart .eq. 0) tpj = tcj
      nje = nje + 1
      ierpj = 0
      jcur = 1
      hl0 = h*el0
      lenp=2*(nqs+nps)+3*nn
      do 110 i = 1,lenp
 110    wm(i+2) = 0.0e0
      call probl(tn,y(neqn1),rarr, iarr)
      call doblks(nn,nnq,nnp,y(neqn1),y(neqn11),y(neqn12),
     *            y(neqn21),y(neqn22))
      call jac(nnq,nnp,y,y(neqn11),y(neqn12),y(neqn22),icontr)
      con = -hl0
c we have the "needed" blocks of Sylvester Equation in y(neqn11),y(neqn22)
c now, we first decompose them with calls to the Sylvester (Lyapunov)
c solver, then after we can multiply them by con and get the right
c signs.  Then, add identity (1/2 on each block)
      do 115 i=1,nps
 115     wm(2+i)=y(neqn22+i-1)
      do 120 i=1,nqs
 120     wm(2+2*nps+i)=y(neqn11+i-1)
c here we perform the orthogonal decomposition of the appropriate
c blocks of the Sylvester equation.  We differentiate depending on
c the possibility of symmetric DRE.
      if (icontr.eq.0) goto 275
      call axpxbd(wm(3),wm(np3),nnp,nnp,nnp,wm(np23),wm(np23q),nnq,
     *            nnq,nnq,uround,uround,wm(ntt),wm(nttp),wm(ntt2p),
     *            wm(ntt1),wm(ntt1q),wm(ntt2q),iwm(21),iwm(nnp+21),ier)
      if (ier.eq.4) then
         ierpj=1
         return
      endif
      do 125 i=1,nps
 125     wm(2+i)=wm(2+i)*con
      do 130 i=1,nqs
 130     wm(2+2*nps+i)=-wm(2+2*nps+i)*con
      j=3
      npp1=nnp+1
      do 135 i=1,nnp
         wm(j)=wm(j)+.5e0
 135     j=j+npp1
      j=3+2*nps
      nnqp1=nnq+1
      do 140 i=1,nnq
         wm(j)=wm(j)+.5e0
 140     j=j+nnqp1
      return
 275  continue
c load the correct matrix to be decomposed for the Lyapunov equation
      do 277 i=1,nps
 277     wm(2+i)=wm(2+2*nps+i)
      call atxxad(wm(3),wm(np3),nnp,nnp,nnp,uround,wm(ntt),wm(nttp),
     *            wm(ntt2p),iwm(21),ier)
      if (ier.eq.4) then
         ierpj=1
         return
      endif
c now complete the transformation by multiplying by con, getting the
c right sign and adding the identity (1/2, for us)
      do 280 i=1,nps
 280     wm(2+i)=-wm(2+i)*con
      j=3
      npp1=nnp+1
      do 285 i=1,nnp
         wm(j)=wm(j)+.5e0
 285     j=j+npp1
      return
c----------------------- end of subroutine prepj -----------------------
      end
      subroutine solsy (wm, x)
c
c revision date: July 10 1990
c author: Luca Dieci, School of Mathematics, GaTech
c
c history: a routine by the same name is in the original lsode of 
c          Alan Hindmarsh (LLL), but this is basically all new
c
      real wm, x
      dimension wm(*), x(*)
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer nnq, nnp, icontr, nonlin, nsuper, nn, nns, nqs, nps,
     1   neqn1, neqn11, neqn12, neqn21, neqn22, np3, np23, np23q,
     2   ntt, nttp, ntt2p, ntt1, ntt1q, ntt2q, nflevc, nflevp
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real reuvct, reuvpt, tcj, tpj
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /intdre/ reuvct, reuvpt, tcj, tpj, nnq, nnp, icontr, 
     1   nonlin, nsuper, nn, nns, nqs, nps, neqn1, neqn11, neqn12, 
     2   neqn21, neqn22, np3, np23, np23q, ntt, nttp, ntt2p, ntt1, 
     3   ntt1q, ntt2q, nflevc, nflevp
c-----------------------------------------------------------------------
c this routine manages the solution of the linear system arising from
c a chord or quasi-Newton iteration.  it is called with miter=1.
c it calls the sylvester-lyapunov equation solver to do this
c
c communication with solsy uses the following variables..
c
c wm    = real work space containing the inverse diagonal matrix if
c         miter = 3 and the lu decomposition of the matrix otherwise.
c         storage of matrix elements starts at wm(3).
c x     = the right-hand side vector on input, and the solution vector
c         on output, of length n.
c iersl = output flag (in common).  iersl = 0 if no trouble occurred.
c *** we set iersl=-1 if troubles arose when solving Sylvester equation ***
c this routine also uses some of the common variables of /intdre/.
c-----------------------------------------------------------------------
      iersl = 0
      if (icontr.eq.0) goto 150
      call axpxbs(wm(3),wm(np3),nnp,nnp,nnp,wm(np23),wm(np23q),nnq,
     *            nnq,nnq,x,nnp,wm(ntt),wm(ntt1),iersl)
      if (iersl.eq.4) iersl=-1
      return
 150  call atxxas(wm(3),wm(np3),x,nnp,nnp,nnp,nnp,wm(ntt),iersl)
      if (iersl.eq.4) iersl=-1
      return
c----------------------- end of subroutine solsy -----------------------
      end
c
      subroutine jac(nnq,nnp,x,a011,a012,a022,icontr)
c
c revision date: July 10 1990
c author: Luca Dieci, School of Mathematics, GaTech
c
      integer nnq, nnp, icontr
      real x, a011, a012, a022
      dimension a011(nnq,nnq),a012(nnq,nnp),a022(nnp,nnp),x(nnp,nnq)
c
c here we form the blocks of the Newton iteration (Sylvester eq.n)
c
      integer i, j, k
      real sum
      do 10 j=1,nnq
         do 10 i=1,nnq
              sum=0.e0
            do 5 k=1,nnp
  5            sum=sum+a012(i,k)*x(k,j)
 10         a011(i,j)=a011(i,j)+sum
      if (icontr.eq.0) goto 25
      do 20 j=1,nnp
          do 20 i=1,nnp
            sum=0.e0
            do 15 k=1,nnq
 15            sum=sum+x(i,k)*a012(k,j)
 20         a022(i,j)=a022(i,j)-sum
      return
 25   do 30 j=1,nnp
           do 30 i=1,nnp
 30           a022(i,j)=-a011(j,i)
      return
      end
c
      subroutine axpxbd(a,u,m,na,nu,b,v,n,nb,nv,epsa,
     1                 epsb,orta,era,eia,ortb,erb,eib,itypea,
     2                 itypeb,iflag)
c
c last modification:  august 12, 1989 
c
c history: based on the software by Bartels & Stewart, Algorithm
c          432, Comm ACM 15, 1972.
c
      integer m, na, nu, n, nb, nv, itypea, itypeb, iflag
      real a, u, b, v, epsa, epsb, orta, era, eia, ortb, erb, eib
      dimension itypea(m), itypeb(n),
     1   a(na,m),u(nu,m),b(nb,n),v(nv,n),orta(m),ortb(n),era(m),
     2   eia(m),erb(n),eib(n)
c
c axpxbd is a fortran 77 subroutine to perform the decomposition of a and b
c for the equation ax + xb = c.  the matrices a and b are transformed
c into real schur form.
c the program requires the auxiliary
c subroutines orthes, ortran, hqr3.
c the parameters in the argument list are
c
c         a         a doubly subscripted array containing
c                   the matrix a. On return, the lower triangle
c                   and superdiagonal of the array a contain
c                   a lower real schur form of a.  the array
c                   a must be dimensioned at least m by m.
c         u         a doubly subscripted array that, on
c                   return, contains the orthogonal matrix
c                   that reduces a to real schur form with
c                   eigenvalues ordered from smallest to largest
c                   real parts.
c         m         the order of the matrix a.
c         na        the first dimension of the array a.
c         nu        the first dimension of the array u.
c         b         a doubly subscripted array containing
c                   the matrix b.  on return, the upper triangle
c                   and subdiagonal of the array b contain an
c                   upper real schur form of b.  the array b
c                   must be dimensioned at least n by n.
c         v         a doubly subscripted array that, on
c                   return, contains the orthogonal matrix
c                   that reduces b to real schur form with
c                   eigenvalues ordered from smallest to largest
c                   real parts.
c         n         the order of the matrix b.
c         nb        the first dimension of the array b.
c         nv        the first dimension of the array v.
c         epsa      a convergence criterion for the reduction
c                   of a to schur form.  epsa should be set
c                   to the unit roundoff if a is inputted
c                   with as many significant digits as possible.
c         epsb      a convergence criterrion for the reduction
c                   of b to real schur form.
c         iflag     an error indicator that, on return,
c                   contains an error signal.  if iflag is
c                   equal to 4 then the program was
c                   unable to reduce a (b) to real schur form.
c                   otherwise, iflag is untouched.
c         orta      a singly subscripted real array of dimension m
c                   needed as work space.
c         era       a singly subscripted real array of dimension m
c                   that on return contains the real parts of the
c                   eigenvalues of a.
c         eia       a singly subscripted real array of dimension m
c                   that on return contains the imaginary parts of
c                   the eigenvalues of a.
c         itypea    a singly subscripted integer array of dimension m
c                   whose i-th entry on return is
c                   0  i-th a-e.value is real
c                   1  i-th a-e.value is complex with positive imaginary part
c                   2  i-th a-e.value is complex with negative imaginary part
c                  -1  i-th a-e.value not successfully computed
c         ortb      same as orta, but for matrix b.
c         erb       same as era, but for the matrix b
c         eib       same as eia, but for matrix b
c         itypeb    same as itypea, but for the matrix b
c
c
c when epsa is negative the reduction of a to real schur form
c is skipped and the arrays a and u are assumed to contain the
c schur form and accompanying orthogonal matrix. likewise, if epsb
c is negative, the reduction of b to real schur form is skipped.
c this is ignored for our purposes of dresol.
c
      integer i, j, nlow, nup, mlow, mup
      real temp
c
      iflag=0
      mlow=1
      mup=m
      nlow=1
      nup=n
c
c if required,  reduce a(transpose) to upper real schur form.
c
      if (epsa .lt. 0.) go to 35
      do 10 i=1,m
         do 10 j=i,m
            temp = a(i,j)
            a(i,j) = a(j,i)
            a(j,i) = temp
   10 continue
      call orthes(na,m,mlow,mup,a,orta)
      call ortran(na,m,mlow,mup,a,orta,u)
      call hqr3(a,u,m,mlow,mup,epsa,era,eia,itypea,na,nu)
      do 20 i=mlow,mup
         if (itypea(i).lt.0) then
            iflag=4
            return
         endif
 20   continue
      do 30 i=1,m
         do 30 j=i,m
            temp = a(i,j)
            a(i,j) = a(j,i)
            a(j,i) = temp
   30 continue
c
c if required,  reduce b to upper real schur form.
c
   35 if (epsb .lt. 0.) goto 45
      call orthes(nb,n,nlow,nup,b,ortb)
      call ortran(nb,n,nlow,nup,b,ortb,v)
      call hqr3(b,v,n,nlow,nup,epsb,erb,eib,itypeb,nb,nv)
      do 40 i=nlow,nup
         if (itypeb(i).lt.0) then
            iflag=4
            return
         endif
 40   continue
 45   return
c end of axpxbd
      end
c
      subroutine axpxbs(a,u,m,na,nu,b,v,n,nb,nv,c,nc,orta,ortb,iflag)
c
c last modification:  august 12, 1989 
c
c history: based on the software by Bartels & Stewart.
c
c axpxbs is a fortran 77 subroutine to solve the sylvester
c equation ax + xb = c.  the matrices a and b are already in
c real schur form.
c the program requires the auxiliary
c subroutines shrslv and sysslv.
c the parameters in the argument list are
c
c         a         a doubly subscripted array as follow.
c                   the lower triangle
c                   and superdiagonal of the array a contain
c                   a lower real schur form of a.  the array
c                   a must be dimensioned at least m by m.
c         u         a doubly subscripted array that,
c                   contains the orthogonal matrix
c                   that has reduced a to real schur form with
c                   eigenvalues ordered from smallest to
c                   largest real parts.
c         m         the order of the matrix a.
c         na        the first dimension of the array a.
c         nu        the first dimension of the array u.
c         b         a doubly subscripted array as follow.
c                   the upper triangle
c                   and subdiagonal of the array b contain an
c                   upper real schur form of b.  the array b
c                   must be dimensioned at least n by n.
c         v         a doubly subscripted array that
c                   contains the orthogonal matrix
c                   that has reduced b to real schur form with
c                   eigenvalues ordered from smallest to
c                   largest real parts.
c         n         the order of the matrix b.
c         nb        the first dimension of the array b.
c         nv        the first dimension of the array v.
c         c         a doubly subscripted array that contains
c                   the right-hand-side matrix c.  on return
c                   it contains the solution matrix x.
c         nc        the first dimension of c.
c         iflag     an error indicator that, on return,
c                   contains an error signal.  if iflag is
c                   equal to 4 then the program detected a division
c                   by zero.  this means sylvester equation not
c                   solvable.  otherwise, iflag is set to 0. 
c         orta      a singly subscripted real array of dimension m
c                   needed as work space.
c         ortb      same as orta, but of dimension n.
c
c
      integer
     *na,nu,nb,nv,nc,m,n,i,j,k,iflag
      real
     *a(na,m),u(nu,m),b(nb,n),v(nv,n),orta(m),ortb(n),c(nc,n)
      iflag=0
c transform c.
      do 20 j=1,n
         do 10 i=1,m
            orta(i)=0.
            do 10 k=1,m
               orta(i) = orta(i) + u(k,i)*c(k,j)
   10 continue
      do 20 i=1,m
         c(i,j) = orta(i)
   20 continue
      do 40 i=1,m
         do 30 j=1,n
            ortb(j) = 0.
            do 30 k=1,n
               ortb(j) = ortb(j) + c(i,k)*v(k,j)
   30    continue
         do 40 j=1,n
            c(i,j) = ortb(j)
   40 continue
c
c solve the transformed system.
c
      call shrslv(a,b,c,m,n,na,nb,nc,iflag)
      if (iflag.eq.4) return
c
c transform c back to the solution.
c
      do 100 j=1,n
         do 90 i=1,m
            orta(i) = 0.
            do 90 k=1,m
               orta(i) = orta(i) + u(i,k)*c(k,j)
   90 continue
      do 100 i=1,m
         c(i,j) = orta(i)
  100 continue
      do 120 i=1,m
         do 110 j=1,n
            ortb(j) = 0.
            do 110 k=1,n
               ortb(j) = ortb(j) + c(i,k)*v(j,k)
  110 continue
      do 120 j=1,n
         c(i,j) = ortb(j)
  120 continue
      return
      end
c
      subroutine shrslv(a,b,c,m,n,na,nb,nc,iflag)
c
c last modification:  august 12, 1989 
c
c history: based on the software by Bartels & Stewart.
c
c shrslv is a fortran iv subroutine to solve the real matrix
c equation ax + bx = c,where a is in lower real schur form
c and b is in upper real schur form.  shrslv uses the aux-
c iliary subroutine sysslv,which it communicates with
c through the common block slvblk.  the parameters in the
c argument list are
c           a      a doubly subscripted array containing the
c                  matrix a in lower real schur form.
c           b      a doubly subscripted array containing the
c                  matrix b in upper real schur form.
c           c      a doubly subscripted array containing rhe
c                  matrix c.
c           m      the order of the matrix a.
c           n      the order of the matrix b.
c           na     the first demension of the array a.
c           nb     the first dimension of the array b.
c           nc     the first dimension of the array c.
c
      integer
     1m,n,na,nb,nc,k,km1,dk,kk,l,lm1,dl,ll,i,ib,j,ja,nsys,iflag
      real
     1a(na,m),b(nb,n),c(nc,n),t,p
      common/slvblk/t(5,5),p(5),nsys
      l = 1
   10   lm1 = l-1
        dl = 1
        if (l .eq. n) go to 15
        if (b(l+1,l) .ne. 0.) dl = 2
   15   ll = l+dl-1
        if (l .eq. 1) go to 30
        do 20 j=l,ll
           do 20 i=1,m
              do 20 ib=1,lm1
                 c(i,j) = c(i,j) - c(i,ib)*b(ib,j)
   20   continue
   30   k = 1
   40     km1 = k-1
          dk = 1
          if (k .eq. m) go to 45
          if (a(k,k+1) .ne. 0.) dk = 2
   45     kk = k+dk-1
          if (k .eq. 1) go to 60
          do 50 i=k,kk
             do 50 j=l,ll
                do 50 ja=1,km1
                   c(i,j) = c(i,j) - a(i,ja)*c(ja,j)
   50     continue
   60     if (dl .eq. 2) go to 80
          if (dk .eq. 2) go to 70
          t(1,1) = a(k,k) + b(l,l)
          if (t(1,1) .eq. 0.) then
             iflag=4
             return
          endif
          c(k,l) = c(k,l)/t(1,1)
          go to 100
   70     t(1,1) = a(k,k) + b(l,l)
          t(1,2) = a(k,kk)
          t(2,1) = a(kk,k)
          t(2,2) = a(kk,kk) + b(l,l)
          p(1) = c(k,l)
          p(2) = c(kk,l)
          nsys = 2
          call sysslv
          c(k,l) = p(1)
          c(kk,l) = p(2)
          go to 100
   80     if (dk .eq. 2) go to 90
          t(1,1) = a(k,k) + b(l,l)
          t(1,2) = b(ll,l)
          t(2,1) = b(l,ll)
          t(2,2) = a(k,k) + b(ll,ll)
          p(1) = c(k,l)
          p(2) = c(k,ll)
          nsys = 2
          call sysslv
          c(k,l) = p(1)
          c(k,ll) = p(2)
          go to 100
   90     t(1,1) = a(k,k) + b(l,l)
          t(1,2) = a(k,kk)
          t(1,3) = b(ll,l)
          t(1,4) = 0.
          t(2,1) = a(kk,k)
          t(2,2) = a(kk,kk) + b(l,l)
          t(2,3) = 0.
          t(2,4) = t(1,3)
          t(3,1) = b(l,ll)
          t(3,2) = 0.
          t(3,3) = a(k,k) + b(ll,ll)
          t(3,4) = t(1,2)
          t(4,1) = 0.
          t(4,2) = t(3,1)
          t(4,3) = t(2,1)
          t(4,4) = a(kk,kk) + b(ll,ll)
          p(1) = c(k,l)
          p(2) = c(kk,l)
          p(3) = c(k,ll)
          p(4) = c(kk,ll)
          nsys = 4
          call sysslv
          c(k,l) = p(1)
          c(kk,l) = p(2)
          c(k,ll) = p(3)
          c(kk,ll) = p(4)
  100   k = k + dk
        if (k .le. m) go to 40
      l = l + dl
      if (l .le. n) go to 10
      return
      end
c
c
      subroutine atxxad(a,u,n,na,nu,eps,ort,er,ei,itype,iflag)
c
c last modification:  august 12, 1989 
c
c history: based on the software by Bartels & Stewart.
c
c
c atxxad is a fortran 77 subroutine to decompose a prior to solving the
c real lyapunov equation trans(a)*x + x*a = c, where c is symmetric and
c trans(a) denotes the transpose of a. the equation is
c transformed so that a is in upper real schur form.
c the program requires the auxiliary subroutines orthes,
c ortan, hqr3 and related. the parameters in the argument
c list are
c           a     a doubly subscripted array containing the
c                 matrix a. on return, the upper triangle
c                 and the first subdiagonal of the array a
c                 contain an upper real schur form of a.
c                 the array a must be dimensioned ar least
c                 n by n.
c           u     a doubly subscripted array that, on
c                 return, contains the orthogonal matrix
c                 that reduces a to upper real schur form with
c                 eigenvalues ordered in descending absolute value.
c           n     the order of the matrix a.
c           na    the first dimension of the the array a.
c           nu    the first dimension of the array u.
c           eps   a convergence criterion for the reduction
c                 of a to real schur form. eps should be
c                 set to the machine precision.
c                 when a is correct within the precision itself.
c           ort   singly subscripted real array of dimension n
c                 needed as work space
c           er    singly subscripted real array of dimension n
c                 that on return contains the real parts of the e.values
c           ei    singly subscripted real array that on return
c                 contains the imaginary parts of the e.values
c         itype   integer array whose i-th entry on return is
c                 0  i-th e.value is real
c                 1  i-th e.value is complex with positive imaginary part
c                 2  i-th e.value is complex with negative imaginary part
c                -1  i-th e.value not successfully computed
c          iflag  an error indicator.  if iflag is set to 4, then
c                 program was unable to reduce to a to real schur
c                 form. it is set to 0 otherwise.
c
c if eps is negative, the reduction of a to real schur
c form is skipped and the arrays a and u are assumed to
c contain the schur form and accompanying orthogonal matrix.
c this makes sense just to check the input.
c
      integer
     *i,n,na,nu,itype(n),nlow,nup,iflag
      real
     *a(na,n),u(nu,n),eps,ort(n),er(n),ei(n)
      iflag=0
      nlow=1
      nup=n
c
c if required,  reduce a to upper real schur form.
c
      if (eps .lt. 0.) go to 15
      call orthes(na,n,nlow,nup,a,ort)
      call ortran(na,n,nlow,nup,a,ort,u)
      call hqr3(a,u,n,nlow,nup,eps,er,ei,itype,na,nu)
      do 10 i=nlow,nup
         if (itype(i).lt.0)then
            iflag=4
            return
         endif
 10   continue
 15   return
c end of atxxad
      end
c
      subroutine atxxas(a,u,c,n,na,nu,nc,ort,iflag)
c
c last modification:  august 12, 1989 
c
c history: based on the software by Bartels & Stewart.
c
c atxxas is a fortran iv subroutine to solve the real matrix
c equation trans(a)*x + x*a = c, where c is symmetric and
c trans(a) denotes the transpose of a.
c a is in upper real schur form, and the
c transformed equation is solved by a recursive procedure.
c the program requires the auxiliary subroutines symslv and sysslv
c
c parameters in the calling list are>
c           a     a doubly subscripted array as follows.
c                 the upper triangle
c                 and the first subdiagonal of the array a
c                 contain an upper real schur form of a.
c                 the array a must be dimensioned ar least
c                 n by n.
c           u     a doubly subscripted array that
c                 contains the orthogonal matrix
c                 that has reduced a to upper real schur form with
c                 eigenvalues ordered in descending absolute value.
c           c     a doubly subscripted array containing
c                 the matrix c. on return, c contains
c                 the solution matrix x.
c           n     the order of the matrix a.
c           na    the first dimension of the the array a.
c           nu    the first dimension of the array u.
c           nc    the first dimension of the array c.
c           ort   singly subscripted real array of dimension n
c                 needed as work space
c          iflag  an error indicator.  if iflag is set to 4, then
c                 program detected a zero eigenvalue of a, thus
c                 cannot solve the lyapunov equation. set to 0 otherwise.
c
      integer
     *i,j,k,n,na,nu,nc,iflag
      real
     *a(na,n),u(nu,n),c(nc,n),ort(n)
      iflag=0
c transform c.
      do 20 i=1,n
         c(i,i) = c(i,i)/2.
   20 continue
      do 40 i=1,n
         do 30 j=1,n
            ort(j)=0.
            do 30 k=i,n
               ort(j) = ort(j) + c(i,k)*u(k,j)
   30 continue
      do 40 j=1,n
         c(i,j) = ort(j)
   40 continue
      do 60 j=1,n
         do 50 i=1,n
            ort(i) = 0.
            do 50 k=1,n
               ort(i) = ort(i) + u(k,i)*c(k,j)
   50 continue
      do 60 i=1,n
         c(i,j) = ort(i)
   60 continue
      do 70 i=1,n
         do 70 j=i,n
            c(i,j)=c(i,j)+c(j,i)
            c(j,i)=c(i,j)
   70 continue
c
c solve the transformed system.
c
      call symslv(a,c,n,na,nc,iflag)
      if (iflag.eq.4) return
c
c transform c back to the solution.
c
      do 80 i=1,n
         c(i,i)=c(i,i)/2.
   80 continue
         do 100 i=1,n
            do 90 j=1,n
            ort(j) = 0.
            do 90 k=i,n
               ort(j) = ort(j) + c(i,k)*u(j,k)
   90 continue
      do 100 j=1,n
         c(i,j) = ort(j)
  100 continue
      do 120 j=1,n
         do 110 i=1,n
            ort(i) = 0.
            do 110 k=1,n
               ort(i) = ort(i) + u(i,k)*c(k,j)
  110 continue
      do 120 i=1,n
         c(i,j) = ort(i)
  120 continue
      do 130 i=1,n
         do 130 j=i,n
            c(i,j)=c(i,j)+c(j,i)
            c(j,i)=c(i,j)
  130 continue
      return
      end
c
c
      subroutine symslv(a,c,n,na,nc,iflag)
c
c last modification:  august 12, 1989 
c
c history: based on the software by Bartels & Stewart.
c
c
c symslv is fortran 77 subroutine to solve the real matrix
c equation trans(a)*x + xa = c, where c is symmetric, a is
c in upper real schur form, and trans(a) denotes the transpose
c of a. symslv uses the auxiliary subroutine sysslv, which it
c communicates with through the common block slvblk. the parameters
c in the argument list are
c        a      a doubly subscripted array containing the matrix a
c               in upper real schur form.
c        c      a doubly subscripted array containing the matrix c.
c        n      the order of the matrix a.
c        na     the first dimension of the array a.
c        nc     the first dimension of the array c.
c
      integer
     *n,na,nc,k,km1,dk,kk,l,ldl,dl,ll,i,ia,j,nsys,iflag
      real
     *a(na,n),c(nc,n),t,p
      common/slvblk/t(5,5),p(5),nsys
      l = 1
   10   dl = 1
        if (l .eq. n) go to 20
        if (a(l+1,l) .ne. 0.) dl = 2
   20   ll = l+dl-1
        k = l
   30     km1 = k-1
          dk = 1
          if (k .eq. n) go to 35
          if (a(k+1,k) .ne. 0.) dk = 2
   35     kk = k+dk-1
          if (k .eq. l) go to 45
          do 40 i=k,kk
             do 40 j=l,ll
                do 40 ia=l,km1
                   c(i,j) = c(i,j) - a(ia,i)*c(ia,j)
   40     continue
   45     if (dl .eq. 2) go to 60
          if (dk .eq. 2) go to 50
          t(1,1) = a(k,k) + a(l,l)
          if (t(1,1) .eq. 0.) then
             iflag=4
             return
          endif
          c(k,l) = c(k,l)/t(1,1)
          go to 90
   50     t(1,1) = a(k,k) + a(l,l)
          t(1,2) = a(kk,k)
          t(2,1) = a(k,kk)
          t(2,2) = a(kk,kk) + a(l,l)
          p(1) = c(k,l)
          p(2) = c(kk,l)
          nsys = 2
          call sysslv
          c(k,l) = p(1)
          c(kk,l) = p(2)
          go to 90
   60     if (dk .eq. 2) go to 70
          t(1,1) = a(k,k) + a(l,l)
          t(1,2) = a(ll,l)
          t(2,1) = a(l,ll)
          t(2,2) = a(k,k) + a(ll,ll)
          p(1) = c(k,l)
          p(2) = c(k,ll)
          nsys = 2
          call sysslv
          c(k,l) = p(1)
          c(k,ll) = p(2)
          go to 90
   70     if (k .ne. l) go to 80
          t(1,1) = a(l,l)
          t(1,2) = a(ll,l)
          t(1,3) = 0.
          t(2,1) = a(l,ll)
          t(2,2) = a(l,l)+a(ll,ll)
          t(2,3) = t(1,2)
          t(3,1) = 0.
          t(3,2) = t(2,1)
          t(3,3) = a(ll,ll)
          p(1) = c(l,l)/2.
          p(2) = c(ll,l)
          p(3) = c(ll,ll)/2.
          nsys = 3
          call sysslv
          c(l,l) = p(1)
          c(ll,l) = p(2)
          c(l,ll) = p(2)
          c(ll,ll) = p(3)
          go to 90
   80     t(1,1) = a(k,k) + a(l,l)
          t(1,2) = a(kk,k)
          t(1,3) = a(ll,l)
          t(1,4) = 0.
          t(2,1) = a(k,kk)
          t(2,2) = a(kk,kk) + a(l,l)
          t(2,3) = 0.
          t(2,4) = t(1,3)
          t(3,1) = a(l,ll)
          t(3,2) = 0.
          t(3,3) = a(k,k) + a(ll,ll)
          t(3,4) = t(1,2)
          t(4,1) = 0.
          t(4,2) = t(3,1)
          t(4,3) = t(2,1)
          t(4,4) = a(kk,kk) + a(ll,ll)
          p(1) = c(k,l)
          p(2) = c(kk,l)
          p(3) = c(k,ll)
          p(4) = c(kk,ll)
          nsys = 4
          call sysslv
          c(k,l) = p(1)
          c(kk,l) = p(2)
          c(k,ll) = p(3)
          c(kk,ll) = p(4)
   90   k = k + dk
        if (k .le. n) go to 30
        ldl = l + dl
        if (ldl .gt. n) return
        do 120 j=ldl,n
           do 100 i=l,ll
              c(i,j)=c(j,i)
  100      continue
           do 120 i=j,n
              do 110 k=l,ll
                 c(i,j)=c(i,j) - c(i,k)*a(k,j) - a(k,i)*c(k,j)
  110         continue
              c(j,i)=c(i,j)
  120      continue
      l=ldl
      go to 10
      end
c
      subroutine sysslv
c
c last modification:  august 12, 1989 
c
c history: taken from the software by Bartels & Stewart.
c
c sysslv is a fortran 77 subroutine that solves the linear
c system ax = b of order n less than 5 by crout reduction
c followed by back substitution.  the matrix a,  the vector
c b,  and the order n are contained on the arrays a, b, and
c the variable n of the common block slvblk.  the solution
c is returned on the array b.
c
      integer
     *n,nm1,n1,k,km1,i,j,kp1,intr,im1,i1,ii
      real
     *rmax,a,b,aa,temp
      common/slvblk/a(5,5),b(5),n
      nm1 = n-1
      n1 = n+1
c
c compute the lu factorization of a.
c
      do 80 k=1,n
         km1 = k-1
         if (k .eq. 1) go to 20
         do 10 i=k,n
            do 10 j=1,km1
               a(i,k) = a(i,k) - a(i,j)*a(j,k)
   10    continue
   20    if (k .eq. n) go to 100
         kp1 = k+1
         rmax =abs(a(k,k))
         intr = k
         do 30 i=kp1,n
            aa = abs(a(i,k))
            if (aa .le. rmax) go to 30
            rmax = aa
            intr = i
   30    continue
         if (rmax .eq. 0.) stop
         a(n1,k) = intr
         if (intr .eq. k) go to 50
         do 40 j=1,n
            temp = a(k,j)
            a(k,j)  = a(intr,j)
            a(intr,j) = temp
   40    continue
   50    do 80 j=kp1,n
            if (k .eq. 1) go to 70
            do 60 i=1,km1
               a(k,j) = a(k,j) - a(k,i)*a(i,j)
   60       continue
   70       a(k,j)  = a(k,j)/a(k,k)
   80 continue
c
c interchange the components of b.
c
  100 do 110 j=1,nm1
         intr = a(n1,j)
         if (intr .eq. j) go to 110
         temp = b(j)
         b(j) = b(intr)
         b(intr) = temp
  110 continue
c
c solve lx = b.
c
  200 b(1) = b(1)/a(1,1)
      do 220 i=2,n
         im1 = i-1
         do 210 j=1,im1
            b(i) = b(i) - a(i,j)*b(j)
  210    continue
         b(i) = b(i)/a(i,i)
  220 continue
c
c solve ux = b.
c
  300 do 310 ii=1,nm1
         i = nm1-ii+1
         i1 = i+1
         do 310 j=i1,n
            b(i) = b(i) - a(i,j)*b(j)
  310 continue
      return
      end
c
c
      subroutine orthes(nm,n,low,igh,a,ort)
c 
c this is from the EISPACK collection
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      real a(nm,n),ort(igh)
      real f,g,h,scale
c *** functions
      real sqrt,abs,sign
c
c     this subroutine is a translation of the algol procedure orthes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a real general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     orthogonal similarity transformations.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrix,
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n,
c
c        a contains the input matrix.
c
c     on output-
c
c        a contains the hessenberg matrix.  information about
c          the orthogonal transformations used in the reduction
c          is stored in the remaining triangle under the
c          hessenberg matrix,
c
c        ort contains further information about the transformations.
c          only elements low through igh are used.
c
c     questions and comments should be directed to B. S. Garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0
         ort(m) = 0.0
         scale = 0.0
c     ********** scale column (algol tol then not needed) **********
         do 90 i = m, igh
   90    scale = scale + abs(a(i,m-1))
c
         if (scale .eq. 0.0) go to 180
         mp = m + igh
c     ********** for i=igh step -1 until m do -- **********
         do 100 ii = m, igh
            i = mp - ii
            ort(i) = a(i,m-1) / scale
            h = h + ort(i) * ort(i)
  100    continue
c
         g = -sign(sqrt(h),ort(m))
         h = h - ort(m) * g
         ort(m) = ort(m) - g
c     ********** form (i-(u*ut)/h) * a **********
         do 130 j = m, n
            f = 0.0
c     ********** for i=igh step -1 until m do -- **********
            do 110 ii = m, igh
               i = mp - ii
               f = f + ort(i) * a(i,j)
  110       continue
c
            f = f / h
c
            do 120 i = m, igh
  120       a(i,j) = a(i,j) - f * ort(i)
c
  130    continue
c     ********** form (i-(u*ut)/h)*a*(i-(u*ut)/h) **********
         do 160 i = 1, igh
            f = 0.0
c     ********** for j=igh step -1 until m do -- **********
            do 140 jj = m, igh
               j = mp - jj
               f = f + ort(j) * a(i,j)
  140       continue
c
            f = f / h
c
            do 150 j = m, igh
  150       a(i,j) = a(i,j) - f * ort(j)
c
  160    continue
c
         ort(m) = scale * ort(m)
         a(m,m-1) = scale * g
  180 continue
c
  200 return
c     ********** last card of orthes **********
      end
c
      subroutine ortran(nm,n,low,igh,a,ort,z)
c 
c this is from the EISPACK collection
c
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      real a(nm,igh),ort(igh),z(nm,n)
      real g
c
c     this subroutine is a translation of the algol procedure ortrans,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine accumulates the orthogonal similarity
c     transformations used in the reduction of a real general
c     matrix to upper hessenberg form by  orthes.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrix,
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n,
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction by  orthes
c          in its strict lower triangle,
c
c        ort contains further information about the trans-
c          formations used in the reduction by  orthes.
c          only elements low through igh are used.
c
c     on output-
c
c        z contains the transformation matrix produced in the
c          reduction by  orthes,
c
c        ort has been altered.
c
c     questions and comments should be directed to B. S. Garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c     ********** initialize z to identity matrix **********
      do 80 i = 1, n
c
         do 60 j = 1, n
   60    z(i,j) = 0.0
c
         z(i,i) = 1.0
   80 continue
c
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
c     ********** for mp=igh-1 step -1 until low+1 do -- **********
      do 140 mm = 1, kl
         mp = igh - mm
         if (a(mp,mp-1) .eq. 0.0) go to 140
         mp1 = mp + 1
c
         do 100 i = mp1, igh
  100    ort(i) = a(i,mp-1)
c
         do 130 j = mp, igh
            g = 0.0
c
            do 110 i = mp, igh
  110       g = g + ort(i) * z(i,j)
c     ********** divisor below is negative of h formed in orthes.
c                double division avoids possible underflow **********
            g = (g / ort(mp)) / a(mp,mp-1)
c
            do 120 i = mp, igh
  120       z(i,j) = z(i,j) + g * ort(i)
c
  130    continue
c
  140 continue
c
  200 return
c     ********** last card of ortran **********
      end
c
c below, the routines hqr3, exchng, qrstep and split 
c are all from the HQR3 software of G.W.Stewart.
c minor modifications have been made to this software, i.e. 
c translating it to single precision and little else 
c
      subroutine hqr3 (a,v,n,nlow,nup,eps,er,ei,itype,na,nv)
c
c this routine is from the HQR3 software by G.W.Stewart
c
c modified on march 1988 by luca dieci
c so that it performs the reordering from
c smallest to largest real parts of the eigenvalues.
c
c     *****parameters:
      integer n,na,nlow,nup,nv,itype(n)
      real a(na,n),ei(n),er(n),eps,v(nv,n)
c
c     *****local variables:
      logical fail
      integer i,it,l,mu,nl,nu
      real e1,e2,p,q,r,s,t,w,x,y,z
c
c     *****functions:
      real abs
c
c     *****subroutines called:
c     exchng,qrstep,split
c
c     ------------------------------------------------------------------
c
c     *****purpose:
c     this subroutine reduces the upper hessenberg matrix a to quasi-
c     triangular form by unitary similarity transformations. the
c     eigenvalues of a, which are contained in the 1 x 1 and 2 x 2
c     diagonal blocks of the reduced matrix, are ordered in descending
c     order of magnitude along the diagonal. the transformations are
c     accumulated in the array v.
c
c     *****parameter description:
c     on input:
c        na,nv            row dimensions of the arrays containing a and
c                         v, respectively, as declared in the calling
c                         program dimension statement;
c
c        a                n x n array containing the upper hessenberg
c                         matrix to be reduced;
c
c        n                order of the matrices a and v;
c
c        nlow,nup         a(nlow,nlow-1) and a(nup,1+nup) are assumed
c                         to be zero, and only rows nlow through nup
c                         and columns nlow through nup are transformed,
c                         resulting in the calculation of eigenvalues
c                         nlow through nup;
c
c        eps              a convergence criterion used to determine when
c                         a subdiagonal element of a is negligible.
c                         specifically, a(i+1,i) is regarded as
c                         negligible if abs(a(i+1),i)).le.eps*
c                         (abs(a(i+1,i+1))). this means that the final
c                         matrix returned by the program will be exactly
c                         similar to a + e where e is of order
c                         eps*norm(a), for any reasonably balanced norm
c                         such as the row-sum norm;
c
c        itype            an integer vector of length n whose
c                         i-th entry is
c                         0  if the i-th eigenvalue is real,
c                         1  if the i-th eigenvalue is complex with
c                            positive imaginary part,
c                         2  if the i-th eigenvalue is complex with
c                            negative imaginary part,
c                         -1 if the i-th eigenvalue was not calculated
c                            successfully.
c
c        on output:
c
c        a                n x n array containing the reduced, quasi-
c                         triangular matrix;
c
c        v                n x n array containing the reducing
c                         transformations to be multiplied;
c
c        er,ei            real scratch vectors of length n which on
c                         return contain the real and imaginary parts,
c                         respectively, of the eigenvalues.
c
c     ------------------------------------------------------------------
c
      do 10 i=nlow,nup
         itype(i)=-1
10    continue
      t=0.0
      nu=nup
20    if (nu.lt.nlow) go to 240
      it=0
30    continue
      l=nu
40    continue
      if (l.eq.nlow) go to 50
      if (abs(a(l,l-1)).le.eps*(abs(a(l-1,l-1))+abs(a(l,l))))
     +    go to 50
      l=l-1
      go to 40
50    continue
      x=a(nu,nu)
      if (l.eq.nu) go to 160
      y=a(nu-1,nu-1)
      w=a(nu,nu-1)*a(nu-1,nu)
      if (l.eq.nu-1) go to 100
      if (it.eq.30) go to 240
      if (it.ne.10 .and. it.ne.20) go to 70
      t=t+x
      do 60 i=nlow,nu
         a(i,i)=a(i,i)-x
60    continue
      s=abs(a(nu,nu-1))+abs(a(nu-1,nu-2))
      x=0.75*s
      y=x
      w=-0.4375*s**2
70    continue
      it=it+1
      nl=nu-2
80    continue
      z=a(nl,nl)
      r=x-z
      s=y-z
      p=(r*s-w)/a(nl+1,nl)+a(nl,nl+1)
      q=a(nl+1,nl+1)-z-r-s
      r=a(nl+2,nl+1)
      s=abs(p)+abs(q)+abs(r)
      p=p/s
      q=q/s
      r=r/s
      if (nl.eq.l) go to 90
      if (abs(a(nl,nl-1))*(abs(q)+abs(r)).le.eps*abs(p)*
     +   (abs(a(nl-1,nl-1))+abs(z)+abs(a(nl+1,nl+1)))) go to 90
      nl=nl-1
      go to 80
90    continue
      call qrstep (a,v,p,q,r,nl,nu,n,na,nv)
      go to 30
100   if (nu.ne.nlow+1) a(nu-1,nu-2)=0.0
      a(nu,nu)=a(nu,nu)+t
      a(nu-1,nu-1)=a(nu-1,nu-1)+t
      itype(nu)=0
      itype(nu-1)=0
      mu=nu
110   continue
      nl=mu-1
      call split (a,v,n,nl,e1,e2,na,nv)
      if (a(mu,mu-1).eq.0.0) go to 170
      if (mu.eq.nup) go to 230
      if (mu.eq.nup-1) go to 130
      if (a(mu+2,mu+1).eq.0.0) go to 130
      if (a(mu-1,mu-1)+a(mu,mu).le.a(mu+1,mu+1)+
     +   a(mu+2,mu+2)) go to 230
      call exchng (a,v,n,nl,2,2,eps,fail,na,nv)
      if (.not.fail) go to 120
      itype(nl)=-1
      itype(nl+1)=-1
      itype(nl+2)=-1
      itype(nl+3)=-1
      go to 240
120   continue
      mu=mu+2
      go to 150
130   continue
      if (a(mu-1,mu-1)+a(mu,mu).le.
     +   2.*a(mu+1,mu+1)) go to 230
      call exchng (a,v,n,nl,2,1,eps,fail,na,nv)
      if (.not.fail) go to 140
      itype(nl)=-1
      itype(nl+1)=-1
      itype(nl+2)=-1
      go to 240
140   continue
      mu=mu+1
150   continue
      go to 110
160   nl=0
      a(nu,nu)=a(nu,nu)+t
      if (nu.ne.nlow) a(nu,nu-1)=0.0
      itype(nu)=0
      mu=nu
170   continue
180   continue
      if (mu.eq.nup) go to 220
      if (mu.eq.nup-1) go to 200
      if (a(mu+2,mu+1).eq.0.0) go to 200
      if (2.*a(mu,mu).le.a(mu+1,mu+1)+a(mu+2,mu+2))
     +    go to 220
      call exchng (a,v,n,mu,1,2,eps,fail,na,nv)
      if (.not.fail) go to 190
      itype(mu)=-1
      itype(mu+1)=-1
      itype(mu+2)=-1
      go to 240
190   continue
      mu=mu+2
      go to 210
200   continue
      if (a(mu,mu).le.a(mu+1,mu+1)) go to 220
      call exchng (a,v,n,mu,1,1,eps,fail,na,nv)
      mu=mu+1
210   continue
      go to 180
220   continue
      mu=nl
      nl=0
      if (mu.ne.0) go to 170
230   continue
      nu=l-1
      go to 20
240   if (nu.lt.nlow) go to 260
      do 250 i=nlow,nu
         a(i,i)=a(i,i)+t
250   continue
260   continue
      nu=nup
270   continue
      if (itype(nu).ne.(-1)) go to 280
      nu=nu-1
      go to 310
280   continue
      if (nu.eq.nlow) go to 290
      if (a(nu,nu-1).eq.0.0) go to 290
      call split (a,v,n,nu-1,e1,e2,na,nv)
      if (a(nu,nu-1).eq.0.0) go to 290
      er(nu)=e1
      ei(nu-1)=e2
      er(nu-1)=er(nu)
      ei(nu)=-ei(nu-1)
      itype(nu-1)=1
      itype(nu)=2
      nu=nu-2
      go to 300
290   continue
      er(nu)=a(nu,nu)
      ei(nu)=0.0
      nu=nu-1
300   continue
310   continue
      if (nu.ge.nlow) go to 270
      return
c
c     last line of hqr3
c
      end
      subroutine exchng (a,v,n,l,b1,b2,eps,fail,na,nv)
c
c
c     based on stewart's hqr3 software.
c     (hammarling's correction added by ajl 01/23/86)
c
c     *****parameters
      integer n,b1,b2,l,na,nv
      real a(na,n),eps,v(nv,n)
      logical fail
c
c     *****local variables:
      integer i,it,j,l1,m
      real p,q,r,s,w,x,y,z
c
c     *****functions:
      real abs,sqrt,amax1
c
c     *****subroutines called:
c     qrstep
c
c     ------------------------------------------------------------------
c
c     *****purpose:
c     given the upper schur matrix a with consecutive b1 x b1 and
c     b2 x b2 diagonal blocks (b1, b2.le.2) starting at a(l,l), this
c     subroutine produces a unitary similarity transformation that
c     exchanges the blocks along with their eigenvalues. the
c     transformation is accumulated in v.
c
c     *****parameter description:
c     on input:
c        na,nv            row dimensions of the arrays containing a
c                         and v, respectively, as declared in the
c                         calling program dimension statement;
c
c        a                n x n matrix whose blocks are to be
c                         interchanged;
c
c        n                order of the matrix a;
c
c        l                position of the blocks;
c
c        b1               an integer containing the size of the first
c                         block;
c
c        b2               an integer containing the size of the second
c                         block;
c
c        eps              a convergence criterion (cf. hqr3).
c
c     on output:
c
c        fail             a logical variable which is .false. on a
c                         normal return. if thirty iterations were
c                         performed without convergence, fail is set to
c                         .true. and the element a(l+b2,l+b2-1) cannot
c                         be assumed zero.
c
c     ------------------------------------------------------------------
c
      fail=.false.
      if (b1.eq.2) go to 70
      if (b2.eq.2) go to 40
      l1=l+1
      q=a(l+1,l+1)-a(l,l)
      if (q .eq. 0.0) return
      p=a(l,l+1)
      r=amax1(abs(p),abs(q))
      if (r.eq.0.0) r = q
      p=p/r
      q=q/r
      r=sqrt(p**2+q**2)
      p=p/r
      q=q/r
      do 10 j=l,n
         s=p*a(l,j)+q*a(l+1,j)
         a(l+1,j)=p*a(l+1,j)-q*a(l,j)
         a(l,j)=s
10    continue
      do 20 i=1,l1
         s=p*a(i,l)+q*a(i,l+1)
         a(i,l+1)=p*a(i,l+1)-q*a(i,l)
         a(i,l)=s
20    continue
      do 30 i=1,n
         s=p*v(i,l)+q*v(i,l+1)
         v(i,l+1)=p*v(i,l+1)-q*v(i,l)
         v(i,l)=s
30    continue
      a(l+1,l)=0.0
      return
40    continue
      x=a(l,l)
      p=1.0
      q=1.0
      r=1.0
      call qrstep (a,v,p,q,r,l,l+2,n,na,nv)
      it=0
50    it=it+1
      if (it.le.30) go to 60
      fail=.true.
      return
60    continue
      p=a(l,l)-x
      q=a(l+1,l)
      r=0.0
      call qrstep (a,v,p,q,r,l,l+2,n,na,nv)
      if (abs(a(l+2,l+1)).gt.eps*(abs(a(l+1,l+1))+abs(a(l+2,l+2))))
     +    go to 50
      a(l+2,l+1)=0.0
      return
70    continue
      m=l+2
      if (b2.eq.2) m=m+1
      x=a(l+1,l+1)
      y=a(l,l)
      w=a(l+1,l)*a(l,l+1)
      p=1.0
      q=1.0
      r=1.0
      call qrstep (a,v,p,q,r,l,m,n,na,nv)
      it=0
80    it=it+1
      if (it.le.30) go to 90
      fail=.true.
      return
90    continue
      z=a(l,l)
      r=x-z
      s=y-z
      p=(r*s-w)/a(l+1,l)+a(l,l+1)
      q=a(l+1,l+1)-z-r-s
      r=a(l+2,l+1)
      s=abs(p)+abs(q)+abs(r)
      p=p/s
      q=q/s
      r=r/s
      call qrstep (a,v,p,q,r,l,m,n,na,nv)
      if (abs(a(m-1,m-2)).gt.eps*(abs(a(m-1,m-1))+abs(a(m-2,m-2))))
     +    go to 80
      a(m-1,m-2)=0.0
      return
c
c     last line of exchng
c
      end
      subroutine qrstep (a,v,p,q,r,nl,nu,n,na,nv)
c
c     based on stewart's hqr3 software.
c
c     *****parameters:
      integer n,na,nl,nu,nv
      real a(na,n),p,q,r,v(nv,n)
c
c     *****local variables:
      logical last
      integer i,j,k,nl2,nl3,num1
      real s,x,y,z
c
c     *****functions:
      integer min0
      real abs,sqrt
c
c     *****subroutines called:
c     none
c
c     ------------------------------------------------------------------
c
c     *****purpose:
c     this subroutine performs one implicit qr step on the upper
c     hessenberg matrix a. the shift is determined by the numbers p,q,
c     and r, and the step is applied to rows and columns nl throuth nu.
c     the transformations are accumulated in the array v.
c
c     *****parameter description:
c     on input:
c        na,nv            row dimensions of the arrays containing a
c                         and v, respectively, as declared in the
c                         calling program dimension statement;
c
c        a                n x n upper hessenberg matrix on which the qr
c                         step is to be performed;
c
c        p,q,r            parameters which determine the shift;
c
c        nl               the lower limit of the step;
c
c        nu               the upper limit of the step;
c
c        n                order of the matrix a.
c
c     on output:
c
c        v                n x n real scratch array containing the
c                         accumulated transformations.
c     ------------------------------------------------------------------
c
      nl2=nl+2
      do 10 i=nl2,nu
         a(i,i-2)=0.0
10    continue
      if (nl2.eq.nu) go to 30
      nl3=nl+3
      do 20 i=nl3,nu
         a(i,i-3)=0.0
20    continue
30    continue
      num1=nu-1
      do 130 k=nl,num1
         last=k.eq.num1
         if (k.eq.nl) go to 40
         p=a(k,k-1)
         q=a(k+1,k-1)
         r=0.0
         if (.not.last) r=a(k+2,k-1)
         x=abs(p)+abs(q)+abs(r)
         if (x.eq.0.0) go to 130
         p=p/x
         q=q/x
         r=r/x
40       continue
         s=sqrt(p**2+q**2+r**2)
         if (p.lt.0.0) s=-s
         if (k.eq.nl) go to 50
         a(k,k-1)=-s*x
         go to 60
50       continue
         if (nl.ne.1) a(k,k-1)=-a(k,k-1)
60       continue
         p=p+s
         x=p/s
         y=q/s
         z=r/s
         q=q/p
         r=r/p
         do 80 j=k,n
            p=a(k,j)+q*a(k+1,j)
            if (last) go to 70
            p=p+r*a(k+2,j)
            a(k+2,j)=a(k+2,j)-p*z
70          continue
            a(k+1,j)=a(k+1,j)-p*y
            a(k,j)=a(k,j)-p*x
80       continue
         j=min0(k+3,nu)
         do 100 i=1,j
            p=x*a(i,k)+y*a(i,k+1)
            if (last) go to 90
            p=p+z*a(i,k+2)
            a(i,k+2)=a(i,k+2)-p*r
90          continue
            a(i,k+1)=a(i,k+1)-p*q
            a(i,k)=a(i,k)-p
100      continue
         do 120 i=1,n
            p=x*v(i,k)+y*v(i,k+1)
            if (last) go to 110
            p=p+z*v(i,k+2)
            v(i,k+2)=v(i,k+2)-p*r
110         continue
            v(i,k+1)=v(i,k+1)-p*q
            v(i,k)=v(i,k)-p
120      continue
130   continue
      return
c
c     last line of qrstep
c
      end
      subroutine split (a,v,n,l,e1,e2,na,nv)
c
c     based on stewart's hqr3 software.
c
c     *****parameters:
      integer l,n,na,nv
      real a(na,n),v(nv,n),e1,e2
c
c     *****local variables:
      integer i,j,l1
      real p,q,r,t,u,w,x,y,z
c
c     *****functions:
      real abs,sqrt
c
c     ****subroutines called:
c     none
c
c     ------------------------------------------------------------------
c
c     *****purpose:
c     given the upper-hessenberg matrix a with a 2 x 2 block starting at
c     a(l,l), this program determines if the corresponding eigenvalues
c     are real or complex. if they are real, a rotation is determined
c     that reduces the block to upper-triangular form with the
c     eigenvalue of largest absolute value appearing first. the
c     rotation is accumulated in the array v.
c
c     *****parameter description:
c     on input:
c        na,nv            row dimensions of the arrays containing
c                         a and v, respectively, as declared in the
c                         calling program dimension statement;
c
c        a                the upper hessenberg matrix whose 2 x 2 block
c                         is to be split;
c
c        n                order of the matrix a;
c
c        l                position of the 2 x 2 block.
c
c     on output:
c
c        v                an n x n array containing the accumulated
c                         splitting transformation;
c
c        e1,e2            real scalars. if the eigenvalues are complex,
c                         e1 and e2 contain their common real part and
c                         positive imaginary part (respectively).
c                         if the eigenvalues are real, e1 contains the
c                         one largest in absolute value and e2 contains
c                         the other one.
c
c     ------------------------------------------------------------------
c
      x=a(l+1,l+1)
      y=a(l,l)
      w=a(l,l+1)*a(l+1,l)
      p=(y-x)/2.0
      q=p**2+w
      if (q.ge.0.0) go to 10
      e1=p+x
      e2=sqrt(-q)
      return
10    continue
      z=sqrt(q)
      if (p.lt.0.0) go to 20
      z=p+z
      go to 30
20    continue
      z=p-z
30    continue
      if (z.eq.0.0) go to 40
      r=-w/z
      go to 50
40    continue
      r=0.0
50    continue
      if (abs(x+z).ge.abs(x+r)) z=r
      y=y-x-z
      x=-z
      t=a(l,l+1)
      u=a(l+1,l)
      if (abs(y)+abs(u).le.abs(t)+abs(x)) go to 60
      q=u
      p=y
      go to 70
60    continue
      q=x
      p=t
70    continue
      r=sqrt(p**2+q**2)
      if (r.gt.0.0) go to 80
      e1=a(l,l)
      e2=a(l+1,l+1)
      a(l+1,l)=0.0
      return
80    continue
      p=p/r
      q=q/r
      do 90 j=l,n
         z=a(l,j)
         a(l,j)=p*z+q*a(l+1,j)
         a(l+1,j)=p*a(l+1,j)-q*z
90    continue
      l1=l+1
      do 100 i=1,l1
         z=a(i,l)
         a(i,l)=p*z+q*a(i,l+1)
         a(i,l+1)=p*a(i,l+1)-q*z
100   continue
      do 110 i=1,n
         z=v(i,l)
         v(i,l)=p*z+q*v(i,l+1)
         v(i,l+1)=p*v(i,l+1)-q*z
110   continue
      a(l+1,l)=0.0
      e1=a(l,l)
      e2=a(l+1,l+1)
      return
c
c     last line of split
c
      end
c
      subroutine ewset (n, itol, rtol, atol, ycur, ewt)
c
c author: Alan Hindmarsh, LLL
c
c-----------------------------------------------------------------------
c this subroutine sets the error weight vector ewt according to
c     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
c with the subscript on rtol and/or atol possibly replaced by 1 above,
c depending on the value of itol.
c-----------------------------------------------------------------------
      integer n, itol
      integer i
      real rtol, atol, ycur, ewt
      dimension rtol(*), atol(*), ycur(n), ewt(n)
c
      go to (10, 20, 30, 40), itol
 10   continue
      do 15 i = 1,n
 15     ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
      return
 20   continue
      do 25 i = 1,n
 25     ewt(i) = rtol(1)*abs(ycur(i)) + atol(i)
      return
 30   continue
      do 35 i = 1,n
 35     ewt(i) = rtol(i)*abs(ycur(i)) + atol(1)
      return
 40   continue
      do 45 i = 1,n
 45     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i)
      return
c----------------------- end of subroutine ewset -----------------------
      end
      real function vnorm (n, v, w)
c 
c revision date: July 10 1990                                            
c author: Luca Dieci, School of Mathematics, GaTech
c
c history: this is a rewriting of function vnorm of lsode
c
c-----------------------------------------------------------------------
c this function routine computes the weighted max-norm
c of the vector of length n contained in the array v, with weights
c contained in the array w of length n..
c   vnorm = amax1( v(i)*w(i), i=1,n )
c-----------------------------------------------------------------------
      integer n, i
      real v, w, sum
c *** functions
      real amax1, abs
      dimension v(n), w(n)
      sum = 0.0e0
      do 10 i = 1,n
 10     sum = amax1(sum, abs(v(i)*w(i)))
      vnorm = sum
      return
c----------------------- end of function vnorm -------------------------
      end
      subroutine srcom (rsav, isav, job)
c
c revision: July 10 1990
c modifications' author: Luca Dieci, School of Mathematics, GaTech
c
c history: identical to the routine srcom of lsode by A.Hindmarsh,
c          but we added an extra block to save-restore
c
c-----------------------------------------------------------------------
c this routine saves or restores (depending on job) the contents of
c the common blocks ls0001, eh0001 and intdre which are used internally
c by one or more dresol solvers.
c
c rsav = real array of length 222 or more.
c isav = integer array of length 66 or more.
c job  = flag indicating to save or restore the common blocks..
c        job  = 1 if common is to be saved (written to rsav/isav)
c        job  = 2 if common is to be restored (read from rsav/isav)
c        a call with job = 2 presumes a prior call with job = 1.
c-----------------------------------------------------------------------
      integer isav, job
      integer ieh, ils, iint
      integer i, lenils, lenrls, lenrnt, lenint
      real rsav, rls, rint
      dimension rsav(*), isav(*)
      common /ls0001/ rls(218), ils(39)
      common /eh0001/ ieh(2)
      common /intdre/ rint(4), iint(25)
      data lenrls/218/, lenils/39/, lenrnt/4/, lenint/25/
c
      if (job .eq. 2) go to 100
c
      do 10 i = 1,lenrls
 10     rsav(i) = rls(i)
      do 15 i = 1,lenrnt
 15     rsav(lenrls+i) = rint(i)
      do 20 i = 1,lenils
 20     isav(i) = ils(i)
      isav(lenils+1) = ieh(1)
      isav(lenils+2) = ieh(2)
      do 25 i = 1,lenint
 25      isav(lenils+2+i) = iint(i)
      return
c
 100  continue
      do 110 i = 1,lenrls
 110     rls(i) = rsav(i)
      do 115 i = 1,lenrnt
 115     rint(i) = rsav(lenrls+i)
      do 120 i = 1,lenils
 120     ils(i) = isav(i)
      ieh(1) = isav(lenils+1)
      ieh(2) = isav(lenils+2)
      do 125 i = 1,lenint
 125     iint(i) = isav(lenils+2+i)
      return
c----------------------- end of subroutine srcom -----------------------
      end
      real function r1mach (idum)
      integer idum
c 
c author: Alan Hindmarsh, LLL
c
c-----------------------------------------------------------------------
c this routine computes the unit roundoff of the machine.
c this is defined as the smallest positive machine number
c u such that  1.0 + u .ne. 1.0
c-----------------------------------------------------------------------
      real u, comp
      u = 1.0e0
 10   u = u*0.5e0
      comp = 1.0e0 + u
      if (comp .ne. 1.0e0) go to 10
      r1mach = u*2.0e0
      return
c----------------------- end of function r1mach ------------------------
      end
      subroutine xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
      integer msg, nmes, nerr, level, ni, i1, i2, nr,
     1   i, lun, lunit, mesflg, ncpw, nch, nwds
      real r1, r2
      dimension msg(nmes)
c author: Alan Hindmarsh, LLL
c-----------------------------------------------------------------------
c subroutines xerrwv, xsetf, and xsetun, as given here, constitute
c a simplified version of the slatec error handling package.
c written by a. c. hindmarsh at llnl.  version of march 30, 1987.
c ----------------------------------------------------------------------
c
c all arguments are input arguments.
c
c msg    = the message (hollerith literal or integer array).
c nmes   = the length of msg (number of characters).
c nerr   = the error number (not used).
c level  = the error level..
c          0 or 1 means recoverable (control returns to caller).
c          2 means fatal (run is aborted--see note below).
c ni     = number of integers (0, 1, or 2) to be printed with message.
c i1,i2  = integers to be printed, depending on ni.
c nr     = number of reals (0, 1, or 2) to be printed with message.
c r1,r2  = reals to be printed, depending on nr.
c
c note..  this routine is machine-dependent and specialized for use
c in limited context, in the following ways..
c 1. the number of hollerith characters stored per word, denoted
c    by ncpw below, is a data-loaded constant.
c 2. the value of nmes is assumed to be at most 60.
c    (multi-line messages are generated by repeated calls.)
c 3. if level = 2, control passes to the statement   stop
c    to abort the run.  this statement may be machine-dependent.
c 4. r1 and r2 are assumed to be in single precision and are printed
c    in e21.13 format.
c-----------------------------------------------------------------------
c the following are instructions for installing this routine
c in different machine environments.
c
c to change the default output unit, change the data statement in the 
c block data atlan2 below
c
c for a different number of characters per word, change the
c data statement setting ncpw below, and format 10.  alternatives for
c various computers are shown in comment cards.
c
c for a different run-abort command, change the statement following
c statement 100 at the end.
c-----------------------------------------------------------------------
      common /eh0001/ mesflg, lunit
c-----------------------------------------------------------------------
c the following data-loaded value of ncpw is valid for the cdc-6600
c and cdc-7600 computers.
c     data ncpw/10/
c the following is valid for the cray-1 and on the cyber 855 and 990.
c     data ncpw/8/
c the following is valid for the burroughs 6700 and 7800 computers.
c     data ncpw/6/
c the following is valid for the pdp-10 computer.
c     data ncpw/5/
c the following is valid for the vax computer with 4 bytes per integer,
c and for the ibm-360, ibm-370, ibm-303x, and ibm-43xx computers.
c and for the sun workstations sun3, sun4 and sparc
      data ncpw/4/
c the following is valid for the pdp-11, or vax with 2-byte integers.
c     data ncpw/2/
c-----------------------------------------------------------------------
c
      if (mesflg .eq. 0) go to 100
c get logical unit number. ---------------------------------------------
      lun = lunit
c get number of words in message. --------------------------------------
      nch = min0(nmes,60)
      nwds = nch/ncpw
      if (nch .ne. nwds*ncpw) nwds = nwds + 1
c write the message. ---------------------------------------------------
      write (lun, 10) (msg(i),i=1,nwds)
c-----------------------------------------------------------------------
c the following format statement is to have the form
c 10  format(1x,mmann)
c where nn = ncpw and mm is the smallest integer .ge. 60/ncpw.
c the following is valid for ncpw = 10.
c 10  format(1x,6a10)
c the following is valid for ncpw = 8.
c 10  format(1x,8a8)
c the following is valid for ncpw = 6.
c 10  format(1x,10a6)
c the following is valid for ncpw = 5.
c 10  format(1x,12a5)
c the following is valid for ncpw = 4.
  10  format(1x,15a4)
c the following is valid for ncpw = 2.
c 10  format(1x,30a2)
c-----------------------------------------------------------------------
      if (ni .eq. 1) write (lun, 20) i1
 20   format(6x,23hin above message,  i1 =,i10)
      if (ni .eq. 2) write (lun, 30) i1,i2
 30   format(6x,23hin above message,  i1 =,i10,3x,4hi2 =,i10)
      if (nr .eq. 1) write (lun, 40) r1
 40   format(6x,23hin above message,  r1 =,e21.13)
      if (nr .eq. 2) write (lun, 50) r1,r2
 50   format(6x,15hin above,  r1 =,e21.13,3x,4hr2 =,e21.13)
c abort the run if level = 2. ------------------------------------------
 100  if (level .ne. 2) return
      stop
c----------------------- end of subroutine xerrwv ----------------------
      end
      block data atlan2
c
c revision date: July 10 1990
c author: Luca Dieci, School of Mathematics, GaTech
c
c this block-data subprograms is needed to load the variables msflag and
c lunit into the common block eh0001.  
c
c    the variables are as follows..
c       mesflg = print control flag..
c                1 means print all messages (the default).
c                0 means no printing.
c       lunit  = logical unit number for messages.
c                the default is 6 (machine-dependent).
      common /eh0001/ mesflg, lunit
      data mesflg/1/, lunit/6/
c----------------------- end of bloch data atlan2 ----------------------
      end
      subroutine xsetf (mflag)
c 
c author: Alan Hindmarsh, LLL
c
c this routine resets the print control flag mflag.
c
      integer mflag, mesflg, lunit
      common /eh0001/ mesflg, lunit
c
      if (mflag .eq. 0 .or. mflag .eq. 1) mesflg = mflag
      return
c----------------------- end of subroutine xsetf -----------------------
      end
      subroutine xsetun (lun)
c
c author: Alan Hindmarsh, LLL
c
c this routine resets the logical unit number for messages.
c
      integer lun, mesflg, lunit
      common /eh0001/ mesflg, lunit
c
      if (lun .gt. 0) lunit = lun
      return
c----------------------- end of subroutine xsetun ----------------------
      end
