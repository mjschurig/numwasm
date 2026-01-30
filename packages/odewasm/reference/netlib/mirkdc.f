      subroutine mirkdc(method,tol,neqns,Nsub,MxNsub,
     +            mesh_input,mesh,sol_input,Y,ldY,output_control,info,
     +            iwork,work,init_de,init_Y,fsub,gsub,dfsub,dgsub)
c
c MIRKDC_4, Dec. 1999.
c
c Authors: Wayne Enright (enright@cs.toronto.edu), 
c          Paul Muir (muir@stmarys.ca).
c Additional programming assistance: Mark Adams, Nicole DeYoung.
c Thanks to Beta-Testers: Ian Gladwell, Larry Shampine, Richard Pancer.
c
c Note: This code makes use of the Level 1 BLAS routines. These must
c be linked with the MIRKDC package in order for it to work.
c
c***declaration of constants***
      integer Mxs
      parameter(Mxs=10)
c
c     `Mxs' Maximum number of stages of the Runge-Kutta method.
c
c***declaration of parameters***
c   imports:
      integer               method
      double precision      tol 
      integer               neqns, Nsub, MxNsub
      integer               mesh_input
      double precision      mesh(0:MxNsub)
      integer               sol_input
      double precision      Y(ldY*(Nsub+1))
      integer               ldY, output_control
c
c     `method' Defines the MIRK method to be used. Possible values 
c              for `method' and the corresponding MIRK schemes are:
c                    method  | MIRK formula
c                    ---------------------------------------------------------
c                    Externally documented methods.
c                     2      | MIRK221 scheme (trapezoidal rule, second order)
c                     4      | MIRK343 scheme (Lobatto scheme, fourth order)
c                     6      | MIRK563 scheme (6th order)
c                    Internally available methods.
c                    121     | MIRK121 scheme (midpoint rule, second order)
c                    221     | MIRK221 scheme (trapezoidal rule, second order)
c                    232     | MIRK232 scheme (one-sided, abscissa 0, 2/3)
c                    2321    | MIRK232 scheme (one-sided, abscissa 1, 1/3)
c                    343     | MIRK343 scheme (Lobatto scheme, fourth order)
c                    453     | MIRK453 scheme (one-sided, 5th order)
c                    563     | MIRK563 scheme (6th order)
c
c     `tol' Tolerance applied to an estimate of the defect of the approximate
c           solution; the defect is the amount by which the continuous 
c           approximate solution, u(t), fails to satisfy the ODE system. 
c           We require |def(t)|/(|f(t,u(t))|+1) to be less than `tol', where 
c           y'(t)=f(t,y(t)) is the ODE and def(t)= u'(t)-f(t,u(t)). 
c           The same tolerance is applied to all defect components.
c
c     `neqns' Number of differential equations. Also number of boundary
c                                                                  conditions.
c     `Nsub' On input, number of subintervals into which the problem
c            interval is initially partitioned. 
c     `MxNsub' Maximum number of subintervals.
c     `mesh_input' a) A value of 0 indicates that a uniform initial mesh 
c                     should be set up by `mirkdc'.
c                  b) A value of 1 indicates that an initial mesh has 
c                     already been provided by the user in the array `mesh'.
c                     This option should be chosen if the user has specific
c                     knowledge of the solution behavior.
c                  c) A value of 2 indicates that the user is employing  
c                     continuation to solve the problem and that the mesh 
c                     from a previous run will be stored in the array `mesh'.
c     `mesh' On input, the set of points that partition the problem
c            interval; initially empty if `mesh_input' = 0.
c     `sol_input' a) Normal use of `mirkdc' requires that `sol_input'=1,
c                    indicating that an initial solution estimate is to be 
c                    provided by the user through the `init_Y' routine. It 
c                    is recommended that `Y' not be set to a constant.
c                 b) The option `sol_input' = 2 should only be chosen when 
c                    using continuation to solve a problem, where `Y' takes 
c                    its value from the solution of a previous run. A value 
c                    of 2 indicates that an initial solution estimate at each
c                    mesh point is available in the array `Y'.
c     `Y' On input, if `sol_input'=2, this is the initial discrete 
c         approximation to the solution, evaluated at each meshpoint. 
c         Externally, in the interfaces with the user main program and
c         subroutines, `Y' is a two-dimensional array; the i-th column of `Y' 
c         stores the solution approximation for the i-th mesh point. Each 
c         solution approximation is of length `neqns'; the resulting 
c         dimensions of `Y' are `neqns x (Nsub+1)'. Inside MIRKDC it is 
c         convenient for `Y' to be stored as a column vector; in this 
c         one-dimensional version the solutions at each meshpoint come one 
c         after the other rather than side by side as they are in the 
c         two-dimensional array. 
c     `ldY' The leading dimension of Y. (`ldY' >= `neqns').
c     `output_control' Controls the amount of information output. Provides 
c                      profiling of code execution for activities such as the
c                      Newton iteration, mesh selection, or defect estimation.
c                             0 - None. 1 - Intermediate. 2 - Full.
c   exports:
c     double precision      mesh(0:Nsub)
c     double precision      Y(ldY*(Nsub+1))
      integer               info
c
c     `mesh' On output, the set of points which define the final mesh.
c     `Y' On output, if `info'=0, the converged discrete solution, evaluated
c         on the final mesh. 
c     `info' Communication flag: 
c                0 - implies a successful termination - solution obtained such 
c                    that the defect is within user specified tolerance.
c               -1 - implies an unsuccessful termination - the required size for
c                    the new mesh is be too large.
c            Internal flag:
c              The Newton iteration can fail for one of several reasons, each
c              signaled by a non-zero `info' value. Even if the Newton 
c              iteration converges, the resultant solution may be 
c              unsuitable because the defect is too large. During the 
c              execution of `mirkdc', `info' may take on the values given 
c              below. None of these conditions lead directly to an 
c              unsuccessful termination. Rather, the code will react by 
c              attempting to double the size of the mesh and try again. 
c                1 - too many iterations were taken without convergence being 
c                    obtained.
c                2 - a singular coefficient matrix was encountered during the 
c                    attempted solution of the Newton system.
c                3 - (a) it was impossible to obtain a suitable damping factor 
c                    for the Newton update (indicative of an "effectively 
c                    singular" Jacobian) or (b) evaluation of natural criterion
c                    function overflowed, indicative of a divergent iteration.
c                4 - the defect was too large and the solution from the Newton
c                    system could not be trusted.
c   work space:
      integer               iwork(neqns*(MxNsub+1))
      double precision      work(2*MxNsub*Mxs*neqns + 
     +                            5*neqns*MxNsub + 2*neqns**2*MxNsub +
     +                            6*neqns + 4*neqns**2 + 2*neqns**2*Mxs)      
c
c     `iwork' Integer work array provided to the `newiter' routine.
c     `work' Double precision work array used within `mirkdc' and provided to
c            `NewIter', `defect_estimate', `mesh_selector', and `interp_eval'.
c            The largest amount of workspace used is during the call to 
c            `Newiter'. The space required then is at most  
c                    2 * MxNsub*Mxs*neqns + 5*neqns*MxNsub + 
c                       2*neqns^2*MxNsub + 6*neqns + 4*neqns^2 + 2*neqns^2*Mxs.
c
c   user-supplied subroutines:
      external    init_de,init_Y,fsub,gsub,dfsub,dgsub
c
c     `init_de' User-supplied subroutine to initialize the problem dependent
c               variables, `leftbc', `a', and `b'.
c
c                  call init_de(leftbc, a, b)
c                    exports:
c                      integer leftbc
c                      double precision a,b 
c                      `leftbc' Number of boundary conditions applied at the
c                               left endpoint of the problem interval.
c                      `a' Left endpoint of the problem interval.
c                      `b' Right endpoint of the problem interval.
c
c     `init_Y' User-supplied subroutine to initialize the discrete solution
c              approximation, `Y'. (Needed if `sol_input'=1).
c
c                  call init_Y(Nsub, neqns, mesh, Y, ldY)
c                    imports:
c                      integer Nsub, neqns
c                      double precision mesh(0:Nsub) 
c                      integer ldY
c                      `Nsub' Number of subintervals in the initial mesh.
c                      `neqns' Number of differential equations in the
c                              first order system.
c                      `mesh' Points making up the partition of the problem
c                             interval.
c                      `ldY' Leading dimension of Y. (`ldY' >= `neqns').
c                    exports:
c                      double precision Y(ldY,Nsub+1)
c                      `Y' Initial approximation to the discrete
c                          solution at the points in the initial mesh.
c                          The i-th column of `Y' is the solution 
c                          approximation vector at mesh point, `mesh(i-1)'.
c               
c     `fsub' User-supplied subroutine which defines f(t,y) for the first order
c            system of differential equations, y' = f(t,y). 
c
c                  call fsub(neqns, t, y, f)
c                    imports:
c                      integer neqns
c                      double precision t, y(neqns) 
c                      `neqns' Number of differential equations in the
c                              first order system.
c                      `t' A point in the problem interval. 
c                      `y' Current solution approximation at t.
c                    exports:
c                      double precision f(neqns)
c                      `f' Value of f(t,y).
c
c     `gsub' User-supplied subroutine which defines the boundary condition
c            equations.
c
c                  call gsub(neqns, ya, yb, bc)
c                    imports:
c                      integer neqns
c                      double precision ya(neqns), yb(neqns)
c                      `neqns' Number of differential equations in the
c                              first order system. Also the total number
c                              of boundary conditions.
c                      `ya' Current solution approximation at left endpoint.
c                      `yb' Current solution approximation at right endpoint.
c                    exports:
c                      double precision bc(neqns)
c                      `bc' Boundary condition equations evaluated at
c                           `ya' and `yb'. The first `leftbc' components
c                           of `bc' are `g_a(ya)'; the remaining
c                           `neqns-leftbc' components of `bc' are `g_b(yb)'.
c
c     `dfsub' User-supplied subroutine which defines the Jacobian, df/dy,
c             of the system of differential equations.  Since the array 
c             `JJ' has already been set to zero, it is only necessary to 
c             set the non-zero elements.
c
c                  call dfsub(neqns, t, y, JJ)
c                    imports:
c                      integer neqns
c                      double precision t, y(neqns)
c                      `neqns' Number of differential equations in the
c                              first order system.
c                      `t' A point in the problem interval. 
c                      `y' Current solution approximation at t.
c                    exports:
c                      double precision JJ(neqns,neqns)
c                      `JJ' Contains the Jacobian, df/dy.
c
c     `dgsub' User-supplied subroutine which defines the Jacobian of
c             the boundary conditions, that is d(bc)/dy. The array `dbc' 
c             has been set to zero prior to the call to `dgsub', so it
c             is only necessary for the user to set the non-zero elements.
c
c                  call dgsub(neqns, ya, yb, dbc)
c                    imports:
c                      integer neqns
c                      double precision ya(neqns), yb(neqns)
c                      `neqns' Number of differential equations in the
c                              first order system. Also the total number
c                              of boundary conditions.
c                      `ya' Current solution approximation at left endpoint.
c                      `yb' Current solution approximation at right endpoint.
c                    exports:
c                      double precision dbc(neqns,neqns)
c                      `dbc' Contains the Jacobian of the boundary conditions.
c
c***declaration of local variables***
      integer              leftbc
      double precision     a,b,h
      integer              maxiter
      double precision     newton_tol
      double precision     defect_norm    
      integer              Nsub_star
      integer              i,j
c
c     `leftbc' The number of boundary conditions at the left end of the 
c              problem interval.
c     `a' Left endpoint.
c     `b' Right endpoint.
c     `h' Initial subinterval size, when a uniform mesh is employed.
c     `maxiter' Maximum number of Newton iterations.
c     `newton_tol' Tolerance for each Newton iteration.
c     `defect_norm' Max norm of the defect estimate.
c     `Nsub_star' Number of subintervals of the new mesh.
c     `i,j' Loop indices.
c
c***declaration of index related info for the work array***
c
c   Arrays local to `mirkdc' to be represented within the `work' array:
c          `k_discrete', `k_interp', `defect', `mesh_new', `Y_new'. 
c   Sizes:
      integer s_k_discrete, s_k_interp, s_defect
      integer s_mesh_new, s_Y_new
c
c     The actual array declarations will be replaced by assignments giving
c     the amount of space need within `work'.
c     double precision     k_discrete(Mxs*neqns*Nsub) => 
c                                       `s_k_discrete' = `Mxs*neqns*Nsub'
c     double precision     k_interp(Mxs*neqns*Nsub) =>
c                                        `s_k_interp' = `Mxs*neqns*Nsub'
c     double precision     defect(Nsub*neqns) => `s_defect' = `Nsub*neqns'
c     double precision     mesh_new(0:MxNsub) => `s_mesh_new' = `MxNsub+1'
c     double precision     Y_new((Nsub_star+1)*neqns) => 
c                                             `s_Y_new' = `(Nsub_star+1)*neqns'
c
c     `k_discrete' Stage information associated with the discrete 
c                  Runge-Kutta method: the i-th set of `s*neqns' locations
c                  contain the `s' vectors of length `neqns' associated
c                  with the Runge-Kutta method on the i-th subinterval.
c                  (`Mxs' is the maximum number of stages; `s' is the
c                  actual number of stages.)
c     `k_interp' Stage information associated with the continuous 
c                Runge-Kutta method; format is the same as `k_discrete'.
c     `defect'   Estimate of the (scaled) defect 
c                                (u'(t) - f(t,u(t)))/(|f(t,u(t))|+1)
c                on each subinterval (where u(t) is the continuous 
c                approximate solution). The defect estimate is a vector with
c                `neqns' components and there is one such vector for each
c                subinterval.
c     `mesh_new' Points defining the new mesh.
c     `Y_new'    Discrete approximation to the solution on `mesh_new'
c
c   Locations within `work':
      integer i_k_discrete, i_k_interp
      integer i_defect, i_mesh_new, i_Y_new
c
c     `i_k_discrete' = 1; 
c     `i_k_interp' = `s_k_discrete' + 1
c     `i_defect' = `s_k_discrete' + `s_k_interp' + 1
c     `i_mesh_new' = `s_k_discrete' + `s_k_interp' + `s_defect' + 1
c     `i_Y_new' = `s_k_discrete'+`s_k_interp'+`s_defect'+`s_mesh_new'+1
c
c   Index related info for the representation within `work' of arrays local 
c   to `newiter'. See `newiter' for explanations of the use of the arrays.
c   Sizes:
      integer s_PHI, s_delta, s_top, s_bot, s_blocks
c
c     `s_PHI' = `neqns * (Nsub+1)'
c     `s_delta' = `neqns * (Nsub+1)'
c     `s_top' = `neqns**2'
c     `s_bot' = `neqns**2'
c     `s_blocks' = `2 * neqns**2 * Nsub'
c
c   Locations within `work' during call to `newiter'.
      integer i_PHI, i_delta, i_top, i_bot, i_blocks
c
c     `i_PHI' = `s_k_discrete + s_k_interp + 1'
c     `i_delta' = `s_k_discrete + s_k_interp + s_PHI + 1'
c     `i_top' = `s_k_discrete + s_k_interp + s_PHI + s_delta + 1'
c     `i_bot' = `s_k_discrete + s_k_interp + s_PHI +  s_delta + s_top + 1'
c     `i_blocks' = `s_k_discrete+s_k_interp+s_PHI+s_delta+s_top+s_bot+1'
c
c   Index related info for the work array for arrays local to `defect_estimate'.c   See `defect_estimate' for explanations of the use of these arrays.
c   Sizes:
c     (All these arrays are of length `neqns'.)
c
c   Locations within `work' during the call to `defect_estimate'.
      integer i_z, i_z_prime, i_def_1, i_def_2
      integer i_f_1, i_f_2, i_temp_1, i_temp_2
c
c     `i_z = i_shift + 1; i_z_prime = i_shift + neqns + 1'
c     `i_def_1 = i_shift + 2*neqns + 1; i_def_2 = i_shift + 3*neqns + 1'
c     `i_f_1 = i_shift + 4*neqns + 1; i_f_2 = i_shift + 5*neqns + 1'
c     `i_temp_1 = i_shift + 6*neqns + 1; i_temp_2 = i_shift + 7*neqns + 1'
c     where `i_shift = s_k_discrete + s_k_interp + s_defect'
c
c   Other work indexes.
      integer i_work, i_shift
c
c     `i_work' Index to free part of work array passed to various routines
c              for use as a work array within those routines.
c     `i_shift' Temporary variable used as a shift to simplify `work' index
c               calculations.
c
c***declaration of variables for the /IOCNTRL/ common block***
c   imports:
      integer    print_level_0, print_level_1, print_level_2
      parameter  (print_level_0 = 0, print_level_1 = 1,
     *                print_level_2 = 2)
      integer    profile         
      common /IOCNTRL/ profile
c
c     `print_level_0' No output.
c     `print_level_1' Intermediate output.
c     `print_level_2' Full output.
c     `profile' Controls output, to standard output, of profiling infor-
c               mation such as Newton iteration counts, mesh selection,
c               relative defect estimate. This variable is identical to
c               the `output_control' parameter.
c
c
c     The following two interfaces are included only to allow convenient loading
c     of the integer work array upon successful exit. This work array is needed
c     when the user calls the `sol_eval' routine, subsequent to the return from
c     `mirkdc'.
c
c***declaration of variables from /RK_s/ common block***
c   imports:
      integer s
      common /RK_s/ s
c     `s' Number of discrete stages.
c
c***declaration of variables from /interp_s_star/ common block***
c   imports:
      integer s_star
      common /interp_s_star/ s_star
c     `s_star' Total number of stages required to form the
c              interpolant on each subinterval. It includes all the
c              stages of the discrete formula plus the additional 
c              stages required for the interpolant.
c----------------------------------------------------------------------------
c     Called by : main
c     Calls to  : `rk_tableau',`interp_tableau',`init_de',`init_Y',`NewIter',
c                 `defect_estimate',`mesh_selector',`interp_eval',`half_mesh'
c----------------------------------------------------------------------------
c***initialization***
c
c     Transfer output_control to common block variable `profile'.
      profile = output_control
c
      if (profile .GT. print_level_0) then
c       Output required size of `work' and `iwork' based on given
c       values for `neqns' and `MxNsub'.
        write(6,*)
        write(6,*) 'The input value for the order of the Runge-Kutta'  
        write(6,80) 'method is ',method,'.'
        write(6,90) 'The input value for the number of ODEs is '
     +                                              ,neqns,'.'
        write(6,*) 'The input value for the maximum number of '
        write(6,100) 'subintervals is ', MxNsub,'.'
        write(6,*) 
        write(6,*) 'Based on these values the size of the double'
        write(6,110) 'precision work array should be at least',
     +   MxNsub*Mxs*neqns +  5*neqns*MxNsub + 2*neqns**2*MxNsub +
     +                    6*neqns + 4*neqns**2 + 2*neqns**2*Mxs,'.'
        write(6,*) 'and the size of the integer work array should'
        write(6,120) 'be at least', neqns*(MxNsub+1),'.'
        write(6,*)
  80    format(1X,a10,i2,a1)
  90    format(1x,a41,i3,a1)
 100    format(1x,a15,i5,a1)
 110    format(1x,a39,i7,a1)
 120    format(1x,a11,i6,a1) 
      endif
c
c     Set up formula dependent coefficients.
      call rk_tableau(method)
      call interp_tableau(method)
c
c     Initialize remaining ODE dependent variables.
      call init_de(leftbc,a,b)
c
c     Setup initial mesh.
c     If `mesh_input'=1 or 2 then there is already a mesh stored in `mesh'.
c     If `mesh_input'=0 then setup uniform initial mesh of `Nsub' subintervals. 
      if (mesh_input.eq.0) then
            mesh(0) = a
            mesh(Nsub) = b
            h = (b - a)/Nsub
            do 10 i = 1, Nsub-1
               mesh(i) = i*h + a
 10         continue
      endif
c
      if (profile .GT. print_level_0) then
          write(6,130) 'The size of the initial mesh is '
     +                 ,Nsub,' subintervals.' 
          write(6,*)
 130      format(1x,a31,i5,a14)
      end if
      if (profile .GT. print_level_1) then
          write(6,140)'The initial mesh of',Nsub,' subintervals:'
          write(6,150)(mesh(i),i=0,Nsub)
          write(6,*)
  140     format(1x,a19,i5,a14)
  150     format(7F10.6)
      endif
c
c     Initialize approximate solution. 
c     If `sol_input' = 2 then an approximate solution is already stored in `Y'. 
c     If `sol_input' = 1 then call `init_Y' to get an initial solution.
      if (sol_input.eq.1) then
          call init_Y(Nsub,neqns,mesh,Y,ldY)
      endif
c
c     Transfer the 2-D array, `Y', to its 1-D equivalent for use within MIRKDC.
      do 15 i = 1,Nsub
           do 20 j = 1,neqns
             Y(i*neqns + j) = Y(i*ldY + j)
 20        continue
 15   continue  
c 
c     Define the Newton tolerance and maximum iteration count
      newton_tol = tol/100.0d0
      maxiter = 40
c
c***main outer loop to control discrete problem sequence***
c   REPEAT-UNTIL LOOP  (Implemented here as a GOTO-IF pair)
 1                                                   continue
c
c     Compute array sizes and indexes to provide access to the work
c     array during the call to `newiter'. `s_defect' is also computed
c     here since its value is required in more than one section of the code.
      s_k_discrete = Mxs * neqns * Nsub
      s_k_interp = Mxs * neqns * Nsub
      i_k_discrete = 1
      i_k_interp = s_k_discrete + 1
      s_PHI = neqns * (Nsub+1)
      s_delta = neqns * (Nsub+1)
      s_top = neqns**2
      s_bot = neqns**2
      s_blocks = 2 * neqns**2 * Nsub
      s_defect = Nsub*neqns
      i_PHI = s_k_discrete + s_k_interp + 1
      i_delta = s_k_discrete + s_k_interp + s_PHI + 1
      i_top = s_k_discrete + s_k_interp + s_PHI + 
     +          s_delta + 1
      i_bot = s_k_discrete + s_k_interp + s_PHI + 
     +          s_delta + s_top + 1
      i_blocks = s_k_discrete + s_k_interp + s_PHI + 
     +             s_delta + s_top + s_bot + 1
      i_work = s_k_discrete + s_k_interp + s_PHI + 
     +           s_delta + s_top + s_bot + s_blocks + 1
c          
c     Call `NewIter' to compute a solution for the discrete problem 
c     based on the current mesh.
c
      if (profile .GT. print_level_0) then
          write(6,*) 'Begin the Newton iteration:' 
          write(6,*)
      end if
c
      call NewIter(neqns,leftbc,Nsub,mesh,Y,
     +       newton_tol,maxiter,info,work(i_PHI),work(i_delta),
     +       work(i_top),work(i_bot),work(i_blocks),iwork,
     +       work(i_k_discrete),work(i_work),fsub,gsub,dfsub,dgsub)
c
      if ((info .NE. 0) .AND. (profile .EQ. print_level_1)) then 
              write(*,*)'the Newton iteration has failed.'
      endif
      if ((info .NE. 0 ) .AND. (profile .GT. print_level_1)) then
          if (info .EQ. 1) then
                write(6,*)'the Newton iteration failed because the',
     +          ' maximum number of allowed iterations was exceeded.'
          endif
          if (info .EQ. 2) then
                write(6,*)'the Newton iteration failed because a ',
     +          'singular matrix was encountered in the computation.'
          endif
          if (info .EQ. 3) then
                write(6,*)'the Newton iteration failed because a ',
     +          'sufficiently large damping factor could not be found.' 
          endif
      endif    
c
c
c     If the Newton iteration has converged, compute the defect estimate.
c
      if (info .EQ. 0) then
c
c         Set indexes to work array.
          i_defect = s_k_discrete + s_k_interp + 1
          i_shift = s_k_discrete + s_k_interp + s_defect
          i_z = i_shift + 1
          i_z_prime = i_shift + neqns + 1
          i_def_1 = i_shift + 2*neqns + 1
          i_def_2 = i_shift + 3*neqns + 1
          i_f_1 = i_shift + 4*neqns + 1
          i_f_2 = i_shift + 5*neqns + 1
          i_temp_1 = i_shift + 6*neqns + 1
          i_temp_2 = i_shift + 7*neqns + 1
          i_work = i_shift + 8*neqns + 1
c
          call defect_estimate(neqns,Nsub,mesh,Y,
     +         work(i_defect), defect_norm,info, 
     +         work(i_z),work(i_z_prime),
     +         work(i_def_1),work(i_def_2), work(i_f_1),
     +         work(i_f_2),work(i_temp_1),work(i_temp_2),
     +         work(i_k_discrete),work(i_k_interp),
     +         work(i_work),fsub)
c
c         The `defect_estimate' routine returns the defect estimate
c         vector along with its max norm. It also checks the size of 
c         the defect norm to make sure it is less than a certain 
c         threshold value (see the `defect_estimate' routine for 
c         further details.) If the defect norm is too large the 
c         `defect_estimate' routine sets `info' to 4 to indicate that 
c         the defect is too large and that the solution should not be trusted.
c
          if (profile .GT. print_level_0) then
                write(6,*)'the Newton iteration was successful. ',
     +          'Build a continuous approximation to the solution and',
     +          ' sample the defect.'
                write(6,*)'The norm of defect is ',defect_norm
                if (info .EQ. 4) then
                  write(6,*)'Since the defect is greater than ',
     +            '10% the solution is not acceptable.'  
                endif
          endif           
c
      endif
c     End of `if-endif (info .EQ. 0)'
c
c 
c     If the Newton iteration failed or the defect was not acceptable, 
c     `info' will not be 0 and we will proceed to the `else' clause of the 
c     following `if-then-else' statement where mesh halving will be 
c     attempted. Otherwise, we will proceed to the `then' clause below 
c     and attempt mesh redistribution.
c
c     Set work array indexes for `mesh_selector' and `half mesh'.
      s_mesh_new = MxNsub + 1
      i_mesh_new = s_k_discrete+s_k_interp
     +                                          +s_defect+1
      i_work = s_k_discrete+s_k_interp+s_defect
     +                                          +s_mesh_new+1
      i_Y_new = s_k_discrete+s_k_interp
     +                           +s_defect+s_mesh_new+1
c
c
      if (info.EQ.0) then                                              
c
c         If we have not yet satisfied the tolerance, we must try again
c         on a new mesh, if possible.
c
          if (defect_norm .GT. tol) then
c
            if (profile .GT. print_level_0) then
                  write(6,*) 'User defined tolerance',tol,' has not '
                  write(6,*) 'been satisfied.'
                  write(6,*) 'Construct a new mesh which ',
     +            'equidistributes the defect.'
            endif
c
c           Attempt to select a new mesh by calling `mesh_selector'.
            call mesh_selector(neqns,Nsub,mesh,work(i_defect),
     +         tol,Nsub_star,work(i_mesh_new),
     +           info,MxNsub,work(i_work))
c
c           If we were not successful in obtaining a new mesh, `info'
c           will not equal 0 and we should terminate the computation. 
c           Otherwise we can compute a new solution estimate on the new
c           mesh, and then update `Nsub', `mesh' and `Y'.
c
            if (info .EQ. 0) then
              if (profile .GT. print_level_0) then
                  write(6,*)
                  write(6,*)
                  write(6,160)'New mesh will be of size ',
     +                         Nsub_star,' subintervals.'
 160              format(1x,a25,i5,a14)
              end if
              if (profile .GT. print_level_1) then
                   write(6,150)(work(i_mesh_new+i), i=0, Nsub_star)
                   write(6,*)
              end if
c
c             Use current computed solution to generate initial guess for next 
c             problem. Compute solution approximations at new mesh points.
c
c             Set work array indexes.
              s_Y_new = neqns*(Nsub_star+1)
              i_Y_new = s_k_discrete+s_k_interp
     +                              +s_defect+s_mesh_new+1
              i_work = s_k_discrete+s_k_interp
     +                        +s_defect+s_mesh_new+s_Y_new+1
c
              do 30 i = 0, Nsub_star
                  call interp_eval(neqns,Nsub,mesh,Y,
     +              work(i_mesh_new+i),
     +              work(i_Y_new+i*neqns), work(i_work), 
     +              work(i_k_discrete),work(i_k_interp))
 30           continue
c
c             Copy `mesh_new' to `mesh', `Y_new' to `Y'; update `Nsub'.
              call dcopy(Nsub_star+1,work(i_mesh_new),1,mesh,1)
              call dcopy((Nsub_star+1)*neqns,work(i_Y_new),1,Y,1)
              Nsub = Nsub_star
c
            endif
c           End of `if-endif (info .EQ. 0)'
c
          endif
c         End of `if-endif (defect_norm .GT. tol)'
c
          if ( (info. EQ. 0) .AND. (defect_norm .LT. tol) .AND.
     +         (profile .GT. print_level_0) ) then
                  write(6,*)'The user defined tolerance',tol,' has ', 
     +            'been satisfied.'
                  write(6,*)'Successful computation. '
                  write(6,*)
          endif
c
      else
c     Else clause associated with `if-then-else (info .EQ. 0)'
c
          if (profile .GT. print_level_0) then
            write(6,*)'Cannot obtain a solution for the current mesh.'
          endif
c
c         Simply half the current mesh, if possible, and try again.
          if (2*Nsub .GT. MxNsub) then
c              New mesh would be too large. Terminate the computation.
               if (profile .GT. print_level_0) then
                  write(6,*)'New mesh would be too large.'
               end if
               info = -1
          else
c           Else clause associated with if_then_else `2*Nsub .GT. MxNsub'
c
            call half_mesh(Nsub,mesh,work(i_mesh_new))
            Nsub_star = 2*Nsub
c
c           Copy `mesh_new' to `mesh'
            call dcopy(Nsub_star+1,work(i_mesh_new),1,mesh,1)
            Nsub = Nsub_star
c            
            if (profile .GT. print_level_0) then
                  write(6,*)'Double the mesh and try again.'
                  write(6,*)
                  write(6,*)
                  write(6,170)'New mesh will be of size', Nsub,'.'
  170             format(1x,a24,i5,a1)
            end if
            if (profile .GT. print_level_1) then
                   write(6,150)(mesh(i), i = 0, Nsub)
                   write(6,*)
            end if 
c
c           Initialize the approximate solution - store in `Y_new'.
            if (sol_input .EQ. 1) then
                  call init_Y(Nsub,neqns,mesh,Y,ldY)
c                 Transferring the 2-D array, Y, to its 1-D equivalent for
c                 use with in MIRKDC.
                  do 35 i = 1,Nsub
                    do 40 j = 1,neqns
                      Y(i*neqns + j) = Y(i*ldY + j)
  40                continue
  35              continue 
            endif
c           If sol_input = 2 then the initial guess on the initial mesh
c           was provided through `Y' and `mesh', respectively. Since this
c           initial guess is no longer usable, with the new mesh, we
c           continue by simply setting `Y_new' to zero.
            if (sol_input .EQ. 2) then
                  if (profile .GT. print_level_0) then
                    write(6,*)'Iteration failed and there is no ',
     +                    'Init_Y routine; Y will be set to zero.'
                    write(6,*)'It would be better to retry the ',
     +                 'continuation sequence with smaller steps'
                  endif
                  call dcopy((Nsub+1)*neqns,0.0d0,0,Y,1)
            endif
c  
c           Reset `info' to 0 and `defect_norm' > `tol' 
c           to force a restart. (See UNTIL conditions below.)
            info = 0            
            defect_norm = 2* tol
c
          endif
c         End of `if-then-else (2*Nsub .GT. MxNsub)'
c
      endif
c     End of `if-then-else (info .EQ. 0)', `Acceptable solution obtained'
c
c   UNTIL( (info .NE. 0) .OR. (defect_norm .LE. tol))
                                       if((info .EQ. 0) .AND.
     +                                 (defect_norm .GT. tol))goto1
c
c***conclusion***
c
c     Set values of integer work array to allow user convenient 
c     use of `sol_eval' routine, after a successful termination.
c     See the `sol_eval' routine for further details.
      if (info .EQ. 0) then
        iwork(1) = neqns
        iwork(2) = Nsub
        iwork(3) = s
        iwork(4) = s_star
        iwork(5) = method
c      
c       Set values of double precision work array to allow user 
c       convenient use of `sol_eval' routine after a successful termination.
        call dcopy(Nsub+1,mesh(0),1,work(2*Mxs*neqns*Nsub+1),1)
        call dcopy((Nsub+1)*neqns,Y,1,
     +                 work(2*Mxs*neqns*Nsub+(Nsub+1)+1),1)
c
c       Transfer the 1-D version of `Y' back to its 2-D equivalent.
        do 45 i = 1,(Nsub-1),-1
          do 50 j = 1,neqns
           Y(i*ldY + j) = Y(i*neqns + j) 
 50       continue
 45     continue
      else
        if (profile.GT.print_level_0) then
          write(6,*)
          write(6,*) 'The computation has failed.'
          write(6,*)
          write(6,*) '1) If the profiling information begins with and ',
     +    'simply proceeds through a series of failed Newton ',
     +    'iterations each followed by a doubling of the mesh, then',
     +    'there may be an error in the coding of the problem.',
     +    'It is recommended that the routines fsub, gsub, dfsub and',
     +    ' dgsub be checked.'
          write(6,*) 
          write(6,*) '2) If the profiling information shows a series',
     +    'of successful, converged Newton iterations with a halt',
     +    'simply because the value of MxNsub was exceeded, then it',
     +    'is recommended that the program be run again with a larger',
     +    ' value of MxNsub.'
          write(6,*)
          write(6,*) '3) More information may be obtained by choosing',
     +    'a larger value for "output_control".'
        endif 
      endif
c
      return
      end

       SUBROUTINE COLROW(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,X,IFLAG)
C
C***************************************************************
C
C  THIS PROGRAM SOLVES THE LINEAR SYSTEM  A*X = B  WHERE  A IS
C  AN ALMOST BLOCK DIAGONAL MATRIX OF THE FORM
C
C               TOPBLK
C               ARRAY(1)
C                     ARRAY(2)
C                          .
C                             .
C                                .
C                                   .
C                                    ARRAY(NBLOKS)
C                                           BOTBLK
C
C  WHERE
C           TOPBLK IS  NRWTOP  BY NOVRLP
C           ARRAY(K), K=1,NBLOKS, ARE NRWBLK BY NRWBLK+NOVRLP
C           BOTBLK IS NRWBOT BY NOVRLP,
C  AND
C           NOVRLP = NRWTOP + NRWBOT
C  WITH
C           NOVRLP.LE.NRWBLK .
C
C  THE LINEAR SYSTEM IS OF ORDER  N = NBLOKS*NRWBLK + NOVRLP.
C
C  THE METHOD IMPLEMENTED IS BASED ON GAUSS ELIMINATION WITH
C  ALTERNATE ROW AND COLUMN ELIMINATION WITH PARTIAL PIVOTING,
C  WHICH PRODUCES A STABLE DECOMPOSITION OF THE MATRIX  A
C  WITHOUT INTRODUCING FILL-IN.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  TO OBTAIN A SINGLE PRECISION VERSION OF THIS PACKAGE, REMOVE
C  ALL DOUBLE PRECISION STATEMENTS.  THERE IS ONE SUCH STATEMENT
C  IN C O L R O W, THREE IN C R D C M P, AND TWO IN C R S O L V.
C  IN ADDITION, REFERENCES TO BUILT-IN FUNCTIONS DABS AND DMAX1
C  MUST BE REPLACED BY ABS AND AMAX1, RESPECTIVELY.  DABS OCCURS
C  NINE TIMES, IN C R D C M P.  DMAX1 OCCURS FOUR TIMES, IN
C  C R D C M P.  FINALLY, ZERO IS INITIALISED TO 0.D0 IN A
C  DATA STATEMENT IN C R D C M P.  THIS MUST BE REPLACED BY:
C               DATA ZERO/0.0/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C                    X - DOUBLE PRECISION(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C                    X - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR (IF IFLAG = 0)
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  AUXILIARY PROGRAMS  *****
C
C       CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
C    *     BOTBLK,NRWBOT,PIVOT,IFLAG)
C            - DECOMPOSES THE MATRIX  A  USING MODIFIED
C              ALTERNATE ROW AND COLUMN ELIMINATON WITH
C              PARTIAL PIVOTING, AND IS USED FOR THIS
C              PURPOSE IN C O L R O W.
C              THE ARGUMENTS ARE AS IN C O L R O W.
C
C       CRSLVE(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
C    *     BOTBLK,NRWBOT,PIVOT,B,X)
C            - SOLVES THE SYSTEM A*X = B ONCE A IS DECOMPOSED.
C              THE ARGUMENTS ARE ALLAS IN C O L R O W.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       THE SUBROUTINE  C O L R O W  AUTOMATICALLY SOLVES THE
C  INPUT SYSTEM WHEN IFLAG=0.  C O L R O W  IS CALLED ONLY ONCE
C  FOR A GIVEN SYSTEM. THE SOLUTION FOR A SEQUENCE OF P RIGHT
C  HAND SIDES CAN BE OBTAINED BY ONE CALL TO  C O L R O W  AND
C  P-1 CALLS TO CRSLVE ONLY. SINCE THE ARRAYS TOPBLK,ARRAY,
C  BOTBLK AND PIVOT CONTAIN THE DECOMPOSITION OF THE GIVEN
C  COEFFICIENT MATRIX AND PIVOTING INFORMATION ON RETURN FROM
C  C O L R O W , THEY MUST NOT BE ALTERED BETWEEN SUCCESSIVE
C  CALLS TO CRSLVE WITH THE SAME LEFT HAND SIDES. FOR THE
C  SAME REASON, IF THE USER WISHES TO SAVE THE COEFFICIENT
C  MATRIX, THE ARRAYS TOPBLK,ARRAY,BOTBLK MUST BE COPIED
C  BEFORE A CALL TO  C O L R O W .
C
C*************************************************************************
C
C             *****  SAMPLE CALLING PROGRAM  *****
C
C      THE FOLLOWING PROGRAM WILL EXERCISE COLROW, IN THE
C      CASE WHEN THE COEFFICIENT MATRIX IS NON-SINGULAR.
C
C       DOUBLE PRECISION TOP,AR,BOT,B,X
C       DOUBLE PRECISION ERROR,ERR
C       DIMENSION TOP(2,4),AR(4,8,2),BOT(2,4),B(12),X(12)
C       INTEGER PIVOT(12)
C       DATA N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT/12,2,4,4,8,2,2/
C       DATA TOP(1,1),TOP(1,2),TOP(1,3),TOP(1,4),
C    *       TOP(2,1),TOP(2,2),TOP(2,3),TOP(2,4)/
C    *0.0 D0,-0.98D0,-0.79D0,-0.15D0,
C    *-1.00D0, 0.25D0,-0.87D0, 0.35D0/
C       DATA AR(1,1,1),AR(1,2,1),AR(1,3,1),AR(1,4,1),
C    *       AR(1,5,1),AR(1,6,1),AR(1,7,1),AR(1,8,1)/
C    *0.78D0, 0.31D0,-0.85D0, 0.89D0,-0.69D0,-0.98D0,-0.76D0,-0.82D0/
C       DATA AR(2,1,1),AR(2,2,1),AR(2,3,1),AR(2,4,1),
C    *       AR(2,5,1),AR(2,6,1),AR(2,7,1),AR(2,8,1)/
C    *0.12D0,-0.01D0, 0.75D0, 0.32D0,-1.00D0,-0.53D0,-0.83D0,-0.98D0/
C       DATA AR(3,1,1),AR(3,2,1),AR(3,3,1),AR(3,4,1),
C    *       AR(3,5,1),AR(3,6,1),AR(3,7,1),AR(3,8,1)/
C    *-0.58D0, 0.04D0, 0.87D0, 0.38D0,-1.00D0,-0.21D0,-0.93D0,-0.84D0/
C       DATA AR(4,1,1),AR(4,2,1),AR(4,3,1),AR(4,4,1),
C    *       AR(4,5,1),AR(4,6,1),AR(4,7,1),AR(4,8,1)/
C    *-0.21D0,-0.91D0,-0.09D0,-0.62D0,-1.99D0,-1.12D0,-1.21D0, 0.07D0/
C       DATA AR(1,1,2),AR(1,2,2),AR(1,3,2),AR(1,4,2),
C    *       AR(1,5,2),AR(1,6,2),AR(1,7,2),AR(1,8,2)/
C    *0.78D0,-0.93D0,-0.76D0, 0.48D0,-0.87D0,-0.14D0,-1.00D0,-0.59D0/
C       DATA AR(2,1,2),AR(2,2,2),AR(2,3,2),AR(2,4,2),
C    *       AR(2,5,2),AR(2,6,2),AR(2,7,2),AR(2,8,2)/
C    *-0.99D0, 0.21D0,-0.73D0,-0.48D0,-0.93D0,-0.91D0, 0.10D0,-0.89D0/
C       DATA AR(3,1,2),AR(3,2,2),AR(3,3,2),AR(3,4,2),
C    *       AR(3,5,2),AR(3,6,2),AR(3,7,2),AR(3,8,2)/
C    *-0.68D0,-0.09D0,-0.58D0,-0.21D0, 0.85D0,-0.39D0, 0.79D0,-0.71D0/
C       DATA AR(4,1,2),AR(4,2,2),AR(4,3,2),AR(4,4,2),
C    *       AR(4,5,2),AR(4,6,2),AR(4,7,2),AR(4,8,2)/
C    *0.39D0,-0.99D0,-0.12D0,-0.75D0,-0.68D0,-0.99D0, 0.50D0,-0.88D0/
C       DATA BOT(1,1),BOT(1,2),BOT(1,3),BOT(1,4),
C    *       BOT(2,1),BOT(2,2),BOT(2,3),BOT(2,4)/
C    *0.71D0,-0.64D0, 0.0 D0, 0.48D0,
C    *0.08D0,100.0D0,50.00D0,15.00D0/
C       DATA B(1),B(2),B(3),B(4),B(5),B(6),B(7),B(8),B(9),B(10),B(11),
C    *       B(12)/
C    *-1.92D0,-1.27D0,-2.12D0,-2.16D0,-2.27D0,-6.08D0,-3.03D0,
C    *-4.62D0,-1.02D0,-3.52D0,.55D0,165.08D0/
C
C*************************************************************************
C
C   THE INPUT MATRIX IS GIVEN BY:
C
C  0.0  -0.98 -0.79 -0.15
C -1.00  0.25 -0.87  0.35
C  0.78  0.31 -0.85  0.89 -0.69 -0.98 -0.76 -0.82
C  0.12 -0.01  0.75  0.32 -1.00 -0.53 -0.83 -0.98
C -0.58  0.04  0.87  0.38 -1.00 -0.21 -0.93 -0.84
C -0.21 -0.91 -0.09 -0.62 -1.99 -1.12 -1.21  0.07
C                          0.78 -0.93 -0.76  0.48 -0.87 -0.14 -1.00 -0.59
C                         -0.99  0.21 -0.73 -0.48 -0.93 -0.91  0.10 -0.89
C                         -0.68 -0.09 -0.58 -0.21  0.85 -0.39  0.79 -0.71
C                          0.39 -0.99 -0.12 -0.75 -0.68 -0.99  0.50 -0.88
C                                                  0.71 -0.64  0.0   0.48
C                                                  0.08 100.0 50.00 15.00
C
C       THE RIGHT HAND SIDE IS GIVEN BY:
C
C         B = (-1.92,-1.27,-2.12,-2.16,-2.27,-6.08,-3.03,-4.62,
C              -1.02,-3.52,0.55,165.08)
C
C       THE SOLUTION OF THIS SYSTEM IS GIVEN BY;
C
C          X = (1,1,1,1,1,1,1,1,1,1,1,1)
C
C*************************************************************************
C
C       CALL COLROW(N,TOP,NRWTOP,NOVRLP,AR,NRWBLK,NCLBLK,NBLOKS,
C    *              BOT,NRWBOT,PIVOT,B,X,IFLAG)
C       IF(IFLAG.NE.0)GO TO 1000
C       ERROR = 0.D0
C       DO 10 I=1,N
C          ERR = 1.D0 - X(I)
C          ERROR = DMAX1(ERROR,DABS(ERR))
C          WRITE(6,100)X(I),ERR
C  10   CONTINUE
C       WRITE(6,200)ERROR
C 200   FORMAT(12H MAX ERROR = ,D15.7)
C 100   FORMAT(1H ,F15.7,D15.7)
C       RETURN
C1000   CONTINUE
C       WRITE(6,300)IFLAG
C 300   FORMAT(9H IFLAG =  ,I3)
C       RETURN
C       END
C
C***************************************************************
C
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,B,X
        INTEGER PIVOT(1)
        DIMENSION TOPBLK(NRWTOP,1),ARRAY(NRWBLK,NCLBLK,1),
     *          BOTBLK(NRWBOT,1),B(1),X(1)
        CALL CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,IFLAG)
        IF(IFLAG.NE.0)RETURN
        CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,B,X)
        RETURN
        END
        SUBROUTINE CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,IFLAG)
C
C***************************************************************
C
C  C R D C M P DECOMPOSES THE ALMOST BLOCK DIAGONAL MATRIX A
C  USING MODIFIED ALTERNATE ROW AND COLUMN ELIMINATION WITH
C  PARTIAL PIVOTING.  THE MATRIX  A  IS STORED IN THE ARRAYS
C  TOPBLK, ARRAY, AND BOTBLK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A TO BE DECOMPOSED
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
C***************************************************************
C
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK
        DOUBLE PRECISION ROWMAX,ROWPIV,ROWMLT,COLMAX,COLPIV
        DOUBLE PRECISION SWAP,COLMLT,PIVMAX,ZERO,TEMPIV
        INTEGER PIVOT(1)
        DIMENSION TOPBLK(NRWTOP,1),ARRAY(NRWBLK,NCLBLK,1),
     *          BOTBLK(NRWBOT,1)
        DATA ZERO/0.0D0/
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        IFLAG = 0
        PIVMAX = ZERO
        NRWTP1 = NRWTOP+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
C
C***************************************************************
C
C          ****  CHECK VALIDITY OF THE INPUT PARAMETERS....
C
C               IF PARAMETERS ARE INVALID THEN TERMINATE AT 10;
C                                         ELSE CONTINUE AT 100.
C
C***************************************************************
C
        IF(N.NE.NBLOKS*NRWBLK+NOVRLP)GO TO 10
        IF(NOVRLP.NE.NRWTOP+NRWBOT)GO TO 10
        IF(NCLBLK.NE.NOVRLP+NRWBLK)GO TO 10
        IF(NOVRLP.GT.NRWBLK)GO TO 10
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        GO TO 100
10      CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IFLAG = 1
        RETURN
100     CONTINUE
C
C***************************************************************
C
C               ****  FIRST, IN TOPBLK....
C
C***************************************************************
C
C          ***  APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN
C                 PIVOTING ....
C
C***************************************************************
C
        DO 190 I = 1,NRWTOP
           IPLUS1 = I+1
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IPVT = I
           COLMAX = DABS(TOPBLK(I,I))
           DO 110 J = IPLUS1,NOVRLP
              TEMPIV = DABS(TOPBLK(I,J))
              IF(TEMPIV.LE.COLMAX)GO TO 110
                 IPVT = J
                 COLMAX = TEMPIV
110        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
           PIVMAX = DMAX1(COLMAX,PIVMAX)
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           PIVOT(I) = IPVT
           IF(IPVT.EQ.I)GO TO 140
              DO 120 L = I,NRWTOP
                 SWAP = TOPBLK(L,IPVT)
                 TOPBLK(L,IPVT) = TOPBLK(L,I)
                 TOPBLK(L,I) = SWAP
120           CONTINUE
              DO 130 L = 1,NRWBLK
                 SWAP = ARRAY(L,IPVT,1)
                 ARRAY(L,IPVT,1) = ARRAY(L,I,1)
                 ARRAY(L,I,1) = SWAP
130           CONTINUE
140        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           COLPIV = TOPBLK(I,I)
           DO 180 J = IPLUS1,NOVRLP
              COLMLT = TOPBLK(I,J)/COLPIV
              TOPBLK(I,J) = COLMLT
              IF(IPLUS1.GT.NRWTOP)GO TO 160
                 DO 150 L = IPLUS1,NRWTOP
                    TOPBLK(L,J) = TOPBLK(L,J)-COLMLT*TOPBLK(L,I)
150              CONTINUE
160           CONTINUE
              DO 170 L = 1,NRWBLK
                 ARRAY(L,J,1) = ARRAY(L,J,1)-COLMLT*ARRAY(L,I,1)
170           CONTINUE
180        CONTINUE
190     CONTINUE
C
C***************************************************************
C
C          ****  IN EACH BLOCK ARRAY(,,K)....
C
C***************************************************************
C
        INCR = 0
        DO 395 K = 1,NBLOKS
           KPLUS1 = K+1
C
C          *****************************************************
C
C          ***  FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH
C                       ROW PIVOTING....
C
C          *****************************************************
C
           DO 270 J = NRWTP1,NRWBLK
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = DABS(ARRAY(JMINN,J,K))
              LOOP = JMINN+1
              DO 210 I = LOOP,NRWBLK
                 TEMPIV = DABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.ROWMAX)GO TO 210
                 IPVT = I
                 ROWMAX = TEMPIV
210           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO  1000
              PIVMAX = DMAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 230
                 DO 220 L = J,NCLBLK
                    SWAP = ARRAY(IPVT,L,K)
                    ARRAY(IPVT,L,K) = ARRAY(JMINN,L,K)
                    ARRAY(JMINN,L,K) = SWAP
220              CONTINUE
230           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = ARRAY(JMINN,J,K)
              DO 240 I = LOOP,NRWBLK
                 ARRAY(I,J,K) = ARRAY(I,J,K)/ROWPIV
240           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 260 L = JPLUS1,NCLBLK
                 ROWMLT = ARRAY(JMINN,L,K)
                 DO 250 I = LOOP,NRWBLK
                    ARRAY(I,L,K) = ARRAY(I,L,K)
     *                                -ROWMLT*ARRAY(I,J,K)
250              CONTINUE
260           CONTINUE
270        CONTINUE
C
C          *****************************************************
C
C          ***  NOW APPLY NRWTOP COLUMN ELIMINATIONS WITH
C                      COLUMN PIVOTING....
C
C          *****************************************************
C
           DO 390 I = NRWEL1,NRWBLK
              IPLUSN = I+NRWTOP
              IPLUS1 = I+1
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = IPLUSN
              COLMAX = DABS(ARRAY(I,IPVT,K))
              LOOP = IPLUSN+1
              DO 310 J = LOOP,NCLBLK
                 TEMPIV = DABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.COLMAX)GO TO 310
                 IPVT = J
                 COLMAX = TEMPIV
310           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = DMAX1(COLMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRN = INCR+IPLUSN
              PIVOT(INCRN) = INCR+IPVT
              IRWBLK = IPLUSN-NRWBLK
              IF(IPVT.EQ.IPLUSN)GO TO 340
                 DO 315 L = I,NRWBLK
                    SWAP = ARRAY(L,IPVT,K)
                    ARRAY(L,IPVT,K) = ARRAY(L,IPLUSN,K)
                    ARRAY(L,IPLUSN,K) = SWAP
315              CONTINUE
                 IPVBLK = IPVT-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 330
                    DO 320 L = 1,NRWBLK
                       SWAP = ARRAY(L,IPVBLK,KPLUS1)
                       ARRAY(L,IPVBLK,KPLUS1)
     *                                 = ARRAY(L,IRWBLK,KPLUS1)
                       ARRAY(L,IRWBLK,KPLUS1) = SWAP
320                 CONTINUE
                    GO TO 340
330              CONTINUE
                 DO 335 L = 1,NRWBOT
                    SWAP = BOTBLK(L,IPVBLK)
                    BOTBLK(L,IPVBLK) = BOTBLK(L,IRWBLK)
                    BOTBLK(L,IRWBLK) = SWAP
335              CONTINUE
340           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              COLPIV = ARRAY(I,IPLUSN,K)
              DO 380 J = LOOP,NCLBLK
                 COLMLT = ARRAY(I,J,K)/COLPIV
                 ARRAY(I,J,K) = COLMLT
                 IF(I.EQ.NRWBLK)GO TO 350
                    DO 345 L = IPLUS1,NRWBLK
                       ARRAY(L,J,K) = ARRAY(L,J,K)
     *                                -COLMLT*ARRAY(L,IPLUSN,K)
345                 CONTINUE
350              CONTINUE
                 JRWBLK = J-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 370
                    DO 360 L = 1,NRWBLK
                       ARRAY(L,JRWBLK,KPLUS1) =
     *                                  ARRAY(L,JRWBLK,KPLUS1)
     *                         -COLMLT*ARRAY(L,IRWBLK,KPLUS1)
360                 CONTINUE
                    GO TO 380
370              CONTINUE
                 DO 375 L = 1,NRWBOT
                    BOTBLK(L,JRWBLK) = BOTBLK(L,JRWBLK)
     *                              -COLMLT*BOTBLK(L,IRWBLK)
375              CONTINUE
380           CONTINUE
390        CONTINUE
           INCR = INCR + NRWBLK
395     CONTINUE
C
C***************************************************************
C
C          ****  FINALLY, IN BOTBLK....
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          ***  APPLY NRWBOT ROW ELIMINATIONS WITH ROW
C                  PIVOTING....
C
C               IF BOT HAS JUST ONE ROW GO TO 500
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 500
           DO 470 J = NRWTP1,NVRLP0
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = DABS(BOTBLK(JMINN,J))
              LOOP = JMINN+1
              DO 410 I = LOOP,NRWBOT
                 TEMPIV = DABS(BOTBLK(I,J))
                 IF(TEMPIV.LE.ROWMAX) GO TO 410
                 IPVT = I
                 ROWMAX = TEMPIV
410           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = DMAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 430
                 DO 420 L = J,NOVRLP
                    SWAP = BOTBLK(IPVT,L)
                    BOTBLK(IPVT,L) = BOTBLK(JMINN,L)
                    BOTBLK(JMINN,L) = SWAP
420              CONTINUE
430           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = BOTBLK(JMINN,J)
              DO 440 I = LOOP,NRWBOT
                 BOTBLK(I,J) = BOTBLK(I,J)/ROWPIV
440           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 460 L = JPLUS1,NOVRLP
                 ROWMLT = BOTBLK(JMINN,L)
                 DO 450 I = LOOP,NRWBOT
                    BOTBLK(I,L) = BOTBLK(I,L)-ROWMLT*BOTBLK(I,J)
450              CONTINUE
460           CONTINUE
470        CONTINUE
500     CONTINUE
C
C***************************************************************
C
C          DONE PROVIDED THE LAST ELEMENT IS NOT ZERO
C
C***************************************************************
C
        IF(PIVMAX+DABS(BOTBLK(NRWBOT,NOVRLP)).NE.PIVMAX) RETURN
C
C***************************************************************
C
C       ****  MATRIX IS SINGULAR - SET IFLAG = - 1.
C                                  TERMINATE AT 1000.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
1000    CONTINUE
        IFLAG = -1
        RETURN
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine criterion(neqns, leftbc, Nsub, mesh, Y, top,blocks,
     +           bot, pivot, PHI, delta, g,k_discrete,work,fsub,gsub)
c
c***overview***
c
c        This routine takes as input the current iterate, `Y' and computes
c     a related function known as the natural criterion function, `g' = g(Y) = 
c     0.5*|inv(J)*PHI(Y)|**2, where PHI(Y) is the residual function evaluated
c     at `Y' and J is the Jacobian of the residual function evaluated at `Y',
c     or possibly at a previous value of the iterate if a Fixed Jacobian
c     iteration is employed. See [Ascher, Mattheij, and Russell, 1988, 
c     Pg. 335]. The intermediate quantities that arise during the calculation
c     of `g', namely, `PHI = PHI(Y)' and `delta = inv(J)*PHI(Y)', are returned
c     because they can be useful during the next step of the iteration. 
c     The mesh of `Nsub' subintervals is stored in `mesh'; `neqns' is the 
c     number of differential equations and the factored Jacobian is provided 
c     in the arrays `top', `blocks', `bot', and `pivot'. The number of boundary
c     conditions imposed at the left endpoint is `leftbc'. The stages values 
c     generated by the `resid' routine are stored in the `k_discrete'
c     array are passed out of this routine for use by other parts of the code.
c     `work' is a work array passed to the `resid' routine and used here for
c     temporary storage. The array, `work', is a work array.
c     The subroutines, `fsub', and `gsub', are provided by the user to define 
c     the boundary value ODE and the boundary conditions.
c------------------------------------------------------------------------------
c***declaration of constants***
      integer Mxs 
      parameter(Mxs=10)
c
c     `Mxs' Maximum value for the total number of stages of the RK method.
c
c     We wish to be able to detect when the value of the criterion
c     function will overflow. Thus we set the overflow guard to 
c     sqrt(8.5d37)=9.22d18. 
c
      double precision overflow_guard
      parameter(overflow_guard=9.22d18)
c
c***declaration of parameters***
c   imports:
      integer                 neqns,leftbc,Nsub
      double precision        mesh(0:Nsub)
      double precision        Y(neqns*(Nsub+1))
c
      double precision        top(neqns*leftbc)
      double precision        blocks(2*neqns**2 * Nsub)
      double precision        bot(neqns*(neqns-leftbc))
      integer                 pivot(neqns*(Nsub+1))
c
c     `neqns' Number of differential equations
c     `leftbc' Number of boundary conditions at the left 
c              end of the problem interval
c     `Nsub' Number of subintervals into which the problem
c            interval is partitioned
c     `mesh' Set of points which partition the problem interval
c     `Y'    Current discrete approximation to the solution, 
c            evaluated at each meshpoint. 
c     `top' Derivative information for left boundary conditions, 
c           in factored form.
c     `blocks' Derivative information for residual function, in factored form.
c     `bot' Derivative information for right boundary conditions,
c           in factored form. 
c     `pivot' Contains pivoting information associated with the factorization 
c             of the Newton matrix.  Together `top', `blocks', `bot', and 
c             `pivot' contain the almost block diagonal matrix representing 
c             the coefficient matrix of the Newton system in factored form 
c             subsequent to the call to `COLROW'.
c   exports:
      double precision      PHI(neqns*(Nsub+1)), delta(neqns*(Nsub+1))
      double precision      g, k_discrete(Mxs*neqns*Nsub)
c
c     `PHI' Value of the residual function evaluated at `Y' .
c     `delta' Solution to the linear system with the usual Newton matrix as its
c             coefficient matrix and the righthand side, `PHI'. That is,
c             `delta = inv(J)*PHI', where J is the Jacobian of the residual 
c             function.
c     `g' Value of natural criterion function at `Y'.
c                        i.e. g = 0.5*|inv(J)*PHI(Y)|^2
c     `k_discrete' Vector containing the discrete stages for
c                  all subintervals. The ith set of `s*neqns' 
c                  locations of `k_discrete' contains the `s'
c                  stages, each of length `neqns' corresponding
c                  to the ith subinterval. Passed from `resid'.
c
c***declaration of work space***
      double precision      work(neqns*(Nsub+1))
c
c     `work' Work array passed to `resid' and also used here for 
c            a temporary copy of the contents of `PHI' provided
c            to `crslve' because `crslve' will overwrite the
c            contents of the `PHI' vector, which we need to export.
c
c***user-supplied subroutines***
      external fsub,gsub
c
c     `fsub' Defines f(t,y) for first order system of differential equations.
c     `gsub' Defines the boundary condition equations.
c
c***declaration of variables for /IOCNTRL/ common block***
      integer   print_level_0, print_level_1, print_level_2
      parameter (print_level_0 = 0, print_level_1 = 1,
     *              print_level_2 = 2)
      integer                 profile
      common /IOCNTRL/        profile
c
c     `print_level_0' No output.
c     `print_level_1' Intermediate output.
c     `print_level_2' Full output.
c     `profile' Controls output, to standard output, of profiling infor-
c               mation such as Newton iteration counts, mesh selection,
c               relative defect estimate. This variable is identical to
c               the `output_control' parameter.
c----------------------------------------------------------------------------  
c      called by: `fixed_jacob', `damp_factor'
c      calls to: `resid', `dcopy', `crslve', `idamax' 
       integer                  idamax
c-----------------------------------------------------------------------------
c***evaluate natural criterion function at Y ****
c     First we compute PHI(Y) by calling `resid'; the value will
c     be returned in the vector `PHI'.
c
      call resid(neqns,leftbc,Nsub,mesh,Y,PHI, 
     +                                 k_discrete,work,fsub,gsub)
c
c     We next compute `delta = inv(J)*PHI(Y)' by solving 
c     the corresponding linear system .  J, stored in `top',
c     `blocks', `bot', and `pivot' has already been 
c     factored so we can simply call `crslve' to get the solution 
c     corresponding to the right hand side, PHI(Y), stored in `PHI'.
c     Since `crslve' would overwrite `PHI', we actually provide
c     only a copy of `PHI' to `crslve'.
c
      call dcopy((Nsub+1)*neqns,PHI,1,work,1)
c
      call crslve(top, leftbc, neqns, blocks,
     +     neqns,2*neqns,Nsub,bot,neqns-leftbc,pivot,work,delta)
c
c     We complete this calculation by computing the value of the 
c     natural criterion function at `Y'. In terms of the quantities we
c     have just computed, this is `g = 0.5*|delta|**2'.
c
      g =  dabs(delta(idamax((Nsub+1)*neqns,delta,1)))
c
c     Here we insert a test here to guard against overflow. If the 
c     norm of the Newton correction is large enough that its square 
c     overflows, the iteration is in difficulty anyway. 
c     We signal difficulty by setting g to -1.0. The negative g 
c     condition will be trapped by the calling routine, and the 
c     iteration will be terminated.
c
      if (g .GT. overflow_guard) then
                  g = -1.0d0
                  if (PROFILE .GT. PRINT_LEVEL_1) then
                      write(6,*)'Natural criterion function overflow'
                  end if
      else
                  g = 0.5d0 * g**2
                  if (PROFILE .GT. PRINT_LEVEL_2) then
                        write(6,*)'Natural criterion function value',g
                  end if
      endif
c
      return
      end


        SUBROUTINE CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,X)
C
C***************************************************************
C
C  C R S L V E  SOLVES THE LINEAR SYSTEM
C                       A*X = B
C  USING THE DECOMPOSITION ALREADY GENERATED IN  C R D C M P.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY  ...
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         THE PIVOT VECTOR FROM  C R D C M P
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C                    X - DOUBLE PRECISION(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C          
C	             X - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR
C
C***************************************************************
C
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,X,B
        DOUBLE PRECISION DOTPRD,XJ,XINCRJ,BINCRJ,SWAP
        INTEGER PIVOT(1)
        DIMENSION TOPBLK(NRWTOP,1),ARRAY(NRWBLK,NCLBLK,1),
     *          BOTBLK(NRWBOT,1),B(1),X(1)
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        NRWTP1 = NRWTOP+1
        NRWBK1 = NRWBLK+1
        NVRLP1 = NOVRLP+1
        NRWTP0 = NRWTOP-1
        NRWBT1 = NRWBOT+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
        NBLKS1 = NBLOKS+1
        NBKTOP = NRWBLK+NRWTOP
C
C***************************************************************
C
C               ****  FORWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST, IN TOPBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 130 J = 1,NRWTOP
           X(J) = B(J)/TOPBLK(J,J)
           IF(J.EQ.NRWTOP)GO TO 120
              XJ = -X(J)
              LOOP = J+1
              DO 110 I = LOOP,NRWTOP
                 B(I) = B(I)+TOPBLK(I,J)*XJ
110           CONTINUE
120        CONTINUE
130     CONTINUE
C
C       ********************************************************
C
C          ***  IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCR = 0
        DO 280 K = 1,NBLOKS
           INCRTP = INCR+NRWTOP
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 220 J = 1,NRWTOP
              INCRJ = INCR+J
              XINCRJ = -X(INCRJ)
              DO 210 I = 1,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
210           CONTINUE
220        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 240 J = NRWTP1,NRWBLK
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 225
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
225           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 230 I = LOOP,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BINCRJ
230           CONTINUE
240        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 270 J = NRWBK1,NBKTOP
              INCRJ = INCR+J
              JRWTOP = J -NRWTOP
              X(INCRJ) = B(INCRJ)/ARRAY(JRWTOP,J,K)
              IF(J.EQ.NBKTOP)GO TO 260
                 XINCRJ = -X(INCRJ)
                 LOOP = J-NRWTP0
                 DO 250 I = LOOP,NRWBLK
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
250              CONTINUE
260           CONTINUE
270        CONTINUE
           INCR = INCR+NRWBLK
280     CONTINUE
C
C       ********************************************************
C
C          ***  FINALLY, IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCRTP = INCR+NRWTOP
        DO 320 J = 1,NRWTOP
           INCRJ = INCR+J
           XINCRJ = -X(INCRJ)
           DO 310 I = 1,NRWBOT
              INCRI = INCRTP+I
              B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
310        CONTINUE
320     CONTINUE
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 350
           DO 340 J = NRWTP1,NVRLP0
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 325
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
325           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 330 I = LOOP,NRWBOT
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*BINCRJ
330           CONTINUE
340        CONTINUE
350     CONTINUE
C
C***************************************************************
C
C               ****  BACKWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 430 LL = 1,NRWBOT
           J = NVRLP1-LL
           INCRJ = INCR+J
           NRWBTL = NRWBT1-LL
           X(INCRJ) = B(INCRJ)/BOTBLK(NRWBTL,J)
           IF(LL.EQ.NRWBOT)GO TO 420
              XINCRJ = -X(INCRJ)
              LOOP = NRWBOT-LL
              DO 410 I = 1,LOOP
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
410           CONTINUE
420        CONTINUE
430     CONTINUE
C
C       ********************************************************
C
C          ***  THEN IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 490 L = 1,NBLOKS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           K = NBLKS1-L
           INCR = INCR-NRWBLK
           DO 450 L1 = NRWEL1,NRWBLK
              I = NRWBLK+NRWEL1-L1
              IPLUSN = I+NRWTOP
              LOOP = IPLUSN+1
              INCRN = INCR+IPLUSN
              DOTPRD = X(INCRN)
              DO 440 J = LOOP,NCLBLK
                 INCRJ = INCR+J
                 DOTPRD = DOTPRD-ARRAY(I,J,K)*X(INCRJ)
440           CONTINUE
              X(INCRN) = DOTPRD
              IPVTN = PIVOT(INCRN)
              IF(INCRN.EQ.IPVTN)GO TO 445
                 SWAP = X(INCRN)
                 X(INCRN) = X(IPVTN)
                 X(IPVTN) = SWAP
445           CONTINUE
450        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           INCRTP = INCR+NRWTOP
           DO 460 J = NRWBK1,NCLBLK
              INCRJ = INCR+J
              XINCRJ = -X(INCRJ)
              DO 455 I = 1,NROWEL
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
455           CONTINUE
460        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 480 LL = 1,NROWEL
              J = NRWBK1-LL
              INCRJ = INCR+J
              NRWELL = NRWEL1-LL
              X(INCRJ) = B(INCRJ)/ARRAY(NRWELL,J,K)
              IF(LL.EQ.NROWEL)GO TO 470
                 XINCRJ = -X(INCRJ)
                 LOOP = NROWEL-LL
                 DO 465 I = 1,LOOP
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
465              CONTINUE
470           CONTINUE
480        CONTINUE
490     CONTINUE
C
C       ********************************************************
C
C          ***  IN TOPBLK FINISH WITH....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 520 L = 1,NRWTOP
           I = NRWTP1-L
           LOOP = I+1
           DOTPRD = X(I)
           DO 510 J = LOOP,NOVRLP
              DOTPRD = DOTPRD-TOPBLK(I,J)*X(J)
510        CONTINUE
           X(I) = DOTPRD
           IPVTI = PIVOT(I)
           IF(I.EQ.IPVTI)GO TO 515
                 SWAP = X(I)
                 X(I) = X(IPVTI)
                 X(IPVTI) = SWAP
515        CONTINUE
520     CONTINUE
        RETURN
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine damp_factor(neqns,leftbc,Nsub,mesh,Y,delta_0,
     +           g_0,top,bot,blocks,pivot,lambda_min,lambda,
     +           PHI,delta,g,Fixed_Jacobian,info,
     +           k_discrete,Y_0,work,fsub,gsub)
c
c***overview***
c
c           This routine takes a discrete solution approximation, `Y', provided
c     on a mesh of `Nsub' subintervals (contained in `mesh') and uses 
c     an iteration to produce a suitable value for the damping 
c     factor `lambda', for use in a damped Newton update of `Y'. The number 
c     of differential equations is `neqns', while the number of boundary 
c     conditions at the lefthand endpoint of the problem interval is `leftbc'. 
c
c           With the iterate `Y' stored in `Y_0' and the corresponding Newton
c     correction input in `delta_0', this routine predicts a new value 
c     for `lambda' and then computes a trial iterate value 
c     `Y = Y_0 - lambda*delta_0'. Its acceptability is measured by 
c     whether or not it produces a sufficient reduction in the natural
c     criterion function `g = 0.5*|delta|^2', where `delta = inv(J)*PHI', 
c     where J is evaluated at `Y_0', the previous iterate, and PHI is the 
c     value of the residual function at the trial iterate `Y'. The arrays 
c     `top', `bot', `blocks', contain the derivative information associated 
c     with the Jacobian matrix J and `pivot' contains the pivoting information 
c     for the linear system solution. The new `g' value is compared with the 
c     previous value stored in `g_0'. If the reduction is not sufficient, a new
c     damping factor is chosen, if possible, and then this evaluation
c     step is repeated.  If it is found that it is impossible to find a 
c     suitable damping factor, this is interpreted as a signal that the 
c     entire Newton iteration will not converge. This is indicated by 
c     setting `info' = 3. `lambda_min' is a minimum value for the
c     damping factor. A successful return is indicated by `info'= 0. A large
c     natural criterion function value, `g', signifies difficulty in the
c     Newton iteration and is indicated by `info' = 3.
c
c           At the end of a successful search for a damping factor, it may be
c     appropriate to take one or more fixed_Jacobian steps. The criterion is 
c     that the current damped Newton step was a full Newton step i.e.
c     `lambda' = 1. This is signaled by setting `Fixed_Jacobian' to TRUE.

c           A by-product of the call to `resid' is the calculation of 
c     intermediate quantities associated with the underlying Runge-Kutta method
c     called stages. These are returned in the array called `k_discrete'.
c     This array is exported through the parameter list for later use
c     during the construction and evaluation of the interpolant associated
c     with the Runge-Kutta method. The array, `work', is a work array. 
c     The subroutines, `fsub', and `gsub', are provided by the user to define 
c     the boundary value ODE and the boundary conditions.
c------------------------------------------------------------------------------
c***declaration of constants***
      integer     Mxs
      parameter   (Mxs=10)
c
c     `Mxs' Maximum number of stages of the RK method.
c
c***declaration of parameters***
c   imports:
      integer            neqns, leftbc, Nsub
      double precision   mesh(0:Nsub), Y(neqns*(Nsub+1))
      double precision   delta_0(neqns*(Nsub+1)), g_0
      double precision   top(neqns*neqns), bot(neqns*neqns)
      double precision   blocks(2*neqns**2*Nsub)
      integer            pivot(neqns*(Nsub+1))    
      double precision   lambda_min, lambda
c
c     `neqns' Number of differential equations.
c     `leftbc' Number of boundary conditions at the left 
c              end of the problem interval.
c     `Nsub' Number of subintervals into which the problem
c            interval is partitioned.
c     `mesh' Set of points which partition the problem interval. 
c     `Y' On input, the current discrete approximation to the solution. 
c     `delta_0' Full Newton correction, equal to `inv(J)*PHI', 
c               where J and PHI are evaluated at the input iterate, `Y'.
c     `g_0' Value of the natural criterion function at the input iterate.
c     `top' Storage of derivative info for left boundary conditions.
c     `bot' Storage of derivative information for right boundary conditions.
c     `blocks' Storage of derivative information for residual function.
c     `pivot' Pivoting information during linear system solution.
c             Together `top', `blocks', and `bot' contain
c             the almost block diagonal matrix representing the 
c             coefficient matrix of the Newton system, in factored form.
c     `lambda_min' Minimum value for the damping factor. If the
c                  algorithm prescribes a value for lambda below this
c                  minimum, it is an indication that the Jacobian is
c                  "effectively singular" and that the Newton iteration 
c                  will not converge.
c     `lambda' On input, initial guess for the value of the damping 
c              factor used to perform the damped Newton update. 
c 
c   exports:
c     double precision   Y(neqns*(Nsub+1)), lambda
      double precision   PHI(neqns*(Nsub+1)), delta(neqns*(Nsub+1)), g
      logical            Fixed_Jacobian
      integer            info
      double precision   k_discrete(Mxs*neqns*Nsub)
c
c     `Y' On output, the updated discrete approximation to the 
c         solution, evaluated at each meshpoint.
c     `lambda' On output, the updated damping factor used to perform 
c              the damped Newton update.
c     `PHI' Value of the residual function associated with the updated iterate.
c     `delta' Newton correction vector, equal to inv(J)*PHI,
c             where PHI is evaluated at the new iterate and J is 
c             evaluated at the previous iterate.
c     `g' Value of natural criterion function at the updated iterate.
c     `Fixed_Jacobian' .TRUE. if the next step to be taken should be 
c                      a fixed Jacobian step. .FALSE. if the next step
c                      to be taken should be a damped Newton step.
c     `info' Communication flag: 
c                        0 - successful iteration converged to within user 
c                            specified tolerance.
c                        3 - unsuccessful termination - it was impossible to 
c                            obtain a suitable damping factor for the Newton 
c                            update (indicative of an "effectively singular" 
c                            Jacobian) or evaluation of natural criterion 
c                            function overflowed
c     `k_discrete' Discrete stages for all subintervals. The i-th set of 
c                  `s*neqns' locations of `k_discrete' contains the `s' stages,
c                  each of length `neqns' corresponding to the i-th subinterval.
c                  Computed by the `resid' routine within the call to 
c                  `criterion'.
c   work space:
      double precision   Y_0(neqns*(Nsub+1))
      double precision     work(neqns*(Nsub+1))

c     `Y_0' Current iterate value is saved here during the search
c           for a suitable damping factor and new iterate value.
c     `work' Passed to the `criterion' routine.
c
c***user-supplied subroutines***
      external fsub,gsub
c
c     `fsub' Defines f(t,y) for first order system of differential equations.
c     `gsub' Defines the boundary condition equations.
c 
c***declaration of local variables***
      logical            accept
      double precision   sigma, tau
c
c     `accept' Used to monitor the search for a new damping factor.
c     `sigma, tau'  Parameters used to the control selection
c                   of the damping factor - see descriptions below.
c
c***declaration of variables for /IOCNTRL/ common block***
c   imports:
      integer   print_level_0,print_level_1,print_level_2
      parameter (print_level_0 = 0,print_level_1 = 1,
     +              print_level_2 = 2)
      integer   profile
      common /IOCNTRL/ profile
c
c     `print_level_0' No output.
c     `print_level_1' Intermediate output.
c     `print_level_2' Full output.
c     `profile' Controls output, to standard output, of profiling infor-
c               mation such as Newton iteration counts, mesh selection,
c               relative defect estimate. This variable is identical to
c               the `output_control' parameter.
c------------------------------------------------------------------------------
c     Called by: `damped_newt'
c     Calls to : `criterion', `dcopy', `daxpy'
c------------------------------------------------------------------------------
c***initialization***
c     Two parameters which control the search for a suitable damping
c     factor are defined here:
c
c     (a)  `sigma' ensures that the reduction in the size of the natural 
c           criterion function `g', evaluated at the new iterate, will be 
c           nonnegligible. This keeps the iteration from stalling.
c 
c     (b)  `tau' controls how much the values for `lambda' are allowed to
c           change from one step to the next in the iteration performed by 
c           this routine. 
      sigma = 0.01d0
      tau = 0.1d0
c
c     During the calculation of trial iterates for the determination
c     of a suitable damping factor, `Y' will be overwritten, so  
c     save the current iterate, stored in `Y', in the vector, `Y_0'.
      call dcopy((Nsub+1)*neqns,Y,1,Y_0,1)
c
c***iterative selection of new damping factor***
c
c     Use a repeat-until loop to control the search for the new
c     damping factor. The loop will be controlled by the flags,
c     `accept' and `info': `accept' = .TRUE. or info .NE. 0
c     will terminate the loop. (Implement the loop in FORTRAN 
c     using an if-goto pair.)
c
      accept = .FALSE.
c     REPEAT
 2                                                         continue
c
c     ***compute trial iterate using current damping factor*** 
c
c            Compute an updated iterate value, stored in `Y', using the
c            current value for the damping factor, stored in `lambda'
c            `Y = Y_0 - lambda*delta_0'.
c
             call dcopy((Nsub+1)*neqns,Y_0,1,Y,1)
             call daxpy((Nsub+1)*neqns,-lambda,delta_0,1,Y,1)
c
c     ***evaluate natural criterion function at this new iterate***
c
c            Measure the suitability of the current value of the
c            damping factor by seeing how much the trial iterate
c            reduces the value of the natural criterion function. Compute
c            the value of the criterion function at `Y' by calling the
c            `criterion' routine. The value is returned in `g'; the
c            residual function, evaluated at `Y' is returned in `PHI'.
c            The value `inv(J)*PHI', (with `J' evaluated at the previous
c            iterate), is returned in `delta'.
c
             call criterion(neqns,leftbc,Nsub,mesh,Y,top,
     +         blocks, bot, pivot,PHI,delta,g,k_discrete,work,fsub,gsub)
c
c            If the natural criterion function value, `g', has overflowed
c            the `criterion' routine will trap that condition and set
c            `g = -1.0' ( a value which cannot otherwise arise).
c            Since the natural criterion function is related to the
c            size of the Newton correction, a large `g' value signifies
c            difficulty in the Newton iteration. Set `info'=3 and exit.
c
             if (g .LT. 0.0d0) then
                 info = 3
             else
c               (g >= 0.0d0)
c                Natural criterion function has not overflowed.
c
c            ***accept current damping factor or compute a new one,
c                                    depending on criterion function value***
c
c                Check value of `g' to ensure that the reduction in size is
c                sufficient. If so `lambda' and thus `Y' are acceptable.
c                If not, select a new `lambda', if possible.
c
                 if (g .LE. (1.0d0 - 2.0d0*lambda*sigma)*g_0) then
                   accept = .TRUE.
c
c                  This damped Newton step is successful. The new iterate
c                  value has already been stored in `Y' so no update is
c                  needed at this point.
c
c                  Decide whether or not to take a fixed Jacobian step 
c                  next time.
c
                   if (lambda .EQ. 1.0d0) Fixed_Jacobian = .TRUE.
                 else
c                  `(g > (1.0d0-2.0d0*lambda*sigma)*g_0')
c                  Compute a new damping factor.
                   lambda = max(tau*lambda, 
     +                     (lambda**2 * g_0)/((2*lambda-1.0d0)*g_0+g))
c
                   if (profile .GT. print_level_2) then
                     write(6,*)'Damping factor lambda=', lambda
                   endif
c
c                  Check this new `lambda' value. If it is too small, then 
c                  the Newton iteration will not converge, set `info' = 3
c                  to signal this.
c
                   if (lambda .LT. lambda_min) info = 3
                 endif
c                End of `if (g.LE.(1.0d0-2.0d0*lambda*sigma)*g_0)'
             end if
c            End of `if (g.LT.0.0d0)'
c
c     UNTIL `(accept .OR. (info .GT. 0))'
c
                                  if ((.NOT. accept) .AND.
     +                              (info .LE. 0)) go to 2
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine damped_newt(neqns,leftbc,Nsub,mesh,Y,newtol,
     +    lambda,PHI,top,bot,blocks,pivot,Fixed_Jacobian,
     +    convrg,delta,delta_0_norm,g,info,k_discrete,delta_0,
     +    work,fsub,gsub,dfsub,dgsub)
c
c***overview***
c
c            This routine takes a mesh of `Nsub' subintervals, partitioning
c     the problem interval and a discrete initial approximation provided
c     at the meshpoints, (contained in the vectors `Y' and `mesh',
c     respectively) and  uses a full or damped Newton step to refine a 
c     solution approximation for the associated discrete nonlinear system. 
c     Convergence occurs if the Newton correction is less than the 
c     tolerance `newtol' and is indicated by `convrg' = TRUE. The 
c     number of differential equations is `neqns', while the number of 
c     boundary conditions at the lefthand endpoint of the problem interval 
c     is `leftbc'. A successful return is indicated by `info'= 0. If this 
c     routine encounters a singular matrix, it signals this be setting
c     `info' = 2. A large natural criterion function value, `g', signifies 
c     difficulty and is indicated by `info' = 3. On the other hand, if the 
c     iteration appears to be converging well, this routine will signal that 
c     the iteration can switch over to a fixed Jacobian step (see below).
c
c            At the beginning of each damped Newton step, it is assumed that 
c     the current iterate value, `Y', and the value of the residual function,
c     `PHI', evaluated at `Y' are available. If the previous step was also
c     a damped Newton step, this routine also expects `delta = inv(J)*PHI',
c     where J, the Jacobian, is evaluated at the previous iterate and 
c     `delta_0_norm', the norm of the previous value of `delta_0 = inv(J)*PHI',
c     where both J and `PHI' are associated with the previous iterate. The 
c     step begins with the evaluation of then Jacobian of the residual function
c     at the current iterate, `Y'. The Newton system, consisting of the 
c     Jacobian as the coefficient matrix and `PHI', as the right hand side, is
c     then solved to obtain the full Newton correction, which we store in the 
c     vector `delta_0'. The factored Jacobian is available in the vectors 
c     `top', `bot', `blocks'. Pivoting information for the solution of the 
c     linear system is stored in `pivot'. During the call to `damp_factor' the
c     current iterate is copied into `Y_0', a damping factor, `lambda' is 
c     computed, and then a damped Newton step is taken to produce a new 
c     iterate, `Y = Y_0 - lambda*delta_0'. 
c
c           At the end of a successful damped Newton step, it may be 
c     appropriate to take one or more fixed Jacobian steps . The criterion 
c     is that the current damped Newton step was a full Newton step 
c     (`lambda'=1). This is signaled by setting `Fixed_Jacobian' to TRUE. 
c     The fixed Jacobian step, for reasons associated with efficient data
c     management, expects to receive a number of quantities : (i) the new 
c     iterate, `Y', (ii) the residual function evaluated at this new iterate, 
c     `PHI', (iii) `delta = inv(J)*PHI', where J, the Jacobian, has been
c     evaluated rather at the previous iterate  and (iv) the value of the 
c     natural criterion function, `g = 0.5*|delta|**2'. All these values are
c     to be naturally computed during the current damped Newton step and 
c     thus will be already set up for the fixed Jacobian step if appropriate. 
c
c           A by-product of the evaluation of the right hand side of the
c     Newton system is the calculation of intermediate quantities associated
c     with the underlying Runge-Kutta method called stages. These are
c     returned to the `damped_newt' routine when it calls the `damp_factor'
c     routine and are stored in the array called `k_discrete'.
c     This array is exported through the parameter list for later use
c     during the construction and evaluation of the interpolant associated
c     with the Runge-Kutta method. The array, `work', is a work array provided 
c     to the `damp_factor' and `newmat' routines. The subroutines, `fsub', 
c     `gsub', `dfsub', and `dgsub', are provided by the user to define the 
c     boundary value ODE and the boundary conditions and their derivatives.
c------------------------------------------------------------------------------ 
c***declaration of constants***
      integer    Mxs
      parameter  (Mxs=10)
c
c     `Mxs'    is the maximum number of stages of the RK method.
c
c     We wish to be able to detect when the criterion function value 
c     will overflow. Thus the overflow_guard to 9.22d18.
c
      double precision       overflow_guard
      parameter   (overflow_guard = 9.22d18)
c
c***declaration of parameters***
c   imports:
      integer            neqns, leftbc, Nsub
      double precision   mesh(0:Nsub), Y(neqns*(Nsub+1))
      double precision   newtol, lambda, PHI(neqns*(Nsub+1))
      double precision   delta(neqns*(Nsub+1)), delta_0_norm
c
c     `neqns' Number of differential equations.
c     `leftbc' Number of boundary conditions at the left 
c              end of the problem interval.
c     `Nsub' Number of subintervals into which the problem
c            interval is partitioned.
c     `mesh' Set of points which partition the problem interval.
c     `Y' On input, the current discrete approximation to the solution.
c     `newtol' Value to be applied to the Newton iteration,
c              i.e. the iteration stops if the Newton correction is
c              less than or equal to the user specified tolerance.
c              A blended relative/absolute tolerance is applied.
c              If `delta_j' is the Newton correction for, `Y_j', the jth 
c              component of `Y', then we require
c                    `(|delta_j| / |Y_j|) <= newtol'.
c     `lambda' On input, the damping factor used in the previous damped
c              Newton step. If its value is set to 0 on input, the
c              previous step was not a damped Newton step.
c     `PHI' On input, the value of the residual function on each
c           subinterval evaluated at the current iterate `Y'.
c     `delta' On input, it is `inv(J)*PHI' where J is evaluated at the
c             previous iterate and `PHI' is evaluated at the input
c             iterate, if the previous step was a damped Newton step.
c     `delta_0_norm' On input, is |inv(J)*PHI| where J and `PHI' were
c                    evaluated at the previous iterate if the 
c                    previous step was a damped Newton step.
c                    It is undefined if the previous step was not a 
c                    damped Newton step.
c   exports:
c     double precision   Y(neqns*(Nsub+1)), lambda
c     double precision   PHI(neqns*(Nsub+1))
      double precision   top(neqns*neqns)
      double precision   bot(neqns*neqns)
      double precision   blocks(2*neqns**2*Nsub)
      integer            pivot(neqns*(Nsub+1)) 
      logical            Fixed_Jacobian, convrg
c     double precision   delta(neqns*(Nsub+1))
c     double precision   delta_0_norm
      double precision   g
      integer            info
      double precision   k_discrete(Mxs*neqns*Nsub)
c
c     `Y' On output, the updated discrete approximation to the solution.
c     `lambda' On output, the value for the new damping factor.
c     `PHI' On output, the value of the residual function on each
c           subinterval evaluated at the updated iterate `Y'.
c     `top' Storage of derivative information for left boundary 
c           conditions.
c     `bot' Storage of derivative information for right boundary
c           conditions.
c     `blocks' Storage of derivative information for residual function.
c     `pivot' Pivoting information during linear system solution.
c             Together `top', `blocks', and `bot' contain
c             the almost block diagonal matrix representing the 
c             coefficient matrix of the Newton system, in factored form.
c     `Fixed_Jacobian' .TRUE. if the next step to be taken should be
c                      a fixed Jacobian step. .FALSE. if the next step
c                      to be taken should be a damped Newton step.
c     `convrg' Used to check for convergence of the Newton iteration. 
c     `delta' Newton correction vector, equal to `inv(J)*PHI(Y)', 
c             where J is evaluated at the input iterate and `PHI'
c             and `PHI' is evaluated at the new iterate `Y'. 
c     `delta_0_norm' on output, |inv(J)*PHI| where J and `PHI' are
c                    evaluated at the input iterate.
c     `g' Value of natural criterion function at the updated iterate.
c     `info' Internal communication flag: 
c            0 - successful iteration converged to within user 
c                specified tolerance.
c            2 - unsuccessful termination - a singular coefficient 
c                matrix was encountered during the attempted solution 
c                of the Newton system.
c            3 - unsuccessful termination - it was impossible to obtain 
c                a suitable damping factor for the Newton update
c                (indicative of an "effectively singular" Jacobian) or 
c                evaluation of natural criterion function overflowed.
c     `k_discrete' Vector containing the discrete stages for all 
c                  subintervals. The i-th set of `s*neqns' locations of 
c                  `k_discrete' contains the `s' stages, each of length
c                  `neqns' corresponding to the i-th subinterval.
c   work space:
      double precision   delta_0(neqns*(Nsub+1))
      double precision   
     +                  work(neqns*(2*Nsub+3+2*(Mxs+1)*neqns))
      integer            i_Y_0, i_work
                                                                       
c     `delta_0' Full Newton correction, equal to `inv(J)*PHI', 
c               where J and `PHI' are evaluated at the input iterate.
c     `work' Work array provided to the `damp_factor' and
c            `newmat' routines.
c     `i_Y_0' Index to the start of the part of `work' to be passed 
c             to `damp_factor' for use as the work array `Y_0'.
c     `i_work' Index to the start of the remainder of `work' to be passed
c              to `damp_factor' for use as the work array `work'.
c
c***user-supplied subroutines***
      external fsub,gsub,dfsub,dgsub
c
c     `fsub' Defines f(t,y) for the first order
c            system of differential equations, y' = f(t,y).
c     `gsub' Defines the boundary condition equations.
c     `dfsub' Defines the Jacobian of the system
c             of differential equations.
c     `dgsub' Defines the Jacobian of the boundary conditions. 
c
c***declaration of local variables***
      integer            factor
      double precision   g_0, lambda_min,mu
      integer            j
c
c     `factor' Used by the `colrow' routine to indicate trouble; a
c              return with `factor' equal to -1 indicates that a singular
c              Jacobian was encountered.      
c     `g_0' Value of the natural criterion function at the input iterate.
c     `lambda_min' Parameter used to the control selection of the 
c                  damping factor - see descriptions below.
c     `mu' Predicted value for the new damping factor.
c     `j' Loop index, 1,...,`(Nsub+1)*neqns'.
c
c***declaration of variables for /IOCNTRL/ common block***
c   imports:
      integer   print_level_0,print_level_1,print_level_2
      parameter (print_level_0 = 0,print_level_1 = 1,
     +               print_level_2 = 2)
      integer   profile
      common /IOCNTRL/ profile
c
c     `print_level_0' No output.
c     `print_level_1' Intermediate output.
c     `print_level_2' Full output.
c     `profile' Controls output, to standard output, of profiling infor-
c               mation such as Newton iteration counts, mesh selection,
c               relative defect estimate. This variable is identical to
c               the `output_control' parameter.
c------------------------------------------------------------------------------
c     Called by : `NewIter'
c     Calls to  : `NewMat', `colrow', `damp_factor', `daxpy', `idamax'
               integer   idamax
c------------------------------------------------------------------------------
c***initialization of work array indexes***
      i_Y_0 = 1
      i_work = neqns*(Nsub+1)+1

c***perform a damped newton step
c     The flag `info' is used to monitor for trouble; It is initialized to 0.
      info = 0
c
c     `lambda_min' is a parameter which controls the search for a
c     suitable damping factor; it gives the minimum value
c     for the damping factor. It is used here in the initial guess for the 
c     damping factor. It is also used with the `damp_factor' routine.  
      lambda_min = 0.01d0
c
c     First construct the Newton system at the input iterate.
c     It is assumed that the value of  the residual function
c     at the current iterate was obtained at the end of the previous 
c     iteration step and is available in the vector `PHI'.  The  
c     coefficient  matrix  is  computed  by calling the `NewMat' 
c     routine and is returned in the  vectors  `top', `blocks', `bot'.
c
      call NewMat(neqns,leftbc,Nsub,mesh,Y,top,blocks,bot,
     +                               k_discrete,work,dfsub,dgsub)
c
c     The coefficient matrix constructed by `NewMat' has an 
c     "almost block diagonal" structure. `colrow', a package of
c     subroutines designed to solve such systems is now employed.
c     This routine will return with `factor' = -1 if  the coefficient
c     matrix of the Newton system is singular and `factor'=0 otherwise.
c     The solution of the linear system, which gives the Newton correction, 
c     is returned in the vector `delta_0'.
c
      call colrow((Nsub+1)*neqns,top,leftbc,neqns,blocks,
     +              neqns,2*neqns,Nsub,bot,neqns-leftbc,
     +              pivot,PHI,delta_0,factor)
c
      if ( factor .EQ. -1 ) then
            info = 2
c           This is a signal that a singular Jacobian matrix was detected.
      end if
c
c
      if (info .EQ. 0) then
c
c          `colrow' was successful in solving the linear system to obtain
c          the Newton correction. Attempt to compute the corresponding
c          new iterate using a damped Newton step. At this time, the value
c          of the natural criterion function for the current iterate is
c          also computed and stored in `g_0'. (It can be shown 
c          that, for this case, `g_0 = 0.5 * |delta_0|**2'.)
c
c          Check here to see if the natural criterion function 
c          value will overflow. If so, trap the condition. Since the
c          natural criterion function is related to the size of the 
c          Newton correction, a large `g' value signifies difficulty in
c          the Newton iteration. Set `info' = 3 and exit.
c
           g_0 = dabs(delta_0(idamax((Nsub+1)*neqns, delta_0,1)))
c
           if (g_0 .GT. overflow_guard) then
                 if (profile .GT. print_level_1) then
                     write(6,*)'Natural criterion function overflow'
                     info = 3
                 end if
           else
                 if (profile .GT. print_level_1) then
                     write(6,*)'Norm of Newton correction',g_0
                 endif
                 g_0 = 0.5d0 * g_0**2
           endif
c
c          Test for convergence now. If the norm of the scaled Newton
c          correction, `delta_0', is sufficiently small, then the iteration
c          has converged. 
c
           convrg = .TRUE.
           do 20 j = 1, (Nsub+1)*neqns
             if( dabs(delta_0(j)) .GT. newtol*(dabs(Y(j))+1.0d0))
     +                convrg = .FALSE.
  20       continue
c
           if (convrg) then
c              Update iterate with full Newton correction.
c              The current iterate is still available in `Y' and
c              the Newton correction is in `delta_0'.
               call daxpy((Nsub+1)*neqns,-1.0d0,delta_0,1,Y,1)
c
c              The stages of the RK method are updated in the next call to 
c              the `interp_setup' routine with this new iterate value.
c
           else
c              Attempt to update the Newton iterate with a damped Newton 
c              correction. Compute a prediction for the new damping factor. 
c              The value for `lambda' must be between `lambda_min' and 1.0d0.
c              If the input `lambda' value is 0.0d0, this is a signal that the
c              previous iteration was not a damped Newton iteration and
c              appropriate information for the prediction of a new `lambda'
c              value is not available, in which case set the predicted 
c              value to 1.0d0.
c
               if (lambda .EQ. 0.0d0 ) then
                   lambda = 1.0d0
               else
c                  For temporary storage overwrite `delta = delta - delta_0'.
                   call daxpy((Nsub+1)*neqns,-1.0d0,delta_0,1,delta,1)
c
c                  Use the norm of this difference in the prediction for 
c                  `lambda'.
                   mu = (lambda*delta_0_norm)/
     +                  dabs(delta(idamax((Nsub+1)*neqns,delta,1)))
                   lambda = max(lambda_min, min(1.0d0, mu))
               endif
c
c              Print out current value for `lambda', if appropriate.
               if (profile .GT. print_level_1) then
                 write(6,*)'Initial damping factor lambda =', lambda
               endif
c
c              The previous value of the norm of the Newton correction,
c              stored in `delta_0_norm', was just used in the prediction
c              of the new `lambda' value, and is no longer needed.          
c              Update the value of `delta_0_norm' using the current value
c              of `delta_0' computed by `colrow' above (the norm of `delta_0'
c              was already computed during the calculation of `g_0'). 
c
               delta_0_norm = dsqrt(2.0d0*g_0)
c
c              Compute a new iterate by updating the current iterate with
c              the (possibly damped) Newton correction.
c
               call damp_factor(neqns,leftbc,Nsub,mesh,Y,delta_0,
     +              g_0,top,bot,blocks,pivot,lambda_min,lambda,
     +              PHI,delta,g,Fixed_Jacobian,info,k_discrete,
     +              work(i_Y_0),work(i_work),fsub,gsub)
               
               if ((info .EQ. 0) .AND. 
     +               (profile .GT. print_level_2)) then 
                     write(6,*)'Norm of damped Newton correction',
     +                                         lambda*delta_0_norm
               endif
c
c              If successful, this routine will compute an appropriate
c              damping factor, returned in `lambda', the updated iterate, 
c              returned in `Y', the corresponding residual function value, 
c              returned in `PHI', `delta = inv(J)*PHI', and `g', the natural 
c              criterion function, evaluated at the new iterate.
c
           endif
c          End of `if (convrg)'
c
      endif
c     End of `if (info.EQ.0)'
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine defect_estimate(neqns,Nsub,mesh,Y,defect,
     +    defect_norm,info,z,z_prime,def_1,def_2,f_1,f_2,
     +    temp_1,temp_2,k_discrete,k_interp,work,fsub) 
c
c***overview***
c
c        This routine uses the discrete solution approximation, `Y', 
c     obtained on the mesh stored in `mesh', plus the stages of the underlying 
c     RK method, stored in `k_discrete', plus some new stages, stored in 
c     `k_interp', to construct an interpolant for `Y'. `neqns' is the number of
c     differential equations; `Nsub' is the the number of subintervals into 
c     which the mesh divides the problem interval. This interpolant is then used
c     to compute an estimate of the defect on each subinterval, for each 
c     solution component. The i-th group of `neqns' locations of `defect' is 
c     used to store the relative defect information associated with the i-th 
c     subinterval. The max norm of the relative defect is returned in 
c     `defect_norm'. A successful return is indicated by `info = 0'. If the 
c     `defect_norm' is bigger than the constant parameter `defect_threshold' 
c     then `info' is set to 4.     
c 
c         The primary information required to define the interpolant on the i-th
c     subinterval is the current solution approximation at the (i-1)-st and 
c     i-th mesh points (which can be obtained from `Y'), and some
c     stages associated with the i-th interval.  The first `s' stages come 
c     from the discrete formula, and have already been computed. They are 
c     available in the `s' vector locations making up the i-th set of 
c     `s*neqns' locations of the vector,'k_discrete'. The remaining 
c     `s_star - s' stages are required only for the interpolant and must 
c     now be computed . For each subinterval, this is done through a call 
c     to the `interp_setup' routine, which returns the new stage evaluations 
c     for that subinterval. The resultant stages are then stored in the 
c     appropriate locations of the `k_interp' vector (which has the same 
c     format as the `k_discrete' vector). 
c
c           On the i-th subinterval the evaluation of the interpolant at 
c     `t = mesh(i-1)+tau*hi' is given by `z(t) = yim1 + hi sum br(tau)*kr',
c     where `hi' is the subinterval size, 'yim1' is the solution approximation 
c     at `mesh(i-1)', `kr' is the r-th stage, `br(tau)' is the corresponding 
c     weight polynomial evaluated at `tau' and the sum runs from 1 to `s_star'.
c              
c            Once the interpolant  is available, we can sample 
c     the defect at the point, `mesh(i-1)+tau*hi', by substituting the
c     interpolant into the differential equation and subtracting the right 
c     side from the left side, i.e. the defect on the i-th subinterval is
c     `delta = z'(mesh(i-1)+tau*hi)-f(mesh(i-1)+tau*hi,z(mesh(i-1)+tau*hi))'.
c     A blended relative/absolute defect measure, `defect', equal to 
c     `|delta|/(|f|+1)' is computed, where `f' is the corresponding right hand
c     side of the ODE. 
c        
c            The sample point is available from the /sample_point/ common
c     block and the subsequent evaluation of the weights is done by a call
c     to the `interp_weights' routine. To improve the reliability of the 
c     defect estimate, we also sample at `1-tau' and choose the larger of
c     the two defect values for the estimate, `delta'.
c
c            `z','z_prime,'def_1','def_2','f_1','f_2','temp_1','temp_2'
c     are work arrays used for temporary storage here while `work'
c     is a work array passed to `interp_setup'. The subroutine, `fsub', 
c     is provided by the user to define the boundary value ODEs.
c------------------------------------------------------------------------------
c***declaration of constants***
      integer Mxs 
      parameter(Mxs=10)
c
c     `Mxs'    Maximum value for the total number of stages of the interpolant.
c
      double precision defect_threshold
      parameter (defect_threshold = 1.0d-1)
c   
c     `defect_threshold' Used to monitor the size of the defect.
c
c***declaration of parameters***
c   imports:
      integer neqns, Nsub
      double precision mesh(0:Nsub), Y(neqns*(Nsub+1))
      double precision k_discrete(Mxs*neqns*Nsub)
c
c     `neqns' Number of differential equations.
c     `Nsub' Number of subintervals.
c     `mesh' Array of points defining the subintervals.
c     `Y' Contains the discrete solution at the mesh points
c         the i-th set of `neqns' locations of `Y' contains 
c         the solution approximation at the i-th mesh point.
c     `k_discrete' Contains the `s' discrete stage values for each subinterval.
c
c   exports:
      double precision         defect(Nsub*neqns), defect_norm
      integer                  info
      double precision k_interp(Mxs*neqns*Nsub)
c
c     `defect' Contains the corresponding blended relative/absolute
c              defect measure.
c     `defect_norm' Max norm of the defect measure.
c     `info' Communication flag :
c                        0 - Successful return.
c                        4 - Size of the defect norm was too large and the 
c                            solution from the current discrete system can not 
c                            be trusted.
c     `k_interp' Contains the `s_star-s' extra stage values for each
c                subinterval need to form the interpolant.
c
c   work space:
      double precision z(neqns),z_prime(neqns)
      double precision def_1(neqns),def_2(neqns)
      double precision f_1(neqns),f_2(neqns)
      double precision temp_1(neqns),temp_2(neqns)
      double precision  work(neqns)

c     `z' Used to accumulate argument for fsub evaluation
c     `z_prime' Used to save derivative of interpolant
c     `def_1,def_2' Contain the values of the defect at the sample
c                   points `tau' and `1-tau' for each subinterval.
c     `f_1,f_2' Contain the corresponding fsub value at the
c               the two sample points.
c     `temp_1', `temp_2' Temporary location for storage of scaled defect 
c                        at `tau' and `1-tau'.
c     `work' Work array passed to `interp_setup'.
c
c***user-supplied subroutines***
      external fsub
c
c     `fsub' Defines f(t,y) for the first order system of differential 
c            equations, y' = f(t,y).
c
c***declaration of local variables***
      integer i,j
      double precision hi, weights_1(Mxs), weights_2(Mxs)
      double precision weights_1_prime(Mxs), weights_2_prime(Mxs)
      double precision estimate_1, estimate_2
c
c     `i' Loop index from 1 to Nsub.
c     `j' Loop index from 1 to `neqns'.
c     `hi' Size of the i-th subinterval.
c     `f' Used to store `fsub' evaluation.
c     `weights_1', `weights_2' Contains the weight polynomials of the 
c                              interpolant evaluated at the sample point `tau' 
c                              and `1-tau'.
c     `weights_1_prime', `weights_2_prime' Contains the first derivatives of 
c                                          the weight polynomials of the 
c                                          interpolant evaluated at the sample 
c                                          points `tau' and `1-tau'.
c     `estimate_1', `estimate_2' Maximum scaled defect at `tau' and `1-tau'
c
c***declaration of variables from /RK_s/ common block***
      integer s
      common /RK_s/ s
c     `s' Number of discrete stages.
c
c***declaration of variables from /interp_s_star/ common block***
      integer s_star
      common /interp_s_star/ s_star
c     `s_star' Total number of stages required to form the
c              interpolant on each subinterval. It includes all the
c              stages of the discrete formula plus the additional 
c              stages required for the interpolant.   
c
c***declaration of variables from /sample_point/ common block***
      double precision tau_star 
      common /sample_point/ tau_star
c
c     `tau_star' Relative position of one of the sample points within each 
c                subinterval. The other sample point is `1-tau_star'.
c-----------------------------------------------------------------------------
c     Called by: `mirkdc'
c     Calls to : `interp_weights', `interp_setup', `sum_stages', `fsub'
c                `daxpy', `dcopy'
      integer    idamax
c------------------------------------------------------------------------------
c     Setup the weights, evaluated at the first sample point.
c
      call interp_weights(s_star,tau_star,weights_1,weights_1_prime)
c
c     Setup the weights, evaluated at the second sample point.
c
      call interp_weights(s_star,1.0d0-tau_star,weights_2,
     +                                        weights_2_prime)
c
      do 30 i = 1, Nsub
            hi = mesh(i) - mesh(i-1)
c
c           We must first setup the new stages for this subinterval.
c
            call interp_setup(neqns,mesh(i-1),hi,Y((i-1)*neqns+1),
     +            Y(i*neqns+1),s,k_discrete((i-1)*s*neqns+1),
     +            s_star,k_interp((i-1)*(s_star-s)*neqns+1),work,fsub)
c
c           The new stages for the i-th subinterval have been computed
c           and stored in the ith set of `(s_star-s)*neqns' locations of
c           `k_interp'.
c
c     ***Sample Point 1***
c
c            With all the information needed to form the interpolant on the
c            i-th subinterval already available, we can compute the value of 
c            the interpolant, and its first derivative at the sample point
c            by calling the `sum_stages' routine.
c
             call sum_stages(neqns,hi,Y((i-1)*neqns+1),
     +                        s,k_discrete((i-1)*s*neqns+1),s_star,
     +                        k_interp((i-1)*(s_star-s)*neqns+1),
     +                        weights_1,weights_1_prime,z,z_prime)
c
c            This yields the value of the interpolant, in `z', and its first
c            derivative, in `z_prime'.   We must next evaluate the righthand
c            side of the system of ode's at `z', i.e. compute f(z) and store 
c            in `f_1'.
c
             call fsub(neqns,mesh(i-1)+tau_star*hi,z,f_1)
c
c            Since the calculation of `z_prime' was already complete after 
c            the return from `sum_stages', we can now compute the defect which 
c            we store in `def_1'.
c
             call daxpy(neqns,-1.0d0,f_1,1,z_prime,1)
             call dcopy(neqns,z_prime,1,def_1,1)
c
c            Compute the norm of the scaled defect for Sample Point 1.
c
             do 10 j =1,neqns
                     temp_1(j) = 
     +                        def_1(j)/(dabs(f_1(j))+1.0d0)
  10         continue
c
             estimate_1=dabs(temp_1(idamax(neqns,temp_1,1)))
c
c     ***Sample Point 2*** 
c
c            With all the information needed to form the interpolant on the
c            i-th subinterval already available, we can compute the value of 
c            the interpolant, and its first derivative at the sample point
c            by calling the `sum_stages' routine.
c
             call sum_stages(neqns,hi,Y((i-1)*neqns+1),
     +                       s,k_discrete((i-1)*s*neqns+1),s_star,
     +                       k_interp((i-1)*(s_star-s)*neqns+1),
     +                       weights_2,weights_2_prime,z,z_prime)
c
c            This yields the value of the interpolant, in `z', and its first
c            derivative, in `z_prime'.  We must evaluate the righthand side 
c            of the system of ode's at `z', i.e. compute f(z).
c
             call fsub(neqns,mesh(i-1)+(1.0d0-tau_star)*hi,z,f_2)
c
c            Since the calculation of `z_prime' was already complete after 
c            the return from `sum_stages', we can now compute the defect, 
c            which we store in `def_2'.
c
             call daxpy(neqns,-1.0d0,f_2,1,z_prime,1)
             call dcopy(neqns,z_prime,1,def_2,1)
c
             do 20 j =1,neqns
                  temp_2(j) = 
     +                     def_2(j)/(dabs(f_2(j))+1.0d0)
 20          continue
c
             estimate_2=dabs(temp_2(idamax(neqns,temp_2,1)))
c
c     ***Compare defect estimates from the two sample points***
c
             if (estimate_1 .GT. estimate_2) then
                call dcopy(neqns,temp_1,1,defect((i-1)*neqns+1),1)
             else
                call dcopy(neqns,temp_2,1,defect((i-1)*neqns+1),1)      
             endif
c
c            The scaled defect for the i-th subinterval has now been copied
c            into the ith group of `neqns' locations of the vector `defect'.
c            
 30   continue
c                                                        
c     Compute the max norm of the defect estimate.
      defect_norm = dabs(defect(idamax(Nsub*neqns,defect,1)))
c
c     If this defect_norm is bigger than `defect_threshold' then
c     this is an indication that, the defect is no longer a 
c     small relative perturbation of the original system of ODEs
c     and the solution from the current discrete system should
c     not be trusted. It is probably better to start over as if the
c     Newton iteration had failed. Here, we check the size of `defect_norm', 
c     and set `info' = 4 (to cause a mesh doubling), if it is too large.
      if (defect_norm .GT. defect_threshold) info = 4
c                                                           
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fixed_jacob(neqns,leftbc,Nsub,mesh,Y,newtol,
     +       delta,g,PHI,top,bot,blocks,pivot,Fixed_Jacobian,lambda,
     +       convrg,info,k_discrete,Y_hat,PHI_hat,work,fsub,gsub)
c
c***overview***
c
c           This routine takes a mesh of `Nsub' subintervals, partitioning
c     the problem interval and a discrete initial approximation provided
c     at the meshpoints, (contained in the vectors `Y' and `mesh',
c     respectively) and uses a Fixed Jacobian step, as part of a Newton 
c     iteration, to produce an updated solution approximation. Convergence
c     occurs if the Newton correction is less than the tolerance `newtol' 
c     and is indicated by `convrg' = TRUE. The number of differential equations
c     is `neqns', while the number of boundary conditions at the lefthand 
c     endpoint of the problem interval is `leftbc'. A successful return is 
c     indicated by `info'= 0. A large natural criterion function value,'g', 
c     signifies difficulty in the Newton iteration and is indicated by 
c     `info' = 3.  
c
c          When a fixed Jacobian step is to be taken, the following
c     quantities are assumed to be available from the previous iteration 
c     step  (i) the iterate, `Y', (ii) the residual function evaluated at this
c     iterate, `PHI', (iii) `delta = inv(J)*PHI', where J, the Jacobian ,
c     available in factored form in the arrays `top', `bot', `blocks', and 
c     `pivot', has not been evaluated at `Y', but rather at some previous 
c     iterate and (iv) the value of the natural criterion function, `g', at
c     the iterate, `Y'. 
c
c          The fixed Jacobian step consists of first computing  a (potential)
c     new iterate, `Y_hat = Y - delta'. The acceptability of `Y_hat' is 
c     determined  by evaluating the natural criterion function at `Y_hat', 
c     giving `g_hat'. During the calculation of `g_hat, it is necessary to 
c     evaluate the residual function at `Y_hat' and store it in the vector, 
c     `PHI_hat', and a new value for `delta = inv(J)*PHI_hat', both of which 
c     will then be available for the next step of the iteration. The 
c     acceptability of the new iterate value is determined by comparing `g_hat'
c     with the value of the natural criterion function for the previous iterate,
c     stored in `g'. If `g_hat < rho*g' then the new iterate value is accepted
c     and another fixed Jacobian step is attempted. If `g_hat > rho*g', the
c     code returns to taking damped Newton steps (the value of `rho' is set in
c     this routine). This is signaled by setting `Fixed_Jacobian' to FALSE. 
c     The value of `lambda', the damping factor, is set to 0 as well. This is
c     to indicate appropriate information for the prediction of a new damping
c     factor is not available for the upcoming damped Newton step.
c
c           A by-product of the evaluation of the right hand side of the
c     Newton system is the calculation of intermediate quantities associated
c     with the underlying Runge-Kutta method called stages. These are
c     returned in the array called `k_discrete', which is exported through 
c     the parameter list for later use during the construction and evaluation 
c     of the interpolant associated with the Runge-Kutta method. The array, 
c     `work', is a work array. The subroutines, `fsub', and `gsub', are provided
c     by the user to define the boundary value ODE and the boundary conditions.
c------------------------------------------------------------------------------ 
c***declaration of constants***
      integer  Mxs
      parameter (Mxs=10)
c
c     `Mxs'    is the maximum number of stages of the RK method.
c
c***declaration of parameters***
c   imports:
      integer            neqns, leftbc, Nsub
      double precision   mesh(0:Nsub), Y(neqns*(Nsub+1))
      double precision   newtol, delta(neqns*(Nsub+1)), g
      double precision   PHI(neqns*(Nsub+1)), top(neqns*neqns)
      double precision   bot(neqns*neqns), blocks(2*neqns**2*Nsub)
      integer            pivot(neqns*(Nsub+1))  
c
c     `neqns'   the number of differential equations
c     `leftbc'  the number of boundary conditions at the left 
c               end of the problem interval
c     `Nsub'    the number of subintervals into which the problem
c               interval is partitioned
c     `mesh'    the set of points which partition the problem interval.
c     `Y' on input, the current discrete approximation to the solution. 
c     `newtol' the tolerance value to be applied to the Newton iteration,
c              i.e. the iteration stops if the Newton correction 
c              is less than or equal to the user specified 
c              tolerance. A blended relative/absolute tolerance
c              is applied. If `delta_j' is the Newton correction
c              for the jth component of Y, Y_j, then we require
c                          |delta_j| / (|Y_j|+1) <= newtol.
c     `delta'   on input, the Newton correction vector, equal to 
c               inv(J)*PHI(Y), where J is evaluated at a previous iterate.
c     `g'       on input, the value of natural criterion function at the 
c               current iterate, `Y'.
c     `PHI'    on input, is the value of the residual function on each
c              subinterval, evaluated at `Y'.
c     `top'     storage of derivative information for left boundary conditions.
c     `bot'     storage of derivative information for right boundary conditions.
c     `blocks'  storage of derivative information for residual function.
c     `pivot'   pivoting information associated with the factorization
c               of J. Together `top', `blocks', and `bot' contain
c               the almost block diagonal matrix representing the 
c               coefficient matrix of the Newton system, in factored form.
c   exports:
c     double precision   Y(neqns*(Nsub+1)),delta(neqns*(Nsub+1))
c     double precision   g, PHI(neqns*(Nsub+1))              
      logical            Fixed_Jacobian
      double precision   lambda
      logical            convrg
      integer            info
      double precision   k_discrete(Mxs*neqns*Nsub)
c
c     `Y' On output, the updated discrete approximation to the solution.
c     `delta' On output, equal to inv(J)*PHI where the
c             Jacobian J is evaluated at the previous iterate 
c             and PHI is evaluated at the updated iterate `Y_hat'
c     `g' On output, value of the criterion function evaluated
c         at the updated iterate.
c     `PHI' On output, is value of the residual function
c           on each subinterval evaluated at the updated iterate.
c     `Fixed_Jacobian' .TRUE. if the next step to be taken should be a
c                      fixed Jacobian step. .FALSE. if the next step
c                      to be taken should be a damped Newton step.
c     `lambda' Set to 0 to indicate appropriate information for the
c              prediction of a new lambda is not available during the 
c              upcoming damped Newton step
c     `convrg' Used to check for convergence of the Newton iteration. 
c     `info' Communication flag: 
c                 0 - successful iteration converged to
c                     within user specified tolerance.
c                 3 - unsuccessful termination - A large natural
c                     criterion function value, g, signifies 
c                     difficulty in the Newton iteration.
c     `k_discrete' Discrete stages for all subintervals. The ith set of 
c                  `s*neqns' locations of `k_discrete' contains the `s'
c                  stages, each of length `neqns', corresponding
c                  to the ith subinterval. Computed by `resid'
c                  from within `criterion`.
c   work space:
      double precision    Y_hat(neqns*(Nsub+1))
      double precision    PHI_hat(neqns*(Nsub+1))
      double precision    work(neqns*(Nsub+1))

c     `Y_hat'  Trial iterate value. 
c     `PHI_hat' Residual function value associated with the trial iterate.
c     `work' Work space passed to `criterion'.
c
c***user - supplied subroutines***
      external fsub,gsub
c     `fsub' Defines f(t,y) for the first order
c            system of differential equations, y' = f(t,y).
c     `gsub' Defines the boundary condition equations.
c 
c***declaration of local variables***
      double precision   g_hat, rho
      integer            j
c
c     `g_hat' Value of the natural criterion function at the trial iterate. 
c     `rho' Used in measuring performance of fixed Jacobian iteration.
c     `j' Loop index, 1,...,(Nsub+1)*neqns.
c
c***declaration of variables for /IOCNTRL/ common block***
c   imports:
      integer   print_level_0,print_level_1,print_level_2
      parameter (print_level_0 = 0,print_level_1 = 1,
     +              print_level_2 = 2)
      integer   profile
      common /IOCNTRL/ profile
c
c     `print_level_0' No output.
c     `print_level_1' Intermediate output.
c     `print_level_2' Full output.
c     `profile' Controls output, to standard output, of profiling infor-
c               mation such as Newton iteration counts, mesh selection,
c               relative defect estimate. This variable is identical to
c               the `output_control' parameter.
c------------------------------------------------------------------------------
c     Called by : `Newiter'
c     Calls to  : `criterion', `dcopy', `daxpy', `idamax'
      integer       idamax         
c------------------------------------------------------------------------------
c     The flag `info' is used to monitor for trouble; it is initialized to 0.
      info = 0
c
c     When a fixed Jacobian step has been taken, the decision about
c     whether or not to follow with another fixed Jacobian step or not 
c     will be determined by the size of the decrease in the value of the
c     natural criterion function. This decrease must be at least a factor 
c     of `rho'. Here `rho' is initialized.
      rho = 0.5d0
c
c     In addition to the new iterate `Y' being available, we also 
c     have the corresponding evaluation of the residual function,
c     `PHI', the Newton correction, `delta', for the current fixed 
c     Jacobian step, and the evaluation of the criterion function
c     for the this iterate, `g = 0.5*|inv(J)*PHI|**2 = 0.5*|delta|**2'.
c
c     If the Newton correction, `delta', is less than the tolerance, the 
c     iteration has converged. A blended relative/absolute tolerance is 
c     employed.
c
      if (profile .GT. print_level_1) then
           write(6,*)'Norm of Newton correction',
     +       dabs(delta(idamax((Nsub+1)*neqns,delta,1)))
      endif
c
      convrg = .TRUE.
      do 10 j = 1, (Nsub+1)*neqns
        if(dabs(delta(j)).GT.newtol*(dabs(Y(j))+1.0d0)) convrg=.FALSE.
  10  continue
c
      if (convrg) then
c          Update current iterate before exiting.
           call daxpy((Nsub+1)*neqns,-1.0d0,delta,1,Y,1)
c
c          The stages of the RK method are updated in the next 
c          call to the `interp_setup' routine with this new
c          iterate value.
c
      else              
c          If the iteration has not converged, some care must be taken
c          in using the fixed Newton step `delta' to  update the
c          current iterate. The correction will not be useful if it
c          does not lead to an iterate, `Y', which provides a sufficient
c          reduction in the size of the natural criterion function.
c          The acceptability of the new iterate is determined by
c          examining the resultant reduction in the size of the natural 
c          criterion function. The following calls compute the trial iterate 
c          `Y_hat = Y - delta'.
c
           call dcopy((Nsub+1)*neqns,Y,1,Y_hat,1)
           call daxpy((Nsub+1)*neqns,-1.0d0,delta,1,Y_hat,1)
c                  
c          The natural criterion function is now evaluated to provide a
c          measure of how good this trial iterate is. By-products
c          of this calculation, are the quantities `PHI_hat', the evaluation
c          of the residual function at `Y_hat', and `delta = 0.5*
c          |inv(J)*PHI_hat|' for the NEXT iterate. These MAY be useful for the
c          next iteration step. If `Y_hat' turns out to be not acceptable `PHI',
c          the residual function at `Y', will be needed, so here it is not 
c          overwritten. Also, `g', the value of the natural criterion function 
c          at `Y', will be needed for comparison with `g_hat', the value of
c          the natural criterion function at `Y_hat', so `g' is not overwritten
c          either. The previous value of `delta = inv(J)*PHI' is no longer 
c          needed in either case, so it is overwritten with the new value 
c          `inv(J)*PHI_hat'.
c
           call criterion(neqns,leftbc,Nsub,mesh,Y_hat,top,blocks,bot,
     +          pivot,PHI_hat,delta,g_hat,k_discrete,work,fsub,gsub)
c
c          If the natural criterion function value is about to overflow
c          the `criterion' routine will trap that condition and set
c          g_hat = -1.0 ( a value which cannot otherwise arise).
c          Since the natural criterion function is related to the
c          size of the Newton correction, a large value signifies
c          difficulty in the Newton iteration. Set `info'=3 and exit
c          the loop.
c
           if (g_hat .LT. 0.0d0) then
               info = 3
           else
c             The Natural criterion function value has not overflowed.  
c             If the new iterate results in a sufficient reduction in
c             the size of the natural criterion function, accept it.
c
              if ( g_hat .LE. rho * g ) then
c                 The new iterate is acceptable.
c                 Update the appropriate values `Y','PHI','g' in 
c                 preparation for another fixed Jacobian step.
c                 (The flag `Fixed_Jacobian' remains TRUE.)
c
                  call dcopy((Nsub+1)*neqns,Y_hat,1,Y,1)
                  call dcopy((Nsub+1)*neqns,PHI_hat,1,PHI,1)
                  g = g_hat
              else
c                 The new iterate is not acceptable. Set the flag so that
c                 next iteration step is a damped Newton step. This also
c                 requires setting `lambda=0'.
                  Fixed_Jacobian = .FALSE.
                  lambda = 0.0d0
c
c                 If `g_hat' is greater than `rho' * `g' but less than `g' 
c                 itself, it may be reasonable to use the new iterate 
c                 and corresponding residual function value to begin the 
c                 subsequent full Newton step. These new values, `Y_hat' 
c                 and `PHI_hat', are copied to `Y' and `PHI'. 
c
                  if (g_hat .LT. g) then
                      call dcopy((Nsub+1)*neqns,Y_hat,1,Y,1)
                      call dcopy((Nsub+1)*neqns,PHI_hat,1,PHI,1)
                  else
c                     Begin the full Newton step using the values in
c                     `Y' and `PHI' from the beginning of this step.
c                     A call to `resid' is made here to obtain the `k_discrete' 
c                     values, evaluated at `Y', as these must be up-to-date
c                     with the `Y' and `PHI' values before the damped Newton
c                     step begins.
c
                      call resid(neqns,leftbc,Nsub,mesh,Y,PHI,
     +                                 k_discrete,work,fsub,gsub)
                  endif
c                 `if (g_hat .LT. g)'
              endif
c             `if (g_hat .LE. rho*g)'
           endif
c          `if (g_hat .LT. 0.0d0)'
      endif
c     `if (convrg)'
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine half_mesh(Nsub,mesh_current,mesh_new)
c
c***overview***
c
c     This routine takes as input a mesh of `Nsub' subintervals which 
c     partitions the problem interval, stored in the vector `mesh_current', 
c     The routine returns a new mesh, stored in `mesh_new', which is exactly 
c     twice as fine as the given one, i.e. has `2*Nsub' subintervals, obtained
c     by splitting each subinterval of the given mesh in half. 
c------------------------------------------------------------------------------
c***declaration of parameters***
c   imports:
      integer                Nsub
      double precision       mesh_current(0:Nsub)
c
c     `Nsub' Number of subintervals in the input mesh.
c     `mesh_current' Input meshpoints.
c
c   exports:
      double precision       mesh_new(0:2*Nsub)
c
c     `mesh_new' Meshpoints making up the new mesh.
c
c***declaration of local variables***
      integer       i
c
c     `i' loop index from 1 to `Nsub'.
c
c------------------------------------------------------------------------------
c     Called by: `mirkdc', `mesh_selector'.
c------------------------------------------------------------------------------
c
      mesh_new(0) = mesh_current(0)
c
      do 10 i = 1, Nsub
c
           mesh_new(2*i) = mesh_current(i)
c
c          Generate a new mesh point at location 2*i-1 using the average of
c          the adjacent meshpoint values.
c
           mesh_new(2*i-1) = (mesh_current(i) + mesh_current(i-1))/2.0d0
c
 10   continue
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine interp_eval(neqns,Nsub,mesh,Y,t,z,z_prime,
     +                                     k_discrete,k_interp)
c***overview***
c
c          This routine uses the current discrete solution approximation, `Y',
c     obtained on a mesh of `Nsub' subintervals, stored in `mesh', plus  
c     additional stage information associated with the underlying discrete RK 
c     method, to compute an evaluation of an interpolant to the discrete
c     solution.  This interpolant is evaluated at `t', with the value returned 
c     in `z', and the value of the first derivative returned in `z_prime'. 
c
c          The interpolant is computed as follows: Given `t', assumed to 
c     be somewhere in the problem interval, the subinterval of `mesh' that 
c     contains `t' is determined first, i.e. we find i such that 
c                `mesh(i-1) =< t < mesh(i) for t < mesh(Nsub)' or 
c                `i = Nsub for t = mesh(Nsub)'.
c     The fraction of the way `t' is through the i-th subinterval, a quantity
c     which we call `tau', is also determined during this step. Then the 
c     interpolant on this i-th subinterval is evaluated - the interpolant has
c     the form `z(mesh(i-1) + tau*hi) = yim1 + hi sum br(tau)*kr'. 
c     The weights, `br(tau)', can be computed with a call to `interp_weights'. 
c
c          The remaining information required to define the interpolant 
c     on the i-th subinterval is the current solution approximation at the 
c     (i-1)st mesh point, `yim1', which can be obtained from `Y', and the
c     stages associated with the i-th subinterval. The first `s' stages,
c     contained in `k_discrete', come from the discrete formula and
c     remaining `s_star - s' stages, contained in `k_interp', are necessary
c     for the interpolant. Each stage is a vector of length `neqns' where 
c     `neqns' is the number of differential equations. 
c
c          `k_discrete' is the vector containing the discrete stages for
c     all subintervals. The i-th set of `s*neqns' locations of `k_discrete' 
c     contains the `s' stages, each of length `neqns' corresponding to the 
c     i-th subinterval. `k_interp' plays a similar role for the `s_star-s'
c     stages associated with the interpolant.
c------------------------------------------------------------------------------
c***declaration of constants***
      integer Mxs 
      parameter(Mxs=10)
c
c     `Mxs' Maximum value for the total number of stages of the interpolant.
c
c***declaration of parameters***
c   imports:
      integer neqns, Nsub
      double precision mesh(0:Nsub), Y(neqns*(Nsub+1)), t
      double precision k_discrete(Mxs*neqns*Nsub)
      double precision k_interp(Mxs*neqns*Nsub)
c
c     `neqns' Number of differential equations.
c     `Nsub' Number of subintervals.
c     `mesh' Array of points defining the subintervals.
c     `Y' Contains the discrete solution at the mesh points;
c         the i-th set of `neqns' locations of `Y' contains 
c         the solution approximation at the i-th mesh point.
c     `t' Point at which the interpolant is to be evaluated.
c     `k_discrete' Contains the `s' discrete stage values for each subinterval.
c     `k_interp' Contains the `s_star-s' extra stage values for each
c                subinterval need to form the interpolant.
c   exports:
      double precision z(neqns), z_prime(neqns)
c
c     `z' Value of the interpolant at `t'.
c     `z_prime' Value of the first derivative of the interpolant at `t'.
c
c***declaration of local variables***
      integer i
      double precision hi
      double precision tau, weights(Mxs), weights_prime(Mxs)
c
c     `i' Index for subinterval containing `t', i ranges from 1 to `Nsub'.
c     `hi' Size of the i-th subinterval.
c     `tau' Equal to (t-mesh(i-1))/hi.
c     `weights' Contains the weight polynomials of the interpolant
c               evaluated at `tau'.
c     `weights_prime' Contains the first derivatives of the weight 
c                     polynomials of the interpolant evaluated at `tau'.
c
c***declaration of variables from /RK_s/ common block***
      integer s
      common /RK_s/ s
c     `s' Number of discrete stages.
c
c***declaration of variables from /interp_s_star/ common block***
      integer s_star
      common /interp_s_star/ s_star
c     `s_star' Total number of stages required to form the
c              interpolant on each subinterval. It includes all the
c              stages of the discrete formula plus the additional 
c              stages required for the interpolant.
c-----------------------------------------------------------------------------
c     Called by: `mirkdc'.
c     Calls to: `interval', `interp_weights', `sum_stages'.
c------------------------------------------------------------------------------
c     Call `interval' to determine the subinterval containing t.
c     The index of the subinterval is returned in `i'.
c
      call interval(Nsub,mesh,t,i)
c                                                                        
c     Compute `tau' for use in the evaluation of the weight polynomials.
c
      hi = mesh(i) - mesh(i-1)
      tau = (t - mesh(i-1))/hi
c
c     Setup the weights, evaluated at `tau'.
c
      call interp_weights(s_star, tau, weights,weights_prime)
c
c     We need to evaluate the interpolant and its derivative at `tau'.
c     This is done by taking weighted sums of the stages, using the weights
c     we have just computed, and the stages corresponding to the i-th
c     subinterval, and calling the `sum_stages' routine.
c
      call sum_stages(neqns,hi,Y((i-1)*neqns+1),s,
     +                    k_discrete((i-1)*s*neqns+1),s_star,
     +                    k_interp((i-1)*(s_star-s)*neqns+1),
     +                    weights, weights_prime, z, z_prime) 
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp_setup(neqns,tim1,hi,yim1,yi,s,ki_discrete,
     +                                     s_star,ki_interp,y,fsub)
c***overview***
c
c         The purpose of this routine is to perform some preliminary 
c     calculations to set up stages for later use in the evaluation of the 
c     interpolant. The interpolant is defined on the i-th subinterval by, 
c          `Z_i(tim1 + hi*tau) = yim1 + hi sum b_r(tau)*k_r',
c     where the stages, `k_r', depend on `yim1' and `yi', the solution 
c     approximations at the left and right endpoints, `hi', the size of the 
c     subinterval, 'tim1' the left endpoint, and the coefficients of the
c     underlying Runge-Kutta method. The weight polynomials, `b_r(tau)', can be
c     evaluated through a call to the `interp_weights' routine.
c     The first `s' stages of each subinterval were already obtained during the
c     computation of the discrete solution, so only the last `s_star - s' stages
c     need to be computed now.
c
c         This routine will construct the `s_star-s' extra stages needed by 
c     the interpolant for a single subinterval. This routine takes as input
c     the vector, `ki_discrete', containing the `s' stages associated with
c     the discrete formula on the i-th subinterval, and then uses the 
c     Runge-Kutta coefficients for the extra stages, available from the 
c     `/interp_coeffs/' common block, to construct the `s_star' new stages. 
c     These are returned in the vector `ki_interp'. Each stage is a vector of 
c     length `neqns', where `neqns' is the number of differential equations. 
c     `y' is a work array. The subroutine, `fsub', is provided by the user to 
c     define the boundary value ODE.
c-------------------------------------------------------------------------------
c***declaration of constants***
      integer  Mxs 
      parameter (Mxs=10)
c
c     `Mxs' Maximum value for the total number of stages of the interpolant.
c
c***declaration of parameters***
c   imports:
      integer            neqns
      double precision   tim1, hi, yim1(neqns), yi(neqns)
      integer            s, s_star
      double precision   ki_discrete(s*neqns)
c
c     `neqns' Number of ordinary differential equations.
c     `tim1' Left endpoint of the ith subinterval.
c     `hi' Size of the current subinterval.
c     `yim1' Solution approximation at the left endpoint of the current 
c                                                             subinterval.
c     `yi' Solution approximation at the right endpoint of the current 
c                                                            subinterval.
c     `s' Number of discrete stages.
c     `s_star' Total number of stages required to form the
c              interpolant on each subinterval. It includes all the
c              stages of the discrete formula plus the additional 
c              stages required for the interpolant.
c     `ki_discrete' Vector of `s*neqns' components. The j-th 
c                   component is of length `neqns' and corresponds
c                   to the j-th of the `s' stage values associated
c                   with the discrete formula on the i-th subinterval.
c   exports:
      double precision ki_interp((s_star-s)*neqns)
c
c     `ki_interp' Vector with same structure as `ki_discrete'.
c                 The j-th vector component corresponds to the j-th 
c                 of the new stages needed to form the interpolant.
c
c   work space:
      double precision y(neqns)
c
c     `y' Used to hold the argument for the function evaluation
c         for each new stage as it is being computed.
c
c***user-supplied subroutines***
      external fsub
c
c     `fsub' Defines f(t,y) for the first order system
c            of differential equations, y' = f(t,y).
c
c***declaration of local variables***
      integer          r, j

c     `r,j' Loop indexes from 1 to `s_star-s'.
c
c***declaration of variables from common block /interp_coeff/***
c   imports:
      double precision c_star(Mxs), v_star(Mxs), x_star(Mxs*Mxs)
      common /interp_coeff/ c_star,v_star,x_star
c
c     `c_star' First component contains the abscissa for the first new stage. 
c              The second component contains the abscissa for the second new 
c              stage, and so on.
c
c     `v_star' First component contains the blending coefficient for the first 
c              new stage. The second component contains the blending coefficient
c              for the second new stage, and so on.
c
c     `x_star' Contains the coupling coefficients for each of the new stages.
c              The Runge-Kutta coefficient matrix, X_star, is related to the 
c              data structure, `x_star' as follows (the matrix is stored 
c              columnwise in the vector):
c
c                     X_star(i,j) = x_star( (j-1)*s + i )
c
c------------------------------------------------------------------------------
c    Called by: `defect_estimate'.
c    Calls to : `fsub', `dcopy', `daxpy'.
c------------------------------------------------------------------------------
c***construction of the new stages***
c
      do 20 r = 1 , s_star-s
c
c          initialize `y' to zero
           call dcopy(neqns,0.0d0,0,y,1)
c
c          Accumulate the argument for the rth new stage in `y'.
c
c          The contributions from the discrete stages.
           do 5 j = 1 , s
c
c                The contribution from the j-th discrete stage.
                 call daxpy(neqns,x_star((j-1)*(s_star-s)+r),
     +                          ki_discrete((j-1)*neqns+1),1,y,1)
 5         continue
c
c          The contributions from the previous new stages.
c
            do 10 j = 1 , r-1
c
c                 The contribution from the j-th new stage.
                  call daxpy(neqns,x_star((j+s-1)*(s_star-s)+r),
     +                           ki_interp((j-1)*neqns+1),1,y,1)
 10         continue
c
c           Multiply the sum of the stages by `hi'.
            call dscal(neqns,hi,y,1)
c
c           Add on the contribution from `yim1'.
            call daxpy(neqns,(1-v_star(r)),yim1,1,y,1)
c
c           Add on the contribution from `yi'.
            call daxpy(neqns,v_star(r),yi,1,y,1)
c
c           All contributions to the argument for r-th new stage have
c           been made. Simply evaluate `fsub' at this argument to get the
c           new stage. The result is copied into the appropriate locations
c           of `ki_interp'.
c
            call fsub(neqns,tim1+c_star(r)*hi,y,
     +                      ki_interp((r-1)*neqns+1))
 20   continue
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp_tableau (method)
c 
c***overview***
c
c        This routine defines all the extra coefficients needed to augment
c      the discrete RK formula with an interpolant. `method' defines the MIRK 
c      method to be used. The possible values for `method' and the corresponding
c      MIRK scheme associated with each of these values are listed below.
c
c             Value of `method' |   MIRK formula
c            -------------------------------------
c                  121          |   MIRK121 scheme
c                  221          |   MIRK221 scheme
c                  232          |   MIRK232 scheme (abscissa = 0,2/3)
c                  2321         |   MIRK232 scheme (abscissa = 1,1/3)
c                  343          |   MIRK343 scheme
c                  453          |   MIRK453 scheme
c                  563          |   MIRK563 scheme
c                  5631         |   An improved MIRK563 scheme
c                         
c      The interpolant is defined on a subinterval by subinterval basis, as 
c      follows: When t can be expressed as t_i + hi*tau for some i and tau,
c      with hi = t_i+1 - t_i then the interpolant, Z(t) = Zi(t_i+ hi*tau) =
c      yi + hi sum br(tau)*kr, where the kr depend on information from the ith
c      subinterval as well as certain blending and coupling coefficients that
c      are specific to the MIRK formula and interpolant. br(tau) is the rth
c      weight polynomial of the interpolant evaluated at tau. 
c
c        The main purpose of this routine is to define the parameter which
c      gives,`s_star', the total number of stages for the interpolant and the 
c      coefficients for the `s_star-s' extra stages needed by the interpolant.
c      Once defined by this routine, the coefficients are then
c      available from the /interp_coeff/ common block, while `s_star' is
c      available from the /interp_s_star/ common block. The `order' of the
c      interpolant is available through the /interp_order/ common block.
c 
c         This routine also defines the sample point, relative to the size 
c      of a subinterval, at which the defect should be sampled. This is the 
c      point at which the interpolant must be evaluated and is thus specific 
c      to the particular interpolant being implemented.  The defect sample 
c      point, which we call `tau_star', is made available through the 
c      /sample_point/ common block.
c
c            The final portion of formula dependent information involves
c      the expressions for the weight polynomials, arising as part of the
c      definition of the interpolant for the discrete formula. This information
c      is provided through the `interp_weights' routine which evaluates all
c      the weight polynomials and their derivatives at a given point. This
c      routine can be called by other parts of the code whenever the weight 
c      polynomials and derivatives are required - i.e. whenever an evaluation 
c      of the interpolant is required.
c-----------------------------------------------------------------------------
c***declaration of constants***
      integer    Mxs
      parameter  (Mxs=10) 
c
c     `Mxs'     is the maximum number of stages of the RK method.
c
c***declaration of parameters***
c   imports:  
      integer       method
c
c     `method'   defines the method to be used
c
c***declaration of variables from common block /interp_s_star/***
      integer s_star
      common /interp_s_star/ s_star
c
c     `s_star'  is the total number of stages required to form the
c               interpolant on each subinterval. It includes all the
c               stages of the discrete formula plus the additional 
c               stages required for the interpolant.
c
c***declaration of variables from common block /interp_coeff/***
c     exports:
        double precision c_star(Mxs), v_star(Mxs), x_star(Mxs*Mxs)
        common /interp_coeff/ c_star,v_star,x_star
c                                                                        
c       The first component of `c_star' contains the abscissa
c       for the first new stage. The second component of `c_star' contains
c       the abscissa for the second new stage, and so on.
c
c       The first component of `v_star' contains the blending coefficient
c       for the first new stage. The second component of `v_star' contains
c       the blending coefficient for the second new stage, and so on.
c
c       The matrix, X_star, is related to the data structure, `x_star'
c       as follows (the matrix is stored columnwise in the vector):
c
c                      X_star(i,j) = x_star( (j-1)*s + i )
c          
c       `x_star' contains the coupling coefficients for each of the new
c       stages.
c
c***declaration of variables for the /sample_point/ common block***
c     exports:
        double precision tau_star
        common /sample_point/ tau_star
c
c       `tau_star'  is the point within each subinterval at which the 
c                       defect is to be sampled.
c
c***declaration of variables for the /rk_method/ common block***
c     exports:
        integer   amethod
        common   /rk_method/ amethod 
c
c       `amethod'   the same as `method', defines the method to be
c       used. It is defined twice here so that it may be
c       passed into a common block.
c
c***declaration of variables for the /interp_order/ common block***
c     exports:
        integer p
        common /interp_order/ p
c
c       `p' is the order of the interpolant => local error is O(h**(p+1)).
c------------------------------------------------------------------------------
c      called by: `mirkdc'
c      calls to: `dcopy'
c------------------------------------------------------------------------------
c      Translate short form of method values.
       if (method .EQ. 2) method = 221
       if (method .EQ. 4) method = 343
       if (method .EQ. 6) method = 563

c      Assign the value of `method' to `amethod' so that it may be
c      passed into a common block.
       amethod = method

c      Beginning method definition based on value of `method'.
c
       if (method. EQ. 121) then
c
c        SETUP OF THE COEFFICIENTS FOR THE NEW STAGES
c
c        Here we give coefficients defining the interpolant for a 1-stage,
c        2nd order MIRK, stage order 1, with a 2nd order interpolant.
c   
         s_star = 3
c 
         c_star(1) = 0.0d0
         c_star(2) = 1.0d0
c  
         v_star(1) = 0.0d0
         v_star(2) = 1.0d0
c
c        Column 1 of x_star            
c  
         x_star(1) = 0.0d0
         x_star(2) = 0.0d0
c
c        Column 2 of x_star
c  
         x_star(3) = 0.0d0
         x_star(4) = 0.0d0
c
c        Column 3 of x_star
c
         x_star(5) = 0.0d0
         x_star(6) = 0.0d0
c
c        Define the sample point for the interpolant of the 1-stage, 2nd 
c        order MIRK discrete formula. 
c  
         tau_star = 0.25d0
c
c        Define the value of the order of the interpolant, `p'.
c  
         p = 2
c
      elseif (method. EQ. 221) then
c
c        SETUP OF THE COEFFICIENTS FOR THE NEW STAGES
c
c        Here we give coefficients defining the interpolant for a 2-stage,
c        2nd order MIRK, stage order 1, with a 2nd order interpolant. 
c
         s_star = 2
c
c        Define the sample point for the interpolant of the 2-stage, 2nd 
c        order MIRK discrete formula. 
c          
         tau_star = 0.25d0
c
c        Define the value of the order of the interpolant, `p'.
c  
         p = 2
c         
      elseif (method. EQ. 232) then
c
c        SETUP OF THE COEFFICIENTS FOR THE NEW STAGES
c
c        Here we give coefficients defining the interpolant for a 2-stage,
c        3rd order MIRK, stage order 2, with a 3rd order interpolant.
c
         s_star = 3
c
         c_star(1) = 1.0d0
c
         v_star(1) = 1.0d0
c
c        Column 1 of x_star
c
         x_star(1) = 0.0d0
c
c        Column 2 of x_star
c
         x_star(2) = 0.0d0
c           
c        Column 3 of x_star
c
         x_star(3) = 0.0d0
c
c        Define the sample point for the interpolant of the 2-stage, 3rd 
c        order MIRK discrete formula. 
c
         tau_star  = 0.25d0
c
c        Define the value of the order of the interpolant, `p'.
c
         p = 3
c
      elseif (method. EQ. 2321) then
c  
c        SETUP OF THE COEFFICIENTS FOR THE NEW STAGES
c
c        Here we give coefficients defining the interpolant for the 
c        reflection of the above 2-stage, 3rd order MIRK, stage order 2, 
c        with a 3rd order interpolant.
c 
         s_star = 3
c
         c_star(1) = 0.0d0
c
         v_star(1) = 0.0d0
c
c        Column 1 of x_star
c
         x_star(1) = 0.0d0
c
c        Column 2 of x_star
c
         x_star(2) = 0.0d0
c            
c        Column 3 of x_star
c
         x_star(3) = 0.0d0
c
c        Define the sample point for the interpolant of the 2-stage, 3rd 
c        order MIRK discrete formula. 
c
         tau_star = 0.25d0
c
c        Define the value of the order of the interpolant, `p'.
c
         p = 3
c          
      elseif (method. EQ. 343) then
c
c        SETUP OF THE COEFFICIENTS FOR THE NEW STAGES
c
c        Here we give coefficients defining the interpolant for a 3-stage,
c        4th order MIRK, stage order 3, with a 4th order interpolant.
c
c        There is one extra stage. The corresponding row of the RK tableau
c        is    3/4 | 27/32 |  3/64  -9/64  0   0
c
         s_star = 4
c
         c_star(1) = 3.0d0/4.0d0
c
         v_star(1) = 27.0d0/32.0d0
c
c        Column 1 of x_star
c
         x_star(1) = 3.0d0/64.0d0
c
c        Column 2 of x_star
c
         x_star(2) = -9.0d0/64.0d0
c
c        Column 3 of x_star
c
         x_star(3) = 0.0d0
c
c        Column 4 of x_star
c
         x_star(4) = 0.0d0
c
c        Define the sample point for the 3-stage, 4th order MIRK 
c        discrete formula. 
c
         tau_star = 0.226d0
c
c        Define the value of the order of the interpolant, `p'.
c
         p = 4
c
      elseif (method .EQ. 453) then
c
c        SETUP OF THE COEFFICIENTS FOR THE NEW STAGES
c
c        Here we give coefficients defining the interpolant for a 4-stage,
c        5th order MIRK, stage order 3, with a 5th order interpolant.
c
         s_star = 6
c
         c_star(1) = 4.0d0/5.0d0
         c_star(2) = 13.0d0/23.0d0
c
         v_star(1) = 4.0d0/5.0d0
         v_star(2) = 13.0d0/23.0d0
c
c        Column 1 of x_star
c
         x_star(1) = 14.0d0/1125.0d0
         x_star(2) = 1.0d0/2.0d0
c
c        Column 2 of x_star
c
         x_star(3) = -74.0d0/875.0d0
         x_star(4) = 4508233.0d0/1958887.0d0
c
c        Column 3 of x_star
c
         x_star(5) = -128.0d0/3375.0d0
         x_star(6) = 48720832.0d0/2518569.0d0
c
c        Column 4 of x_star
c
         x_star(7) = 104.0d0/945.0d0
         x_star(8) = -27646420.0d0/17629983.0d0
c
c        Column 5 of x_star
c
         x_star(9)  = 0.0d0
         x_star(10) = -11517095.0d0/559682.0d0
c
c        Column 6 of x_star
c
         x_star(11) = 0.0d0
         x_star(12) = 0.0d0
c
c        Define the sample point for the interpolant of the 4-stage, 5th 
c        order MIRK discrete formula. 
c
         tau_star = 0.30d0
c
c        Define the value of the order of the interpolant, `p'.
c
         p = 5
c
      elseif (method .EQ. 4531) then
c
c        SETUP OF THE COEFFICIENTS FOR THE NEW STAGES
c
c        Here we give coefficients defining the interpolant for a 4-stage,
c        5th order MIRK, stage order 3, with a 5th order interpolant.
c
         s_star = 6
         c_star(1) = 7.0d0/13.0d0
         c_star(2) = 13.0d0/23.0d0
c
         v_star(1) = 4.0d0/5.0d0
         v_star(2) = 23.0d0/33.0d0
c
c        Column 1 of x_star
c
         x_star(1) = 36445.0d0/1542294.0d0
         x_star(2) = 1.0d0/2.0d0
c
c        Column 2 of x_star
c
         x_star(3) = -125759.0d0/1999270.0d0
         x_star(4) = 1499810527.0d0/2327157756.0d0
c
c        Column 3 of x_star
c
         x_star(5) = -3224608.0d0/11567205.0d0
         x_star(6) = -9339357584.0d0/2742721641.0d0
c
c        Column 4 of x_star
c
         x_star(7) = 915050.0d0/16194087.0d0
         x_star(8) = -11845932250.0d0/4918765257.0d0
c
c        Column 5 of x_star
c
         x_star(9)  = 0.0d0
         x_star(10) = 514365992257.0d0/113365827828.0d0
c
c        Column 6 of x_star
c
         x_star(11) = 0.0d0
         x_star(12) = 0.0d0
c
c        Define the sample point for the interpolant of the 4-stage, 5th 
c        order MIRK discrete formula. 
c
         tau_star = 0.30d0
c
c        Define the value of the order of the interpolant, `p'.
c
         p = 5
c
      elseif (method .EQ. 563) then
c
c        SETUP OF THE COEFFICIENTS FOR THE NEW STAGES 
c
c        Here we give coefficients defining the interpolant for a 5-stage,
c        6th order MIRK, stage order 3, with a 6th order interpolant.
c       
c        There are four extra stages. The corresponding rows of the RK tableau:
c        7/16 | 7/16 |  1547/32768  -1225/32768  749/4096   -287/2048 -861/16384
c        3/8  | 3/8  |  83/1536      -13/384     283/1536   -167/1536  -49/512
c        9/16 | 9/16 |  1225/32768  -1547/32768  287/2048 -749/4096    861/16384
c        1/8  | 1/8  |  233/3456     -19/1152       0         0           0
c                                            ...   -5/72      7/72      -17/216
c
         s_star = 9
c
         c_star(1) = 7.0d0/16.0d0
         c_star(2) = 3.0d0/8.0d0
         c_star(3) = 9.0d0/16.0d0
         c_star(4) = 1.0d0/8.0d0
c
         v_star(1) = 7.0d0/16.0d0
         v_star(2) = 3.0d0/8.0d0
         v_star(3) = 9.0d0/16.0d0
         v_star(4) = 1.0d0/8.0d0
c
c        Column 1 of x_star
c
         x_star(1) = 1547.0d0/32768.0d0  
         x_star(2) = 83.0d0/1536.0d0
         x_star(3) = 1225.0d0/32768.0d0
         x_star(4) = 233.0d0/3456.0d0
c        
c        Column 2 of x_star
c
         x_star(5) = -1225.0d0/32768.0d0
         x_star(6) = -13.0d0/384.0d0
         x_star(7) = -1547.0d0/32768.0d0
         x_star(8) = -19.0d0/1152.0d0
c
c        Column 3 of x_star
c
         x_star(9)  = 749.0d0/4096.0d0
         x_star(10) = 283.0d0/1536.0d0
         x_star(11) = 287.0d0/2048.0d0
         x_star(12) = 0.0d0
c
c        Column 4 of x_star
c
         x_star(13) = -287.0d0/2048.0d0
         x_star(14) = -167.0d0/1536.0d0
         x_star(15) = -749.0d0/4096.0d0
         x_star(16) = 0.0d0
c
c        Column 5 of x_star
c
         x_star(17) = -861.0d0/16384.0d0
         x_star(18) = -49.0d0/512.0d0
         x_star(19) = 861.0d0/16384.0d0
         x_star(20) = 0.0d0
c
c        Column 6 of x_star
c
         x_star(21) = 0.0d0
         x_star(22) = 0.0d0
         x_star(23) = 0.0d0
         x_star(24) = -5.0d0/72.0d0
c
c        Column 7 of x_star
c
         x_star(25) = 0.0d0
         x_star(26) = 0.0d0
         x_star(27) = 0.0d0
         x_star(28) = 7.0d0/72.0d0
c
c        Column 8 of x_star
c
         x_star(29) = 0.0d0
         x_star(30) = 0.0d0
         x_star(31) = 0.0d0
         x_star(32) = -17.0d0/216.0d0
c
c        Column 9 of x_star
c
         x_star(33) = 0.0d0
         x_star(34) = 0.0d0
         x_star(35) = 0.0d0
         x_star(36) = 0.0d0
c
c        Define the sample point for the interpolant of the 5-stage, 6th 
c        order MIRK discrete formula. 
c
         tau_star = 0.7156d0
c
c        Define the value of the order of the interpolant, `p'.
c
         p = 6
c
      elseif (method .EQ. 5631) then
c
c        SETUP OF THE COEFFICIENTS FOR THE NEW STAGES 
c
c        Here we give coefficients defining the interpolant for a 5-stage,
c        6th order MIRK, stage order 3, with a 6th order interpolant.
c
c        There are four extra stages. The corresponding rows of the RK tableau:
c        7/16 | 7/16 |  1547/32768  -1225/32768  749/4096   -287/2048 -861/16384
c        3/8  | 3/8  |  83/1536      -13/384     283/1536   -167/1536  -49/512
c        9/16 | 9/16 |  1225/32768  -1547/32768  287/2048 -749/4096    861/16384
c        1/8  | 1/8  |  233/3456     -19/1152       0         0           0
c                                            ...   -5/72      7/72      -17/216
c
         s_star = 8
c
         c_star(1) = 1.0d0/2.0d0
         c_star(2) = 1.0d0/2.0d0-sqrt(7.0d0)/14.0d0
         c_star(3) = 87.0d0/100.0d0
c
         v_star(1) = 1.0d0/2.0d0
         v_star(2) = 1.0d0/2.0d0-sqrt(7.0d0)/14.0d0
         v_star(3) = 87.0d0/100.0d0
c
c        (The elements of the matrix containing the stage coefficients are
c        stored columnwise in the vector `x_star'.
c        Column 1 of x_star
c
      x_star(1) = 1.D0/64.D0
      x_star(2) = 3.D0/112.D0+9.D0/1960.D0*dsqrt(7.D0)
      x_star(3) = 2707592511.D0/1000000000000.D0-1006699707.D0/
     #1000000000000.D0*dsqrt(7.D0)
c        
c        Column 2 of x_star
c
      x_star(4) = -1.D0/64.D0
      x_star(5) = -3.D0/112.D0+9.D0/1960.D0*dsqrt(7.D0)
      x_star(6) = -51527976591.D0/1000000000000.D0-1006699707.D0/
     #1000000000000.D0*dsqrt(7.D0)
c
c        Column 3 of x_star
c
      x_star(7) = 7.D0/192.D0*dsqrt(21.D0)
      x_star(8) = 11.D0/840.D0*dsqrt(7.D0)+3.D0/112.D0*dsqrt(7.D0)
     #*dsqrt(3.D0)
      x_star(9) = -610366393.D0/75000000000.D0+7046897949.D0/
     #1000000000000.D0*dsqrt(7.D0)+14508670449.D0/1000000000000.D0
     #*dsqrt(7.D0)*dsqrt(3.D0)
c
c        Column 4 of x_star
c
      x_star(10) = -7.D0/192.D0*dsqrt(21.D0)
      x_star(11) = 11.D0/840.D0*dsqrt(7.D0)-3.D0/112.D0*dsqrt(7.D0)
     #*dsqrt(3.D0)
      x_star(12) = -610366393.D0/75000000000.D0+7046897949.D0/
     #1000000000000.D0*dsqrt(7.D0)-14508670449.D0/1000000000000.D0
     #*dsqrt(7.D0)*dsqrt(3.D0)
c
c        Column 5 of x_star
c
      x_star(13) = 0
      x_star(14) = 88.D0/5145.D0*dsqrt(7.D0)
      x_star(15) = -12456457.D0/1171875000.D0+1006699707.D0/
     #109375000000.D0*dsqrt(7.D0)
c
c        Column 6 of x_star
c
      x_star(16) = 0
      x_star(17) = -18.D0/343.D0*dsqrt(7.D0)
      x_star(18) = 3020099121.D0/437500000000.D0*dsqrt(7.D0)+
     #47328957.D0/625000000.D0
c
c        Column 7 of x_star
c
      x_star(19) = 0
      x_star(20) = 0
      x_star(21) = -7046897949.D0/250000000000.D0*dsqrt(7.D0)
c
c        Column 8 of x_star
c
      x_star(22) = 0
      x_star(23) = 0
      x_star(24) = 0
c
c        Define the sample point for the interpolant of the 5-stage, 6th 
c        order MIRK discrete formula. 
c
         tau_star = 0.4d0
c
c        Define the value of the order of the interpolant, `p'.
c
         p = 6
c
      endif
c
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          subroutine interp_weights(s_star,tau,w,wp)
c
c***overview***
c
c          This routine evaluates the weight polynomials (and their first 
c      derivatives) corresponding to the `s_star' stages used to form the 
c      interpolant for the discrete formula defined by the `RK_tableau' 
c      routine. The evaluation takes place at `tau'. The results are returned 
c      in the `s_star' locations of the vectors `w' and `wp'. The 
c      interp_weights routine will be dependent on the parameter, `method'. 
c      `method' defines the method to be used and will be accessed via the
c      /rk_method/ common block. The possible values for `method' and the 
c      corresponding MIRK method associated with each of these values are
c      listed below.
c
c               Value of `method'   |   MIRK formula
c               -------------------------------------
c                      121          |   MIRK121 scheme
c                      221          |   MIRK221 scheme
c                      232          |   MIRK232 scheme (abscissa = 0,2/3)
c                      2321         |   MIRK232 scheme (abscissa = 1,1/3)
c                      343          |   MIRK343 scheme
c                      453          |   MIRK453 scheme
c                      563          |   MIRK563 scheme
c                      5631         |   An improved MIRK563 scheme
c                            
c          The presence of the interp_weights routine demonstrates the use 
c      of an "information hiding" principle. The more common approach would be
c      to simply provide the coefficients of the polynomial weights in
c      an array made available to the rest of the software package through
c      the /interp_coeffs/ common block. In the approach employed here we
c      "hide" the polynomial expression for each weight inside the  
c      `interp_weights' routine whose function is to provide an evaluation
c      of the weights at a given point. An advantage of this approach is 
c      that changes can be made in the method by which the weights are 
c      determined with no subsequent change needed in the rest of the 
c      software. 
c----------------------------------------------------------------------------
c***declaration of parameters***
c      imports:
         integer s_star
         double precision tau
c
c        `s_star' is the number of weight polynomials.
c        `tau'    is the point at which the weight polynomials 
c                 are to be evaluated. 
c      exports:
         double precision w(s_star), wp(s_star)
c
c        `w'    contains the values of the weight polynomials
c               of the interpolant, evaluated at `tau'.
c
c        `wp'   contains the values of the first derivative
c               of the weight polynomials of the interpolant,
c               evaluated  at `tau'.
c
c***declaration of local variables***
         double precision t2, tm1, t4m3, t2m1
c
c        Factors appearing in the expressions for the weight polynomials
c        and derivatives, defined below.
c
c***declaration of variables for the /rk_method/ common block***
c      exports:
         integer   method
         common   /rk_method/ method 
c
c        `method'  defines the method to be used
c------------------------------------------------------------------------------
c     called by: `defect_estimate', `interp_eval'
c------------------------------------------------------------------------------
c      The polynomials in `tau' that define the interpolant are very
c      dependent on the underlying discrete formula and the  interpolant
c      derived. The following part of this routine will be different for 
c      each formula implemented.
c
       if (method. EQ. 121) then
c
c         WEIGHT POLYNOMIAL IMPLEMENTATION for a 1-stage, 2nd order MIRK formula
c
          w(1) = 0.0d0
c
          w(2) = tau*(1.0d0-tau/2.0d0)
c   
          w(3) = tau**2/2.0d0
c
c         Derivative polynomials.
c
          wp(1) = 0.0d0
c  
          wp(2) = 1.0d0 - tau
c  
          wp(3) = tau
c
       elseif (method. EQ. 221) then
c
c         WEIGHT POLYNOMIAL IMPLEMENTATION for a 2-stage, 2nd order MIRK formula
c
          w(1) = tau*(1.0d0-tau/2.0d0)
c  
          w(2) = tau**2/2.0d0  
c
c         Derivative polynomials.
c
          wp(1) = 1.0d0 - tau
c  
          wp(2) = tau
c 
       elseif (method .EQ. 232) then
c
c         WEIGHT POLYNOMIAL IMPLEMENTATION for a 2-stage, 3rd order MIRK formula
c
          w(1) = tau*(2.0d0*tau**2 - 5.0d0*tau + 4.0d0)/4.0d0
c
          w(2) = -3.0d0*tau**2*(2.0d0*tau - 3.0d0)/4.0d0
c
          w(3) = tau**2*(tau - 1.0d0)
c
c         Derivative polynomials.
c
          wp(1) = (tau - 2.0d0/3.0d0)*(tau - 1.0d0)/(2.0d0/3.0d0)       
c
          wp(2) = tau*(tau - 1.0d0)/(-2.0d0/9.0d0)
c
          wp(3) = tau*(tau - 2.0d0/3.0d0)/(1.0d0/3.0d0)
c
       elseif (method .EQ. 2321) then
c
c         WEIGHT POLYNOMIAL IMPLEMENTATION for the reflection of the above
c         2-stage, 3rd order MIRK formula
c
          w(1) = 1.0d0/2.0d0*tau**3 - 1.0d0/4.0d0*tau**2
c
          w(2) = -3.0d0/2.0d0*tau**3 + 9.0d0/4.0d0*tau**2
c
          w(3) = tau**3 - 2.0d0*tau**2 + tau
c
c         Derivative polynomials.
c
          wp(1) = 3.0d0*tau*(tau-1.0d0/3.0d0)/2.0d0
c
          wp(2) = -9.0d0*tau*(tau-1.0d0)/2.0d0
c
          wp(3) = 3.0d0*(tau-1.0d0)*(tau-1.0d0/3.0d0)
c       
       elseif (method .EQ. 343) then
c 
c         WEIGHT POLYNOMIAL IMPLEMENTATION for a 3-stage, 4th order MIRK formula
c
          t2 = tau*tau
          tm1 = tau - 1.0d0
          t4m3 = tau*4.0d0 - 3.0d0
          t2m1 = tau*2.0d0 - 1.0d0
c
          w(1) = -tau*(2.0d0*tau-3.0d0)*(2.0d0*t2-3.0d0*tau+2.0d0)/6.0d0
c                                                                         
          w(2) = t2*(12.0d0*t2-20.0d0*tau + 9.0d0)/6.0d0
c
          w(3) = 2.0d0*t2*(6.0d0*t2-14.0d0*tau+9.0d0)/3.0d0
c
          w(4) = -16.0d0*t2*tm1*tm1/3.0d0
c
c         Derivative polynomials.
c
          wp(1) = -tm1*t4m3*t2m1/3.0d0
c
          wp(2) = tau*t2m1*t4m3
c
          wp(3) = 4.0d0*tau*t4m3*tm1
c
          wp(4) = -32.0d0*tau*t2m1*tm1/3.0d0
c
      elseif (method .EQ. 453) then                                
c                                                                             
      w(1) = tau*(22464.0d0-83910.0d0*tau+143041.0d0*tau**2
     + -113808.0d0*tau**3+33256.0d0*tau**4)/22464.0d0        
c   
      w(2) = tau**2*(-2418.0d0+12303.0d0*tau-19512.0d0*tau**2
     ++10904.0d0*tau**3)/3360.0d0                      
c
      w(3) = -8.0d0/81.0d0*tau**2*(-78.0d0+209.0d0*tau-204.0d0
     +*tau**2+8.0d0*tau**3)                           
c
      w(4) = -25.0d0/1134.0d0*tau**2*(-390.0d0+1045.0d0*tau
     +-1020.0d0*tau**2+328.0d0*tau**3)                   
c
      w(5) = -25.0d0/5184.0d0*tau**2*(390.0d0+255.0d0*tau-
     +1680.0d0*tau**2+2072.0d0*tau**3)                    
c
      w(6) = 279841.0d0/168480.0d0*tau**2*(-6.0d0+21.0d0*tau
     +-24.0d0*tau**2+8.0d0*tau**3)                      
c
c     Derivative polynomials.
c                                                                              
      wp(1) = 1.0d0-13985.0d0/1872.0d0*tau+143041.0d0/7488.0d0
     +*tau**2-2371.0d0/117.0d0*tau**3+20785.0d0/2808.0d0*tau**4
c                                                          
      wp(2) = -403.0d0/280.0d0*tau+12303.0d0/1120.0d0*tau**2
     +-813.0d0/35.0d0*tau**3+1363.0d0/84.0d0*tau**4               
c                                                     
      wp(3) = 416.0d0/27.0d0*tau-1672.0d0/27.0d0*tau**2+2176.0d0
     +/27.0d0*tau**3-320.0d0/81.0d0*tau**4                               
c                                          
      wp(4) = 3250.0d0/189.0d0*tau-26125.0d0/378.0d0*tau**2+
     +17000.0d0/189.0d0*tau**3-20500.0d0/567.0d0*tau**4                        
c                                        
      wp(5) = -1625.0d0/432.0d0*tau-2125.0d0/576.0d0*tau**2
     ++875.0d0/27.0d0*tau**3-32375.0d0/648.0d0*tau**4         
c
      wp(6) = -279841.0d0/14040.0d0*tau+1958887.0d0/18720.0d0*
     +tau**2-279841.0d0/1755.0d0*tau**3+279841.0d0/4212.0d0*tau**4            
c                
      elseif (method .EQ. 4531) then
c                                                                               
      w(1) = tau*(3857490.0d0-14075970.0d0*tau+24780520.0d0
     +*tau**2-20977215.0d0*tau**3+7479584.0d0*tau**4)/3857490.0d0   
c
      w(2) = tau**2*(-959595.0d0+3459675.0d0*tau-4676415.0d0
     +*tau**2+1893368.0d0*tau**3)/989100.0d0            
c
      w(3) = -32.0d0/63585.0d0*tau**2*(-13650.0d0+42100.0d0*tau
     +-47175.0d0*tau**2+10376.0d0*tau**3)           
c 
      w(4) = -250.0d0/89019.0d0*tau**2*(-2730.0d0+8420.0d0*tau
     +-9435.0d0*tau**2+4336.0d0*tau**3)               
c
      w(5) = 28561.0d0/791280.0d0*tau**2*(390.0d0+255.0d0*tau
     +-1680.0d0*tau**2+2072.0d0*tau**3)                
c
      w(6) = -279841.0d0/2449200.0d0*tau**2*(210.0d0-225.0d0*
     +tau-180.0d0*tau**2+536.0d0*tau**3)               
c                                                                            
      wp(1) = 1.0d0-938398.0d0/128583.0d0*tau+2478052.0d0/128583.0d0
     +*tau**2-399566.0d0/18369.0d0*tau**3+534256.0d0/55107.0d0*tau**4 
c                                           
      wp(2) = -9139.0d0/4710.0d0*tau+46129.0d0/4396.0d0*tau**2-
     +311761.0d0/16485.0d0*tau**3+473342.0d0/49455.0d0*tau**4 
c                                                      
      wp(3) = 58240.0d0/4239.0d0*tau-269440.0d0/4239.0d0*tau**2+
     +402560.0d0/4239.0d0*tau**3-332032.0d0/12717.0d0*tau**4 
c                                                      
      wp(4) = 65000.0d0/4239.0d0*tau-2105000.0d0/29673.0d0*tau**2
     ++3145000.0d0/29673.0d0*tau**3-5420000.0d0/89019.0d0*tau**4 
c                                                     
      wp(5) = 371293.0d0/13188.0d0*tau+485537.0d0/17584.0d0*tau**2
     +-114244.0d0/471.0d0*tau**3+1056757.0d0/2826.0d0*tau**4  
c                                                 
      wp(6) = -1958887.0d0/40820.0d0*tau+2518569.0d0/32656.0d0*
     +tau**2+839523.0d0/10205.0d0*tau**3-18749347.0d0/61230.0d0*tau**4  
c                                                                        
      elseif (method .EQ. 563) then
c 
c     WEIGHT POLYNOMIAL IMPLEMENTATION for a 5-stage, 6th order MIRK formula
c
      w(1) =
     +tau-28607.D0/7434.D0*tau**2-166210.D0/33453.D0*tau**3+334780.D0/1  
     +1151.D0*tau**4-1911296.D0/55755.D0*tau**5+406528.D0/33453.D0*tau*
     +*6
c                                                                        
      w(2) =
     +777.D0/590.D0*tau**2-2534158.D0/234171.D0*tau**3+0.208858D7/78057
     +.D0*tau**4-10479104.D0/390285.D0*tau**5+0.11328512D8/0.1170855D7*
     +tau**6
c
      w(3) =
     +-1008.D0/59.D0*tau**2+222176.D0/1593.D0*tau**3-180032.D0/531.D0
     +*tau**4+876544.D0/2655.D0*tau**5-180224.D0/1593.D0*tau**6
c                                                                        
      w(4) =
     +-1008.D0/59.D0*tau**2+222176.D0/1593.D0*tau**3-180032.D0
     +/531.D0*tau**4+876544.D0/2655.D0*tau**5-180224.D0/1593.D0*tau**6
c
      w(5) =
     +-378.D0/59.D0*tau**2+27772.D0/531.D0*tau**3-22504.D0/177.D0*tau*
     +*4+109568.D0/885.D0*tau**5-22528.D0/531.D0*tau**6
c
      w(6) =
     +-95232.D0/413.D0*tau**2+0.62384128D8/33453.D0*tau**3-49429504.D0
     +/11151.D0*tau**4+0.46759936D8/11151.D0*tau**5-46661632.D0/33453.
     +D0*tau**6
c
      w(7) =
     +896.D0/5.D0*tau**2-4352.D0/3.D0*tau**3+3456*tau**4-16384.D0/5.D0
     +*tau**5+16384.D0/15.D0*tau**6
c
      w(8) =
     +50176.D0/531.D0*tau**2-179554304.D0/234171.D0*tau**3+0.143363072D9
     +/78057.D0*tau**4-136675328.D0/78057.D0*tau**5+0.137363456D9
     +/234171.D0*tau**6
c
      w(9) =
     +16384.D0/441.D0*tau**3-16384.D0/147.D0*tau**4+16384.D0/147.D0
     +*tau**5-16384.D0/441.D0*tau**6
c
c     Derivative polynomials.
c
      wp(1) =
     +1-28607.D0/3717.D0*tau-166210.D0/11151.D0*tau**2+0.133912D7/
     +11151.D0*tau**3-1911296.D0/11151.D0*tau**4+813056.D0/11151.
     +D0*tau**5
c
      wp(2) =
     +777.D0/295.D0*tau-2534158.D0/78057.D0*tau**2+0.835432D7
     +/78057.D0*tau**3-10479104.D0/78057.D0*tau**4+0.22657024D8
     +/390285.D0*tau**5
c    
      wp(3) =
     +-2016.D0/59.D0*tau+222176.D0/531.D0*tau**2-720128.D0/531.D0
     +*tau**3+876544.D0/531.D0*tau**4-360448.D0/531.D0*tau**5
c
      wp(4) =
     +-2016.D0/59.D0*tau+222176.D0/531.D0*tau**2-720128.D0/531.D0
     +*tau**3+876544.D0/531.D0*tau**4-360448.D0/531.D0*tau**5
c
      wp(5) =
     +-756.D0/59.D0*tau+27772.D0/177.D0*tau**2-90016.D0/177.D0*
     +tau**3+109568.D0/177.D0*tau**4-45056.D0/177.D0*tau**5
c
      wp(6) =
     +-190464.D0/413.D0*tau+0.62384128D8/11151.D0*tau**2-1977
     +18016.D0/11151.D0*tau**3+0.23379968D9/11151.D0*tau**4-
     +93323264.D0/11151.D0*tau**5
c
      wp(7) =
     +1792.D0/5.D0*tau-4352*tau**2+13824*tau**3-16384*tau**4
     ++32768.D0/5.D0*tau**5 
c
      wp(8) =
     +100352.D0/531.D0*tau-179554304.D0/78057.D0*tau**2+0.57345
     +2288D9/78057.D0*tau**3-683376640.D0/78057.D0*tau**4+0.274
     +726912D9/78057.D0*tau**5
c
      wp(9) =
     +16384.D0/147.D0*tau**2-65536.D0/147.D0*tau**3+81920.D0/147
     +.D0*tau**4-32768.D0/147.D0*tau**5
c
      elseif (method .EQ. 5631) then
c 
c     WEIGHT POLYNOMIAL IMPLEMENTATION for an improved 5-stage, 6th order 
c     MIRK formula
c
      w(1) =
     #-(12233+1450*dsqrt(7.D0))*(800086000*tau**5+63579600*dsqrt
     #(7.D0)*tau**4-2936650584.D0*tau**4+4235152620.D0*tau**3
     #-201404565*dsqrt(7.D0)*tau**3+232506630*dsqrt(7.D0)*tau**2
     #-3033109390.D0*tau**2+1116511695*tau-116253315*dsqrt(7.D0)
     #*tau+22707000*dsqrt(7.D0)-191568780)*tau/2112984835740.D0
c                                                                        
      w(2) =
     #-(-10799+650*dsqrt(7.D0))*(24962000*tau**4
     #+473200*dsqrt(7.
     #D0)*tau**3-67024328*tau**3-751855*dsqrt(7.D0)*tau**2
     #+66629600*tau**2-29507
     #250*tau+236210*dsqrt(7.D0)*tau+5080365+50895*dsqrt(7.D0))
     #*tau**2/29551834260.D0
c
      w(3) =
     #7.D0/1274940.D0*(259+50*dsqrt(7.D0))
     #*(14000*tau**4-48216*tau
     #**3+1200*dsqrt(7.D0)*tau**3-3555*dsqrt(7.D0)
     #*tau**2+62790*tau**2+3610*ds
     #qrt(7.D0)*tau-37450*tau+9135-1305*dsqrt(7.D0))*tau**2
c                                                                        
      w(4) =
     #7.D0/1274940.D0*(259+50*dsqrt(7.D0))*(14000*tau**4
     #-48216*tau
     #**3+1200*dsqrt(7.D0)*tau**3-3555*dsqrt(7.D0)
     #*tau**2+62790*tau**2+3610*ds
     #qrt(7.D0)*tau-37450*tau+9135-1305*dsqrt(7.D0))*tau**2
c
      w(5) =
     #16.D0/2231145.D0*(259+50*dsqrt(7.D0))
     #*(14000*tau**4-48216*
     #tau**3+1200*dsqrt(7.D0)*tau**3-3555*dsqrt(7.D0)*tau**2
     #+62790*tau**2+3610*d
     #sqrt(7.D0)*tau-37450*tau+9135
     #-1305*dsqrt(7.D0))*tau**2
c
      w(6) =
     #4.D0/1227278493.D0*(740*dsqrt(7.D0)-6083)*(1561000*tau**2-
     #2461284*tau-109520*dsqrt(7.D0)*tau+979272+86913*
     #dsqrt(7.D0))*(tau-1)**2*tau**2
c
      w(7) =
     #-49.D0/63747.D0*dsqrt(7.D0)*(20000*tau**2
     #-20000*tau+3393)*(tau-1)**2*tau**2
c
      w(8) =
     #-1250000000.D0/889206903.D0*(28*tau**2-28*tau+9)
     #*(tau-1)**2*tau**2
c
c     Derivative polynomials.
c
      wp(1) =
     #(1450*dsqrt(7.D0)+12233)*(14*tau-7
     #+dsqrt(7.D0))*(tau-1)*(-4
     #00043*tau+75481+2083*dsqrt(7.D0))*(100*tau-87)
     #*(2*tau-1)/493029795006.D0
c
      wp(2) =
     # -(650*dsqrt(7.D0)-10799)*(14*tau-7+dsqrt(7.D0))*(37443*tau-
     #13762-2083*dsqrt(7.D0))*(100*tau-87)
     #*(2*tau-1)*tau/20686283982.D0
c    
      wp(3) =
     # 7.D0/42498.D0*(259+50*dsqrt(7.D0))*(14*tau-7+dsqrt(7.D0))
     #*(tau-1)*(100*tau-87)*(2*tau-1)*tau
c
      wp(4) =
     # 7.D0/42498.D0*(259+50*dsqrt(7.D0))*(14*tau-7+dsqrt(7.D0))
     #*(tau-1)*(100*tau-87)*(2*tau-1)*tau
c
      wp(5) =
     # 32.D0/148743.D0*(259+50*dsqrt(7.D0))*(14*tau-7+dsqrt(7.D0
     #))*(tau-1)*(100*tau-87)*(2*tau-1)*tau
c
      wp(6) =
     # 4.D0/1227278493.D0*(740*dsqrt(7.D0)-6083)*(14*tau-7+dsqrt
     #(7.D0))*(tau-1)*(100*tau-87)*(6690*tau-4085
     #-869*dsqrt(7.D0))*tau
c
      wp(7) =
     # -98.D0/21249.D0*dsqrt(7.D0)*(tau-1)*(100*tau-13)
     #*(100*tau-87)*(2*tau-1)*tau
c
      wp(8) =
     # -1250000000.D0/2074816107.D0*(14*tau-7+dsqrt(7.D0))
     #*(tau-1)*(14*tau-7-dsqrt(7.D0))*(2*tau-1)*tau
c
      endif
c
      return
      end




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interval(Nsub, mesh, t, i)
c      
c***overview***
c
c     This routine finds `i' such that `mesh(i-1) =< t < mesh(i)' for
c     `t < mesh(Nsub)'. The value of `Nsub' is returned when `t = mesh(Nsub)'.
c     `Nsub' is the number of subintervals and `mesh' is the list of points
c     defining the subintervals. 
c----------------------------------------------------------------------------
c***declaration of parameters***
c   imports:
      integer Nsub
      double precision mesh(0:Nsub), t                  
c
c     `Nsub' Number of subintervals.
c     `mesh' Array of points defining the subintervals.
c     `t' Point to be located within `mesh'.
c
c   exports:
      integer i
c
c     `i' Mesh index such that `t' satisfies `mesh(i-1) =< t < mesh(i)'  
c         for `t < mesh(Nsub)' or `i = Nsub for t = mesh(Nsub)'.
c
c***declaration of local variables***
      double precision     meshhalfnew
      integer              inew, iold, iright, ileft
c
c     `meshhalfnew' New mesh value resulting from halving the current interval.
c     `inew' Current subinterval containing `meshhalfnew'.
c     `iold' Previous subinterval containing `meshhalfnew'.
c     `ileft' Left hand mesh point bounding `meshhalfnew'.
c     `iright' Right hand mesh point bounding `meshhalfnew'.
c--------------------------------------------------------------------------- 
c      Called by: `interp_eval' or `sol_eval'.
c---------------------------------------------------------------------------
c     Handle special cases for the location of `t' first
c
      if ((t .LT. mesh(0)) .OR. (t .GT. mesh(Nsub))) then
         write(6,*)'Error-t is outside of problem interval'
      endif
c      
      if (t .EQ. mesh(Nsub)) then
                  i = Nsub
                  return
      end if
c
      if (t .EQ. mesh(0)) then
                  i = 1
                  return
      end if
c
c     At this point we know that `mesh(0) < t < mesh(Nsub)'.
c     Binary search to determine subinterval containing `t'.
c             
      inew = (Nsub/2)
      meshhalfnew = mesh(inew)
c
c     Determine if `t' is equal to `meshhalfnew'.
c
      if (t .EQ. meshhalfnew) then
         i = inew + 1
         return
      endif
c
c     `t' is not equal to the meshpoint.
c
      if (t .LT. meshhalfnew) then
c
c          `t' is on the left hand side of the interval.
c
           dowhile (t .LT. meshhalfnew)
              iold = inew
              inew = iold/2
              meshhalfnew = mesh(inew)
           enddo         
c  
           iright = iold
           ileft = inew
c
      else
c
c          `t' is on the right hand side of the interval.
c
           dowhile (t .GT. meshhalfnew)
              iold = inew
              inew = nint(real((iold + (Nsub - iold)/2)) + 0.5)
              meshhalfnew = mesh(inew)                                   
           enddo
c
           ileft  = iold
           iright = inew
c
      endif 
c     
      i = iright
c
c     Determine if `t' is equal to `meshhalfnew'.
c
      if (t .EQ. meshhalfnew) then
         i = inew + 1
         return
      endif
c       
c     `t' is not equal to the meshpoint.
c
c     `t' is bounded by `ileft' and `iright'.
c
      dowhile (abs(iright - ileft) .NE. 1)
           inew = (ileft + iright)/2
           if(mesh(inew) .LT. t) then
                 ileft = inew
           elseif (mesh(inew) .GT. t) then
                 iright = inew
           else
                 i = inew + 1
                 return
           endif
      enddo
c
      i = iright
c
c     `t' is in the `i'-th subinterval.
c 
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine jacblk(neqns,hi,ti,yi,yip1,k,LLi,RRi,
     +                            yr,dfdy,ur,dkdyi,dkdyip1,dfsub)
c
c***overview***
c 
c       ****Evaluation of the Jacobian for one subinterval****
c
c           This routine computes the left and right blocks (of the Jacobian of
c       the residual function) associated with the i-th subinterval. The left
c       and right block values are returned in the parameters `LLi' and `RRi'.
c       Each block is of size `neqns**2', where `neqns' is the number of
c       differential equations. We have,
c
c                         s                               s
c                         --                              --
c       LLi =  -I - hi * >  br * dkr/dyi,  RRi = I - hi * >  br * dkr/dyi+1
c                         --                              --
c                         r=1                             r=1
c
c      Each partial derivative, dkr/dyi, is computed as follows,
c
c                                          r-1
c                                          --
c      dkr/dyi = df/dy(tr,yr)*((1-vr)*I+hi*>  xrj * dkj/dyi)
c                                          --
c                                          j=1,
c
c      with a similar expression for each dkr/dyi+1.
c
c         The solution value, yr, at which df/dy is evaluated depends on
c      the input parameters `yi' and `yip1', the solution approximations
c      at each end of the ith subinterval, and the previously computed
c      stages provided by the input parameter `k'. The size of the ith
c      subinterval is `hi' and the endpoint is `ti'. The computation of
c      `yr' also depends on the coefficients of the RK method provided
c      through the /RK_s/,/RK_coeff/ common blocks. 
c      
c            (All matrices employed in this routine are actually stored
c      columnwise within vectors.)
c--------------------------------------------------------------------------------
c***declaration of constants***
            integer Mxs 
            parameter(Mxs=10)
c  
c            `Mxs' is the maximum value for the total number of stages
c             of the interpolant.
c
c***declaration of parameters***
c         imports:
             integer                   neqns
             double precision          hi, ti, yi(neqns),yip1(neqns)
             double precision          k(Mxs*neqns)
c
c            `neqns' is the number of differential equations.
c
c            `hi'    is the size of the current subinterval.
c
c            `ti'    is the value of the lefthand endpoint of the 
c                       ith subinterval.
c
c            `yi'    is the current solution approximation at `ti'.
c
c            `yip1'  is the current solution approximation at `ti+h', i.e.
c                    at the right hand endpoint of the subinterval.
c
c            `k'     contains the `s' stage evaluations, each of size 
c                       `neqns'associated with the `ith' subinterval.
c         exports:
              double precision          LLi(neqns**2),RRi(neqns**2)
c
c            `LLi' is the partial derivative of the ith component of
c                  the residual function with respect to yi
c
c            `RRi' is the partial derivative of the ith component of the
c                  residual function with respect to yi+1
c
c***declaration of work arrays***
            double precision        yr(neqns)
            double precision        dfdy(neqns**2)
            double precision        ur(neqns**2)
            double precision        dkdyi(neqns**2*Mxs)
            double precision        dkdyip1(neqns**2*Mxs)

c            `yr'      the argument for the evaluation of the `rth' stage.
c
c            `dfdy'    during the `rth' time through the main loop below, it
c                      contains the derivative of the right hand side of the 
c                      system of ODE's, evaluated at the arguments of the 
c                      'rth' stage of the RK scheme.
c
c            `ur'      used to store the intermediate values arising in the
c                      computation of each partial derivative.
c
c            `dkdyi'   The `rth' set of `neqns**2' locations of this vector
c                      will contain the value of the rth partial derivative
c                      with respect to `yi'.
c            
c            `dkdyip1' The `rth' set of `neqns**2' locations of this vector
c                      will contain the value of the rth partial derivative
c                      with respect to `yip1'.            
c
c***user - supplied subroutines***
            external dfsub
c
c           dfsub - This routine defines the Jacobian of the system of
c                   differential equations with respect to y, i.e.
c                   df/dy.
c
c***declaration of local variables***
            integer                 l,i,m,r, j,i1,i2
            double precision        sum
            integer                 neqnssq 
            integer                 alpha,beta,sigma          
            double precision        tr
            logical                 dkdyi_zero(Mxs)
            logical                 dkdyip1_zero(Mxs)
            logical                 matrix_mult_dkdyi
            logical                 matrix_mult_dkdyip1
c
c            `l,i,m,r,j'     are loop index variables
c
c            `sum'     is an accumulator
c
c            `i1',`i2' are indexes used in calculating ur and yr
c
c            `neqnssq' is the value of neqns squared
c
c            `alpha',`beta',`sigma' are indexes used in the
c                     multiplication of the two matrices `ur' and `dfdy'
c
c            `tr'     the abscissa for the evaluation of the `rth' stage.
c
c            `dkdyi_zero' is a logical array. If the computation of the rth
c                         stage has been completed then a value of `.TRUE.' in
c                         `dkdyi_zero(r)' indicates that every component of
c                         the rth partial derivative with respect to `yi'
c                         is zero. If the computation of the rth stage has
c                         not been completed then a value of `.TRUE.'
c                         indicates that the `ur' components of the rth
c                         partial derivative are currently zero. A value of
c                         `.FALSE.' in `dkdyi_zero(r)' implies the rth
c                         partial derivative is non-zero in at least one of
c                         its components.
c
c            `dkdyip1_zero' is a logical array. If the computation of the rth
c                           stage has been completed then a value of `.TRUE.'
c                           in `dkdyip1_zero(r)' indicates that every component
c                           of the rth partial derivative with respect to
c                           `yip1' is zero. If the computation of the rth
c                           stage has not been completed then a value of
c                           `.TRUE.' indicates that the `ur' components of
c                           the rth partial derivative are currently zero.
c                           A value of `.FALSE.' in `dkdyip1_zero(r)'
c                           implies the rth partial derivative is non-zero
c                           in at least one of its components. `dkdyi_zero'
c                           and `dkdyip1_zero' are both used to avoid
c                           multiplications by zero.
c
c            `matrix_mult_dkdyi' is a logical variable. For the compuation
c                                of a given stage, a value of `.FALSE.'
c                                indicates that  the partial derivative
c                                with respect to `yi'
c                                has the form `(1-v(r))*dfdy' and that
c                                multiplying the matrix `ur' by the matrix
c                                `dfdy' is unnecessary. A value of `.TRUE.'
c                                indicates that the partial derivative
c                                does not have the form `(1-v(r))*dfdy' and
c                                that multiplying the matrix `ur' by the
c                                matrix `dfdy' may be necessary.
c
c            `matrix_mult_dkdyip1' is a logical variable. For the compuation
c                                  of a given stage, a value of `.FALSE.'
c                                  indicates that  the partial derivative
c                                  with respect to `yip1'
c                                  has the form `v(r)*dfdy' and that
c                                  multiplying the  matrix `ur' by the matrix
c                                  `dfdy' is unnecessary. A value of `.TRUE.'
c                                  indicates that the partial derivative
c                                  does not have the form `v(r)*dfdy' and
c                                  that multiplying the matrix `ur' by the
c                                  matrix `dfdy' may be necessary.
c
c***declaration of variables in common block /RK_s/***
            integer                 s
            common /RK_s/ s
c
c           `s' is the number of discrete stages of the Runge-Kutta formula.
c
c***declaration of variables in common block /RK_coeff/***
            double precision        v(Mxs), b(Mxs)
            double precision        c(Mxs), x(Mxs**2)
            common /RK_coeff/ v, b, c, x
c         
c           `v', `b', `c',and `x' are coefficients that define 
c           the discrete Runge-Kutta formula.  
c
c           The first component of `c' contains the abscissa for the first 
c           new stage. The second component of `c' contains the abscissa for 
c           the second new stage, and so on.
c
c           The first component of `v' contains the blending coefficient
c           for the first new stage. The second component of `v' contains
c           the blending coefficient for the second new stage, and so on.
c
c           The matrix, X, is related to  the data structure, `x'
c           as follows (the matrix is stored columnwise in the vector):
c
c                 X(i,j) = x( (j-1)*s + i )
c          
c           `x' contains the coupling coefficients for each of the new
c           stages
c
c------------------------------------------------------------------------------
c            called by: `newmat'
c            calls to: `dfsub'
c------------------------------------------------------------------------------
c***COMPUTATION OF THE PARTIAL DERIVATIVES***
c
c     We compute the partial derivatives in pairs: dkr/dyi and dkr/dyi+1 
c     during the `rth' pass of the following loop. The results are stored 
c     in the appropriate locations of the vectors `dkdyi' and `dkdyip1', 
c     respectively.
c
c     Initialize the components `dkdyi_zero' and `dkdyi_zero' to `.TRUE.'
c     indicating that all of the partial derivatives are currently zero.
c
      do 1 j = 1,s
        dkdyi_zero(j)   = .TRUE.
        dkdyip1_zero(j) = .TRUE.
 1    continue
c
      neqnssq = neqns**2
c
      do 40 r = 1, s
c
c         Initialize `matrix_mult_dkdyi' to `.TRUE.' indicating that the
c         partial derivative is not of the form `(1-v(r))*dfdy' and that
c         multiplying the matrix `ur' by `dfdy' may still be necessary.
c
          matrix_mult_dkdyi = .TRUE.
c
c         Construct the arguments for the `rth' stage, needed here for the 
c         evaluation of `df/dy' 
c
c         Compute the abscissa
c
          tr = ti + c(r)*hi
c
c         Begin the construction of `yr' by setting it to zero.
c
          do 2 j=1, neqns
            yr(j) = 0.0d0
 2        continue
c
c          Add on weighted contributions from first `r-1' stages.
c          The `rth' stage of the current subinterval is stored in
c          `k((r-1)*neqns+1)'.
c
           do 5 j = 1, r-1
c
            i1 = (j-1)*s + r
            i2 = (j-1)*neqns
c
            if ( x(i1) .NE. 0.0d0 ) then
c
c             If `x(i1) = 0' then there is no need to perform the following
c             multiplication.

              do 3 l = 1, neqns
                yr(l) = yr(l) + x(i1)*k(i2+l)
 3            continue
c
             endif
c
 5        continue
c 
c          Multiply yr by hi.
c          
           do 7 l =1, neqns
              yr(l) = hi*yr(l)
 7         continue
c
c          Add on weighted contributions from `yi' and `yip1'.
c
           do 10 l=1,neqns
              yr(l) = (1.0d0 - v(r))*yi(l) + v(r)*yip1(l) + yr(l)
 10        continue
c
c          The arguments for the `rth' stage, `tr' and `yr' have been 
c          constructed. Evaluate df/dy, the derivative of the system of ODE's 
c          by first setting it equal to zero then calling the user defined 
c          subroutine `dfsub' at these arguments to fill in the non-zero
c          entries.
c
            do 12 i = 1,neqnssq
               dfdy(i) = 0.0d0
 12         continue
c          
            call dfsub(neqns,tr,yr,dfdy)
c
c          The values of the `rth' partial derivatives with respect to `yi'
c          and `yip1' are obtained by multiplying `dfdy' by a weighted linear
c          combination of the partial derivatives of the previous `r-1'
c          partial derivatives. We will store these linear combinations, as we
c          construct them, in the matrix, `ur'.
c          
c          Construct the `ur' matrix for dkr/dyi
c
           do 14 j = 1, neqnssq
                 ur(j) = 0.0d0
 14        continue
c
c          Add on weighted multiples of the previous partial derivatives.
c          The `jth' previous partial derivative is stored in
c          `dkdyi((j-1)*neqns**2+1)'. This can be done through a call to
c          `daxpy' since we are storing all matricies in vector form.
c
           do 15 j = 1, r-1
c     
            i1 = (j-1)*s + r
            i2 = (j-1)*neqnssq
c
            if ( x(i1) .NE. 0  .AND. ( .NOT. dkdyi_zero(j) ) ) then
c
c             If `x(i1) = 0' or the jth partial derivative equals zero then
c             there is no need to perform the following multiplication.
c
              dkdyi_zero(r) = .FALSE.

c             Indicates the matrix `ur' is non zero in at least one of the
c             components and we will require a matrix multiply later.
c             It also indicates that the rth partial derivative of
c             `dkdyi' will be non zero in one of the components.
c
              do 16 l= 1, neqnssq
                 ur(l) = ur(l) + x(i1)*dkdyi(i2+l)
 16           continue
c
            endif
c
 15        continue
c
c          Multiply ur by hi
c
           if ( .NOT. dkdyi_zero(r) ) then
c
c            If the `ur' matrix for the rth partial derivative equals zero
c            then there is no need to perform the following multiplication.

             do 18 j = 1, neqnssq
               ur(j) = hi*ur(j)
 18          continue
c
           endif
c
c          Add on weighted contribution (1.0-v(r))*I.
c
           alpha = (r-1)*neqns**2

           if ( v(r) .NE. 1.0d0 ) then
c
c            If `1-v(r) = 0' then there is no need to perform the following
c            multiplication.
c
             if ( dkdyi_zero(r) ) then

               matrix_mult_dkdyi = .FALSE.

c              Signifies that multiplying the matrix `ur' by the
c              matrix `dfdy' is not necessary later on
c              because the matrix dkdyi can be obtained simply
c              by multiplying the scalar (1-v(r)) times the matrix dfdy.
c              We proceed to perform the scalar multiplication here.
c
               do j = 1,neqnssq
                 dkdyi(alpha + j) = ( 1.0d0-v(r) )*dfdy(j)
               enddo
c
             else
    
               do 20 j =1, neqnssq, neqns+1
                 ur(j) = 1.0d0-v(r) + ur(j)
 20            continue
c
             endif
c
             dkdyi_zero(r) = .FALSE.
c
c            Indicates the matrix `ur' is non zero in at least one of the
c            components and we will require a matrix multiply later.
c            It also indicates that the rth partial derivative of
c            `dkdyi' will be non zero in one of the components.
c
           endif
c
c           `ur' now contains the proper linear combination of previously
c           computed partial derivatives with respect to `yi'. 
c           The value of dkr/dyi is obtained simply by multiplying `ur' and
c           `dfdy' together. The resultant matrix value is stored in
c           the `neqns*neqns' locations of `dkdyi((r-1)*neqns**2+1)'.
c
c           Multiply `ur' and `dfdy'
c
           if ( matrix_mult_dkdyi ) then
c
            if ( .NOT. dkdyi_zero(r) ) then
c
c            If the `ur' matrix for the rth partial derivative equals zero
c            or rth partial derivative could be obtained by multiplying
c            the scalar `(1-v(r)' by `dfdy' then there is no need to perform
c            the following matrix multiplication.
c
             do 26 i =1, neqns
                 beta = alpha + i
c                       
                 do 27 j = 1, neqns
                   sigma = (j-1)*neqns
                   sum = 0.0d0
c
                   do 28 m = 1,neqns
                      sum = sum + dfdy((m-1)*neqns+i)*ur(sigma+m)        
  28               continue
c
                   dkdyi(sigma+beta) = sum
  27             continue
c
  26         continue
c
            else

c             If the `ur' matrix is equal to zero then the corresponding rth
c             partial derivative is set to zero. As well, the value of
c             `dkdyi_zero(r)' will left as `.TRUE.' to avoid
c             multiplying by zero in subsequent stage computations.
c
              do 25 j = 1,neqnssq
               dkdyi(alpha + j) = 0.0d0
 25           continue
c
            endif
c                 `(.NOT. dkdyi_zero(r) )'
c
           endif
c               `( matrix_mult_dkdyi) '
c
c          Initialize `matrix_mult_dkdyip1' to `.TRUE.' indicating that the
c          partial derivative is currently not of the form `v(r)*dfdy' and
c          that multiplying the matrix `ur' by `dfdy' may still be necessary.
c
            matrix_mult_dkdyip1 = .TRUE.
c                 
c           Construct the `ur' matrix for dkr/dyip1
c
            do 24 j = 1, neqnssq
                  ur(j) = 0.0d0
 24         continue
c 
c           Add on weighted multiples of the previous partial derivatives.
c           The `jth' previous partial derivative is stored in
c           `dkdyip1((j-1)*neqns**2+1)'.
c
            do 30 j = 1, r-1
c            
            i1 = (j-1)*s+ r
            i2 = (j-1)*neqnssq
c    
            if ( x(i1) .NE. 0  .AND.( .NOT. dkdyip1_zero(j)) ) then

              dkdyip1_zero(r) = .FALSE.

c             Indicates the matrix `ur' is non zero in at least one of the
c             components and we will require a matrix multiply later.
c             It also indicates that the rth partial derivative of
c             `dkdyip1' will be non zero in one of the components.
c
              do 31 l = 1, neqnssq
                ur(l) = ur(l) + x(i1)*dkdyip1(i2+l)
 31           continue
c
             endif
c
 30         continue
c
c           Multiply ur by hi
c
            if ( .NOT. dkdyip1_zero(r) ) then
c
c            If the `ur' matrix for the rth partial derivative equals zero
c            then there is no need to perform the following multiplication.

             do 35 j = 1, neqnssq
               ur(j) = hi*ur(j)
 35          continue
c
            endif
c
c           Add on weighted contribution v(r)*I
c
            if ( v(r) .NE. 0.0d0 ) then
c
c            If `v(r) = 0' then there is no need to perform the following
c            multiplication.
c
               if ( dkdyip1_zero(r) ) then
c
                matrix_mult_dkdyip1 = .FALSE.

c               Signifies that multiplying the matrix `ur' by the
c               matrix `dfdy' is not necessary later on
c               because the matrix dkdyi can be obtained simply
c               by multiplying the scalar v(r) times the matrix dfdy.
c               We proceed to perform the scalar multiplication here.
c
                do j = 1,neqnssq
                  dkdyip1(alpha + j) = v(r) * dfdy(j)
                enddo
c
               else

                do 36 j = 1, neqnssq, neqns+1
                  ur(j) = v(r) + ur(j)
 36             continue
c
               endif
c
               dkdyip1_zero(r) = .FALSE.
c
c              Indicates the matrix `ur' is non zero in at least one of the
c              components and we will require a matrix multiply later.
c              It also indicates that the rth partial derivative of
c              `dkdyip1' will be non zero in one of the components.
c
            endif
c
c           `ur' now contains the proper linear combination of previously
c           computed partial derivatives with respect to `yip1'. 
c           The value of dkr/dyi+1 is obtained simply by multiplying `ur' 
c           and `dfdy' together. The resultant matrix value is stored in
c           the `neqns*neqns' locations of `dkdyip1((r-1)*neqns**2+1)'.
c        
c          Multiply `ur' by `dfdy'
c   
          if ( matrix_mult_dkdyip1 ) then

           if ( .NOT. dkdyip1_zero(r) ) then

c            If the `ur' matrix for the rth partial derivative equals zero
c            or rth partial derivative could be obtained by multiplying
c            the scalar `v(r)' by `dfdy' then there is no need to perform
c            the following matrix multiplication.
c
             do 39 i = 1, neqns
                   beta = alpha + i
c
                   do 38 j = 1 ,neqns
                    sum=0.0d0
                    sigma = (j-1)*neqns
c
                    do 37 m = 1, neqns
                     sum = sum + dfdy((m-1)*neqns+i)*ur(sigma+m)        
 37                 continue
c
                    dkdyip1(sigma+beta) = sum
 38                continue
c
 39          continue
c               
            else
c
c            If the `ur' matrix is equal to zero then the corresponding rth
c            partial derivative is set to zero. As well, the value of
c            `dkdyip1_zero(r)' will left as `.TRUE.' to avoid
c            multiplying by zero in subsequent stage computations.
c
             do 41 j = 1,neqnssq
              dkdyip1(alpha + j) = 0.0d0
 41          continue
c
            endif
c               `( .NOT. dkdyip1_zero(r) )'
c
          endif
c               `( matrix_mult_dkdyip1 )'
c
 40     continue
c
c      All partial derivatives have now been computed.
c
c      COMPUTATION OF THE BLOCKS, LLi AND RRi.
c
c      The left and right blocks are now constructed from linear combinations
c      of the partial derivatives.
c
c      Initialize LLi and RRi to 0
c
            do 44 j = 1, neqnssq
                  LLi(j) = 0.0d0
                  RRi(j) = 0.0d0
 44         continue
c
c      Add appropriately weighted partial derivative values to LLi and RRi.
c
       do 50 r = 1, s
c
        i1 = (r-1)*neqnssq
c
        do 51 l = 1, neqnssq
         LLi(l) = LLi(l) - b(r)*dkdyi(i1+l)
         RRi(l) = RRi(l) - b(r)*dkdyip1(i1+l)
 51     continue
c
 50    continue
c
c      Multiply LLi and RRi by hi
c
       do 53 j = 1, neqnssq
         LLi(j) = hi*LLi(j)
         RRi(j) = hi*RRi(j)
 53    continue
c      
c      Add on -I and I, respectively
c
       do 55 j= 1, neqnssq, neqns+1
         LLi(j) = -1.0d0 + LLi(j)
         RRi(j) =  1.0d0 + RRi(j)
 55    continue
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mesh_selector(neqns,Nsub,mesh_current,defect,tol,
     *                Nsub_star,mesh_new,info,MxNsub,s_hat)
c
c***overview***
c
c     The purpose of this routine is to examine the relative defect 
c     associated with a continuous solution approximation in order to decide
c     upon a good mesh. The `Nsub+1' points describing the current mesh are 
c     passed in through the vector `mesh_current'. An estimate of the blended 
c     relative/absolute defect, is passed in through the vector `defect', 
c     which has `Nsub' components, each of size `neqns', where `neqns' is
c     the number of differential equations. The user defined tolerance, which
c     will be used to assess the relative defect, is passed in through the 
c     parameter, `tol'.
c 
c       After examining the blended relative/absolute defect, this routine 
c     will call either the `half_mesh' routine, to generate a new mesh by
c     halving each subinterval of the current mesh, or the routine 
c     `redistribute'. `redistribute' generates a new mesh with a different 
c     number of points than the current mesh and with a different distribution
c     of points, depending how the relative defect is distributed over the
c     current mesh.
c
c       The redistribution process involves monitoring the mesh function
c     values which are based on the defect. A mesh redistribution is performed 
c     if the maximum mesh function value, is sufficiently greater than
c     the average value. If so, a new mesh is determined that approximately
c     equidistributes the mesh function over the subintervals of the new mesh. 
c     The number of points in the new mesh is determined from accuracy
c     considerations.
c
c          The algorithm implemented here is based on [Ascher, Mattheij, and 
c     Russell, Chapter 9]. The number of subintervals in the new mesh is 
c     returned in `Nsub_star', while the points of the new mesh are returned 
c     in the vector `mesh_new'. This routine will fail if the prescribed number
c     of new points is greater than the maximum number allowed, `MxNsub', in 
c     which case we set `info = -1' to signal this.
c
c     `s_hat' is used here for temporary storage.
c-------------------------------------------------------------------------------
c***declaration of constants***
         integer Mxs
         parameter (Mxs=10)

c        `Mxs'    is the maximum value for the total number of stages
c                 of the interpolant.
c
c***declaration of factors controlling mesh redistribution***
         double precision     safety_factor
         double precision     rho
         double precision     upper_new_mesh
         double precision     lower_new_mesh
c
         parameter(safety_factor = 1.3d0)
         parameter(rho = 1.0d0)
         parameter(upper_new_mesh = 4.0d0)
         parameter(lower_new_mesh = 0.5d0)
c                                                                    
c        `safety_factor' a factor to increase the prediction of the appropriate 
c                        number of subintervals needed to achieve the given
c                        tolerance in the new mesh.
c
c        `rho'   a factor which determines how much bigger the maximum mesh
c                function value must be, compared to the average value, in 
c                order for a mesh redistribution to occur.
c 
c                Setting `rho=1' implies that a mesh redistribution will take
c                place everytime.
c
c        `upper_new_mesh', `lower_new_mesh'
c
c                are upper and lower bounds on the ratio between the number
c                of subintervals in current mesh and the number of subintervals
c                in the new mesh.
c
c***declaration of parameters***
c      imports:
         integer              neqns
         integer              Nsub
         double precision     mesh_current(0:Nsub)
         double precision     defect(neqns*Nsub)
         double precision     tol
         integer              MxNsub
c
c        `neqns'        the number of differential equations.
c
c        `Nsub'         the number of subintervals in the current mesh.
c 
c        `mesh_current' a vector of points defining the current mesh.
c
c        `defect'       a vector giving an estimate of the blended 
c                       relative/absolute defect, with `Nsub' components,
c                       each of size `neqns'.
c 
c        `tol'          the user defined tolerance. 
c
c        `MxNsub'       is the maximum number of subintervals.
c      exports:
         integer Nsub_star
         double precision mesh_new(0:MxNsub)
         integer info
c
c        `Nsub_star' the number of subintervals of the new mesh.
c
c        `mesh_new' the points defining the new mesh.
c
c        `info' communication flag
c
c                       info = -1    size of the new mesh would be too
c                                      large
c
c                       info =  0    successful generation of a new mesh
c
c***declaration of work array***
         double precision    s_hat(Nsub)
c
c        `s_hat' a vector of mesh function values, with one component per
c                subinterval. The mesh function value is based on the
c                value of the norm of the defect, the tolerance `tol',
c                the size of the subinterval, and the order of the Runge-Kutta
c                scheme, `p'.

c***declaration of local variables***
         integer             i
         double precision    hi
         double precision    r1, r2, r3 
         integer             offset
         integer             Nsub_pred
         double precision    norm
c
c        `i'     loop index from 1 to `Nsub'.
c
c        `hi'    size of ith subinterval of current mesh.
c
c                                    1/(p+1)
c                s_hat(i) = (|| defect(i) ||/tol )
c                           ----------------------
c                                    hi
c
c        `r1'    maximum size of weighted mesh function values for any 
c                subinterval of current mesh.
c
c        `r2'    weighted sum of mesh function values, `s_hat(i)', with
c                weights, `hi'. `r2' is the predicted number of points 
c                that should be in the new mesh.
c
c        `r3'    average of weighted mesh function values.
c
c        `offset'     offset used in the vector `defect'  
c                     for the `neqns' solution components of the ith 
c                     subinterval.
c
c        `Nsub_pred'  predicted value for number of subintervals in
c                     new mesh if a redistribution is needed.
c
c        `norm'       max norm of the defect estimate for the ith subinterval.
c
c***declaration of variables in common block /interp_order/***
c      imports:
         integer p
         common /interp_order/ p
c
c        `p'      the order of the interpolant associated with the underlying
c                 discrete RK formula.
c
c***declaration of variables in common block /IOCNTRL/***
c      imports:
         integer   print_level_0, print_level_1, print_level_2
         parameter (print_level_0 = 0, print_level_1 = 1,
     *              print_level_2 = 2)
         integer                 profile
         common /IOCNTRL/  profile
c
c        `print_level_0' No output.
c        `print_level_1' Intermediate output.
c        `print_level_2' Full output.        
c        `profile' Controls output, to standard output, of profiling infor-
c                  mation such as Newton iteration counts, mesh selection,
c                  relative defect estimate. This variable is identical to
c                  the `output_control' parameter.
c-----------------------------------------------------------------------------
c        called by: `mirkdc'
c        calls to : `half_mesh', `redistribute', `idamax'
                    integer idamax
c------------------------------------------------------------------------------
c     Compute weighted mesh function values as well as r1 and r2.
c     `r1' will be the maximum value of `s_hat(i)*hi'; `r2' will be the sum
c     of these same values. Initialize each to zero. Also set `info' to 0.
c
      info = 0
      r1 = 0.0d0
      r2 = 0.0d0
c
      do 10 i = 1 , Nsub
c
           hi = mesh_current(i)-mesh_current(i-1)
c
c          Compute offset to beginning of defect estimate vector for
c          ith subinterval. This vector is of length `neqns'.
c
           offset = (i-1)*neqns
c
c          Compute norm of the relative defect estimate vector for the 
c          ith subinterval.
c
           norm = dabs(defect(offset+idamax(neqns,defect(offset+1),1)))
c
c          Compute corresponding mesh function value, `s_hat(i)'.
c
           s_hat(i) = (norm/tol)**(1.0d0/(p+1))/hi
c
c          Update maximum weighted mesh function value, if necessary.
c
           if ((s_hat(i)*hi) .GT. r1) r1 = s_hat(i)*hi
c
c          Accumulate ith weighted mesh function value in `r2'
c
           r2 = r2 + s_hat(i)*hi
c
 10   continue
c
c     `r3' is the average of the weighted mesh function values
c
      r3 = r2/Nsub
c
c     `r2' gives an estimate of the number of subintervals that the
c     new mesh should have. We employ a safety factor to increase
c     the number of points slightly to increase the probability
c     of achieving the tolerance.
c  
      Nsub_pred = (safety_factor*r2)+1
c
c     To avoid cycling, make sure new mesh has a few more points than
c     the previous one.
c
      if ( dabs( (Nsub_pred - Nsub)/(1.0d0*Nsub)) .lt. 0.1d0)
     +                             Nsub_pred = 1.1d0 * Nsub
c
c     Decide whether to half the current mesh or produce a new mesh by 
c     redistributing the points of the current mesh, while adding or deleting
c     a few. Redistribute if the maximum weighted mesh function value
c     is more than `rho' times the size of the average value.
c
      if (r1 .LE. rho*r3) then
c
c          The maximum value is not sufficiently large. That is, the
c          distribution of mesh function values is sufficiently uniform
c          and mesh redistribution is unwarrented. In this case, the current
c          mesh is halved.
c
           Nsub_star = 2*Nsub
c
           if (Nsub_star .GT. MxNsub) then
c
c              mesh halving fails; signal this by setting info = -1
               info = -1
c
               if (profile .GT. print_level_0) then
                 write(6,*)' New mesh would be too large'
               endif
c
           else 
c
               if (profile .GT. print_level_0) then
                 write(6,*)' Half the current mesh'
               end if
c
               call half_mesh(Nsub,mesh_current,mesh_new)
c              A new mesh twice as fine as the current one is returned
c              in `mesh_new'.
c
           endif
c
      else
c
c          The distribution of mesh function values is sufficiently
c          non-uniform to warrent a mesh redistribution.
c
c          The value of `Nsub_pred' gives the required number of
c          subintervals for the new mesh. However, if this is 
c          too big or too small we will not use it. 
c
           Nsub_star = Nsub_pred
c
           if ( Nsub_star .GT. upper_new_mesh * Nsub ) 
     +               Nsub_star = upper_new_mesh * Nsub
           if ( Nsub_star .LT. lower_new_mesh * Nsub ) 
     +               Nsub_star = lower_new_mesh * Nsub
c
c          Also, check that `Nsub_star' is less than the maximum value.
c
           if (Nsub_star .GT. MxNsub) then
c
c              mesh redistribution fails; signal this by setting info = -1
               info = -1
c
               if (profile .GT. print_level_0) then
                  write(6,*)' New mesh would be too large'
               endif
c
           else 
c
c               if (profile .GT. print_level_0) then
c                  write(6,*)'Call to redistibute with Nsub_star=',
c    +                                                Nsub_star          
c               end if                                                         
c                                                                        
               call redistribute(Nsub,mesh_current,s_hat,
     +                                       Nsub_star,mesh_new)          
c              The new mesh, having `Nsub_star' subintervals, is returned in
c              `mesh_new'
c
           endif
c
      endif
c          `r1. LE. rho*r3'
c
      return
      end









cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine NewIter(neqns,leftbc,Nsub,mesh,Y,newtol,maxiter,
     +             info,PHI,delta,top,bot,blocks,pivot,k_discrete,
     +             work,fsub,gsub,dfsub,dgsub)
c
c***overview***
c
c            This routine takes a mesh of `Nsub' subintervals, partitioning
c     the problem interval and a discrete initial approximation provided
c     at the meshpoints, (contained in the vectors `Y' and `mesh',
c     respectively) and uses a modified Newton method to produce a 
c     solution approximation to within a tolerance of `newtol', for
c     the associated discrete nonlinear system `PHI(Y) = 0', where `PHI(Y)' is
c     the residual function. The number of differential equations 
c     is `neqns', while the number of boundary conditions at the lefthand 
c     endpoint of the problem interval is `leftbc'. The iteration completes
c     successfully if it produces a converged solution within `maxiter' 
c     iterations. A successful return is indicated by `info'= 0. The 
c     updated solution is returned in `Y'. The arrays `PHI', `delta', `top', 
c     `bot', `blocks', `pivot', and `work' provide working space for the 
c     `resid', `damped_newt', and `fixed_jacob' routines. This routine 
c     provides, as a by-product, the stage values associated with the
c     application of a Runge-Kutta method on each subinterval of the mesh.
c     The array `k_discrete' returns these stages. The subroutines,
c     `fsub', `gsub', `dfsub', and `dgsub', are provided by the user to
c     define the boundary value ODE and the boundary conditions and their
c     derivatives.
c
c            The modified Newton iteration implemented here proceeds
c     through a sequence of steps, each of which must be either a damped 
c     Newton step or a fixed Jacobian step. The code decides, depending 
c     on how well the iterates appear to be converging, which type of step
c     to take each time. This routine is organized as follows. After some 
c     initialization, we begin a loop which executes the modified Newton 
c     iteration. The kind of step to be taken is controlled by the value of 
c     the flag, `Fixed_Jacobian', which is initially set to FALSE since the 
c     first step must be a full Newton step.  
c
c            During a damped Newton step an iterative process is used to 
c     choose a damping factor, `lambda', which is used to provide a damped 
c     update to the current iterate `Y = Y - lambda*delta' where `delta' is 
c     the Newton correction; `delta = inv(J)*PHI' where `PHI' is the residual
c     function evaluated at `Y' and J is the corresponding Jacobian also 
c     evaluated at `Y'. The use of a damping factor can improve the 
c     convergence properties of the Newton iteration as well as allowing for 
c     the early detection of a divergent iteration. The update is accepted
c     if it produces a sufficient reduction in the value of the natural
c     criterion function (which is related to the size of the residual
c     function). At the end of a successful damped Newton step, the routine
c     takes a fixed Jacobian step if the step was a full Newton step 
c     (i.e. `lambda' = 1). 
c 
c          The fixed Jacobian step consists of first computing a potential
c     new iterate, `Y_hat = Y - delta' where `delta' is the Newton correction
c     for this fixed Jacobian step. Here `delta = inv(J)*PHI', where `PHI' is 
c     the residual function evaluated at `Y' and where J, the Jacobian, has 
c     not been evaluated at `Y', but rather at some previous iterate. 
c     The new iterate, `Y_hat' will be acceptable if it leads to a sufficient 
c     reduction in the size of the natural criterion function. If the new 
c     iterate value is accepted then another fixed Jacobian step is attempted.
c     Otherwise, the code returns to taking full Newton steps. 
c      
c            This routine can fail in three ways. If the Newton iteration takes
c     too many steps, the iteration halts with `info'=1. If a singular matrix is
c     encountered during the factorization of the Newton matrix (during the call
c     to `damped_newt'), the iteration halts with `info'=2. If it is impossible 
c     for the `damped_newt' routine to obtain a suitable damping factor for the
c     Newton update (indicative of an ``effectively singular" Jacobian) the 
c     iteration halts with `info'=3. 
c-----------------------------------------------------------------------------
c***declaration of constants***
      integer  Mxs
      parameter(Mxs=10)
c
c     `Mxs' Maximum number of stages of the Runge-Kutta method.
c
c***declaration of parameters***
c   imports:
      integer                 neqns, leftbc, Nsub
      double precision        mesh(0:Nsub)
      double precision        Y(neqns*(Nsub+1))
      double precision        newtol
      integer                 maxiter
c
c     `neqns' Number of differential equations.
c     `leftbc' Number of boundary conditions at the left 
c              end of the problem interval.
c     `Nsub' Number of subintervals into which the problem
c            interval is partitioned.
c     `mesh' Set of points which partition the problem interval.
c      `Y' On input, the initial approximation for the Newton iteration.
c     `newtol' Tolerance value to be applied to the Newton iteration,
c              i.e. the iteration stops if the Newton correction 
c              is less than or equal to the user specified 
c              tolerance. A blended relative/absolute tolerance
c              is applied. If `delta_j' is the Newton correction
c              for, `Y_j', the jth component of `Y', then we require
c                     `(|delta_j| / (|Y_j|+1)) <= newtol'.
c     `maxiter' Maximum number of damped Newton iterations 
c               to be allowed during each call to this routine.
c   exports:
c     double precision      Y(neqns*(Nsub+1))
      integer               info
      double precision      k_discrete(Mxs*neqns*Nsub)
c
c     `Y' On output, when `info' = 0, the converged solution to 
c         the Newton iteration.
c     `info' Internal communication flag: 
c            0 - successful iteration converged to
c                within user specified tolerance.
c            1 - unsuccessful termination - more than `maxiter' iterations 
c                were taken without convergence being achieved.
c            2 - unsuccessful termination - a singular coefficient matrix 
c                was encountered during the attempted solution of the
c                Newton system.
c            3 - unsuccessful termination - it was impossible to obtain a 
c                suitable damping factor for the Newton update (indicative 
c                of an ``effectively singular" Jacobian) or evaluation
c                of natural criterion function overflowed.
c     `k_discrete' The discrete stages for all subintervals. The i-th set of
c                  `s*neqns' locations of `k_discrete' contains the `s
c                  stages, each of length `neqns' corresponding to the i-th
c                  subinterval. 
c   work space:
      double precision      PHI(neqns*(Nsub+1))
      double precision      delta(neqns*(Nsub+1))
      double precision      top(neqns*neqns)
      double precision      bot(neqns*neqns)
      double precision      blocks(2*neqns**2 * Nsub)
      integer               pivot(neqns*(Nsub+1))
      double precision      
     +            work(3*neqns*(Nsub+1)+neqns+2*neqns**2*(Mxs+1))
c
c     `PHI' Residual function; it has one set of `neqns' locations per 
c          subinterval. The first `leftbc' locations correspond to the 
c          left boundary conditions; the final `neqns-leftbc' locations
c          correspond to the right boundary conditions. Used as the 
c          right hand side for the Newton systems.      
c     `delta' Newton correction vector, equal to `inv(J)*PHI(Y)',
c            where J is the Jacobian evaluated at the current iterate, 
c            `Y', when a damped Newton step is taken.
c     `top' Storage of derivative information for left boundary 
c          conditions.
c     `bot' Storage of derivative information for right boundary
c          conditions. 
c     `blocks' Storage of derivative information for residual function.
c     `pivot' Pivoting information associated with the factorization
c              of J. Together `top', `blocks', and `bot' contain
c              the almost block diagonal matrix representing the 
c              coefficient matrix of the Newton system.
c     `work' Work array passed to `resid' and `fixed_jacob'.
c
c***user-supplied subroutines***
      external fsub,gsub,dfsub,dgsub
c
c     `fsub' Defines f(t,y) for the first order system
c            of differential equations, y' = f(t,y).
c     `gsub' Defines the boundary condition equations.
c     `dfsub' Defines the Jacobian of the system of
c             differential equations.
c     `dgsub' Defines the Jacobian of the boundary conditions. 
c       
c***declaration of local variables***
c
      logical            convrg
      integer            count
      double precision   delta_0_norm, lambda, g
      logical            Fixed_Jacobian
      integer            i_delta_0, i_work_d
      integer            i_Y_hat, i_PHI_hat, i_work_f
c
c     convrg' Used to check for convergence of the Newton iteration.
c     count'  Counts number of damped Newton iterations.
c     delta_0_norm' Norm of the Newton correction vector `delta_0'.
c     lambda' Damping factor used to perform a damped Newton update.
c     g' Value of natural criterion function.
c     Fixed_Jacobian' .TRUE. if the next step to be taken should be
c                     a fixed Jacobian step. .FALSE. if the next step
c                     to be taken should be a damped Newton step.
c     i_delta_0' Index into `work' for the start of the part of
c                `work' that will correspond to the temporary
c                storage array called `delta_0' within `damped_newt'.
c     i_work_d' Index to `work' for the start of the part of
c               `work' that will correspond to the work space used
c                within `damped_newt'.
c     i_Y_hat' Index to `work' for the start of the part of 
c              `work' that will correspond to the temporary storage
c              array `Y_hat' used within `fixed_jacob'.
c     i_PHI_hat' Index to `work' for the start of the part of
c                `work' that will correspond to the temporary storage
c                array `PHI_hat' used within `fixed_jacob'.
c     i_work_f' Index to `work' for the start of the part of
c               `work' that will correspond to the work space used
c                within `fixed_jacob'.
c
c***declaration of variables for /IOCNTRL/ common block***
c   imports:
      integer   print_level_0, print_level_1, print_level_2
      parameter (print_level_0 = 0, print_level_1 = 1,
     *              print_level_2 = 2)
      integer                 profile
      common /IOCNTRL/        profile
c
c     `print_level_0' No output.
c     `print_level_1' Intermediate output.
c     `print_level_2' Full output.
c     `profile' Controls output, to standard output, of profiling infor-
c               mation such as Newton iteration counts, mesh selection,
c               and defect estimates.
c-----------------------------------------------------------------------------
c      Called by: `mirkdc'
c      Calls to: `fixed_jacob',`damped_newt', `resid'         
c-----------------------------------------------------------------------------
c***initialization of work array indexes
       i_delta_0 = 1
       i_work_d    = neqns*(Nsub+1) + 1
       i_Y_hat = 1
       i_PHI_hat    = neqns*(Nsub+1) + 1
       i_work_f    = 2*neqns*(Nsub+1) + 1
c
c***initialization for Newton iteration*** 
c      The number of damped Newton iterations (equal to the number of 
c      Jacobian evaluations and factorizations) will be counted using the 
c      variable `count' which is initialized to 0. The variable `convrg', 
c      monitors whether the iteration has converged or not and is 
c      initialized to .FALSE..  
c
       count = 0
       convrg = .FALSE.
c
c      The first iteration step must be a damped Newton step.
c      This is handled by setting `Fixed_Jacobian' to .FALSE. 
c
       Fixed_Jacobian = .FALSE.
c
c      When a damped Newton step is used, it is assumed that the current 
c      iterate, `Y', and also the value of, `PHI', the residual function
c      at `Y' will be available from the previous iteration step. For the 
c      first step, which will be a damped Newton step, `Y' will be
c      available, but the residual function must be computed here in
c      preparation for the first time through the damped Newton iteration 
c      step. This is done by calling `resid' which will return with the
c      residual function stored in the vector `PHI'.
c
       call resid(neqns, leftbc, Nsub, mesh, Y, PHI, 
     +                                   k_discrete, work, fsub, gsub)
c
c      In order to predict the new damping factor, each damped Newton step
c      also assumes several other quantities are available from the previous
c      iteration step. For the first step, these will not be available. In 
c      order for the prediction algorithm to detect this situation we set 
c      `lambda=0' (The value of `lambda' can never reach 0 during the  
c      usual execution of the algorithm).
c
       lambda = 0.0d0
c
c***Newton iteration loop***
c   REPEAT-UNTIL LOOP (Implemented here using an if-goto pair)
 1                                                        continue
       if (Fixed_Jacobian) then     
c
              if (profile .GT. print_level_1) then
                write(6,*)'Fixed Jacobian step.' 
              end if
c
              call fixed_jacob(neqns,leftbc,Nsub,mesh,Y,newtol,
     +             delta,g,PHI,top,bot,blocks,pivot,Fixed_Jacobian,
     +             lambda, convrg,info, k_discrete,work(i_Y_hat),
     +             work(i_PHI_hat),work(i_work_f), fsub,gsub)
c
       else
c
c             Update `count' to monitor the number of damped Newton iterations
              count = count + 1
c
              if (profile .GT. print_level_1) then
                write(6,100)'Damped Newton step, count =',count,'.'
 100            format(1x,a27,i3,a1)
              end if
c 
              call damped_newt(neqns,leftbc,Nsub,mesh,Y,newtol,
     +               lambda,PHI,top,bot,blocks,pivot,Fixed_Jacobian,
     +               convrg,delta,delta_0_norm,g,info,k_discrete,
     +               work(i_delta_0),work(i_work_d),fsub,gsub,
     +               dfsub,dgsub)
       endif
c      End of `If (Fixed Jacobian)'
c
       if ((.NOT. convrg) .AND. (info .LE. 0)) then
c                Check to see that the maximum number of iterations has
c                not been exceeded.
                 if ( count .GT. maxiter) info = 1
       endif
c
c      Repeat the loop unless convergence has been obtained, the allowed
c      number of iterations has been exceeded, a singular matrix was
c      encountered during the Newton iteration, or a suitable damping 
c      factor has not been found (the latter three signaled by `info' > 0).
c
c   UNTIL ( convrg .OR. (info .GT. 0))
                                    if ((.NOT. convrg) .AND. 
     +                              (info.LE. 0)) go to 1
c***end of Newton iteration loop***
c
c***conclusion***
      if (profile .GT. print_level_0) then
           if (count .GT. 1) then 
             write(6,200)'After',count,' Newton steps'
           else
             write(6,300)'After',count,' Newton step'
  200        format(1x,a5,i3,a13)
  300        format(1x,a5,i3,a12)
           endif
      end if
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine newmat(neqns,leftbc,Nsub,mesh,Y,top,blocks,bot,
     +                                      k_discrete,work,dfsub,dgsub)
c
c***overview***
c
c            This routine computes the Newton matrix for the discrete system 
c     modeling the BVP.  The Newton matrix is the Jacobian of the residual
c     function.  It is evaluated at the current iterate `Y'. The Jacobian is 
c     returned in the three vectors `top',`bot', and `blocks'; `top' contains
c     the elements of the Jacobian associated with the boundary  conditions 
c     evaluated at the left endpoint - it consists of `leftbc' rows, each
c     of size `neqns', where `leftbc' is the number of boundary conditions 
c     imposed at the left endpoint, and `neqns' is the number of differential
c     equations; `bot'  contains  the  corresponding information for the right
c     endpoint. It consists of `neqns-leftbc' rows, each of length `neqns'. 
c     The vector `blocks' contains, for the i-th subinterval of the mesh which
c     partitions the problem interval, as prescribed by the vector, `mesh', the
c     blocks L_i and R_i, i=1,...,Nsub where `Nsub' is the number of subinter-
c     vals defined by the mesh.
c
c           L_i is the derivative of the i-th component of the residual function
c     with respect to y_i; R_i is the derivative of the i-th component of the
c     residual function with respect to y_i+1. These blocks, stored column 
c     by column, are placed in `blocks'  one  after the  other  in  the order
c     L_i, R_i, i=1,...,Nsub. Each such block is `neqns' rows, each of length
c     `neqns'.
c
c            This routine accesses the discrete stages for each subinterval,
c     which are computed by the `subcom' routine, and saved in the array
c     `k_discrete' provided to this routine through the parameter list. Thus 
c     a precondition for the use of this routine is that the `resid' routine 
c     has already been called, with the same `Y' input.  `work' is a work
c     array provided to `jacblk' and also used here.
c-------------------------------------------------------------------------------
c***declaration of constants***
      integer Mxs 
      parameter(Mxs=10)
c
c     `Mxs' Maximum value for the number of stages of the RK method.
c
c***declaration of parameters***
c   imports:
      integer                 neqns,leftbc,Nsub
      double precision        mesh(0:Nsub), Y(neqns*(Nsub+1))
      double precision        k_discrete(Mxs*neqns*Nsub)
c
c     `neqns' Number of differential equations.
c     `leftbc' Number of boundary conditions at the left end of the problem 
c                                                                     interval.
c     `Nsub' Number of subintervals into which the problem
c                                                 interval is partitioned.
c     `mesh' Set of points which partition the problem interval.
c     `Y' Current discrete approximation to the solution, 
c                                          evaluated at each meshpoint.
c     `k_discrete' Vector containing the discrete stages for
c                  all subintervals. The i-th set of `s*neqns' 
c                  locations of `k_discrete' contains the `s' 
c                  stages, each of length `neqns' corresponding
c                  to the i-th subinterval.
c   exports:
      double precision      top(leftbc*neqns)
      double precision      blocks(2*neqns**2*(Nsub+1))
      double precision      bot((neqns-leftbc)*neqns)
c
c     `top' Derivative information for left boundary conditions.
c     `blocks' Derivative information for residual function.
c     `bot' Derivative information for right boundary conditions.
c
c***declaration of work arrays***
      double precision         work(neqns+2*(Mxs+1)*neqns**2)
c
c     `work' Passed for working storage to the `jacblk' routine and
c            used here as the vector where the values of the 
c            derivatives of the boundary conditions are returned from
c            the user defined routine `dgsub'.
c
c***user-supplied subroutines***
      external dfsub, dgsub
c
c     `dfsub' Defines the Jacobian of the system of differential equations.
c     `dgsub' Defines the Jacobian of the boundary condition equations.
c
c***declaration of local variables***
      integer                  i, j, k, neqnsq
      integer                  i1,i2,i3,i4,n1
      integer                  toploc, botloc, workloc
      double precision         h, t
c
c     `i,j,k' Loop variables; `i' ranges from 1 to `Nsub'; `j,k' range from 1 
c             to `neqns'.
c     `i1',`i2',`i3',`i4',`n1' Indexes for the variables `Y', `k_discrete', and
c                              `blocks' .
c     `toploc', `botloc', `workloc' Offset locations in the vectors `top', 
c                                   `bot', and `work'.
c     `neqnsq' Equals `neqns*neqns'.
c
c     `h' Size of the i-th subinterval, `mesh(i) - mesh(i-1)'.
c     `t' Value of the (i-1)st mesh point, `mesh(i-1)'.
c
c***declarations of indices for work array passing to `jacblk'***
      integer        i_yr, i_dfdy, i_ur, i_dkdyi, i_dkdyip1

c     `i_yr' Starting location in `work': work(1), for `y_r', the argument for 
c            the evaluation of the `rth' stage.
c     `i_dfdy' Starting location in `work' :  work(neqns+1), for `dfdy', which
c              during the `rth' time through the main loop below, contains the 
c              derivative of the right hand side of the system of ODEs, 
c              evaluated at the arguments of the `rth' stage of the RK scheme.
c     `i_ur' Starting location in `work' : work(neqns+neqns^2+1), for `ur',
c            which is used to store the intermediate values arising in the
c            computation of each partial derivative.
c     `i_dkdyi' Starting location in `work' : work(neqns+2*neqns^2+1), for
c               `dkdyi'; the `rth' set of `neqns**2' locations of this vector
c               will contain the value of the `rth' partial derivative
c               with respect to `y_i'.
c     `i_dkdyip1' Starting location in `work'  : work(neqns+2*neqns^2+
c                                                              Mxs*neqns^2+1),
c                 for `dkdyip1'; the `rth' set of `neqns**2' locations of this 
c                 vector will contain the value of the `rth' partial derivative
c                 with respect to `y_ip1'.            
c
c***declaration of variables in common block /RK_s/***
      integer                 s
      common /RK_s/ s
c
c     `s' Number of discrete stages of the Runge-Kutta formula.
c----------------------------------------------------------------------------
c      Called by: `damped_newt'.
c      Calls to: `jacblk', `dgsub'.
c-----------------------------------------------------------------------------
c***computation of `blocks'*** 
c
      neqnsq     = neqns**2
      i_yr       =  1
      i_dfdy     =  neqns+1
      i_ur       =  neqns+neqnsq+1
      i_dkdyi    =  neqns+2*neqnsq+1
      i_dkdyip1  =  neqns+2*neqnsq+Mxs*neqnsq+1
c
c     We need to compute the Jacobian blocks L and R for each subinterval.
c
c     For each subinterval ...
c
      i1 = s*neqns
      n1= Nsub*neqns
c
      do 5 i=1, Nsub
           i2 = (i-1)*neqns
           i3 = (i-1)*i1
           i4 = 2*(i-1)*neqnsq
           h = mesh(i) - mesh(i-1)
           t = mesh(i-1)
c
c          y_i and y_i+1, the current solution approximations at the left and 
c          right endpoints of the i-th subinterval, respectively, are required 
c          for the computation of the i-th Jacobian blocks. Y((i-1)*neqns+1) 
c          is y_i; Y(i*neqns+1) is y_i+1. We will also need the `s' stage
c          evaluations associated with the `ith' subinterval. These are already
c          available in the `ith' set of `s*neqns' locations of `k_discrete',
c          `k_discrete((i-1)*s*neqns+1)'.
c
c          The Jacobian blocks are passed to the appropriate places 
c          in `blocks.' L is stored in the neqnsq locations starting at
c          location 2*(i-1)*neqnsq+1; R is stored in the neqnsq locations
c          (2*i-1)*neqnsq+1.
c
           call jacblk(neqns,h,t,Y(i2+1),Y(i2+neqns+1), 
     +          k_discrete(i3+1), blocks(i4+1), blocks(i4+neqnsq+1),
     +          work(i_yr),work(i_dfdy),work(i_ur),
     +          work(i_dkdyi), work(i_dkdyip1),dfsub)
c
 5      continue
c
c***computation of `top' and `bot'***
c
c       We now compute the parts of the Jacobian of the residual
c       function associated with the boundary conditions.
c
c       `Y(1)' is y_0, `Y(Nsub*neqns+1)' is y_Nsub. We call `dgsub' to get
c       the boundary condition derivative information. The information is
c       returned in the vector `work'. Although, `work' is viewed inside
c       `dgsub' as a square matrix of size `neqns' by `neqns', it is more
c       convenient to view it here in the usual columnwise fashion, i.e.
c       with the columns of the square matrix laid out one after the
c       other in the column vector, `work'.
c
c       Set the Jacobian to zero then call `dgsub' to obtain the 
c       non-zero entries.
c
        do 7 i = 1,neqns*neqns
          work(i) = 0.0d0
 7      continue
c
        call dgsub(neqns,Y(1),Y(n1+1),work)
c
c       The vector `top' must be assigned the values corresponding to 
c       the first `leftbc' rows of the square version of the boundary condition
c       matrix. Since `work' contains the elements of the boundary condition
c       matrix stored in a columnwise fashion, this corresponds to selecting
c       the first `leftbc' elements from each column stored within `work', and
c       copying them into the appropriate locations of the vector `top'.
c
c       The vector `top' will also store the first `leftbc' rows of the
c       boundary condition Jacobian in a columnwise fashion. Thus it will
c       consist of `neqns' vectors, each of length `leftbc', representing the
c       columns of the `leftbc' by `neqns' matrix corresponding to the 
c       derivatives of the boundary conditions at the left end of the problem
c       interval.
c       
c       A similar situation holds for the vector `bot', which will
c       contain, in a columnwise fashion, the `neqns-leftbc' by `neqns'
c       matrix, corresponding to the derivatives of the boundary conditions
c       at the right end of the problem interval. The variables `toploc',
c       `workloc', and `botloc' are used to identify, within their respective
c       vectors, the appropriate locations for the jth column, for j=1 to 
c       `neqns'.
c
        do 25 j = 1, neqns
c
           toploc = (j-1)*leftbc
           workloc = (j-1)*neqns
           botloc = (j-1)*(neqns-leftbc)
c
c          Copy the first leftbc locations of the jth "column" of `work'
c          into the leftbc locations corresponding to the jth "column" of `top'.
c
           do 10 k = 1, leftbc
            top(toploc+k) = work(workloc+k)
10         continue
c
c          Copy the remaining `neqns-leftbc' locations of the jth "column" 
c          of `work' into the `neqns-leftbc' locations corresponding to the 
c          jth "column" of `bot'.
c
           do 20 k = leftbc+1, neqns
            bot(botloc+k-leftbc) = work(workloc+k)
20         continue
c
25      continue
c
      return
      end










cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine redistribute(Nsub,mesh_current,s_hat,Nsub_star,
     +                         mesh_new)
c
c***overview***
c
c         This routine takes as input a mesh of `Nsub' subintervals which 
c      partitions the problem interval, stored in the vector `mesh_current', 
c      and, for each subinterval of this mesh, the value of the mesh function,
c      stored in the vector, `s_hat'. The mesh function, which is piecewise
c      constant, i.e. has a constant value on each subinterval of the current
c      mesh, associated with some measure of the relative defect on each 
c      subinterval (see `mesh_selector' routine for further details). The 
c      routine returns a new mesh, of `Nsub_star' points, which equidistributes
c      the integral of the mesh function over the subintervals of the new mesh
c      in the vector `mesh_new'. 
c
c         The equidistribution is performed as follows. The idea is to select 
c      new mesh points so that the area under the mesh function is equally 
c      distributed over the subintervals of the new mesh. The total area under
c      the mesh function is first computed and then divided by the total number
c      of subintervals in the new mesh, giving a quantity we call `zeta'.
c      The value of the integral of the mesh function over each subinterval of
c      the new mesh must equal `zeta'. The area under the mesh function for
c      each subinterval of the old mesh is computed and added to a sum, 
c      proceeding from left to right. Whenever the sum reaches a value beyond
c      `zeta', only the area of that fraction of the subinterval needed to 
c      reach the value of `zeta' is added. This fraction is used to decide 
c      where to place the new mesh point. The procedure beginning with the new
c      mesh point is repeated and continued until every subinterval of the
c      current mesh is processed. See [Ascher, Mattheij, Russell, Chap.9] for
c      details
c
c         This routine assumes that `Nsub_star' is less than or equal to the
c      maximum allowable number of subintervals.
c--------------------------------------------------------------------------------
c***declaration of parameters***
c      imports:
          integer                Nsub
          double precision       mesh_current(0:Nsub)
          double precision       s_hat(Nsub)
          integer                Nsub_star
c
c         `Nsub'         the number of subintervals in the input mesh
c
c         `mesh_current' the vector containing the current mesh
c
c         `s_hat'        the vector containing the mesh function 
c                        value for each subinterval.
c
c         `Nsub_star'    the number of points in the new mesh.
c      exports:
          double precision      mesh_new(0:Nsub_star)
c
c         `mesh_new'   the meshpoints making up the new mesh
c
c***declaration of local variables***
          integer i,k
          double precision sum
          double precision zeta
          double precision t
          double precision Integral
          double precision Next_Piece
c
c         `i'          loop index from 1 to Nsub_star.
c
c         `k'          loop index from 1 to Nsub.
c
c         `sum'        temporary variable used to compute `zeta'.
c
c         `zeta'       weighted sum of mesh function values, averaged over
c                      number of subintervals in the new mesh, `Nsub_star'.
c
c         `t'          used as marker to indicate the point within the 
c                      problem interval to which the integration for the
c                      determination of new mesh points has proceeded so far.
c
c         `Integral'   used to store accumulation of contributions 
c                      to the area under the mesh function for each
c                      subinterval of the new mesh.
c
c         `Next_Piece' used to store the value of each contribution
c                      to `Integral' as it arises during the determination
c                      of the position of each new mesh point.
c------------------------------------------------------------------------------
c         called by: `mesh_selector'
c------------------------------------------------------------------------------
c      Compute equidistribution constant, `zeta'. It equals the
c      weighted sum of the mesh function values divided by the desired
c      number of mesh points in the new mesh, `Nsub_star'.
c
       sum = 0.0d0
c
       do 10 k = 1, Nsub
            sum = sum + s_hat(k)*(mesh_current(k)-mesh_current(k-1))
 10    continue
c
       zeta = sum/Nsub_star
c
c      Loop to determine points of new mesh.
c
c      k gives index of the subinterval of the current mesh whose mesh
c      function contribution is being added to the integral for the ith
c      new mesh point.
c
       k = 1
       i = 0
       mesh_new(0) = mesh_current(0)
       t = mesh_current(0)
       Integral = 0.0d0
c
c      REPEAT
c                                                 --------
 1                                                continue
c                                                 --------
c
           Next_Piece = s_hat(k)*(mesh_current(k)-t)
c
           if ((Integral+Next_Piece) .GT. zeta) then
c
c              The next piece is too big. Set the next new meshpoint
c              so that it allows us to include only the appropriate
c              fraction of the current piece necessary to make the
c              integral for this new subinterval equal to `zeta'.
c
               mesh_new(i+1) = (zeta - Integral)/s_hat(k) + t
c
c              Set the "point of integration" marker, `t', to this new mesh
c              point value, set Integral to zero, and increment `i' in 
c              preparation for the determination of the next new mesh point.
c
               t = mesh_new(i+1)
               i = i + 1
               Integral = 0
c
           else
c
c              The next piece to be added to the integral is still not
c              large enough to make the value of the integral `zeta'
c              or bigger, so just add it on to the current integral
c              value and adjust the "point of integration" marker, `t'.
c
               Integral = Integral + Next_Piece
               t = mesh_current(k)
c
c              The contribution from the kth subinterval of the current
c              mesh is complete. Increment k.
c
               k = k + 1
c
           endif
c
c      UNTIL ( k .GT. Nsub)
c                                         ----------------------
                                          if (k. LE. Nsub) goto1
c                                         ----------------------
c
c      Finish up the process by setting the last new mesh point value
c      to the value of the old last mesh point.
c
       mesh_new(Nsub_star) = mesh_current(Nsub)
c                                                                               
       return
       end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine resid(neqns,leftbc,Nsub,mesh,Y,PHI,
     +                                    k_discrete,work,fsub,gsub)
c
c***overview***
c
c          This routine evaluates the residual function at the current 
c     iterate, `Y', on the current mesh, `mesh'. The residual function value 
c     is  returned  in  the vector `PHI'. The mesh has `Nsub' subintervals. 
c
c          The evaluation  of the residual function is done in two phases:
c     one involving the calculation of the components  of `PHI' associated 
c     with the boundary conditions (`leftbc' conditions at the left boundary
c     and `neqns-leftbc' at the right boundary, where `neqns' is the number of 
c     differential equations), and the other associated  with  the 
c     other `Nsub' components, each of size `neqns', corresponding  to each 
c     subinterval of the current mesh. The  components  associated with
c     the boundary  conditions  are  obtained  directly  from  the  subroutine
c     `gsub'.  Each  of  the other components  is obtained through a call
c     to the `subcom' routine.
c
c          In  the  first  part of  this  routine we perform the  necessary
c     initializations.  We next treat the components of `PHI'  related  to the
c     subintervals, and lastly deal with the components  corresponding  to the
c     boundary conditions.
c
c          A by-product of the call to `subcom' routine is the return
c     of the `s' stages that arise as intermediate quantities during the
c     computation of the part of `PHI' associated with the i-th subinterval.
c     These stage values will be used later during the call to `Newmat' and the
c     defect estimation phases of this code, and are stored here in the 
c     `k_discrete' vector which is returned from this routine through the 
c     parameter list. `work' is a work array provided to the `subcom' routine.
c     The subroutines, `fsub' and `gsub' are provided by the user to define the 
c     boundary value ODE and the boundary conditions.
c-------------------------------------------------------------------------------
c***declaration of constants***
      integer Mxs 
      parameter (Mxs=10)
c
c     `Mxs'    is the maximum number of stages.
c
c***declaration of parameters***
c   imports:
      integer                 neqns,leftbc,Nsub
      double precision        mesh(0:Nsub), Y(neqns*(Nsub+1))
c
c     `neqns' Number of differential equations.
c     `leftbc' Number of boundary conditions at the left end of the problem 
c                                                                    interval.
c     `Nsub' Number of subintervals into which the problem interval is 
c                                                               partitioned.
c     `mesh' Set of points which partition the problem interval.
c     `Y' Current discrete approximation to the solution, 
c                                        evaluated at each meshpoint.
c
c   exports:
      double precision     PHI(neqns*(Nsub+1))
      double precision     k_discrete(Mxs*neqns*Nsub)
c
c     `PHI' Value of the residual function on each subinterval.
c           It has one set of `neqns' locations per subinterval. 
c           In addition, the first `leftbc' locations correspond to the
c           left boundary conditions; the final `neqns-leftbc' locations
c           correspond to the right boundary conditions.
c     `k_discrete' Vector containing the discrete stages for
c                  all subintervals. The i-th set of `s*neqns' 
c                  locations of `k_discrete' contains the `s'
c                  stages, each of length `neqns' corresponding
c                  to the i-th subinterval.
c   work arrays:
      double precision      work(neqns)
c
c     `work' Temporary work array required by the `subcom' routine
c            and also used here to receive the values of the residual 
c            function corresponding to the boundary conditions as
c            returned from the user defined routine `gsub'.
c
c***declaration of variables in common block /RK_s/***
      integer                 s
      common /RK_s/ s
c
c     `s' Number of discrete stages of the Runge-Kutta formula.
c
c***user-supplied subroutines***
      external fsub, gsub
c
c     `fsub' Defines f(t,y) for the first order
c                      system of differential equations, y' = f(t,y).
c     `gsub' Defines the boundary condition equations.
c
c***declaration of local variables***
      integer                 i, j, i1, n1
      double precision        h, t
c       
c     `i,j'  Loop variables; `i' ranges from 1 to `Nsub';
c            `j' ranges from 1 to `neqns'.
c     `i1','n1' Indexes used for the variables `Y', `PHI', and `k_discrete'.
c     `h' Size of the i-th subinterval, `mesh(i) - mesh(i-1)'.
c     `t' Value of the i-th mesh point, `mesh(i-1)'.
c-------------------------------------------------------------------------------
c     Called by: `newiter','criterion', `fixed_jacob'.
c     Calls: `subcom', `gsub'.
c-------------------------------------------------------------------------------
c***computation of the main part of `PHI' associated with the subintervals***
c
c     For each subinterval ...
c
      do 5 i=1, Nsub
           i1 = (i-1)*neqns
           h = mesh(i) - mesh(i-1)
           t = mesh(i-1)
c
c          y_i and y_ip1, the current solution approximations at the endpoints 
c          of the i-th subinterval, are required for the computation of the i-th
c          residual component. Y((i-1)*neqns+1) is y_i, Y(i*neqns+1) is y_ip1.
c
c          We now compute the component of `PHI' associated with the i-th sub-
c          interval. Note that `PHI(leftbc+(i-1)*neqns+1)' is simply
c          the place in `PHI' where this information is stored. The stages
c          of the discrete formula, for the i-th subinterval, are returned
c          by the `subcom' routine to the appropriate part of the `k_discrete'
c          vector. 
c
           call subcom(neqns,h,t,Y(i1+1),Y(i1+neqns+1),
     +          PHI(leftbc+i1+1), k_discrete(i1*s+1),work,fsub)
c
 5    continue
c
c***boundary condition case***
c
      n1=Nsub*neqns
c
c     Y(1) is y_0, Y(Nsub*neqns+1) is y_Nsub. We call `gsub' to get the
c     boundary condition information.
c
      call gsub(neqns,Y(1),Y(n1+1),work)
c
c     We next assign the boundary condition info, stored in the array 
c     `work' to the appropriate part of `PHI'.
c
      do 10 j = 1, leftbc
           PHI(j) = work(j)
 10   continue
c
      do 15 j = leftbc+1, neqns
           PHI(n1 + j) = work(j)
 15   continue
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine RK_tableau (method)
c
c***overview****
c
c          This routine sets up all coefficients related to the discrete
c      RK formula. `method' defines the MIRK method to be used. The possible
c      values for `method' and the corresponding MIRK scheme associated with
c      each of these values are listed below.
c
c             Value of `method' |   MIRK formula
c            -------------------------------------
c                  121          |   MIRK121 scheme
c                  221          |   MIRK221 scheme
c                  232          |   MIRK232 scheme (abscissa = 0,2/3)
c                  2321         |   MIRK232 scheme (abscissa = 1,1/3)
c                  343          |   MIRK343 scheme
c                  453          |   MIRK453 scheme
c                  563          |   MIRK563 scheme
c                  5631         |   An improved MIRK563 scheme
c
c      This routine assigns values to the parameters `v',`b',`c',`x' and `s'; 
c      thus defining the MIRK formula to be used. The tableau 
c      of the method is
c
c            c(1) | v(1) |   0        0      ...     0
c            c(2) | v(2) | x(2,1)     0      ...     0
c            c(3) | v(3) | x(3,1)  x(3,2)    ...     0
c              .       .      .      .           .
c              .       .      .      .         .
c              .       .      .      .           .
c            c(s) | v(s) | x(s,1)  x(s,2)    ...     0
c            --------------------------------------------------------
c                          b(1)    b(2)      ...    b(s)
c
c          These coefficients are then available through the 
c          /RK_s/,/RK_coeff/ common blocks. 
c------------------------------------------------------------------------------
c***declaration of constants***
            integer Mxs  
            parameter (Mxs=10)
c   
c           `Mxs'    maximum number of stages in the discrete formula.

c***declaration of parameters***
c        imports:
            integer        method
c
c           `method'  defines the method to be used
c
c***declaration of variables in common block /RK_s/***
            integer                 s
            common /RK_s/ s
c
c           `s' is the number of discrete stages of the Runge-Kutta formula.

c***declaration of variables in common block /RK_coeff/***
            double precision        v(Mxs), b(Mxs)
            double precision        c(Mxs), x(Mxs**2)
            common /RK_coeff/ v, b, c, x
c         
c           `v', `b', `c',and `x' are coefficients that define 
c           the discrete Runge-Kutta formula.  See tableau above.
c
c           The first component of `c' contains the abscissa for the first 
c           new stage. The second component of `c' contains the abscissa for 
c           the second new stage, and so on.
c
c           The first component of `v' contains the blending coefficient
c           for the first new stage. The second component of `v' contains
c           the blending coefficient for the second new stage, and so on.
c
c           The matrix, X, is related to  the data structure, `x'
c           as follows (the matrix is stored columnwise in the vector):
c
c                 X(i,j) = x( (j-1)*s + i )
c          
c           `x' contains the coupling coefficients for each of the new
c           stages
c----------------------------------------------------------------------------- 
c     called by: `mirkdc'
c------------------------------------------------------------------------------
c       Translate short form of method values.
        if (method .EQ. 2) method = 221
        if (method .EQ. 4) method = 343
        if (method .EQ. 6) method = 563

c      Beginning method definition based on value of `method'.

       if (method .EQ. 121) then
c
c          Define the 1-stage, 2nd order, stage order 1 MIRK formula
c
           s = 1
c
           c(1) = 1.0d0/2.0d0
c   
           v(1) = 1.0d0/2.0d0
c  
           b(1) = 1.0d0
c
c          Column 1 of x
c 
           x(1) = 0.0d0
c
       elseif (method .EQ. 221) then
c 
c          Define the 2-stage, 2nd order, stage order 1 MIRK formula
c
           s = 2
c  
           c(1) = 0.0d0
           c(2) = 1.0d0
c 
           v(1) = 0.0d0
           v(2) = 1.0d0
c 
           b(1) = 1.0d0/2.0d0
           b(2) = 1.0d0/2.0d0
c 
c          Column 1 of x
c
           x(1) = 0.0d0
           x(2) = 0.0d0
c
c          Column 2 of x
c
           x(3) = 0.0d0
           x(4) = 0.0d0
c
       elseif (method .EQ. 232) then
c
c          Define the 2-stage, 3rd order, stage order 2 MIRK formula
c 
           s = 2
c
           c(1) = 0.0d0
           c(2) = 2.0d0/3.0d0
c
           v(1) = 0.0d0
           v(2) = 4.0d0/9.0d0
c
           b(1) = 1.0d0/4.0d0
           b(2) = 3.0d0/4.0d0
c
c          Column 1 of x
c
           x(1) = 0.0d0
           x(2) = 2.0d0/9.0d0
c
c          Column 2 of x
c
           x(3) = 0.0d0
           x(4) = 0.0d0
c
       elseif (method. EQ. 2321) then
c
c          Define the reflection of the above 2-stage, 3rd order, stage order 2
c          MIRK formula
c
           s = 2
c
           c(1) = 1.0d0
           c(2) = 1.0d0/3.0d0
 
           v(1) = 1.0d0
           v(2) = 5.0d0/9.0d0
 
           b(1) = 1.0d0/4.0d0
           b(2) = 3.0d0/4.0d0
c 
c          Column 1 of x
c
           x(1) = 0.0d0
           x(2) = -2.0d0/9.0d0
c
c          Column 2 of x
c
           x(3) = 0.0d0
           x(4) = 0.0d0
c 
       elseif (method .EQ. 343) then
c
c          Define the 3-stage, 4th order, stage order 3 MIRK formula
c
           s = 3
c      
           c(1) = 0.0d0
           c(2) = 1.0d0
           c(3) = 1.0d0/2.0d0
c
           v(1) = 0.0d0
           v(2) = 1.0d0
           v(3) = 1.0d0 / 2.0d0
c
           b(1) = 1.0d0 / 6.0d0
           b(2) = 1.0d0 / 6.0d0
           b(3) = 2.0d0 / 3.0d0
c
c          Column 1 of x
c
           x(1) = 0.0d0
           x(2) = 0.0d0
           x(3) = 1.0d0 / 8.0d0
c
c          Column 2 of x
c
           x(4) = 0.0d0
           x(5) = 0.0d0
           x(6) = -1.0d0 / 8.0d0
c
c          Column 3 of x
c
           x(7) = 0.0d0
           x(8) = 0.0d0
           x(9) = 0.0d0
c
       elseif (method. EQ. 453) then
c
c          Define the 4-stage, 5th order, stage order 3 MIRK formula
c
           s = 4
c
           c(1) = 0.0d0
           c(2) = 1.0d0
           c(3) = 3.0d0/4.0d0
           c(4) = 3.0d0/10.0d0
c
           v(1) = 0.0d0
           v(2) = 1.0d0
           v(3) = 27.0d0/32.0d0
           v(4) = 837.0d0/1250.0d0
c
           b(1) = 5.0d0/54.0d0
           b(2) = 1.0d0/14.0d0
           b(3) = 32.0d0/81.0d0
           b(4) = 250.0d0/567.0d0
c 
c          Column 1 of x
c
           x(1) = 0.0d0
           x(2) = 0.0d0
           x(3) = 3.0d0/64.0d0
           x(4) = 21.0d0/1000.0d0
c
c          Column 2 of x
c
           x(5) = 0.0d0
           x(6) = 0.0d0
           x(7) = -9.0d0/64.0d0
           x(8) = 63.0d0/5000.0d0
c
c          Column 3 of x
c
           x(9)  = 0.0d0
           x(10) = 0.0d0
           x(11) = 0.0d0
           x(12) = -252.0d0/625.0d0
c
c          Column 4 of x
c
           x(13) = 0.0d0
           x(14) = 0.0d0
           x(15) = 0.0d0
           x(16) = 0.0d0
c
       elseif (method. EQ. 4531) then
c
           s = 4
c
           c(1) = 0.0d0
           c(2) = 1.0d0
           c(3) = 3.0d0/4.0d0
           c(4) = 3.0d0/10.0d0
c
           v(1) = 0.0d0
           v(2) = 1.0d0
           v(3) = 27.0d0/32.0d0
           v(4) = 837.0d0/1250.0d0
c
           b(1) = 5.0d0/54.0d0
           b(2) = 1.0d0/14.0d0
           b(3) = 32.0d0/81.0d0
           b(4) = 250.0d0/567.0d0
c 
c          Column 1 of x
c
           x(1) = 0.0d0
           x(2) = 0.0d0
           x(3) = 3.0d0/64.0d0
           x(4) = 21.0d0/1000.0d0
c
c          Column 2 of x
c
           x(5) = 0.0d0
           x(6) = 0.0d0
           x(7) = -9.0d0/64.0d0
           x(8) = 63.0d0/5000.0d0
c
c          Column 3 of x
c
           x(9)  = 0.0d0
           x(10) = 0.0d0
           x(11) = 0.0d0
           x(12) = -252.0d0/625.0d0
c
c          Column 4 of x
c
           x(13) = 0.0d0
           x(14) = 0.0d0
           x(15) = 0.0d0
           x(16) = 0.0d0
c
       elseif (method. EQ. 563) then
c
c          Define the 5-stage, 6th order, stage order 3 MIRK formula
c
           s = 5
c
           c(1) = 0.0d0
           c(2) = 1.0d0
           c(3) = 1.0d0/4.0d0
           c(4) = 3.0d0/4.0d0
           c(5) = 1.0d0/2.0d0
c
           v(1) = 0.0d0
           v(2) = 1.0d0
           v(3) = 5.0d0 / 32.0d0
           v(4) = 27.0d0 / 32.0d0
           v(5) = 1.0d0 / 2.0d0
c
           b(1) = 7.0d0 / 90.0d0
           b(2) = 7.0d0 / 90.0d0
           b(3) = 32.0d0 / 90.0d0
           b(4) = 32.0d0 / 90.0d0
           b(5) = 12.0d0 / 90.0d0
c
c          Column 1 of x          
c
           x(1) = 0.0d0
           x(2) = 0.0d0
           x(3) = 9.0d0 / 64.0d0
           x(4) = 3.0d0 / 64.0d0
           x(5) = -5.0d0 / 24.0d0
c
c          Column 2 of x
c
           x(6) = 0.0d0
           x(7) = 0.0d0
           x(8) = -3.0d0 / 64.0d0
           x(9) = -9.0d0 / 64.0d0
           x(10) = 5.0d0 / 24.0d0
c
c          Column 3 of x
c
           x(11) = 0.0d0
           x(12) = 0.0d0
           x(13) = 0.0d0
           x(14) = 0.0d0
           x(15) = 2.0d0 / 3.0d0
c
c          Column 4 of x
c
           x(16) = 0.0d0
           x(17) = 0.0d0
           x(18) = 0.0d0
           x(19) = 0.0d0
           x(20) = -2.0d0 / 3.0d0
c         
c          Column 5 of x
c
           x(21) = 0.0d0
           x(22) = 0.0d0
           x(23) = 0.0d0
           x(24) = 0.0d0
           x(25) = 0.0d0
c
       elseif (method. EQ. 5631) then
c
c          Define the improved 5-stage, 6th order, stage order 3 MIRK formula
c
           s = 5
c
           c(1) = 0.0d0
           c(2) = 1.0d0
           c(3) = 1.0d0/2.0d0-sqrt(21.0d0)/14.0d0
           c(4) = 1.0d0/2.0d0+sqrt(21.0d0)/14.0d0
           c(5) = 1.0d0/2.0d0
c
           v(1) = 0.0d0
           v(2) = 1.0d0
           v(3) = 1.0d0/2.0d0-9.0d0*sqrt(21.0d0)/98.0d0
           v(4) = 1.0d0/2.0d0+9.0d0*sqrt(21.0d0)/98.0d0
           v(5) = 1.0d0 / 2.0d0
c
           b(1) = 1.0d0 / 20.0d0
           b(2) = 1.0d0 / 20.0d0
           b(3) = 49.0d0 / 180.0d0
           b(4) = 49.0d0 / 180.0d0
           b(5) = 16.0d0 / 45.0d0
c
c          Column 1 of x          
c
           x(1) = 0.0d0
           x(2) = 0.0d0
           x(3) = 1.0d0/14.0d0+sqrt(21.0d0)/98.0d0
           x(4) = 1.0d0/14.0d0-sqrt(21.0d0)/98.0d0
           x(5) = -5.0d0 / 128.0d0
c
c          Column 2 of x
c
           x(6) = 0.0d0
           x(7) = 0.0d0
           x(8) = -1.0d0/14.0d0+sqrt(21.0d0)/98.0d0
           x(9) = -1.0d0/14.0d0-sqrt(21.0d0)/98.0d0
           x(10) = 5.0d0 / 128.0d0
c
c          Column 3 of x
c
           x(11) = 0.0d0
           x(12) = 0.0d0
           x(13) = 0.0d0
           x(14) = 0.0d0
           x(15) = 7.0d0*sqrt(21.0d0)/128.0d0
c
c          Column 4 of x
c
           x(16) = 0.0d0
           x(17) = 0.0d0
           x(18) = 0.0d0
           x(19) = 0.0d0
           x(20) = -7.0d0*sqrt(21.0d0)/128.0d0
c         
c          Column 5 of x
c
           x(21) = 0.0d0
           x(22) = 0.0d0
           x(23) = 0.0d0
           x(24) = 0.0d0
           x(25) = 0.0d0
c
       endif
c   
       return
       end








cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sol_eval(t,z,z_prime,fail,iwork,work)
c
c***overview***
c                
c          This routine uses the current solution approximation, `Y',
c       obtained on the mesh stored in `mesh', plus the stages of the 
c       underlying discrete RK method plus some new stages associated with
c       the interpolant in order to evaluate the interpolant.  This interpolant 
c       is evaluated at `t', with the value returned in `z', and the value of 
c       the first derivative returned in `z_prime'.
c
c          The interpolant is computed as follows: Given `t', assumed to 
c       be somewhere in the problem interval, the subinterval of `mesh' that 
c       contains `t' is determined first, i.e. find i such that mesh(i-1) 
c       =< t < mesh(i) for t < mesh(Nsub) or i = Nsub for t = mesh(Nsub).
c       The fraction of the way `t' is through the ith subinterval, a quantity
c       which we call `tau',  is also determined during this step. Then the 
c       interpolant on this ith subinterval is evaluated - the interpolant has
c       the form z(mesh(i-1) + tau*hi) = yim1 + hi sum br(tau)*kr. The weights,
c       br(tau), can be computed with a call to `interp_weights'. `Nsub' is 
c       the number of subintervals.
c
c          The remaining information required to define the interpolant 
c       on the ith subinterval is the current solution approximation at the 
c       (i-1)st mesh point, `yim1', which can be obtained from `Y', and the
c       stages associated with the ith subinterval. The first `s' stages,
c       contained in `k_discrete', come from the discrete formula and
c       remaining `s_star - s' stages, contained in `k_interp', are necessary
c       for the interpolant. Each stage is a vector of length `neqns' where 
c       `neqns' is the number of differential equations. 
c
c          `k_discrete' is the vector containing the discrete stages for
c       all subintervals. The ith set of `s*neqns' locations of `k_discrete' 
c       contains the `s' stages, each of length `neqns' corresponding to the 
c       ith subinterval. `k_interp' plays a similar role for the `s_star-s'
c       stages associated with the interpolant.
c------------------------------------------------------------------------------
c***declaration of constants***
        integer Mxs 
        parameter(Mxs=10)
c
c       `Mxs'    is the maximum value for the total number of stages
c                of the interpolant.
c
c***declaration of parameters***
c      imports:
        double precision t
c       `t'     the point at which the interpolant is to be evaluated.
c      exports:
c       double precision z(neqns)
        double precision z(*)
c
c       double precision z_prime(neqns)
        double precision z_prime(*)
c
        integer fail
c
c       `z'       is the value of the interpolant at x.
c
c       `z_prime' is the value of the first derivative of the
c                 interpolant at x.
c
c       `fail'    if fail = 1 then the point `t', at which the interpolant
c                 is to be evaluated, is outside the problem interval and 
c                 therefore `sol_eval' will not produce a solution
c                 approximation. If fail = 0 then `t' is a valid meshpoint.
c
c***LAYOUT OF WORK ARRAYS***
c
c        All of the input info to this routine except for the evaluation 
c        point `t' is contained within the `iwork' and `work' arrays.
c       
         integer iwork(5)
c       
c        iwork(1) = neqns
         integer neqns
c        `neqns' is the number of differential equations.
c
c        iwork(2) = Nsub
         integer Nsub
c        `Nsub' is the number of subintervals.
c
c        iwork(3) = s
         integer s
c        `s' is the number of discrete stages.
c
c        iwork(4) = s_star
         integer s_star
c        `s_star'  is the total number of stages required to form the
c                  interpolant on each subinterval. It includes all the
c                  stages of the discrete formula plus the additional 
c                  stages required for the interpolant.
c
c        iwork(5) = method
         integer method
c        `method' defines the Runge-Kutta method that was used to compute the
c                 solution. It is passed via the /rk_method/ common block to
c                 the `interp_weights' routine.
c
c
c        let MnN = Mxs*neqns*Nsub.
c        
c        double precision work(2*Mxs*neqns*Nsub+(Nsub+1)+(Nsub+1)*neqns)
         double precision work(*)
c
c        work(1..MnN) = k_discrete(Mxs*neqns*Nsub)
c        `k_discrete' contains the `s' discrete stage values for each 
c                                                                subinterval.
c        work(MnN+1..2*MnN) = k_interp(Mxs*neqns*Nsub)
c        `k_interp' contains the `s_star-s' extra stage values for each
c                                 subinterval need to form the interpolant.
c
c        work(2*MnN+1..2*MnN+(Nsub+1)) = mesh(0:Nsub)
c        `mesh'  is the list of points defining the subintervals
c
c        work(2*MnN+(Nsub+1)+1..2*MnN+(Nsub+1)+(Nsub+1)*neqns)) = 
c                                                   Y(neqns*(Nsub+1))
c        `Y'     contains the discrete solution at the mesh points;
c                the ith set of `neqns' locations of `Y' contains 
c                the solution approximation at the ith mesh point.
c         
c***declaration of local variables***
        integer i
        double precision hi
        double precision tau
        double precision weights(Mxs)
        double precision weights_prime(Mxs)
        double precision a, b
c
c       `i'             index for subinterval containing `t', i ranges
c                       from 1 to `Nsub'
c
c       `hi'            the size of the ith subinterval.
c
c       `tau'           equal to (t-mesh(i-1))/hi
c
c       `weights'       contains the weight polynomials of the interpolant
c                       evaluated at `tau'.
c
c       `weights_prime' contains the first derivatives of the weight 
c                       polynomials of the interpolant evaluated at `tau'.
c 
c       `a'             the left hand endpoint of the interval.
c       `b'             the right hand endpoint of the interval.
c
c***declaration of work array indexes***
        integer i_k_discrete 
        integer i_k_interp
        integer i_mesh, i_Y

c       `i_k_discrete' = 1; `i_k_interp' = Mxs*neqns*Nsub+1;
c       `i_mesh' = 2*Mxs*neqns*Nsub+1; `i_Y' = 2*Mxs*neqns*Nsub+(Nsub+1)+1

c***declaration of variables in /rk_method/ common block***
c      exports:
c       integer   method
        common   /rk_method/ method

c       `method'  defines the method to be used
c
c       This value is passed through the rk_method common block
c       to the `interp_weights' routine. 
c-----------------------------------------------------------------------------
c      called by:`main',    user's main program subsequent to the
c                           completion of the call to mirkdc.
c
c      calls to: `interval', `interp_weights', `sum_stages'
c------------------------------------------------------------------------------
c         Setup access to info in `iwork' and `work'

          neqns = iwork(1)
          Nsub = iwork(2)
          s = iwork(3)
          s_star = iwork(4)
          method = iwork(5)

          i_k_discrete = 1
          i_k_interp = Mxs*neqns*Nsub+1
          i_mesh = 2*Mxs*neqns*Nsub+1
          i_Y = 2*Mxs*neqns*Nsub+(Nsub+1)+1
          
          a = work(i_mesh)
          b = work(i_mesh+(Nsub))
c
          fail = 0
          if ((t.LT.a).or.(t.GT.b)) then
           fail = 1
          else
c          Call `interval' to determine the subinterval containing t.
c          The index of the subinterval is returned in `i'.
c
           call interval(Nsub,work(i_mesh),t,i)
c                                                                        
c          Compute `tau' for use in the evaluation of the weight polynomials.
c
           hi = work(i_mesh+i) - work(i_mesh+(i-1))
           tau = (t - work(i_mesh+(i-1)))/hi
c
c          Setup the weights, evaluated at `tau'.
c
           call interp_weights(s_star, tau, weights,weights_prime)
c
c          We need to evaluate the interpolant and its derivative at `tau'.
c   
c          This is done by taking weighted sums of the stages, using the weights
c          we have just computed, and the stages corresponding to the ith
c          subinterval, and calling the `sum_stages' routine.
c
           call sum_stages(neqns,hi,work(i_Y+(i-1)*neqns),s,
     +               work(i_k_discrete+(i-1)*s*neqns),s_star,
     +               work(i_k_interp+(i-1)*(s_star-s)*neqns),
     +                     weights, weights_prime, z, z_prime) 
c
          end if
          end









ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine subcom(neqns,hi,ti,yi,yip1,phii,k,y,fsub)
c
c***overview***
c
c  ***Subcom Routine(Calculation of subcomponents of residual function)***
c
c            The purpose of this routine is to evaluate the component
c       of PHI, the residual function, associated with the i-th sub-
c       interval. This component of PHI is returned through the parameter 
c       `phii'.
c                                          s
c                                          --
c                     phii = yip1-yi- hi * <  br*kr
c                                          --
c                                          r=1
c
c       The evaluation takes place for the current iterate `Y'.
c       However, the i-th component of PHI depends only on the solution 
c       information directly associated with the i-th subinterval, namely
c       yi and yi+1, which have been assigned to `yi' and `yip1'.  Also    
c       the size of the i-th subinterval is assigned to `hi' and the left
c       hand endpoint of the ith subinterval is assigned to `ti'.
c 
c            The calculation uses the s-stage MIRK formula defined by the 
c       values of `s', `c',`v', `b' and `x' (as passed in from the 
c       /RK_s/,/RK_coeff/ commom blocks), The solution approximation vectors, 
c       `yi' and `yip1', as well as `phii', are of size `neqns', the number 
c       of differential equations.
c
c            During the computation of `phii', we will compute `s' 
c       intermediate vectors, kr, r=1..s, each of size `neqns', and called 
c       the stages of the formula. The stages are stored in the  vector `k'  
c       in order, and are returned by this routine for later use during the 
c       call to `newmat' and during the interpolation phase.
c       `y' provides a work array of size `neqns'.
c-----------------------------------------------------------------------------
c***declaration of constants***
         integer Mxs 
         parameter (Mxs=10)
c
c         `Mxs' maximum number of stages. 
c
c***declaration of parameters***
c        imports :
           integer                    neqns 
           double precision           hi, ti, yi(neqns), yip1(neqns)
c
c         `neqns' is the number of differential equations
c
c         `hi'    is the size of the ith subinterval
c
c         `ti'    is the lefthand endpoint of the ith subinterval
c
c         `yi'    is the solution approximation value at the left endpoint
c
c         `yip1'  is the solution approxiamtion at the right endpoint
c        exports :
           double precision           phii(neqns)
           double precision           k(Mxs*neqns)
c
c         `phii' is the value of the component of the residual function
c                associated with the ith subinterval
c
c         `k'    stores the `s' stage values k1, k2, ..., ks, each of size
c                `neqns'
c
c***declaration of work array***
           double precision           y(neqns)

c          `y'     a work array for temporary storage.
c
c***user-supplied subroutines***
           external fsub
c
c          `fsub' - This routine defines f(t,y) for the first order system 
c                   of differential equations, y' = f(t,y).
c
c***declaration of local variables***
           integer                    j, r, l
           integer                    i1,i2
           double precision           tr
c
c          `j,r,l' are loop index variables; `j,r' range from 1 to `s'
c                 `l' ranges from 1 to `neqns'.
c
c          `i1',`i2' are indexes used in calculating ur.
c
c          `tr'    point at which rth stage of Runge-Kutta method is 
c                  evaluated.
c
c***declaration of variables in common block /RK_s/***
            integer                 s
            common /RK_s/ s
c
c           `s' is the number of discrete stages of the Runge-Kutta formula.

c***declaration of variables in common block /RK_coeff/***
            double precision        v(Mxs), b(Mxs)
            double precision        c(Mxs), x(Mxs**2)
            common /RK_coeff/ v, b, c, x
c         
c           `v', `b', `c',and `x' are coefficients that define 
c           the discrete Runge-Kutta formula.  
c
c           The first component of `c' contains the abscissa for the first 
c           new stage. The second component of `c' contains the abscissa for 
c           the second new stage, and so on.
c
c           The first component of `v' contains the blending coefficient
c           for the first new stage. The second component of `v' contains
c           the blending coefficient for the second new stage, and so on.
c
c           The matrix, X, is related to  the data structure, `x'
c           as follows (the matrix is stored columnwise in the vector):
c
c                 X(i,j) = x( (j-1)*s + i )
c          
c           `x' contains the coupling coefficients for each of the new
c           stages
c------------------------------------------------------------------------------
c      called by: `resid'
c      calls to: `fsub'
c------------------------------------------------------------------------------
c***COMPUTATION OF `phii'***
c
c       For each stage of the Runge-Kutta method, we compute the
c       argument of the function f, the right hand side of the ODE system,
c       and then evaluate f there by calling `fsub', with the results being
c       stored in the appropiate part of the vector `k'. The vector `y' 
c       provides  storage for accumulation of the argument for f.
c
        do  7  r = 1, s
c
           do 2 l = 1, neqns
               y(l) = 0.0d0
 2         continue 
c
c          Add the weighted sum of the previous stages 
c 
           do 4 j = 1, r-1
            i1 = (j-1)*s + r
            i2 = (j-1)*neqns
c
             do 3 l = 1, neqns
                 y(l) = y(l) + x(i1)*k(i2+l)
 3           continue
c            
 4         continue
c
c         Multiply by hi
c
          do 5 l = 1, neqns
            y(l) = hi*y(l)
 5        continue

c          Add (1-v(r))*yi + v(r)*yip1 
c
           do  6  l = 1, neqns
              y(l) = (1.0d0 - v(r))*yi(l) + v(r)*yip1(l) + y(l)
 6         continue
c
c             The second argument to f has been accumulated in the vector `y'.
c             Compute the abscissa for the argument to fsub.
c
              tr = ti + c(r)*hi
c
c            Evaluate the function f to get kr, the rth stage.
c            The user routine which evaluates the right hand side of the
c            system of ODE's must be called `fsub'
c
              call fsub(neqns, tr, y, k((r-1)*neqns+1))
c
c            The rth stage has now been computed and stored in the rth 
c            vector location of k
c
 7      continue
c
c      All `s' stages of the formula have now been computed. We can now
c      take their weighted sum in order to complete the calculation of `phii'.
c
        do  8 l = 1, neqns
            phii(l) = 0.0d0
 8      continue
c
        do 9 r = 1, s
         i1 = (r-1)*neqns
         do 10 l = 1, neqns
            phii(l) = phii(l) - b(r)*k(i1+l)
 10         continue
 9      continue
c
c       Multiply phii by hi
c
        do 11 l = 1, neqns
          phii(l) = hi*phii(l)
 11     continue
c
c       Add  (yip1 - yi) to phii
c
        do  12  l = 1, neqns
           phii(l) = yip1(l) - yi(l) + phii(l)
 12     continue
c
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sum_stages(neqns,hi,yim1,s,ki_discrete,
     +            s_star,ki_interp,weights,weights_prime,z,z_prime)
c
c***overview***
c
c          Given the stages and corresponding weights, this routine
c      computes the value of the interpolant, `z', and its first derivative, 
c      `z_prime', as follows: 
c            z = yim1 + hi sum weights(r)*ki(r),
c            z_prime = sum weights_prime(r)*ki(r),
c      where the ki(r)'s are the stages. The first `s' stages, contained
c      in `ki_discrete', come from the discrete formula; the remaining 
c      `s_star - s' stages, contained in `ki_interp', are necessary for
c      the interpolant only. Each stage is a vector of length `neqns' where
c      `neqns' is the number of differential equations. `hi' is the size of
c      the ith subinterval and `yim1'  is the discrete solution at the left
c      endpoint of the ith subinterval. `weights' contains the values of the 
c      weight polynomials of the interpolant and `weights_prime' contains the
c      values of the first derivative of the weight polynomials.
c------------------------------------------------------------------------------
c***declaration of parameters***
c      imports:
            integer neqns
            double precision hi
            double precision yim1(neqns)
            integer s
            double precision ki_discrete(s*neqns)
            integer s_star
            double precision ki_interp((s_star-s)*neqns)
            double precision weights(s_star)
            double precision weights_prime(s_star)
c
c           `neqns'  is the number of differential equations.
c
c           `hi'     is the size of the current (ith) subinterval.
c
c           `yim1'   is the discrete solution at the left endpoint
c                    of the ith subinterval.
c
c           `s'      is the number of discrete stages.
c
c           `ki_discrete'   the set of `s' discrete stages associated with 
c                           the ith subinterval. 
c
c           `s_star'        is the total number of stages.
c
c           `ki_interp'     the set of `s_star-s' interpolation stages 
c                           associated with the ith subinterval. 
c
c           `weights'       the weights associated with the stages needed to
c                           evaluate the interpolant at the appropriate
c                           point. The weights have aleady been evaluated at
c                           that point.
c           
c           `weights_prime' the weights associated with the stages needed to
c                           form the first derivative of the interpolant
c                           at the appropriate point. The weights have already
c                           been evaluated at this point.
c      exports:
            double precision z(neqns)
            double precision z_prime(neqns)
c
c           `z'       is the value of the interpolant.
c
c           `z_prime' is the value of the first derivative of the 
c                     interpolant.
c
c***declaration of local variables***
            integer j
c
c           `j'       loop index from 1 to 's', or 1 to 's_star-s'
c---------------------------------------------------------------------------
c      Called by: `defect_estimate', `interp_eval'
c      Calls to : `dcopy', `daxpy'
c---------------------------------------------------------------------------
c          Initialize all vectors to be used in this computation to zero.
c
           call dcopy(neqns,0.0d0,0,z,1)
           call dcopy(neqns,0.0d0,0,z_prime,1)
c
c          We can now sum the appropriately weighted stages to compute
c          the interpolant `z' and the first derivative of the 
c          interpolant, `z_prime'.
c
c          Accumulate weighted sum of discrete stages for `z' and `z_prime'.
c
           do 10 j = 1, s
c
                call daxpy(neqns,weights(j), 
     +                      ki_discrete((j-1)*neqns+1),1,z,1)
                call daxpy(neqns,weights_prime(j), 
     +                      ki_discrete((j-1)*neqns+1),1,z_prime,1)
 10        continue
c
c          Accumulate weighted sum of new interpolant stages for `z' 
c          and `z_prime'.
c
           do 20 j = 1,(s_star-s)
c
            call daxpy(neqns,weights(s+j), 
     +                      ki_interp((j-1)*neqns+1),1,z,1)

            call daxpy(neqns,weights_prime(s+j), 
     +                      ki_interp((j-1)*neqns+1),1,z_prime,1)

 20        continue
c
c          The calculation of `z_prime' is now complete. 
c          Scale the weighted sum of stages by hi, and then add
c          on yim1. This will complete the calculation of `z'. 
c 
           call dscal(neqns,hi,z,1)
           call daxpy(neqns,1.0d0,yim1,1,z,1)
c            
           return
           end




