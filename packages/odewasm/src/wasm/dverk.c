/* ../reference/netlib/dverk.f -- translated by f2c (version 20240504).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b32 = .16666666666666666;
static integer c__1 = 1;

/* Subroutine */ int dverk_(integer *n, S_fp fcn, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *tol, integer *ind, doublereal *c__, 
	integer *nw, doublereal *w)
{
    /* Format strings */
    static char fmt_505[] = "(///\0020\002,\002computation stopped in dverk "
	    "with the following values -  \002/\0020\002,\002ind =\002,i4,5x"
	    ",\002tol  =\002,1pd13.6,5x,\002x         =\002,1pd22.15/\002 "
	    "\002,\002n   =\002,i4,5x,\002hmin =\002,1pd13.6,5x,\002xend     "
	    " =\002,1pd22.15/\002 \002,\002nw  =\002,i4,5x,\002hmax =\002,1pd"
	    "13.6,5x,\002prev xend =\002,1pd22.15/\0020\002,14x,\002no of suc"
	    "cessful steps    =\002,0pf8.0/\002 \002,14x,\002no of successive"
	    " failures =\002,0pf8.0/\002 \002,14x,\002no of function evals   "
	    "   =\002,0pf8.0/\0020\002,\002the components of y are\002//(\002 "
	    "\002,1p5d24.15))";

    /* System generated locals */
    integer w_dim1, w_offset, i__1;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), d_sign(doublereal *, 
	    doublereal *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer k;
    doublereal temp;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, fmt_505, 0 };



/* *********************************************************************** */
/*                                                                      * */
/* note added 11/14/85.                                                 * */
/*                                                                      * */
/* if you discover any errors in this subroutine, please contact        * */
/*                                                                      * */
/*        kenneth r. jackson                                            * */
/*        department of computer science                                * */
/*        university of toronto                                         * */
/*        toronto, ontario,                                             * */
/*        canada   m5s 1a4                                              * */
/*                                                                      * */
/*        phone: 416-978-7075                                           * */
/*                                                                      * */
/*        electronic mail:                                              * */
/*        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     * */
/*        csnet:  krj@toronto                                           * */
/*        arpa:   krj.toronto@csnet-relay                               * */
/*        bitnet: krj%toronto@csnet-relay.arpa                          * */
/*                                                                      * */
/* dverk is written in fortran 66.                                      * */
/*                                                                      * */
/* the constants dwarf and rreb -- c(10) and c(11), respectively -- are * */
/* set for a  vax  in  double  precision.  they  should  be  reset,  as * */
/* described below, if this program is run on another machine.          * */
/*                                                                      * */
/* the c array is declared in this subroutine to have one element only, * */
/* although  more  elements  are  referenced  in this subroutine.  this * */
/* causes some compilers to issue warning messages.  there is,  though, * */
/* no  error  provided  c is declared sufficiently large in the calling * */
/* program, as described below.                                         * */
/*                                                                      * */
/* the following external statement  for  fcn  was  added  to  avoid  a * */
/* warning  message  from  the  unix  f77 compiler.  the original dverk * */
/* comments and code follow it.                                         * */
/*                                                                      * */
/* *********************************************************************** */


/* *********************************************************************** */
/*                                                                      * */
/*     purpose - this is a runge-kutta  subroutine  based  on  verner's * */
/* fifth and sixth order pair of formulas for finding approximations to * */
/* the solution of  a  system  of  first  order  ordinary  differential * */
/* equations  with  initial  conditions. it attempts to keep the global * */
/* error proportional to  a  tolerance  specified  by  the  user.  (the * */
/* proportionality  depends  on the kind of error control that is used, * */
/* as well as the differential equation and the range of integration.)  * */
/*                                                                      * */
/*     various options are available to the user,  including  different * */
/* kinds  of  error control, restrictions on step sizes, and interrupts * */
/* which permit the user to examine the state of the  calculation  (and * */
/* perhaps make modifications) during intermediate stages.              * */
/*                                                                      * */
/*     the program is efficient for non-stiff systems.  however, a good * */
/* variable-order-adams  method  will probably be more efficient if the * */
/* function evaluations are very costly.  such a method would  also  be * */
/* more suitable if one wanted to obtain a large number of intermediate * */
/* solution values by interpolation, as might be the case  for  example * */
/* with graphical output.                                               * */
/*                                                                      * */
/*                                    hull-enright-jackson   1/10/76    * */
/*                                                                      * */
/* *********************************************************************** */
/*                                                                      * */
/*     use - the user must specify each of the following                * */
/*                                                                      * */
/*     n  number of equations                                           * */
/*                                                                      * */
/*   fcn  name of subroutine for evaluating functions - the  subroutine * */
/*           itself must also be provided by the user - it should be of * */
/*           the following form                                         * */
/*              subroutine fcn(n, x, y, yprime)                         * */
/*              integer n                                               * */
/*              double precision x, y(n), yprime(n)                     * */
/*                      *** etc ***                                     * */
/*           and it should evaluate yprime, given n, x and y            * */
/*                                                                      * */
/*     x  independent variable - initial value supplied by user         * */
/*                                                                      * */
/*     y  dependent variable - initial values of components y(1), y(2), * */
/*           ..., y(n) supplied by user                                 * */
/*                                                                      * */
/*  xend  value of x to which integration is to be carried out - it may * */
/*           be less than the initial value of x                        * */
/*                                                                      * */
/*   tol  tolerance - the subroutine attempts to control a norm of  the * */
/*           local  error  in  such  a  way  that  the  global error is * */
/*           proportional to tol. in some problems there will be enough * */
/*           damping  of  errors, as well as some cancellation, so that * */
/*           the global error will be less than tol. alternatively, the * */
/*           control   can   be  viewed  as  attempting  to  provide  a * */
/*           calculated value of y at xend which is the exact  solution * */
/*           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) * */
/*           is proportional to tol.  (the norm  is  a  max  norm  with * */
/*           weights  that  depend on the error control strategy chosen * */
/*           by the user.  the default weight for the k-th component is * */
/*           1/max(1,abs(y(k))),  which therefore provides a mixture of * */
/*           absolute and relative error control.)                      * */
/*                                                                      * */
/*   ind  indicator - on initial entry ind must be set equal to  either * */
/*           1  or  2. if the user does not wish to use any options, he * */
/*           should set ind to 1 - all that remains for the user to  do * */
/*           then  is  to  declare c and w, and to specify nw. the user * */
/*           may also  select  various  options  on  initial  entry  by * */
/*           setting ind = 2 and initializing the first 9 components of * */
/*           c as described in the next section.  he may also  re-enter * */
/*           the  subroutine  with ind = 3 as mentioned again below. in * */
/*           any event, the subroutine returns with ind equal to        * */
/*              3 after a normal return                                 * */
/*              4, 5, or 6 after an interrupt (see options c(8), c(9))  * */
/*              -1, -2, or -3 after an error condition (see below)      * */
/*                                                                      * */
/*     c  communications vector - the dimension must be greater than or * */
/*           equal to 24, unless option c(1) = 4 or 5 is used, in which * */
/*           case the dimension must be greater than or equal to n+30   * */
/*                                                                      * */
/*    nw  first dimension of workspace w -  must  be  greater  than  or * */
/*           equal to n                                                 * */
/*                                                                      * */
/*     w  workspace matrix - first dimension must be nw and second must * */
/*           be greater than or equal to 9                              * */
/*                                                                      * */
/*     the subroutine  will  normally  return  with  ind  =  3,  having * */
/* replaced the initial values of x and y with, respectively, the value * */
/* of xend and an approximation to y at xend.  the  subroutine  can  be * */
/* called  repeatedly  with new values of xend without having to change * */
/* any other argument.  however, changes in tol, or any of the  options * */
/* described below, may also be made on such a re-entry if desired.     * */
/*                                                                      * */
/*     three error returns are also possible, in which  case  x  and  y * */
/* will be the most recently accepted values -                          * */
/*     with ind = -3 the subroutine was unable  to  satisfy  the  error * */
/*        requirement  with a particular step-size that is less than or * */
/*        equal to hmin, which may mean that tol is too small           * */
/*     with ind = -2 the value of hmin  is  greater  than  hmax,  which * */
/*        probably  means  that the requested tol (which is used in the * */
/*        calculation of hmin) is too small                             * */
/*     with ind = -1 the allowed maximum number of fcn evaluations  has * */
/*        been  exceeded,  but  this  can only occur if option c(7), as * */
/*        described in the next section, has been used                  * */
/*                                                                      * */
/*     there are several circumstances that will cause the calculations * */
/* to  be  terminated,  along with output of information that will help * */
/* the user determine the cause of  the  trouble.  these  circumstances * */
/* involve  entry with illegal or inconsistent values of the arguments, * */
/* such as attempting a normal  re-entry  without  first  changing  the * */
/* value of xend, or attempting to re-enter with ind less than zero.    * */
/*                                                                      * */
/* *********************************************************************** */
/*                                                                      * */
/*     options - if the subroutine is entered with ind = 1, the first 9 * */
/* components of the communications vector are initialized to zero, and * */
/* the subroutine uses only default values  for  each  option.  if  the * */
/* subroutine  is  entered  with ind = 2, the user must specify each of * */
/* these 9 components - normally he would first set them all  to  zero, * */
/* and  then  make  non-zero  those  that  correspond to the particular * */
/* options he wishes to select. in any event, options may be changed on * */
/* re-entry  to  the  subroutine  -  but if the user changes any of the * */
/* options, or tol, in the course of a calculation he should be careful * */
/* about  how  such changes affect the subroutine - it may be better to * */
/* restart with ind = 1 or 2. (components 10 to 24 of c are used by the * */
/* program  -  the information is available to the user, but should not * */
/* normally be changed by him.)                                         * */
/*                                                                      * */
/*  c(1)  error control indicator - the norm of the local error is  the * */
/*           max  norm  of  the  weighted  error  estimate  vector, the * */
/*           weights being determined according to the value of c(1) -  * */
/*              if c(1)=1 the weights are 1 (absolute error control)    * */
/*              if c(1)=2 the weights are 1/abs(y(k))  (relative  error * */
/*                 control)                                             * */
/*              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) * */
/*                 (relative  error  control,  unless abs(y(k)) is less * */
/*                 than the floor value, abs(c(2)) )                    * */
/*              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) * */
/*                 (here individual floor values are used)              * */
/*              if c(1)=5 the weights are 1/abs(c(k+30))                * */
/*              for all other values of c(1), including  c(1) = 0,  the * */
/*                 default  values  of  the  weights  are  taken  to be * */
/*                 1/max(1,abs(y(k))), as mentioned earlier             * */
/*           (in the two cases c(1) = 4 or 5 the user must declare  the * */
/*           dimension of c to be at least n+30 and must initialize the * */
/*           components c(31), c(32), ..., c(n+30).)                    * */
/*                                                                      * */
/*  c(2)  floor value - used when the indicator c(1) has the value 3    * */
/*                                                                      * */
/*  c(3)  hmin specification - if not zero, the subroutine chooses hmin * */
/*           to be abs(c(3)) - otherwise it uses the default value      * */
/*              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     * */
/*           where dwarf is a very small positive  machine  number  and * */
/*           rreb is the relative roundoff error bound                  * */
/*                                                                      * */
/*  c(4)  hstart specification - if not zero, the subroutine  will  use * */
/*           an  initial  hmag equal to abs(c(4)), except of course for * */
/*           the restrictions imposed by hmin and hmax  -  otherwise it * */
/*           uses the default value of hmax*(tol)**(1/6)                * */
/*                                                                      * */
/*  c(5)  scale specification - this is intended to be a measure of the * */
/*           scale of the problem - larger values of scale tend to make * */
/*           the method more reliable, first  by  possibly  restricting * */
/*           hmax  (as  described  below) and second, by tightening the * */
/*           acceptance requirement - if c(5) is zero, a default  value * */
/*           of  1  is  used.  for  linear  homogeneous  problems  with * */
/*           constant coefficients, an appropriate value for scale is a * */
/*           norm  of  the  associated  matrix.  for other problems, an * */
/*           approximation to  an  average  value  of  a  norm  of  the * */
/*           jacobian along the trajectory may be appropriate           * */
/*                                                                      * */
/*  c(6)  hmax specification - four cases are possible                  * */
/*           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            * */
/*              min(abs(c(6)),2/abs(c(5)))                              * */
/*           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) * */
/*           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            * */
/*              2/abs(c(5))                                             * */
/*           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value * */
/*              of 2                                                    * */
/*                                                                      * */
/*  c(7)  maximum number of function evaluations  -  if  not  zero,  an * */
/*           error  return with ind = -1 will be caused when the number * */
/*           of function evaluations exceeds abs(c(7))                  * */
/*                                                                      * */
/*  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will * */
/*           interrupt   the  calculations  after  it  has  chosen  its * */
/*           preliminary value of hmag, and just before choosing htrial * */
/*           and  xtrial  in  preparation for taking a step (htrial may * */
/*           differ from hmag in sign, and may  require  adjustment  if * */
/*           xend  is  near) - the subroutine returns with ind = 4, and * */
/*           will resume calculation at the point  of  interruption  if * */
/*           re-entered with ind = 4                                    * */
/*                                                                      * */
/*  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will * */
/*           interrupt   the  calculations  immediately  after  it  has * */
/*           decided whether or not to accept the result  of  the  most * */
/*           recent  trial step, with ind = 5 if it plans to accept, or * */
/*           ind = 6 if it plans to reject -  y(*)  is  the  previously * */
/*           accepted  result, while w(*,9) is the newly computed trial * */
/*           value, and w(*,2) is the unweighted error estimate vector. * */
/*           the  subroutine  will  resume calculations at the point of * */
/*           interruption on re-entry with ind = 5 or 6. (the user  may * */
/*           change ind in this case if he wishes, for example to force * */
/*           acceptance of a step that would otherwise be rejected,  or * */
/*           vice versa. he can also restart with ind = 1 or 2.)        * */
/*                                                                      * */
/* *********************************************************************** */
/*                                                                      * */
/*  summary of the components of the communications vector              * */
/*                                                                      * */
/*     prescribed at the option       determined by the program         * */
/*           of the user                                                * */
/*                                                                      * */
/*                                    c(10) rreb(rel roundoff err bnd)  * */
/*     c(1) error control indicator   c(11) dwarf (very small mach no)  * */
/*     c(2) floor value               c(12) weighted norm y             * */
/*     c(3) hmin specification        c(13) hmin                        * */
/*     c(4) hstart specification      c(14) hmag                        * */
/*     c(5) scale specification       c(15) scale                       * */
/*     c(6) hmax specification        c(16) hmax                        * */
/*     c(7) max no of fcn evals       c(17) xtrial                      * */
/*     c(8) interrupt no 1            c(18) htrial                      * */
/*     c(9) interrupt no 2            c(19) est                         * */
/*                                    c(20) previous xend               * */
/*                                    c(21) flag for xend               * */
/*                                    c(22) no of successful steps      * */
/*                                    c(23) no of successive failures   * */
/*                                    c(24) no of fcn evals             * */
/*                                                                      * */
/*  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        * */
/*                                                                      * */
/* *********************************************************************** */
/*                                                                      * */
/*  an overview of the program                                          * */
/*                                                                      * */
/*     begin initialization, parameter checking, interrupt re-entries   * */
/*  ......abort if ind out of range 1 to 6                              * */
/*  .     cases - initial entry, normal re-entry, interrupt re-entries  * */
/*  .     case 1 - initial entry (ind .eq. 1 or 2)                      * */
/*  v........abort if n.gt.nw or tol.le.0                               * */
/*  .        if initial entry without options (ind .eq. 1)              * */
/*  .           set c(1) to c(9) equal to zero                          * */
/*  .        else initial entry with options (ind .eq. 2)               * */
/*  .           make c(1) to c(9) non-negative                          * */
/*  .           make floor values non-negative if they are to be used   * */
/*  .        end if                                                     * */
/*  .        initialize rreb, dwarf, prev xend, flag, counts            * */
/*  .     case 2 - normal re-entry (ind .eq. 3)                         * */
/*  .........abort if xend reached, and either x changed or xend not    * */
/*  .        re-initialize flag                                         * */
/*  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    * */
/*  v        transfer control to the appropriate re-entry point.......  * */
/*  .     end cases                                                  .  * */
/*  .  end initialization, etc.                                      .  * */
/*  .                                                                v  * */
/*  .  loop through the following 4 stages, once for each trial step .  * */
/*  .     stage 1 - prepare                                          .  * */
/* ***********error return (with ind=-1) if no of fcn evals too great .  * */
/*  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  * */
/*  .        calc hmin, scale, hmax                                  .  * */
/* ***********error return (with ind=-2) if hmin .gt. hmax            .  * */
/*  .        calc preliminary hmag                                   .  * */
/* ***********interrupt no 1 (with ind=4) if requested.......re-entry.v  * */
/*  .        calc hmag, xtrial and htrial                            .  * */
/*  .     end stage 1                                                .  * */
/*  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  * */
/*  .     stage 3 - calc the error estimate                          .  * */
/*  .     stage 4 - make decisions                                   .  * */
/*  .        set ind=5 if step acceptable, else set ind=6            .  * */
/* ***********interrupt no 2 if requested....................re-entry.v  * */
/*  .        if step accepted (ind .eq. 5)                              * */
/*  .           update x, y from xtrial, ytrial                         * */
/*  .           add 1 to no of successful steps                         * */
/*  .           set no of successive failures to zero                   * */
/* **************return(with ind=3, xend saved, flag set) if x .eq. xend * */
/*  .        else step not accepted (ind .eq. 6)                        * */
/*  .           add 1 to no of successive failures                      * */
/* **************error return (with ind=-3) if hmag .le. hmin            * */
/*  .        end if                                                     * */
/*  .     end stage 4                                                   * */
/*  .  end loop                                                         * */
/*  .                                                                   * */
/*  begin abort action                                                  * */
/*     output appropriate  message  about  stopping  the  calculations, * */
/*        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, * */
/*        previous xend,  no of  successful  steps,  no  of  successive * */
/*        failures, no of fcn evals, and the components of y            * */
/*     stop                                                             * */
/*  end abort action                                                    * */
/*                                                                      * */
/* *********************************************************************** */

/*     ****************************************************************** */
/*     * begin initialization, parameter checking, interrupt re-entries * */
/*     ****************************************************************** */

/*  ......abort if ind out of range 1 to 6 */
    /* Parameter adjustments */
    --y;
    --c__;
    w_dim1 = *nw;
    w_offset = 1 + w_dim1;
    w -= w_offset;

    /* Function Body */
    if (*ind < 1 || *ind > 6) {
	goto L500;
    }

/*        cases - initial entry, normal re-entry, interrupt re-entries */
    switch (*ind) {
	case 1:  goto L5;
	case 2:  goto L5;
	case 3:  goto L45;
	case 4:  goto L1111;
	case 5:  goto L2222;
	case 6:  goto L2222;
    }
/*        case 1 - initial entry (ind .eq. 1 or 2) */
/*  .........abort if n.gt.nw or tol.le.0 */
L5:
    if (*n > *nw || *tol <= 0.) {
	goto L500;
    }
    if (*ind == 2) {
	goto L15;
    }
/*              initial entry without options (ind .eq. 1) */
/*              set c(1) to c(9) equal to 0 */
    for (k = 1; k <= 9; ++k) {
	c__[k] = 0.;
/* L10: */
    }
    goto L35;
L15:
/*              initial entry with options (ind .eq. 2) */
/*              make c(1) to c(9) non-negative */
    for (k = 1; k <= 9; ++k) {
	c__[k] = (d__1 = c__[k], abs(d__1));
/* L20: */
    }
/*              make floor values non-negative if they are to be used */
    if (c__[1] != 4. && c__[1] != 5.) {
	goto L30;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	c__[k + 30] = (d__1 = c__[k + 30], abs(d__1));
/* L25: */
    }
L30:
L35:
/*           initialize rreb, dwarf, prev xend, flag, counts */
    c__[10] = 1.3877787807814457e-17;
    c__[11] = 1e-35;
/*           set previous xend initially to initial value of x */
    c__[20] = *x;
    for (k = 21; k <= 24; ++k) {
	c__[k] = 0.;
/* L40: */
    }
    goto L50;
/*        case 2 - normal re-entry (ind .eq. 3) */
/*  .........abort if xend reached, and either x changed or xend not */
L45:
    if (c__[21] != 0. && (*x != c__[20] || *xend == c__[20])) {
	goto L500;
    }
/*           re-initialize flag */
    c__[21] = 0.;
    goto L50;
/*        case 3 - re-entry following an interrupt (ind .eq. 4 to 6) */
/*           transfer control to the appropriate re-entry point.......... */
/*           this has already been handled by the computed go to        . */
/*        end cases                                                     v */
L50:

/*     end initialization, etc. */

/*     ****************************************************************** */
/*     * loop through the following 4 stages, once for each trial  step * */
/*     * until the occurrence of one of the following                   * */
/*     *    (a) the normal return (with ind .eq. 3) on reaching xend in * */
/*     *        stage 4                                                 * */
/*     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 * */
/*     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if * */
/*     *        requested, in stage 1 or stage 4                        * */
/*     ****************************************************************** */

L99999:

/*        *************************************************************** */
/*        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., * */
/*        * and some parameter  checking,  and  end  up  with  suitable * */
/*        * values of hmag, xtrial and htrial in preparation for taking * */
/*        * an integration step.                                        * */
/*        *************************************************************** */

/* ***********error return (with ind=-1) if no of fcn evals too great */
    if (c__[7] == 0. || c__[24] < c__[7]) {
	goto L100;
    }
    *ind = -1;
    return 0;
L100:

/*           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6 */
    if (*ind == 6) {
	goto L105;
    }
    (*fcn)(n, x, &y[1], &w[w_dim1 + 1]);
    c__[24] += 1.;
L105:

/*           calculate hmin - use default unless value prescribed */
    c__[13] = c__[3];
    if (c__[3] != 0.) {
	goto L165;
    }
/*              calculate default value of hmin */
/*              first calculate weighted norm y - c(12) - as specified */
/*              by the error control indicator c(1) */
    temp = 0.;
    if (c__[1] != 1.) {
	goto L115;
    }
/*                 absolute error control - weights are 1 */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = temp, d__3 = (d__1 = y[k], abs(d__1));
	temp = max(d__2,d__3);
/* L110: */
    }
    c__[12] = temp;
    goto L160;
L115:
    if (c__[1] != 2.) {
	goto L120;
    }
/*                 relative error control - weights are 1/dabs(y(k)) so */
/*                 weighted norm y is 1 */
    c__[12] = 1.;
    goto L160;
L120:
    if (c__[1] != 3.) {
	goto L130;
    }
/*                 weights are 1/max(c(2),abs(y(k))) */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = temp, d__3 = (d__1 = y[k], abs(d__1)) / c__[2];
	temp = max(d__2,d__3);
/* L125: */
    }
    c__[12] = min(temp,1.);
    goto L160;
L130:
    if (c__[1] != 4.) {
	goto L140;
    }
/*                 weights are 1/max(c(k+30),abs(y(k))) */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = temp, d__3 = (d__1 = y[k], abs(d__1)) / c__[k + 30];
	temp = max(d__2,d__3);
/* L135: */
    }
    c__[12] = min(temp,1.);
    goto L160;
L140:
    if (c__[1] != 5.) {
	goto L150;
    }
/*                 weights are 1/c(k+30) */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = temp, d__3 = (d__1 = y[k], abs(d__1)) / c__[k + 30];
	temp = max(d__2,d__3);
/* L145: */
    }
    c__[12] = temp;
    goto L160;
L150:
/*                 default case - weights are 1/max(1,abs(y(k))) */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = temp, d__3 = (d__1 = y[k], abs(d__1));
	temp = max(d__2,d__3);
/* L155: */
    }
    c__[12] = min(temp,1.);
L160:
/* Computing MAX */
/* Computing MAX */
    d__3 = c__[12] / *tol, d__4 = abs(*x);
    d__1 = c__[11], d__2 = c__[10] * max(d__3,d__4);
    c__[13] = max(d__1,d__2) * 10.;
L165:

/*           calculate scale - use default unless value prescribed */
    c__[15] = c__[5];
    if (c__[5] == 0.) {
	c__[15] = 1.;
    }

/*           calculate hmax - consider 4 cases */
/*           case 1 both hmax and scale prescribed */
    if (c__[6] != 0. && c__[5] != 0.) {
/* Computing MIN */
	d__1 = c__[6], d__2 = 2. / c__[5];
	c__[16] = min(d__1,d__2);
    }
/*           case 2 - hmax prescribed, but scale not */
    if (c__[6] != 0. && c__[5] == 0.) {
	c__[16] = c__[6];
    }
/*           case 3 - hmax not prescribed, but scale is */
    if (c__[6] == 0. && c__[5] != 0.) {
	c__[16] = 2. / c__[5];
    }
/*           case 4 - neither hmax nor scale is provided */
    if (c__[6] == 0. && c__[5] == 0.) {
	c__[16] = 2.;
    }

/* ***********error return (with ind=-2) if hmin .gt. hmax */
    if (c__[13] <= c__[16]) {
	goto L170;
    }
    *ind = -2;
    return 0;
L170:

/*           calculate preliminary hmag - consider 3 cases */
    if (*ind > 2) {
	goto L175;
    }
/*           case 1 - initial entry - use prescribed value of hstart, if */
/*              any, else default */
    c__[14] = c__[4];
    if (c__[4] == 0.) {
	c__[14] = c__[16] * pow_dd(tol, &c_b32);
    }
    goto L185;
L175:
    if (c__[23] > 1.) {
	goto L180;
    }
/*           case 2 - after a successful step, or at most  one  failure, */
/*              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible */
/*              overflow. then avoid reduction by more than half. */
    temp = c__[14] * 2.;
    if (*tol < c__[19] * 120.42729108217097) {
	d__1 = *tol / c__[19];
	temp = pow_dd(&d__1, &c_b32) * .9 * c__[14];
    }
/* Computing MAX */
    d__1 = temp, d__2 = c__[14] * .5;
    c__[14] = max(d__1,d__2);
    goto L185;
L180:
/*           case 3 - after two or more successive failures */
    c__[14] *= .5;
L185:

/*           check against hmax */
    c__[14] = min(c__[14],c__[16]);

/*           check against hmin */
    c__[14] = max(c__[14],c__[13]);

/* ***********interrupt no 1 (with ind=4) if requested */
    if (c__[8] == 0.) {
	goto L1111;
    }
    *ind = 4;
    return 0;
/*           resume here on re-entry with ind .eq. 4   ........re-entry.. */
L1111:

/*           calculate hmag, xtrial - depending on preliminary hmag, xend */
    if (c__[14] >= (d__1 = *xend - *x, abs(d__1))) {
	goto L190;
    }
/*              do not step more than half way to xend */
/* Computing MIN */
    d__2 = c__[14], d__3 = (d__1 = *xend - *x, abs(d__1)) * .5;
    c__[14] = min(d__2,d__3);
    d__1 = *xend - *x;
    c__[17] = *x + d_sign(&c__[14], &d__1);
    goto L195;
L190:
/*              hit xend exactly */
    c__[14] = (d__1 = *xend - *x, abs(d__1));
    c__[17] = *xend;
L195:

/*           calculate htrial */
    c__[18] = c__[17] - *x;

/*        end stage 1 */

/*        *************************************************************** */
/*        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). * */
/*        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in * */
/*        * stage 3. w(*,9) is temporary storage until finally it holds * */
/*        * ytrial.                                                     * */
/*        *************************************************************** */

    temp = c__[18] / 1.39816908e12;

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w[k + w_dim1 * 9] = y[k] + temp * w[k + w_dim1] * 2.3302818e11;
/* L200: */
    }
    d__1 = *x + c__[18] / 6.;
    (*fcn)(n, &d__1, &w[w_dim1 * 9 + 1], &w[(w_dim1 << 1) + 1]);

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w[k + w_dim1 * 9] = y[k] + temp * (w[k + w_dim1] * 74569017600. + w[k 
		+ (w_dim1 << 1)] * 298276070400.);
/* L205: */
    }
    d__1 = *x + c__[18] * .26666666666666666;
    (*fcn)(n, &d__1, &w[w_dim1 * 9 + 1], &w[w_dim1 * 3 + 1]);

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w[k + w_dim1 * 9] = y[k] + temp * (w[k + w_dim1] * 1.1651409e12 - w[k 
		+ (w_dim1 << 1)] * 3.72845088e12 + w[k + w_dim1 * 3] * 
		3.4954227e12);
/* L210: */
    }
    d__1 = *x + c__[18] * .66666666666666663;
    (*fcn)(n, &d__1, &w[w_dim1 * 9 + 1], &w[(w_dim1 << 2) + 1]);

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w[k + w_dim1 * 9] = y[k] + temp * (-w[k + w_dim1] * 3604654659375. + 
		w[k + (w_dim1 << 1)] * 1.28165499e13 - w[k + w_dim1 * 3] * 
		9284716546875. + w[k + (w_dim1 << 2)] * 1237962206250.);
/* L215: */
    }
    d__1 = *x + c__[18] * .83333333333333337;
    (*fcn)(n, &d__1, &w[w_dim1 * 9 + 1], &w[w_dim1 * 5 + 1]);

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w[k + w_dim1 * 9] = y[k] + temp * (w[k + w_dim1] * 3.355605792e12 - w[
		k + (w_dim1 << 1)] * 1.118535264e13 + w[k + w_dim1 * 3] * 
		9.17262885e12 - w[k + (w_dim1 << 2)] * 4.2721833e11 + w[k + 
		w_dim1 * 5] * 4.82505408e11);
/* L220: */
    }
    d__1 = *x + c__[18];
    (*fcn)(n, &d__1, &w[w_dim1 * 9 + 1], &w[w_dim1 * 6 + 1]);

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w[k + w_dim1 * 9] = y[k] + temp * (-w[k + w_dim1] * 770204740536. + w[
		k + (w_dim1 << 1)] * 2311639545600. - w[k + w_dim1 * 3] * 
		1.322092233e12 - w[k + (w_dim1 << 2)] * 453006781920. + w[k + 
		w_dim1 * 5] * 326875481856.);
/* L225: */
    }
    d__1 = *x + c__[18] / 15.;
    (*fcn)(n, &d__1, &w[w_dim1 * 9 + 1], &w[w_dim1 * 7 + 1]);

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w[k + w_dim1 * 9] = y[k] + temp * (w[k + w_dim1] * 2.845924389e12 - w[
		k + (w_dim1 << 1)] * 9.754668e12 + w[k + w_dim1 * 3] * 
		7.897110375e12 - w[k + (w_dim1 << 2)] * 1.9208266e11 + w[k + 
		w_dim1 * 5] * 4.00298976e11 + w[k + w_dim1 * 7] * 2.01586e11);
/* L230: */
    }
    d__1 = *x + c__[18];
    (*fcn)(n, &d__1, &w[w_dim1 * 9 + 1], &w[(w_dim1 << 3) + 1]);

/*           calculate ytrial, the extrapolated approximation and store */
/*              in w(*,9) */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w[k + w_dim1 * 9] = y[k] + temp * (w[k + w_dim1] * 1.04862681e11 + w[
		k + w_dim1 * 3] * 5.4518625e11 + w[k + (w_dim1 << 2)] * 
		4.46637345e11 + w[k + w_dim1 * 5] * 1.88806464e11 + w[k + 
		w_dim1 * 7] * 1.5076875e10 + w[k + (w_dim1 << 3)] * 
		9.7599465e10);
/* L235: */
    }

/*           add 7 to the no of fcn evals */
    c__[24] += 7.;

/*        end stage 2 */

/*        *************************************************************** */
/*        * stage 3 - calculate the error estimate est. first calculate * */
/*        * the  unweighted  absolute  error  estimate vector (per unit * */
/*        * step) for the unextrapolated approximation and store it  in * */
/*        * w(*,2).  then  calculate the weighted max norm of w(*,2) as * */
/*        * specified by the error  control  indicator  c(1).  finally, * */
/*        * modify  this result to produce est, the error estimate (per * */
/*        * unit step) for the extrapolated approximation ytrial.       * */
/*        *************************************************************** */

/*           calculate the unweighted absolute error estimate vector */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	w[k + (w_dim1 << 1)] = (w[k + w_dim1] * 8738556750. + w[k + w_dim1 * 
		3] * 9735468750. - w[k + (w_dim1 << 2)] * 9709507500. + w[k + 
		w_dim1 * 5] * 8.582112e9 + w[k + w_dim1 * 6] * 9.532971e10 - 
		w[k + w_dim1 * 7] * 1.5076875e10 - w[k + (w_dim1 << 3)] * 
		9.7599465e10) / 1.39816908e12;
/* L300: */
    }

/*           calculate the weighted max norm of w(*,2) as specified by */
/*           the error control indicator c(1) */
    temp = 0.;
    if (c__[1] != 1.) {
	goto L310;
    }
/*              absolute error control */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = temp, d__3 = (d__1 = w[k + (w_dim1 << 1)], abs(d__1));
	temp = max(d__2,d__3);
/* L305: */
    }
    goto L360;
L310:
    if (c__[1] != 2.) {
	goto L320;
    }
/*              relative error control */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = temp, d__3 = (d__1 = w[k + (w_dim1 << 1)] / y[k], abs(d__1));
	temp = max(d__2,d__3);
/* L315: */
    }
    goto L360;
L320:
    if (c__[1] != 3.) {
	goto L330;
    }
/*              weights are 1/max(c(2),abs(y(k))) */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
/* Computing MAX */
	d__5 = c__[2], d__6 = (d__2 = y[k], abs(d__2));
	d__3 = temp, d__4 = (d__1 = w[k + (w_dim1 << 1)], abs(d__1)) / max(
		d__5,d__6);
	temp = max(d__3,d__4);
/* L325: */
    }
    goto L360;
L330:
    if (c__[1] != 4.) {
	goto L340;
    }
/*              weights are 1/max(c(k+30),abs(y(k))) */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
/* Computing MAX */
	d__5 = c__[k + 30], d__6 = (d__2 = y[k], abs(d__2));
	d__3 = temp, d__4 = (d__1 = w[k + (w_dim1 << 1)], abs(d__1)) / max(
		d__5,d__6);
	temp = max(d__3,d__4);
/* L335: */
    }
    goto L360;
L340:
    if (c__[1] != 5.) {
	goto L350;
    }
/*              weights are 1/c(k+30) */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = temp, d__3 = (d__1 = w[k + (w_dim1 << 1)] / c__[k + 30], abs(
		d__1));
	temp = max(d__2,d__3);
/* L345: */
    }
    goto L360;
L350:
/*              default case - weights are 1/max(1,abs(y(k))) */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
/* Computing MAX */
	d__5 = 1., d__6 = (d__2 = y[k], abs(d__2));
	d__3 = temp, d__4 = (d__1 = w[k + (w_dim1 << 1)], abs(d__1)) / max(
		d__5,d__6);
	temp = max(d__3,d__4);
/* L355: */
    }
L360:

/*           calculate est - (the weighted max norm of w(*,2))*hmag*scale */
/*              - est is intended to be a measure of the error  per  unit */
/*              step in ytrial */
    c__[19] = temp * c__[14] * c__[15];

/*        end stage 3 */

/*        *************************************************************** */
/*        * stage 4 - make decisions.                                   * */
/*        *************************************************************** */

/*           set ind=5 if step acceptable, else set ind=6 */
    *ind = 5;
    if (c__[19] > *tol) {
	*ind = 6;
    }

/* ***********interrupt no 2 if requested */
    if (c__[9] == 0.) {
	goto L2222;
    }
    return 0;
/*           resume here on re-entry with ind .eq. 5 or 6   ...re-entry.. */
L2222:

    if (*ind == 6) {
	goto L410;
    }
/*              step accepted (ind .eq. 5), so update x, y from xtrial, */
/*                 ytrial, add 1 to the no of successful steps, and set */
/*                 the no of successive failures to zero */
    *x = c__[17];
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	y[k] = w[k + w_dim1 * 9];
/* L400: */
    }
    c__[22] += 1.;
    c__[23] = 0.;
/* **************return(with ind=3, xend saved, flag set) if x .eq. xend */
    if (*x != *xend) {
	goto L405;
    }
    *ind = 3;
    c__[20] = *xend;
    c__[21] = 1.;
    return 0;
L405:
    goto L420;
L410:
/*              step not accepted (ind .eq. 6), so add 1 to the no of */
/*                 successive failures */
    c__[23] += 1.;
/* **************error return (with ind=-3) if hmag .le. hmin */
    if (c__[14] > c__[13]) {
	goto L415;
    }
    *ind = -3;
    return 0;
L415:
L420:

/*        end stage 4 */

    goto L99999;
/*     end loop */

/*  begin abort action */
L500:

    s_wsfe(&io___3);
    do_fio(&c__1, (char *)&(*ind), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*tol), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&c__[13], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*xend), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*nw), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&c__[16], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&c__[20], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&c__[22], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&c__[23], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&c__[24], (ftnlen)sizeof(doublereal));
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&y[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();

    s_stop("", (ftnlen)0);

/*  end abort action */

    return 0;
} /* dverk_ */

