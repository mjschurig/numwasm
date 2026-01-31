/* ../reference/netlib/rksuite.f -- translated by f2c (version 20240504).
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

/* Common Block Declarations */

struct {
    doublereal tstrt, tnd, dir, hstrt, tolr;
    integer neqn;
} rkcom1_;

#define rkcom1_1 rkcom1_

struct {
    doublereal t, h__, told, hold;
    integer nfcn, svnfcn, okstp, flstp;
    logical first, last;
} rkcom2_;

#define rkcom2_1 rkcom2_

struct {
    integer prthrs, prerst, prwt, pryold, prscr, pry, pryp, prstgs, printp, 
	    lnintp;
} rkcom3_;

#define rkcom3_1 rkcom3_

struct {
    doublereal a[169]	/* was [13][13] */, b[13], c__[13], bhat[13], r__[66]	
	    /* was [11][6] */, e[7];
    integer ptr[13], nstage, methd, mintp;
    logical intp;
} rkcom4_;

#define rkcom4_1 rkcom4_

struct {
    doublereal toosml, cost, safety, expon, stbrad, tanang, rs, rs1, rs2, rs3,
	     rs4;
    integer order, lststg, maxtry, nsec;
    logical fsal;
} rkcom5_;

#define rkcom5_1 rkcom5_

struct {
    doublereal maxerr, locmax;
    integer gnfcn, przstg, przy, przyp, przers, przerr, przynu;
    logical erason, erasfl;
} rkcom6_;

#define rkcom6_1 rkcom6_

struct {
    doublereal mcheps, dwarf, rndoff, sqrrmc, cubrmc, tiny;
    integer outch;
} rkcom7_;

#define rkcom7_1 rkcom7_

struct {
    logical msg, utask;
} rkcom8_;

#define rkcom8_1 rkcom8_

struct {
    char rec[800];
} rkcom9_;

#define rkcom9_1 rkcom9_

/* Table of constant values */

static logical c_false = FALSE_;
static integer c__1 = 1;
static doublereal c_b65 = 1.;
static logical c_true = TRUE_;
static integer c__5000 = 5000;
static doublereal c_b517 = .33333333333333331;
static integer c__9 = 9;

/* Subroutine */ int setup_(integer *neq, doublereal *tstart, doublereal *
	ystart, doublereal *tend, doublereal *tol, doublereal *thres, integer 
	*method, char *task, logical *errass, doublereal *hstart, doublereal *
	work, integer *lenwrk, logical *mesage, ftnlen task_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    integer l, ier, flag__, nrec;
    doublereal hmin;
    integer lreq;
    char task1[1];
    extern /* Subroutine */ int rkmsg_(integer *, char *, integer *, integer *
	    , ftnlen), const_(integer *, integer *, logical *, integer *), 
	    rksit_(logical *, char *, integer *, ftnlen);
    logical legalt;
    integer freepr, vecstg, lintpl;
    logical reqstg;
    extern /* Subroutine */ int mconst_(integer *);

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

/*  If you are not familiar with the code SETUP and how it is used in */
/*  conjunction with UT or CT to solve initial value problems, you should study */
/*  the document file rksuite.doc carefully before attempting to use the code. */
/*  The following "Brief Reminder" is intended only to remind you of the */
/*  meaning, type, and size requirements of the arguments. */

/*  The environmental parameters OUTCH, MCHEPS, and DWARF are used in the */
/*  following description.  To find out their values */

/*       CALL ENVIRN(OUTCH,MCHEPS,DWARF) */

/*  INPUT VARIABLES */

/*     NEQ       - INTEGER */
/*                 The number of differential equations in the system. */
/*                 Constraint: NEQ >= 1 */
/*     TSTART    - DOUBLE PRECISION */
/*                 The initial value of the independent variable. */
/*     YSTART(*) - DOUBLE PRECISION array of length NEQ */
/*                 The vector of initial values of the solution components. */
/*     TEND      - DOUBLE PRECISION */
/*                 The integration proceeds from TSTART in the direction of */
/*                 TEND. You cannot go past TEND. */
/*                 Constraint: TEND must be clearly distinguishable from TSTART */
/*                 in the precision available. */
/*     TOL       - DOUBLE PRECISION */
/*                 The relative error tolerance. */
/*                 Constraint: 0.01D0 >= TOL >= 10*MCHEPS */
/*     THRES(*)  - DOUBLE PRECISION array of length NEQ */
/*                 THRES(L) is the threshold for the Ith solution component. */
/*                 Constraint: THRES(L) >= SQRT(DWARF) */
/*     METHOD    - INTEGER */
/*                 Specifies which Runge-Kutta pair is to be used. */
/*                  = 1 - use the (2,3) pair */
/*                  = 2 - use the (4,5) pair */
/*                  = 3 - use the (7,8) pair */
/*     TASK      - CHARACTER*(*) */
/*                 Only the first character of TASK is significant. */
/*                 TASK(1:1) = `U' or `u' - UT is to be used */
/*                           = `C' or `c' - CT is to be used */
/*                 Constraint: TASK(1:1) = `U'or `u' or`C' or `c' */
/*     ERRASS    - LOGICAL */
/*                 = .FALSE. - do not attempt to assess the true error. */
/*                 = .TRUE.  - assess the true error. Costs roughly twice */
/*                             as much as the integration with METHODs 2 and */
/*                             3, and three times with METHOD = 1. */
/*     HSTART    - DOUBLE PRECISION */
/*                 0.0D0     - select automatically the first step size. */
/*                 non-zero  - try HSTART for the first step. */

/*  WORKSPACE */

/*     WORK(*) - DOUBLE PRECISION array of length LENWRK */
/*               Do not alter the contents of this array after calling SETUP. */

/*  INPUT VARIABLES */

/*     LENWRK  - INTEGER */
/*               Length of WORK(*): How big LENWRK must be depends */
/*               on the task and how it is to be solved. */

/*               LENWRK = 32*NEQ is sufficient for all cases. */

/*               If storage is a problem, the least storage possible */
/*               in the various cases is: */

/*                 If TASK = `U' or `u', then */
/*                   if ERRASS = .FALSE. and */
/*                     METHOD = 1, LENWRK must be at least 10*NEQ */
/*                            = 2                          20*NEQ */
/*                            = 3                          16*NEQ */
/*                   if ERRASS = .TRUE. and */
/*                     METHOD = 1, LENWRK must be at least 15*NEQ */
/*                            = 2                          32*NEQ */
/*                            = 3                          21*NEQ */

/*                 If TASK = `C' or `c', then */
/*                   if ERRASS = .FALSE. and */
/*                     METHOD = 1, LENWRK must be at least 10*NEQ */
/*                            = 2                          14*NEQ */
/*                            = 3                          16*NEQ */
/*                   if ERRASS = .TRUE. and */
/*                     METHOD = 1, LENWRK must be at least 15*NEQ */
/*                            = 2                          26*NEQ */
/*                            = 3                          21*NEQ */

/*                 Warning:  To exploit the interpolation capability */
/*                 of METHODs 1 and 2, you have to call INTRP.  This */
/*                 subroutine requires working storage in addition to */
/*                 that specified here. */

/*     MESAGE    - LOGICAL */
/*                 Specifies whether you want informative messages written to */
/*                 the standard output channel OUTCH. */
/*                 = .TRUE.   - provide messages */
/*                 = .FALSE.  - do not provide messages */

/*  In the event of a "catastrophic" failure to call SETUP correctly, the */
/*  nature of the catastrophe is reported on the standard output channel, */
/*  regardless of the value of MESAGE.  Unless special provision was made */
/*  in advance (see rksuite.doc), the computation then comes to a STOP. */

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block for General Workspace Pointers .. */
/*     .. Common Block to hold Formula Definitions .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Global Error Assessment .. */
/*     .. Common Block for Environment Parameters .. */
/*     .. Common Block for Integrator Options .. */
/*     .. Common Block for Error Message .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*  Clear previous flag values of subprograms in the suite. */

    /* Parameter adjustments */
    --work;
    --thres;
    --ystart;

    /* Function Body */
    ier = -1;
    rksit_(&c_false, "SETUP ", &ier, (ftnlen)6);

    ier = 1;
    nrec = 0;

/*  Fetch output channel and machine constants; initialise common */
/*  block /RKCOM7/ */

    mconst_(method);

/*  Check for valid input of trivial arguments */
    *(unsigned char *)task1 = *(unsigned char *)task;
    legalt = *(unsigned char *)task1 == 'U' || *(unsigned char *)task1 == 'u' 
	    || *(unsigned char *)task1 == 'C' || *(unsigned char *)task1 == 
	    'c';
    if (! legalt) {
	ier = 911;
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A,A,A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have set the first character of ", (ftnlen)40);
	do_fio(&c__1, " ** TASK to be '", (ftnlen)16);
	do_fio(&c__1, task1, (ftnlen)1);
	do_fio(&c__1, "'. It must be one of ", (ftnlen)21);
	do_fio(&c__1, " ** 'U','u','C' or 'c'.", (ftnlen)23);
	e_wsfi();
    } else if (*neq < 1) {
	ier = 911;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A,I6,A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have set NEQ = ", (ftnlen)23);
	do_fio(&c__1, (char *)&(*neq), (ftnlen)sizeof(integer));
	do_fio(&c__1, " which is less than 1.", (ftnlen)22);
	e_wsfi();
    } else if (*method < 1 || *method > 3) {
	ier = 911;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A,I6,A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have set METHOD = ", (ftnlen)26);
	do_fio(&c__1, (char *)&(*method), (ftnlen)sizeof(integer));
	do_fio(&c__1, " which is not 1, 2, or 3.", (ftnlen)25);
	e_wsfi();
    } else if (*tstart == *tend) {
	ier = 911;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A,D13.5,A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have set TSTART = TEND = ", (ftnlen)33);
	do_fio(&c__1, (char *)&(*tstart), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, ".", (ftnlen)1);
	e_wsfi();
    } else if (*tol > .01 || *tol < rkcom7_1.rndoff) {
	ier = 911;
	nrec = 2;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A,D13.5,A/A,D13.5,A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have set TOL = ", (ftnlen)23);
	do_fio(&c__1, (char *)&(*tol), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " which is not permitted. The", (ftnlen)28);
	do_fio(&c__1, " ** range of permitted values is (", (ftnlen)34);
	do_fio(&c__1, (char *)&rkcom7_1.rndoff, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, ",0.01D0).", (ftnlen)9);
	e_wsfi();
    } else {
	l = 1;
L20:
	if (thres[l] < rkcom7_1.tiny) {
	    ier = 911;
	    nrec = 2;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A,I6,A,D13.5,A/A,D13.5,A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have set THRES(", (ftnlen)23);
	    do_fio(&c__1, (char *)&l, (ftnlen)sizeof(integer));
	    do_fio(&c__1, ") to be ", (ftnlen)8);
	    do_fio(&c__1, (char *)&thres[l], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, " which is ", (ftnlen)10);
	    do_fio(&c__1, " ** less than the permitted minimum,", (ftnlen)36);
	    do_fio(&c__1, (char *)&rkcom7_1.tiny, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ".", (ftnlen)1);
	    e_wsfi();
	}
	++l;
	if (ier != 911 && l <= *neq) {
	    goto L20;
	}
    }

/*  Return if error detected */

    if (ier != 1) {
	goto L80;
    }

/*  Set formula definitions and characteristics by means of arguments */
/*  in the call list and COMMON blocks /RKCOM4/ and /RKCOM5/ */

    const_(method, &vecstg, &reqstg, &lintpl);

/*  Set options in /RKCOM8/ */
    rkcom8_1.utask = *(unsigned char *)task1 == 'U' || *(unsigned char *)
	    task1 == 'u';
    rkcom8_1.msg = *mesage;

/*  Initialise problem status in /RKCOM1/ and /RKCOM2/ */
    rkcom1_1.neqn = *neq;
    rkcom1_1.tstrt = *tstart;
    rkcom1_1.tnd = *tend;
    rkcom2_1.t = *tstart;
    rkcom2_1.told = *tstart;
    d__1 = *tend - *tstart;
    rkcom1_1.dir = d_sign(&c_b65, &d__1);

/*  In CT the first step taken will have magnitude H.  If HSTRT = ABS(HSTART) */
/*  is not equal to zero, H = HSTRT.  If HSTRT is equal to zero, the code is */
/*  to find an on-scale initial step size H.  To start this process, H is set */
/*  here to an upper bound on the first step size that reflects the scale of */
/*  the independent variable.  UT has some additional information, namely the */
/*  first output point, that is used to refine this bound in UT when UTASK */
/*  is .TRUE..  If HSTRT is not zero, but it is either too big or too small, */
/*  the input HSTART is ignored and HSTRT is set to zero to activate the */
/*  automatic determination of an on-scale initial step size. */

    rkcom1_1.hstrt = abs(*hstart);
/* Computing MAX */
/* Computing MAX */
    d__3 = abs(*tstart), d__4 = abs(*tend);
    d__1 = rkcom7_1.tiny, d__2 = rkcom5_1.toosml * max(d__3,d__4);
    hmin = max(d__1,d__2);
    if (rkcom1_1.hstrt > (d__1 = *tend - *tstart, abs(d__1)) || 
	    rkcom1_1.hstrt < hmin) {
	rkcom1_1.hstrt = 0.;
    }
    if (rkcom1_1.hstrt == 0.) {
/* Computing MAX */
	d__2 = (d__1 = *tend - *tstart, abs(d__1)) / rkcom5_1.rs3;
	rkcom2_1.h__ = max(d__2,hmin);
    } else {
	rkcom2_1.h__ = rkcom1_1.hstrt;
    }
    rkcom2_1.hold = 0.;
    rkcom1_1.tolr = *tol;
    rkcom2_1.nfcn = 0;
    rkcom2_1.svnfcn = 0;
    rkcom2_1.okstp = 0;
    rkcom2_1.flstp = 0;
    rkcom2_1.first = TRUE_;
    rkcom2_1.last = FALSE_;

/*  WORK(*) is partioned into a number of arrays using pointers. These */
/*  pointers are set in /RKCOM3/. */
    rkcom3_1.prthrs = 1;
/*                           the threshold values */
    rkcom3_1.prerst = rkcom3_1.prthrs + *neq;
/*                           the error estimates */
    rkcom3_1.prwt = rkcom3_1.prerst + *neq;
/*                           the weights used in the local error test */
    rkcom3_1.pryold = rkcom3_1.prwt + *neq;
/*                           the previous value of the solution */
    rkcom3_1.prscr = rkcom3_1.pryold + *neq;
/*                           scratch array used for the higher order */
/*                           approximate solution and for the previous */
/*                           value of the derivative of the solution */
    rkcom3_1.pry = rkcom3_1.prscr + *neq;
/*                           the dependent variables */
    rkcom3_1.pryp = rkcom3_1.pry + *neq;
/*                           the derivatives */
    rkcom3_1.prstgs = rkcom3_1.pryp + *neq;
/*                           intermediate stages held in an internal */
/*                           array STAGES(NEQ,VECSTG) */

    freepr = rkcom3_1.prstgs + vecstg * *neq;

/*  Allocate storage for interpolation if the TASK = `U' or `u' was */
/*  specified. INTP and LINTPL returned by CONST indicate whether there */
/*  is an interpolation scheme associated with the pair and how much */
/*  storage is required. */

    rkcom3_1.printp = 1;
    rkcom3_1.lnintp = 1;
    if (rkcom8_1.utask) {
	if (rkcom4_1.intp) {
	    rkcom3_1.lnintp = lintpl * *neq;
	    if (reqstg) {
		rkcom3_1.printp = freepr;
		freepr = rkcom3_1.printp + rkcom3_1.lnintp;
	    } else {
		rkcom3_1.printp = rkcom3_1.prstgs;
/* Computing MAX */
		i__1 = rkcom3_1.printp + vecstg * *neq, i__2 = 
			rkcom3_1.printp + rkcom3_1.lnintp;
		freepr = max(i__1,i__2);
	    }
	}
    }

/*  Initialise state and allocate storage for global error assessment */
/*  using /RKCOM6/ */
    rkcom6_1.gnfcn = 0;
    rkcom6_1.maxerr = 0.;
    rkcom6_1.locmax = *tstart;
    rkcom6_1.erason = *errass;
    rkcom6_1.erasfl = FALSE_;
    if (*errass) {

/*  Storage is required for the stages of a secondary integration. The */
/*  stages of the primary intergration can only be overwritten in the */
/*  cases where there is no interpolant or the interpolant does not */
/*  require information about the stages (e.g. METHOD 3 and METHOD 1, */
/*  respectively). */
	if (! reqstg) {
	    rkcom6_1.przstg = rkcom3_1.prstgs;
	} else {
	    rkcom6_1.przstg = freepr;
	    freepr = rkcom6_1.przstg + vecstg * *neq;
	}
	rkcom6_1.przy = freepr;
	rkcom6_1.przyp = rkcom6_1.przy + *neq;
	rkcom6_1.przers = rkcom6_1.przyp + *neq;
	rkcom6_1.przerr = rkcom6_1.przers + *neq;
	rkcom6_1.przynu = rkcom6_1.przerr + *neq;
	freepr = rkcom6_1.przynu + *neq;
    } else {
	rkcom6_1.przstg = 1;
	rkcom6_1.przy = 1;
	rkcom6_1.przyp = 1;
	rkcom6_1.przers = 1;
	rkcom6_1.przerr = 1;
	rkcom6_1.przynu = 1;
    }

    lreq = freepr - 1;

/*  Check for enough workspace and suitable range of integration */

    if (*lenwrk < lreq) {
	ier = 911;
	nrec = 2;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A,I6,A,I6,A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have not supplied enough workspace. You gave "
		"LENWRK ", (ftnlen)60);
	do_fio(&c__1, " ** as", (ftnlen)6);
	do_fio(&c__1, (char *)&(*lenwrk), (ftnlen)sizeof(integer));
	do_fio(&c__1, ", but it must be at least ", (ftnlen)26);
	do_fio(&c__1, (char *)&lreq, (ftnlen)sizeof(integer));
	do_fio(&c__1, ".", (ftnlen)1);
	e_wsfi();
    } else {
/* Computing MAX */
/* Computing MAX */
	d__3 = abs(*tstart), d__4 = abs(*tend);
	d__1 = rkcom7_1.tiny, d__2 = rkcom5_1.toosml * max(d__3,d__4);
	hmin = max(d__1,d__2);
	if ((d__1 = *tend - *tstart, abs(d__1)) < hmin) {
	    ier = 911;
	    nrec = 4;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A/A/A,D13.5/A,D13.5,A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have set values for TEND and TSTART that "
		    "are not ", (ftnlen)57);
	    do_fio(&c__1, " ** clearly distinguishable for the method and th"
		    "e precision ", (ftnlen)61);
	    do_fio(&c__1, " ** of the computer being used. ABS(TEND-TSTART) "
		    "is ", (ftnlen)52);
	    d__2 = (d__1 = *tend - *tstart, abs(d__1));
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, " ** but should be at least ", (ftnlen)27);
	    do_fio(&c__1, (char *)&hmin, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ".", (ftnlen)1);
	    e_wsfi();
	}
    }

/*  Return if error detected */

    if (ier != 1) {
	goto L80;
    }

/*  Initialize elements of the workspace */
    i__1 = *neq;
    for (l = 1; l <= i__1; ++l) {
	work[rkcom3_1.prthrs - 1 + l] = thres[l];
	work[rkcom3_1.pry - 1 + l] = ystart[l];
/* L40: */
    }

/*  Initialize the global error to zero when ERRASS = .TRUE. */
    if (*errass) {
	i__1 = *neq;
	for (l = 1; l <= i__1; ++l) {
	    work[rkcom6_1.przerr - 1 + l] = 0.;
/* L60: */
	}
    }

L80:

    rkmsg_(&ier, "SETUP ", &nrec, &flag__, (ftnlen)6);

    return 0;
} /* setup_ */

/* Subroutine */ int ut_(U_fp f, doublereal *twant, doublereal *tgot, 
	doublereal *ygot, doublereal *ypgot, doublereal *ymax, doublereal *
	work, integer *uflag)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    integer l;
    extern /* Subroutine */ int ct_(U_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    integer ier, nrec;
    doublereal hmin, tnow;
    integer cflag;
    extern /* Subroutine */ int chkfl_(logical *, logical *);
    static doublereal utend;
    integer state;
    extern /* Subroutine */ int reset_(doublereal *), rkmsg_(integer *, char *
	    , integer *, integer *, ftnlen);
    static doublereal tlast;
    extern /* Subroutine */ int intrp_(doublereal *, char *, integer *, 
	    doublereal *, doublereal *, U_fp, doublereal *, doublereal *, 
	    integer *, ftnlen), rksit_(logical *, char *, integer *, ftnlen);
    logical goback, baderr;

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

/*  If you are not familiar with the code UT and how it is used in */
/*  conjunction with SETUP to solve initial value problems, you should study */
/*  the document file rksuite.doc carefully before proceeding further.  The */
/*  following "Brief Reminder" is intended only to remind you of the meaning, */
/*  type, and size requirements of the arguments. */

/*  NAME DECLARED IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM: */

/*     F         - name of the subroutine for evaluating the differential */
/*                 equations. */

/*  The subroutine F must have the form */

/*  SUBROUTINE F(T,Y,YP) */
/*  DOUBLE PRECISION T,Y(*),YP(*) */
/*     Given input values of the independent variable T and the solution */
/*     components Y(*), for each L = 1,2,...,NEQ evaluate the differential */
/*     equation for the derivative of the Ith solution component and place the */
/*     value in YP(L).  Do not alter the input values of T and Y(*). */
/*  RETURN */
/*  END */

/*  INPUT VARIABLE */

/*     TWANT     - DOUBLE PRECISION */
/*                 The next value of the independent variable where a */
/*                 solution is desired. */

/*                 Constraints: TWANT must lie between the previous value */
/*                 of TGOT (TSTART on the first call) and TEND. TWANT can be */
/*                 equal to TEND, but it must be clearly distinguishable from */
/*                 the previous value of TGOT (TSTART on the first call) in */
/*                 the precision available. */

/*  OUTPUT VARIABLES */

/*     TGOT      - DOUBLE PRECISION */
/*                 A solution has been computed at this value of the */
/*                 independent variable. */
/*     YGOT(*)   - DOUBLE PRECISION array of length NEQ */
/*                 Approximation to the true solution at TGOT. Do not alter */
/*                 the contents of this array */
/*     YPGOT(*)  - DOUBLE PRECISION array of length NEQ */
/*                 Approximation to the first derivative of the true */
/*                 solution at TGOT. */
/*     YMAX(*)   - DOUBLE PRECISION array of length NEQ */
/*                 YMAX(L) is the largest magnitude of YGOT(L) computed at any */
/*                 time in the integration from TSTART to TGOT. Do not alter */
/*                 the contents of this array. */

/*  WORKSPACE */

/*     WORK(*)   - DOUBLE PRECISION array as used in SETUP */
/*                 Do not alter the contents of this array. */

/*  OUTPUT VARIABLE */

/*     UFLAG     - INTEGER */

/*                       SUCCESS.  TGOT = TWANT. */
/*                 = 1 - Complete success. */

/*                       "SOFT" FAILURES */
/*                 = 2 - Warning:  You are using METHOD = 3 inefficiently */
/*                       by computing answers at many values of TWANT.  If */
/*                       you really need answers at so many specific points, */
/*                       it would be more efficient to compute them with */
/*                       METHOD = 2.  To do this you would need to restart */
/*                       from TGOT, YGOT(*) by a call to SETUP.  If you wish */
/*                       to continue as you are, you may. */
/*                 = 3 - Warning:  A considerable amount of work has been */
/*                       expended.  If you wish to continue on to TWANT, just */
/*                       call UT again. */
/*                 = 4 - Warning:  It appears that this problem is "stiff". */
/*                       You really should change to another code that is */
/*                       intended for such problems, but if you insist, you can */
/*                       continue with UT by calling it again. */

/*                       "HARD" FAILURES */
/*                 = 5 - You are asking for too much accuracy. You cannot */
/*                       continue integrating this problem. */
/*                 = 6 - The global error assessment may not be reliable beyond */
/*                       the current point in the integration.  You cannot */
/*                       continue integrating this problem. */

/*                       "CATASTROPHIC" FAILURES */
/*                 = 911 - The nature of the catastrophe is reported on */
/*                         the standard output channel. Unless special */
/*                         provision was made in advance (see rksuite.doc), */
/*                         the computation then comes to a STOP. */

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block for General Workspace Pointers .. */
/*     .. Common Block to hold Formula Definitions .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Environment Parameters .. */
/*     .. Common Block for Integrator Options .. */
/*     .. Common Block for Error Message .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Save statement .. */
/*     .. Executable Statements .. */
    /* Parameter adjustments */
    --work;
    --ymax;
    --ypgot;
    --ygot;

    /* Function Body */
    ier = 1;
    nrec = 0;
    goback = FALSE_;
    baderr = FALSE_;

/*  Is it permissible to call UT? */

    rksit_(&c_true, "SETUP", &state, (ftnlen)5);
    if (state == 911) {
	ier = 912;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** A catastrophic error has already been detected el"
		"sewhere.", (ftnlen)61);
	e_wsfi();
	goto L100;
    }
    if (state == -1) {
	ier = 911;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have not called SETUP, so you cannot use UT.", 
		(ftnlen)52);
	e_wsfi();
	goto L100;
    }
    if (! rkcom8_1.utask) {
	ier = 911;
	nrec = 2;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have called UT after you specified in SETUP t"
		"hat ", (ftnlen)57);
	do_fio(&c__1, " ** you were going to use CT. This is not permitted.", 
		(ftnlen)52);
	e_wsfi();
	goto L100;
    }
    rksit_(&c_true, "UT    ", &state, (ftnlen)6);
    if (state == 5 || state == 6) {
	ier = 911;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** This routine has already returned with a hard fai"
		"lure.", (ftnlen)58);
	do_fio(&c__1, " ** You must call SETUP to start another problem.", (
		ftnlen)49);
	e_wsfi();
	goto L100;
    }
    state = -2;
    rksit_(&c_false, "UT    ", &state, (ftnlen)6);

    if (rkcom2_1.first) {

/*  First call. */

/*  A value of TND is specified in SETUP. When INTP = .FALSE., as with */
/*  METHD = 3, output is obtained at the specified TWANT by resetting TND */
/*  to TWANT.  At this point, before the integration gets started, this can */
/*  be done with a simple assignment.  Later it is done with a call to RESET. */
/*  The original TND is SAVEd as a local variable UTEND. */

	utend = rkcom1_1.tnd;
	if (! rkcom4_1.intp) {
	    rkcom1_1.tnd = *twant;
	}

/*  The last TGOT returned is SAVEd in the variable TLAST.  T (a variable */
/*  passed through the common block RKCOM2) records how far the integration */
/*  has advanced towards the specified TND.  When output is obtained by */
/*  interpolation, the integration goes past the TGOT returned (T is closer */
/*  to the specified TND than TGOT).  Initialize these variables and YMAX(*). */
	tlast = rkcom1_1.tstrt;
	*tgot = rkcom1_1.tstrt;
	i__1 = rkcom1_1.neqn;
	for (l = 1; l <= i__1; ++l) {
	    ymax[l] = (d__1 = work[rkcom3_1.pry - 1 + l], abs(d__1));
/* L20: */
	}

/*  If the code is to find an on-scale initial step size H, a bound was placed */
/*  on H in SETUP.  Here the first output point is used to refine this bound. */
	if (rkcom1_1.hstrt == 0.) {
/* Computing MIN */
	    d__2 = abs(rkcom2_1.h__), d__3 = (d__1 = *twant - rkcom1_1.tstrt, 
		    abs(d__1));
	    rkcom2_1.h__ = min(d__2,d__3);
/* Computing MAX */
/* Computing MAX */
	    d__3 = abs(rkcom1_1.tstrt), d__4 = abs(rkcom1_1.tnd);
	    d__1 = rkcom7_1.tiny, d__2 = rkcom5_1.toosml * max(d__3,d__4);
	    hmin = max(d__1,d__2);
	    rkcom2_1.h__ = max(rkcom2_1.h__,hmin);
	}

    } else {

/*  Subsequent call. */

	if (tlast == utend) {
	    ier = 911;
	    nrec = 3;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A/A/A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have called UT after reaching TEND. (Your"
		    " last    ", (ftnlen)58);
	    do_fio(&c__1, " ** call to UT resulted in TGOT = TEND.)  To star"
		    "t a new ", (ftnlen)57);
	    do_fio(&c__1, " ** problem, you will need to call SETUP.", (
		    ftnlen)41);
	    e_wsfi();
	    goto L100;
	}

    }

/*  Check for valid TWANT. */

    if (rkcom1_1.dir * (*twant - tlast) <= 0.) {
	ier = 911;
	nrec = 4;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A/A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have made a call to UT with a TWANT that does"
		"   ", (ftnlen)56);
	do_fio(&c__1, " ** not lie between the previous value of TGOT (TSTAR"
		"T  ", (ftnlen)56);
	do_fio(&c__1, " ** on the first call) and TEND. This is not permitte"
		"d. ", (ftnlen)56);
	do_fio(&c__1, " ** Check your program carefully.", (ftnlen)33);
	e_wsfi();
	goto L100;
    }
    if (rkcom1_1.dir * (*twant - utend) > 0.) {
/* Computing MAX */
/* Computing MAX */
	d__3 = abs(*twant), d__4 = abs(utend);
	d__1 = rkcom7_1.tiny, d__2 = rkcom5_1.toosml * max(d__3,d__4);
	hmin = max(d__1,d__2);
	if ((d__1 = *twant - utend, abs(d__1)) < hmin) {
	    ier = 911;
	    nrec = 5;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A/A/A/A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have made a call to UT with a TWANT that "
		    "does      ", (ftnlen)59);
	    do_fio(&c__1, " ** not lie between the previous value of TGOT (T"
		    "START on  ", (ftnlen)59);
	    do_fio(&c__1, " ** the first call) and TEND. This is not permitt"
		    "ed. TWANT ", (ftnlen)59);
	    do_fio(&c__1, " ** is very close to TEND, so you may have meant "
		    "to set    ", (ftnlen)59);
	    do_fio(&c__1, " ** it to be TEND exactly.  Check your program ca"
		    "refully.  ", (ftnlen)59);
	    e_wsfi();
	} else {
	    ier = 911;
	    nrec = 4;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A/A/A/A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have made a call to UT with a TWANT that "
		    "does   ", (ftnlen)56);
	    do_fio(&c__1, " ** not lie between the previous value of TGOT (T"
		    "START  ", (ftnlen)56);
	    do_fio(&c__1, " ** on the first call) and TEND. This is not perm"
		    "itted. ", (ftnlen)56);
	    do_fio(&c__1, " ** Check your program carefully.", (ftnlen)33);
	    e_wsfi();
	}
	goto L100;
    }
    if (! rkcom4_1.intp) {
/* Computing MAX */
/* Computing MAX */
	d__3 = abs(tlast), d__4 = abs(*twant);
	d__1 = rkcom7_1.tiny, d__2 = rkcom5_1.toosml * max(d__3,d__4);
	hmin = max(d__1,d__2);
	if ((d__1 = *twant - tlast, abs(d__1)) < hmin) {
	    ier = 911;
	    nrec = 4;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A/A/A/A,D13.5,A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have made a call to UT with a TWANT that "
		    "is not ", (ftnlen)56);
	    do_fio(&c__1, " ** sufficiently different from the last value of"
		    " TGOT  ", (ftnlen)56);
	    do_fio(&c__1, " ** (TSTART on the first call).  When using METHO"
		    "D = 3, ", (ftnlen)56);
	    do_fio(&c__1, " ** it must differ by at least ", (ftnlen)31);
	    do_fio(&c__1, (char *)&hmin, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ".", (ftnlen)1);
	    e_wsfi();
	    goto L100;
	}

/*  We have a valid TWANT. There is no interpolation with this METHD and */
/*  therefore we step to TWANT exactly by resetting TND with a call to RESET. */
/*  On the first step this matter is handled differently as explained above. */

	if (! rkcom2_1.first) {
	    reset_(twant);
	    chkfl_(&c_true, &baderr);
	    if (baderr) {
		goto L100;
	    }
	}
    }

/*  Process output, decide whether to take another step. */

L40:

    if (rkcom4_1.intp) {

/*  Interpolation is possible with this METHD.  The integration has */
/*  already reached T. If this is past TWANT, GOBACK is set .TRUE. and */
/*  the answers are obtained by interpolation. */

	goback = rkcom1_1.dir * (rkcom2_1.t - *twant) >= 0.;
	if (goback) {
	    intrp_(twant, "Both solution and derivative", &rkcom1_1.neqn, &
		    ygot[1], &ypgot[1], (U_fp)f, &work[1], &work[
		    rkcom3_1.printp], &rkcom3_1.lnintp, (ftnlen)28);
	    chkfl_(&c_true, &baderr);
	    if (baderr) {
		goto L100;
	    }
	    *tgot = *twant;
	}
    } else {

/*  Interpolation is not possible with this METHD, so output is obtained */
/*  by integrating to TWANT = TND.  Both YGOT(*) and YPGOT(*) are then */
/*  already loaded with the solution at TWANT by CT. */

	goback = rkcom2_1.t == *twant;
	if (goback) {
	    *tgot = *twant;
	}
    }

/*  Updating of YMAX(*) is done here to account for the fact that when */
/*  interpolation is done, the integration goes past TGOT.  Note that YGOT(*) */
/*  is not defined until CT is called.  YMAX(*) was initialized at TSTRT */
/*  from values stored in WORK(*), so only needs to be updated for T */
/*  different from TSTRT. */
    if (rkcom2_1.t != rkcom1_1.tstrt) {
	i__1 = rkcom1_1.neqn;
	for (l = 1; l <= i__1; ++l) {
/* Computing MAX */
	    d__2 = ymax[l], d__3 = (d__1 = ygot[l], abs(d__1));
	    ymax[l] = max(d__2,d__3);
/* L60: */
	}
    }

/*  If done, go to the exit point. */
    if (goback) {
	goto L100;
    }

/*  Take a step with CT in the direction of TND.  On exit, the solution is */
/*  advanced to TNOW.  The way CT is written, the approximate solution at */
/*  TNOW is available in both YGOT(*) and in WORK(*).  If output is obtained by */
/*  stepping to the end (TNOW = TWANT = TND), YGOT(*) can be returned directly. */
/*  If output is obtained by interpolation, the subroutine INTRP that does this */
/*  uses the values in WORK(*) for its computations and places the approximate */
/*  solution at TWANT in the array YGOT(*) for return to the calling program. */
/*  The approximate derivative is handled in the same way. TNOW is output from */
/*  CT and is actually a copy of T declared above in a common block. */

    ct_((U_fp)f, &tnow, &ygot[1], &ypgot[1], &work[1], &cflag);
    ier = cflag;

/*  A successful step by CT is indicated by CFLAG = 1 or = 2. */
    if (cflag == 1) {
	goto L40;
    } else if (cflag == 2) {

/*  Supplement the warning message written in CT. */
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** The last message was produced on a call to CT fro"
		"m UT.  ", (ftnlen)60);
	do_fio(&c__1, " ** In UT the appropriate action is to change to METH"
		"OD = 2,", (ftnlen)60);
	do_fio(&c__1, " ** or, if insufficient memory is available, to METHO"
		"D = 1. ", (ftnlen)60);
	e_wsfi();
    } else if (cflag <= 6) {
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** The last message was produced on a call to CT fro"
		"m UT.", (ftnlen)58);
	e_wsfi();
    } else {
	baderr = TRUE_;
    }
    *tgot = rkcom2_1.t;

/*  Update YMAX(*) before the return. */
    i__1 = rkcom1_1.neqn;
    for (l = 1; l <= i__1; ++l) {
/* Computing MAX */
	d__2 = ymax[l], d__3 = (d__1 = ygot[l], abs(d__1));
	ymax[l] = max(d__2,d__3);
/* L80: */
    }

/*  Exit point for UT. */

L100:

    if (baderr) {
	ier = 911;
	nrec = 4;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A/A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** An internal call by UT to a subroutine resulted i"
		"n an  ", (ftnlen)59);
	do_fio(&c__1, " ** error that should not happen.  Check your program"
		"      ", (ftnlen)59);
	do_fio(&c__1, " ** carefully for array sizes, correct number of argu"
		"ments,", (ftnlen)59);
	do_fio(&c__1, " ** type mismatches ... .", (ftnlen)25);
	e_wsfi();
    }

    tlast = *tgot;

/*  All exits are done here after a call to RKMSG to report */
/*  what happened and set UFLAG. */

    rkmsg_(&ier, "UT    ", &nrec, uflag, (ftnlen)6);

    return 0;
} /* ut_ */

/* Subroutine */ int stat_(integer *totfcn, integer *stpcst, doublereal *
	waste, integer *stpsok, doublereal *hnext)
{
    /* System generated locals */
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    integer ier, flag__, nrec, state;
    extern /* Subroutine */ int rkmsg_(integer *, char *, integer *, integer *
	    , ftnlen), rksit_(logical *, char *, integer *, ftnlen);

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

/*  If you are not familiar with the code STAT and how it is used in */
/*  conjunction with the integrators CT and UT, you should study the */
/*  document file rksuite.doc carefully before attempting to use the code. */
/*  The following "Brief Reminder" is intended only to remind you of the */
/*  meaning, type, and size requirements of the arguments. */

/*  STAT is called to obtain some details about the integration. */

/*  OUTPUT VARIABLES */

/*     TOTFCN    - INTEGER */
/*                 Total number of calls to F in the integration so far -- */
/*                 a measure of the cost of the integration. */
/*     STPCST    - INTEGER */
/*                 Cost of a typical step with this METHOD measured in */
/*                 calls to F. */
/*     WASTE     - DOUBLE PRECISION */
/*                 The number of attempted steps that failed to meet the */
/*                 local error requirement divided by the total number of */
/*                 steps attempted so far in the integration. */
/*     STPSOK    - INTEGER */
/*                 The number of accepted steps. */
/*     HNEXT     - DOUBLE PRECISION */
/*                 The step size the integrator plans to use for the next step. */

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
/*     .. Scalar Arguments .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Global Error Assessment .. */
/*     .. Common Block for Integrator Options .. */
/*     .. Common Block for Error Message .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    ier = 1;
    nrec = 0;

/*  Is it permissible to call STAT? */

    rksit_(&c_true, "STAT  ", &state, (ftnlen)6);
    if (state == 911) {
	ier = 912;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** A catastrophic error has already been detected el"
		"sewhere.", (ftnlen)61);
	e_wsfi();
	goto L20;
    }
    if (state == -2) {
	ier = 911;
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have already made a call to STAT after a hard"
		"   ", (ftnlen)56);
	do_fio(&c__1, " ** failure was reported from the integrator. You can"
		"not", (ftnlen)56);
	do_fio(&c__1, " ** call STAT again.", (ftnlen)20);
	e_wsfi();
	goto L20;
    }
    rksit_(&c_true, "CT", &state, (ftnlen)2);
    if (state == -1) {
	ier = 911;
	nrec = 1;
	if (rkcom8_1.utask) {
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have not called UT, so you cannot use STA"
		    "T.", (ftnlen)51);
	    e_wsfi();
	} else {
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have not called CT, so you cannot use STA"
		    "T.", (ftnlen)51);
	    e_wsfi();
	}
	goto L20;
    }

/*  Set flag so that the routine can only be called once after a hard */
/*  failure from the integrator. */
    if (state == 5 || state == 6) {
	ier = -2;
    }

    *totfcn = rkcom2_1.svnfcn + rkcom2_1.nfcn;
    if (rkcom6_1.erason) {
	*totfcn += rkcom6_1.gnfcn;
    }
    *stpcst = (integer) rkcom5_1.cost;
    *stpsok = rkcom2_1.okstp;
    if (rkcom2_1.okstp <= 1) {
	*waste = 0.;
    } else {
	*waste = (doublereal) rkcom2_1.flstp / (doublereal) (rkcom2_1.flstp + 
		rkcom2_1.okstp);
    }
    *hnext = rkcom2_1.h__;

L20:

    rkmsg_(&ier, "STAT  ", &nrec, &flag__, (ftnlen)6);

    return 0;
} /* stat_ */

/* Subroutine */ int glberr_(doublereal *rmserr, doublereal *errmax, 
	doublereal *terrmx, doublereal *work)
{
    /* System generated locals */
    integer i__1;
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    double sqrt(doublereal);

    /* Local variables */
    integer l, ier, flag__, nrec, state;
    extern /* Subroutine */ int rkmsg_(integer *, char *, integer *, integer *
	    , ftnlen), rksit_(logical *, char *, integer *, ftnlen);

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

/*  If you are not familiar with the code GLBERR and how it is used in */
/*  conjunction with UT and CT to solve initial value problems, you should */
/*  study the document file rksuite.doc carefully before attempting to use */
/*  the code.  The following "Brief Reminder" is intended only to remind you */
/*  of the meaning, type, and size requirements of the arguments. */

/*  If ERRASS was set .TRUE. in the call to SETUP, then after any call to UT */
/*  or CT to advance the integration to TNOW or TWANT, the subroutine GLBERR */
/*  may be called to obtain an assessment of the true error of the integration. */
/*  At each step and for each solution component Y(L), a more accurate "true" */
/*  solution YT(L), an average magnitude "size(L)" of its size, and its error */
/*                abs(Y(L) - YT(L))/max("size(L)",THRES(L)) */
/*  are computed.  The assessment returned in RMSERR(L) is the RMS (root-mean- */
/*  square) average of the error in the Lth solution component over all steps */
/*  of the integration from TSTART through TNOW. */

/*  OUTPUT VARIABLES */

/*     RMSERR(*) - DOUBLE PRECISION array of length NEQ */
/*                 RMSERR(L) approximates the RMS average of the true error */
/*                 of the numerical solution for the Ith solution component, */
/*                 L = 1,2,...,NEQ.  The average is taken over all steps from */
/*                 TSTART to TNOW. */
/*     ERRMAX    - DOUBLE PRECISION */
/*                 The maximum (approximate) true error taken over all */
/*                 solution components and all steps from TSTART to TNOW. */
/*     TERRMX    - DOUBLE PRECISION */
/*                 First value of the independent variable where the */
/*                 (approximate) true error attains the maximum value ERRMAX. */

/*  WORKSPACE */

/*     WORK(*)   - DOUBLE PRECISION array as used in SETUP and UT or CT */
/*                 Do not alter the contents of this array. */

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block for Global Error Assessment .. */
/*     .. Common Block for Integrator Options .. */
/*     .. Common Block for Error Message .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --work;
    --rmserr;

    /* Function Body */
    ier = 1;
    nrec = 0;

/*  Is it permissible to call GLBERR? */

    rksit_(&c_true, "GLBERR", &state, (ftnlen)6);
    if (state == 911) {
	ier = 912;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** A catastrophic error has already been detected el"
		"sewhere.", (ftnlen)61);
	e_wsfi();
	goto L40;
    }
    if (state == -2) {
	ier = 911;
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have already made a call to GLBERR after a ha"
		"rd ", (ftnlen)56);
	do_fio(&c__1, " ** failure was reported from the integrator. You can"
		"not", (ftnlen)56);
	do_fio(&c__1, " ** call GLBERR again.", (ftnlen)22);
	e_wsfi();
	goto L40;
    }
    rksit_(&c_true, "CT", &state, (ftnlen)2);
    if (state == -1) {
	ier = 911;
	nrec = 1;
	if (rkcom8_1.utask) {
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have not yet called UT, so you cannot cal"
		    "l GLBERR.", (ftnlen)58);
	    e_wsfi();
	} else {
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have not yet called CT, so you cannot cal"
		    "l GLBERR.", (ftnlen)58);
	    e_wsfi();
	}
	goto L40;
    }

/*  Set flag so that the routine can only be called once after a hard */
/*  failure from the integrator. */
    if (state == 5 || state == 6) {
	ier = -2;
    }

/*  Check that ERRASS was set properly for error assessment in SETUP. */

    if (! rkcom6_1.erason) {
	ier = 911;
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** No error assessment is available since you did no"
		"t ", (ftnlen)55);
	do_fio(&c__1, " ** ask for it in your call to the routine SETUP.", (
		ftnlen)49);
	do_fio(&c__1, " ** Check your program carefully.", (ftnlen)33);
	e_wsfi();
	goto L40;
    }

/*  Check to see if the integrator has not actually taken a step. */

    if (rkcom2_1.okstp == 0) {
	ier = 911;
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** The integrator has not actually taken any success"
		"ful ", (ftnlen)57);
	do_fio(&c__1, " ** steps.  You cannot call GLBERR in this circumstan"
		"ce. ", (ftnlen)57);
	do_fio(&c__1, " ** Check your program carefully.", (ftnlen)33);
	e_wsfi();
	goto L40;
    }

/*  Compute RMS error and set output variables. */

    *errmax = rkcom6_1.maxerr;
    *terrmx = rkcom6_1.locmax;
    i__1 = rkcom1_1.neqn;
    for (l = 1; l <= i__1; ++l) {
	rmserr[l] = sqrt(work[rkcom6_1.przerr - 1 + l] / rkcom2_1.okstp);
/* L20: */
    }

L40:

    rkmsg_(&ier, "GLBERR", &nrec, &flag__, (ftnlen)6);

    return 0;
} /* glberr_ */

/* Subroutine */ int ct_(S_fp f, doublereal *tnow, doublereal *ynow, 
	doublereal *ypnow, doublereal *work, integer *cflag)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    double pow_dd(doublereal *, doublereal *), d_sign(doublereal *, 
	    doublereal *);

    /* Local variables */
    integer l, ier;
    doublereal err, tau, beta;
    static doublereal havg;
    integer nrec;
    logical main;
    doublereal hmin;
    extern /* Subroutine */ int step_(S_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, doublereal *, doublereal *, logical *);
    static integer ynew;
    doublereal htry, temp1, temp2, alpha;
    static integer ntend;
    extern /* Subroutine */ int stiff_(S_fp, doublereal *, integer *, logical 
	    *, integer *, doublereal *, integer *, integer *);
    integer state;
    extern /* Subroutine */ int rkmsg_(integer *, char *, integer *, integer *
	    , ftnlen);
    integer point;
    static integer ypold;
    extern /* Subroutine */ int rksit_(logical *, char *, integer *, ftnlen);
    logical phase1;
    static logical phase2;
    logical phase3, failed;
    static logical chkeff;
    static doublereal errold;
    logical toomch;
    static integer jflstp;
    doublereal ypnorm;
    extern /* Subroutine */ int truerr_(S_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

/*  If you are not familiar with the code CT and how it is used in */
/*  conjunction with SETUP to solve initial value problems, you should study */
/*  the document file rksuite.doc carefully before attempting to use the code. */
/*  The following "Brief Reminder" is intended only to remind you of the */
/*  meaning, type, and size requirements of the arguments. */

/*  NAME DECLARED IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM: */

/*     F         - name of the subroutine for evaluating the differential */
/*                 equations. */

/*  The subroutine F must have the form */

/*  SUBROUTINE F(T,Y,YP) */
/*  DOUBLE PRECISION T,Y(*),YP(*) */
/*     Using the input values of the independent variable T and the solution */
/*     components Y(*), for each L = 1,2,...,NEQ evaluate the differential */
/*     equation for the derivative of the Lth solution component and place the */
/*     value in YP(L).  Do not alter the input values of T and Y(*). */
/*  RETURN */
/*  END */

/*  OUTPUT VARIABLES */

/*     TNOW      - DOUBLE PRECISION */
/*                 Current value of the independent variable. */
/*     YNOW(*)   - DOUBLE PRECISION array of length NEQ */
/*                 Approximation to the true solution at TNOW. */
/*     YPNOW(*)  - DOUBLE PRECISION array of length NEQ */
/*                 Approximation to the first derivative of the */
/*                 true solution at TNOW. */

/*  WORKSPACE */

/*     WORK(*)   - DOUBLE PRECISION array as used in SETUP */
/*                 Do not alter the contents of this array. */

/*  OUTPUT VARIABLE */

/*     CFLAG     - INTEGER */

/*                       SUCCESS.  A STEP WAS TAKEN TO TNOW. */
/*                 = 1 - Complete success. */

/*                       "SOFT" FAILURES */
/*                 = 2 - Warning:  You have obtained an answer by integrating */
/*                       to TEND (TNOW = TEND).  You have done this at least */
/*                       100 times, and monitoring of the computation reveals */
/*                       that this way of getting output has degraded the */
/*                       efficiency of the code. If you really need answers at */
/*                       so many specific points, it would be more efficient to */
/*                       get them with INTRP.  (If METHOD = 3, you would need */
/*                       to change METHOD and restart from TNOW, YNOW(*) by a */
/*                       call to SETUP.)  If you wish to continue as you are, */
/*                       you may. */
/*                 = 3 - Warning:  A considerable amount of work has been */
/*                       expended. To continue the integration, just call */
/*                       CT again. */
/*                 = 4 - Warning:  It appears that this problem is "stiff". */
/*                       You really should change to another code that is */
/*                       intended for such problems, but if you insist, you */
/*                       can continue with CT by calling it again. */

/*                       "HARD" FAILURES */
/*                 = 5 - You are asking for too much accuracy. You cannot */
/*                       continue integrating this problem. */
/*                 = 6 - The global error assessment may not be reliable beyond */
/*                       the current point in the integration.  You cannot */
/*                       continue integrating this problem. */

/*                       "CATASTROPHIC" FAILURES */
/*                 = 911 - The nature of the catastrophe is reported on */
/*                         the standard output channel. Unless special */
/*                         provision was made in advance (see rksuite.doc), */
/*                         the computation then comes to a STOP. */

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block for General Workspace Pointers .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Global Error Assessment .. */
/*     .. Common Block for Environment Parameters .. */
/*     .. Common Block for Integrator Options .. */
/*     .. Common Block for Error Message .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Save statement .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --work;
    --ypnow;
    --ynow;

    /* Function Body */
    ier = 1;
    nrec = 0;

/*  Is it permissible to call CT? */

    rksit_(&c_true, "SETUP", &state, (ftnlen)5);
    if (state == 911) {
	ier = 912;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** A catastrophic error has already been detected el"
		"sewhere.", (ftnlen)61);
	e_wsfi();
	goto L180;
    }
    if (state == -1) {
	ier = 911;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have not called SETUP, so you cannot use CT.", 
		(ftnlen)52);
	e_wsfi();
	goto L180;
    }
    if (rkcom8_1.utask) {
	rksit_(&c_true, "UT", &state, (ftnlen)2);
	if (state != -2) {
	    ier = 911;
	    nrec = 2;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A/A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have called CT after you specified in SET"
		    "UP that ", (ftnlen)57);
	    do_fio(&c__1, " ** you were going to use UT. This is not permitt"
		    "ed.", (ftnlen)52);
	    e_wsfi();
	    rkcom8_1.utask = FALSE_;
	    goto L180;
	}
    }
    rksit_(&c_true, "CT    ", &state, (ftnlen)6);
    if (state == 5 || state == 6) {
	ier = 911;
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** CT has already returned with a flag value of 5 or"
		" 6.", (ftnlen)56);
	do_fio(&c__1, " ** You cannot continue integrating this problem. You"
		" must ", (ftnlen)59);
	do_fio(&c__1, " ** call SETUP to start another problem.", (ftnlen)40);
	e_wsfi();
	goto L180;
    }

    if (rkcom2_1.first) {

/*  First call in an integration -- initialize everything. */

	chkeff = FALSE_;
	ntend = 0;
	jflstp = 0;

/*  A scratch area of WORK(*) starting at PRSCR is used to hold two */
/*  arrays in this subroutine: the higher order approximate solution at */
/*  the end of a step and the approximate derivative of the solution at */
/*  the end of the last step. To make this clear, local pointers YNEW and */
/*  YPOLD are used. */
	ynew = rkcom3_1.prscr;
	ypold = rkcom3_1.prscr;

/*  For this first step T was initialized to TSTRT in SETUP and the */
/*  starting values YSTART(*) were loaded into the area of WORK(*) reserved */
/*  for the current solution approximation starting at location PRY. The */
/*  derivative is now computed and stored in WORK(*) starting at PRYP. */
/*  Subsequently these arrays are copied to the output vectors YNOW(*) */
/*  and YPNOW(*). */
	(*f)(&rkcom2_1.t, &work[rkcom3_1.pry], &work[rkcom3_1.pryp]);
	++rkcom2_1.nfcn;
	i__1 = rkcom1_1.neqn;
	for (l = 1; l <= i__1; ++l) {
	    ynow[l] = work[rkcom3_1.pry - 1 + l];
	    ypnow[l] = work[rkcom3_1.pryp - 1 + l];
/* L20: */
	}

/*  Set dependent variables for error assessment. */
	if (rkcom6_1.erason) {
	    i__1 = rkcom1_1.neqn;
	    for (l = 1; l <= i__1; ++l) {
		work[rkcom6_1.przy - 1 + l] = ynow[l];
		work[rkcom6_1.przyp - 1 + l] = ypnow[l];
/* L40: */
	    }
	}

/*  The weights for the control of the error depend on the size of the */
/*  solution at the beginning and at the end of the step. On the first */
/*  step we do not have all this information. Whilst determining the */
/*  initial step size we initialize the weight vector to the larger of */
/*  abs(Y(L)) and the threshold for this component. */
	i__1 = rkcom1_1.neqn;
	for (l = 1; l <= i__1; ++l) {
/* Computing MAX */
	    d__2 = (d__1 = ynow[l], abs(d__1)), d__3 = work[rkcom3_1.prthrs - 
		    1 + l];
	    work[rkcom3_1.prwt - 1 + l] = max(d__2,d__3);
/* L60: */
	}

/*  If HSTRT is equal to zero, the code is to find an on-scale initial step */
/*  size H.  CT has an elaborate scheme of three phases for finding such an H, */
/*  and some preparations are made earlier.  In SETUP an upper bound is placed */
/*  on H that reflects the scale of the independent variable. When UTASK is */
/*  .TRUE., UT refines this bound using the first output point.  Here in CT */
/*  PHASE1 applies a rule of thumb based on the error control, the order of the */
/*  the formula, and the size of the initial slope to get a crude approximation */
/*  to an on-scale H.  PHASE2 may reduce H in the course of taking the first */
/*  step.  PHASE3 repeatedly adjusts H and retakes the first step until H is */
/*  on scale. */

/*  A guess for the magnitude of the first step size H can be provided to SETUP */
/*  as HSTART.  If it is too big or too small, it is ignored and the automatic */
/*  determination of an on-scale initial step size is activated.  If it is */
/*  acceptable, H is set to HSTART in SETUP.  Even when H is supplied to CT, */
/*  PHASE3 of the scheme for finding an on-scale initial step size is made */
/*  active so that the code can deal with a bad guess. */

	phase1 = rkcom1_1.hstrt == 0.;
	phase2 = phase1;
	phase3 = TRUE_;
	if (phase1) {
	    rkcom2_1.h__ = abs(rkcom2_1.h__);
	    ypnorm = 0.;
	    i__1 = rkcom1_1.neqn;
	    for (l = 1; l <= i__1; ++l) {
		if ((d__1 = ynow[l], abs(d__1)) != 0.) {
/* Computing MAX */
		    d__2 = ypnorm, d__3 = (d__1 = ypnow[l], abs(d__1)) / work[
			    rkcom3_1.prwt - 1 + l];
		    ypnorm = max(d__2,d__3);
		}
/* L80: */
	    }
	    tau = pow_dd(&rkcom1_1.tolr, &rkcom5_1.expon);
	    if (rkcom2_1.h__ * ypnorm > tau) {
		rkcom2_1.h__ = tau / ypnorm;
	    }
/* Computing MAX */
/* Computing MAX */
	    d__3 = abs(rkcom1_1.tstrt), d__4 = abs(rkcom1_1.tnd);
	    d__1 = rkcom7_1.tiny, d__2 = rkcom5_1.toosml * max(d__3,d__4);
	    hmin = max(d__1,d__2);
	    rkcom2_1.h__ = rkcom1_1.dir * max(rkcom2_1.h__,hmin);
	    phase1 = FALSE_;
	}

    } else {

/* Continuation call */

	if (rkcom2_1.last) {
	    ier = 911;
	    nrec = 3;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A,D13.5,A/A/A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have already reached TEND ( = ", (ftnlen)
		    38);
	    do_fio(&c__1, (char *)&rkcom1_1.tnd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ").", (ftnlen)2);
	    do_fio(&c__1, " ** To integrate further with the same problem yo"
		    "u must ", (ftnlen)56);
	    do_fio(&c__1, " ** call the routine RESET with a new value of TE"
		    "ND.", (ftnlen)52);
	    e_wsfi();
	    goto L180;
	}
    }

/*  Begin computation of a step here. */

    failed = FALSE_;

L100:
    d__1 = abs(rkcom2_1.h__);
    rkcom2_1.h__ = d_sign(&d__1, &rkcom1_1.dir);

/*  Reduce the step size if necessary so that the code will not step */
/*  past TND.  "Look ahead" to prevent unnecessarily small step sizes. */
    rkcom2_1.last = rkcom1_1.dir * (rkcom2_1.t + rkcom2_1.h__ - rkcom1_1.tnd) 
	    >= 0.;
    if (rkcom2_1.last) {
	rkcom2_1.h__ = rkcom1_1.tnd - rkcom2_1.t;
    } else if (rkcom1_1.dir * (rkcom2_1.t + rkcom2_1.h__ * 2. - rkcom1_1.tnd) 
	    >= 0.) {
	rkcom2_1.h__ = (rkcom1_1.tnd - rkcom2_1.t) / 2.;
    }

/*  When the integrator is at T and attempts a step of H, the function */
/*  defining the differential equations will be evaluated at a number of */
/*  arguments between T and T+H.  If H is too small, these arguments cannot */
/*  be clearly distinguished in the precision available. */

/* Computing MAX */
/* Computing MAX */
    d__4 = abs(rkcom2_1.t), d__5 = (d__1 = rkcom2_1.t + rkcom2_1.h__, abs(
	    d__1));
    d__2 = rkcom7_1.tiny, d__3 = rkcom5_1.toosml * max(d__4,d__5);
    hmin = max(d__2,d__3);
    if (abs(rkcom2_1.h__) < hmin) {
	ier = 5;
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A,D13.5,A,D13.5/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** In order to satisfy your error requirements CT wo"
		"uld ", (ftnlen)57);
	do_fio(&c__1, " ** have to use a step size of ", (ftnlen)31);
	do_fio(&c__1, (char *)&rkcom2_1.h__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " at TNOW = ", (ftnlen)11);
	do_fio(&c__1, (char *)&rkcom2_1.t, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " ** This is too small for the machine precision.", (
		ftnlen)48);
	e_wsfi();
	goto L180;
    }

/*  Monitor the impact of output on the efficiency of the integration. */

    if (chkeff) {
	++ntend;
	if (ntend >= 100 && ntend >= rkcom2_1.okstp / 3) {
	    ier = 2;
	    nrec = 6;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A/A/A/A/A/A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** More than 100 output points have been obtaine"
		    "d by ", (ftnlen)54);
	    do_fio(&c__1, " ** integrating to TEND.  They have been sufficie"
		    "ntly close ", (ftnlen)60);
	    do_fio(&c__1, " ** to one another that the efficiency of the int"
		    "egration has ", (ftnlen)62);
	    do_fio(&c__1, " ** been degraded. It would probably be (much) mo"
		    "re efficient ", (ftnlen)62);
	    do_fio(&c__1, " ** to obtain output by interpolating with INTRP "
		    "(after ", (ftnlen)56);
	    do_fio(&c__1, " ** changing to METHOD=2 if you are using METHOD "
		    "= 3).", (ftnlen)54);
	    e_wsfi();
	    ntend = 0;
	    goto L180;
	}
    }

/*  Check for stiffness and for too much work.  Stiffness can be */
/*  checked only after a successful step. */

    if (! failed) {

/*  Check for too much work. */
	toomch = rkcom2_1.nfcn > 5000;
	if (toomch) {
	    ier = 3;
	    nrec = 3;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A,I6,A/A/A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** Approximately ", (ftnlen)18);
	    do_fio(&c__1, (char *)&c__5000, (ftnlen)sizeof(integer));
	    do_fio(&c__1, " function evaluations have been ", (ftnlen)32);
	    do_fio(&c__1, " ** used to compute the solution since the integr"
		    "ation ", (ftnlen)55);
	    do_fio(&c__1, " ** started or since this message was last printe"
		    "d.", (ftnlen)51);
	    e_wsfi();

/*  After this warning message, NFCN is reset to permit the integration */
/*  to continue.  The total number of function evaluations in the primary */
/*  integration is SVNFCN + NFCN. */
	    rkcom2_1.svnfcn += rkcom2_1.nfcn;
	    rkcom2_1.nfcn = 0;
	}

/*  Check for stiffness.  NREC is passed on to STIFF because when */
/*  TOOMCH = .TRUE. and stiffness is diagnosed, the message about too */
/*  much work is augmented inside STIFF to explain that it is due to */
/*  stiffness. */
	stiff_((S_fp)f, &havg, &jflstp, &toomch, &c__5000, &work[1], &ier, &
		nrec);

	if (ier != 1) {
	    goto L180;
	}
    }

/*  Take a step.  Whilst finding an on-scale H (PHASE2 = .TRUE.), the input */
/*  value of H might be reduced (repeatedly), but it will not be reduced */
/*  below HMIN.  The local error is estimated, a weight vector is formed, */
/*  and a weighted maximum norm, ERR, of the local error is returned. */
/*  The variable MAIN is input as .TRUE. to tell STEP that this is the */
/*  primary, or "main", integration. */

/*  H resides in the common block /RKCOM2/ which is used by both CT and STEP; */
/*  since it may be changed inside STEP, a local copy is made to ensure */
/*  portability of the code. */

    main = TRUE_;
    htry = rkcom2_1.h__;
    step_((S_fp)f, &rkcom1_1.neqn, &rkcom2_1.t, &work[rkcom3_1.pry], &work[
	    rkcom3_1.pryp], &work[rkcom3_1.prstgs], &rkcom1_1.tolr, &htry, &
	    work[rkcom3_1.prwt], &work[ynew], &work[rkcom3_1.prerst], &err, &
	    main, &hmin, &work[rkcom3_1.prthrs], &phase2);
    rkcom2_1.h__ = htry;

/*  Compare the norm of the local error to the tolerance. */

    if (err > rkcom1_1.tolr) {

/*  Failed step.  Reduce the step size and try again. */

/*  First step:  Terminate PHASE3 of the search for an on-scale step size. */
/*               The step size is not on scale, so ERR may not be accurate; */
/*               reduce H by a fixed factor.  Failed attempts to take the */
/*               first step are not counted. */
/*  Later step:  Use ERR to compute an "optimal" reduction of H.  More than */
/*               one failure indicates a difficulty with the problem and an */
/*               ERR that may not be accurate, so reduce H by a fixed factor. */

	if (rkcom2_1.first) {
	    phase3 = FALSE_;
	    alpha = rkcom5_1.rs1;
	} else {
	    ++rkcom2_1.flstp;
	    ++jflstp;
	    if (failed) {
		alpha = rkcom5_1.rs1;
	    } else {
		d__1 = rkcom1_1.tolr / err;
		alpha = rkcom5_1.safety * pow_dd(&d__1, &rkcom5_1.expon);
		alpha = max(alpha,rkcom5_1.rs1);
	    }
	}
	rkcom2_1.h__ = alpha * rkcom2_1.h__;
	failed = TRUE_;
	goto L100;
    }

/*  Successful step. */

/*  Predict a step size appropriate for the next step.  After the first */
/*  step the prediction can be refined using an idea of H.A. Watts that */
/*  takes account of how well the prediction worked on the previous step. */
    d__1 = err / rkcom1_1.tolr;
    beta = pow_dd(&d__1, &rkcom5_1.expon);
    if (! rkcom2_1.first) {
	temp1 = pow_dd(&err, &rkcom5_1.expon) / rkcom2_1.h__;
	temp2 = pow_dd(&errold, &rkcom5_1.expon) / rkcom2_1.hold;
	if (temp1 < temp2 * 100. && temp2 < temp1 * 100.) {
	    beta *= temp1 / temp2;
	}
    }
    alpha = rkcom5_1.rs3;
    if (rkcom5_1.safety < beta * alpha) {
	alpha = rkcom5_1.safety / beta;
    }

/*  On the first step a search is made for an on-scale step size.  PHASE2 */
/*  of the scheme comes to an end here because a step size has been found */
/*  that is both successful and has a credible local error estimate. Except */
/*  in the special case that the first step is also the last, the step is */
/*  repeated in PHASE3 as long as an increase greater than RS2 appears */
/*  possible.  An increase as big as RS3 is permitted.  A step failure */
/*  terminates PHASE3. */

    if (rkcom2_1.first) {
	phase2 = FALSE_;
	phase3 = phase3 && ! rkcom2_1.last && alpha > rkcom5_1.rs2;
	if (phase3) {
	    rkcom2_1.h__ = alpha * rkcom2_1.h__;
	    goto L100;
	}
    }

/*  After getting on scale, step size changes are more restricted. */
    alpha = min(alpha,rkcom5_1.rs);
    if (failed) {
	alpha = min(alpha,1.);
    }
    alpha = max(alpha,rkcom5_1.rs1);
    rkcom2_1.hold = rkcom2_1.h__;
    rkcom2_1.h__ = alpha * rkcom2_1.h__;

/*  For the diagnosis of stiffness, an average accepted step size, HAVG, */
/*  must be computed and SAVEd. */
    if (rkcom2_1.first) {
	havg = rkcom2_1.hold;
    } else {
	havg = havg * .9 + rkcom2_1.hold * .1;
    }

    rkcom2_1.first = FALSE_;
    errold = err;
    rkcom2_1.told = rkcom2_1.t;

/*  Take care that T is set to precisely TND when the end of the */
/*  integration is reached. */
    if (rkcom2_1.last) {
	rkcom2_1.t = rkcom1_1.tnd;
    } else {
	rkcom2_1.t += rkcom2_1.hold;
    }

/*  Increment counter on accepted steps.  Note that successful steps */
/*  that are repeated whilst getting on scale are not counted. */
    ++rkcom2_1.okstp;

/*  Advance the current solution and its derivative.  (Stored in WORK(*) */
/*  with the first location being PRY and PRYP, respectively.)  Update the */
/*  previous solution and its derivative.  (Stored in WORK(*) with the first */
/*  location being PRYOLD and YPOLD, respectively.)  Note that the previous */
/*  derivative will overwrite YNEW(*). */

    i__1 = rkcom1_1.neqn;
    for (l = 1; l <= i__1; ++l) {
	work[rkcom3_1.pryold - 1 + l] = work[rkcom3_1.pry - 1 + l];
	work[rkcom3_1.pry - 1 + l] = work[ynew - 1 + l];
	work[ypold - 1 + l] = work[rkcom3_1.pryp - 1 + l];
/* L120: */
    }

    if (rkcom5_1.fsal) {

/*  When FSAL = .TRUE., YP(*) is the last stage of the step. */
	point = rkcom3_1.prstgs + (rkcom5_1.lststg - 1) * rkcom1_1.neqn;
	i__1 = rkcom1_1.neqn;
	for (l = 1; l <= i__1; ++l) {
	    work[rkcom3_1.pryp - 1 + l] = work[point - 1 + l];
/* L140: */
	}
    } else {

/*  Call F to evaluate YP(*). */
	(*f)(&rkcom2_1.t, &work[rkcom3_1.pry], &work[rkcom3_1.pryp]);
	++rkcom2_1.nfcn;
    }

/*  If global error assessment is desired, advance the secondary */
/*  integration from TOLD to T. */

    if (rkcom6_1.erason) {
	truerr_((S_fp)f, &rkcom1_1.neqn, &work[rkcom3_1.pry], &rkcom1_1.tolr, 
		&work[rkcom3_1.prwt], &work[rkcom6_1.przy], &work[
		rkcom6_1.przyp], &work[rkcom6_1.przerr], &work[
		rkcom6_1.przynu], &work[rkcom6_1.przers], &work[
		rkcom6_1.przstg], &ier);
	if (ier == 6) {

/*  The global error estimating procedure has broken down. Treat it as a */
/*  failed step. The solution and derivative are reset to their values at */
/*  the beginning of the step since the last valid error assessment refers */
/*  to them. */
	    --rkcom2_1.okstp;
	    rkcom6_1.erasfl = TRUE_;
	    rkcom2_1.last = FALSE_;
	    rkcom2_1.t = rkcom2_1.told;
	    rkcom2_1.h__ = rkcom2_1.hold;
	    i__1 = rkcom1_1.neqn;
	    for (l = 1; l <= i__1; ++l) {
		work[rkcom3_1.pry - 1 + l] = work[rkcom3_1.pryold - 1 + l];
		work[rkcom3_1.pryp - 1 + l] = work[ypold - 1 + l];
/* L160: */
	    }
	    if (rkcom2_1.okstp > 1) {
		nrec = 2;
		ici__1.icierr = 0;
		ici__1.icirnum = 10;
		ici__1.icirlen = 80;
		ici__1.iciunit = rkcom9_1.rec;
		ici__1.icifmt = "(A/A,D13.5,A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, " ** The global error assessment may not be re"
			"liable for T past ", (ftnlen)63);
		do_fio(&c__1, " ** TNOW = ", (ftnlen)11);
		do_fio(&c__1, (char *)&rkcom2_1.t, (ftnlen)sizeof(doublereal))
			;
		do_fio(&c__1, ".  The integration is being terminated.", (
			ftnlen)39);
		e_wsfi();
	    } else {
		nrec = 2;
		ici__1.icierr = 0;
		ici__1.icirnum = 10;
		ici__1.icirlen = 80;
		ici__1.iciunit = rkcom9_1.rec;
		ici__1.icifmt = "(A/A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, " ** The global error assessment algorithm fai"
			"led at the start", (ftnlen)61);
		do_fio(&c__1, " ** the integration.  The integration is bein"
			"g terminated.", (ftnlen)58);
		e_wsfi();
	    }
	    goto L180;
	}
    }


/*  Exit point for CT */

L180:

/*  Set the output variables and flag that interpolation is permitted */

    if (ier < 911) {
	*tnow = rkcom2_1.t;
	rkcom2_1.last = *tnow == rkcom1_1.tnd;
	chkeff = rkcom2_1.last;
	i__1 = rkcom1_1.neqn;
	for (l = 1; l <= i__1; ++l) {
	    ynow[l] = work[rkcom3_1.pry - 1 + l];
	    ypnow[l] = work[rkcom3_1.pryp - 1 + l];
/* L200: */
	}
	if (ier == 1) {
	    state = -2;
	    rksit_(&c_false, "INTRP", &state, (ftnlen)5);
	}
    }

/*  Call RKMSG to report what happened and set CFLAG. */

    rkmsg_(&ier, "CT    ", &nrec, cflag, (ftnlen)6);

    return 0;
} /* ct_ */

/* Subroutine */ int intrp_(doublereal *twant, char *reqest, integer *nwant, 
	doublereal *ywant, doublereal *ypwant, U_fp f, doublereal *work, 
	doublereal *wrkint, integer *lenint, ftnlen reqest_len)
{
    /* Initialized data */

    static integer nwntsv = -1;

    /* System generated locals */
    integer i__1, i__2;
    logical L__1;
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    integer ier, flag__, ichk, nrec;
    extern /* Subroutine */ int evali_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, char *, integer *, doublereal *, doublereal *, 
	    ftnlen), formi_(U_fp, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     doublereal *, doublereal *);
    integer state;
    extern /* Subroutine */ int rkmsg_(integer *, char *, integer *, integer *
	    , ftnlen), rksit_(logical *, char *, integer *, ftnlen);
    integer state1;
    char reqst1[1];
    logical legalr;
    static logical inintp;
    static integer startp;

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

/*  If you are not familiar with the code INTRP and how it is used in */
/*  conjunction with CT to solve initial value problems, you should study the */
/*  document file rksuite.doc carefully before attempting to use the code. The */
/*  following "Brief Reminder" is intended only to remind you of the meaning, */
/*  type, and size requirements of the arguments. */

/*  When integrating with METHOD = 1 or 2, answers may be obtained inexpensively */
/*  between steps by interpolation. INTRP is called after a step by CT from a */
/*  previous value of T, called TOLD below, to the current value of T to get */
/*  an answer at TWANT. You can specify any value of TWANT you wish, but */
/*  specifying a value outside the interval [TOLD,T] is likely to yield */
/*  answers with unsatisfactory accuracy. */

/*  INPUT VARIABLE */

/*     TWANT     - DOUBLE PRECISION */
/*                 The value of the independent variable where a solution */
/*                 is desired. */

/*  The interpolant is to be evaluated at TWANT to approximate the solution */
/*  and/or its first derivative there.  There are three cases: */

/*  INPUT VARIABLE */

/*     REQEST    - CHARACTER*(*) */
/*                 Only the first character of REQEST is significant. */
/*                 REQEST(1:1) = `S' or `s'- compute approximate `S'olution */
/*                                           only. */
/*                             = `D' or `d'- compute approximate first */
/*                                           `D'erivative of the solution only. */
/*                             = `B' or `b'- compute `B'oth approximate solution */
/*                                           and first derivative. */
/*                 Constraint: REQEST(1:1) must be `S',`s',`D',`d',`B' or `b'. */

/*  If you intend to interpolate at many points, you should arrange for the */
/*  the interesting components to be the first ones because the code */
/*  approximates only the first NWANT components. */

/*  INPUT VARIABLE */

/*     NWANT     - INTEGER */
/*                 Only the first NWANT components of the answer are to be */
/*                 computed. */
/*                 Constraint:  NEQ >= NWANT >= 1 */

/*  OUTPUT VARIABLES */

/*     YWANT(*)  - DOUBLE PRECISION array of length NWANT */
/*                 Approximation to the first NWANT components of the true */
/*                 solution at TWANT when REQESTed. */
/*     YPWANT(*) - DOUBLE PRECISION array of length NWANT */
/*                 Approximation to the first NWANT components of the first */
/*                 derivative of the true solution at TWANT when REQESTed. */

/*  NAME DECLARED IN AN EXTERNAL STATEMENT IN THE PROGRAM CALLING INTRP: */

/*     F         - name of the subroutine for evaluating the differential */
/*                 equations as provided to CT. */

/*  WORKSPACE */

/*     WORK(*)   - DOUBLE PRECISION array as used in SETUP and CT */
/*                 Do not alter the contents of this array. */

/*     WRKINT(*) - DOUBLE PRECISION array of length LENINT */
/*                 Do not alter the contents of this array. */

/*     LENINT    - INTEGER */
/*                 Length of WRKINT. If */
/*                 METHOD = 1, LENINT must be at least 1 */
/*                        = 2, LENINT must be at least NEQ+MAX(NEQ,5*NWANT) */
/*                        = 3--CANNOT BE USED WITH THIS SUBROUTINE */

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block for General Workspace Pointers .. */
/*     .. Common Block to hold Formula Definitions .. */
/*     .. Common Block for Integrator Options .. */
/*     .. Common Block for Error Message .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Save statement .. */
/*     .. Data statements .. */
    /* Parameter adjustments */
    --wrkint;
    --work;
    --ypwant;
    --ywant;

    /* Function Body */
/*     .. Executable Statements .. */

    ier = 1;
    nrec = 0;

/*  Is it permissible to call INTRP? */

    rksit_(&c_true, "CT", &state, (ftnlen)2);
    if (state == 911) {
	ier = 912;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** A catastrophic error has already been detected el"
		"sewhere.", (ftnlen)61);
	e_wsfi();
	goto L20;
    }
    if (rkcom8_1.utask) {
	rksit_(&c_true, "UT", &state1, (ftnlen)2);
	if (state1 != -2) {
	    ier = 911;
	    nrec = 2;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A/A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have called INTRP after you specified to "
		    "SETUP ", (ftnlen)55);
	    do_fio(&c__1, " ** that you were going to use UT. This is not pe"
		    "rmitted.", (ftnlen)57);
	    e_wsfi();
	    goto L20;
	}
    }
    if (state == -1) {
	ier = 911;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have not called CT, so you cannot use INTRP.", 
		(ftnlen)52);
	e_wsfi();
	goto L20;
    }
    if (state > 1) {
	ier = 911;
	nrec = 2;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** CT has returned with a flag value greater than 1.",
		 (ftnlen)53);
	do_fio(&c__1, " ** You cannot call INTRP in this circumstance.", (
		ftnlen)47);
	e_wsfi();
	goto L20;
    }

/*  Check input */

    *(unsigned char *)reqst1 = *(unsigned char *)reqest;
    legalr = *(unsigned char *)reqst1 == 'S' || *(unsigned char *)reqst1 == 
	    's' || *(unsigned char *)reqst1 == 'D' || *(unsigned char *)
	    reqst1 == 'd' || *(unsigned char *)reqst1 == 'B' || *(unsigned 
	    char *)reqst1 == 'b';
    if (! legalr) {
	ier = 911;
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A,A,A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have set the first character of ", (ftnlen)40);
	do_fio(&c__1, " ** REQEST to be '", (ftnlen)18);
	do_fio(&c__1, reqst1, (ftnlen)1);
	do_fio(&c__1, "'. It must be one of ", (ftnlen)21);
	do_fio(&c__1, " ** 'S','s','D','d','B' or 'b'.", (ftnlen)31);
	e_wsfi();
	goto L20;
    }

    if (*nwant > rkcom1_1.neqn) {
	ier = 911;
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A,I6,A/A,I6,A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have specified the value of NWANT to be ", (
		ftnlen)48);
	do_fio(&c__1, (char *)&(*nwant), (ftnlen)sizeof(integer));
	do_fio(&c__1, ". This", (ftnlen)6);
	do_fio(&c__1, " ** is greater than ", (ftnlen)20);
	do_fio(&c__1, (char *)&rkcom1_1.neqn, (ftnlen)sizeof(integer));
	do_fio(&c__1, ", which is the number of equations ", (ftnlen)35);
	do_fio(&c__1, " ** in the system being integrated.", (ftnlen)35);
	e_wsfi();
	goto L20;
    } else if (*nwant < 1) {
	ier = 911;
	nrec = 3;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A,I6,A/A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have specified the value of NWANT to be ", (
		ftnlen)48);
	do_fio(&c__1, (char *)&(*nwant), (ftnlen)sizeof(integer));
	do_fio(&c__1, ", but ", (ftnlen)6);
	do_fio(&c__1, " ** this is less than 1. You cannot interpolate a zer"
		"o or ", (ftnlen)58);
	do_fio(&c__1, " ** negative number of components.", (ftnlen)34);
	e_wsfi();
	goto L20;
    }

    if (rkcom4_1.methd == 1) {
	if (*lenint < 1) {
	    ier = 911;
	    nrec = 2;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A,I6,A/A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have specified LENINT to be ", (ftnlen)36);
	    do_fio(&c__1, (char *)&(*lenint), (ftnlen)sizeof(integer));
	    do_fio(&c__1, ".", (ftnlen)1);
	    do_fio(&c__1, " ** This is too small. LENINT must be at least 1.",
		     (ftnlen)49);
	    e_wsfi();
	    goto L20;
	}
	startp = 1;
    } else if (rkcom4_1.methd == 2) {
/* Computing MAX */
	i__1 = rkcom1_1.neqn, i__2 = *nwant * 5;
	ichk = rkcom1_1.neqn + max(i__1,i__2);
	if (*lenint < ichk) {
	    ier = 911;
	    nrec = 3;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A,I6,A/A/A,I6,A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have specified LENINT to be ", (ftnlen)36);
	    do_fio(&c__1, (char *)&(*lenint), (ftnlen)sizeof(integer));
	    do_fio(&c__1, ". This is too", (ftnlen)13);
	    do_fio(&c__1, " ** small. NINT must be at least NEQ + MAX(NEQ, 5"
		    "*NWANT) ", (ftnlen)57);
	    do_fio(&c__1, " ** which is ", (ftnlen)13);
	    do_fio(&c__1, (char *)&ichk, (ftnlen)sizeof(integer));
	    do_fio(&c__1, ".", (ftnlen)1);
	    e_wsfi();
	    goto L20;
	}
	startp = rkcom1_1.neqn + 1;
    } else if (rkcom4_1.methd == 3) {
	ier = 911;
	nrec = 5;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A/A/A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have been using CT with METHOD = 3 to integra"
		"te your  ", (ftnlen)62);
	do_fio(&c__1, " ** equations. You have just called INTRP, but interp"
		"olation  ", (ftnlen)62);
	do_fio(&c__1, " ** is not available for this METHOD. Either use METH"
		"OD = 2,  ", (ftnlen)62);
	do_fio(&c__1, " ** for which interpolation is available, or use RESE"
		"T to make", (ftnlen)62);
	do_fio(&c__1, " ** CT step exactly to the points where you want outp"
		"ut.", (ftnlen)56);
	e_wsfi();
	goto L20;
    }

/*  Has the interpolant been initialised for this step? */

    rksit_(&c_true, "INTRP ", &state, (ftnlen)6);
    inintp = state != -2;

/*  Some initialization must be done before interpolation is possible. */
/*  To reduce the overhead, the interpolating polynomial is formed for */
/*  the first NWANT components.  In the unusual circumstance that NWANT */
/*  is changed while still interpolating within the span of the current */
/*  step, the scheme must be reinitialized to accomodate the additional */
/*  components. */

    if (! inintp || *nwant != nwntsv) {

/*  At present the derivative of the solution at the previous step, YPOLD(*), */
/*  is stored in the scratch array area starting at PRSCR. In the case of */
/*  METHD = 1 we can overwrite the stages. */

	if (rkcom4_1.methd == 1) {
	    L__1 = ! inintp;
	    formi_((U_fp)f, &rkcom1_1.neqn, nwant, &work[rkcom3_1.pry], &work[
		    rkcom3_1.pryp], &work[rkcom3_1.pryold], &work[
		    rkcom3_1.prscr], &work[rkcom3_1.prstgs], &L__1, &work[
		    rkcom3_1.prstgs], &work[rkcom3_1.prstgs]);
	} else {
	    L__1 = ! inintp;
	    formi_((U_fp)f, &rkcom1_1.neqn, nwant, &work[rkcom3_1.pry], &work[
		    rkcom3_1.pryp], &work[rkcom3_1.pryold], &work[
		    rkcom3_1.prscr], &work[rkcom3_1.prstgs], &L__1, &wrkint[1]
		    , &wrkint[startp]);
	}

/*  Set markers to show that interpolation has been initialized for */
/*  NWANT components. */
	nwntsv = *nwant;
	inintp = TRUE_;
    }

/*  The actual evaluation of the interpolating polynomial and/or its first */
/*  derivative is done in EVALI. */

    if (rkcom4_1.methd == 1) {
	evali_(&work[rkcom3_1.pry], &work[rkcom3_1.pryp], &work[
		rkcom3_1.prstgs], twant, reqest, nwant, &ywant[1], &ypwant[1],
		 reqest_len);
    } else {
	evali_(&work[rkcom3_1.pry], &work[rkcom3_1.pryp], &wrkint[startp], 
		twant, reqest, nwant, &ywant[1], &ypwant[1], reqest_len);
    }

L20:

    rkmsg_(&ier, "INTRP ", &nrec, &flag__, (ftnlen)6);

    return 0;
} /* intrp_ */

/* Subroutine */ int reset_(doublereal *tendnu)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    integer ier, flag__, nrec;
    doublereal hmin, tdiff;
    integer state;
    extern /* Subroutine */ int rkmsg_(integer *, char *, integer *, integer *
	    , ftnlen), rksit_(logical *, char *, integer *, ftnlen);
    integer state1;

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

/*  If you are not familiar with the code RESET and how it is used in */
/*  conjunction with CT to solve initial value problems, you should study the */
/*  document file rksuite.doc carefully before attempting to use the code. The */
/*  following "Brief Reminder" is intended only to remind you of the meaning, */
/*  type, and size requirements of the arguments. */

/*  The integration using CT proceeds from TSTART in the direction of TEND, and */
/*  is now at TNOW.  To reset TEND to a new value TENDNU, just call RESET with */
/*  TENDNU as the argument.  You must continue integrating in the same */
/*  direction, so the sign of (TENDNU - TNOW) must be the same as the sign of */
/*  (TEND - TSTART). To change direction you must restart by a call to SETUP. */

/*  INPUT VARIABLE */

/*     TENDNU    - DOUBLE PRECISION */
/*                 The new value of TEND. */
/*                 Constraint: TENDNU and TNOW must be clearly distinguishable */
/*                 in the precision used.  The sign of (TENDNU - TNOW) must be */
/*                 the same as the sign of (TEND - TSTART). */

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
/*     .. Scalar Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Environment Parameters .. */
/*     .. Common Block for Integrator Options .. */
/*     .. Common Block for Error Message .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */
    ier = 1;
    nrec = 0;

/*  Is it permissible to call RESET? */

    rksit_(&c_true, "CT", &state, (ftnlen)2);
    if (state == 911) {
	ier = 912;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** A catastrophic error has already been detected el"
		"sewhere.", (ftnlen)61);
	e_wsfi();
	goto L20;
    }
    if (rkcom8_1.utask) {
	rksit_(&c_true, "UT", &state1, (ftnlen)2);
	if (state1 != -2) {
	    ier = 911;
	    nrec = 2;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A/A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** You have called RESET after you specified to "
		    "SETUP that ", (ftnlen)60);
	    do_fio(&c__1, " ** you were going to use UT. This is not permitt"
		    "ed.", (ftnlen)52);
	    e_wsfi();
	    goto L20;
	}
    }
    if (state == -1) {
	ier = 911;
	nrec = 1;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** You have not called CT, so you cannot use RESET.", 
		(ftnlen)52);
	e_wsfi();
	goto L20;
    }
    if (state == 5 || state == 6) {
	ier = 911;
	nrec = 2;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A,I1,A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** CT has returned with CFLAG =  ", (ftnlen)34);
	do_fio(&c__1, (char *)&state, (ftnlen)sizeof(integer));
	do_fio(&c__1, ".", (ftnlen)1);
	do_fio(&c__1, " ** You cannot call RESET in this circumstance.", (
		ftnlen)47);
	e_wsfi();
	goto L20;
    }

/*  Check value of TENDNU */

    if (rkcom1_1.dir > 0. && *tendnu <= rkcom2_1.t) {
	ier = 911;
	nrec = 4;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A,D13.5/A,D13.5,A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** Integration is proceeding in the positive directi"
		"on. The ", (ftnlen)61);
	do_fio(&c__1, " ** current value for the independent variable is ", (
		ftnlen)50);
	do_fio(&c__1, (char *)&rkcom2_1.t, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " ** and you have set TENDNU = ", (ftnlen)30);
	do_fio(&c__1, (char *)&(*tendnu), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, ".  TENDNU must be ", (ftnlen)18);
	do_fio(&c__1, " ** greater than T.", (ftnlen)19);
	e_wsfi();
    } else if (rkcom1_1.dir < 0. && *tendnu >= rkcom2_1.t) {
	ier = 911;
	nrec = 4;
	ici__1.icierr = 0;
	ici__1.icirnum = 10;
	ici__1.icirlen = 80;
	ici__1.iciunit = rkcom9_1.rec;
	ici__1.icifmt = "(A/A,D13.5/A,D13.5,A/A)";
	s_wsfi(&ici__1);
	do_fio(&c__1, " ** Integration is proceeding in the negative directi"
		"on. The ", (ftnlen)61);
	do_fio(&c__1, " ** current value for the independent variable is ", (
		ftnlen)50);
	do_fio(&c__1, (char *)&rkcom2_1.t, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " ** and you have set TENDNU = ", (ftnlen)30);
	do_fio(&c__1, (char *)&(*tendnu), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, ".  TENDNU must be ", (ftnlen)18);
	do_fio(&c__1, " ** less than T.", (ftnlen)16);
	e_wsfi();
    } else {
/* Computing MAX */
/* Computing MAX */
	d__3 = abs(rkcom2_1.t), d__4 = abs(*tendnu);
	d__1 = rkcom7_1.tiny, d__2 = rkcom5_1.toosml * max(d__3,d__4);
	hmin = max(d__1,d__2);
	tdiff = (d__1 = *tendnu - rkcom2_1.t, abs(d__1));
	if (tdiff < hmin) {
	    ier = 911;
	    nrec = 4;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 10;
	    ici__1.icirlen = 80;
	    ici__1.iciunit = rkcom9_1.rec;
	    ici__1.icifmt = "(A,D13.5,A/A,D13.5,A/A/A,D13.5,A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, " ** The current value of the independent variable"
		    " T is ", (ftnlen)55);
	    do_fio(&c__1, (char *)&rkcom2_1.t, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ".", (ftnlen)1);
	    do_fio(&c__1, " ** The TENDNU you supplied has ABS(TENDNU-T) = ", 
		    (ftnlen)48);
	    do_fio(&c__1, (char *)&tdiff, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ".", (ftnlen)1);
	    do_fio(&c__1, " ** For the METHOD and the precision of the compu"
		    "ter being ", (ftnlen)59);
	    do_fio(&c__1, " ** used, this difference must be at least ", (
		    ftnlen)43);
	    do_fio(&c__1, (char *)&hmin, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, ".", (ftnlen)1);
	    e_wsfi();
	}
    }
    if (ier == 911) {
	goto L20;
    }

    rkcom1_1.tnd = *tendnu;
    rkcom2_1.last = FALSE_;

L20:

    rkmsg_(&ier, "RESET ", &nrec, &flag__, (ftnlen)6);

    return 0;
} /* reset_ */

/* Subroutine */ int mconst_(integer *method)
{
    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern /* Subroutine */ int envirn_(integer *, doublereal *, doublereal *)
	    ;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:   Sets machine-dependent global quantities */

/*  Common:    Initializes:    /RKCOM7/ OUTCH, MCHEPS, DWARF, RNDOFF, */
/*                                      SQRRMC, CUBRMC, TINY */
/*             Reads:          none */
/*             Alters:         none */

/*  Comments: */
/*  ========= */
/*  OUTCH, MCHEPS, DWARF are pure environmental parameters with values */
/*  obtained from a call to ENVIRN. The other quantities set depend on */
/*  the environmental parameters, the implementation, and, possibly, */
/*  METHOD. At present the METHODs implemented in the RK suite do not */
/*  influence the values of these quantities. */
/*  OUTCH  - Standard output channel */
/*  MCHEPS - Largest positive number such that 1.0D0 + MCHEPS = 1.0D0 */
/*  DWARF  - Smallest positive number */
/*  RNDOFF - 10 times MCHEPS */
/*  SQRRMC - Square root of MCHEPS */
/*  CUBRMC - Cube root of MCHEPS */
/*  TINY   - Square root of DWARF */

/*     .. Scalar Arguments .. */
/*     .. Common Block for Environment Parameters .. */
/*     .. Parameters .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    envirn_(&rkcom7_1.outch, &rkcom7_1.mcheps, &rkcom7_1.dwarf);

    rkcom7_1.rndoff = rkcom7_1.mcheps * 10.;
    rkcom7_1.sqrrmc = sqrt(rkcom7_1.mcheps);
    rkcom7_1.cubrmc = pow_dd(&rkcom7_1.mcheps, &c_b517);
    rkcom7_1.tiny = sqrt(rkcom7_1.dwarf);

    return 0;
} /* mconst_ */

/* Subroutine */ int envirn_(integer *outch, doublereal *mcheps, doublereal *
	dwarf)
{
    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___78 = { 0, 6, 0, 0, 0 };
    static cilist io___79 = { 0, 6, 0, 0, 0 };
    static cilist io___80 = { 0, 6, 0, 0, 0 };
    static cilist io___81 = { 0, 6, 0, 0, 0 };
    static cilist io___82 = { 0, 6, 0, 0, 0 };


/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

/*  The RK suite requires some environmental parameters that are provided by */
/*  this subroutine.  The values provided with the distribution codes are those */
/*  appropriate to the IEEE standard.  They must be altered, if necessary, to */
/*  those appropriate to the computing system you are using before calling the */
/*  codes of the suite. */

/*        ================================================================ */
/*        ================================================================ */
/*        TO MAKE SURE THAT THESE MACHINE AND INSTALLATION DEPENDENT */
/*        QUANTITIES ARE SPECIFIED PROPERLY, THE DISTRIBUTION VERSION */
/*        WRITES A MESSAGE ABOUT THE MATTER TO THE STANDARD OUTPUT CHANNEL */
/*        AND TERMINATES THE RUN.  THE VALUES PROVIDED IN THE DISTRIBUTION */
/*        VERSION SHOULD BE ALTERED, IF NECESSARY, AND THE "WRITE" AND */
/*        "STOP" STATEMENTS COMMENTED OUT. */
/*        ================================================================ */
/*        ================================================================ */

/*  OUTPUT VARIABLES */

/*     OUTCH     - INTEGER */
/*                 Standard output channel */
/*     MCHEPS    - DOUBLE PRECISION */
/*                 MCHEPS is the largest positive number such that */
/*                 1.0D0 + MCHEPS = 1.0D0. */
/*     DWARF     - DOUBLE PRECISION */
/*                 DWARF is the smallest positive number. */

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

/*     .. Scalar Arguments .. */
/*     .. Executable Statements .. */

/*  The following six statements are to be Commented out after verification that */
/*  the machine and installation dependent quantities are specified correctly. */
/*  If you pass copies of RKSUITE on to others, please give them the whole */
/*  distribution version of RKSUITE, and in particular, give them a version */
/*  of ENVIRN that does not have the following six statements Commented out. */
    s_wsle(&io___78);
    do_lio(&c__9, &c__1, " Before using RKSUITE, you must verify that the  ", 
	    (ftnlen)49);
    e_wsle();
    s_wsle(&io___79);
    do_lio(&c__9, &c__1, " machine- and installation-dependent quantities  ", 
	    (ftnlen)49);
    e_wsle();
    s_wsle(&io___80);
    do_lio(&c__9, &c__1, " specified in the subroutine ENVIRN are correct, ", 
	    (ftnlen)49);
    e_wsle();
    s_wsle(&io___81);
    do_lio(&c__9, &c__1, " and then Comment these WRITE statements and the ", 
	    (ftnlen)49);
    e_wsle();
    s_wsle(&io___82);
    do_lio(&c__9, &c__1, " STOP statement out of ENVIRN.                   ", 
	    (ftnlen)49);
    e_wsle();
    s_stop("", (ftnlen)0);

/*  The following values are appropriate to IEEE arithmetic with the typical */
/*  standard output channel. */

    *outch = 6;
    *mcheps = 1.11e-16;
    *dwarf = 2.23e-308;

/* ------------------------------------------------------------------------------ */
/*  If you have the routines D1MACH and I1MACH on your system, you could */
/*  replace the preceding statements by the following ones to obtain the */
/*  appropriate machine dependent numbers. The routines D1MACH and I1MACH */
/*  are public domain software.  They are available from NETLIB. */
/*      .. Scalar Arguments .. */
/*      INTEGER           OUTCH */
/*      DOUBLE PRECISION  DWARF, MCHEPS */
/*      .. External Functions .. */
/*      INTEGER           I1MACH */
/*      DOUBLE PRECISION  D1MACH */
/*      .. Executable Statements .. */

/*      OUTCH = I1MACH(2) */
/*      MCHEPS = D1MACH(3) */
/*      DWARF = D1MACH(1) */

/*  If you have the NAG Fortran Library available on your system, you could */
/*  replace the preceding statements by the following ones to obtain the */
/*  appropriate machine dependent numbers. */

/*      .. Scalar Arguments .. */
/*      INTEGER           OUTCH */
/*      DOUBLE PRECISION  DWARF, MCHEPS */
/*      .. External Functions .. */
/*      DOUBLE PRECISION  X02AJF, X02AMF */
/*      .. Executable Statements .. */

/*      CALL X04AAF(0,OUTCH) */
/*      MCHEPS = X02AJF() */
/*      DWARF = X02AMF() */

/*  If you have the IMSL MATH/LIBRARY available on your system, you could */
/*  replace the preceding statements by the following ones to obtain the */
/*  appropriate machine dependent numbers. */

/*      .. Scalar Arguments .. */
/*      INTEGER           OUTCH */
/*      DOUBLE PRECISION  DWARF, MCHEPS */
/*      .. External Functions .. */
/*      DOUBLE PRECISION  DMACH */
/*      .. Executable Statements .. */

/*      CALL UMACH(2,OUTCH) */
/*      MCHEPS = DMACH(4) */
/*      DWARF = DMACH(1) */
/* ------------------------------------------------------------------------------ */

    return 0;
} /* envirn_ */

/* Subroutine */ int const_(integer *method, integer *vecstg, logical *reqstg,
	 integer *lintpl)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    integer i__, j;
    doublereal diff, cdiff;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:   Set formula definitions and formula characteristics for */
/*             selected method. Return storage requirements for the */
/*             selected method. */

/*  Input:     METHOD */
/*  Output:    VECSTG, REQSTG, LINTPL */

/*  Common:    Initializes:    /RKCOM4/ A(*,*), B(*), C(*), BHAT(*), R(*), */
/*                                      E(*), PTR(*), NSTAGE, METHD, INTP, MINTP */
/*                             /RKCOM5/ TOOSML, COST, SAFETY, EXPON, STBRAD, */
/*                                      TANANG, RS, RS1, RS2, RS3, RS4, ORDER, */
/*                                      LSTSTG, MAXTRY, NSEC, FSAL */
/*             Reads:          /RKCOM7/ RNDOFF */
/*             Alters:         none */

/*  Comments: */
/*  ========= */
/*  Runge-Kutta formula pairs are described by a set of coefficients */
/*  and by the setting of a number of parameters that describe the */
/*  characteristics of the pair.  The input variable METHD indicates */
/*  which of the three pairs available is to be set. In the case of */
/*  METHD = 2 additional coefficients are defined that make interpolation */
/*  of the results possible and provide an additional error estimator. */
/*  VECSTG is the number of columns of workspace required to compute the */
/*  stages of a METHD. For interpolation purposes the routine returns via */
/*  COMMON the logical variable INTP indicating whether interpolation is */
/*  possible and via the call list: */
/*  REQSTG - whether the stages are required to form the */
/*           interpolant (set .FALSE. if INTP=.FALSE.) */
/*  LINTPL - the number of extra columns of storage required for use */
/*           with UT (set 0 if INTP=.FALSE.) */

/*  Quantities set in common blocks: */
/*  METHD - copy of METHOD */
/*  A, B, C, BHAT - coefficients of the selected method */
/*  R      - extra coefficents for interpolation with METHD = 2 */
/*  E      - extra coefficients for additional local error estimate */
/*           with METHD = 2 */
/*  PTR    - vector of pointers indicating how individual stages are to */
/*           be stored.  With it zero coefficients of the formulas can */
/*           be exploited to reduce the storage required */
/*  NSTAGE - number of stages for the specified METHD */
/*  INTP   - indicates whether there is an associated interpolant */
/*           (depending on the method, the user may have to supply */
/*           extra workspace) */
/*  MINTP  - the degree of the interpolating polynomial, if one exists */
/*  FSAL   - indicates whether the last stage of a step can be used as */
/*           the first stage of the following step */
/*  LSTSTG - pointer to location of last stage for use with FSAL=.TRUE. */
/*  ORDER  - the lower order of the pair of Runge-Kutta formulas that */
/*           constitute a METHD */
/*  TANANG, */
/*  STBRAD - the stability region of the formula used to advance */
/*           the integration is approximated by a sector in the left half */
/*           complex plane.  TANANG is the tangent of the interior angle */
/*           of the sector and STBRAD is the radius of the sector. */
/*  COST   - cost of a successful step in function evaluations */
/*  MAXTRY - limit on the number of iterations in the stiffness check. As */
/*           set, no more than 24 function evaluations are made in the check. */
/*  NSEC   - each step of size H in the primary integration corresponds to */
/*           NSEC steps of size H/NSEC in the secondary integration when */
/*           global error assessment is done. */
/*  EXPON  - used to adjust the step size; this code implements an error */
/*           per step control for which EXPON = 1/(ORDER + 1). */
/*  SAFETY - quantity used in selecting the step size */
/*  TOOSML - quantity used to determine when a step size is too small for */
/*           the precision available */
/*  RS, RS1, */
/*  RS2, RS3, */
/*  RS4    - quantities used in determining the maximum and minimum change */
/*           change in step size (set independently of METHD) */

/*  Further comments on SAFETY: */
/* ============================ */
/*  The code estimates the largest step size that will yield the specified */
/*  accuracy.  About half the time this step size would result in a local */
/*  error that is a little too large, and the step would be rejected.  For */
/*  this reason a SAFETY factor is used to reduce the "optimal" value to one */
/*  more likely to succeed.  Unfortunately, there is some art in choosing this */
/*  value. The more expensive a failed step, the smaller one is inclined to */
/*  take this factor. However, a small factor means that if the prediction were */
/*  good, more accuracy than desired would be obtained and the behavior of the */
/*  error would then be irregular.  The more stringent the tolerance, the better */
/*  the prediction, except near limiting precision. Thus the general range of */
/*  tolerances expected influences the choice of SAFETY. */

/*     .. Scalar Arguments .. */
/*     .. Common Block to hold Formula Definitions .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Environment Parameters .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    rkcom4_1.methd = *method;

    switch (rkcom4_1.methd) {
	case 1:  goto L20;
	case 2:  goto L40;
	case 3:  goto L100;
    }

/*  METHD = 1. */
/*    This pair is from "A 3(2) Pair of Runge-Kutta Formulas" by P. Bogacki */
/*    and L.F. Shampine, Appl. Math. Lett., 2, pp. 321-325, 1989.  The authors */
/*    are grateful to P. Bogacki for his assistance in implementing the pair. */

L20:
    rkcom4_1.nstage = 4;
    rkcom5_1.fsal = TRUE_;
    rkcom5_1.order = 2;
    rkcom5_1.tanang = 8.9;
    rkcom5_1.stbrad = 2.3;
    rkcom5_1.safety = .8;
    rkcom4_1.intp = TRUE_;
    rkcom4_1.mintp = 3;
    *reqstg = FALSE_;
    *lintpl = 2;
    rkcom5_1.nsec = 3;

    rkcom4_1.ptr[0] = 0;
    rkcom4_1.ptr[1] = 1;
    rkcom4_1.ptr[2] = 2;
    rkcom4_1.ptr[3] = 3;

    rkcom4_1.a[1] = .5;
    rkcom4_1.a[2] = 0.;
    rkcom4_1.a[15] = .75;
    rkcom4_1.a[3] = .22222222222222221;
    rkcom4_1.a[16] = .33333333333333331;
    rkcom4_1.a[29] = .44444444444444442;

/*  The coefficients BHAT(*) refer to the formula used to advance the */
/*  integration, here the one of order 3.  The coefficients B(*) refer */
/*  to the other formula, here the one of order 2. For this pair, BHAT(*) */
/*  is not needed since FSAL = .TRUE. */

    rkcom4_1.b[0] = .29166666666666669;
    rkcom4_1.b[1] = .25;
    rkcom4_1.b[2] = .33333333333333331;
    rkcom4_1.b[3] = .125;

    rkcom4_1.c__[0] = 0.;
    rkcom4_1.c__[1] = .5;
    rkcom4_1.c__[2] = .75;
    rkcom4_1.c__[3] = 1.;

    goto L120;

/*  METHD = 2 */
/*    This pair is from "An Efficient Runge-Kutta (4,5) Pair" by P. Bogacki */
/*    and L.F. Shampine, Rept. 89-20, Math. Dept., Southern Methodist */
/*    University, Dallas, Texas, USA, 1989.  The authors are grateful to */
/*    P. Bogacki for his assistance in implementing the pair.  Shampine and */
/*    Bogacki subsequently modified the formula to enhance the reliability of */
/*    the pair.  The original fourth order formula is used in an estimate of */
/*    the local error.  If the step fails, the computation is broken off.  If */
/*    the step is acceptable, the first evaluation of the next step is done, */
/*    i.e., the pair is implemented as FSAL and the local error of the step */
/*    is again estimated with a fourth order formula using the additional data. */
/*    The step must succeed with both estimators to be accepted.  When the */
/*    second estimate is formed, it is used for the subsequent adjustment of */
/*    the step size because it is of higher quality.  The two fourth order */
/*    formulas are well matched to leading order, and only exceptionally do */
/*    the estimators disagree -- problems with discontinuous coefficients are */
/*    handled more reliably by using two estimators as is global error */
/*    estimation. */

L40:
    rkcom4_1.nstage = 8;
    rkcom5_1.fsal = TRUE_;
    rkcom5_1.order = 4;
    rkcom5_1.tanang = 5.2;
    rkcom5_1.stbrad = 3.9;
    rkcom5_1.safety = .8;
    rkcom4_1.intp = TRUE_;
    *reqstg = TRUE_;
    rkcom4_1.mintp = 6;
    *lintpl = 6;
    rkcom5_1.nsec = 2;

    rkcom4_1.ptr[0] = 0;
    rkcom4_1.ptr[1] = 1;
    rkcom4_1.ptr[2] = 2;
    rkcom4_1.ptr[3] = 3;
    rkcom4_1.ptr[4] = 4;
    rkcom4_1.ptr[5] = 5;
    rkcom4_1.ptr[6] = 6;
    rkcom4_1.ptr[7] = 7;

    rkcom4_1.a[1] = .16666666666666666;
    rkcom4_1.a[2] = .07407407407407407;
    rkcom4_1.a[15] = .14814814814814814;
    rkcom4_1.a[3] = .13338192419825073;
    rkcom4_1.a[16] = -.47230320699708456;
    rkcom4_1.a[29] = .76749271137026243;
    rkcom4_1.a[4] = .22895622895622897;
    rkcom4_1.a[17] = -.36363636363636365;
    rkcom4_1.a[30] = .2937062937062937;
    rkcom4_1.a[43] = .50764050764050761;
    rkcom4_1.a[5] = .026500355113636364;
    rkcom4_1.a[18] = .23011363636363635;
    rkcom4_1.a[31] = .10772747760052448;
    rkcom4_1.a[44] = .16021907779720279;
    rkcom4_1.a[57] = .22543945312499999;
    rkcom4_1.a[6] = .18159821692916506;
    rkcom4_1.a[19] = -.38707982536247293;
    rkcom4_1.a[32] = .41288268560074054;
    rkcom4_1.a[45] = .64099131912534968;
    rkcom4_1.a[58] = -1.0121439383514517;
    rkcom4_1.a[71] = 1.1637515420586695;
    rkcom4_1.a[7] = .072792658730158735;
    rkcom4_1.a[20] = 0.;
    rkcom4_1.a[33] = .28662437773692473;
    rkcom4_1.a[46] = .19513621794871794;
    rkcom4_1.a[59] = .0086383928571428566;
    rkcom4_1.a[72] = .3595655806182122;
    rkcom4_1.a[85] = .077242772108843533;

/*  The coefficients B(*) refer to the formula of order 4. */

    rkcom4_1.b[0] = .07084476451760402;
    rkcom4_1.b[1] = 0.;
    rkcom4_1.b[2] = .29567307692307693;
    rkcom4_1.b[3] = .17965747482208388;
    rkcom4_1.b[4] = .029861111111111113;
    rkcom4_1.b[5] = .3462886755067825;
    rkcom4_1.b[6] = .071762401338705387;
    rkcom4_1.b[7] = .0059124957806361723;

/*  The coefficients E(*) refer to an estimate of the local error based on */
/*  the first formula of order 4.  It is the difference of the fifth order */
/*  result, here located in A(8,*), and the fourth order result.  By */
/*  construction both E(2) and E(7) are zero. */

    rkcom4_1.e[0] = -.0023437499999999999;
    rkcom4_1.e[1] = 0.;
    rkcom4_1.e[2] = .0103760754048583;
    rkcom4_1.e[3] = -.016490384615384615;
    rkcom4_1.e[4] = .018984375000000001;
    rkcom4_1.e[5] = -.010526315789473684;
    rkcom4_1.e[6] = 0.;

    rkcom4_1.c__[0] = 0.;
    rkcom4_1.c__[1] = .16666666666666666;
    rkcom4_1.c__[2] = .22222222222222221;
    rkcom4_1.c__[3] = .42857142857142855;
    rkcom4_1.c__[4] = .66666666666666663;
    rkcom4_1.c__[5] = .75;
    rkcom4_1.c__[6] = 1.;
    rkcom4_1.c__[7] = 1.;

/*  To do interpolation with this pair, some extra stages have to be computed. */
/*  The following additional A(*,*) and C(*) coefficients are for this purpose. */
/*  In addition there is an array R(*,*) that plays a role for interpolation */
/*  analogous to that of BHAT(*) for the basic step. */

    rkcom4_1.c__[8] = .5;
    rkcom4_1.c__[9] = .83333333333333337;
    rkcom4_1.c__[10] = .1111111111111111;

    rkcom4_1.a[8] = .074055989583333329;
    rkcom4_1.a[9] = -.063587240361623443;
    rkcom4_1.a[10] = .063609077240098705;
    rkcom4_1.a[21] = 0.;
    rkcom4_1.a[22] = .57424619248188691;
    rkcom4_1.a[23] = .010578541828541829;
    rkcom4_1.a[34] = .28964485093442743;
    rkcom4_1.a[35] = -.063650630072499534;
    rkcom4_1.a[36] = .066001009456705312;
    rkcom4_1.a[47] = .12839214966168092;
    rkcom4_1.a[48] = .043159777438314964;
    rkcom4_1.a[49] = .020483915553584021;
    rkcom4_1.a[60] = -.0037792968750000001;
    rkcom4_1.a[61] = .83701128838987326;
    rkcom4_1.a[62] = .0036822703302195488;
    rkcom4_1.a[73] = .014230019493177388;
    rkcom4_1.a[74] = -.34045447246719235;
    rkcom4_1.a[75] = .15525863227100201;
    rkcom4_1.a[86] = -.033793712797619051;
    rkcom4_1.a[87] = .049265038183349222;
    rkcom4_1.a[88] = -.085097025138180266;
    rkcom4_1.a[99] = .03125;
    rkcom4_1.a[100] = -.0068826776691659668;
    rkcom4_1.a[101] = .10000000000000001;
    rkcom4_1.a[113] = -.19577394258960973;
    rkcom4_1.a[114] = -.10000000000000001;
    rkcom4_1.a[127] = -.12340531043086005;

    for (i__ = 1; i__ <= 11; ++i__) {
	rkcom4_1.r__[i__ - 1] = 0.;
/* L60: */
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	rkcom4_1.r__[i__ * 11 - 10] = 0.;
/* L80: */
    }
    rkcom4_1.r__[55] = -11.547607240534195;
    rkcom4_1.r__[44] = -32.389845531126397;
    rkcom4_1.r__[33] = -34.205462244802526;
    rkcom4_1.r__[22] = -16.577402175822957;
    rkcom4_1.r__[11] = -3.2869708803427939;

    rkcom4_1.r__[57] = -27.245925284262768;
    rkcom4_1.r__[46] = -88.616760918474498;
    rkcom4_1.r__[35] = -106.67089628427365;
    rkcom4_1.r__[24] = -57.048459705648732;
    rkcom4_1.r__[13] = -12.035023433323735;

    rkcom4_1.r__[58] = -14.307769126529058;
    rkcom4_1.r__[47] = -47.606576610356406;
    rkcom4_1.r__[36] = -58.886200999259962;
    rkcom4_1.r__[25] = -32.574021109464375;
    rkcom4_1.r__[14] = -7.1817638119804785;

    rkcom4_1.r__[59] = -.72967792238625573;
    rkcom4_1.r__[48] = -2.3963551957301958;
    rkcom4_1.r__[37] = -3.000847537878788;
    rkcom4_1.r__[26] = -1.7486179638262971;
    rkcom4_1.r__[15] = -.42308609214859216;

    rkcom4_1.r__[60] = 19.121088881250603;
    rkcom4_1.r__[49] = 48.733692708914717;
    rkcom4_1.r__[38] = 34.32116842799104;
    rkcom4_1.r__[27] = -1.7935169069960781;
    rkcom4_1.r__[16] = -6.8616470879412184;

    rkcom4_1.r__[61] = -14.676126700680273;
    rkcom4_1.r__[50] = -45.882206632653059;
    rkcom4_1.r__[39] = -51.675414540816327;
    rkcom4_1.r__[28] = -24.563201530612243;
    rkcom4_1.r__[17] = -4.1711096938775514;

    rkcom4_1.r__[62] = 19.453125;
    rkcom4_1.r__[51] = 62.359375;
    rkcom4_1.r__[40] = 73.671875;
    rkcom4_1.r__[29] = 39.078125;
    rkcom4_1.r__[18] = 9.3125;
    rkcom4_1.r__[7] = 1.;

    rkcom4_1.r__[63] = 18.333333333333332;
    rkcom4_1.r__[52] = 71.;
    rkcom4_1.r__[41] = 103.;
    rkcom4_1.r__[30] = 66.333333333333329;
    rkcom4_1.r__[19] = 16.;

    rkcom4_1.r__[64] = -23.400440940191388;
    rkcom4_1.r__[53] = -70.201322820574163;
    rkcom4_1.r__[42] = -73.554221820959782;
    rkcom4_1.r__[31] = -30.106238940962648;
    rkcom4_1.r__[20] = -3.3528990003856314;

    rkcom4_1.r__[65] = 35.;
    rkcom4_1.r__[54] = 105.;
    rkcom4_1.r__[43] = 117.;
    rkcom4_1.r__[32] = 59.;
    rkcom4_1.r__[21] = 12.;

    goto L120;

/*  METHD = 3 */
/*    This pair is from "High Order Embedded Runge-Kutta Formulae" by P.J. */
/*    Prince and J.R. Dormand, J. Comp. Appl. Math.,7, pp. 67-75, 1981.  The */
/*    authors are grateful to P. Prince and J. Dormand for their assistance in */
/*    implementing the pair. */

L100:
    rkcom4_1.nstage = 13;
    rkcom5_1.fsal = FALSE_;
    rkcom5_1.order = 7;
    rkcom5_1.tanang = 11.;
    rkcom5_1.stbrad = 5.2;
    rkcom5_1.safety = .8;
    rkcom4_1.intp = FALSE_;
    *reqstg = FALSE_;
    rkcom4_1.mintp = 0;
    *lintpl = 0;
    rkcom5_1.nsec = 2;

    rkcom4_1.ptr[0] = 0;
    rkcom4_1.ptr[1] = 1;
    rkcom4_1.ptr[2] = 2;
    rkcom4_1.ptr[3] = 1;
    rkcom4_1.ptr[4] = 3;
    rkcom4_1.ptr[5] = 2;
    rkcom4_1.ptr[6] = 4;
    rkcom4_1.ptr[7] = 5;
    rkcom4_1.ptr[8] = 6;
    rkcom4_1.ptr[9] = 7;
    rkcom4_1.ptr[10] = 8;
    rkcom4_1.ptr[11] = 9;
    rkcom4_1.ptr[12] = 1;

    rkcom4_1.a[1] = .0555555555555555555555555555556;
    rkcom4_1.a[2] = .0208333333333333333333333333333;
    rkcom4_1.a[15] = .0625;
    rkcom4_1.a[3] = .03125;
    rkcom4_1.a[16] = 0.;
    rkcom4_1.a[29] = .09375;
    rkcom4_1.a[4] = .3125;
    rkcom4_1.a[17] = 0.;
    rkcom4_1.a[30] = -1.171875;
    rkcom4_1.a[43] = 1.171875;
    rkcom4_1.a[5] = .0375;
    rkcom4_1.a[18] = 0.;
    rkcom4_1.a[31] = 0.;
    rkcom4_1.a[44] = .1875;
    rkcom4_1.a[57] = .15;
    rkcom4_1.a[6] = .0479101371111111111111111111111;
    rkcom4_1.a[19] = 0.;
    rkcom4_1.a[32] = 0.;
    rkcom4_1.a[45] = .112248712777777777777777777778;
    rkcom4_1.a[58] = -.0255056737777777777777777777778;
    rkcom4_1.a[71] = .0128468238888888888888888888889;
    rkcom4_1.a[7] = .016917989787292281181431107136;
    rkcom4_1.a[20] = 0.;
    rkcom4_1.a[33] = 0.;
    rkcom4_1.a[46] = .387848278486043169526545744159;
    rkcom4_1.a[59] = .0359773698515003278967008896348;
    rkcom4_1.a[72] = .196970214215666060156715256072;
    rkcom4_1.a[85] = -.172713852340501838761392997002;
    rkcom4_1.a[8] = .0690957533591923006485645489846;
    rkcom4_1.a[21] = 0.;
    rkcom4_1.a[34] = 0.;
    rkcom4_1.a[47] = -.634247976728854151882807874972;
    rkcom4_1.a[60] = -.161197575224604080366876923982;
    rkcom4_1.a[73] = .138650309458825255419866950133;
    rkcom4_1.a[86] = .94092861403575626972423968413;
    rkcom4_1.a[99] = .211636326481943981855372117132;
    rkcom4_1.a[9] = .183556996839045385489806023537;
    rkcom4_1.a[22] = 0.;
    rkcom4_1.a[35] = 0.;
    rkcom4_1.a[48] = -2.46876808431559245274431575997;
    rkcom4_1.a[61] = -.291286887816300456388002572804;
    rkcom4_1.a[74] = -.026473020233117375688439799466;
    rkcom4_1.a[87] = 2.84783876419280044916451825422;
    rkcom4_1.a[100] = .281387331469849792539403641827;
    rkcom4_1.a[113] = .123744899863314657627030212664;
    rkcom4_1.a[10] = -1.21542481739588805916051052503;
    rkcom4_1.a[23] = 0.;
    rkcom4_1.a[36] = 0.;
    rkcom4_1.a[49] = 16.6726086659457724322804132886;
    rkcom4_1.a[62] = .915741828416817960595718650451;
    rkcom4_1.a[75] = -6.05660580435747094755450554309;
    rkcom4_1.a[88] = -16.0035735941561781118417064101;
    rkcom4_1.a[101] = 14.849303086297662557545391898;
    rkcom4_1.a[114] = -13.3715757352898493182930413962;
    rkcom4_1.a[127] = 5.13418264817963793317325361166;
    rkcom4_1.a[11] = .258860916438264283815730932232;
    rkcom4_1.a[24] = 0.;
    rkcom4_1.a[37] = 0.;
    rkcom4_1.a[50] = -4.77448578548920511231011750971;
    rkcom4_1.a[63] = -.43509301377703250944070041181;
    rkcom4_1.a[76] = -3.04948333207224150956051286631;
    rkcom4_1.a[89] = 5.57792003993609911742367663447;
    rkcom4_1.a[102] = 6.15583158986104009733868912669;
    rkcom4_1.a[115] = -5.06210458673693837007740643391;
    rkcom4_1.a[128] = 2.19392617318067906127491429047;
    rkcom4_1.a[141] = .134627998659334941535726237887;
    rkcom4_1.a[12] = .822427599626507477963168204773;
    rkcom4_1.a[25] = 0.;
    rkcom4_1.a[38] = 0.;
    rkcom4_1.a[51] = -11.6586732572776642839765530355;
    rkcom4_1.a[64] = -.757622116690936195881116154088;
    rkcom4_1.a[77] = .713973588159581527978269282765;
    rkcom4_1.a[90] = 12.0757749868900567395661704486;
    rkcom4_1.a[103] = -2.12765911392040265639082085897;
    rkcom4_1.a[116] = 1.99016620704895541832807169835;
    rkcom4_1.a[129] = -.234286471544040292660294691857;
    rkcom4_1.a[142] = .17589857770794226507310510589;
    rkcom4_1.a[155] = 0.;

/*  The coefficients BHAT(*) refer to the formula used to advance the */
/*  integration, here the one of order 8.  The coefficients B(*) refer */
/*  to the other formula, here the one of order 7. */

    rkcom4_1.bhat[0] = .0417474911415302462220859284685;
    rkcom4_1.bhat[1] = 0.;
    rkcom4_1.bhat[2] = 0.;
    rkcom4_1.bhat[3] = 0.;
    rkcom4_1.bhat[4] = 0.;
    rkcom4_1.bhat[5] = -.0554523286112393089615218946547;
    rkcom4_1.bhat[6] = .239312807201180097046747354249;
    rkcom4_1.bhat[7] = .70351066940344302305804641089;
    rkcom4_1.bhat[8] = -.759759613814460929884487677085;
    rkcom4_1.bhat[9] = .660563030922286341461378594838;
    rkcom4_1.bhat[10] = .158187482510123335529614838601;
    rkcom4_1.bhat[11] = -.238109538752862804471863555306;
    rkcom4_1.bhat[12] = .25;

    rkcom4_1.b[0] = .029553213676353496981964883112;
    rkcom4_1.b[1] = 0.;
    rkcom4_1.b[2] = 0.;
    rkcom4_1.b[3] = 0.;
    rkcom4_1.b[4] = 0.;
    rkcom4_1.b[5] = -.828606276487797039766805612689;
    rkcom4_1.b[6] = .311240900051118327929913751627;
    rkcom4_1.b[7] = 2.46734519059988698196468570407;
    rkcom4_1.b[8] = -2.54694165184190873912738007542;
    rkcom4_1.b[9] = 1.44354858367677524030187495069;
    rkcom4_1.b[10] = .0794155958811272872713019541622;
    rkcom4_1.b[11] = .0444444444444444444444444444445;
    rkcom4_1.b[12] = 0.;

    rkcom4_1.c__[0] = 0.;
    rkcom4_1.c__[1] = .0555555555555555555555555555556;
    rkcom4_1.c__[2] = .0833333333333333333333333333334;
    rkcom4_1.c__[3] = .125;
    rkcom4_1.c__[4] = .3125;
    rkcom4_1.c__[5] = .375;
    rkcom4_1.c__[6] = .1475;
    rkcom4_1.c__[7] = .465;
    rkcom4_1.c__[8] = .564865451382259575398358501426;
    rkcom4_1.c__[9] = .65;
    rkcom4_1.c__[10] = .924656277640504446745013574318;
    rkcom4_1.c__[11] = 1.;
    rkcom4_1.c__[12] = rkcom4_1.c__[11];

    goto L120;

/*  The definitions of all pairs come here for the calculation of */
/*  LSTSTG, RS1, RS2, RS3, RS4, COST, MAXTRY, EXPON, TOOSML, and VECSTG. */

L120:
    rkcom5_1.lststg = rkcom4_1.ptr[rkcom4_1.nstage - 1];
    if (rkcom5_1.fsal) {
	rkcom5_1.cost = (doublereal) (rkcom4_1.nstage - 1);
    } else {
	rkcom5_1.cost = (doublereal) rkcom4_1.nstage;
    }

/*  MAXTRY - limit on the number of iterations of a computation made in */
/*  diagnosing stiffness.  There are at most Q = 3 function calls per */
/*  iteration. MAXTRY is determined so that  Q*MAXTRY <= 5% of the cost of */
/*  50 steps and 1 <= MAXTRY <= 8. This limits the number of calls to FCN */
/*  in each diagnosis of stiffness to 24 calls. */

/* Computing MIN */
/* Computing MAX */
    i__3 = 1, i__4 = (integer) (rkcom5_1.cost * .05 * 50.);
    i__1 = 8, i__2 = max(i__3,i__4);
    rkcom5_1.maxtry = min(i__1,i__2);

    rkcom5_1.expon = 1. / (rkcom5_1.order + 1.);

/*     In calculating CDIFF it is assumed that there will be a non-zero */
/*     difference |C(I) - C(J)| less than one. If C(I) = C(J) for any I not */
/*     equal to J, they should be made precisely equal by assignment. */

    cdiff = 1.;
    i__1 = rkcom4_1.nstage - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = rkcom4_1.nstage;
	for (j = i__ + 1; j <= i__2; ++j) {
	    diff = (d__1 = rkcom4_1.c__[i__ - 1] - rkcom4_1.c__[j - 1], abs(
		    d__1));
	    if (diff != 0.) {
		cdiff = min(cdiff,diff);
	    }
/* L140: */
	}
/* L160: */
    }
    rkcom5_1.toosml = rkcom7_1.rndoff / cdiff;

/*  Determine the number of columns needed in STAGES(1:NEQ,*) (this will be */
/*  at most NSTAGE-1 since the first stage is held in a separate array). */
/*  The PTR array contains the column positions of the stages. */

    *vecstg = 0;
    i__1 = rkcom4_1.nstage;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MAX */
	i__2 = rkcom4_1.ptr[i__ - 1];
	*vecstg = max(i__2,*vecstg);
/* L180: */
    }

    rkcom5_1.rs = 2.;
    rkcom5_1.rs1 = 1. / rkcom5_1.rs;
/* Computing 2nd power */
    d__1 = rkcom5_1.rs;
    rkcom5_1.rs2 = d__1 * d__1;
/* Computing 3rd power */
    d__1 = rkcom5_1.rs;
    rkcom5_1.rs3 = d__1 * (d__1 * d__1);
    rkcom5_1.rs4 = 1. / rkcom5_1.rs3;

    return 0;
} /* const_ */

/* Subroutine */ int formi_(S_fp f, integer *neq, integer *nwant, doublereal *
	y, doublereal *yp, doublereal *yold, doublereal *ypold, doublereal *
	stages, logical *calstg, doublereal *xstage, doublereal *p)
{
    /* System generated locals */
    integer stages_dim1, stages_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer i__, j, k, l;
    doublereal d1, d2, d3, d4, hyp, hypold;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:    Forms an interpolating polynomial for use with */
/*              METHDs 1 or 2. */

/*  Input:      NEQ, NWANT, T, Y(*), YP(*), HOLD, YOLD(*), YPOLD(*), */
/*              STAGES(NEQ,*), CALSTG */
/*  Output:     P(*), XSTAGE(NEQ) */
/*  External:   F */

/*  Common:     Initializes:    none */
/*              Reads:          /RKCOM4/ A(*,*), C(*), R(*), METHD, MINTP */
/*                              /RKCOM2/ T, TOLD, HOLD */
/*              Alters:         /RKCOM2/ NFCN */

/*  Comments: */
/*  ========= */
/*  The integration has reached T with a step HOLD from TOLD = T-HOLD. */
/*  Y(*),YP(*) and YOLD(*),YPOLD(*) approximate the solution and its */
/*  derivative at T and TOLD respectively.  STAGES(NEQ,*) holds the stages */
/*  computed in taking this step. In the case of METHD = 2 it is necessary */
/*  to compute some more stages in this subroutine. CALSTG indicates whether */
/*  or not the extra stages need to be computed. A(*,*) and C(*) are used in */
/*  computing these stages. The extra stages are stored in STAGES(NEQ,*) and */
/*  XSTAGE(*).  The coefficients of the interpolating polynomials for the first */
/*  NWANT components of the solution are returned in the array P(*). The */
/*  polynomial is of degree MINTP = 3 for METHD = 1 and of degree MINTP = 6 */
/*  for METHD = 2. The vector R(*) is used for workspace when METHD = 2. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block to hold Formula Definitions .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    stages_dim1 = *neq;
    stages_offset = 1 + stages_dim1;
    stages -= stages_offset;
    --y;
    --yp;
    --yold;
    --ypold;
    --xstage;
    --p;

    /* Function Body */
    if (rkcom4_1.methd == 1) {

/*  METHD = 1.  Use the cubic Hermite interpolant that is is fully */
/*  specified by the values and slopes at the two ends of the step. */

	i__1 = *nwant;
	for (l = 1; l <= i__1; ++l) {
	    d1 = y[l] - yold[l];
	    hyp = rkcom2_1.hold * yp[l];
	    hypold = rkcom2_1.hold * ypold[l];
	    d2 = hyp - d1;
	    d3 = d1 - hypold;
	    d4 = d2 - d3;
	    p[l] = d2 + d4;
	    p[*nwant + l] = d4;
/* L20: */
	}

    } else {

/*  METHD = 2. */

	if (*calstg) {

/*  Compute the extra stages needed for interpolation using the facts that */
/*       1. Stage 1 is YPOLD(*). */
/*       2. Stage i (i>1) is stored in STAGES(1:NEQ,i). */
/*       3. This pair is FSAL, i.e. STAGES(1:NEQ,7)=YP(1:NEQ), which frees */
/*          up STAGES(1:NEQ,7) for use by stage 9. */
/*       4. XSTAGE(1:NEQ) is used for stage 10. */
/*       5. The coefficient of stage 2 in the interpolant is always 0, so */
/*          STAGES(1:NEQ,1) is used for stage 11. */
/*  The vector P(1:NEQ) is used as workspace for computing the stages. */

	    for (i__ = 9; i__ <= 11; ++i__) {
		i__1 = *neq;
		for (l = 1; l <= i__1; ++l) {
		    p[l] = rkcom4_1.a[i__ - 1] * ypold[l];
/* L40: */
		}
		i__1 = i__ - 1;
		for (j = 2; j <= i__1; ++j) {
		    if (j <= 7) {
			i__2 = *neq;
			for (l = 1; l <= i__2; ++l) {
			    p[l] += rkcom4_1.a[i__ + j * 13 - 14] * stages[l 
				    + (j - 1) * stages_dim1];
/* L60: */
			}
		    } else if (j == 8) {
			i__2 = *neq;
			for (l = 1; l <= i__2; ++l) {
			    p[l] += rkcom4_1.a[i__ + j * 13 - 14] * yp[l];
/* L80: */
			}
		    } else if (j == 9) {
			i__2 = *neq;
			for (l = 1; l <= i__2; ++l) {
			    p[l] += rkcom4_1.a[i__ + j * 13 - 14] * stages[l 
				    + stages_dim1 * 7];
/* L100: */
			}
		    } else if (j == 10) {
			i__2 = *neq;
			for (l = 1; l <= i__2; ++l) {
			    p[l] += rkcom4_1.a[i__ + j * 13 - 14] * xstage[l];
/* L120: */
			}
		    }
/* L140: */
		}
		i__1 = *neq;
		for (l = 1; l <= i__1; ++l) {
		    p[l] = yold[l] + rkcom2_1.hold * p[l];
/* L160: */
		}
		if (i__ == 9) {
		    d__1 = rkcom2_1.told + rkcom4_1.c__[i__ - 1] * 
			    rkcom2_1.hold;
		    (*f)(&d__1, &p[1], &stages[stages_dim1 * 7 + 1]);
		    ++rkcom2_1.nfcn;
		} else if (i__ == 10) {
		    d__1 = rkcom2_1.told + rkcom4_1.c__[i__ - 1] * 
			    rkcom2_1.hold;
		    (*f)(&d__1, &p[1], &xstage[1]);
		    ++rkcom2_1.nfcn;
		} else {
		    d__1 = rkcom2_1.told + rkcom4_1.c__[i__ - 1] * 
			    rkcom2_1.hold;
		    (*f)(&d__1, &p[1], &stages[stages_dim1 + 1]);
		    ++rkcom2_1.nfcn;
		}
/* L180: */
	    }
	}

/*  Form the coefficients of the interpolating polynomial in its shifted */
/*  and scaled form.  The transformation from the form in which the */
/*  polynomial is derived can be somewhat ill-conditioned.  The terms */
/*  are grouped so as to minimize the errors of the transformation. */

/*  Coefficient of SIGMA**6 */
	k = *nwant << 2;
	i__1 = *nwant;
	for (l = 1; l <= i__1; ++l) {
	    p[k + l] = rkcom4_1.r__[59] * stages[l + (stages_dim1 << 2)] + (
		    rkcom4_1.r__[64] * xstage[l] + rkcom4_1.r__[62] * yp[l] + 
		    (rkcom4_1.r__[61] * stages[l + stages_dim1 * 6] + 
		    rkcom4_1.r__[60] * stages[l + stages_dim1 * 5])) + (
		    rkcom4_1.r__[58] * stages[l + stages_dim1 * 3] + 
		    rkcom4_1.r__[63] * stages[l + stages_dim1 * 7] + (
		    rkcom4_1.r__[57] * stages[l + (stages_dim1 << 1)] + 
		    rkcom4_1.r__[65] * stages[l + stages_dim1]) + 
		    rkcom4_1.r__[55] * ypold[l]);
/* L200: */
	}

/*  Coefficient of SIGMA**5 */
	k = *nwant * 3;
	i__1 = *nwant;
	for (l = 1; l <= i__1; ++l) {
	    p[k + l] = rkcom4_1.r__[53] * xstage[l] + rkcom4_1.r__[52] * 
		    stages[l + stages_dim1 * 7] + (rkcom4_1.r__[50] * stages[
		    l + stages_dim1 * 6] + rkcom4_1.r__[49] * stages[l + 
		    stages_dim1 * 5] + rkcom4_1.r__[48] * stages[l + (
		    stages_dim1 << 2)]) + (rkcom4_1.r__[47] * stages[l + 
		    stages_dim1 * 3] + rkcom4_1.r__[51] * yp[l] + (
		    rkcom4_1.r__[46] * stages[l + (stages_dim1 << 1)] + 
		    rkcom4_1.r__[54] * stages[l + stages_dim1]) + 
		    rkcom4_1.r__[44] * ypold[l]);
/* L220: */
	}

/*  Coefficient of SIGMA**4 */
	k = *nwant << 1;
	i__1 = *nwant;
	for (l = 1; l <= i__1; ++l) {
	    p[k + l] = rkcom4_1.r__[36] * stages[l + stages_dim1 * 3] + 
		    rkcom4_1.r__[40] * yp[l] + (rkcom4_1.r__[39] * stages[l + 
		    stages_dim1 * 6] + rkcom4_1.r__[38] * stages[l + 
		    stages_dim1 * 5]) + rkcom4_1.r__[37] * stages[l + (
		    stages_dim1 << 2)] + (rkcom4_1.r__[42] * xstage[l] + 
		    rkcom4_1.r__[41] * stages[l + stages_dim1 * 7] + (
		    rkcom4_1.r__[35] * stages[l + (stages_dim1 << 1)] + 
		    rkcom4_1.r__[43] * stages[l + stages_dim1]) + 
		    rkcom4_1.r__[33] * ypold[l]);
/* L240: */
	}

/*  Coefficient of SIGMA**3 */
	k = *nwant;
	i__1 = *nwant;
	for (l = 1; l <= i__1; ++l) {
	    p[k + l] = rkcom4_1.r__[26] * stages[l + (stages_dim1 << 2)] + 
		    rkcom4_1.r__[27] * stages[l + stages_dim1 * 5] + (
		    rkcom4_1.r__[24] * stages[l + (stages_dim1 << 1)] + 
		    rkcom4_1.r__[30] * stages[l + stages_dim1 * 7] + (
		    rkcom4_1.r__[31] * xstage[l] + rkcom4_1.r__[29] * yp[l]) 
		    + rkcom4_1.r__[22] * ypold[l]) + (rkcom4_1.r__[25] * 
		    stages[l + stages_dim1 * 3] + rkcom4_1.r__[32] * stages[l 
		    + stages_dim1] + rkcom4_1.r__[28] * stages[l + 
		    stages_dim1 * 6]);
/* L260: */
	}

/*  Coefficient of SIGMA**2 */

	i__1 = *nwant;
	for (l = 1; l <= i__1; ++l) {
	    p[l] = rkcom4_1.r__[15] * stages[l + (stages_dim1 << 2)] + (
		    rkcom4_1.r__[16] * stages[l + stages_dim1 * 5] + 
		    rkcom4_1.r__[18] * yp[l] + rkcom4_1.r__[11] * ypold[l]) + 
		    (rkcom4_1.r__[13] * stages[l + (stages_dim1 << 1)] + 
		    rkcom4_1.r__[19] * stages[l + stages_dim1 * 7] + 
		    rkcom4_1.r__[20] * xstage[l]) + (rkcom4_1.r__[14] * 
		    stages[l + stages_dim1 * 3] + rkcom4_1.r__[21] * stages[l 
		    + stages_dim1] + rkcom4_1.r__[17] * stages[l + 
		    stages_dim1 * 6]);
/* L280: */
	}

/*  Scale all the coefficients by the step size. */

	i__1 = *nwant * (rkcom4_1.mintp - 1);
	for (l = 1; l <= i__1; ++l) {
	    p[l] = rkcom2_1.hold * p[l];
/* L300: */
	}

    }

    return 0;
} /* formi_ */

/* Subroutine */ int evali_(doublereal *y, doublereal *yp, doublereal *p, 
	doublereal *twant, char *reqest, integer *nwant, doublereal *ywant, 
	doublereal *ypwant, ftnlen reqest_len)
{
    /* System generated locals */
    integer p_dim1, p_offset, i__1;

    /* Local variables */
    integer k, l;
    doublereal sigma;
    char reqst1[1];

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:    Evaluation of an interpolating polynomial and/or its */
/*              first derivative. */

/*  Input:      Y(*), YP(*), P(NWANT,*), TWANT, REQEST, NWANT */
/*  Output:     YWANT(*), YPWANT(*) */

/*  Common:     Initializes:    none */
/*              Reads:          /RKCOM2/ HOLD, T */
/*                              /RKCOM4/ MINTP */
/*              Alters:         none */

/*  Comments: */
/*  ========= */
/*  The interpolant is evaluated at TWANT to approximate the solution, */
/*  YWANT, and/or its first derivative there, YPWANT. Only the first */
/*  NWANT components of the answer are computed. There are three cases */
/*  that are indicated by the first character of REQEST: */
/*    REQEST(1:1) = `S' or `s'- compute approximate `S'olution only. */
/*                = `D' or `d'- compute approximate first `D'erivative */
/*                              of the solution only. */
/*                = `B' or `b'- compute `B'oth approximate solution and */
/*                              first derivative. */
/*  The coefficents of the polynomial are contained in Y(*), YP(*) and */
/*  P(NWANT,*). */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block to hold Formula Definitions .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

/*  Evaluate the interpolating polynomial of degree MINTP in terms of the */
/*  shifted and scaled independent variable SIGMA. */

    /* Parameter adjustments */
    --y;
    --yp;
    p_dim1 = *nwant;
    p_offset = 1 + p_dim1;
    p -= p_offset;
    --ywant;
    --ypwant;

    /* Function Body */
    sigma = (*twant - rkcom2_1.t) / rkcom2_1.hold;

    *(unsigned char *)reqst1 = *(unsigned char *)reqest;
    if (*(unsigned char *)reqst1 == 'S' || *(unsigned char *)reqst1 == 's' || 
	    *(unsigned char *)reqst1 == 'B' || *(unsigned char *)reqst1 == 
	    'b') {

	i__1 = *nwant;
	for (l = 1; l <= i__1; ++l) {
	    ywant[l] = p[l + (rkcom4_1.mintp - 1) * p_dim1] * sigma;
/* L20: */
	}
	for (k = rkcom4_1.mintp - 2; k >= 1; --k) {
	    i__1 = *nwant;
	    for (l = 1; l <= i__1; ++l) {
		ywant[l] = (ywant[l] + p[l + k * p_dim1]) * sigma;
/* L40: */
	    }
/* L60: */
	}
	i__1 = *nwant;
	for (l = 1; l <= i__1; ++l) {
	    ywant[l] = (ywant[l] + rkcom2_1.hold * yp[l]) * sigma + y[l];
/* L80: */
	}
    }

/*  Evaluate the derivative of the interpolating polynomial. */

    if (*(unsigned char *)reqst1 == 'D' || *(unsigned char *)reqst1 == 'd' || 
	    *(unsigned char *)reqst1 == 'B' || *(unsigned char *)reqst1 == 
	    'b') {

/*  The derivative of the interpolating polynomial with respect to TWANT */
/*  is the derivative with respect to S divided by HOLD. */

	i__1 = *nwant;
	for (l = 1; l <= i__1; ++l) {
	    ypwant[l] = rkcom4_1.mintp * p[l + (rkcom4_1.mintp - 1) * p_dim1] 
		    * sigma;
/* L100: */
	}
	for (k = rkcom4_1.mintp - 1; k >= 2; --k) {
	    i__1 = *nwant;
	    for (l = 1; l <= i__1; ++l) {
		ypwant[l] = (ypwant[l] + k * p[l + (k - 1) * p_dim1]) * sigma;
/* L120: */
	    }
/* L140: */
	}
	i__1 = *nwant;
	for (l = 1; l <= i__1; ++l) {
	    ypwant[l] = (ypwant[l] + rkcom2_1.hold * yp[l]) / rkcom2_1.hold;
/* L160: */
	}
    }

    return 0;
} /* evali_ */

/* Subroutine */ int rkmsg_(integer *ier, char *srname, integer *nrec, 
	integer *flag__, ftnlen srname_len)
{
    /* System generated locals */
    integer i__1;
    cilist ci__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), do_fio(
	    integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen), s_copy(char *, char *, 
	    ftnlen, ftnlen);

    /* Local variables */
    integer i__;
    logical ok, on;
    extern /* Subroutine */ int chkfl_(logical *, logical *), rksit_(logical *
	    , char *, integer *, ftnlen);
    logical baderr, utcall;
    extern /* Subroutine */ int softfl_(logical *, logical *);

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:      To process error messages and terminate the program */
/*                in the event of a "catastrophic" failure. */

/*  Input:        IER, SRNAME, NREC */
/*  Output:       FLAG */

/*  Common:       Initializes:    none */
/*                Reads:          /RKCOM7/ OUTCH */
/*                                /RKCOM8/ MSG, UTASK */
/*                                /RKCOM9/ REC */
/*                Alters:         none */

/*  Comments: */
/*  ========= */
/*  The output variable FLAG is assigned the value of the input variable IER. */

/*  IER = -2  reports a successful call of the subroutine SRNAME and */
/*            indicates that special action is to be taken elsewhere */
/*            in the suite.  FLAG is set and a return is effected. */

/*  IER = 1   reports a successful call of the subroutine SRNAME.  FLAG */
/*            is set and a return is effected. */

/*  1 < IER < 911 and MSG = .TRUE.: a message of NREC records contained in */
/*            the array REC(*) is written to the standard output channel, */
/*            OUTCH.  FLAG is set and a return is effected. */

/*  IER = 911 reports a "catastrophic" error was detected in SRNAME.  A */
/*            message is written to OUTCH regardless of the value of MSG and */
/*            normally the execution of the program is terminated.  The */
/*            execution is not terminated if the error is the result of an */
/*            indirect call to CT, RESET, or INTRP through UT (UTASK = .TRUE.). */
/*            Termination can be prevented by using the subroutine SOFTFL. */

/*  IER = 912 reports that a "catastrophic" error was detected earlier and */
/*            termination was prevented, but the user has failed to take */
/*            appropriate remedial action.  Execution is terminated. */

/*     .. Scalar Arguments .. */
/*     .. Common Block for Environment Parameters .. */
/*     .. Common Block for Integrator Options .. */
/*     .. Common Block for Error Message .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

/*  Check where the call came from - if it is an indirect call from UT, */
/*  the run is not STOPped. */
    utcall = (s_cmp(srname, "RESET", srname_len, (ftnlen)5) == 0 || s_cmp(
	    srname, "CT", srname_len, (ftnlen)2) == 0 || s_cmp(srname, "INTRP"
	    , srname_len, (ftnlen)5) == 0) && rkcom8_1.utask;

/*  Check if can continue with integrator. */
    ok = (s_cmp(srname, "CT", srname_len, (ftnlen)2) == 0 || s_cmp(srname, 
	    "UT", srname_len, (ftnlen)2) == 0) && (*ier == 2 || *ier == 3 || *
	    ier == 4);

/*  Check if program termination has been overridden. */
    softfl_(&c_true, &on);

    if (rkcom8_1.msg && *ier > 1 || *ier >= 911) {
	ci__1.cierr = 0;
	ci__1.ciunit = rkcom7_1.outch;
	ci__1.cifmt = "(/A)";
	s_wsfe(&ci__1);
	do_fio(&c__1, " **", (ftnlen)3);
	e_wsfe();
	ci__1.cierr = 0;
	ci__1.ciunit = rkcom7_1.outch;
	ci__1.cifmt = "(A)";
	s_wsfe(&ci__1);
	i__1 = *nrec;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, rkcom9_1.rec + (i__ - 1) * 80, (ftnlen)80);
	}
	e_wsfe();
	if (*ier >= 911) {
	    ci__1.cierr = 0;
	    ci__1.ciunit = rkcom7_1.outch;
	    ci__1.cifmt = "(A/A,A,A/A/)";
	    s_wsfe(&ci__1);
	    do_fio(&c__1, " **", (ftnlen)3);
	    do_fio(&c__1, " ** Catastrophic error detected in ", (ftnlen)35);
	    do_fio(&c__1, srname, srname_len);
	    do_fio(&c__1, ".", (ftnlen)1);
	    do_fio(&c__1, " **", (ftnlen)3);
	    e_wsfe();
	    if (! (utcall || on) && *ier == 911 || *ier == 912) {
		ci__1.cierr = 0;
		ci__1.ciunit = rkcom7_1.outch;
		ci__1.cifmt = "(A/A/A)";
		s_wsfe(&ci__1);
		do_fio(&c__1, " **", (ftnlen)3);
		do_fio(&c__1, " ** Execution of your program is being termin"
			"ated.", (ftnlen)50);
		do_fio(&c__1, " **", (ftnlen)3);
		e_wsfe();
		s_stop("", (ftnlen)0);
	    }
	} else if (ok) {
	    ci__1.cierr = 0;
	    ci__1.ciunit = rkcom7_1.outch;
	    ci__1.cifmt = "(A/A,A,A,I2,A/A/A)";
	    s_wsfe(&ci__1);
	    do_fio(&c__1, " **", (ftnlen)3);
	    do_fio(&c__1, " ** Warning from routine ", (ftnlen)25);
	    do_fio(&c__1, srname, srname_len);
	    do_fio(&c__1, " with flag set ", (ftnlen)15);
	    do_fio(&c__1, (char *)&(*ier), (ftnlen)sizeof(integer));
	    do_fio(&c__1, ".", (ftnlen)1);
	    do_fio(&c__1, " ** You can continue integrating this problem.", (
		    ftnlen)46);
	    do_fio(&c__1, " **", (ftnlen)3);
	    e_wsfe();
	} else {
	    ci__1.cierr = 0;
	    ci__1.ciunit = rkcom7_1.outch;
	    ci__1.cifmt = "(A/A,A,A,I2,A/A/A)";
	    s_wsfe(&ci__1);
	    do_fio(&c__1, " **", (ftnlen)3);
	    do_fio(&c__1, " ** Warning from routine ", (ftnlen)25);
	    do_fio(&c__1, srname, srname_len);
	    do_fio(&c__1, " with flag set ", (ftnlen)15);
	    do_fio(&c__1, (char *)&(*ier), (ftnlen)sizeof(integer));
	    do_fio(&c__1, ".", (ftnlen)1);
	    do_fio(&c__1, " ** You cannot continue integrating this problem.",
		     (ftnlen)49);
	    do_fio(&c__1, " **", (ftnlen)3);
	    e_wsfe();
	}
    }
    for (i__ = *nrec + 1; i__ <= 10; ++i__) {
	s_copy(rkcom9_1.rec + (i__ - 1) * 80, " ", (ftnlen)80, (ftnlen)1);
/* L20: */
    }
    *flag__ = *ier;

/*  TELL RKSIT the status of the routine associated with SRNAME */
    rksit_(&c_false, srname, flag__, srname_len);

/*  Indicate that a catastrophic error has been detected */
    baderr = *flag__ >= 911;
    chkfl_(&c_false, &baderr);

    return 0;

} /* rkmsg_ */

/* Subroutine */ int rksit_(logical *ask, char *srname, integer *state, 
	ftnlen srname_len)
{
    /* Initialized data */

    static integer svsta[7] = { -1,-1,-1,-1,-1,-1,-1 };

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__, name__;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:      To save or enquire about the status of each */
/*                subprogram in the suite. */

/*  Input:        ASK, SRNAME */
/*  Input/output: STATE */


/*  Comments: */
/*  ========= */
/*  SRNAME indicates which routine is of interest in the call to RKSIT. */

/*  If ASK=.FALSE., then the value of STATE (which as used is the error */
/*  flag) for the routine SRNAME is saved internally. This value of STATE */
/*  is usually positive.  There are two exceptions: */
/*     1. SRNAME='SETUP' and STATE=-1 indicates a completely new problem, */
/*        so all SAVEd states are cleared. */
/*     2. STATE=-2 is used by some routines in the suite to indicate */
/*        circumstances which require special action. */

/*  If ASK=.TRUE., then RKSIT first checks to see if there were any */
/*  catastrophic errors, that is, a SAVEd state has value 911. This should */
/*  happen only when the user has overridden program termination in the event */
/*  of catastrophic failures from routines in the package but has failed to */
/*  take appropriate action. If this is the case, then RKSIT returns a value */
/*  of STATE = 911 which forces a termination of execution inside RKMSG. If */
/*  no catastrophic errors are flagged, then STATE returns the saved state */
/*  value for the routine specified by SRNAME. */

/*     .. Scalar Arguments .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. Save statement .. */
/*     .. Data statements .. */
/*     .. Executable Statements .. */

    if (s_cmp(srname, "SETUP", srname_len, (ftnlen)5) == 0) {
	name__ = 1;
    } else if (s_cmp(srname, "UT", srname_len, (ftnlen)2) == 0) {
	name__ = 2;
    } else if (s_cmp(srname, "STAT", srname_len, (ftnlen)4) == 0) {
	name__ = 3;
    } else if (s_cmp(srname, "GLBERR", srname_len, (ftnlen)6) == 0) {
	name__ = 4;
    } else if (s_cmp(srname, "CT", srname_len, (ftnlen)2) == 0) {
	name__ = 5;
    } else if (s_cmp(srname, "INTRP", srname_len, (ftnlen)5) == 0) {
	name__ = 6;
    } else if (s_cmp(srname, "RESET", srname_len, (ftnlen)5) == 0) {
	name__ = 7;
    } else {
	name__ = 0;
    }

/*  (Re)initialize if SETUP is telling RKSIT to do so. */
    if (! (*ask) && name__ == 1 && *state == -1) {
	for (i__ = 1; i__ <= 7; ++i__) {
	    svsta[i__ - 1] = -1;
/* L20: */
	}
	goto L60;
    }

/*  Check for 911 on exit from a previous call. */
    if (*ask) {
	for (i__ = 1; i__ <= 7; ++i__) {
	    if (svsta[i__ - 1] == 911) {
		*state = 911;
		goto L60;
	    }
/* L40: */
	}
    }

    if (*ask) {
	*state = svsta[name__ - 1];
    } else {
	svsta[name__ - 1] = *state;
    }

L60:

    return 0;
} /* rksit_ */

/* Subroutine */ int truerr_(S_fp f, integer *neq, doublereal *y, doublereal *
	tol, doublereal *weight, doublereal *zy, doublereal *zyp, doublereal *
	zerror, doublereal *zynew, doublereal *zerres, doublereal *zstage, 
	integer *ier)
{
    /* System generated locals */
    integer zstage_dim1, zstage_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    integer l;
    doublereal diff, hsec;
    logical main;
    doublereal hmin, tsec;
    extern /* Subroutine */ int step_(S_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, doublereal *, doublereal *, logical *);
    integer level, istep;
    doublereal zlerr, ztest1, ztest2, mxerlc, dumarr[1], errmax;
    logical ldummy;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:      Compute a running RMS measure of the true (global) error */
/*                for a general Runge-Kutta pair. */


/*  Input:        NEQ, Y(*), TOL, WEIGHT(*), */
/*  Input/output: ZY(*), ZYP(*), ZERROR(*) */
/*  Workspace:    ZYNEW(*), ZERRES(*), ZSTAGE(NEQ,*) */
/*  Output:       IER */
/*  External:     F */

/*  Common:       Initializes:    none */
/*                Reads:          /RKCOM2/ T, HOLD */
/*                                /RKCOM5/ TOOSML, ORDER, NSEC */
/*                                /RKCOM7/ TINY */
/*                Alters:         /RKCOM6/ MAXERR, LOCMAX, GNFCN */

/*  Comments: */
/*  ========= */
/*  A secondary integration is performed using a fraction of the step size */
/*  of the primary integration. ZY(*) and ZYP(*) are the approximate solution */
/*  and first derivative of this secondary integration. ZERRES(*) contains the */
/*  error estimates for the secondary integration. ZYNEW(*) and ZSTAGE(*,*) are */
/*  workspace for taking a step. The error assessment is computed using the */
/*  difference of the primary and secondary solutions at the primary */
/*  integration points as an estimate of the true error there.  The weights */
/*  used are those of the error test of the primary integration. This error */
/*  assessment is maintained in the vector ZERROR(*).  MAXERR and LOCMAX */
/*  contain the maximum contribution to the assessment and its location, */
/*  respectively.  The number of calls to F is counted by GNFCN. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Global Error Assessment .. */
/*     .. Common Block for Environment Parameters .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */
    /* Parameter adjustments */
    zstage_dim1 = *neq;
    zstage_offset = 1 + zstage_dim1;
    zstage -= zstage_offset;
    --y;
    --weight;
    --zy;
    --zyp;
    --zerror;
    --zynew;
    --zerres;

    /* Function Body */
    tsec = rkcom2_1.t - rkcom2_1.hold;
    hsec = rkcom2_1.hold / (doublereal) rkcom5_1.nsec;
/* Computing MAX */
/* Computing MAX */
    d__3 = abs(tsec), d__4 = abs(rkcom2_1.t);
    d__1 = rkcom7_1.tiny, d__2 = rkcom5_1.toosml * max(d__3,d__4);
    hmin = max(d__1,d__2);
    if (abs(hsec) < hmin) {
	*ier = 6;
	goto L120;
    }
    ztest1 = *tol / (doublereal) rkcom5_1.nsec;
    ztest2 = *tol / 10.;
    level = 0;

/*  The subroutine STEP is used to take a step.  In its use in the primary */
/*  integration provision is made for getting on scale in the first step. */
/*  In this situation only the subroutine might reduce the step size.  By */
/*  setting MAIN = .FALSE., the subroutine will take a step of the size input. */
/*  In this use of the subroutine, all items of the call list appearing after */
/*  MAIN are dummy variables. */

/*  Perform secondary integration. */
    main = FALSE_;
    ldummy = FALSE_;
    i__1 = rkcom5_1.nsec;
    for (istep = 1; istep <= i__1; ++istep) {

/*  Take a step. */
	step_((S_fp)f, neq, &tsec, &zy[1], &zyp[1], &zstage[zstage_offset], &
		ztest1, &hsec, &weight[1], &zynew[1], &zerres[1], &zlerr, &
		main, &c_b65, dumarr, &ldummy);

/*  The primary integration is using a step size of HUSED and the secondary */
/*  integration is using the smaller step size HSEC = HUSED/NSEC.  If steps */
/*  of this size were taken from the same starting point and the asymptotic */
/*  behavior were evident, the smaller step size would result in a local error */
/*  that is considerably smaller, namely by a factor of 1/(NSEC**(ORDER+1)). */
/*  If the two approximate solutions are close and TOLR is neither too large nor */
/*  too small, this should be approximately true.  The step size is chosen in */
/*  the primary integration so that the local error ERR is no larger than TOLR. */
/*  The local error, ZLERR, of the secondary integration is compared to TOLR in */
/*  an attempt to diagnose a secondary integration that is not rather more */
/*  accurate than the primary integration. */

	if (zlerr >= ztest1) {
	    level = 2;
	} else if (zlerr > ztest2) {
	    ++level;
	}
	if (level >= 2) {
	    *ier = 6;
	    goto L120;
	}

/*  Advance TSEC and the dependent variables ZY(*) and ZYP(*). */
	tsec = rkcom2_1.t - (doublereal) (rkcom5_1.nsec - istep) * hsec;
	i__2 = *neq;
	for (l = 1; l <= i__2; ++l) {
	    zy[l] = zynew[l];
/* L20: */
	}

	if (rkcom5_1.fsal) {

/*  When FSAL = .TRUE., the derivative ZYP(*) is the last stage of the step. */
	    i__2 = *neq;
	    for (l = 1; l <= i__2; ++l) {
		zyp[l] = zstage[l + rkcom5_1.lststg * zstage_dim1];
/* L40: */
	    }
	} else {

/*  Call F to evaluate ZYP(*). */
	    (*f)(&tsec, &zy[1], &zyp[1]);
	    ++rkcom6_1.gnfcn;
	}

/* L60: */
    }

/*  Update the maximum error seen, MAXERR, and its location, LOCMAX. */
/*  Use local variables ERRMAX and MXERLC. */

    errmax = rkcom6_1.maxerr;
    mxerlc = rkcom6_1.locmax;
    i__1 = *neq;
    for (l = 1; l <= i__1; ++l) {
	diff = (d__1 = zy[l] - y[l], abs(d__1)) / weight[l];
	if (diff > errmax) {
	    errmax = diff;
	    mxerlc = rkcom2_1.t;
	}
/* L80: */
    }

/*  If the global error is greater than 0.1D0, the solutions have diverged so */
/*  far that comparing them may not provide a reliable estimate of the global */
/*  error. The test is made before ZERROR(*) and MAXERR, LCMXER are updated so */
/*  that on a failure, they refer to the last reliable results. */

    if (errmax > .1) {
	*ier = 6;
	goto L120;
    } else {
	rkcom6_1.maxerr = errmax;
	rkcom6_1.locmax = mxerlc;
	i__1 = *neq;
	for (l = 1; l <= i__1; ++l) {
	    diff = (d__1 = zy[l] - y[l], abs(d__1)) / weight[l];
/* Computing 2nd power */
	    d__1 = diff;
	    zerror[l] += d__1 * d__1;
/* L100: */
	}
	*ier = 1;
    }

/*  Exit point for TRUERR */
L120:

    return 0;
} /* truerr_ */

/* Subroutine */ int step_(S_fp f, integer *neq, doublereal *tnow, doublereal 
	*y, doublereal *yp, doublereal *stages, doublereal *tol, doublereal *
	htry, doublereal *weight, doublereal *ynew, doublereal *errest, 
	doublereal *err, logical *main, doublereal *hmin, doublereal *thres, 
	logical *phase2)
{
    /* System generated locals */
    integer stages_dim1, stages_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    integer i__, j, l;
    doublereal avgy, tstg;
    extern /* Subroutine */ int stepa_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *), stepb_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, doublereal *);
    logical cutbak;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:      To compute a step of an explicit Runge-Kutta */
/*                method and estimate the local error of the step. */

/*  Input:        NEQ, TNOW, Y(*), YP(*), TOL, MAIN, HMIN, THRES(*) */
/*  Input/output: HTRY, PHASE2, LAST, WEIGHT(*) */
/*  Output:       STAGES(NEQ,*), YNEW(*), ERREST(*), ERR */

/*  Common:       Initializes:    none */
/*                Reads:          /RKCOM1/ TND */
/*                                /RKCOM2/ LAST */
/*                                /RKCOM4/ A, B, C, BHAT, PTR, NSTAGE, METHD */
/*                                /RKCOM5/ FSAL */
/*                Alters:         /RKCOM2/ NFCN, LAST */
/*                                /RKCOM6/ GNFCN */

/*  Comments: */
/*  ========= */
/*  From an approximate solution Y(*) at TNOW and first derivative there, */
/*  YP(*) = F(TNOW,Y,YP), a step is taken to get an approximation YNEW(*) */
/*  at TNOW + HTRY. The Runge-Kutta method and how it is used are defined */
/*  by A, B, C, BHAT, PTR, NSTAGE, METHD and FSAL. Intermediate stages */
/*  of the method are stored in the array STAGES(NEQ,*). The error in */
/*  each solution component is estimated and returned in ERREST(*). A */
/*  weighted maximum norm of the local error, ERR, is formed. For some */
/*  methods an intermediate error estimate can be computed before completion */
/*  of the step (see routine STEPB); if the estimate is greater than the */
/*  specified tolerance TOL, the computation of the step is terminated. */

/*  When global error estimation is desired, two integrations are done. */
/*  The usual integration is referred to as the "primary", or "main", */
/*  integration (MAIN=.TRUE.).  For global error estimation another, */
/*  "secondary" integration (MAIN=.FALSE.) is carried out with a smaller */
/*  step size.  The weight vector WEIGHT(*) used in computing ERR is */
/*  determined by the main integration.  Thus this argument is output when */
/*  MAIN = .TRUE. and input when MAIN = .FALSE.. */

/*  When taking the first step in an integration, the logical variable */
/*  PHASE2 may be input as .TRUE. and if the first step is the whole of */
/*  the range of integration, then LAST will be .TRUE.. When PHASE2=.TRUE., */
/*  the first three stages are monitored to help assure that the step */
/*  size H is small enough for the integration to be stable and for the */
/*  estimate of the error of the step to be credible. Calls are made to */
/*  the subroutine STEPA for this purpose. If necessary, H will be */
/*  reduced in STEPA (and LAST altered accordingly) and the step retried */
/*  in STEP until an acceptable value is found. */

/*  In the primary integration the number of calls to F is counted by */
/*  NFCN, and in the secondary integration, by GNFCN. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block to hold Formula Definitions .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Global Error Assessment .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*  Many of the following loops over L = 1, NEQ have constant array values */
/*  inside. The code is written with clarity in mind.  Any optimizing */
/*  compiler will identify these occurrences and take appropriate action. */
/*  A check for zero multipliers has been included so as to prevent */
/*  needless computation resulting from the storing of zero coefficients */
/*  in the arrays for the sake of clarity.  The array ERREST(*) is used */
/*  for working storage in this computation. */

    /* Parameter adjustments */
    stages_dim1 = *neq;
    stages_offset = 1 + stages_dim1;
    stages -= stages_offset;
    --y;
    --yp;
    --weight;
    --ynew;
    --errest;
    --thres;

    /* Function Body */
L20:
    if (*main) {
	if (*phase2) {

/*  Initialize weights for measuring the local error. */
	    i__1 = *neq;
	    for (l = 1; l <= i__1; ++l) {
/* Computing MAX */
		d__2 = thres[l], d__3 = (d__1 = y[l], abs(d__1));
		weight[l] = max(d__2,d__3);
/* L40: */
	    }
	}
    }

    i__1 = rkcom4_1.nstage;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    if (j == 1) {
		i__3 = *neq;
		for (l = 1; l <= i__3; ++l) {
		    errest[l] = rkcom4_1.a[i__ - 1] * yp[l];
/* L60: */
		}
	    } else {
		if (rkcom4_1.a[i__ + j * 13 - 14] != 0.) {
		    i__3 = *neq;
		    for (l = 1; l <= i__3; ++l) {
			errest[l] += rkcom4_1.a[i__ + j * 13 - 14] * stages[l 
				+ rkcom4_1.ptr[j - 1] * stages_dim1];
/* L80: */
		    }
		}
	    }
/* L100: */
	}
	i__2 = *neq;
	for (l = 1; l <= i__2; ++l) {
	    ynew[l] = y[l] + *htry * errest[l];
/* L120: */
	}

/*  METHD = 2 is special in that an estimate of the local error can be */
/*  formed before the step is completed.  If the step is a failure, */
/*  return immediately.  Otherwise, complete the step and compute a more */
/*  accurate error estimate. */
	if (rkcom4_1.methd == 2 && i__ == 7) {
	    stepb_(neq, &y[1], &yp[1], htry, &ynew[1], &stages[stages_offset],
		     &thres[1], err, main, &weight[1]);
	    if (*err > *tol) {
		return 0;
	    }
	}

	tstg = *tnow + rkcom4_1.c__[i__ - 1] * *htry;
	if (*main && rkcom2_1.last && rkcom4_1.c__[i__ - 1] == 1.) {
	    tstg = rkcom1_1.tnd;
	}
	(*f)(&tstg, &ynew[1], &stages[rkcom4_1.ptr[i__ - 1] * stages_dim1 + 1]
		);

/*  Increment the counter for the number of function evaluations */
/*  depending on whether the primary or secondary integration is taking */
/*  place. */
	if (*main) {
	    ++rkcom2_1.nfcn;
	} else {
	    ++rkcom6_1.gnfcn;
	}

/* ---------------------------------------------------------------------- */
/*  When PHASE2 is .TRUE. we are in the second phase of the automatic */
/*  selection of the initial step size.  The results of the first three */
/*  stages are monitored in the subroutine STEPA for evidence that H is */
/*  too large -- instability and/or an unreliable estimate of the error */
/*  of the step is then possible.  When the subroutine believes H to be */
/*  too large, it returns CUTBAK = .TRUE. and a suitably reduced H for */
/*  another try. */

	if (*main) {
	    if (*phase2) {
		if (i__ <= 3 && abs(*htry) > *hmin) {
		    stepa_(tnow, &y[1], &yp[1], &tstg, &ynew[1], &stages[
			    rkcom4_1.ptr[i__ - 1] * stages_dim1 + 1], htry, &
			    weight[1], &cutbak);
		    if (cutbak) {
			rkcom2_1.last = FALSE_;

/*  Make sure that STEPA does not reduce the step size below the */
/*  minimum. If it does, reset H to HMIN and deactivate PHASE2. */
			if (abs(*htry) <= *hmin) {
			    *htry = d_sign(hmin, htry);
			    *phase2 = FALSE_;
			}
			goto L20;
		    }
		}
	    }
	}
/* ---------------------------------------------------------------------- */

/* L140: */
    }

/*  Some formulas are constructed so that the last stage represents */
/*  the result of the step (FSAL=.TRUE.), hence if the step is acceptable, */
/*  it will be the first stage for the next step. When FSAL=.FALSE., we */
/*  have to complete the computation of the step. */

    if (! rkcom5_1.fsal) {
	i__1 = rkcom4_1.nstage;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (i__ == 1) {
		i__2 = *neq;
		for (l = 1; l <= i__2; ++l) {
		    errest[l] = rkcom4_1.bhat[0] * yp[l];
/* L160: */
		}
	    } else {
		if (rkcom4_1.bhat[i__ - 1] != 0.) {
		    i__2 = *neq;
		    for (l = 1; l <= i__2; ++l) {
			errest[l] += rkcom4_1.bhat[i__ - 1] * stages[l + 
				rkcom4_1.ptr[i__ - 1] * stages_dim1];
/* L180: */
		    }
		}
	    }
/* L200: */
	}
	i__1 = *neq;
	for (l = 1; l <= i__1; ++l) {
	    ynew[l] = y[l] + *htry * errest[l];
/* L220: */
	}
    }

/*  Form an estimate of the error in the lower order formula by comparing */
/*  it to the higher order formula of the pair. ERREST(*) has been used */
/*  as working storage above.  The higher order approximation has been */
/*  formed as YNEW(*) = Y(*) + HTRY*ERREST(*) where ERREST(*) is a linear */
/*  combination of the stages of the formula. The lower order result also */
/*  has the form Y(*) plus HTRY times a different linear combination of */
/*  the stages. Hence, this different linear combination of stages for */
/*  the lower order formula can just be subtracted from the combination */
/*  stored in ERREST(*) to produce the errors. The result is then */
/*  multiplied by HTRY to obtain the error estimate. */

    i__1 = rkcom4_1.nstage;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == 1 && rkcom4_1.b[0] != 0.) {
	    i__2 = *neq;
	    for (l = 1; l <= i__2; ++l) {
		errest[l] -= rkcom4_1.b[0] * yp[l];
/* L240: */
	    }
	} else {
	    if (rkcom4_1.b[i__ - 1] != 0.) {
		i__2 = *neq;
		for (l = 1; l <= i__2; ++l) {
		    errest[l] -= rkcom4_1.b[i__ - 1] * stages[l + 
			    rkcom4_1.ptr[i__ - 1] * stages_dim1];
/* L260: */
		}
	    }
	}
/* L280: */
    }
    i__1 = *neq;
    for (l = 1; l <= i__1; ++l) {
	errest[l] = *htry * errest[l];
/* L300: */
    }

/*  The error in a solution component is measured relative to a weight */
/*  that is the larger of a threshold and the size of the solution over */
/*  the step.  Using the magnitude of a solution component at both ends */
/*  of the step in the definition of "size" increases the robustness of */
/*  the test. When global error estimation is specified, the weight */
/*  vector WEIGHT(*) is defined by the primary integration and is then */
/*  used in the secondary integration. */

    if (*main) {
	i__1 = *neq;
	for (l = 1; l <= i__1; ++l) {
	    avgy = ((d__1 = y[l], abs(d__1)) + (d__2 = ynew[l], abs(d__2))) * 
		    .5;
/* Computing MAX */
	    d__1 = avgy, d__2 = thres[l];
	    weight[l] = max(d__1,d__2);
/* L320: */
	}
    }

    *err = 0.;
    i__1 = *neq;
    for (l = 1; l <= i__1; ++l) {
/* Computing MAX */
	d__2 = *err, d__3 = (d__1 = errest[l] / weight[l], abs(d__1));
	*err = max(d__2,d__3);
/* L340: */
    }

    return 0;
} /* step_ */

/* Subroutine */ int stepa_(doublereal *tnow, doublereal *y, doublereal *yp, 
	doublereal *tstg, doublereal *ystg, doublereal *ypstg, doublereal *
	htry, doublereal *weight, logical *cutbak)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer l;
    doublereal wt, scl, twt, ynrm, fdiff, tdiff, argdif, ystgnm;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:      To calculate an "on-scale" step size for phase 2 of */
/*                the initial step size computation. */

/*  Input:        TNOW, Y(*), YP(*), TSTG, YSTG(*), YPSTG(*) */
/*  Input/output: HTRY, WEIGHT */
/*  Output:       CUTBAK */

/*  Common:       Initializes:    none */
/*                Reads:          /RKCOM1/ TND, NEQ */
/*                                /RKCOM5/ STBRAD, RS1, RS4 */
/*                                /RKCOM7/ RNDOFF */
/*                Alters:         none */

/*  Comments: */
/*  ========= */
/*  This subroutine is used during the first three stages of the first step. */
/*  A Lipschitz constant L for the differential equation in autonomous form */
/*  is approximated, and the product abs(HTRY)*L is compared to an approximate */
/*  radius, STBRAD, of the stability region of the method. The step size is */
/*  reduced as necessary, within a range specified by the step size control */
/*  parameters RS1 and RS4, to assure stability and give some confidence in */
/*  the error estimator.  If HTRY is reduced, CUTBAK is set .TRUE.. */

/*  Y(*) and YP(*) contain the solution and its derivative at TNOW and */
/*  similarly YSTG(*) and YPSTG(*) contain approximations at TSTG. */

/*  Normally the weights used in the control of the error depend on the */
/*  size of the solution at the beginning and at the end of the step, but */
/*  at this time we do not have a solution at the end of the step.  Each */
/*  stage YSTG(*) of the Runge - Kutta process represents a low order */
/*  approximation to the solution at TSTG.  Because the initial value of */
/*  WEIGHT(*) provided in the first phase of the scheme is based only on */
/*  the solution at T and THRES(*), it is continually updated in STEPA to */
/*  account for the size of the solution throughout the step as revealed */
/*  by the intermediate stages YSTG(*). Inside this subroutine only, the */
/*  differential equation is converted to autonomous form. After the */
/*  conversion, the end of the interval of integration, TND, is used */
/*  to define a suitable weight for the independent variable. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Environment Parameters .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*  Update the weights to account for the current intermediate solution */
/*  approximation YSTG(*).  Compute the sizes of Y(*) and YSTG(*) in the */
/*  new norm.  The size of the Lipschitz constant is assessed by a difference */
/*  in the arguments Y(*), YSTG(*) and a difference in the function evaluated */
/*  at these arguments. */

    /* Parameter adjustments */
    --weight;
    --ypstg;
    --ystg;
    --yp;
    --y;

    /* Function Body */
    ynrm = 0.;
    ystgnm = 0.;
    argdif = 0.;
    fdiff = 0.;
    i__1 = rkcom1_1.neqn;
    for (l = 1; l <= i__1; ++l) {
/* Computing MAX */
	d__2 = weight[l], d__3 = (d__1 = ystg[l], abs(d__1));
	wt = max(d__2,d__3);
	weight[l] = wt;
/* Computing MAX */
	d__2 = ynrm, d__3 = (d__1 = y[l], abs(d__1)) / wt;
	ynrm = max(d__2,d__3);
/* Computing MAX */
	d__2 = ystgnm, d__3 = (d__1 = ystg[l], abs(d__1)) / wt;
	ystgnm = max(d__2,d__3);
/* Computing MAX */
	d__2 = argdif, d__3 = (d__1 = ystg[l] - y[l], abs(d__1)) / wt;
	argdif = max(d__2,d__3);
/* Computing MAX */
	d__2 = fdiff, d__3 = (d__1 = ypstg[l] - yp[l], abs(d__1)) / wt;
	fdiff = max(d__2,d__3);
/* L20: */
    }

/*  The transformation of the equation to autonomous form is done */
/*  implicitly.  The difference of the arguments must take into account */
/*  the difference between the values of the independent variable T and */
/*  TSTG. The difference of the corresponding component of the function */
/*  is zero because of the way the standard transformation is done. */

    tdiff = *tstg - *tnow;
    twt = (d__1 = rkcom1_1.tnd - *tnow, abs(d__1));
/* Computing MAX */
    d__1 = ynrm, d__2 = abs(*tnow) / twt;
    ynrm = max(d__1,d__2);
/* Computing MAX */
    d__1 = ystgnm, d__2 = abs(*tstg) / twt;
    ystgnm = max(d__1,d__2);
/* Computing MAX */
    d__1 = argdif, d__2 = abs(tdiff) / twt;
    argdif = max(d__1,d__2);

/*  The ratio FDIFF/ARGDIF is a lower bound for, and an approximation to, a */
/*  Lipschitz constant L for the differential equation written in autonomous */
/*  form.  First we must ask if the difference ARGDIF is significant in the */
/*  precision available.  If it appears to be, we insist that abs(HTRY)*L be */
/*  less than an approximate radius, STBRAD, of the stability region of the */
/*  method.  This is more stringent than necessary for stability, possibly a */
/*  lot more stringent, but the aim is to get an HTRY small enough that the */
/*  error estimate for the step is credible.  The reduction is required to be */
/*  at least as much as the step control parameter RS1. It is necessary to */
/*  limit the reduction of HTRY at any one time because we may be misled in */
/*  the size of the reduction that is appropriate due to nonlinearity of the */
/*  differential equation and to inaccurate weights caused by HTRY much too */
/*  large.  The reduction is not permitted to be more than the step control */
/*  parameter RS4. */

    *cutbak = FALSE_;
    if (argdif > rkcom7_1.rndoff * max(ynrm,ystgnm)) {
	if (abs(*htry) * fdiff > rkcom5_1.stbrad * argdif) {
	    scl = rkcom5_1.stbrad * argdif / (abs(*htry) * fdiff);
	    scl = min(scl,rkcom5_1.rs1);
	    scl = max(scl,rkcom5_1.rs4);
	    *htry = scl * *htry;
	    *cutbak = TRUE_;
	}
    }

    return 0;
} /* stepa_ */

/* Subroutine */ int stepb_(integer *neq, doublereal *y, doublereal *yp, 
	doublereal *h__, doublereal *ynew, doublereal *stages, doublereal *
	thres, doublereal *err, logical *main, doublereal *weight)
{
    /* System generated locals */
    integer stages_dim1, stages_offset, i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer l;
    doublereal wt, sum, avgy;
    integer index;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:      To compute an error estimate for METHD = 2 prior */
/*                to completing the step. */

/*  Input:        NEQ, Y(*), YP(*), H, STAGES(NEQ,*), THRES(*), MAIN, */
/*                WEIGHT(*) */
/*  Output:       ERR */

/*  Common:       Initializes:    none */
/*                Reads:          /RKCOM4/ E, PTR */
/*                Alters:         none */

/*  Comments: */
/*  ========= */
/*  If global error assessment is taking place, then MAIN = .FALSE. and */
/*  the weight vector generated by the primary integration is used.  The */
/*  error estimate is a linear combination (with coefficients in E(*)) */
/*  of the stages stored in STAGES(*,*) (located by PTR(*)). */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Common Block to hold Formula Definitions .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    stages_dim1 = *neq;
    stages_offset = 1 + stages_dim1;
    stages -= stages_offset;
    --y;
    --yp;
    --ynew;
    --thres;
    --weight;

    /* Function Body */
    *err = 0.;
    i__1 = *neq;
    for (l = 1; l <= i__1; ++l) {

/*  Estimate the local error of component L. The coding makes use of */
/*  E(2) = 0.0D0 and E(7) = 0.0D0. */

	sum = rkcom4_1.e[0] * yp[l];
	for (index = 3; index <= 6; ++index) {
	    sum += rkcom4_1.e[index - 1] * stages[l + rkcom4_1.ptr[index - 1] 
		    * stages_dim1];
/* L20: */
	}

/*  The local error is H*SUM.  A weighted maximum norm of SUM is formed */
/*  and then the factor of H is taken into account. */

	if (*main) {
	    avgy = ((d__1 = y[l], abs(d__1)) + (d__2 = ynew[l], abs(d__2))) * 
		    .5;
/* Computing MAX */
	    d__1 = avgy, d__2 = thres[l];
	    wt = max(d__1,d__2);
	} else {
	    wt = weight[l];
	}

/* Computing MAX */
	d__2 = *err, d__3 = (d__1 = sum / wt, abs(d__1));
	*err = max(d__2,d__3);
/* L40: */
    }
    *err = abs(*h__) * *err;

    return 0;
} /* stepb_ */

/* Subroutine */ int stiff_(S_fp f, doublereal *havg, integer *jflstp, 
	logical *toomch, integer *maxfcn, doublereal *work, integer *ier, 
	integer *nrec)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    icilist ici__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;

    /* Local variables */
    integer l;
    logical stif;
    doublereal avgy;
    extern /* Subroutine */ int stiffa_(S_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, logical *, logical *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    logical lotsfl;
    doublereal xtrawk;
    logical unsure;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:      Diagnose stiffness.  This depends on two things: whether */
/*                the step size is being restricted on grounds of stability */
/*                and whether the integration to TND can be completed in no */
/*                more than MAXFCN function evaluations. */

/*  Input:        HAVG, TOOMCH, MAXFCN, WORK(*) */
/*  Input/output: JFLSTP */
/*  Output:       IER, NREC */
/*  Workspace:    WORK(*) */
/*  External:     F */

/*  Common:       Initializes:    /RKCOM9/ REC */
/*                Reads:          /RKCOM1/ TND, NEQN */
/*                                /RKCOM2/ T, H, NFCN, SVNFCN, OKSTP */
/*                                /RKCOM3/ PRY, PRYP, PRTHRS, PRWT, PRSCR, */
/*                                         PRSTGS, PRYOLD */
/*                                /RKCOM5/ COST */
/*                Alters:         /RKCOM2/ NFCN */
/*                                /RKCOM9/ REC */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Common Block for General Workspace Pointers .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Error Message .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --work;

    /* Function Body */
    if ((rkcom2_1.okstp - 10) % 40 == 0) {
	lotsfl = *jflstp >= 10;
	*jflstp = 0;
    } else {
	lotsfl = FALSE_;
    }

/*  If either too much work has been done or there are lots of failed steps, */
/*  test for stiffness. */

    if (*toomch || lotsfl) {

/*  Regenerate weight vector */
	i__1 = rkcom1_1.neqn;
	for (l = 1; l <= i__1; ++l) {
	    avgy = ((d__1 = work[rkcom3_1.pry - 1 + l], abs(d__1)) + (d__2 = 
		    work[rkcom3_1.pryold - 1 + l], abs(d__2))) * .5;
/* Computing MAX */
	    d__1 = avgy, d__2 = work[rkcom3_1.prthrs - 1 + l];
	    work[rkcom3_1.prwt - 1 + l] = max(d__1,d__2);
/* L20: */
	}

/*  STIFFA determines whether the problem is STIFF. In some circumstances it */
/*  is UNSURE.  The decision depends on two things: whether the step size is */
/*  being restricted on grounds of stability and whether the integration to */
/*  TND can be completed in no more than MAXFCN function evaluations.  The */
/*  last four arguments of STIFFA are vectors of length NEQN used for working */
/*  storage.  Some storage in WORK(*) reserved for the stages (there are a */
/*  minimum of three such vectors reserved for the METHDs implemented) and */
/*  the scratch vector starting at PRSCR are used for this purpose. */

	stiffa_((S_fp)f, &rkcom2_1.t, &work[rkcom3_1.pry], &rkcom2_1.h__, 
		havg, &rkcom1_1.tnd, maxfcn, &work[rkcom3_1.prwt], &work[
		rkcom3_1.pryp], &work[rkcom3_1.prerst], &unsure, &stif, &work[
		rkcom3_1.prstgs], &work[rkcom3_1.prstgs + rkcom1_1.neqn], &
		work[rkcom3_1.prstgs + (rkcom1_1.neqn << 1)], &work[
		rkcom3_1.prscr]);
	if (! unsure) {
	    if (stif) {

/*  Predict how much eXTRA WorK will be needed to reach TND. */
		xtrawk = rkcom5_1.cost * (d__1 = (rkcom1_1.tnd - rkcom2_1.t) /
			 *havg, abs(d__1)) / (doublereal) (rkcom2_1.svnfcn + 
			rkcom2_1.nfcn);
		*ier = 4;
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 80;
		ici__1.iciunit = rkcom9_1.rec + *nrec * 80;
		ici__1.icifmt = "(A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, " ** Your problem has been diagnosed as stiff."
			"  If the ", (ftnlen)54);
		e_wsfi();
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 80;
		ici__1.iciunit = rkcom9_1.rec + (*nrec + 1) * 80;
		ici__1.icifmt = "(A,D13.5)";
		s_wsfi(&ici__1);
		do_fio(&c__1, " ** situation persists, it will cost roughly ",
			 (ftnlen)45);
		do_fio(&c__1, (char *)&xtrawk, (ftnlen)sizeof(doublereal));
		e_wsfi();
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 80;
		ici__1.iciunit = rkcom9_1.rec + (*nrec + 2) * 80;
		ici__1.icifmt = "(A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, " ** times as much to reach TEND as it has cos"
			"t to reach TNOW.", (ftnlen)61);
		e_wsfi();
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 80;
		ici__1.iciunit = rkcom9_1.rec + (*nrec + 3) * 80;
		ici__1.icifmt = "(A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, " ** You should probably change to a code inte"
			"nded for ", (ftnlen)54);
		e_wsfi();
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 80;
		ici__1.iciunit = rkcom9_1.rec + (*nrec + 4) * 80;
		ici__1.icifmt = "(A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, " ** stiff problems. ", (ftnlen)20);
		e_wsfi();
		*nrec += 5;
	    }
	}
    }

    return 0;
} /* stiff_ */

/* Subroutine */ int stiffa_(S_fp f, doublereal *x, doublereal *y, doublereal 
	*hnow, doublereal *havg, doublereal *xend, integer *maxfcn, 
	doublereal *wt, doublereal *fxy, doublereal *v0, logical *unsure, 
	logical *stif, doublereal *v1, doublereal *v2, doublereal *v3, 
	doublereal *vtemp)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer l;
    doublereal d1, d2, r1[2], r2[2], rho, v0v0, v0v1, v0v2, v1v1, v1v2, v1v3, 
	    v2v2, v2v3, v3v3, det1, det2, rho2, res2, rold, dist, ynrm, beta1,
	     beta2;
    integer ntry;
    doublereal v0nrm, root1[2], v3nrm, root2[2], scale, alpha1, alpha2;
    extern /* Subroutine */ int stiffb_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *), stiffc_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), stiffd_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, S_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal dotprd_(doublereal *, doublereal *, doublereal *, 
	    integer *);
    doublereal xtrfcn;
    logical rootre;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  External:     F */
/*  Input:        X, Y(*), HNOW, HAVG, XEND, MAXFCN, WT(*), FXY(*) */
/*  Input/Output  V0(*) */
/*  Output:       UNSURE, STIF */
/*  Workspace:    V1(*), V2(*), V3(*), VTEMP(*) */

/*  Common:       Initializes:    none */
/*                Reads:          /RKCOM1/ TND, NEQN */
/*                                /RKCOM5/ COST, STBRAD, TANANG */
/*                                /RKCOM7/ SQRRMC, CUBRMC */
/*                Alters:         none */

/*  STIFFA diagnoses stiffness for an explicit Runge-Kutta code.  When it */
/*  is called, either many step failures have been observed, or a lot of */
/*  work has been done. */

/*  The NEQ equations of the problem are defined by the subroutine F(X,Y,YP). */
/*  When STIFFA is called, the integration has reached X where the approximate */
/*  solution is Y(*).  The vector FXY(*) is defined by a call of F(X,Y,FXY). */
/*  It is an input argument because it is usually available from the integrator. */

/*  The last successful step was of size HNOW, and an average step size is */
/*  HAVG.  A weighted norm is used to measure the local error with the error */
/*  in solution component L divided by the positive weight WT(L) provided in */
/*  the vector WT(*). */

/*  Explicit Runge - Kutta codes estimate the local error of Y(*) by */
/*  forming the difference of two approximate solutions.  This difference */
/*  must be provided in the vector V0(*).  When this difference is too */
/*  small to be significant, STIFFA will replace it with a "random" vector. */

/*  STIF is set .TRUE. when the average step size appears to be restricted */
/*  on grounds of stability.  In certain cases the variable UNSURE is set */
/*  .TRUE.; the value of STIF is then not defined. */

/*  The stability region of the explicit Runge-Kutta formula is described */
/*  by quantities TANANG and STBRAD that are communicated by the setup routine */
/*  via COMMON.  Stability regions often change sharply near the imaginary */
/*  axis so that it is difficult to classify the stiffness of a problem with */
/*  eigenvalues of a local Jacobian that are "near" the imaginary axis.  For */
/*  this reason,  we consider only points Z in the upper left half complex */
/*  plane for which TAN( IMAG(Z)/( - RE(Z))) <= TANANG. Eigenvalues outside */
/*  this region are one reason for the code being UNSURE.  The stability */
/*  region is approximated by the intersection of a disk with this sector. */
/*  The radius of this disk is called STBRAD. */

/*  Working storage must be provided via the four vectors V1(*),V2(*), */
/*  V3(*),VTEMP(*).  These vectors must be of length at least NEQ. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block to hold Formula Characterisitcs .. */
/*     .. Common Block for Environment Parameters .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*  If the current step size differs substantially from the average, */
/*  the problem is not stiff. */

    /* Parameter adjustments */
    --vtemp;
    --v3;
    --v2;
    --v1;
    --v0;
    --fxy;
    --wt;
    --y;

    /* Function Body */
    if ((d__1 = *hnow / *havg, abs(d__1)) > 5. || (d__2 = *hnow / *havg, abs(
	    d__2)) < .2) {
	*stif = FALSE_;
	*unsure = FALSE_;
	return 0;
    } else {
	*unsure = TRUE_;
    }

/*  The average step size is used to predict the cost in function evaluations */
/*  of finishing the integration to XEND.  If this cost is no more than MAXFCN, */
/*  the problem is declared not stiff: If the step size is being restricted on */
/*  grounds of stability, it will stay close to HAVG.  The prediction will */
/*  then be good, but the cost is too low to consider the problem stiff.  If */
/*  the step size is not close to HAVG, the problem is not stiff.  Either way */
/*  there is no point to testing for a step size restriction due to stability. */

    xtrfcn = rkcom5_1.cost * (d__1 = (*xend - *x) / *havg, abs(d__1));
    if (xtrfcn <= (doublereal) (*maxfcn)) {
	*stif = FALSE_;
	*unsure = FALSE_;
	return 0;
    } else {
	*unsure = TRUE_;
    }

/*  There have been many step failures or a lot of work has been done.  Now */
/*  we must determine if this is due to the stability characteristics of the */
/*  formula.  This is done by calculating the dominant eigenvalues of the */
/*  local Jacobian and then testing whether HAVG corresponds to being on the */
/*  boundary of the stability region. */

/*  The size of Y(*) provides scale information needed to approximate */
/*  the Jacobian by differences. */

    ynrm = sqrt(dotprd_(&y[1], &y[1], &wt[1], &rkcom1_1.neqn));
    scale = ynrm * rkcom7_1.sqrrmc;
    if (scale == 0.) {

/*  Degenerate case.  Y(*) is (almost) the zero vector so the scale is not */
/*  defined.  The input vector V0(*) is the difference between Y(*) and a */
/*  lower order approximation to the solution that is within the error */
/*  tolerance.  When Y(*) vanishes, V0(*) is itself an acceptable approximate */
/*  solution, so we take SCALE from it, if this is possible. */

	ynrm = sqrt(dotprd_(&v0[1], &v0[1], &wt[1], &rkcom1_1.neqn));
	scale = ynrm * rkcom7_1.sqrrmc;
	if (scale == 0.) {
	    *unsure = TRUE_;
	    return 0;
	}
    }

    v0v0 = dotprd_(&v0[1], &v0[1], &wt[1], &rkcom1_1.neqn);
    if (v0v0 == 0.) {

/*  Degenerate case.  V0(*) is (almost) the zero vector so cannot */
/*  be used to define a direction for an increment to Y(*).  Try a */
/*  "random" direction. */

	i__1 = rkcom1_1.neqn;
	for (l = 1; l <= i__1; ++l) {
	    v0[l] = 1.;
/* L20: */
	}
	v0v0 = dotprd_(&v0[1], &v0[1], &wt[1], &rkcom1_1.neqn);
    }
    v0nrm = sqrt(v0v0);
    i__1 = rkcom1_1.neqn;
    for (l = 1; l <= i__1; ++l) {
	v0[l] /= v0nrm;
/* L40: */
    }
    v0v0 = 1.;

/*  Use a nonlinear power method to estimate the two dominant eigenvalues. */
/*  V0(*) is often very rich in the two associated eigenvectors.  For this */
/*  reason the computation is organized with the expectation that a minimal */
/*  number of iterations will suffice.  Indeed, it is necessary to recognize */
/*  a kind of degeneracy when there is a dominant real eigenvalue.  The */
/*  subroutine STIFFB does this.  In the first try, NTRY = 1, a Rayleigh */
/*  quotient for such an eigenvalue is initialized as ROLD.  After each */
/*  iteration, REROOT computes a new Rayleigh quotient and tests whether the */
/*  two approximations agree to one tenth of one per cent and the eigenvalue, */
/*  eigenvector pair satisfy a stringent test on the residual.  ROOTRE = .TRUE. */
/*  signals that a single dominant real root has been found. */

    ntry = 1;
L60:

    stiffd_(&v0[1], havg, x, &y[1], (S_fp)f, &fxy[1], &wt[1], &scale, &v0v0, &
	    v1[1], &v1v1, &vtemp[1]);

/*  The quantity SQRT(V1V1/V0V0) is a lower bound for the product of HAVG */
/*  and a Lipschitz constant.  If it should be LARGE, stiffness is not */
/*  restricting the step size to the stability region.  The principle is */
/*  clear enough, but the real reason for this test is to recognize an */
/*  extremely inaccurate computation of V1V1 due to finite precision */
/*  arithmetic in certain degenerate circumstances. */

    if (sqrt(v1v1) > sqrt(v0v0) * 1e10) {
	*unsure = TRUE_;
	return 0;
    }

    v0v1 = dotprd_(&v0[1], &v1[1], &wt[1], &rkcom1_1.neqn);
    if (ntry == 1) {
	rold = v0v1 / v0v0;

/*  This is the first Rayleigh quotient approximating the product of HAVG */
/*  and a dominant real eigenvalue.  If it should be very small, the */
/*  problem is not stiff.  It is important to test for this possibility so */
/*  as to prevent underflow and degeneracies in the subsequent iteration. */

	if (abs(rold) < rkcom7_1.cubrmc) {
	    *unsure = FALSE_;
	    *stif = FALSE_;
	    return 0;
	}
    } else {
	stiffb_(&v1v1, &v0v1, &v0v0, &rold, &rho, root1, root2, &rootre);
	if (rootre) {
	    goto L100;
	}
    }
    stiffd_(&v1[1], havg, x, &y[1], (S_fp)f, &fxy[1], &wt[1], &scale, &v1v1, &
	    v2[1], &v2v2, &vtemp[1]);
    v0v2 = dotprd_(&v0[1], &v2[1], &wt[1], &rkcom1_1.neqn);
    v1v2 = dotprd_(&v1[1], &v2[1], &wt[1], &rkcom1_1.neqn);
    stiffb_(&v2v2, &v1v2, &v1v1, &rold, &rho, root1, root2, &rootre);
    if (rootre) {
	goto L100;
    }

/*  Fit a quadratic in the eigenvalue to the three successive iterates */
/*  V0(*),V1(*),V2(*) of the power method to get a first approximation to */
/*  a pair of eigenvalues.  A test made earlier in STIFFB implies that */
/*  the quantity DET1 here will not be too small. */

/* Computing 2nd power */
    d__1 = v0v1;
    det1 = v0v0 * v1v1 - d__1 * d__1;
    alpha1 = (-v0v0 * v1v2 + v0v1 * v0v2) / det1;
    beta1 = (v0v1 * v1v2 - v1v1 * v0v2) / det1;

/*  Iterate again to get V3, test again for degeneracy, and then fit a */
/*  quadratic to V1(*),V2(*),V3(*) to get a second approximation to a pair */
/*  of eigenvalues. */

    stiffd_(&v2[1], havg, x, &y[1], (S_fp)f, &fxy[1], &wt[1], &scale, &v2v2, &
	    v3[1], &v3v3, &vtemp[1]);
    v1v3 = dotprd_(&v1[1], &v3[1], &wt[1], &rkcom1_1.neqn);
    v2v3 = dotprd_(&v2[1], &v3[1], &wt[1], &rkcom1_1.neqn);
    stiffb_(&v3v3, &v2v3, &v2v2, &rold, &rho, root1, root2, &rootre);
    if (rootre) {
	goto L100;
    }
/* Computing 2nd power */
    d__1 = v1v2;
    det2 = v1v1 * v2v2 - d__1 * d__1;
    alpha2 = (-v1v1 * v2v3 + v1v2 * v1v3) / det2;
    beta2 = (v1v2 * v2v3 - v2v2 * v1v3) / det2;

/*  First test the residual of the quadratic fit to see if we might */
/*  have determined a pair of eigenvalues. */

/* Computing 2nd power */
    d__2 = alpha2;
/* Computing 2nd power */
    d__3 = beta2;
    res2 = (d__1 = v3v3 + v2v2 * (d__2 * d__2) + v1v1 * (d__3 * d__3) + v2v3 *
	     2. * alpha2 + v1v3 * 2. * beta2 + v1v2 * 2. * alpha2 * beta2, 
	    abs(d__1));
    if (res2 <= v3v3 * 9.9999999999999995e-7) {

/*  Calculate the two approximate pairs of eigenvalues. */

	stiffc_(&alpha1, &beta1, r1, r2);
	stiffc_(&alpha2, &beta2, root1, root2);

/*  The test for convergence is done on the larger root of the second */
/*  approximation.  It is complicated by the fact that one pair of roots */
/*  might be real and the other complex.  First calculate the spectral */
/*  radius RHO of HAVG*J as the magnitude of ROOT1.  Then see if one of */
/*  the roots R1,R2 is within one per cent of ROOT1.  A subdominant root */
/*  may be very poorly approximated if its magnitude is much smaller than */
/*  RHO -- this does not matter in our use of these eigenvalues. */

/* Computing 2nd power */
	d__1 = root1[0];
/* Computing 2nd power */
	d__2 = root1[1];
	rho = sqrt(d__1 * d__1 + d__2 * d__2);
/* Computing 2nd power */
	d__1 = root1[0] - r1[0];
/* Computing 2nd power */
	d__2 = root1[1] - r1[1];
	d1 = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
	d__1 = root1[0] - r2[0];
/* Computing 2nd power */
	d__2 = root1[1] - r2[1];
	d2 = d__1 * d__1 + d__2 * d__2;
	dist = sqrt((min(d1,d2)));
	if (dist <= rho * .001) {
	    goto L100;
	}
    }

/*  Do not have convergence yet.  Because the iterations are cheap, and */
/*  because the convergence criterion is stringent, we are willing to try */
/*  a few iterations. */

    if (ntry < rkcom5_1.maxtry) {
	++ntry;
	v3nrm = sqrt(v3v3);
	i__1 = rkcom1_1.neqn;
	for (l = 1; l <= i__1; ++l) {
	    v0[l] = v3[l] / v3nrm;
/* L80: */
	}
	v0v0 = 1.;
	goto L60;
    } else {
	*unsure = TRUE_;
	return 0;
    }

/*                        ************** */

/*  We now have the dominant eigenvalues.  Decide if the average step */
/*  size is being restricted on grounds of stability.  Check the real */
/*  parts of the eigenvalues.  First see if the dominant eigenvalue is */
/*  in the left half plane -- there won't be a stability restriction */
/*  unless it is. If there is another eigenvalue of comparable magnitude */
/*  with a positive real part, the problem is not stiff. If the dominant */
/*  eigenvalue is too close to the imaginary axis, we cannot diagnose */
/*  stiffness. */

L100:
    if (root1[0] > 0.) {
	*stif = FALSE_;
	*unsure = FALSE_;
	return 0;
    }
/* Computing 2nd power */
    d__1 = root2[0];
/* Computing 2nd power */
    d__2 = root2[1];
    rho2 = sqrt(d__1 * d__1 + d__2 * d__2);
    if (rho2 >= rho * .9 && root2[0] > 0.) {
	*stif = FALSE_;
	*unsure = FALSE_;
	return 0;
    }
    if (abs(root1[1]) > abs(root1[0]) * rkcom5_1.tanang) {
	*unsure = TRUE_;
	return 0;
    }

/*  If the average step size corresponds to being well within the */
/*  stability region, the step size is not being restricted because */
/*  of stability. */

    *stif = rho >= rkcom5_1.stbrad * .9;
    *unsure = FALSE_;
    return 0;
} /* stiffa_ */

/* Subroutine */ int stiffb_(doublereal *v1v1, doublereal *v0v1, doublereal *
	v0v0, doublereal *rold, doublereal *rho, doublereal *root1, 
	doublereal *root2, logical *rootre)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    doublereal r__, det, res;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Input:        V1V1, V0V1, V0V0 */
/*  Input/output: ROLD */
/*  Output:       RHO, ROOT1(*),ROOT2(*),ROOTRE */

/*  Decide if the iteration has degenerated because of a strongly */
/*  dominant real eigenvalue.  Have just computed the latest iterate. */
/*  V1V1 is its dot product with itself, V0V1 is the dot product */
/*  of the previous iterate with the current one, and V0V0 is the */
/*  dot product of the previous iterate with itself.  ROLD is a */
/*  previous Rayleigh quotient approximating a dominant real */
/*  eigenvalue.  It must be computed directly the first time the */
/*  subroutine is called.  It is updated each call to STIFFB, hence */
/*  is available for subsequent calls. */

/*  If there is a strongly dominant real eigenvalue, ROOTRE is set */
/*  .TRUE., ROOT1(*) returns the eigenvalue, RHO returns the magnitude */
/*  of the eigenvalue, and ROOT2(*) is set to zero. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --root2;
    --root1;

    /* Function Body */
    r__ = *v0v1 / *v0v0;
    *rho = abs(r__);
/* Computing 2nd power */
    d__1 = *v0v1;
    det = *v0v0 * *v1v1 - d__1 * d__1;
    res = (d__1 = det / *v0v0, abs(d__1));
    *rootre = det == 0. || res <= *v1v1 * 9.9999999999999995e-7 && (d__1 = 
	    r__ - *rold, abs(d__1)) <= *rho * .001;
    if (*rootre) {
	root1[1] = r__;
	root1[2] = 0.;
	root2[1] = 0.;
	root2[2] = 0.;
    }
    *rold = r__;

    return 0;
} /* stiffb_ */

/* Subroutine */ int stiffc_(doublereal *alpha, doublereal *beta, doublereal *
	r1, doublereal *r2)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal disc, temp, sqdisc;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Input:  ALPHA, BETA */
/*  Output: R1(*), R2(*) */

/*  This subroutine computes the two complex roots R1 and R2 of */
/*  the quadratic equation X**2 + ALPHA*X + BETA = 0.  The magnitude */
/*  of R1 is greater than or equal to the magnitude of R2. R1 and R2 are */
/*  returned as vectors of two components with the first being the real */
/*  part of the complex number and the second being the imaginary part. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */
    /* Parameter adjustments */
    --r2;
    --r1;

    /* Function Body */
    temp = *alpha / 2.;
/* Computing 2nd power */
    d__1 = temp;
    disc = d__1 * d__1 - *beta;
    if (disc == 0.) {

/*  Double root. */

	r1[1] = -temp;
	r1[2] = 0.;
	r2[1] = r1[1];
	r2[2] = r1[2];
	return 0;
    }

    sqdisc = sqrt((abs(disc)));
    if (disc < 0.) {

/*  Complex conjugate roots. */

	r1[1] = -temp;
	r1[2] = sqdisc;
	r2[1] = r1[1];
	r2[2] = -r1[2];
    } else {

/*  Real pair of roots.  Calculate the bigger one in R1(1). */

	if (temp > 0.) {
	    r1[1] = -temp - sqdisc;
	} else {
	    r1[1] = -temp + sqdisc;
	}
	r1[2] = 0.;
	r2[1] = *beta / r1[1];
	r2[2] = 0.;
    }

    return 0;
} /* stiffc_ */

/* Subroutine */ int stiffd_(doublereal *v, doublereal *havg, doublereal *x, 
	doublereal *y, S_fp f, doublereal *fxy, doublereal *wt, doublereal *
	scale, doublereal *vdotv, doublereal *z__, doublereal *zdotz, 
	doublereal *vtemp)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer l;
    doublereal temp1, temp2;
    extern doublereal dotprd_(doublereal *, doublereal *, doublereal *, 
	    integer *);

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  External:     F */
/*  Input:        V(*), HAVG, X, Y(*), FXY(*), WT(*), SCALE, VDOTV, */
/*  Output:       Z(*), ZDOTZ */
/*  Workspace:    VTEMP(*) */

/*  For an input vector V(*) of length NEQ, this subroutine computes a vector */
/*  Z(*) that approximates the product HAVG*J*V where HAVG is an input scalar */
/*  and J is the Jacobian matrix of a function F evaluated at the input */
/*  arguments (X,Y(*)).  This function is defined by a subroutine of the form */
/*  F(T,U,F) that when given T and U(*), returns the value of the function in */
/*  F(*).  The input vector FXY(*) is defined by F(X,Y,FXY).  Scaling is a */
/*  delicate matter.  A weighted Euclidean norm is used with the (positive) */
/*  weights provided in WT(*).  The input scalar SCALE is the square root of */
/*  the unit roundoff times the norm of Y(*).  The square of the norm of the */
/*  input vector V(*) is input as VDOTV.  The routine outputs the square of */
/*  the norm of the output vector Z(*) as ZDOTZ.  The subroutine calls the */
/*  DOUBLE PRECISION FUNCTION DOTPRD(U,V,WT,NEQ) to compute the dot (inner) */
/*  product.  The vector VTEMP(*) is used for working storage. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Common Block for Problem Definition .. */
/*     .. Common Block to hold Problem Status .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*  Scale V(*) so that it can be used as an increment to Y(*) */
/*  for an accurate difference approximation to the Jacobian. */

    /* Parameter adjustments */
    --vtemp;
    --z__;
    --wt;
    --fxy;
    --y;
    --v;

    /* Function Body */
    temp1 = *scale / sqrt(*vdotv);
    i__1 = rkcom1_1.neqn;
    for (l = 1; l <= i__1; ++l) {
	vtemp[l] = y[l] + temp1 * v[l];
/* L20: */
    }

    (*f)(x, &vtemp[1], &z__[1]);
    ++rkcom2_1.nfcn;

/*  Form the difference approximation.  At the same time undo */
/*  the scaling of V(*) and introduce the factor of HAVG. */

    temp2 = *havg / temp1;
    i__1 = rkcom1_1.neqn;
    for (l = 1; l <= i__1; ++l) {
	z__[l] = temp2 * (z__[l] - fxy[l]);
/* L40: */
    }

    *zdotz = dotprd_(&z__[1], &z__[1], &wt[1], &rkcom1_1.neqn);

    return 0;
} /* stiffd_ */

doublereal dotprd_(doublereal *u, doublereal *v, doublereal *wt, integer *neq)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    integer l;
    doublereal sum;

/* ************************************************ */
/* **** NOT A DESIGNATED USER-CALLABLE ROUTINE **** */
/* ************************************************ */

/*  Purpose:   To compute a weighted Euclidean dot (inner) product of */
/*             two vectors. */

/*  Input:     U(*), V(*), WT(*), NEQ */
/*  Output:    the result DOTPRD is returned via the subprogram name */

/*  Comments: */
/*  ========= */
/*  The vectors U(*), V(*), and WT(*) are of length NEQ. The components */
/*  of WT(*) are weights that must be non-zero. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --wt;
    --v;
    --u;

    /* Function Body */
    sum = 0.;
    i__1 = *neq;
    for (l = 1; l <= i__1; ++l) {
	sum += u[l] / wt[l] * (v[l] / wt[l]);
/* L20: */
    }

    ret_val = sum;

    return ret_val;
} /* dotprd_ */

/* Subroutine */ int softfl_(logical *ask, logical *on)
{
    /* Initialized data */

    static logical soft = FALSE_;


/*  Purpose:      To prevent a program STOP after a "catastrophic" */
/*                failure when using a routine from RKSUITE. */

/*  Input:        ASK */
/*  Input/output: ON */

/*  Comments: */
/*  ========= */
/*  When a "catastrophic" failure is detected, the default action of */
/*  RKSUITE is to write an explanation to the standard output channel, */
/*  OUTCH, and STOP.  This subroutine can be used to prevent the STOP and */
/*  so allow the main program to continue.  To do this, you call SOFTFL with */
/*  ASK = .FALSE. and ON = .TRUE.  You must then call the subroutine CHKFL */
/*  after every call to a user-callable routine in RKSUITE to check whether */
/*  a catastrophic error occurred and take appropriate action if it did.  Of */
/*  course, you may call SETUP at any time to start a new problem, but calling */
/*  any other user-callable routine in RKSUITE after a catastrophic error will */
/*  lead to a STOP (even when "soft failure" has been set "on"). */

/*  When ON is set by a call to SOFTFL with ASK = .FALSE., the value of ON */
/*  is SAVEd.  The subroutine RKMSG in RKSUITE calls SOFTFL with ASK = .TRUE. */
/*  to find out the SAVEd value of ON. */

/*     .. Scalar Arguments .. */
/*     .. Local Scalars .. */
/*     .. Save statement .. */
/*     .. Data statements .. */
/*     .. Executable Statements .. */

    if (*ask) {
	*on = soft;
    } else {
	soft = *on;
    }

    return 0;
} /* softfl_ */

/* Subroutine */ int chkfl_(logical *ask, logical *error)
{
    /* Initialized data */

    static logical saverr = FALSE_;


/*  Purpose:      Enquiry routine used in conjunction with SOFTFL. */
/*                Reports whether a "catastrophic" error was detected. */

/*  Input:        ASK */
/*  Input/output: ERROR */

/*  Comments: */
/*  ========= */
/*  When a "catastrophic" failure is detected, the default action of */
/*  RKSUITE is to write an explanation to the standard output channel, */
/*  OUTCH, and STOP.  SOFTFL can be used to prevent the STOP and so */
/*  allow the main program to continue.  It is then necessary to call */
/*  CHKFL with ASK = .TRUE. after every call to a user-callable routine */
/*  in RKSUITE to check whether a catastrophic error occurred and take */
/*  appropriate action if it did.  If there was a catastrophic error, */
/*  ERROR is returned .TRUE.  Of course, you may call SETUP at any time */
/*  to start a new problem, but calling any other user-callable routine */
/*  in RKSUITE after a catastrophic error will lead to a STOP (even when */
/*  "soft failure" has been set "on"). */

/*  When a catastrophic failure (IER = 911) is detected in one of */
/*  the routines in RKSUITE, it calls CHKFL with ASK = .FALSE. and */
/*  ERROR = .TRUE.  This value of ERROR is SAVEd. */

/*     .. Scalar Arguments .. */
/*     .. Local Scalars .. */
/*     .. Save statement .. */
/*     .. Data statements .. */
/*     .. Executable Statements .. */

    if (*ask) {
	*error = saverr;
    } else {
	saverr = *error;
    }

    return 0;
} /* chkfl_ */

