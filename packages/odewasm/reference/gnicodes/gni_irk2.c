// -----------------------------------------------------------------------
//                 VERSION OF SEPTEMBER 4,2002
//  E-MAIL CONTACT ADDRESS : Ernst.Hairer@math.unige.ch
// -----------------------------------------------------------------------
//  SOLVES SECOND ORDER ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM
//                       Q'' = F(X,Q)
//  BASED ON THE SYMPLECTIC AND SYMMETRIC GAUSS (IRK) METHODS
//  DESCRIBED IN SECTIONS II.1, VIII.6 OF THE BOOK:
//
//      E. HAIRER, C. LUBICH, G. WANNER, GEOMETRIC NUMERICAL INTEGRATION,
//         STRUCTURE-PRESERVING ALGORITHMS FOR ODES.
//         SPRINGER SERIES IN COMPUT. MATH. 31, SPRINGER 2002.
//
//  AND IN THE PUBLICATION
//
//      E. HAIRER, M. HAIRER, GNI-CODES - MATLAB PROGRAMS FOR
//         GEOMETRIC NUMERICAL INTEGRATION.
//
//  INPUT..
//     ndim        DIMENSION OF Q AND F(X,Q)
//
//     task        An instance of class task_item.
//                 As such, it has methods including the
//                 function that evaluates the equation of motion
//                 and the function that prints the results.

//     NSTEP       NUMBER OF INTEGRATION STEPS
//                    CONSTANT STEP SIZE, H=(XEND-X)/NSTEP

//     X           INITIAL X-VALUE
//     P(ndim)     INITIAL VELOCITY VECTOR
//     Q(ndim)     INITIAL POSITION VECTOR
//     XEND        FINAL X-VALUE
//
//     METH        NUMBER OF STAGES OF THE GAUSS METHOD
//                    FOR THE MOMENT ONLY POSSIBLE VALUES: 2,4,6.
//     RPAR(LR)    REAL PARAMETER ARRAY; LR MUST BE AT LEAST LR=10
//                    RPAR(1),...,RPAR(10) SERVE AS PARAMETERS FOR
//                    THE CODE. FURTHER VALUES CAN BE USED FOR DEFINING
//                    PARAMETERS IN THE PROBLEM but you should never do
//                    this;  use member variables in the task_item class
//     IPAR(LI)    INTEGER PARAMETER ARRAY; LI MUST BE AT LEAST LI=10
//                    IPAR(1),...,IPAR(10) SERVE AS PARAMETERS FOR
//                    PARAMETERS IN THE PROBLEM but you should never do
//                    this;  use member variables in the task_item class

//  OUTPUT..
//     P(ndim)     SOLUTION (VELOCITY) AT XEND
//     Q(ndim)     SOLUTION (POSITION) AT XEND
//     function return value: nonzero if integration was interrupted
// -----------------------------------------------------------------------
//     SOPHISTICATED SETTING OF PARAMETERS
// -----------------------------------------------------------------------
//    RPAR(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
//    IPAR(1)   NITMAX, MAXIMAL NUMER OF FIXED POINT ITERAT., DEFAULT 50
// -----------------------------------------------------------------------
// -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0 */

using namespace std;
#include <iostream>
#include "task.h"
#include "gni_irk2.h"

/* cfeh dr_irk2 */
/* ----------------------------------------------------------------------- */

INT coefg_(INT *ns,
       FLP *cv,
       FLP *b,
       FLP *bc,
       const INT& nsd,
       FLP *aa,
       FLP *e,
       INT *nm,
       FLP *sm,
       const INT& nmd,
       FLP *am,
       FLP *hstep);

INT gni_irk2(const INT ndim, task_item* task, INT *nstep,
       const FLP *tbegin, FLP *p, FLP *q, FLP *tend,
       INT *meth, FLP *rpar, INT *
       ipar)
{
   /* System generated locals */
   INT i__1, i__2, i__3;
   FLP d__1;

   /* Local variables */
   static FLP b[6], cv[6], e[54]	/* was [6][9] */, f[3000],
           h__;
   static INT i__;
   static FLP aa[36]	/* was [6][6] */, bc[6], am[9], ep[500], eq[
           500], fs[500], qi;
   static INT nm, is;
   static FLP ay, pi, sm[3], yh[500], qq[500], ps[500];
   static INT ns;
   static FLP zq[3000]	/* was [500][6] */, fac, epi, eqi, dyno;
   static INT niter, istep;
   static FLP dynold;
   static INT nitmax;
   static FLP uround;

    /* Parameter adjustments */
    FLP* q_X1(q-1);
    FLP* p_X1(p-1);
    FLP* rpar_X1(rpar-1);
    int* ipar_X1(ipar-1);

    /* Function Body */
    if (rpar_X1[1] == 0.) {
	uround = 1e-16;
    } else {
	uround = rpar_X1[1];
    }
/* -------- NITMAX, MAXIMAL NUMER OF FIXED POINT ITERATIONS */
    if (ipar_X1[1] == 0) {
	nitmax = 50;
    } else {
	nitmax = ipar_X1[1];
    }
/* -------- */
    ns = *meth;
    h__ = (*tend - *tbegin) / *nstep;
    coefg_(&ns, cv, b, bc, 6, aa, e, &nm, sm, 3, am, &h__);
    int flag = task->solfix(0, *tbegin, *tbegin, &p_X1[1], &q_X1[1], ndim,
        &rpar_X1[1], &ipar_X1[1]);
    if (flag) return flag;
    (task->EoM)(ndim, *tbegin, &q_X1[1], fs, &rpar_X1[1], &ipar_X1[1]);
    i__1 = ns;
    for (is = 1; is <= i__1; ++is) {
/* Computing 2nd power */
	d__1 = cv[is - 1];
	fac = d__1 * d__1 / 2;
	i__2 = ndim;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    zq[i__ + is * 500 - 501] = cv[is - 1] * p_X1[i__] + fac * fs[i__ -
		    1];
	}
    }
    i__1 = ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ps[i__ - 1] = p_X1[i__];
	eq[i__ - 1] = 0.;
	ep[i__ - 1] = 0.;
    }
/* --- LOOP FOR THE ITERATIONS */
    i__1 = *nstep;
    FLP tnow(*tbegin);
    for (istep = 1; istep <= i__1; ++istep) {
	if (istep > 1) {
	    startb(task, ndim, tnow,
                &p_X1[1], &q_X1[1], &ns, 500, fs, ps, zq,
		    6, e, yh, &nm, sm, am, f, cv, &rpar_X1[1], &ipar_X1[1]);
	}
/* --- FIXED POINT ITERATION */
	niter = 0;
	dynold = 0.;
L40:
	rknite(task, ndim, &ns, tnow, &q_X1[1], &p_X1[1], 6, aa, cv, 500,
		qq, zq, f, &dyno, &rpar_X1[1], &ipar_X1[1]);
	++niter;
	if (dynold < dyno && dyno < uround * 10) {
	    goto L50;
	}
	if (niter >= nitmax) {
            cout << "no convergence;  iteration: " << niter
                << "  dyno: " << dyno
                << endl;
	    if (dyno > uround * 1e6) {
		return 0;
	    }
	}
	if (dyno > uround * .1) {
	    goto L40;
	}
L50:
/* --- UPDATE OF THE SOLUTION */
	tnow += h__;
	i__2 = ndim;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    eqi = eq[i__ - 1];
	    qi = q_X1[i__];
	    ay = qi;
	    eqi += h__ * p_X1[i__];
	    qi = ay + eqi;
	    eqi += ay - qi;
	    i__3 = ns;
	    for (is = 1; is <= i__3; ++is) {
		ay = qi;
		eqi += f[i__ + (is - 1) * ndim - 1] * bc[is - 1];
		qi = ay + eqi;
		eqi += ay - qi;
	    }
	    ay = q_X1[i__];
	    q_X1[i__] = qi;
	    eq[i__ - 1] = eqi;
	    epi = ep[i__ - 1];
	    pi = p_X1[i__];
	    i__3 = ns / 2;
	    for (is = 1; is <= i__3; ++is) {
		ay = pi;
		epi += (f[i__ + (is - 1) * ndim - 1] + f[i__ + (ns - is) * ndim -
			1]) * b[is - 1];
		pi = ay + epi;
		epi += ay - pi;
	    }
	    p_X1[i__] = pi;
	    ep[i__ - 1] = epi;
	}
        FLP told = tnow - h__;
        flag = task->solfix(istep, told, tnow, &p_X1[1], &q_X1[1], ndim,
                &rpar_X1[1], &ipar_X1[1]);
        if (flag) return flag;
    }
    return 0;
} /* gni_irk2__ */


INT startb(task_item* task,
        const INT ndim, const FLP tnow, FLP *p,
        FLP *q, INT *ns, INT ndgl, FLP *fs,
	FLP *ps, FLP *zq, INT nsd, FLP *e,
	FLP *yh, INT *nm, FLP *sm, FLP *am,
	FLP *f, FLP *cv, FLP *rpar, INT *ipar)
{
    /* System generated locals */
    INT zq_dim1, zq_offset, e_dim1, e_offset, i__1, i__2, i__3;
    FLP d__1;

    /* Local variables */
    static FLP h__;
    static INT i__, is, js, ns1, ns2;
    static FLP sav;
    static INT nsm;
    static FLP zqiis;

/* ---------------------------------------------------------- */
    /* Parameter adjustments */
    FLP* ps_X1(ps-1);
    FLP* q_X1(q-1);
    FLP* p_X1(p-1);
    FLP* cv_X1(cv-1);
    FLP* f_X1(f-1);
    FLP* yh_X1(yh-1);
    zq_dim1 = ndgl;
    zq_offset = 1 + zq_dim1;
    FLP* zq_X1(zq - zq_offset);
    FLP* fs_X1(fs-1);
    FLP* am_X1(am-1);
    FLP* sm_X1(sm-1);
    e_dim1 = nsd;
    e_offset = 1 + e_dim1;
    FLP* e_X1(e - e_offset);
    FLP* rpar_X1(rpar-1);
    int* ipar_X1(ipar-1);

    /* Function Body */
    ns1 = *ns + 1;
    ns2 = *ns + 2;
    nsm = *ns + *nm;
    i__1 = ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sav = 0.;
	i__2 = *ns;
	for (js = 1; js <= i__2; ++js) {
	    sav += am_X1[js] * zq_X1[i__ + js * zq_dim1];
	}
	yh_X1[i__] = sav + am_X1[ns1] * ps_X1[i__] + am_X1[ns2] * p_X1[i__] + q_X1[i__];
	i__2 = *ns;
	for (is = 1; is <= i__2; ++is) {
	    sav = 0.;
	    i__3 = *ns;
	    for (js = 1; js <= i__3; ++js) {
		sav += e_X1[is + js * e_dim1] * f_X1[i__ + (js - 1) * ndim];
	    }
	    zq_X1[i__ + is * zq_dim1] = sav + e_X1[is + ns1 * e_dim1] * fs_X1[i__];
	}
    }

    d__1 = tnow + h__;
    task->EoM(ndim, d__1, &q_X1[1], &fs_X1[1], &rpar_X1[1], &ipar_X1[1]);
    d__1 = tnow + h__ * sm_X1[*nm];
    task->EoM(ndim, d__1, &yh_X1[1], &f_X1[1], &rpar_X1[1], &ipar_X1[1]);
    i__1 = ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ps_X1[i__] = p_X1[i__];
	i__2 = *ns;
	for (is = 1; is <= i__2; ++is) {
	    zqiis = zq_X1[i__ + is * zq_dim1] + e_X1[is + ns2 * e_dim1] * fs_X1[i__] +
		    e_X1[is + nsm * e_dim1] * f_X1[i__];
	    zq_X1[i__ + is * zq_dim1] = zqiis + cv_X1[is] * p_X1[i__];
	}
    }
    return 0;
} /* startb */


INT rknite(task_item* task, const INT ndim,
        INT *ns, const FLP tnow,
        FLP *q, FLP *p, INT nsd, FLP *aa,
	FLP *c_, INT ndgl, FLP *qq, FLP *zq,
	FLP *f, FLP *dyno, FLP *rpar, INT *ipar)
{
    /* System generated locals */
    INT zq_dim1, zq_offset, aa_dim1, aa_offset, i__1, i__2, i__3;
    FLP d__1, d__2, d__3;

    /* Local variables */
    static INT i__, j, is, js;
    static FLP sum, dnom;

/* ---------------------------------------------------------- */
/* --- */
    /* Parameter adjustments */
    FLP* p_X1(p-1);
    FLP* q_X1(q-1);
    FLP* cv_X1(c_-1);
    aa_dim1 = nsd;
    aa_offset = 1 + aa_dim1;
    FLP* aa_X1(aa - aa_offset);
    FLP* f_X1(f-1);
    zq_dim1 = ndgl;
    zq_offset = 1 + zq_dim1;
    FLP* zq_X1(zq - zq_offset);
    FLP* qq_X1(qq-1);
    FLP* rpar_X1(rpar-1);
    int* ipar_X1(ipar-1);

    /* Function Body */
    i__1 = *ns;
    for (js = 1; js <= i__1; ++js) {
	i__2 = ndim;
	for (j = 1; j <= i__2; ++j) {
	    qq_X1[j] = q_X1[j] + zq_X1[j + js * zq_dim1];
	}
	d__1 = tnow + cv_X1[js];
	task->EoM(ndim, d__1, &qq_X1[1], &f_X1[(js - 1) * ndim + 1], &rpar_X1[1], &ipar_X1[1]);
    }
/* --- */
    *dyno = 0.;
    i__1 = ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = .1, d__3 = (d__1 = q_X1[i__], fabs(d__1));
	dnom = max(d__2,d__3);
	i__2 = *ns;
	for (is = 1; is <= i__2; ++is) {
	    sum = cv_X1[is] * p_X1[i__];
	    i__3 = *ns;
	    for (js = 1; js <= i__3; ++js) {
		sum += aa_X1[is + js * aa_dim1] * f_X1[i__ + (js - 1) * ndim];
	    }
/* Computing 2nd power */
	    d__1 = (sum - zq_X1[i__ + is * zq_dim1]) / dnom;
	    *dyno += d__1 * d__1;
	    zq_X1[i__ + is * zq_dim1] = sum;
	}
    }
    *dyno = sqrt(*dyno / (*ns * ndim));
/* --- */
    return 0;
} /* rknite */


/* Subroutine */ INT coefg_(INT *ns, FLP *cv, FLP *bv,
	FLP *bc,
        const INT& nsd, FLP *aa, FLP *e,
        INT *nm,
        FLP *sm,
        const INT& nmd, FLP *am, FLP *hstep)
{
    /* System generated locals */
    INT aa_dim1, aa_offset, e_dim1, e_offset, i__1, i__2;

    /* Local variables */
    static INT im, is, js;
    static FLP hstep2;

    /* Parameter adjustments */
    FLP* bc_X1(bc-1);
    FLP* bv_X1(bv-1);
    FLP* cv_X1(cv-1);
    aa_dim1 = nsd;
    aa_offset = 1 + aa_dim1;
    FLP* aa_X1(aa - aa_offset);
    FLP* am_X1(am-1);
    FLP* sm_X1(sm-1);
    e_dim1 = nsd;
    e_offset = 1 + e_dim1;
    FLP* e_X1(e - e_offset);
/* Function Body */
    *nm = 3;
    if (*ns == 2) {     // emacs: esc-1 esc | calc_new -ns 2
         cv_X1[1]=  .21132486540518710671499036E+00L;
         cv_X1[2]=  .78867513459481286552943402E+00L;
         bv_X1[1]=  .50000000000000000000000000E+00L;
         bv_X1[2]=  .50000000000000000000000000E+00L;
         bc_X1[1]=  .39433756729740643276471701E+00L;
         bc_X1[2]=  .10566243270259355335749518E+00L;
         aa_X1[1 + 1*aa_dim1]=  .41666666666666664353702032E-01;
         aa_X1[1 + 2*aa_dim1]= -.19337567297406439703610914E-01;
         aa_X1[2 + 1*aa_dim1]=  .26933756729740643276471701E+00;
         aa_X1[2 + 2*aa_dim1]=  .41666666666666664353702032E-01;
         e_X1[1 + 1*e_dim1]= -.76474568177964663551193780E-01;
         e_X1[1 + 2*e_dim1]=  .17615393673685059150457732E+00;
         e_X1[1 + 3*e_dim1]=  .29533739145129821179747154E-01;
         e_X1[1 + 4*e_dim1]= -.10688400833475553142193348E+00;
         e_X1[1 + 5*e_dim1]=  .00000000000000000000000000E+00;
         e_X1[2 + 1*e_dim1]=  .83748242689951313266760735E+00;
         e_X1[2 + 2*e_dim1]= -.16038284621250658013025259E+01;
         e_X1[2 + 3*e_dim1]= -.41135192096331174571588463E+00;
         e_X1[2 + 4*e_dim1]=  .14887021901529375877970551E+01;
         e_X1[2 + 5*e_dim1]=  .00000000000000000000000000E+00;
         sm_X1[1]=  .00000000000000000000000000E+00L;
         sm_X1[2]=  .10000000000000000000000000E+01L;
         sm_X1[3]=  .16000000000000000888178420E+01L;
         am_X1[1]=  .25279583039343449968328059E+02L;
         am_X1[2]= -.86907830393434419846698802E+01L;
         am_X1[3]= -.80640000000000033875124927E+00L;
         am_X1[4]=  .29184000000000009933387446E+01L;
         am_X1[5]=  .00000000000000000000000000E+00L;
    }
    if (*ns == 4) { // emacs: esc-1 esc | calc_new -ns 4
         cv_X1[1]=  .69431844202973713731097405E-01L;
         cv_X1[2]=  .33000947820757187134432797E+00L;
         cv_X1[3]=  .66999052179242812865567203E+00L;
         cv_X1[4]=  .93056815579702634178005383E+00L;
         bv_X1[1]=  .17392742256872692485636378E+00L;
         bv_X1[2]=  .32607257743127304738806060E+00L;
         bv_X1[3]=  .32607257743127304738806060E+00L;
         bv_X1[4]=  .17392742256872692485636378E+00L;
         bc_X1[1]=  .16185132086231029946432614E+00L;
         bc_X1[2]=  .21846553629538056906511656E+00L;
         bc_X1[3]=  .10760704113589250607851966E+00L;
         bc_X1[4]=  .12076101706416621922590693E-01L;
         aa_X1[1 + 1*aa_dim1]=  .40381914508467314783857205E-02;
         aa_X1[1 + 2*aa_dim1]= -.32958609449446961241203535E-02;
         aa_X1[1 + 3*aa_dim1]=  .26447829520668537754690686E-02;
         aa_X1[1 + 4*aa_dim1]= -.97672296325588170602671756E-03;
         aa_X1[2 + 1*aa_dim1]=  .43563580902396259464381956E-01;
         aa_X1[2 + 2*aa_dim1]=  .13818951406296126407924341E-01;
         aa_X1[2 + 3*aa_dim1]= -.43401341944349957138737928E-02;
         aa_X1[2 + 4*aa_dim1]=  .14107297391595338166558893E-02;
         aa_X1[3 + 1*aa_dim1]=  .10586435263357640845782726E+00;
         aa_X1[3 + 2*aa_dim1]=  .10651836096505307160953180E+00;
         aa_X1[3 + 3*aa_dim1]=  .13818951406296126407924341E-01;
         aa_X1[3 + 4*aa_dim1]= -.17580153590805495077142862E-02;
         aa_X1[4 + 1*aa_dim1]=  .14879849619263779691991090E+00;
         aa_X1[4 + 2*aa_dim1]=  .19847049885237719180075544E+00;
         aa_X1[4 + 3*aa_dim1]=  .81671359795877571108313475E-01;
         aa_X1[4 + 4*aa_dim1]=  .40381914508467314783857205E-02;
         e_X1[1 + 1*e_dim1]= -.25294291374367464109162995E-01;
         e_X1[1 + 2*e_dim1]=  .27329964075078944790719859E-01;
         e_X1[1 + 3*e_dim1]= -.40345491602452029922254439E-01;
         e_X1[1 + 4*e_dim1]=  .68940879518524886271357843E-01;
         e_X1[1 + 5*e_dim1]=  .13699161458581856146143885E-01;
         e_X1[1 + 6*e_dim1]= -.41919831580653182934170786E-01;
         e_X1[1 + 7*e_dim1]=  .00000000000000000000000000E+00;
         e_X1[2 + 1*e_dim1]=  .27598778734680451396599210E+00;
         e_X1[2 + 2*e_dim1]= -.23490311281159778622473766E+00;
         e_X1[2 + 3*e_dim1]=  .30637582212469710585267535E+00;
         e_X1[2 + 4*e_dim1]= -.62857823183857075566294270E+00;
         e_X1[2 + 5*e_dim1]= -.16289618260612598321657174E+00;
         e_X1[2 + 6*e_dim1]=  .49846704563820987132771734E+00;
         e_X1[2 + 7*e_dim1]=  .00000000000000000000000000E+00;
         e_X1[3 + 1*e_dim1]=  .31925410269310541266918335E+01;
         e_X1[3 + 2*e_dim1]= -.28107548361348548837668204E+01;
         e_X1[3 + 3*e_dim1]=  .38816354537247685918543993E+01;
         e_X1[3 + 4*e_dim1]= -.78928482234580394916179102E+01;
         e_X1[3 + 5*e_dim1]= -.18707844382649314596278600E+01;
         e_X1[3 + 6*e_dim1]=  .57246546668478481834085869E+01;
         e_X1[3 + 7*e_dim1]=  .00000000000000000000000000E+00;
         e_X1[4 + 1*e_dim1]=  .12717659277813599061346395E+02;
         e_X1[4 + 2*e_dim1]= -.11634766434694126857607444E+02;
         e_X1[4 + 3*e_dim1]=  .16462584365145023213017339E+02;
         e_X1[4 + 4*e_dim1]= -.32284611573813158713619487E+02;
         e_X1[4 + 5*e_dim1]= -.73649996101467776199456239E+01;
         e_X1[4 + 6*e_dim1]=  .22537112521987179292182191E+02;
         e_X1[4 + 7*e_dim1]=  .00000000000000000000000000E+00;
         sm_X1[1]=  .00000000000000000000000000E+00L;
         sm_X1[2]=  .10000000000000000000000000E+01L;
         sm_X1[3]=  .16499999999999999111821580E+01L;
         am_X1[1]=  .10806374869243995817669202E+04L;
         am_X1[2]= -.66008818661284658446675166E+03L;
         am_X1[3]=  .61810154357557496496156091E+03L;
         am_X1[4]= -.31341427826212844820474857E+03L;
         am_X1[5]= -.10187174765624995131929609E+02L;
         am_X1[6]=  .31173050390624986505372362E+02L;
         am_X1[7]=  .00000000000000000000000000E+00L;
    }
    if (*ns == 6) {    // emacs: esc-1 esc | calc_new
         cv_X1[1]=  .33765242898423988848755073E-01L;
         cv_X1[2]=  .16939530676686773147388010E+00L;
         cv_X1[3]=  .38069040695840156152129907E+00L;
         cv_X1[4]=  .61930959304159849398985216E+00L;
         cv_X1[5]=  .83060469323313224077054429E+00L;
         cv_X1[6]=  .96623475710157602502903273E+00L;
         bv_X1[1]=  .85662246189585178335335058E-01L;
         bv_X1[2]=  .18038078652406930313389921E+00L;
         bv_X1[3]=  .23395696728634551853076573E+00L;
         bv_X1[4]=  .23395696728634551853076573E+00L;
         bv_X1[5]=  .18038078652406930313389921E+00L;
         bv_X1[6]=  .85662246189585178335335058E-01L;
         bc_X1[1]=  .82769839639769235417610105E-01L;
         bc_X1[2]=  .14982512785597568161222171E+00L;
         bc_X1[3]=  .14489179419935321879719936E+00L;
         bc_X1[4]=  .89065173086992313611354177E-01L;
         bc_X1[5]=  .30555658668093604174442746E-01L;
         bc_X1[6]=  .28924065498159377135545256E-02L;
         aa_X1[1 + 1*aa_dim1]=  .90625420195651151897969777E-03;
         aa_X1[1 + 2*aa_dim1]= -.72859711612531399913678110E-03;
         aa_X1[1 + 3*aa_dim1]=  .79102695861167692057552836E-03;
         aa_X1[1 + 4*aa_dim1]= -.70675390218535383626830004E-03;
         aa_X1[1 + 5*aa_dim1]=  .45647714224056920410965699E-03;
         aa_X1[1 + 6*aa_dim1]= -.14836147050330408076954103E-03;
         aa_X1[2 + 1*aa_dim1]=  .11272367531794365638764255E-01;
         aa_X1[2 + 2*aa_dim1]=  .39083482447840696069607525E-02;
         aa_X1[2 + 3*aa_dim1]= -.14724868010943911907084658E-02;
         aa_X1[2 + 4*aa_dim1]=  .10992669056588430607362961E-02;
         aa_X1[2 + 5*aa_dim1]= -.67689040729401424908301399E-03;
         aa_X1[2 + 6*aa_dim1]=  .21677950347174140864961456E-03;
         aa_X1[3 + 1*aa_dim1]=  .30008019623627545796606952E-01;
         aa_X1[3 + 2*aa_dim1]=  .36978289259468145877551848E-01;
         aa_X1[3 + 3*aa_dim1]=  .65490339168957825005668028E-02;
         aa_X1[3 + 4*aa_dim1]= -.16615098173008262250061051E-02;
         aa_X1[3 + 5*aa_dim1]=  .84753461862041611692836218E-03;
         aa_X1[3 + 6*aa_dim1]= -.25877462623437422510261352E-03;
         aa_X1[4 + 1*aa_dim1]=  .49900269650650898312083115E-01;
         aa_X1[4 + 2*aa_dim1]=  .82003427445271614981692210E-01;
         aa_X1[4 + 3*aa_dim1]=  .54165111295060068552498223E-01;
         aa_X1[4 + 4*aa_dim1]=  .65490339168957816332050648E-02;
         aa_X1[4 + 5*aa_dim1]= -.11352871017627468985117245E-02;
         aa_X1[4 + 6*aa_dim1]=  .28963081055952377392340158E-03;
         aa_X1[5 + 1*aa_dim1]=  .68475836671617246187437900E-01;
         aa_X1[5 + 2*aa_dim1]=  .11859257878058808433063831E+00;
         aa_X1[5 + 3*aa_dim1]=  .10635984886129551396649617E+00;
         aa_X1[5 + 4*aa_dim1]=  .47961474042181380383897960E-01;
         aa_X1[5 + 5*aa_dim1]=  .39083482447840704743224904E-02;
         aa_X1[5 + 6*aa_dim1]= -.34600839001342452178749953E-03;
         aa_X1[6 + 1*aa_dim1]=  .79729071619449989882788543E-01;
         aa_X1[6 + 2*aa_dim1]=  .14419100392702230428731980E+00;
         aa_X1[6 + 3*aa_dim1]=  .13628542646896577017479046E+00;
         aa_X1[6 + 4*aa_dim1]=  .81956586217401899974177581E-01;
         aa_X1[6 + 5*aa_dim1]=  .23736460480774324022235078E-01;
         aa_X1[6 + 6*aa_dim1]=  .90625420195651141055948052E-03;
         e_X1[1 + 1*e_dim1]= -.14683890443713711396189758E-01;
         e_X1[1 + 2*e_dim1]=  .12267719252911186630927709E-01;
         e_X1[1 + 3*e_dim1]= -.13721292099138041042927227E-01;
         e_X1[1 + 4*e_dim1]=  .17284063971557486055763420E-01;
         e_X1[1 + 5*e_dim1]= -.23467545192348750887623154E-01;
         e_X1[1 + 6*e_dim1]=  .37038410388557908847761979E-01;
         e_X1[1 + 7*e_dim1]=  .86345757899393606876747853E-02;
         e_X1[1 + 8*e_dim1]= -.22781995853770650728442604E-01;
         e_X1[1 + 9*e_dim1]=  .00000000000000000000000000E+00;
         e_X1[2 + 1*e_dim1]=  .11319881604039469524813910E+00;
         e_X1[2 + 2*e_dim1]= -.83351663787741198108349749E-01;
         e_X1[2 + 3*e_dim1]=  .80335070059241703299690585E-01;
         e_X1[2 + 4*e_dim1]= -.91591316463245198153941828E-01;
         e_X1[2 + 5*e_dim1]=  .12617084620765961489041729E+00;
         e_X1[2 + 6*e_dim1]= -.24345393346076149576617809E+00;
         e_X1[2 + 7*e_dim1]= -.68991285957545139817526092E-01;
         e_X1[2 + 8*e_dim1]=  .18203085233931762254577791E+00;
         e_X1[2 + 9*e_dim1]=  .00000000000000000000000000E+00;
         e_X1[3 + 1*e_dim1]=  .22426502882937966454335310E+01;
         e_X1[3 + 2*e_dim1]= -.16410163756541664703547667E+01;
         e_X1[3 + 3*e_dim1]=  .15869185051017999921185719E+01;
         e_X1[3 + 4*e_dim1]= -.18677442282417928698379228E+01;
         e_X1[3 + 5*e_dim1]=  .27233402774076758134924603E+01;
         e_X1[3 + 6*e_dim1]= -.52164345820792457075754101E+01;
         e_X1[3 + 7*e_dim1]= -.13700344488544737853175093E+01;
         e_X1[3 + 8*e_dim1]=  .36147831570014838931115264E+01;
         e_X1[3 + 9*e_dim1]=  .00000000000000000000000000E+00;
         e_X1[4 + 1*e_dim1]=  .23382153299762862275201769E+02;
         e_X1[4 + 2*e_dim1]= -.17521634244141424119334260E+02;
         e_X1[4 + 3*e_dim1]=  .17575690145495602934033741E+02;
         e_X1[4 + 4*e_dim1]= -.21391501223389877139879900E+02;
         e_X1[4 + 5*e_dim1]=  .31180056211340744454219021E+02;
         e_X1[4 + 6*e_dim1]= -.56300296535982262469133275E+02;
         e_X1[4 + 7*e_dim1]= -.14200702567022359446013979E+02;
         e_X1[4 + 8*e_dim1]=  .37468007099953389626989519E+02;
         e_X1[4 + 9*e_dim1]=  .00000000000000000000000000E+00;
         e_X1[5 + 1*e_dim1]=  .11139364767813630407999881E+03;
         e_X1[5 + 2*e_dim1]= -.84723572571360307392751565E+02;
         e_X1[5 + 3*e_dim1]=  .86716574160386812764045317E+02;
         e_X1[5 + 4*e_dim1]= -.10700093080334335127190570E+03;
         e_X1[5 + 5*e_dim1]=  .15482887797769487292498525E+03;
         e_X1[5 + 6*e_dim1]= -.27129109870701967111017439E+03;
         e_X1[5 + 7*e_dim1]= -.67393377171550454818316211E+02;
         e_X1[5 + 8*e_dim1]=  .17781483151526623487370671E+03;
         e_X1[5 + 9*e_dim1]=  .00000000000000000000000000E+00;
         e_X1[6 + 1*e_dim1]=  .26074561352854505003051599E+03;
         e_X1[6 + 2*e_dim1]= -.19989352725426564916233474E+03;
         e_X1[6 + 3*e_dim1]=  .20661428715509001108330267E+03;
         e_X1[6 + 4*e_dim1]= -.25628367171165763238604995E+03;
         e_X1[6 + 5*e_dim1]=  .36890868592525151825611829E+03;
         e_X1[6 + 6*e_dim1]= -.63754672856348258846992394E+03;
         e_X1[6 + 7*e_dim1]= -.15741727502997960641550890E+03;
         e_X1[6 + 8*e_dim1]=  .41533942075341445843150723E+03;
         e_X1[6 + 9*e_dim1]=  .00000000000000000000000000E+00;
         sm_X1[1]=  .00000000000000000000000000E+00L;
         sm_X1[2]=  .10000000000000000000000000E+01L;
         sm_X1[3]=  .17500000000000000000000000E+01L;
         am_X1[1]=  .58080578375796241743955761E+05L;
         am_X1[2]= -.33214989339522820955608040E+05L;
         am_X1[3]=  .28376088288311992073431611E+05L;
         am_X1[4]= -.27923430684614970232360065E+05L;
         am_X1[5]=  .29743005589491014688974246E+05L;
         am_X1[6]= -.15525927919158812073874287E+05L;
         am_X1[7]= -.27700591278076012713427190E+03L;
         am_X1[8]=  .73086943817138603662897367E+03L;
         am_X1[9]=  .00000000000000000000000000E+00L;
    }

    hstep2 = *hstep * *hstep;
    i__1 = *ns;
    for (is = 1; is <= i__1; ++is) {
	bv_X1[is] = *hstep * bv_X1[is];
	bc_X1[is] = hstep2 * bc_X1[is];
	cv_X1[is] = *hstep * cv_X1[is];
	i__2 = *ns;
	for (js = 1; js <= i__2; ++js) {
	    aa_X1[is + js * aa_dim1] = hstep2 * aa_X1[is + js * aa_dim1];
	    e_X1[is + js * e_dim1] = hstep2 * e_X1[is + js * e_dim1];
	}
    }
    i__1 = *nm;
    for (im = 1; im <= i__1; ++im) {
	i__2 = *ns;
	for (is = 1; is <= i__2; ++is) {
	    e_X1[is + (*ns + im) * e_dim1] = hstep2 * e_X1[is + (*ns + im) * e_dim1]
		    ;
	}
	am_X1[*ns + im] = *hstep * am_X1[*ns + im];
    }

    return 0;
} /* coefg_ */