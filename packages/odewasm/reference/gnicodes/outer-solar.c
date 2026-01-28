using namespace std;
#include <iostream>
#include <math.h>

#include "task.h"

class outer_item : public tasklist_item {
public:
  outer_item(const string _name) :
        tasklist_item(_name) {}
  virtual informer info;
  virtual setupper setup;
  virtual EoM_eval EoM;
  virtual hamiltoner hamil;
};

outer_item my_outer_item("outer-solar-system");

void outer_item::info(FLP *tbegin, FLP *tend,
        FLP* stepsize, INT *ndim){
    *ndim = 18;
    *tbegin = 0.;
    *tend = 5e5;
    *stepsize = 0.1;
    odt = (*tend - *tbegin + 0.5) / 1000.;
}


///////////////////////////
// set up the initial conditions
//
// qex, pex are the values to be *expected* at the end of the run
//
void outer_item::setup(const FLP /* tbegin not used */,
        const FLP /* tend not used */,
	INT /* ndim not used */,
        FLP *q, FLP *p, FLP *qex, FLP
	*pex, FLP *rpar, INT *ipar)
{
    /* Parameter adjustments */
    FLP* pex_X1(pex-1);
    FLP* qex_X1(qex-1);
    FLP* p_X1(p-1);
    FLP* q_X1(q-1);

    /* Function Body */
    q_X1[1] = -3.5023653;
    q_X1[2] = -3.8169847;
    q_X1[3] = -1.5507963;
    q_X1[4] = 9.0755314;
    q_X1[5] = -3.0458353;
    q_X1[6] = -1.6483708;
    q_X1[7] = 8.310142;
    q_X1[8] = -16.2901086;
    q_X1[9] = -7.2521278;
    q_X1[10] = 11.4707666;
    q_X1[11] = -25.7294829;
    q_X1[12] = -10.8169456;
    q_X1[13] = -15.5387357;
    q_X1[14] = -25.2225594;
    q_X1[15] = -3.1902382;
    q_X1[16] = 0.;
    q_X1[17] = 0.;
    q_X1[18] = 0.;
    p_X1[1] = .00565429;
    p_X1[2] = -.0041249;
    p_X1[3] = -.00190589;
    p_X1[4] = .00168318;
    p_X1[5] = .00483525;
    p_X1[6] = .00192462;
    p_X1[7] = .00354178;
    p_X1[8] = .00137102;
    p_X1[9] = 5.5029e-4;
    p_X1[10] = .0028893;
    p_X1[11] = .00114527;
    p_X1[12] = 3.9677e-4;
    p_X1[13] = .00276725;
    p_X1[14] = -.00170702;
    p_X1[15] = -.00136504;
    p_X1[16] = 0.;
    p_X1[17] = 0.;
    p_X1[18] = 0.;
    qex_X1[1] = 7.766584086800482;
    qex_X1[2] = .2531065754551048;
    qex_X1[3] = -.09410571402013185;
    qex_X1[4] = -5.564967162844037;
    qex_X1[5] = 1.674849740822012;
    qex_X1[6] = .9767232069533176;
    qex_X1[7] = 19.63899572895227;
    qex_X1[8] = 8.95850455228646;
    qex_X1[9] = 3.611839157057347;
    qex_X1[10] = 24.93570870305177;
    qex_X1[11] = 17.69518676153705;
    qex_X1[12] = 6.583785164549242;
    qex_X1[13] = 31.78592511375764;
    qex_X1[14] = 38.63618958160644;
    qex_X1[15] = 3.192794169732889;
    qex_X1[16] = 3.084118473380683;
    qex_X1[17] = -1.227726356581642;
    qex_X1[18] = -.6162537634647217;
    pex_X1[1] = -.002495503201917009;
    pex_X1[2] = .006896467194473328;
    pex_X1[3] = .003007950247474123;
    pex_X1[4] = -.002255335935351989;
    pex_X1[5] = -.004905913854771086;
    pex_X1[6] = -.001938473641716708;
    pex_X1[7] = -.002186170231167942;
    pex_X1[8] = .002817177012110666;
    pex_X1[9] = .001262882639181183;
    pex_X1[10] = -.002148728705895163;
    pex_X1[11] = .002128650077635786;
    pex_X1[12] = 9.248501411662923e-4;
    pex_X1[13] = -.001675173186229401;
    pex_X1[14] = .001011833320388655;
    pex_X1[15] = 8.23180003857652e-4;
    pex_X1[16] = 9.417379703028725e-6;
    pex_X1[17] = -7.855256238249194e-6;
    pex_X1[18] = -3.646926313230521e-6;
} /* outer_item::setup */


//////////////////////
// evaluate the equation of motion
//
void outer_item::EoM(const INT /* ndim not used */,
        const FLP /* tnow not used */, FLP *q,
	FLP *f, FLP *rpar, INT *ipar)
{
    /* System generated locals */
    FLP d__1, d__2, d__3, d__4;

    /* Local variables */
    static FLP d__[36]	/* was [6][6] */;
    static INT i__, j;
    static FLP k, m[6];
    static INT i1, j1;

    /* Parameter adjustments */
    FLP* f_X1(f-1);
    FLP* q_X1(q-1);

    /* Function Body */
    ++eom_evals;
    k = 2.95912208286e-4;
    m[0] = 9.54786104043e-4;
    m[1] = 2.85583733151e-4;
    m[2] = 4.37273164546e-5;
    m[3] = 5.17759138449e-5;
    m[4] = 7.6923076923076926e-9;
    m[5] = 1.00000597682;
    for (i__ = 1; i__ <= 5; ++i__) {
        i1 = (i__ - 1) * 3 + 1;
        for (j = i__ + 1; j <= 6; ++j) {
            j1 = (j - 1) * 3 + 1;
/* Computing 2nd power */
            d__2 = q_X1[i1] - q_X1[j1];
/* Computing 2nd power */
            d__3 = q_X1[i1 + 1] - q_X1[j1 + 1];
/* Computing 2nd power */
            d__4 = q_X1[i1 + 2] - q_X1[j1 + 2];
/* Computing 3rd power */
            d__1 = sqrt(d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
            d__[i__ + j * 6 - 7] = d__1 * (d__1 * d__1);
            d__[j + i__ * 6 - 7] = d__[i__ + j * 6 - 7];
        }
    }
    for (i__ = 1; i__ <= 6; ++i__) {
        i1 = (i__ - 1) * 3 + 1;
        f_X1[i1] = 0.;
        f_X1[i1 + 1] = 0.;
        f_X1[i1 + 2] = 0.;
        for (j = 1; j <= 6; ++j) {
            if (j != i__) {
                j1 = (j - 1) * 3 + 1;
                f_X1[i1] += m[j - 1] * (q_X1[j1] - q_X1[i1]) / d__[i__ + j * 6 - 7]
                        ;
                f_X1[i1 + 1] += m[j - 1] * (q_X1[j1 + 1] - q_X1[i1 + 1]) / d__[i__
                        + j * 6 - 7];
                f_X1[i1 + 2] += m[j - 1] * (q_X1[j1 + 2] - q_X1[i1 + 2]) / d__[i__
                        + j * 6 - 7];
            }
        }
        f_X1[i1] = k * f_X1[i1];
        f_X1[i1 + 1] = k * f_X1[i1 + 1];
        f_X1[i1 + 2] = k * f_X1[i1 + 2];
    }
} /* outer_item::EoM */


////////////////////////////////
// evaluate the Hamiltonian
//
FLP outer_item::hamil(const INT ndim,
        const FLP tnow, FLP *q,
	FLP *p, FLP *rpar, INT *ipar)
{
    /* System generated locals */
    INT i__1;
    FLP d__1, d__2, d__3;

    /* Local variables */
    static FLP d__[36]	/* was [6][6] */;
    static INT i__, j;
    static FLP k, m[6], y[100];
    static INT i1, j1;
    static FLP pot;
    FLP rslt(0);

    /* Parameter adjustments */
    FLP* p_X1(p-1);
    FLP* q_X1(q-1);

    /* Function Body */
    i__1 = ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
        y[i__ - 1] = q_X1[i__];
        y[i__ + ndim - 1] = p_X1[i__];
    }
    k = 2.95912208286e-4;
    m[0] = 9.54786104043e-4;
    m[1] = 2.85583733151e-4;
    m[2] = 4.37273164546e-5;
    m[3] = 5.17759138449e-5;
    m[4] = 7.6923076923076926e-9;
    m[5] = 1.00000597682;
    for (i__ = 1; i__ <= 5; ++i__) {
        i1 = (i__ - 1) * 3 + 1;
        for (j = i__ + 1; j <= 6; ++j) {
            j1 = (j - 1) * 3 + 1;
/* Computing 2nd power */
            d__1 = y[i1 - 1] - y[j1 - 1];
/* Computing 2nd power */
            d__2 = y[i1] - y[j1];
/* Computing 2nd power */
            d__3 = y[i1 + 1] - y[j1 + 1];
            d__[i__ + j * 6 - 7] = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 *
                     d__3);
            d__[j + i__ * 6 - 7] = d__[i__ + j * 6 - 7];
        }
    }

    for (i__ = 1; i__ <= 6; ++i__) {
        i1 = ndim + (i__ - 1) * 3 + 1;
/* Computing 2nd power */
        d__1 = y[i1 - 1];
/* Computing 2nd power */
        d__2 = y[i1];
/* Computing 2nd power */
        d__3 = y[i1 + 1];
        rslt += m[i__ - 1] * (d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    }
    rslt /= 2.;
    pot = 0.;
    for (i__ = 2; i__ <= 6; ++i__) {
        i__1 = i__ - 1;
        for (j = 1; j <= i__1; ++j) {
            pot += m[i__ - 1] * m[j - 1] / d__[i__ + j * 6 - 7];
        }
    }
    rslt -= k * pot;
    return rslt;
} /* outer_item::hamiltonian */