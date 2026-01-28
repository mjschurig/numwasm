using namespace std;
#include <iostream>
#include <iomanip>
#include <math.h>

#include "task.h"

const FLP s43(sqrt(4.0L/3.0L));
const FLP asymp(exp(M_PIl * sqrt(2.0L)));

FLP partx(const FLP ttt){
  return ttt*sqrt(ttt * 4.L / 3.L);
}

FLP partv(const FLP ttt){
  return sqrt(3. * ttt);
}

class ratio_item : public tasklist_item{
public:
  ratio_item(const string _name) :
        tasklist_item(_name){}
  virtual solfixer solfix;
  virtual informer info;
  virtual setupper setup;
  virtual EoM_eval EoM;
  virtual hamiltoner hamil;
  virtual statter statistics;

/////////////////////////////////////////////
// some variables to support the output routine:

FLP last_t;
FLP last_p;
FLP last_q;
FLP last_z;
FLP last_cross_t;

double maxerr_p;
double maxerr_q;

int partmode;

};

ratio_item my_ratio_item("ratio-gun");

void ratio_item::info(FLP *tbegin, FLP *tend,
        FLP* stepsize, INT *ndim){
    *ndim = 1;
    *tbegin = 0;
    *tend = *tbegin * 10.;
    if (! *tend) *tend = 0.001;
// 100 steps per decade is more than enough;
// 1000 steps per decade is waaay more than enough,
// but may be helpful until we get smarter about interpolation
// in the zero-crossing-finder.
    *stepsize = *tend / 1000.;
    aux = -1.;                  // initial velocity
    odt = 0.1;
    verbosity = 2;
}

///////////////////////////
// set up the initial conditions
//
// qex, pex are the values to be *expected* at the end of the run
//
void ratio_item::setup(const FLP tbegin,
        const FLP /* tend not used */,
	INT /* ndim not used */,
        FLP *q, FLP *p, FLP *qex, FLP
	*pex, FLP *rpar, INT *ipar)
{
 last_t = 0;
 last_p = 0;
 last_q = 0;
 last_z = 0;
 last_cross_t = 0;
 partmode = aux == "part";
 maxerr_p = 0;
 maxerr_q = 0;

  if (partmode){
    q[0] = partx(tbegin);
    p[0] = partv(tbegin);
    qex[0] = 0.;
    pex[0] = 0.;
  } else {
    FLP v0(atof(aux.c_str()));
    q[0] = 1. + v0 * tbegin;
    p[0] = v0;
    qex[0] = 0.;
    pex[0] = 0.;
  }
} /* ratio_gun_setup */

//////////////////////
// evaluate the equation of motion
//
void ratio_item::EoM(const INT /* ndim not used */,
        const FLP tnow, FLP *q,
	FLP *f, FLP *rpar, INT *ipar)
{
    ++eom_evals;              // count number of times we've been called
    f[0] = tnow / q[0];
} /* ratio_gun_EoM */


////////////////////////////////
// evaluate the energy
//
FLP ratio_item::hamil(const INT ndim,
        const FLP tnow, FLP *q,
	FLP *p, FLP *rpar, INT *ipar)
{
  return 0.5*p[0]*p[0] - tnow * log(q[0]);
} /* ratio_gun_hamiltonian */

FLP zfunc(const FLP ttt, const FLP qqq){
  return (qqq - partx(ttt)) / sqrt(ttt);
}

/////////////////
// a little output utility:
//
void dump_record(const FLP intfact,
        const FLP t0, const FLP p0, const FLP q0,
        const FLP t1, const FLP p1, const FLP q1,
        const int flag = 0){

  FLP z0 = zfunc(t0, q0);
  FLP z1 = zfunc(t1, q1);
  FLP ttt = (1.-intfact)*t0 + intfact*t1;
  FLP ppp = (1.-intfact)*p0 + intfact*p1;
  FLP qqq = (1.-intfact)*q0 + intfact*q1;
  FLP zzz = (1.-intfact)*z0 + intfact*z1;

  cout << setw(16) << ttt << ", "
       << setw(16) << qqq << ",,"
       << setw(16) << ppp << ",,"
       << setw(16) << ttt << ", "
       << setw(16) << zzz << ","
       << setw(16) << qqq/partx(ttt) - 1. << ","
       << setw(16) << ppp/partv(ttt) - 1. << ","
       << (flag ? " T" : "");
  cout << endl;
}

INT ratio_item::solfix(INT nr,
        const FLP t_old, const FLP tnow,
	FLP *p, FLP *q, const INT ndim,
	FLP * rpar, INT *ipar)
{

// Collect and store statistics
#if 0
    FLP energy;
    energy = hamil(ndim, tnow, q, p, rpar, ipar);
    if (!raw_points) {
	e0 = energy;       // save the initial value
	delta_e_max = 0.;
    } else {
        delta_e_max = max(delta_e_max, fabs(e0-energy));
    }
#endif

    maxerr_q = max(maxerr_q, fabs(double(q[0] / partx(tnow) -1.)));
    maxerr_p = max(maxerr_p, fabs(double(p[0] / partv(tnow) -1.)));

// That's all the statistics.

    if (!raw_points && verbosity > 0) {
      // print column headers:
      cout << "t, x,, v,, t, z, crossing" << endl;
    }
    raw_points++;               // count # of times we have been called
    cout.precision(10);

    FLP zval = zfunc(tnow, q[0]);
    FLP cross_t;
    if (!partmode && last_t && (last_z * zval < 0.)) {
      FLP intfact = (0. - last_z) / (zval - last_z);
      cross_t = (1.-intfact)*last_t + intfact*tnow;
      if ((verbosity & 1) && last_cross_t) {
        // print analysis of zero crossing:
        FLP ratio = cross_t / last_cross_t;
        /*{*/
        cout << "X   ,,,,,,,}, "
             << setw(12) << cross_t
             << " / " << setw(12) << last_cross_t
             << " = " << setw(12) << ratio
             << "  =>  " << (ratio / asymp - 1.0) * 1e6 << " ppm"
             << endl;
        // same as above, but rounded to 6 figures:
        int save = cout.precision(6);
        /*{*/
        cout << "Y   ,,,,,,,}, "
             << setw(12) << cross_t
             << " / " << setw(12) << last_cross_t
             << " = " << setw(12) << ratio
             << "  =>  " << (ratio / asymp - 1.0) * 1e6 << " ppm"
             << endl;
        cout.precision(save);
      }
      if ((verbosity & 3)){
        // special interpolated record, showing the zero-crossing
        // printed without regard to odt
        dump_record(intfact, last_t, last_p, last_q, tnow, p[0], q[0], 1);
      }
      last_cross_t = cross_t;
    }

    if ((verbosity & 2) && tnow > good_out_t) {
      // normal output record, not interpolated
      // spacing controlled by odt
      dump_record(0, tnow, p[0], q[0], tnow, p[0], q[0]);
      good_out_t = tnow * (1. + odt);
    }

    last_t = tnow;
    last_z = zval;
    last_p = p[0];
    last_q = q[0];
    return 0;
} /* ratio_item::solfix */


FLP ratio_item::statistics(
        const INT ndim, const INT nstep,
        const FLP tbegin, const FLP tend,
        FLP *q, FLP *p,
        FLP *qex, FLP *pex,
        const FLP wall_time
){
    cout.precision(10);

    if (partmode) {
#if 1 // KLUDGE :  should probably replace this with a call to solfix(...)
        const int wid(16);
        cout << setw(wid) << last_t
             << ", " << setw(wid) << last_q
             << ", " << setw(wid) << last_p
             << endl;
        cout << setw(wid) << tend
             << ", " << setw(wid) << q[0]
             << ", " << setw(wid) << p[0]
             << endl;
#endif
      double fake(1e-18);
      cout << "** t-t/nstep:, " << 1000. * (tend-tbegin)/nstep
           << ", x/x-1:, " << max(fake, maxerr_q)
           << ", v/v-1:, " << max(fake, maxerr_p)
           << endl;
    }

// calculate error by comparing to expected values
    FLP err(0);
    for (int ii = 1; ii <= ndim; ++ii) {
	FLP dq = q[ii - 1] - qex[ii - 1];
	FLP dp = p[ii - 1] - pex[ii - 1];
	err = err + dq * dq + dp * dp;
    }
    err = sqrt(err);
    cout << "**"
         << "  eom_evals: " << eom_evals
         << "  wall_time: " << setw(14) << wall_time
         << endl;

    return 0;
}

/*

New Data:

Long doubles:

** 1/nstep:, 0.25, x/x-1:, 2.253606959e-08, v/v-1:, 6.759167681e-09
** 1/nstep:, 0.1666666667, x/x-1:, 1.23943084e-09, v/v-1:, 2.940423616e-10
** 1/nstep:, 0.125, x/x-1:, 1.203233638e-10, v/v-1:, 2.537343903e-11
** 1/nstep:, 0.08333333333, x/x-1:, 3.114395623e-12, v/v-1:, 6.271194501e-13
** 1/nstep:, 0.05882352941, x/x-1:, 2.068737976e-13, v/v-1:, 3.492506848e-13
** 1/nstep:, 0.03448275862, x/x-1:, 5.527096792e-13, v/v-1:, 9.192545813e-13
** 1/nstep:, 0.02272727273, x/x-1:, 6.146754113e-14, v/v-1:, 1.004025422e-13
** 1/nstep:, 0.01694915254, x/x-1:, 2.184999686e-13, v/v-1:, 3.928498694e-13
** 1/nstep:, 0.01123595506, x/x-1:, 3.362611838e-13, v/v-1:, 5.07689214e-13
** 1/nstep:, 0.0078125, x/x-1:, 7.053161223e-13, v/v-1:, 1.070319939e-12
** 1/nstep:, 0.005586592179, x/x-1:, 9.563866084e-13, v/v-1:, 1.452146346e-12
** 1/nstep:, 0.003344481605, x/x-1:, 3.785350939e-14, v/v-1:, 1.307076192e-13
** 1/nstep:, 0.002222222222, x/x-1:, 1.272978034e-14, v/v-1:, 1.820510973e-14
** 1/nstep:, 0.001669449082, x/x-1:, 2.190576279e-15, v/v-1:, 3.026767205e-15
** 1/nstep:, 0.001111111111, x/x-1:, 5.485520892e-16, v/v-1:, 6.819089564e-16
** 1/nstep:, 0.0007782101167, x/x-1:, 2.919756451e-16, v/v-1:, 2.125036258e-16
** 1/nstep:, 0.0005555555556, x/x-1:, 9.622185861e-15, v/v-1:, 1.826696346e-14
** 1/nstep:, 0.0003333333333, x/x-1:, 4.034352062e-13, v/v-1:, 5.664641933e-13
** 1/nstep:, 0.0002222222222, x/x-1:, 1.367180457e-12, v/v-1:, 1.910144135e-12
** 1/nstep:, 0.0001666666667, x/x-1:, 2.121085967e-12, v/v-1:, 2.963764759e-12
** 1/nstep:, 0.0001111111111, x/x-1:, 6.262440653e-13, v/v-1:, 8.747994733e-13
** 1/nstep:, 7.777864198e-05, x/x-1:, 2.154965659e-13, v/v-1:, 3.009750652e-13
** 1/nstep:, 5.555555556e-05, x/x-1:, 9.242812678e-14, v/v-1:, 1.178092996e-13

plain doubles:

** 1/nstep:, 0.25, x/x-1:, 2.253606968e-08, v/v-1:, 6.75916767e-09
** 1/nstep:, 0.1666666667, x/x-1:, 1.23943078e-09, v/v-1:, 2.94042346e-10
** 1/nstep:, 0.1111111111, x/x-1:, 4.343270188e-11, v/v-1:, 8.757439218e-12
** 1/nstep:, 0.08333333333, x/x-1:, 3.114508651e-12, v/v-1:, 6.272760089e-13
** 1/nstep:, 0.05555555556, x/x-1:, 1.094235813e-12, v/v-1:, 1.951327988e-12
** 1/nstep:, 0.03333333333, x/x-1:, 4.654054919e-13, v/v-1:, 7.736034036e-13
** 1/nstep:, 0.02272727273, x/x-1:, 6.150635556e-14, v/v-1:, 1.003641614e-13
** 1/nstep:, 0.01666666667, x/x-1:, 2.081668171e-13, v/v-1:, 3.470557175e-13
** 1/nstep:, 0.01123595506, x/x-1:, 3.360645096e-13, v/v-1:, 5.077049892e-13
** 1/nstep:, 0.0078125, x/x-1:, 7.043254868e-13, v/v-1:, 1.070588063e-12
** 1/nstep:, 0.005586592179, x/x-1:, 9.570122472e-13, v/v-1:, 1.451949672e-12
** 1/nstep:, 0.003333333333, x/x-1:, 3.996802889e-14, v/v-1:, 1.28119737e-13
** 1/nstep:, 0.002222222222, x/x-1:, 2.076117056e-14, v/v-1:, 1.765254609e-14
** 1/nstep:, 0.001666666667, x/x-1:, 2.575717417e-14, v/v-1:, 5.995204333e-15
** 1/nstep:, 0.001111111111, x/x-1:, 3.552713679e-14, v/v-1:, 8.215650382e-15
** 1/nstep:, 0.0007782101167, x/x-1:, 5.018208071e-14, v/v-1:, 1.421085472e-14
** 1/nstep:, 0.0005555555556, x/x-1:, 8.015810238e-14, v/v-1:, 3.17523785e-14
** 1/nstep:, 0.0003333333333, x/x-1:, 5.013767179e-13, v/v-1:, 5.824229987e-13
** 1/nstep:, 0.0002222222222, x/x-1:, 1.341593503e-12, v/v-1:, 1.906030889e-12
** 1/nstep:, 0.0001666666667, x/x-1:, 2.185140957e-12, v/v-1:, 2.960076628e-12
** 1/nstep:, 0.0001111111111, x/x-1:, 7.618350395e-13, v/v-1:, 8.708589405e-13
** 1/nstep:, 7.777864198e-05, x/x-1:, 5.280220705e-13, v/v-1:, 3.212985433e-13
** 1/nstep:, 5.555555556e-05, x/x-1:, 5.648814749e-13, v/v-1:, 1.663114091e-13

*/

/*

For taking data:

for step in $( ratio-step 2e-4 1e-9 6 ) ; do
 ./dr_irk2 -task rat -aux part -tb 0.0001 -step $step | fgrep /nstep
done

*/