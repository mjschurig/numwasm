using namespace std;
#include <iostream>
#include <math.h>

#include "task.h"

class harmonic_item : public tasklist_item {
public:
  harmonic_item(const string _name) :
        tasklist_item(_name) {}
  virtual informer info;
  virtual setupper setup;
  virtual EoM_eval EoM;
  virtual hamiltoner hamil;
};

harmonic_item my_harmonic_item("harmonic-oscillator");

void harmonic_item::info(FLP *tbegin, FLP *tend,
        FLP* stepsize, INT *ndim){
    *ndim = 1;
    *tbegin = 0.;
    *tend = M_PIl * 2;
    *stepsize = 0.01;
}

///////////////////////////
// set up the initial conditions
//
// qex, pex are the values to be *expected* at the end of the run
//
void harmonic_item::setup(const FLP /* tbegin not used */,
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
    q_X1[1] = 0.;
    p_X1[1] = 1.;
    qex_X1[1] = 0.;
    pex_X1[1] = 1.;
} /* harmonic_setup */

//////////////////////
// evaluate the equation of motion
//
void harmonic_item::EoM(const INT  /* ndim not used */,
        const FLP /* tnow not used */, FLP *q,
	FLP *f, FLP *rpar, INT *ipar)
{
    ++eom_evals;
    f[0] = -q[0];
} /* harmonic_EoM */

////////////////////////////////
// evaluate the energy
//
FLP harmonic_item::hamil(const INT ndim,
        const FLP tnow, FLP *q,
	FLP *p, FLP *rpar, INT *ipar)
{
    FLP d__1, d__2;

    d__1 = p[0];
    d__2 = q[0];
    return (d__1 * d__1 + d__2 * d__2) * .5;
} /* harmonic_hamiltonian */