using namespace std;
#include <iostream>
#include <math.h>

#include "task.h"

class pendulum_item : public tasklist_item {
public:
  pendulum_item(const string _name) :
        tasklist_item(_name) {}
  virtual informer info;
  virtual setupper setup;
  virtual EoM_eval EoM;
  virtual hamiltoner hamil;
};

pendulum_item my_pendulum_item("pendulum");

void pendulum_item::info(FLP *tbegin, FLP *tend,
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
void pendulum_item::setup(const FLP /* tbegin not used */,
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
    qex_X1[1] = -.443944662290259;
    pex_X1[1] = .897846803312944;
} /* pendulum_setup */

//////////////////////
// evaluate the equation of motion
//
void pendulum_item::EoM(const INT /* ndim not used */,
        FLP /* tnow not used */, FLP *q,
	FLP *f, FLP *rpar, INT *ipar)
{
    ++eom_evals;
    f[0] = -sin(q[0]);
} /* pendulum_eom */


////////////////////////////////
// evaluate the hamiltonian
//
FLP pendulum_item::hamil(const INT ndim,
        const FLP /* tnow not used */, FLP *q,
	FLP *p, FLP *rpar, INT *ipar)
{
    FLP d__1 = p[0];
    return d__1 * d__1 * .5 - cos(q[0]);

} /* pendulum_hamiltonian */