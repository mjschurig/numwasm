using namespace std;
#include <iostream>
#include <math.h>

#include "task.h"

class kepler_item : public tasklist_item {
public:
  kepler_item(const string _name) :
        tasklist_item(_name){}
  virtual informer info;
  virtual setupper setup;
  virtual EoM_eval EoM;
  virtual hamiltoner hamil;
};

kepler_item my_kepler_item("kepler-orbit");

void kepler_item::info(FLP *tbegin, FLP *tend,
           FLP* stepsize, INT *ndim){
    *ndim = 2;
    *tbegin = 0.;
    *tend = M_PIl * 2;
    *stepsize = .01;
    aux = "0.6";                 // eccentricity
}

///////////////////////////
// set up the initial conditions
//
// qex, pex are the values to be *expected* at the end of the run
//
void kepler_item::setup(const FLP /* tbegin not used */,
        const FLP /* tend not used */,
	INT /* ndim not used */,
        FLP *q, FLP *p, FLP *qex, FLP
	*pex, FLP *rpar, INT *ipar)
{
    FLP eccent(atof(aux.c_str()));

    q[0] = 1 - eccent;
    q[1] = 0.;
    p[0] = 0.;
    p[1] = sqrt((eccent + 1) / (1 - eccent));
    qex[0] = 1 - eccent;
    qex[1] = 0.;
    pex[0] = 0.;
    pex[1] = sqrt((eccent + 1) / (1 - eccent));
} /* kepler_setup */

//////////////////////
// evaluate the equation of motion
//
void kepler_item::EoM(const INT /* ndim not used */,
        FLP /* tnow not used */, FLP *q,
	FLP *f, FLP *rpar, INT *ipar)
{
    ++eom_evals;
    FLP radius = q[0] * q[0] + q[1] * q[1];
    radius *= sqrt(radius);
    f[0] = -q[0] / radius;
    f[1] = -q[1] / radius;
} /* kepler_EoM */


////////////////////////////////
// evaluate the energy
//
FLP kepler_item::hamil(const INT ndim, const FLP tnow, FLP *q,
	FLP *p, FLP *rpar, INT *ipar)
{
    /* System generated locals */
    FLP d__1, d__2, d__3, d__4;

    /* Parameter adjustments */
    FLP* p_X1(p-1);
    FLP* q_X1(q-1);

    /* Function Bdy */
    d__1 = p_X1[1];
    d__2 = p_X1[2];
    d__3 = q_X1[1];
    d__4 = q_X1[2];
    return (d__1 * d__1 + d__2 * d__2) * .5
        - 1 / sqrt(d__3 * d__3 + d__4 * d__4);
} /* kepler_hamiltonian */