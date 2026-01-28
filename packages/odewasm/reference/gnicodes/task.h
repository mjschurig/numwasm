using namespace std;
#include <sys/types.h>
#include "mytypes.h"
#include <string>
#include <stdlib.h>     /* for exit */
#include <math.h>       /* for sqrtl */

typedef void informer(FLP *tbegin, FLP *tend,
        FLP* stepsize, INT *ndim);

typedef void setupper(const FLP tbegin, const FLP tend,
	const INT ndim,
        FLP *q, FLP *p,
        FLP *qex, FLP *pex,
        FLP *rpar, INT *ipar);

typedef void EoM_eval(const INT ndim,
        const FLP tnow, FLP *q,
	FLP *f, FLP *rpar, INT *ipar);

typedef FLP hamiltoner(const INT ndim,
        const FLP tnow, FLP *q,
	FLP *p, FLP *rpar, INT *ipar);

typedef INT solfixer(INT nr,
        const FLP t_old, const FLP tnow,
	FLP *p, FLP *q, const INT ndim,
        FLP *rpar, INT *ipar);

typedef FLP statter(
        const INT ndim, const INT nstep,
        const FLP tbegin, const FLP tend,
        FLP *q, FLP *p,
        FLP *qex, FLP *pex,
        const FLP wall_time);

class task_item{
public:
  int verbosity;
// Collect some statistics
  FLP e0;            // initial energy (used by solfix)
  double delta_e_max;   // max deviation from its initial value
  int eom_evals;        // number of calls to equation-of-motion function
// generic adjustable parameter
// ... eccentricity (for kepler task)
// ... initial velocity (for ratio-gun task)
  string aux;
// for controlling the density of output points:
  FLP odt;
  FLP good_out_t;
  int raw_points;       // number of times solfix has been called
              // (possibly less than the number of recorts output)

// constructor:
  task_item() : verbosity(1), e0(0), delta_e_max(0),
        eom_evals(0), aux(), odt(0), good_out_t(0),
        raw_points(0) {}

// methods:
  virtual solfixer solfix;
  virtual statter statistics;
// You don't have to write a hamil method, but if you don't,
// the default just returns zero:
  virtual hamiltoner hamil;
  virtual informer info = 0;
  virtual setupper setup = 0;
  virtual EoM_eval EoM = 0;
};

// Pull in ye olde linked-list class:
#include "list-item.h"
extern list_item* menu_root;

class tasklist_item : public task_item, public list_item{
public:
  tasklist_item(const string _name) :
    list_item(menu_root, _name)
  {}
};

inline long double sqrt(long double foo) {
  return sqrtl(foo);
}

inline long double sin(long double foo) {
  return sinl(foo);
}

inline long double cos(long double foo) {
  return cosl(foo);
}