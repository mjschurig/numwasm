//     SOLFIX      routine that prints
//                 numerical solution during integration.
//
//                    int solfix (nr,t_old,tnow,p,q,ndim,rpar,ipar)
//
//                 solfix furnishes the solution "q,p" at the nr-th
//                    grid-point "tnow" (initial value for nr=0).
//                 "t_old" is the previous grid-point (if any).
//                 if the return value is <0, the integration will
//                 be aborted.

using namespace std;
#include <stdlib.h>             /* for exit() */

#include <iostream>
#include <sstream>
#include "task.h"

#include <string>
#include <iomanip>

#include <sys/timeb.h>
#include <math.h>

//////////////////////
// Print an output record
//
// we use task_item::e0 to store the initial value of the hamiltonian
// we use task_iten::delta_e_max to keep a running maximum of the
//  excursion of the hamiltonian from its initial value.
INT task_item::solfix(INT nr,
        const FLP t_old, const FLP tnow,
	FLP *p, FLP *q, const INT ndim,
        FLP *rpar, INT *ipar)
{

    int age = raw_points++;     // count how many times we've been called
// Collect and store statistics on how the Hamiltonian is doing:
    FLP energy;
    energy = hamil(ndim, tnow, q, p, rpar, ipar);
    if (nr == 0) {
	e0 = energy;       // save the initial value
	delta_e_max = 0.;
    } else {
/* Computing MAX */
        delta_e_max = max(delta_e_max, fabs(double(e0-energy)));
    }
// That's all the statistics.

// Now do some output, if desired:
    if (verbosity <= 0) return 0;
    int fwid(16);
    string sep(", ");

    if (!age) {
      ostringstream fmt;
      // print column headers:
      cout << setw(fwid) << "t";
      for (INT ii = 0; ii < ndim; ii++) {
        fmt << "q[" << ii << "]";
        cout << sep << setw(fwid) << fmt.str();
        fmt.str("");
      }
      for (INT ii = 0; ii < ndim; ii++) {
        fmt << "p[" << ii << "]";
        cout << sep << setw(fwid) << fmt.str();
        fmt.str("");
      }
      cout << endl;
    }

    if (verbosity > 1 || tnow >= good_out_t) {
      good_out_t = tnow + odt;
      cout.precision(10);
      cout << setw(fwid) << tnow;
      for (INT ii = 0; ii < ndim; ii++) {
        cout << sep << setw(fwid) << q[ii];
      }
      for (INT ii = 0; ii < ndim; ii++) {
        cout << sep << setw(fwid) << p[ii];
      }
      cout << endl;
    }
    return 0;
} /* solfix */


/////////////////////////////////
// trivial place-holder
//
FLP task_item::hamil(const INT ndim,
        const FLP tnow, FLP *q,
	FLP *p, FLP *rpar, INT *ipar){
  return 0;
}


FLP task_item::statistics(
        const INT ndim, const INT nstep,
        const FLP /* tbegin not used */, const FLP /* tend not used */,
        FLP *q, FLP *p,
        FLP *qex, FLP *pex,
        const FLP wall_time
){
    cout.precision(10);

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
         << "  err: "  << err
         << "  max_hamil_dev: " << delta_e_max
         << endl;
    return 0;
}