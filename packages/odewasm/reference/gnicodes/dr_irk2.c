using namespace std;
#include <stdlib.h>             /* for exit() */

#include <iostream>
#include "task.h"
#include "gni_irk2.h"

#include <string>
#include <iomanip>

#include <sys/timeb.h>

void dump_menu();       // forward reference

void usage() {
    cout <<
"dr_irk2 [options]\n"
"\n"
"Options include\n"
"  -task tttt   Use the equations of motion for task tttt.\n"
"               See below for a list.\n"
"  -step dd.dd  Set stepsize.\n"
"  -tbegin tt   Set initial time.\n"
"  -method d    Choose Gaussian method of order 2*d.\n"
"  -verbosity v Determine verbosity of output; -v 1 = normal.\n"
"  -quiet       Same as -verb 0;  print only summary.\n"
"  -iterate N   Iterate N times, advancing the end-time 10X each round\n"
"  -odt tt.tt   Set minimum spacing of output points.\n"
"  -aux aa.aa   Set auxiliary parameter:\n"
"                ... eccentricity for kepler-orbit\n"
"                ... initial velocity for ratio-gun\n"
"  -h           Print this message (and exit).\n"
"\n"
"The following tasks are implemented in this version:\n"
;
  dump_menu();
  cout <<
"\n"
"Task names can be abbreviated as much as you dare.\n"
;
}

list_item* menu_root(0);
static list_guard whatever(menu_root, "guard");

void dump_menu(){
    tasklist_item* item((tasklist_item*)menu_root->next);
    FLP tbegin, tend, stepsize;
    int ndim;
    cout << "   " << setw(22) <<  "default:";
    cout << "  tbegin     tend stepsize   nsteps ndim      aux     odt"
         << endl;
    while (item) {
      cout << "   " << setw(22) << left << item->name;
      cout << right;
      item->info(&tbegin, &tend, &stepsize, &ndim);
      INT nstep = (tend - tbegin + stepsize/2.) / stepsize;
      cout << setw(8) << tbegin << " "
        << setw(8) << tend << " "
        << setw(8) << stepsize << " "
        << setw(8) << nstep << " "
        << setw(4) << ndim << " "
        << setw(8) << item->aux
        << setw(8) << item->odt
        << endl;
      item = (tasklist_item*)item->next;
    }
}

//////////////////////////////////////////////////////////////////////
// return time with millisecond resolution
// seconds since the epoch
FLP jtime() {
  struct timeb mytime;
  ftime(&mytime);
  return mytime.time + mytime.millitm / 1000.;
}

////////////////
// little utility to help with argument parsing:
//
int prefix(const string shorter, const string longer){
  return shorter == longer.substr(0, shorter.length());
}

///////////////////////////
// Main program
//
int main(const INT _argc, const char* _argv[])
{
    tasklist_item* thisItem((tasklist_item*)menu_root->next);
    if (!thisItem) {
      cerr << "No tasks implemented in this version!?!?" << endl;
      exit(1);
    }
    INT meth(6), ndim;
    FLP stepsize(0);
    string task_name("kepler");
    FLP tbegin, tend;
    INT maxiter(1);

// argument processing:
    INT argc(_argc);
    const char **argv(_argv);
    string progname(*argv); argv++; argc--;
    int did_info(0);
    while (argc) {
      string opt(*argv); argv++; argc--;
      if (opt.substr(0,2) == "--") opt = opt.substr(1);
      if (prefix(opt, "-help")) {
        usage();
        exit(0);
      }
      if (prefix(opt, "-task")) {
        if (did_info > 0) {
          cerr << "The -task option should appear only once." << endl;
          exit(1);
        }
        if (did_info < 0) {
          cerr << "The -task option should come first." << endl;
          exit(1);
        }
        if (!argc) {
          cerr << "Option -task requires an argument" << endl;
          exit(1);
        }
        task_name = *argv; argv++; argc--;
        while (thisItem) {
          if (prefix(task_name, thisItem->name)) goto got1;
          thisItem = (tasklist_item*)thisItem->next;
        }
// Here if we fell off the end of the task-menu:
        cerr << "Task '" << task_name
            << "' is not on the menu;" << endl;
        cerr << "  try 'dr_irk2 -help' for a list of possibilities." << endl;
        exit(1);
got1:;;;
// Here with a valid task.
        thisItem->info(&tbegin, &tend, &stepsize, &ndim);
        did_info = 1;
        continue;
      }
// Here if we have an option other than -h and -task.
// Need a task.  Default to first thing on menu list:
      if (!did_info){
        thisItem->info(&tbegin, &tend, &stepsize, &ndim);
        did_info = -1;
      }
      if (prefix(opt, "-quiet")) {
        thisItem->verbosity = 0;
      } else if (prefix(opt, "-verbosity")) {
        if (!argc) {
          cerr << "Option -verbosity requires an argument" << endl;
          exit(1);
        }
        thisItem->verbosity = atoi(*argv); argv++; argc--;
      } else if (prefix(opt, "-method")) {
        if (!argc) {
          cerr << "Option -method requires an argument" << endl;
          exit(1);
        }
        meth = atoi(*argv); argv++; argc--;
      } else if (prefix(opt, "-iterate")) {
        if (!argc) {
          cerr << "Option -iterate requires an argument" << endl;
          exit(1);
        }
        maxiter = atoi(*argv); argv++; argc--;
      } else if (prefix(opt, "-stepsize")) {
        if (!argc) {
          cerr << "Option -stepsize requires an argument" << endl;
          exit(1);
        }
        stepsize = atof(*argv); argv++; argc--;
      } else if (prefix(opt, "-tbegin")) {
        if (!argc) {
          cerr << "Option -tbegin requires an argument" << endl;
          exit(1);
        }
        tbegin = atof(*argv); argv++; argc--;
      } else if (prefix(opt, "-tend")) {
        if (!argc) {
          cerr << "Option -tend requires an argument" << endl;
          exit(1);
        }
        tend = atof(*argv); argv++; argc--;
      } else if (prefix(opt, "-aux")) {
        if (!argc) {
          cerr << "Option -aux requires an argument" << endl;
          exit(1);
        }
        thisItem->aux = *argv; argv++; argc--;
      } else if (prefix(opt, "-odt")) {
        if (!argc) {
          cerr << "Option -odt requires an argument" << endl;
          exit(1);
        }
        thisItem->odt = atof(*argv); argv++; argc--;
      } else {
        cerr << "Unrecognized option '" << opt << "'" << endl;
        exit(1);
      }
    }

// STILL no task?  Default to first thing on menu list:
    if (!did_info){
      thisItem->info(&tbegin, &tend, &stepsize, &ndim);
    }

    FLP p[ndim], q[ndim], pex[ndim], qex[ndim];
    const int par_needed(10);
    const int parsize(20);      // KLUDGE
    INT ipar[parsize];
    FLP rpar[parsize];
    FLP time0, time1;
    INT nstep;

    thisItem->setup(tbegin, tend, ndim, q, p, qex, pex, rpar, ipar);
    nstep = (INT) ((tend - tbegin) / stepsize);

    cout << "** main:: "
        << "  method: " << meth
        << "  tbegin: " << tbegin
        << "  tend: " << tend
        << "  nstep: " << nstep
        << "  stepsize: " << stepsize
        << "  iter: " << maxiter
        << endl;

    cout << "** " << thisItem->name
        << " ... "
        << "  ndim: " << ndim
        << "  aux: " << thisItem->aux
        << "  odt: " << thisItem->odt
        << endl;

    for (int ii = par_needed; ii < parsize; ii++) {
      ipar[ii] = 10001 * ii;
      rpar[ii] = 100001. * ii;
    }
/* --- CALL OF THE METHOD */
    time0 = jtime();
    for (int ii = 0; ;) {
      for (int kk = 0; kk < par_needed; ++kk) {
          rpar[kk] = 0.;
          ipar[kk] = 0;
      }
      gni_irk2(ndim, thisItem, &nstep,
          &tbegin, p, q, &tend, &meth,
          rpar, ipar);
      if (++ii >= maxiter) break;
      tbegin = tend;
      tend = 10. * tbegin;
      // keep the same number of steps per decade
    }
    time1 = jtime();
    for (int ii = par_needed; ii < parsize; ii++) {
      if (ipar[ii] != 10001 * ii) {
        cerr << "Porridge eaten: " << ii << endl;
      }
      if (rpar[ii] != 100001. * ii) {
        cerr << "Chair crushed: " << ii << endl;
      }
    }
    thisItem->statistics(ndim, nstep, tbegin, tend,
        q, p, qex, pex, time1-time0);
} /* int main() */