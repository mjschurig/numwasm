INT gni_irk2(const INT ndim, task_item* eom, INT *nstep,
         const FLP *tbegin, FLP *p, FLP *q, FLP *tend,
         INT *meth, FLP *rpar, INT *ipar);

INT startb(task_item* task,
        const INT ndim, const FLP tnow, FLP *p,
        FLP *q, INT *ns, INT ndgl, FLP *fs,
	FLP *ps, FLP *zq, INT nsd, FLP *e,
	FLP *yh, INT *nm, FLP *sm, FLP *am,
	FLP *f, FLP *c__, FLP *rpar, INT *ipar);

INT rknite(task_item* task,
        const INT ndim, INT *ns, const FLP tnow,
	FLP *q, FLP *p, INT nsd, FLP *aa,
	FLP *c__, INT ndgl, FLP *qq, FLP *zq,
	FLP *f, FLP *dyno, FLP *rpar, INT *ipar);