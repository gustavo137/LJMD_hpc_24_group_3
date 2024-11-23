#ifndef WALLCLOCK_FUNC
#define WALLCLOCK_FUNC
#include <sys/time.h>
double wallclock();
#endif

#ifndef AZZERO_FUNC
#define AZZERO_FUNC
void azzero(double *d, const int n);
#endif

#ifndef PBC_FUNC
#define PBC_FUNC
double pbc(double x, const double boxby2);
#endif
