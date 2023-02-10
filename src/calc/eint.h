#ifndef EINT_H
#define EINT_H

#include <stdlib.h>

#include "cint.h"

#include "tpdef.h"

//#define eri(i, j, l, m) eri[(i)+((j)+((l)+(m)*sum_di)*sum_di)*sum_di]

typedef CACHE_SIZE_T CINTintf(double *out, FINT *dims, FINT *shls,
                                          FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
                                          CINTOpt *opt, double *cache);
typedef void CINToptf(CINTOpt **opt,
                                   FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);

int get_eint(Molecular *I_atoms, Libcint *I_libcint, Integral *O_inte);

#endif //EINT_H
