#ifndef GUESS_H
#define GUESS_H

#include "tpdef.h"

int init_hcore(const Integral * const I_inte, double *const O_z);

int init_uhf(const int I_sum_di, const int I_len, const int I_homo
        , double * const IO_z_a
        , double * const O_z_b
        );

#endif
