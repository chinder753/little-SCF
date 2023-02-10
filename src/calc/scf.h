#ifndef HF_H
#define HF_H

#include "tpdef.h"

double rhf(const int I_loop, const double I_D_error, const int I_occ
        , const Integral * const I_inte
        , double * const IO_z
        , double * const O_w
);

double uhf(const int I_loop, const int I_diis, double const I_D_error
           , const int I_sum_di, const int I_len, const int I_lenh
           , const int * const I_alpha, const int * const I_beta
           , const double * const I_s, const double * const I_hcore, const double * const I_eri
           , double * const IO_z_a, double * const IO_z_b
           );

double rmp2(const int I_occ
        , const Integral * const I_inte
        , const double * const I_z, const double * const I_w
);

int rcis(const int I_electron
        , const Integral * const I_inte
        , const double * const I_z, const double * const I_w
        , double * const O_excite
);

#endif
