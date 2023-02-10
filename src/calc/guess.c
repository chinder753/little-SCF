#include "guess.h"

#include <math.h>

#include "diyblas.h"

int init_hcore(const Integral * const I_inte, double *const O_z){
    double t_hcore[I_inte->lenh], t_s[I_inte->lenh], w[I_inte->sum_di];

    cblas_dcopy(I_inte->lenh, I_inte->hcore, 1, t_hcore, 1);
    cblas_dcopy(I_inte->lenh, I_inte->overlap, 1, t_s, 1);

    LAPACKE_dspgvd(LAPACK_COL_MAJOR, 1
            , 'V', 'U'
            , I_inte->sum_di
            , t_hcore, t_s
            , w
            , O_z, I_inte->sum_di);
}



int init_ohf(const int I_sum_di, const int I_len, const int I_spin, const int I_electron, const int I_state
            , const double * const I_z_hcore
            , const double * const O_z_c, const double * const O_z_o
            ){

}



int init_uhf(const int I_sum_di, const int I_len, const int I_homo
            , double * const IO_z_a
            , double * const O_z_b
            ){
    double temp[I_sum_di], N;
    int homo = I_homo * I_sum_di, lumo = (I_homo + 1) * I_sum_di;
    // beta = alpha = hcore
    cblas_dcopy(I_len, IO_z_a, 1, O_z_b, 1);
    // alpha_homo = homo + lumo
    cblas_daxpy(I_sum_di, 1, &IO_z_a[lumo], 1, &IO_z_a[homo], 1);
    cblas_dscal(I_sum_di, sqrt(2.0), &IO_z_a[homo], 1);
    // alpha_lumo = lumo - homo
    cblas_daxpy(I_sum_di, -1, &O_z_b[homo], 1, &IO_z_a[lumo], 1);
    cblas_dscal(I_sum_di, sqrt(2.0), &IO_z_a[lumo], 1);
    // beta_homo = -alpha_lumo
    cblas_dcopy(I_sum_di, &IO_z_a[lumo], 1, &O_z_b[homo], 1);
    cblas_dscal(I_sum_di, -1, &O_z_b[homo], 1);
    // beta_lumo = -alpha_homo
    cblas_dcopy(I_sum_di, &IO_z_a[homo], 1, &O_z_b[lumo], 1);
    cblas_dscal(I_sum_di, -1, &O_z_b[lumo], 1);
}