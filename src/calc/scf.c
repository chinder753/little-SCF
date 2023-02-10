#include "scf.h"

#include <stdio.h>
#include <stdlib.h>

#include "fock.h"
#include "diyblas.h"

#define eri(i, j, l, m) I_inte->eri[(((i)*I_inte->sum_di+(j))*I_inte->sum_di+(l))*I_inte->sum_di+(m)]


int rcis(const int I_electron
         , const Integral * const I_inte
        , const double * const I_z, const double * const I_w
        , double * O_excite
         ){
    //index
    int nstate = I_electron * (I_inte->sum_di*2 - I_electron);
    int cis_index[nstate][2];
    int index_occ[I_electron], index_vri[nstate];
    for(int i=0; i<I_electron/2; i++){
        index_occ[i] = i+1;
        index_occ[i+I_electron/2] = -i-1;
    }
    for(int i=0; i<(I_inte->sum_di*2 - I_electron)/2; i++){
        index_vri[i] = i+I_electron/2+1;
        index_vri[i+(I_inte->sum_di*2 - I_electron)/2] = -i-I_electron/2-1;
    }

    int temp_index = 0;
    for(int i=0; i<I_electron; i++){
        for(int j=0; j<I_inte->sum_di*2-I_electron; j++){
            cis_index[temp_index][0] = index_occ[i];
            cis_index[temp_index][1] = index_vri[j];
            temp_index++;
        }
    }

    // hamilton
    double cis_h[nstate][nstate];
    for(int l=0; l < nstate; l++){
        for(int m=0; m < l+1; m++) {
            int a = abs(cis_index[l][1]) - 1, b = abs(cis_index[m][1]) - 1, i = abs(cis_index[l][0]) - 1, j = abs(cis_index[m][0]) - 1;
            cis_h[l][m] = 0;
            if(cis_index[l][1] * cis_index[l][0] > 0 && cis_index[m][1] * cis_index[m][0] > 0){
                for(int lambda=0; lambda < I_inte->sum_di; lambda++) {
                    for(int mu=0; mu < I_inte->sum_di; mu++){
                        for(int sigma=0; sigma < I_inte->sum_di; sigma++) {
                            for(int nu=0; nu < I_inte->sum_di; nu++) {
                                cis_h[l][m] += I_z[a*I_inte->sum_di+lambda] * I_z[i*I_inte->sum_di+mu] * eri(lambda,mu,sigma,nu) * I_z[b*I_inte->sum_di+sigma] * I_z[j*I_inte->sum_di+nu];
                            }
                        }
                    }
                }
            }
            if(cis_index[l][1] * cis_index[m][1] > 0 && cis_index[l][0] * cis_index[m][0] > 0){
                for(int lambda=0; lambda < I_inte->sum_di; lambda++) {
                    for(int mu=0; mu < I_inte->sum_di; mu++){
                        for(int sigma=0; sigma < I_inte->sum_di; sigma++) {
                            for(int nu=0; nu < I_inte->sum_di; nu++) {
                                cis_h[l][m] -= I_z[a*I_inte->sum_di+lambda] * I_z[b*I_inte->sum_di+mu] * eri(lambda,mu,sigma,nu) * I_z[i*I_inte->sum_di+sigma] * I_z[j*I_inte->sum_di+nu];
                            }
                        }
                    }
                }
            }
            if(cis_index[l][0] == cis_index[m][0] && cis_index[l][1] == cis_index[m][1]){
                cis_h[l][m] += I_w[a] - I_w[i];
            }
        }
    }

    LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', nstate, cis_h, nstate, O_excite);
    return 0;
}



double rmp2(const int I_occ
           , const Integral * const I_inte
           , const double * const I_z, const double * const I_w
           ){
    double temp=0;
    for(int i=0; i<I_occ; i++){
        for(int j=0; j<I_occ; j++) {
            for(int a=I_occ; a<I_inte->sum_di; a++){
                for(int b=I_occ; b<I_inte->sum_di; b++) {

                    double in_1=0, in_2=0;
                    for(int mu=0; mu < I_inte->sum_di; mu++){
                        for(int nu=0; nu < I_inte->sum_di; nu++) {
                            for(int lambda=0; lambda < I_inte->sum_di; lambda++) {
                                for(int sigma=0; sigma < I_inte->sum_di; sigma++) {
                                    in_1 += I_z[i*I_inte->sum_di+lambda] * I_z[a*I_inte->sum_di+mu] * eri(lambda,mu,sigma,nu) * I_z[j*I_inte->sum_di+sigma] * I_z[b*I_inte->sum_di+nu];
                                    in_2 += I_z[i*I_inte->sum_di+lambda] * I_z[b*I_inte->sum_di+nu] * eri(lambda,nu,sigma,mu) * I_z[j*I_inte->sum_di+sigma] * I_z[a*I_inte->sum_di+mu];
                                }
                            }
                        }
                    }
                    temp += in_1*(2*in_1-in_2)/(I_w[i] + I_w[j] - I_w[a] - I_w[b]);
                }
            }
        }
    }
    return temp;
}



double uhf(const int I_loop, const int I_diis, double const I_D_error
           , const int I_sum_di, const int I_len, const int I_lenh
           , const int * const I_alpha, const int * const I_beta
           , const double * const I_s, const double * const I_hcore, const double * const I_eri
           , double * const IO_z_a, double * const IO_z_b
           ){
    int count=0;
    double E=0, E_old=1, error_a=1, error_b=1;
    double p_a[I_len], p_b[I_len], delta_p_a[I_len], delta_p_b[I_len];
    double F_a[I_len], F_b[I_len];
    double w_a[I_sum_di], w_b[I_sum_di];

    u_density(I_sum_di
            , I_alpha, I_beta
            , IO_z_a, IO_z_b
            , p_a, p_b
    );
    ufock_build(I_sum_di, I_len, I_lenh
                , I_alpha, I_beta
                , I_hcore, I_eri, p_a, p_b
                , F_a, F_b
                );


    while((error_a>I_D_error || error_b>I_D_error) && (count < I_loop)){
//        E_old = E;
//        E = 0;
        double temp_s[I_lenh];
        cblas_dcopy(I_lenh, I_s, 1, temp_s, 1);
        // solve FC=eSC
        LAPACKE_dspgvd(LAPACK_COL_MAJOR, 1
                , 'V', 'U'
                , I_sum_di
                , F_a, temp_s
                , w_a
                , IO_z_a, I_sum_di
        );
        cblas_dcopy(I_lenh, I_s, 1, temp_s, 1);
        LAPACKE_dspgvd(LAPACK_COL_MAJOR, 1
                , 'V', 'U'
                , I_sum_di
                , F_b, temp_s
                , w_b
                , IO_z_b, I_sum_di
        );
        cblas_dcopy(I_len, p_a, 1, delta_p_a, 1);
        cblas_dcopy(I_len, p_b, 1, delta_p_b, 1);

        u_density(I_sum_di
                , I_alpha, I_beta
                , IO_z_a, IO_z_b
                , p_a, p_b
        );

        ufock_build(I_sum_di, I_len, I_lenh
                    , I_alpha, I_beta
                    , I_hcore, I_eri, p_a, p_b
                    , F_a, F_b
                    );

        cblas_daxpy(I_len, -1, p_a, 1, delta_p_a, 1);
        error_a = cblas_dnrm2(I_len, delta_p_a, 1);

        cblas_daxpy(I_len, -1, p_b, 1, delta_p_b, 1);
        error_b = cblas_dnrm2(I_len, delta_p_b, 1);

        printf("loop %3d:\n  D_alpha_error:%15e\n  D_beta_error:%15e\n", count, error_a, error_b);

        count++;
    }
    if(count==I_loop){
        printf("scf not converged\n");
    }else{
        printf("converged in %d loop\n", count);
    }
    return E;
}



double rhf(const int I_loop, const double I_D_error, const int I_occ
           , const Integral * const I_inte
           , double * const IO_z
           , double * const O_w
           ){
    int count=0;
    double E=0, E_old=1, error=1;
    double p_i[I_inte->len], delta_p_i[I_inte->len], F[I_inte->lenh];

    cblas_dgemm(LAPACK_COL_MAJOR, CblasNoTrans, CblasTrans
            , I_inte->sum_di, I_inte->sum_di, I_occ
            , 2
            , IO_z, I_inte->sum_di
            , IO_z, I_inte->sum_di
            , 0
            , p_i, I_inte->sum_di
    );


    rfock(I_inte->sum_di, I_inte->len, I_inte->lenh
               , I_inte->hcore, I_inte->eri, p_i
               , F
               );

    while(((E-E_old)>1e-5 || (E-E_old)<-1e-5 || error>I_D_error) && (count<I_loop)){
        E_old = E;
        E = 0;
        double temp_s[I_inte->lenh];
        cblas_dcopy(I_inte->lenh, I_inte->overlap, 1, temp_s, 1);
        // solve FC=eSC
        LAPACKE_dspgvd(LAPACK_COL_MAJOR, 1
                , 'V', 'U'
                , I_inte->sum_di
                , F, temp_s
                , O_w
                , IO_z, I_inte->sum_di
                );
        cblas_dcopy(I_inte->len, p_i, 1, delta_p_i, 1);

        count++;
        cblas_dgemm(LAPACK_COL_MAJOR, CblasNoTrans, CblasTrans
                , I_inte->sum_di, I_inte->sum_di, I_occ
                , 2
                , IO_z, I_inte->sum_di
                , IO_z, I_inte->sum_di
                , 0
                , p_i, I_inte->sum_di
        );

        rfock(I_inte->sum_di, I_inte->len, I_inte->lenh
                , I_inte->hcore, I_inte->eri, p_i
                , F
        );

        cblas_daxpy(I_inte->len, -1, p_i, 1, delta_p_i, 1);
        error = cblas_dnrm2(I_inte->len, delta_p_i, 1);
        printf("%d D_error:%15e\n", count, error);

        for(int i=0; i<I_inte->sum_di; i++){
            E += p_i[i*I_inte->sum_di+i] * (I_inte->hcore[i*(i+1)/2+i] + F[i * (i + 1) / 2 + i]);
        }
        for(int i=0; i<I_inte->sum_di; i++){
            for(int j=0; j<i; j++) {
                E += 2 * p_i[i*I_inte->sum_di+j] * (I_inte->hcore[i*(i+1)/2+j] + F[i * (i + 1) / 2 + j]);
            }
        }
        E /= 2;
        printf("  delta E:%15e\n  Energy:%15e\n", E-E_old, E);
    }
    if(count==I_loop){
        printf("scf not converged\n");
    }else{
        printf("converged in %d loop\n", count);
    }
    return E;
}
