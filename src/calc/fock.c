#include "fock.h"

#include "diyblas.h"

#define eri(i, j, l, m) I_eri[(((i)*I_sum_di+(j))*I_sum_di+(l))*I_sum_di+(m)]


//////////////////////////
// restricted hartree fock
//////////////////////////
int rfock(const int I_sum_di, const int I_len, const int I_lenh
        , const double * const I_hcore, const double * const I_eri, const double * const I_p
        , double * const O_fock
        ){
    double JK[I_sum_di][I_sum_di];
    cblas_dcopy(I_lenh, I_hcore, 1, O_fock, 1);
    for (int i = 0; i < I_sum_di; i++) {
        for (int j = 0; j < i+1; j++) {
            cblas_dcopy(I_len, &eri(i, j, 0, 0), 1, JK, 1);
            for (int l = 0; l < I_sum_di; l++) {
                cblas_daxpy(I_sum_di, -0.5, &eri(j, l, i, 0), 1, JK[l], 1);
            }
            O_fock[i * (i + 1) / 2 + j] += cblas_ddot(I_len, I_p, 1, JK, 1);
        }
    }
}


//////////////////////////
// open shell hartree fock
//////////////////////////
int ofock_c(const int I_sum_di, const int I_len, const double I_gamma
          , const double * const I_hcore, const double * const I_jk, const double * const I_j, const double * const I_k
          , double * const O_fock
          ){
    double temp_jk[I_len];
    // temp_jk = I_j + I_k
    cblas_daxpy(I_len, 1, I_j, 1, I_k, 1);
    // O_fock = I_hcore
    cblas_dcopy(I_len, I_hcore, 1, O_fock, 1);
    // O_fock = I_hcore + I_jk
    cblas_daxpy(I_len, 1, I_jk, 1, O_fock, 1);
    // O_fock = I_hcore + I_jk + I_gamma * ( I_j + I_k )
    cblas_daxpy(I_len, I_gamma, temp_jk, 1, O_fock, 1);
}



int ofock_o(const int I_sum_di, const int I_len, const double I_agamma, const double I_bgamma
            , const double * const I_hcore, const double * const I_jk, const double * const I_j, const double * const I_k
            , double * const O_fock
            ){
    // O_fock = I_hcore
    cblas_dcopy(I_len, I_hcore, 1, O_fock, 1);
    // O_fock = I_hcore + I_jk
    cblas_daxpy(I_len, 1, I_jk, 1, O_fock, 1);
    // O_fock = I_hcore + I_jk + I_agamma * I_j
    cblas_daxpy(I_len, I_agamma, I_j, 1, O_fock, 1);
    // O_fock = I_hcore + I_jk + I_agamma * I_j + I_bgamma * I_k
    cblas_daxpy(I_len, I_bgamma, I_k, 1, O_fock, 1);
}



// SPF + FPS
int spf(const int I_sum_di, const int I_len
        , const double * const I_s, const double * const I_p, const double * const I_fa
        , double const * O_spf
        ){
    double temp_s[I_len],temp_fa[I_len], temp_spf[I_len];
    sp2sy(I_sum_di, I_s, temp_s);
    sp2sy(I_sum_di, I_fa, temp_fa);
    // temp_s \cdot I_p
    cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper
                , I_sum_di, I_sum_di
                , 1
                , temp_s, 1
                , I_p, 1
                , 0
                , temp_spf, I_sum_di
                );
    // temp_s \cdot I_p \cdot I_fa
    cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper
                , I_sum_di, I_sum_di
                , 1
                , temp_spf, 1
                , temp_fa, 1
                , 0
                , temp_spf, I_sum_di
                );
    tran(I_sum_di, temp_spf, temp_s);
    cblas_daxpy(I_sum_di, 1, temp_spf, 1, temp_s, 1);
    sy2sp(I_sum_di, temp_s, O_spf);
}



int h(const int I_sum_di, const int I_len, const int I_lenh, const double I_gamma
      , const double * const I_s, const double * const I_p_o, const double * const I_p_c, const double * const I_fo, const double * const I_fc
      , double * const O_h_o, double * const O_h_c
      ){
    double temp_gamma = I_gamma/(1-I_gamma);
    double tr_fa[I_lenh];
    // tr_fa = I_fo
    cblas_dcopy(I_lenh, I_fo, 1, tr_fa, 1);
    // tr_fa = I_fo - I_fc
    cblas_daxpy(I_lenh, -1, I_fc, 1, tr_fa, 1);
    // tr_fa = temp_gamma * ( I_fo - I_fc )
    cblas_dscal(I_lenh, temp_gamma, tr_fa, 1);
    double spf_o[I_lenh], spf_c[I_lenh];
    spf(I_sum_di, I_len
        , I_s, I_p_o, tr_fa
        , spf_o
        );
    spf(I_sum_di, I_len
        , I_s, I_p_c, tr_fa
        , spf_c
        );
    cblas_daxpy(I_lenh, -1, I_fc, 1, spf_o, 1);
    cblas_daxpy(I_lenh, -1/I_gamma, I_fo, 1, spf_c, 1);
}



//////////////////////////
// unrestricted hartree fock
//////////////////////////
int ufock(const int I_sum_di, const int I_len
          , const int I_m, const int I_n
          , const int * const I_electron
          , const double * const I_eri, const double * const I_p, const double * const I_p_t
          , double * const O_fock
          ){
    double J[I_sum_di][I_sum_di], K[I_sum_di][I_sum_di];
    cblas_dcopy(I_len, &eri(I_m, I_n, 0, 0), 1, J, 1);
    for (int l = 0; l < I_sum_di; l++) {
        cblas_dcopy(I_sum_di, &eri(I_n, l, I_m, 0), 1, K[l], 1);
    }
    for(int i=2; i<I_electron[1]; i+=2){
        cblas_dswap(I_sum_di, &J[I_electron[i]][0], 1, &J[I_electron[i+1]][0], 1);
        cblas_dswap(I_sum_di, &K[I_electron[i]][0], 1, &K[I_electron[i+1]][0], 1);
        cblas_dswap(I_sum_di, &J[0][I_electron[i]], I_sum_di, &J[0][I_electron[i+1]], I_sum_di);
        cblas_dswap(I_sum_di, &K[0][I_electron[i]], I_sum_di, &K[0][I_electron[i+1]], I_sum_di);
    }
    *O_fock += cblas_ddot(I_len, I_p_t, 1, J, 1) - cblas_ddot(I_len, I_p, 1, K, 1);
}



int u_density(const int I_sum_di
              , const int * const I_alpha, const int * const I_beta
              , const double * const I_z_alpha, const double * const I_z_beta
              , double * const O_p_alpha, double * const O_p_beta
              ){
    double alpha[I_sum_di*I_alpha[0]], beta[I_sum_di*I_beta[0]];
    cblas_dcopy(I_sum_di*I_alpha[0], I_z_alpha, 1, alpha, 1);
    cblas_dcopy(I_sum_di*I_beta[0], I_z_beta, 1, beta, 1);
    for(int i=2; i<I_alpha[1]; i+=2){
        cblas_dcopy(I_sum_di, &I_z_alpha[I_alpha[2+i]], 1, &alpha[I_alpha[3+i]], 1);
    }
    for(int i=2; i<I_beta[1]; i+=2){
        cblas_dcopy(I_sum_di, &I_z_beta[I_beta[2+i]], 1, &beta[I_beta[3+i]], 1);
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans
            , I_sum_di, I_sum_di, I_alpha[0]
            , 1
            , alpha, I_sum_di
            , alpha, I_sum_di
            , 0
            , O_p_alpha, I_sum_di
    );
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans
            , I_sum_di, I_sum_di, I_beta[0]
            , 1
            , beta, I_sum_di
            , beta, I_sum_di
            , 0
            , O_p_beta, I_sum_di
    );
}



int ufock_build(const int I_sum_di, const int I_len, const int I_lenh
                , const int * const I_alpha, const int * const I_beta
                , const double * const I_hcore, const double * const I_eri, const double * const I_p_a, const double * const I_p_b
                , double * const O_fock_a, double * const O_fock_b
                ){
    double p_t[I_len];
    // P_t = P_alpha
    cblas_dcopy(I_len, I_p_a, 1, p_t, 1);
    // P_t = P_alpha + P_beta
    cblas_daxpy(I_len, 1, I_p_b, 1, p_t, 1);
    // F_alpha = H_core
    cblas_dcopy(I_lenh, I_hcore, 1, O_fock_a, 1);
    // F_beta = H_core
    cblas_dcopy(I_lenh, I_hcore, 1, O_fock_b, 1);
    for (int i = 0; i < I_sum_di; i++) {
        for (int j = 0; j < i+1; j++) {
            ufock(I_sum_di, I_len
                  , i, j
                  , I_alpha
                  , I_eri, I_p_a, p_t
                  , &O_fock_a[i*(i+1)/2+j]
                  );
            ufock(I_sum_di, I_len
                  , i, j
                  , I_beta
                  , I_eri, I_p_b, p_t
                  , &O_fock_b[i*(i+1)/2+j]
                  );
        }
    }
}



//////////////////////////
// CI Hamilton
//////////////////////////
