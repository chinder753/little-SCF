#ifndef FOCK_H
#define FOCK_H

int rfock(const int I_sum_di, const int I_len, const int I_lenh
          , const double * const I_hcore, const double * const I_eri, const double * const I_p
          , double * const O_fock
          );

int ofock_c(const int I_sum_di, const int I_len, const double I_gamma
            , const double * const I_hcore, const double * const I_jk, const double * const I_j, const double * const I_k
            , double * const O_fock
            );

int ofock_o(const int I_sum_di, const int I_len, const double I_agamma, const double I_bgamma
            , const double * const I_hcore, const double * const I_jk, const double * const I_j, const double * const I_k
            , double * const O_fock
            );

int ufock_build(const int I_sum_di, const int I_len, const int I_lenh
                , const int * const I_alpha, const int * const I_beta
                , const double * const I_hcore, const double * const I_eri, const double * const I_p_a, const double * const I_p_b
                , double * const O_fock_a, double * const O_fock_b
                );

int u_density(const int I_sum_di
              , const int * const I_alpha, const int * const I_beta
              , const double * const I_z_alpha, const double * const I_z_beta
              , double * const O_p_alpha, double * const O_p_beta
              );

double cis_Hamilton(const int I_sum_di
                    , const int I_i, const int I_j
                    , const int I_a, const int I_b
                    , const double const * I_F, const double const * I_eri
                    );

#endif
