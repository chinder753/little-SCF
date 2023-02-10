#include "eint.h"

#include <math.h>
#include <stdio.h>

#include "diyblas.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define eri(i, j, l, m) out[(((i)*I_sum_di+(j))*I_sum_di+(l))*I_sum_di+(m)]

typedef CACHE_SIZE_T CINTIntegralFunction(double *out, FINT *dims, FINT *shls,
                                          FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
                                          CINTOpt *opt, double *cache);
typedef void CINTOptimizerFunction(CINTOpt **opt,
                                   FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);

extern CINTOptimizerFunction int2e_optimizer;
extern CINTIntegralFunction int2e_cart;
extern CINTOptimizerFunction int1e_ovlp_optimizer;
extern CINTIntegralFunction int1e_ovlp_cart;
extern CINTOptimizerFunction int1e_nuc_optimizer;
extern CINTIntegralFunction int1e_nuc_cart;
extern CINTOptimizerFunction int1e_kin_optimizer;
extern CINTIntegralFunction int1e_kin_cart;



int trint1e(double * out
        , Libcint *libcint
        , int I_max_di, int I_sum_di, const int * const I_di
        , CINTintf (*intf), CINToptf (*optf)
){

    CINTOpt *opt = NULL;
    optf( &opt, libcint->atm, libcint->natm, libcint->bas, libcint->nbas, libcint->env);

    int shls[2];
    double temp_buf[I_max_di * I_max_di];
    int temp_i=0, temp_j=0, start_j=0;

    for(int i=0; i<libcint->nbas; i++){
        shls[0] = i; shls[1] = i;
        intf(temp_buf
                , NULL, shls
                , libcint->atm, libcint->natm
                , libcint->bas, libcint->nbas
                , libcint->env
                , opt, NULL
        );
        for(int a=0; a < I_di[i]; a++){
            for(int b=a; b < I_di[i]; b++){
                out[(temp_i+b) * I_sum_di + a + temp_i] = temp_buf[a * I_di[i] + b];
            }
        }
        temp_i += I_di[i];
    }

    temp_i = I_di[0];

    for(int i=1; i<libcint->nbas; i++){
        shls[0] = i;
        start_j = 0;
        for(int j=0; j<i; j++){
            shls[1] = j;
            intf(temp_buf
                    , NULL, shls
                    , libcint->atm, libcint->natm
                    , libcint->bas, libcint->nbas
                    , libcint->env
                    , opt, NULL
            );
            int temp_index=0;
            for(int a=0; a < I_di[i]; a++){
                temp_j = start_j;
                for(int b=0; b < I_di[j]; b++){
                    out[(temp_i+a) * I_sum_di + (b + temp_j)] = temp_buf[temp_index];
                    temp_index++;
                }
            }
            start_j += I_di[j];
        }
        temp_i += I_di[i];
    }

    CINTdel_optimizer(&opt);

    return 0;
}



int spint1e(double * out
        , Libcint *libcint
        , const int I_sum_di, const int I_max_di, const int *I_di
        , CINTintf (*intf), CINToptf (*optf)
){

    double temp_int[I_sum_di][I_sum_di];
    trint1e(temp_int, libcint, I_max_di, I_sum_di, I_di, intf, optf);

    for(int i=0; i<I_sum_di; i++){
        for(int j=i; j<I_sum_di; j++){
            out[i+(j+1)*j/2] = temp_int[j][i];
        }
    }

    return 0;
}




int allint2(double * out
        , Libcint *libcint
        , int I_max_di, int I_sum_di, int *I_di
        , CINTintf (*intf), CINToptf (*optf)
){
    int len_bas = libcint->nbas*libcint->nbas;
    int bas_eri[len_bas][2], index_eri[len_bas][2], temp_index=0, temp_i=0, temp_j=0;
    for(int i=0; i<libcint->nbas; i++){
        temp_j = 0;
        for(int j=0; j<libcint->nbas; j++) {
            bas_eri[temp_index][0] = i;
            bas_eri[temp_index][1] = j;
            index_eri[temp_index][0] = temp_i;
            index_eri[temp_index][1] = temp_j;
            temp_j += I_di[j];
            temp_index++;
        }
        temp_i += I_di[i];
    }

    int shls[4];
    CINTOpt *opt=NULL;
    optf(&opt, libcint->atm, libcint->natm, libcint->bas, libcint->nbas, libcint->env);

    for(int i=0; i<len_bas; i++){
        shls[0] = bas_eri[i][0]; shls[1] = bas_eri[i][1];
        for(int j=0; j<len_bas; j++) {
            shls[2] = bas_eri[j][0]; shls[3] = bas_eri[j][1];
            double temp_buf[I_di[shls[3]] * I_di[shls[2]] * I_di[shls[1]] * I_di[shls[0]]];
            intf(temp_buf, NULL, shls, libcint->atm, libcint->natm, libcint->bas, libcint->nbas, libcint->env, opt, NULL);
            temp_index = 0;
            for(int d=0; d < I_di[shls[3]]; d++) {
                for(int c=0; c < I_di[shls[2]]; c++){
                    for(int b=0; b < I_di[shls[1]]; b++){
                        for(int a=0; a < I_di[shls[0]]; a++){
                            eri(index_eri[i][0]+a, index_eri[i][1]+b, index_eri[j][0]+c, index_eri[j][1]+d) = temp_buf[temp_index];
                            temp_index++;
                        }
                    }
                }
            }
        }
    }
    CINTdel_optimizer(&opt);

    return 0;
}



int get_eint(Molecular *I_atoms, Libcint *I_libcint, Integral *O_inte){
    O_inte->di = (int *) malloc(I_libcint->nbas*sizeof(int));
    O_inte->sum_di=0;
    O_inte->max_di=0;
    for(int i=0; i<I_libcint->nbas; i++){
        O_inte->di[i] = CINTcgto_cart(i, I_libcint->bas);
        O_inte->max_di = max(O_inte->max_di, O_inte->di[i]);
        O_inte->sum_di += O_inte->di[i];
    }
    O_inte->len = O_inte->sum_di*O_inte->sum_di;
    O_inte->lenh = (O_inte->sum_di+1)*O_inte->sum_di/2;

    O_inte->nucr = 0;
    for(int i=0; i < I_atoms->sum_atoms; i++){
        for(int j=0; j<i; j++){
            double R2=0;
            for(int k=0; k<3; k++) R2 += (I_atoms->atom[i].coordinates[k] - I_atoms->atom[j].coordinates[k]) * (I_atoms->atom[i].coordinates[k] - I_atoms->atom[j].coordinates[k]);
            O_inte->nucr += (I_atoms->atom[i].nuc + 1) * (I_atoms->atom[j].nuc + 1) / sqrt(R2);
        }
    }

    double nuc[O_inte->lenh];
    O_inte->overlap = (double *) malloc(O_inte->lenh * sizeof(double));
    O_inte->hcore = (double *) malloc(O_inte->lenh * sizeof(double));
    // overlap
    spint1e(O_inte->overlap, I_libcint, O_inte->sum_di, O_inte->max_di, O_inte->di, int1e_ovlp_cart,
            int1e_ovlp_optimizer);
    // nuc
    spint1e(nuc, I_libcint, O_inte->sum_di, O_inte->max_di, O_inte->di, int1e_nuc_cart, int1e_nuc_optimizer);
    // kin
    spint1e(O_inte->hcore, I_libcint, O_inte->sum_di, O_inte->max_di, O_inte->di, int1e_kin_cart, int1e_kin_optimizer);
    // hcore
    cblas_daxpy(O_inte->lenh
            , 1
            , nuc, 1
            , O_inte->hcore, 1
    );
    // eri
    O_inte->eri = (double *) malloc(O_inte->len*O_inte->len*sizeof(double));
    allint2(O_inte->eri, I_libcint, O_inte->max_di, O_inte->sum_di, O_inte->di, int2e_cart, int2e_optimizer);

    printf("Nuclear repulsion: %15e Eh\n", O_inte->nucr);
    spprintf(O_inte->sum_di, "OVERLAP", O_inte->overlap);
    spprintf(O_inte->sum_di, "HCORE", O_inte->hcore);
}