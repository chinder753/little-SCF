#ifndef TYPE_DEF_H
#define TYPE_DEF_H

#include <stdint.h>

typedef struct{
    int nuc;
    double coordinates[3];
} Atom;

typedef struct{
    uint8_t sum_eles;
    uint16_t sum_electron;
    uint16_t sum_atoms;
    uint8_t * ele;
    uint16_t * num_atom;
    Atom * atom;
} Molecular;

typedef struct{
    uint8_t nuc, ang, num_prim, num_contr;
    uint32_t off_exp, off_coeff;
}Shell;

typedef struct{
    uint16_t sum_shells;
    uint8_t * num_shell;
    Shell * shell;
}Shells;

typedef struct{
    int charge_of           ;
    int ptr_coord           ;
    int nuc_mod_of          ;
    int ptr_zeta            ;
    int ptr_frac_charge     ;
    int reserve_atmlot      ;
}Atm;

typedef struct{
    int atom_of       ;
    int ang_of        ;
    int nprim_of      ;
    int nctr_of       ;
    int kappa_of      ;
    int ptr_exp       ;
    int ptr_coeff     ;
    int reserve_baslot;
}Bas;

typedef struct{
    int natm;
    Atm *atm;
    int nbas;
    Bas *bas;
    int nenv;
    double *env;
}Libcint;

typedef struct{
    int sum_di;
    int max_di;
    int len;
    int lenh;
    int *di;
    double nucr;
    double *overlap;
    double *hcore;
    double *eri;
}Integral;

#endif //TYPE_DEF_H
