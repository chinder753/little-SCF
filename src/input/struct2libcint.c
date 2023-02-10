#include "struct2libcint.h"

#include <stdio.h>

#include "cint.h"

#include "byte2shell.h"

uint8_t struct2libcint(const char fileName[], Libcint *libcint, const Molecular * const atoms, const Shells * const shells){
    uint32_t len_bas = 0, len_env = PTR_ENV_START + 3 * atoms->sum_atoms;
    uint32_t temp_index=0;
    uint32_t exp_off[shells->sum_shells], coeff_off[shells->sum_shells];
    for(uint8_t i=1; i<atoms->sum_eles; i++){
        len_bas += atoms->num_atom[i] * shells->num_shell[i];
        for(uint8_t j=0; j<shells->num_shell[i]; j++){
            exp_off[temp_index] = len_env;
            coeff_off[temp_index] = len_env + shells->shell[temp_index].num_prim;
            len_env = coeff_off[temp_index] + shells->shell[temp_index].num_prim * shells->shell[temp_index].num_contr;
            temp_index++;
        }
    }

    libcint->atm = (Atm *) malloc(atoms->sum_atoms * sizeof(Atm));
    libcint->bas = (Bas *) malloc(len_bas * sizeof(Bas));
    libcint->env = (double *) malloc(len_env * sizeof(double));

    for(uint8_t i=0; i<PTR_ENV_START; i++){
        libcint->env[i] = 0.0;
    }

    uint8_t last_ele=0;  // index in atoms->ele
    temp_index = 0;
    for(uint16_t i=0; i<atoms->sum_atoms; i++){
        Atom now_atom = atoms->atom[i];
        libcint->atm[i].charge_of      = now_atom.nuc+1;
        libcint->atm[i].ptr_coord      = PTR_ENV_START + (int)i*3;
        libcint->atm[i].nuc_mod_of     = 0;
        libcint->atm[i].ptr_zeta       = 0;
        libcint->atm[i].ptr_frac_charge= 0;
        libcint->env[PTR_ENV_START + i*3+0] = now_atom.coordinates[0];
        libcint->env[PTR_ENV_START + i*3+1] = now_atom.coordinates[1];
        libcint->env[PTR_ENV_START + i*3+2] = now_atom.coordinates[2];

        if(now_atom.nuc == atoms->ele[last_ele]){
            for(uint16_t j=0; j<shells->num_shell[last_ele]; j++){
                libcint->bas[temp_index+j].atom_of  = (int)i;
                libcint->bas[temp_index+j].ang_of   = libcint->bas[temp_index - shells->num_shell[last_ele] + j].ang_of;
                libcint->bas[temp_index+j].nprim_of = libcint->bas[temp_index - shells->num_shell[last_ele] + j].nprim_of;
                libcint->bas[temp_index+j].nctr_of  = libcint->bas[temp_index - shells->num_shell[last_ele] + j].nctr_of;
                libcint->bas[temp_index+j].kappa_of = libcint->bas[temp_index - shells->num_shell[last_ele] + j].kappa_of;
                libcint->bas[temp_index+j].ptr_exp  = libcint->bas[temp_index - shells->num_shell[last_ele] + j].ptr_exp;
                libcint->bas[temp_index+j].ptr_coeff= libcint->bas[temp_index - shells->num_shell[last_ele] + j].ptr_coeff;
            }
        }else {
            uint16_t shell_index=0;
            for(uint16_t j = 1; j < atoms->sum_eles; j++){
                if(now_atom.nuc == atoms->ele[j]){
                    last_ele = j;
                    break;
                }
                shell_index += shells->num_shell[j];
            }
            for (uint16_t j = 0; j < shells->num_shell[last_ele]; j++) {
                libcint->bas[temp_index+j].atom_of   = (int)i;
                libcint->bas[temp_index+j].ang_of    = shells->shell[shell_index+j].ang;
                libcint->bas[temp_index+j].nprim_of  = shells->shell[shell_index+j].num_prim;
                libcint->bas[temp_index+j].nctr_of   = shells->shell[shell_index+j].num_contr;
                libcint->bas[temp_index+j].kappa_of  = 0;
                libcint->bas[temp_index+j].ptr_exp   = (int)exp_off[shell_index+j];
                libcint->bas[temp_index+j].ptr_coeff = (int)coeff_off[shell_index+j];
            }
        }
        temp_index += shells->num_shell[last_ele];
    }

    FILE *ebyteFile = fopen(fileName, "rb");
    temp_index = PTR_ENV_START + 3 * atoms->sum_atoms;
    for(uint32_t i=0; i<shells->sum_shells; i++){
        Shell now_shell = shells->shell[i];
        double temp_contr[shells->shell[i].num_prim];
        double temp_norm[shells->shell[i].num_prim];

        fseek(ebyteFile, (long int)now_shell.off_exp, SEEK_SET);
        fread(&libcint->env[temp_index], sizeof(double), now_shell.num_prim, ebyteFile);
        for(uint8_t j=0; j<now_shell.num_prim; j++){
            temp_norm[j] = CINTgto_norm(now_shell.ang, libcint->env[temp_index+j]);
        }
        temp_index += now_shell.num_prim;

        for(uint8_t j=0; j<now_shell.num_contr; j++) {
            fread(&temp_contr, sizeof(double), now_shell.num_prim, ebyteFile);
            for(uint8_t k=0; k<now_shell.num_prim; k++){
                libcint->env[temp_index+k] = temp_contr[k] * temp_norm[k];
            }
            temp_index += now_shell.num_prim;
        }
    }

    libcint->nbas = (int) len_bas;
    libcint->natm = (int) atoms->sum_atoms;
    libcint->nenv = (int) len_env;

    return 1;
}

void del_Libcint(Libcint *libcint){
    free(libcint->atm);
    free(libcint->bas);
    free(libcint->env);
}
