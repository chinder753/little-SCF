#include "byte2shell.h"

uint8_t byte2shell(const char I_fileName[], const Molecular * const I_atoms, Shells *O_shells) {
    // open the bas file
    FILE *bbyteFile = fopen(I_fileName, "rb");
    if(bbyteFile==NULL){
        printf("can't open bas file (%s)", I_fileName);
        exit(-1);
    }

    // read offset from file
    O_shells->num_shell = (uint8_t *)malloc(I_atoms->sum_eles * sizeof(uint8_t));
    uint16_t off_shell[I_atoms->sum_eles];
    uint16_t sum_shell = read_shell(bbyteFile, I_atoms, O_shells->num_shell, off_shell);

    uint16_t temp_count_shell = 0;
    O_shells->shell = (Shell *) malloc(sum_shell * sizeof(Shell));

    for(uint16_t i=0; i<(I_atoms->sum_eles); i++){
        fseek(bbyteFile, (int32_t)off_shell[i], SEEK_SET);
        for(uint16_t j=0; j < O_shells->num_shell[i]; j++){
            O_shells->shell[temp_count_shell].nuc = I_atoms->ele[i];
            // read shell
            O_shells->shell[temp_count_shell].ang = fgetc(bbyteFile);
            O_shells->shell[temp_count_shell].num_prim = fgetc(bbyteFile);
            O_shells->shell[temp_count_shell].num_contr = fgetc(bbyteFile);
            fread(&O_shells->shell[temp_count_shell].off_exp, EXP_OFF, 1, bbyteFile);
            fread(&O_shells->shell[temp_count_shell].off_coeff, COEFF_OFF, 1, bbyteFile);
            // read end
            temp_count_shell++;
        }
    }

    fclose(bbyteFile);

    O_shells->sum_shells = sum_shell;

    return 0;
}


uint16_t read_shell(FILE *I_bbyteFile, const Molecular * const I_atoms, uint8_t * const O_num_shell, uint16_t  * const O_off_shell){
    uint16_t sum_shell=0;
    for(uint8_t i=0; i < I_atoms->sum_eles; i++){
        fseek(I_bbyteFile, NUM_ATOMS + ATOM_BITS * I_atoms->ele[i], SEEK_SET);
        O_num_shell[i] = fgetc(I_bbyteFile);
        sum_shell += O_num_shell[i];
        O_off_shell[i]=0;
        fread(&O_off_shell[i], BAS_OFF, 1, I_bbyteFile);
    }
    return sum_shell;
}

void del_Shells(Shells *shells){
    free(shells->shell);
    free(shells->num_shell);
}