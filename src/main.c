#include <string.h>
#include "little_SCF.h"

int main(int argc,char *argv []){
    char method[strlen(argv[2])];
    strcpy(method, argv[2]);
    for(int i=0; i<strlen(method); i++){
        if(method[i]>96 && method[i]<123){
            method[i] -= 32;
        }
    }

    if(argc!=3 && strcmp(method, "-UHF")!=0){
        printf("The number of argument must be 2.\n");
        return -1;
    }else if(argc!=4 && strcmp(method, "-UHF")==0){
        printf("The number of argument must be 3.\n");
        return -1;
    }

    Molecular atoms;
    printf("Raed xyz file (%s)...", argv[1]);
    xyz2atom(argv[1], &atoms);
    printf("OK\n\n");

    char bas_path[100] = PATH_BAS_FILE;
    strncat(bas_path, "sto-3g", 100- strlen(PATH_BAS_FILE));
    Shells shells;
    printf("Raed bas file (%s)...", bas_path);
    byte2shell(bas_path, &atoms, &shells);
    printf("OK\n\n");

    char env_path[100] = PATH_ENV_FILE;
    strncat(env_path, "sto-3g", 100-strlen(PATH_ENV_FILE));
    Libcint libcint;
    printf("Raed env file (%s)...", env_path);
    struct2libcint(env_path, &libcint, &atoms, &shells);
    printf("OK\n\n");

    Integral integral;
    printf("Calculate integral...\n");
    get_eint(&atoms, &libcint, &integral);
    printf("\n");


    if(strcmp(method, "-RHF")==0){
        printf("RHF running...\n");
        double z[integral.len], w[integral.sum_di];
        init_hcore(&integral, z);
        double e = rhf(200, 1e-6, atoms.sum_electron/2, &integral, z, w);
        printf("RHF DOWN!\n\n");

        printf("Molecule orbital coefficient:\n");
        for(int i=0; i<integral.sum_di; i++){
            for(int j=0; j<integral.sum_di; j++){
                printf("    %15e", z[i*integral.sum_di+j]);
            }
            printf("\n");
        }

        printf("Orbital energy:\n    num    occ    E\n");
        for(int i=0; i<atoms.sum_electron/2; i++){
            printf("    %3d    %3d    %15e\n", i, 2, w[i]);
        }
        for(int i=atoms.sum_electron/2; i<integral.sum_di; i++){
            printf("    %3d    %3d    %15e\n", i, 0, w[i]);
        }
        printf("\nMolecule energy:%15e\n", e+integral.nucr);
    }else if(strcmp(method, "-UHF")==0){
        int spin= atoi(argv[3]);
        printf("UHF running...\n\nMolecule spin is %d.\n", spin);

        double z_a[integral.len], z_b[integral.len];
        int alpha[] = {atoms.sum_electron/2+spin, 0}, beta[] = {atoms.sum_electron/2-spin, 0};
        init_hcore(&integral, z_a);
        if(spin==0){
            init_uhf(integral.sum_di, integral.len, atoms.sum_electron/2, z_a, z_b);
        }else{
            memmove(z_b, z_a, integral.len*sizeof(double));
        }
        uhf(200, 1, 1e-6, integral.sum_di, integral.len, integral.lenh
            , alpha, beta
            , integral.overlap, integral.hcore, integral.eri
            , z_a, z_b);
        printf("UHF DOWN!\n\n");

        printf("Molecule alpha orbital coefficient:\n");
        for(int i=0; i<integral.sum_di; i++){
            for(int j=0; j<integral.sum_di; j++){
                printf("    %15e", z_a[i*integral.sum_di+j]);
            }
            printf("\n");
        }
        printf("\n");
        printf("Molecule beta orbital coefficient:\n");
        for(int i=0; i<integral.sum_di; i++){
            for(int j=0; j<integral.sum_di; j++){
                printf("    %15e", z_b[i*integral.sum_di+j]);
            }
            printf("\n");
        }
    }else if(strcmp(method, "-MP2")==0){
        printf("RHF running...\n");
        double z[integral.len], w[integral.sum_di];
        init_hcore(&integral, z);
        double e = rhf(200, 1e-6, atoms.sum_electron/2, &integral, z, w);
        printf("RHF DOWN!\n\n");

        printf("MP2 running...");
        double e_mp2 = rmp2(atoms.sum_electron/2, &integral, z, w);
        printf("DOWN!\n\n");

        printf("rhf:%15e\nmp2:%15e\ntotal:%15e\n", e, e_mp2, e+e_mp2);
    }else if(strcmp(method, "-CIS")==0){
        printf("RHF running...\n");
        double z[integral.len], w[integral.sum_di];
        init_hcore(&integral, z);
        double e = rhf(200, 1e-6, atoms.sum_electron/2, &integral, z, w);
        printf("RHF DOWN!\n\n");

        double excited[atoms.sum_electron*(integral.sum_di*2 - atoms.sum_electron)];
        printf("CIS running...\n");
        rcis(atoms.sum_electron, &integral, z, w, excited);
        printf("CIS DOWN!\n\n");

        printf("Ground state energy:%15e\n", e);
        printf("Excited energy:\n");
        for(int i=0; i<atoms.sum_electron*(integral.sum_di*2 - atoms.sum_electron);i++){
            printf("    %3d:%15e\n", i, excited[i]);
        }
    }else{
        printf("Don't know method %s\n", argv[2]);
    }
}