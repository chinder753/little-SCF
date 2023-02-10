#include "xyz2atoms.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char *sym_element[] = {
          "H", "He"
        , "Li", "Be", "B", "C", "N", "O", "F", "Ne"
        , "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"
        , "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"
        , "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pb", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"};




/**
 * @brief error_xyz_read
 * @param i               <em>[I]</em>  which lines arise error
 *
 * @warning               this funciton will kill the application
 */
void error_xyz_read(uint16_t i){
    printf("read xyz file failed in %d line\n", i);
    exit(-1);
}



/**
 * @brief isINele
 * @param I_target_ele    <em>[I]</em>  which element you need to find
 * @param I_ele           <em>[I]</em>  elements have been read from xyz file
 * @param I_num_ele       <em>[I]</em>  how many atom of one element
 * @param O_index_ele     <em>[O]</em>  if return is 1, it will be the index of target_ele in ele
 *
 * @return 1
 * - target_ele have found in ele
 * @return 0
 * - target_ele have not in ele
 */
_Bool isINele(const uint8_t I_target_ele, const uint8_t * const I_ele, const uint8_t I_num_ele, uint8_t * const O_index_ele){
    for(uint16_t i=1; i < I_num_ele; i++){
        if(I_target_ele == I_ele[i] ){
            *O_index_ele = i;
            return 1;
        }
    }
    return 0;
}



/**
 * @brief count_ele
 * @param I_atoms         <em>[I]</em>  atoms read from xyz file
 * @param I_sum_atoms     <em>[I]</em>  number of atom in atoms
 * @param O_ele           <em>[O]</em>  elements read from xyz file
 * @param O_num_atom      <em>[O]</em>  how many atom of one element
 *
 * @return              number of element in atoms
 * @example
 *   <table>
 *      <caption>C2H6O</caption>
 *      <tr>  <th>          <th>0    <th>1    <th>2    <th>3
 *      <tr>  <td>ele       <td>255  <td>6    <td>1    <td>8
 *      <tr>  <td>num_atom  <td>255  <td>2    <td>6    <td>1
 *   </table>
 */
uint8_t count_ele(const Atom *I_atoms, const uint8_t I_sum_atoms, uint8_t * const O_ele, uint16_t * const O_num_atom){
    uint8_t num_ele=1, index_ele=1, last_ele=255;
    O_ele[0] = 255;
    O_num_atom[0] = 255;
    for(int i=0; i < I_sum_atoms; i++){
        if(I_atoms[i].nuc == last_ele){
            O_num_atom[index_ele]++;
            continue;
        }
        if( isINele(I_atoms[i].nuc, O_ele, num_ele, &index_ele) ){
            O_num_atom[index_ele]++;
            last_ele = I_atoms[i].nuc;
            continue;
        }
        O_num_atom[num_ele] = 1;
        O_ele[num_ele] = I_atoms[i].nuc;
        index_ele = num_ele;
        num_ele++;
    }
    return num_ele;
}



/**
 * @brief parse_xyz
 * @param I_sum_atoms     <em>[I]</em>  number of atom in atoms
 * @param I_xyzFile       <em>[I]</em>  xyz FILE pointer variable
 * @param O_atoms         <em>[O]</em>  atoms read from xyz file
 *
 * @return 0
 *      - read file success
 */
uint16_t parse_xyz(uint16_t I_sum_atoms, FILE *I_xyzFile, Atom *const O_atoms) {
    for(uint16_t index_atom=0; index_atom < I_sum_atoms; index_atom++){
        char line_str[100];
        if(fgets(line_str, 100, I_xyzFile) == NULL) return index_atom + 3;
        char * sym = strtok(line_str, " ");
        if(sym==NULL)  return index_atom+3;
        for(uint8_t j=0; j<54; j++){
            if(strcmp(sym, sym_element[j])==0){
                O_atoms[index_atom].nuc = j;
                break;
            }
        }
        for(uint8_t i=0; i<3; i++){
            if(sym==NULL)  return index_atom;
            sym = strtok(NULL, " ");
            O_atoms[index_atom].coordinates[i] = atof(sym) * 1.88972;
        }
    }
    return 0;
}



/**
 * @brief xyz2atom
 * @param I_fileName      <em>[I]</em>  xyz file name
 * @param O_atoms         <em>[O]</em>  molecular read from xyz file
 *
 * @return 0
 *     - read file success
*/
uint8_t xyz2atom(const char I_fileName[], Molecular *O_atoms) {
    // open the xyz file
    FILE * xyzFile = fopen(I_fileName, "r");
    if(xyzFile==NULL){
        printf("can't open xyz file\n");
        exit(-1);
    }

    // read the number of atom
    char num_atoms_str[100];
    if( fgets(num_atoms_str, 100, xyzFile) == NULL ) error_xyz_read(1);

    // parse str to int
    uint16_t num_atoms = atoi(num_atoms_str);
    if( num_atoms == 0 ) error_xyz_read(1);

    // jump to atom line
    if( fgets(num_atoms_str, 100, xyzFile) == NULL ) error_xyz_read(1);

    // create atom array
    O_atoms->atom = (Atom *)malloc(num_atoms * sizeof(Atom));
    if(O_atoms->atom == NULL ){
        printf("malloc error");
        exit(-1);
    }

    // read and parse atom symbol and coordinates
    uint16_t index_atom = parse_xyz(num_atoms, xyzFile, O_atoms->atom);
    if(index_atom!=0) error_xyz_read(index_atom+3);

    fclose(xyzFile);

    // count how many element in atom
    uint8_t t_ele[50];
    uint16_t t_atom_num[50];
    uint8_t num_ele = count_ele(O_atoms->atom, num_atoms, t_ele, t_atom_num);
    O_atoms->ele = (uint8_t *)malloc(num_ele * sizeof(uint8_t));
    O_atoms->num_atom = (uint16_t *)malloc(num_ele * sizeof(uint16_t));
    for(uint8_t i=0; i<num_ele; i++){
        O_atoms->ele[i] = t_ele[i];
        O_atoms->num_atom[i] = t_atom_num[i];
    }

    O_atoms->sum_atoms = num_atoms;
    O_atoms->sum_eles = num_ele;

    O_atoms->sum_electron = 0;
    for(int i=1; i<num_ele; i++){
        O_atoms->sum_electron += (O_atoms->ele[i]+1)*O_atoms->num_atom[i];
    }

    return 0;
}



/**
 * @brief del_Molecular
 * @param I_atoms         <em>[I]</em>  the atoms you want to delete
 */
void del_Molecular(Molecular * I_atoms){
    free(I_atoms->atom);
    free(I_atoms->ele);
    free(I_atoms->num_atom);
}