#ifndef BYTE2ENV_H
#define BYTE2ENV_H

#include <stdio.h>
#include <stdlib.h>
#include "../tpdef.h"


// this all defines same as libcint defined

#define ATOM_BITS   3
#define SHELL_BITS  11

#define NUM_ATOMS   1
#define NUM_SHELLS  1
#define BAS_OFF     2
#define ANG         1
#define NUM_PRI     1
#define NUM_CONTR   1
#define EXP_OFF     4
#define COEFF_OFF   4

// end define



/**
 * @brief byte2shell
 * @param I_fileName        <em>[I]</em>  bas file name
 * @param I_atoms           <em>[I]</em>  molecular read from xyz file
 * @param O_shells          <em>[O]</em>  shells read from bas file
 *
 * @return 0
 *      - read file success
*/
uint8_t byte2shell(const char I_fileName[], const Molecular * const I_atoms, Shells *O_shells);



/**
 *
 * @param I_bbyteFile         <em>[I]</em>  bas FILE pointer variable
 * @param I_atoms             <em>[I]</em>  molecular read from xyz file
 * @param O_num_shell         <em>[O]</em>  shells number of each atom read from bas file
 * @param O_off_shell         <em>[O]</em>  offset of shell read from bas file
 *
 * @return                    number of shell in atoms
 */
uint16_t read_shell(FILE *I_bbyteFile, const Molecular * const I_atoms, uint8_t * const O_num_shell, uint16_t  * const O_off_shell);



/**
 * @brief del_Shells
 * @param shells              <em>[I]</em>  the shells you want to delete
 */
void del_Shells(Shells *shells);

#endif //BYTE2ENV_H
