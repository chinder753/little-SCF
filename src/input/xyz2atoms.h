/**
 * @file xyz2atoms.h
 * @details define some function for read molecular information from xyz file
 * @author chinder753
 * @daate 2022.12.22
 * @version v1.0
 */

#ifndef XYZ2ATOMS_H
#define XYZ2ATOMS_H

#include "../tpdef.h"

uint8_t xyz2atom(const char I_fileName[], Molecular *O_atoms);

void del_Molecular(Molecular * I_atoms);

#endif //XYZ2ATOMS_H
