#ifndef STRUCT2LIBCINT_H
#define STRUCT2LIBCINT_H

#include "../tpdef.h"

uint8_t struct2libcint(const char fileName[], Libcint *libcint,const Molecular * const atoms, const Shells * const shells);
void del_Libcint(Libcint *libcint);

#endif //STRUCT2LIBCINT_H
