# little SCF

## Introduction

This is a rough program that can calculate molecule energy using RHF, UHF, MP2, and CIS. It uses Libcint to get the atomic integral and BLAS to solve some math problems. It will cost a lot of memory. Because the one-electron integral is stored in a packed symmetric matrix, but the two-electron integral is stored in a 4-dimensional array.

## Installation

There have two-way to get this terrible program.

1. Download the precompiled package that matches your platform.
2. If you want to compile it by yourself, you need to prepare some libraries, such as Pthread, Libcint and BLAS. Then type this command:
``` bash
export LIB_PATH_LIBCINT=<where your libcint install>
mkdir build
cd ./build
cmake ..
make -j
```

## How to use

You can use this command for a test:
``` bash
./little_SCF H2.xyz -rhf
```

It will output:

```
Raed xyz file (../test/H2.xyz)...OK

Raed bas file (/root/little_SCF/data/bas_sto-3g)...OK

Raed env file (/root/little_SCF/data/env_sto-3g)...OK

Calculate integral...
Nuclear repulsion:    5.291789e-01 Eh
OVERLAP
   1.000000e+00
   4.964866e-01   1.000000e+00
HCORE
  -9.794978e-01
  -6.828680e-01  -9.794978e-01

RHF running...
1 D_error:   1.099065e-15
  delta E:  -1.595288e+00
  Energy:  -1.595288e+00
2 D_error:   0.000000e+00
  delta E:   0.000000e+00
  Energy:  -1.595288e+00
converged in 2 loop
RHF DOWN!

Molecule orbital coefficient:
      -5.780276e-01      -5.780276e-01
       9.965050e-01      -9.965050e-01
Orbital energy:
    num    occ    E
      0      2      -4.844427e-01
      1      0       4.575041e-01

Molecule energy:  -1.066110e+00
```

And other methods are the same as this too, but UHF is an exception:

``` bash
./little_SCF O2.xyz -uhf 2

```

It means there are two more alpha electrons than beta electrons.

## Plan

The input module used in this program just can read the basis that defined elements one by one. So I'm rewriting the input module for higher performance and easier usage.