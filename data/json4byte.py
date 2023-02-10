import json, struct

atom_bits = 3
shell_bits = 11

# https://docs.python.org/zh-cn/3/library/struct.html#format-characters
num_atoms       = ">B"
num_shells      = ">B"
bas_off         = "<H"
ang             = ">B"
num_pri         = ">B"
num_contr       = ">B"
exp_off         = "<L"
coeff_off       = "<L"


"""
******************************
**********bas_sto-3g**********
******************************

uchar = unsigned char
usint = unsigned short int
ullint = unsigned long long int

total atoms
    [0,7]           uchar   8bits       m, the number of atoms defined in basic
one atom
    [8,15]          uchar   8bits       how many shells in atoms
    [16,47]         usint   16bits      offset to save bas
    ...
one electron shell
    [0,7]           uchar   8bits       angular momentum
    [8,15]          uchar   8bits       number of primitive GTO
    [16,23]         uchar   8bits       number of contracted GTO
    [24,87]         ullint  32bits      env offset to save exponents of primitive GTOs
    [56,151]        ullint  32bits      env offset to save column-major contraction coefficients
    ...

------------------
|   total atoms  |
------------------
|     atom_1     |  <-- 8
------------------
|     atom_2     |  <-- 8 + atom_bits
------------------
|       ...      |
------------------
| atom_1 shell_1 |  <-- offset_1
------------------
| atom_1 shell_2 |  <-- offset_1 + shell_bits
------------------
|       ...      |
------------------
| atom_2 shell_1 |  <-- offset_2
------------------
|       ...      |
------------------
"""


"""
******************************
**********env_sto-3g**********
******************************

---------------------------
|          ...            |
---------------------------
|   atom_1 exponents_1    |  <-- exp_offset_1
---------------------------
|   atom_1 exponents_2    |  <-- exp_offset_1 + 64
---------------------------
|          ...            |  <-- untill atom_1 exponents_<number of primitive GTO>
---------------------------
| atom_1 coefficients_1_1 |  <-- coeff_offset_1
----------------------
| atom_1 coefficients_1_2 |  <-- coeff_offset_1 + 64
---------------------------
|           ...           |  <-- untill atom_1 exponents_1_<number of contracted GTO>
---------------------------
| atom_1 coefficients_2_1 |  <-- coeff_offset_1 + <number of contracted GTO> * 64
---------------------------
|           ...           |  <-- untill atom_1 coefficients_<number of primitive GTO>_<number of contracted GTO>
---------------------------
|   atom_2 exponents_1    |  <-- exp_offset_2 = coeff_offset_1 
---------------------------
|           ...           |
---------------------------
"""


# read json from file
js_bas = ""
with open("sto-3g.1.json") as file:
    js_bas = json.load(file)['elements']


# initial offset, len(js_bas.keys()) is the number of atoms defined in basic
shls_offset = 1+len(js_bas.keys())*atom_bits
env_offset = 0
# initial byte info
atom_str = struct.pack(num_atoms, len(js_bas.keys()))
shell_str = struct.pack('B', 0)
env_str = struct.pack('B', 0)


for element in js_bas:
    sum_shells = 0
    temp_shls_offset = shls_offset
    for shls in js_bas[element]["electron_shells"]:
        sum_shells += len(shls["angular_momentum"])
        for index_ang in range(len(shls["angular_momentum"])):
            x_ang = shls["angular_momentum"][index_ang]
            shell_str += struct.pack(ang, x_ang)  # 8bit | angular momentum
            shell_str += struct.pack(num_pri, len(shls["exponents"]))  # 8bit | number of primitive GTO
            shell_str += struct.pack(num_contr, 1)  # 8bit | number of contracted GTO
            shell_str += struct.pack(exp_off, env_offset)  # 32bits | env offset to save exponents of primitive GTOs
            for exp in shls["exponents"]:
                env_str += struct.pack('<d', float(exp))  # 64bits | exponents
                env_offset += 8
            shell_str += struct.pack(coeff_off, env_offset)  # 32bits | env offset to save exponents of primitive GTOs
            for coeff in shls["coefficients"][index_ang]:
                env_str += struct.pack('<d', float(coeff))  # 64bits | coefficients
                env_offset += 8
            shls_offset += shell_bits
    atom_str += struct.pack(num_shells, sum_shells)  # 8bit | how many shells in atoms
    atom_str += struct.pack(bas_off, temp_shls_offset)  # 16bit | offset to save bas


with open("bas_sto-3g", "wb") as bas:
    bas.write(atom_str+shell_str[1:])
with open("env_sto-3g", "wb") as env:
    env.write(env_str[1:])