cmake_minimum_required(VERSION 3.16)

# compile setting
project(little_scf C)
set(CMAKE_C_STANDARD 99)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_FLAGS "-Wall -Wextra -O3")

# data file
add_definitions(-DPATH_BAS_FILE=\"${PROJECT_SOURCE_DIR}/data/bas_\")
add_definitions(-DPATH_ENV_FILE=\"${PROJECT_SOURCE_DIR}/data/env_\")


###############
###   lib   ###
###############
include_directories(src)
add_library(lscf
        STATIC
        src/input/xyz2atoms.c
        src/input/byte2shell.c
        src/input/struct2libcint.c
        src/calc/eint.c
        src/calc/scf.c
        src/calc/diyblas.c
        src/calc/fock.c
        src/calc/guess.c
        )
target_link_libraries(lscf PUBLIC m)
########## -lib END- ##########



###############
### 3rd lib ###
###############

# Blas
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
    message("BLAS has found in ${BLAS_LIBRARIES}")
    target_link_libraries(lscf PUBLIC ${BLAS_LIBRARIES})
endif()

# LAPACK
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    message("LAPACK has found in  ${LAPACK_LIBRARIES}")
    target_link_libraries(lscf PUBLIC ${LAPACK_LIBRARIES})
endif()

# LIBCINT
target_include_directories(lscf PUBLIC $ENV{LIB_PATH_LIBCINT}/include)
target_link_libraries(lscf PUBLIC $ENV{LIB_PATH_LIBCINT}/libcint.so)
########## -3rd lib END- ##########



##############
###  main  ###
##############
add_executable(little_scf
        src/main.c
        )

target_include_directories(little_scf PUBLIC src/calc)
target_include_directories(little_scf PUBLIC src/input)

find_package(Threads)
target_link_libraries(lscf PUBLIC ${CMAKE_THREAD_LIBS_INIT})

target_link_libraries(little_scf lscf)
########## -main END- ##########
