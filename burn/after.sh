#!/bin/sh
#
# Compiling of ECHAM5 afterburner:

gcc -g -O2 -o after after.c -lm -DHAVE_LIBNETCDF -I_NETCDF_ROOT_/include -L_NETCDF_ROOT_/lib -lnetcdf
