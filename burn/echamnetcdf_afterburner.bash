#!/bin/bash
#

# Intel-Compiler:

#ifort -o echamnetcdf_afterburner -I/usr/local/netcdf-3.6.3-ifort/include echamnetcdf_afterburner.f90 #interpolate_index.f90 -L/usr/local/netcdf-3.6.3-ifort/lib -lnetcdf

# for beverin

ifort -o echamnetcdf_afterburner -I/opt/netcdf-3.6.3/include echamnetcdf_afterburner.f90 interpolate_index.f90 -L/opt/netcdf-3.6.3/lib -lnetcdf





# Lahey-Compiler:

#F90DEBUG="-Cpp -g --chk[s] --chk[e] --co  --f95 --lst --nsav --xref --warn --verbose --trap --trace"

#lf95 -o echamnetcdf_afterburner "$F90DEBUG" -I/usr/local/netcdf-lf95/include echamnetcdf_afterburner.f90 interpolate_index.f90 -L/usr/local/netcdf-lf95/lib -lnetcdf

