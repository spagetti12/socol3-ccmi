#
# modules to load:
#
# # On brutus:
# module load intel/11.0.081 netcdf open_mpi/1.3.3 
# 
# # On IAC Linux workstation
# module load ifort/10.1.017 netcdf/3.6.3-ifort openmpi/1.3.3-ifort-10.1.017
#

export

SHELL = /bin/sh
 
ARCH  = LINUX

srcdir = .
top_srcdir = .
prefix = .
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
libdir = ${exec_prefix}/lib
includedir = ${prefix}/include
oldincludedir = /usr/include
infodir = ${prefix}/info
mandir = ${prefix}/man

sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var

program_transform_name = s,x,x,

MPIROOT        = $(MPI_ROOT)
MPI_LIB        = -L$(MPI_ROOT)/lib
MPI_INCLUDE    = -I$(MPI_ROOT)/include

NETCDFROOT     = $(NETCDF)
NETCDF_LIB     = -L$(NETCDFROOT)/lib -lnetcdff -lnetcdf 
NETCDF_INCLUDE = -I$(NETCDFROOT)/include 

LIB      = -L../lib -lsupport -llapack -lblas 
LIBS     = $(LIB) $(NETCDF_LIB) $(MPI_LIB)

MODOPT   = -I
MODULES  = ../modules

INCLUDE  = ../include
INCLUDES = $(MODOPT)$(MODULES) -I$(INCLUDE) $(NETCDF_INCLUDE) $(MPI_INCLUDE)

F90      = mpif90
FC       = mpif90
CC       = mpicc
CPP      = mpicc -E

DEFS     = -DHAVE_CONFIG_H

#CFLAGS   = -I../config -O3 -xT -fp-model precise -DNAGf90Fortran
CFLAGS   = -I../config -O3 -DNAGf90Fortran
#CFLAGS   = -I../config -O3
#FFLAGS   = -O3 -xT 
FFLAGS   = -O3
#FFLAGS   = $(INCLUDES) -axP -O3 -fp-model precise -r8
#F90FLAGS = $(INCLUDES) -O3 -xT -fp-model precise -r8 -fpp -DNOMPI 
#F90FLAGS = $(INCLUDES) -O3 -fp-model precise -r8 -fpp -DNOMPI 
#F90FLAGS = $(INCLUDES) -O3 -fp-model precise -r8 -fpp
#F90FLAGS = $(INCLUDES) -g -debug all -check all -inline-debug-info -traceback -fp-model precise -r8 -fpp
F90FLAGS = $(INCLUDES) -O3 -fp-model precise -r8 -fpp -traceback
CPPFLAGS =
ARFLAGS  = crv
#LDFLAGS  = -O3 -xT -fp-model precise -r8 -fpp   
#LDFLAGS  = -O3 -fp-model precise -r8 -fpp   
LDFLAGS  = -O3 -fp-model precise -r8 -fpp  

SRCDIRS = blas lapack support modules src

all:
	@for DIR in $(SRCDIRS) ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR ;\
	    $(MAKE) ; status=$$? ; \
	    if [ $$status != 0 ] ; then \
	      echo "Exit status from make was $$status" ; exit $$status ; \
	    fi ; \
	    cd $$back ; \
	  done 

clean:
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done
	-rm -f config.cache
	-rm -f lib/*.a bin/echam5
	-rm -f html/[a-z]*

tar:
	@tarfile=../echam5.f90.`date +%y%m%d`.taz ; gtar zcvf $$tarfile \
	`ls */*.f90 */*.[fhc] */*inc */Makefile Makefile.in Makefile run/hjob*`

index:
	-rm -f html/[a-z]*
	util/f2html.pl -f util/fgenrc -d html \
          support modules src include
