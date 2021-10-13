#!/bin/bash

# Here the parameters of the simulation are set.

# Martin Schraner, ETH Zuerich, February 2010

###############################################################################


CALL_SOCOL_DEFAULT=$(pwd)/call_socol_default_refc2.bash


# 1. Mandatory parameters:

EXPNO=00020
EXPERIMENTSYEAR=1950
EXPERIMENTEYEAR=2011
NPROCA=8                       # Number of divisions in north-south direction
NPROCB=5                        # Number of divisions in west-east direction
                                # NPROCA*NPROCB = number of CPUs
                                # nlat/NPROCA must be >= 4!

# 2. Optional parameters:

EXPERIMENTEMONTH=1               # (default: 12)
RESTARTYEAR=2000                # (default: $EXPERIMENTSYEAR)
RESTARTMONTH=                   # (default: 1)
RES=42                          # (default: 42)
LEV=                            # (default: 39)
LSOCOL=                         # (default: .TRUE.)
LCHEM=                          # (default: .TRUE.)
CO2FAC=                         # applied at RESTARTMONTH of RESTARTYEAR (d: 1)
BRUTUS_RUNTIME_MODEL=           # (default: 08:00)
BRUTUS_RUNTIME_POST=03:00       # (default: 01:00)
BRUTUS_JOB_NMONTH=              # (default: 6 if LCHEM=.TRUE., =36 else)
BRUTUS_RUNTIME_KILL_FACT=       # (default: 2.0)
NPROMA=                         # (default: 16)
CODENAME=                       # (default: echam5_socol)
EXECUTABLENAME=                 # (default: echam5_${EXPNO})
CTRLLAB_BCOND_LINKS=${EXPNO}            # (default: default)
CTRLLAB_ECHAM_NAMELIST=${EXPNO}         # (default: default)
#CTRLLAB_AFTERLEVELS_CHEM_M_NAMELIST=  # (default: default)
CTRLLAB_AFTER_ECHAM_H_NAMELIST= # (default: default)
CTRLLAB_AFTER_ECHAM_C_NAMELIST= # (default: default)
CTRLLAB_AFTER_CHEM1_H_NAMELIST= # (default: default)
CTRLLAB_AFTER_CHEM1_C_NAMELIST= # (default: default)
CTRLLAB_AFTER_CHEM2_H_NAMELIST= # (default: default)
CTRLLAB_AFTER_CHEM2_C_NAMELIST= # (default: default)
#CTRLLAB_AFTER_O3ORIG_D_NAMELIST= # (default: default)
CTRLLAB_RUNSCRIPT=              # (default: default)
CTRLLAB_POSTSCRIPT=             # (default: default)
LAFTER_ECHAM_H=                 # (default: .TRUE.)
LAFTER_ECHAM_C=                 # (default: .TRUE.)
LAFTER_CHEM1_H=                 # (default: .TRUE.)
LAFTER_CHEM1_C=                 # (default: .TRUE.)
LAFTER_CHEM2_H=                 # (default: .TRUE.)
LAFTER_CHEM2_C=                 # (default: .TRUE.)
#LAFTER_O3ORIG_D=                 # (default: .TRUE.)
CHEMINITIAL_FILE=T${RES}L39_20000restc_197412+tropchem+ods.nc               # (default: T${RES}L${LEV}_20000restc_197412.nc}
#CHEMINITIAL_FILE=T${RES}L39_20000restc_197412+tropchem.nc
RAW_FILETYPE_ECHAM=             # 1: GRIB, 2: NetCdf (default: 1)
RAW_FILETYPE_CHEM=              # 1: GRIB, 2: NetCdf (default: 1)
SMSNOTIFICATION=                # (default: .FALSE.)

# Call make make_run_default:
EXPNO=$EXPNO EXPERIMENTSYEAR=$EXPERIMENTSYEAR \
    EXPERIMENTEYEAR=$EXPERIMENTEYEAR EXPERIMENTEMONTH=$EXPERIMENTEMONTH \
    RESTARTYEAR=$RESTARTYEAR RESTARTMONTH=$RESTARTMONTH RES=$RES LEV=$LEV \
    LSOCOL=$LSOCOL LCHEM=$LCHEM CO2FAC=$CO2FAC NPROCA=$NPROCA NPROCB=$NPROCB \
    BRUTUS_RUNTIME_MODEL=$BRUTUS_RUNTIME_MODEL \
    BRUTUS_RUNTIME_POST=$BRUTUS_RUNTIME_POST \
    BRUTUS_JOB_NMONTH=$BRUTUS_JOB_NMONTH \
    BRUTUS_RUNTIME_KILL_FACT=$BRUTUS_RUNTIME_KILL_FACT NPROMA=$NPROMA \
    CODENAME=$CODENAME EXECUTABLENAME=$EXECUTABLENAME \
    CTRLLAB_BCOND_LINKS=$CTRLLAB_BCOND_LINKS \
    CTRLLAB_ECHAM_NAMELIST=$CTRLLAB_ECHAM_NAMELIST \
    CTRLLAB_AFTERLEVELS_CHEM_M_NAMELIST=$CTRLLAB_AFTERLEVELS_CHEM_M_NAMELIST \
    CTRLLAB_AFTER_ECHAM_H_NAMELIST=$CTRLLAB_AFTER_ECHAM_H_NAMELIST \
    CTRLLAB_AFTER_ECHAM_C_NAMELIST=$CTRLLAB_AFTER_ECHAM_C_NAMELIST \
    CTRLLAB_AFTER_CHEM1_H_NAMELIST=$CTRLLAB_AFTER_CHEM1_H_NAMELIST \
    CTRLLAB_AFTER_CHEM1_C_NAMELIST=$CTRLLAB_AFTER_CHEM1_C_NAMELIST \
    CTRLLAB_AFTER_CHEM2_H_NAMELIST=$CTRLLAB_AFTER_CHEM2_H_NAMELIST \
    CTRLLAB_AFTER_CHEM2_C_NAMELIST=$CTRLLAB_AFTER_CHEM2_C_NAMELIST \
    CTRLLAB_AFTER_O3ORIG_D_NAMELIST=$CTRLLAB_AFTER_O3ORIG_D_NAMELIST \
    CTRLLAB_RUNSCRIPT=$CTRLLAB_RUNSCRIPT \
    CTRLLAB_POSTSCRIPT=$CTRLLAB_POSTSCRIPT \
    LAFTER_ECHAM_H=$LAFTER_ECHAM_H LAFTER_ECHAM_C=$LAFTER_ECHAM_C \
    LAFTER_CHEM1_H=$LAFTER_CHEM1_H CHEMINITIAL_FILE=$CHEMINITIAL_FILE \
    LAFTER_CHEM1_C=$LAFTER_CHEM1_C CHEMINITIAL_FILE=$CHEMINITIAL_FILE \
    LAFTER_CHEM2_H=$LAFTER_CHEM2_H CHEMINITIAL_FILE=$CHEMINITIAL_FILE \
    LAFTER_CHEM2_C=$LAFTER_CHEM2_C CHEMINITIAL_FILE=$CHEMINITIAL_FILE \
    LAFTER_O3ORIG_D=$LAFTER_O3ORIG_D CHEMINITIAL_FILE=$CHEMINITIAL_FILE \
    RAW_FILETYPE_ECHAM=$RAW_FILETYPE_ECHAM \
    RAW_FILETYPE_CHEM=$RAW_FILETYPE_CHEM SMSNOTIFICATION=$SMSNOTIFICATION \
    BRUTUS_YEARLOOP=$BRUTUS_YEARLOOP $CALL_SOCOL_DEFAULT

exit
