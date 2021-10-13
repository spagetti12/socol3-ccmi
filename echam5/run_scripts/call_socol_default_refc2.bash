#!/bin/bash

# Prepares SOCOL simulation on cluster and calls *run_socol_default.bash*.
# *call_socol_default.bash*  is called from *start_run_template.bash*. 

# (to be called as EXPNO=value1 EXPERIMENTSYEAR=value2 ... <dirpath>/call_socol_default.bash) 

# Keyword parameters:

# EXPNO : Experiment number
# EXPERIMENTSYEAR
# EXPERIMENTEYEAR
# EXPERIMENTEMONTH (optional; default: 12)
# RESTARTYEAR (optional; default: $EXPERIMENTSYEAR)
# RESTARTMONTH (optional; default: 1)
# RES : Horizontal resolution (RES) (optional; default: 42)
# LEV : Vertical resoultion (LEV) (optional; default: 39)
# LSOCOL : True for SOCOL (optional; default: .TRUE.)
# LCHEM : True for chemistry module (optional; default: .TRUE.)
# CO2FAC : CO2 factor for creating ensembles; applied at RESTARTMONTH of
#            RESTARTYEAR (optional; default: 1)
# NPROCA
# NPROCB
# BRUTUS_RUNTIME_MODEL : runtime of one BRUTUS job for model given as 
#            HH:MM (optional; default: 08:00)
# BRUTUS_RUNTIME_POST : runtime for one BRUTUS job for postprocessing
#            given as HH:MM (optional; default: 01:00)
# BRUTUS_JOB_NMONTH : Number of months per submitted job
#            (optional; default: 6)
# BRUTUS_RUNTIME_KILL_FACT : If remaining model runtime of current BRUTUS job
#            is less than $BRUTUS_RUNTIME_KILL_FACT*(mean runtime of a month),
#            all BRUTUS jobs are killed and re-started (optional; default: 2.0)
# NPROMA : (optional; default: 16)
# CODENAME : (optional; default: echam5_socol)
# EXECUTABLENAME : (optional; default: echam5_${EXPNO})
# CTRLLAB_BCOND_LINKS : Control number for links of boundary conditions
#            (optional; default: default)
# CTRLLAB_ECHAM_NAMELIST : Control number for ECHAM namelist 
#            (optional; default: default)
# CTRLLAB_AFTERLEVELS_CHEM_M_NAMELIST : Control number for namelist of 
#            afterburner levels of chem_m files (optional; default: default)
# CTRLLAB_AFTER_ECHAM_H_NAMELIST : Control number for namelist of echam_h files
#            (optional; default: default)
# CTRLLAB_AFTER_ECHAM_C_NAMELIST : Control number for namelist of echam_c files
#            (optional; default: default)
# CTRLLAB_AFTER_CHEM1_H_NAMELIST : Control number for namelist of chem1_h files
#            (optional; default: default)
# CTRLLAB_AFTER_CHEM1_C_NAMELIST : Control number for namelist of chem1_c files
#            (optional; default: default)
# CTRLLAB_AFTER_CHEM2_H_NAMELIST : Control number for namelist of chem2_h files
#            (optional; default: default)
# CTRLLAB_AFTER_CHEM2_C_NAMELIST : Control number for namelist of chem2_c files
#            (optional; default: default)
# CTRLLAB_AFTER_O3ORIG_NAMELIST : Control number for namelist of o3orig_d files
#            (optional; default: default)
# CTRLLAB_RUNSCRIPT : Control number for SOCOL run script
#            (optional; default: default)
# CTRLLAB_POSTSCRIPT : Control number for afterburner
#            (optional; default: default)
# LAFTER_ECHAM_H : True for echam_h (optional; default: .TRUE.)
# LAFTER_ECHAM_C : True for echam_c (optional; default: .TRUE.)
# LAFTER_CHEM1_H : True for chem1_h (optional; default: .TRUE.)
# LAFTER_CHEM1_C : True for chem1_c (optional; default: .TRUE.)
# LAFTER_CHEM2_H : True for chem2_h (optional; default: .TRUE.)
# LAFTER_CHEM2_C : True for chem2_c (optional; default: .TRUE.)
# LAFTER_O3ORIG_D : True for o3orig_d afterburner (optional; default: .TRUE.)
# CHEMINITIAL_FILE : Initial file for chemistry 
#            (optional; T${RES}L${LEV}_20000restc_197412.nc}
# RAW_FILETYPE_ECHAM : Output file type for ECHAM-files. 1: GRIB, 2: NetCdf
#            (optional; default: 1)
# RAW_FILETYPE_CHEM : Output file type for chem-files. 1: GRIB, 2: NetCdf
#            (optional; default: 1)
# SMSNOTIFICATION : True for SMS notification in case of killed job
#            (optional; default: .FALSE.)


# Martin Schraner, ETH Zuerich, February 2010

###############################################################################


# Functions:
kill_job(){
    [ -f ${ERR_MES} ] \
	&& echo "job.exp.${EXPNO} killed at Year $1, Month $2" >> ${ERR_MES} \
	|| echo "job.exp.${EXPNO} killed at Year $1, Month $2" > ${ERR_MES}
    #Email-Notification:
    mail -s "job.exp.${EXPNO} killed" $MAIL_NOTIFICATION < $ERR_MES
    #SMS-Notification:
    [ $SMSNOTIFICATION = .TRUE. ] \
	&& mail -s "job.exp.${EXPNO} killed" $SMS_NOTIFICATION < $ERR_MES
    exit
}

export -f kill_job
#################################################


#################################################
# Check input:
#################################################

: ${EXPNO:?'EXPNO is missing'}; export EXPNO
: ${EXPERIMENTSYEAR:?'EXPERIMENTSYEAR is missing'}; export EXPERIMENTSYEAR
: ${EXPERIMENTEYEAR:?'EXPERIMENTEYEAR is missing'}; export EXPERIMENTEYEAR
EXPERIMENTEMONTH=${EXPERIMENTEMONTH:-12}
EXPERIMENTEMONTH=${EXPERIMENTEMONTH#0}  # remove leading zero
export RESTARTYEAR=${RESTARTYEAR:-$EXPERIMENTSYEAR}
export RESTARTMONTH=${RESTARTMONTH:-1}
RESTARTMONTH=${RESTARTMONTH#0}  # remove leading zero
export RES=${RES:-42}
export LEV=${LEV:-39}
export LSOCOL=${LSOCOL:-.TRUE.}
export LCHEM=${LCHEM:-.TRUE.}
export CHEMINITIAL_FILE=${CHEMINITIAL_FILE:-T${RES}L${LEV}_20000restc_197412.nc}
export CO2FAC=${CO2FAC:-1.00000}
: ${NPROCA:?'NPROCA is missing'}; export NPROCA
: ${NPROCB:?'NPROCB is missing'}; export NPROCB
export NCPUS=$((NPROCA*NPROCB))
export BRUTUS_RUNTIME_MODEL=${BRUTUS_RUNTIME_MODEL:-08:00}
export BRUTUS_RUNTIME_POST=${BRUTUS_RUNTIME_POST:-01:00}
BRUTUS_JOB_NMONTH=${BRUTUS_JOB_NMONTH:-default}
if [ ${BRUTUS_JOB_NMONTH} = default ]; then
    if ([ $LSOCOL = .TRUE. ] && [ $LCHEM = .TRUE. ]); then
	[ $RES -eq 31 ] && BRUTUS_JOB_NMONTH=8 || BRUTUS_JOB_NMONTH=6
    else
	BRUTUS_JOB_NMONTH=36
    fi	
fi
export BRUTUS_RUNTIME_KILL_FACT=${BRUTUS_RUNTIME_KILL_FAC:-2.0}
export NPROMA=${NPROMA:-16}
CODENAME=${CODENAME:-echam5_socol}
export EXECUTABLENAME=${EXECUTABLENAME:-echam5_${EXPNO}}
export CTRLLAB_BCOND_LINKS=${CTRLLAB_BCOND_LINKS:-default}
CTRLLAB_AFTER_ECHAM_H_NAMELIST=${CTRLLAB_AFTER_ECHAM_H_NAMELIST:-default}
CTRLLAB_AFTER_ECHAM_C_NAMELIST=${CTRLLAB_AFTER_ECHAM_C_NAMELIST:-default}
CTRLLAB_AFTER_CHEM1_H_NAMELIST=${CTRLLAB_AFTER_CHEM1_H_NAMELIST:-default}
CTRLLAB_AFTER_CHEM1_C_NAMELIST=${CTRLLAB_AFTER_CHEM1_C_NAMELIST:-default}
CTRLLAB_AFTER_CHEM2_H_NAMELIST=${CTRLLAB_AFTER_CHEM2_H_NAMELIST:-default}
CTRLLAB_AFTER_CHEM2_C_NAMELIST=${CTRLLAB_AFTER_CHEM2_C_NAMELIST:-default}
#CTRLLAB_AFTER_O3ORIG_D_NAMELIST=${CTRLLAB_AFTER_O3ORIG_D_NAMELIST:-default}
CTRLLAB_RUNSCRIPT=${CTRLLAB_RUNSCRIPT:-default}
export CTRLLAB_POSTSCRIPT=${CTRLLAB_POSTSCRIPT:-default}
export LAFTER_ECHAM_H=${LAFTER_ECHAM_H:-.TRUE.}
export LAFTER_ECHAM_C=${LAFTER_ECHAM_C:-.TRUE.}
export LAFTER_CHEM1_H=${LAFTER_CHEM1_H:-.TRUE.}
export LAFTER_CHEM1_C=${LAFTER_CHEM1_C:-.TRUE.}
export LAFTER_CHEM2_H=${LAFTER_CHEM2_H:-.TRUE.}
export LAFTER_CHEM2_C=${LAFTER_CHEM2_C:-.TRUE.}
#export LAFTER_O3ORIG_D=${LAFTER_O3ORIG_D:-.TRUE.}
#([ $LSOCOL = .FALSE. ] || [ $LCHEM = .FALSE. ]) && LAFTER_CHEM1_D=.FALSE.
#export LAFTER_CHEM2_D=${LAFTER_CHEM2_D:-.TRUE.}
#([ $LSOCOL = .FALSE. ] || [ $LCHEM = .FALSE. ]) && LAFTER_CHEM2_D=.FALSE.
export RAW_FILETYPE_ECHAM=${RAW_FILETYPE_ECHAM:-1}
export RAW_FILETYPE_CHEM=${RAW_FILETYPE_CHEM:-1}
if [ $RAW_FILETYPE_ECHAM -ne 1 ]; then
    LAFTER_ECHAM_H=.FALSE.
    LAFTER_ECHAM_C=.FALSE.
fi
if [ $RAW_FILETYPE_CHEM -ne 1 ]; then
    LAFTER_CHEM1_H=.FALSE.
    LAFTER_CHEM1_C=.FALSE.
    LAFTER_CHEM2_H=.FALSE.
    LAFTER_CHEM2_C=.FALSE.
fi  
SMSNOTIFICATION=${SMSNOTIFICATION:-.FALSE.}

echo "Submitting Job $EXPNO for Year $RESTARTYEAR month $RESTARTMONTH until Year $EXPERIMENTEYEAR month $EXPERIMENTEMONTH"
echo "  Number of CPUs: $NCPUS"
echo "  Requested wall clock time per job: $BRUTUS_RUNTIME_MODEL hours"


#####################################
# Check time:
#####################################

# Check provided date:
if [ ${EXPERIMENTEYEAR} -lt ${RESTARTYEAR} ]; then
    echo " --> RESTARTYEAR cannot be lower than EXPERIMENTEYEAR! Stop."
    exit
fi 

if ([ ${EXPERIMENTEYEAR} -eq ${RESTARTYEAR} ] \
    && [ ${EXPERIMENTEMONTH} -lt ${RESTARTMONTH} ]); then
    echo " --> EXPERIMENTEMONTH cannot be lower than RESTARTMONTH! Stop."
    exit
fi


###############################################
# Read personnel data of model user:
###############################################

USERDATA=~/echam5/user_data.bash
source ${USERDATA}


#################################################
# Path settings:
#################################################

# Path to code folder:
SOCOLCODESDIR=~/echam5/codes


# Path and filenames for run scripts:
export RUNSCRIPTDIR=~/echam5/run_scripts
export RUNSCRIPT=run_socol_${CTRLLAB_RUNSCRIPT}_refc2.bash
export POSTSCRIPT=post_socol_${CTRLLAB_POSTSCRIPT}_refc2.bash
export CHECKRUNTIMESCRIPT=check_runtime_socol.bash
export KILLANDRESTARTSCRIPT=kill_and_restart_socol.bash

# Path and filenames for afterburner scripts:
export BURNSCRIPTDIR=~/echam5/burn

# Path and filenames for namelists:
NAMELISTDIR=~/echam5/namelists
export ECHAM_NAMELIST=echam_${CTRLLAB_ECHAM_NAMELIST}.nml
#export AFTERLEVELS_CHEM_M_NAMELIST=afterlevels_chem_m_${CTRLLAB_AFTERLEVELS_CHEM_M_NAMELIST}.nml
export AFTER_ECHAM_H_NAMELIST=after_echam_d_${CTRLLAB_AFTER_ECHAM_H_NAMELIST}_hybrid.nml
export AFTER_ECHAM_C_NAMELIST=after_echam_d_${CTRLLAB_AFTER_ECHAM_C_NAMELIST}_ccmi.nml
export AFTER_CHEM1_H_NAMELIST=after_chem1_d_${CTRLLAB_AFTER_CHEM1_H_NAMELIST}_hybrid.nml
export AFTER_CHEM1_C_NAMELIST=after_chem1_d_${CTRLLAB_AFTER_CHEM1_C_NAMELIST}_ccmi.nml
export AFTER_CHEM2_H_NAMELIST=after_chem2_d_${CTRLLAB_AFTER_CHEM2_H_NAMELIST}_hybrid.nml
export AFTER_CHEM2_C_NAMELIST=after_chem2_d_${CTRLLAB_AFTER_CHEM2_C_NAMELIST}_ccmi.nml

# Path and filenames for boundary condition links:
BCONDLINKSDIR=~/echam5/bcond_links
export BCOND_LINKS=bcond_links_${CTRLLAB_BCOND_LINKS}.bash

# Path to ECHAM folder:
export WORKDIR=/cluster/work/uwis/stenkea/echam5
export WORKDIR2=/cluster/scratch_xl/public/${BRUTUS_LOGINNAME}/echam5
export WORKDIR3=~revelll/ccmi_data

# Model executable (model binary):
export MODELBIN=$WORKDIR2/bin/$EXECUTABLENAME

# Path to boundary conditions and chemical initial data files:
export BOUNDARYCOND=$WORKDIR/data
export BOUNDARYCONDAMIP=$WORKDIR/data/T${RES}/amip2
export BOUNDARYCONDSOCOL=$WORKDIR/data/SOCOL/T${RES}
export CHEMINITIALDIR=$WORKDIR/data/SOCOL/T${RES}/chem_initial

# Temporary output before sending to HSM/local PC...

# Path for raw model output and restart files on Cluster HSM system:
#export HSM_DIR=/UWIS/iactrans/HSM/brutus/${BRUTUS_LOGINNAME}/raw/run$EXPNO
export RUNDIR=$WORKDIR2/run/run$EXPNO
# ...grouped in the following subdirectories:
# a) Model output (before afterburner):
export MODELOUTPUTDIR=$RUNDIR/modeloutput
# b) Temporary output of netCDF-files (after afterburner):
export AFTERPATH=$RUNDIR/modeloutput
# c) Restart files:
export RESTARTDIR=$RUNDIR/restartfiles
# d) Log files:
export LOGDIR=$RUNDIR/logfiles
# e) Scripts and namelists:
export MODELSCRIPTSDIR=$RUNDIR/modelscripts
# f) Model code:
export MODELCODEDIR=$RUNDIR/modelcode
# g) Model executable (model binary):
export MODELBINDIR=$RUNDIR/bin
# h) Run folder (rest):
export RUNMODELDIR=$RUNDIR/run

# HSM sub-folders:
#export HSM_RESTART=${HSM_DIR}/restart
#export HSM_DAILY_GRIB=${HSM_DIR}/daily_grib
#export HSM_CHEM_M_BACKUP=${HSM_DIR}/chem_m_backup

# Log files:
export ERR_MES=${LOGDIR}/${EXPNO}_error_message.txt 
export RUN_LOG=${LOGDIR}/${EXPNO}_run.log
export RUN_JOB_NUMBERS_LOG=${LOGDIR}/${EXPNO}_run_job_numbers.log


#####################################
# Make some default directories:
#####################################

# Directories at BRUTUS:
[ -d $RUNDIR ] || mkdir -p $RUNDIR
[ -d $MODELOUTPUTDIR ] || mkdir -p $MODELOUTPUTDIR
[ -d $AFTERPATH ] || mkdir -p $AFTERPATH
[ -d $RESTARTDIR ] || mkdir -p $RESTARTDIR
[ -d $LOGDIR ] || mkdir -p $LOGDIR
[ -d $MODELSCRIPTSDIR ] || mkdir -p $MODELSCRIPTSDIR
#[ -d $MODELCODEDIR ] || mkdir -p $MODELCODEDIR
[ -d $MODELBINDIR ] || mkdir -p $MODELBINDIR
[ -d $RUNMODELDIR ] || mkdir -p $RUNMODELDIR

# Directories at HSM system:
#[ -d $HSM_DIR ] || mkdir -p $HSM_DIR
#[ -d $HSM_RESTART ] || mkdir -p $HSM_RESTART
#[ -d $HSM_DAILY_GRIB ] || mkdir -p $HSM_DAILY_GRIB 
#([ -d $HSM_CHEM_M_BACKUP ] || [ ${LSOCOL} = .FALSE. ] \
#    || [ ${LCHEM} = .FALSE. ]) || mkdir -p $HSM_CHEM_M_BACKUP

# Directories at local PC:
ok=$(/usr/bin/ssh ${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME} [ -d $LOCALPC_DIR ]; \
    echo $?)
[ $ok -ne 0 ] && /usr/bin/ssh ${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME} \
    mkdir -p $LOCALPC_DIR

ok=$(/usr/bin/ssh ${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME} \
    [ -d $LOCALPC_DIR/daily ]; echo $?)
[ $ok -ne 0 ] && /usr/bin/ssh ${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME} \
    mkdir -p $LOCALPC_DIR/daily

ok=$(/usr/bin/ssh ${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME} \
    [ -d $LOCALPC_DIR/daily ]; echo $?)
[ $ok -ne 0 ] && /usr/bin/ssh ${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME} \
    mkdir -p $LOCALPC_DIR/daily_grib

ok=$(/usr/bin/ssh ${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME} \
    [ -d $LOCALPC_DIR/daily ]; echo $?)
[ $ok -ne 0 ] && /usr/bin/ssh ${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME} \
    mkdir -p $LOCALPC_DIR/restart

##############################################
# Check availability of scripts and namelists:
##############################################

if ([ ! -s ${RUNSCRIPTDIR}/${RUNSCRIPT} ] \
    && [ ! -s ${MODELSCRIPTSDIR}/${RUNSCRIPT} ]); then
    error_message=" --> ${RUNSCRIPT} is missing. Stop"
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job $RESTARTYEAR $RESTARTMONTH

fi
if ([ ! -s ${RUNSCRIPTDIR}/${CHECKRUNTIMESCRIPT} ] \
    && [ ! -s ${MODELSCRIPTSDIR}/${CHECKRUNTIMESCRIPT} ]); then
    error_message=" --> ${CHECKRUNTIMESCRIPT} is missing. Stop"
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job $RESTARTYEAR $RESTARTMONTH
fi
if ([ ! -s ${RUNSCRIPTDIR}/${KILLANDRESTARTSCRIPT} ] \
    && [ ! -s ${MODELSCRIPTSDIR}/${KILLANDRESTARTSCRIPT} ]); then
    error_message=" --> ${KILLANDRESTARTSCRIPT} is missing. Stop"
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job $RESTARTYEAR $RESTARTMONTH
fi
if ([ ! -s ${RUNSCRIPTDIR}/${POSTSCRIPT} ] \
    && [ ! -s ${MODELSCRIPTSDIR}/${POSTSCRIPT} ]); then
    error_message=" --> ${POSTSCRIPT} is missing. Stop"
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job $RESTARTYEAR $RESTARTMONTH
fi
if ([ ! -s ${NAMELISTDIR}/${ECHAM_NAMELIST} ] \
    && [ ! -s ${MODELSCRIPTSDIR}/${ECHAM_NAMELIST} ]); then
    error_message=" --> Namelist ${ECHAM_NAMELIST} is missing. Stop"
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job $RESTARTYEAR $RESTARTMONTH
fi
#if ([ ! -s ${NAMELISTDIR}/${AFTERLEVELS_CHEM_M_NAMELIST} ] \
#    && [ ! -s ${MODELSCRIPTSDIR}/${AFTERLEVELS_CHEM_M_NAMELIST} ]); then
#    error_message=" --> Namelist ${AFTERLEVELS_CHEM_M_NAMELIST} is missing. Stop"
#    echo ${error_message} >> ${RUN_LOG}
#    echo ${error_message} > ${ERR_MES}
#    kill_job $RESTARTYEAR $RESTARTMONTH
#fi
if ([ ! -s ${NAMELISTDIR}/${AFTER_ECHAM_H_NAMELIST} ] \
    && [ ! -s ${MODELSCRIPTSDIR}/${AFTER_ECHAM_H_NAMELIST} ]); then
    error_message=" --> Namelist ${AFTER_ECHAM_H_NAMELIST} is missing. Stop"
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job $RESTARTYEAR $RESTARTMONTH
fi
if ([ ! -s ${NAMELISTDIR}/${AFTER_ECHAM_C_NAMELIST} ] \
    && [ ! -s ${MODELSCRIPTSDIR}/${AFTER_ECHAM_C_NAMELIST} ]); then
    error_message=" --> Namelist ${AFTER_ECHAM_C_NAMELIST} is missing. Stop"
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job $RESTARTYEAR $RESTARTMONTH
fi
#if ([ ! -s ${NAMELISTDIR}/${AFTER_CHEM1_H_NAMELIST} ] \
#    && [ ! -s ${MODELSCRIPTSDIR}/${AFTER_CHEM1_H_NAMELIST} ]); then
#    error_message=" --> Namelist ${AFTER_CHEM1_H_NAMELIST} is missing. Stop"
#    echo ${error_message} >> ${RUN_LOG}
#    echo ${error_message} > ${ERR_MES}
#    kill_job $RESTARTYEAR $RESTARTMONTH
#fi
#if ([ ! -${NAMELISTDIR}/${AFTER_CHEM1_C_NAMELIST} ] \
#    && [ ! -s ${MODELSCRIPTSDIR}/${AFTER_CHEM1_C_NAMELIST} ]); then
#    error_message=" --> Namelist ${AFTER_CHEM1_C_NAMELIST} is missing. Stop"
#    echo ${error_message} >> ${RUN_LOG}
#    echo ${error_message} > ${ERR_MES}
#    kill_job $RESTARTYEAR $RESTARTMONTH
#fi
#if ([ ! -s ${NAMELISTDIR}/${AFTER_CHEM2_H_NAMELIST} ] \
#    && [ ! -s ${MODELSCRIPTSDIR}/${AFTER_CHEM2_H_NAMELIST} ]); then
#    error_message=" --> Namelist ${AFTER_CHEM2_H_NAMELIST} is missing. Stop"
#    echo ${error_message} >> ${RUN_LOG}
#    echo ${error_message} > ${ERR_MES}
#    kill_job $RESTARTYEAR $RESTARTMONTH
#fi
#if ([ ! -${NAMELISTDIR}/${AFTER_CHEM2_C_NAMELIST} ] \
#    && [ ! -s ${MODELSCRIPTSDIR}/${AFTER_CHEM2_C_NAMELIST} ]); then
#    error_message=" --> Namelist ${AFTER_CHEM2_C_NAMELIST} is missing. Stop"
#    echo ${error_message} >> ${RUN_LOG}
#    echo ${error_message} > ${ERR_MES}
#    kill_job $RESTARTYEAR $RESTARTMONTH
#fi
if ([ ! -s ${BCONDLINKSDIR}/${BCOND_LINKS} ] \
    && [ ! -s ${MODELSCRIPTSDIR}/${BCOND_LINKS} ]); then
    error_message=" --> ${BCOND_LINKS} missing. Stop"
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job $RESTARTYEAR $RESTARTMONTH
fi
if [ ! -s ${CHEMINITIALDIR}/${CHEMINITIAL_FILE} ]; then
    error_message=" --> ${CHEMINITIAL_FILE} is missing. Stop"
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job $RESTARTYEAR $RESTARTMONTH
fi
if ([ ! -s ${WORKDIR2}/bin/${EXECUTABLENAME} ] \
    && [ ! -s ~/socolv3_ccmi/bin/${EXECUTABLENAME} ] \
    && [ ! -s ${MODELBINDIR}/${EXECUTABLENAME} ]); then
    error_message=" --> bin/${EXECUTABLENAME} is missing. Stop"
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job $RESTARTYEAR $RESTARTMONTH
fi


###############################################
# Copy code, scripts and executable to run folder:
###############################################

# Copy model code (if available) to run folder for documentation:
if ([ -d ${SOCOLCODESDIR}/${CODENAME} ] \
    && [ ! -d ${MODELCODEDIR}/${CODENAME} ]) ; then
    mkdir ${MODELCODEDIR}/${CODENAME}
    cp -rp ${SOCOLCODESDIR}/${CODENAME}/include \
	${SOCOLCODESDIR}/${CODENAME}/modules \
	${SOCOLCODESDIR}/${CODENAME}/src ${MODELCODEDIR}/${CODENAME}
fi

# Copy control-scripts and namelists to run folder for documentation and 
# execution:
/bin/mv -f $(pwd)/*${EXPNO}* $MODELSCRIPTSDIR 2> /dev/null
/bin/mv -f ${NAMELISTDIR}/*${EXPNO}* $MODELSCRIPTSDIR 2> /dev/null
/bin/mv -f ${BCONDLINKSDIR}/*${EXPNO}* $MODELSCRIPTSDIR 2> /dev/null
[ -s ${RUNSCRIPTDIR}/$RUNSCRIPT ] \
    && /bin/cp -p ${RUNSCRIPTDIR}/$RUNSCRIPT $MODELSCRIPTSDIR
[ -s ${RUNSCRIPTDIR}/$POSTSCRIPT ] \
    && /bin/cp -p ${RUNSCRIPTDIR}/$POSTSCRIPT $MODELSCRIPTSDIR
[ -s ${NAMELISTDIR}/$ECHAM_NAMELIST ] \
    && /bin/cp -p ${NAMELISTDIR}/$ECHAM_NAMELIST $MODELSCRIPTSDIR
#[ -s ${NAMELISTDIR}/$AFTERLEVELS_CHEM_M_NAMELIST ] \
#    && /bin/cp -p ${NAMELISTDIR}/$AFTERLEVELS_CHEM_M_NAMELIST $MODELSCRIPTSDIR
([ $LAFTER_ECHAM_H = .TRUE. ] \
    && [ -s ${NAMELISTDIR}/$AFTER_ECHAM_H_NAMELIST ]) \
    && /bin/cp -p ${NAMELISTDIR}/$AFTER_ECHAM_H_NAMELIST $MODELSCRIPTSDIR
([ $LAFTER_ECHAM_C = .TRUE. ] \
    && [ -s ${NAMELISTDIR}/$AFTER_ECHAM_C_NAMELIST ]) \
    && /bin/cp -p ${NAMELISTDIR}/$AFTER_ECHAM_C_NAMELIST $MODELSCRIPTSDIR
#([ $LAFTER_CHEM1_H = .TRUE. ] \
#    && [ -s ${NAMELISTDIR}/$AFTER_CHEM1_H_NAMELIST ]) \
#    && /bin/cp -p ${NAMELISTDIR}/$AFTER_CHEM1_H_NAMELIST $MODELSCRIPTSDIR
#([ $LAFTER_CHEM2_C = .TRUE. ] \
#    && [ -s ${NAMELISTDIR}/$AFTER_CHEM2_D_NAMELIST ]) \
#    && /bin/cp -p ${NAMELISTDIR}/$AFTER_CHEM2_D_NAMELIST $MODELSCRIPTSDIR
#([ $LAFTER_O3ORIG_D = .TRUE. ] \
#    && [ -s ${NAMELISTDIR}/$AFTER_O3ORIG_D_NAMELIST ]) \
#    && /bin/cp -p ${NAMELISTDIR}/$AFTER_O3ORIG_D_NAMELIST $MODELSCRIPTSDIR
[ -s ${BCONDLINKSDIR}/$BCOND_LINKS ] \
    && /bin/cp -p ${BCONDLINKSDIR}/$BCOND_LINKS $MODELSCRIPTSDIR

# Move model executable (binary file) to run folder for executation:
[ -s $WORKDIR2/bin/$EXECUTABLENAME ] \
    && /bin/cp -p $WORKDIR2/bin/$EXECUTABLENAME $MODELBINDIR

[ -s ~/socolv3_ccmi/bin/$EXECUTABLENAME ] \
    && /bin/cp -p ~/socolv3_ccmi/bin/$EXECUTABLENAME $MODELBINDIR

####################################
# Clean up RUNSCRIPTDIR, BURNSCRIPTDIR, NAMELISTDIR and BCONDLINKSDIR:
####################################

/bin/rm -f "${BURNSCRIPTDIR}/*_${EXPNO}*" 2> /dev/null
/bin/rm -f "${NAMELISTDIR}/*_${EXPNO}*" 2> /dev/null
/bin/rm -f "${BCONDLINKSDIR}/*_${EXPNO}*" 2> /dev/null
/bin/rm -f "${WORKDIR2}/bin/*_${EXPNO}*" 2> /dev/null


#####################################
# Submit job:
#####################################

YYYY=$RESTARTYEAR
M=$RESTARTMONTH

while ([ $YYYY -lt $EXPERIMENTEYEAR ] || \
    ([ $YYYY -eq $EXPERIMENTEYEAR ] && [ $M -le $EXPERIMENTEMONTH ])); do

    export Y_S=$YYYY      # Start year of loop in SOCOL script
    export M_S=$M         # Start month of start year of loop in SOCOL script

    export Y_E=$YYYY      # End year of loop in SOCOL script
    export M_E=$((M+BRUTUS_JOB_NMONTH-1))  # End month of end year

    # Correct M_E/Y_E if necessary:
    while [ $M_E -gt 12 ]; do
	M_E=$((M_E-12))
	Y_E=$((Y_E+1))
    done
    if [ $Y_E -gt $EXPERIMENTEYEAR ]; then
	Y_E=$EXPERIMENTEYEAR
	M_E=$EXPERIMENTEMONTH
    fi
    ([ $Y_E -eq $EXPERIMENTEYEAR ] && [ $M_E -gt $EXPERIMENTEMONTH ]) \
	&& M_E=$EXPERIMENTEMONTH

    [ $M_S -le 9 ] && MM_S=0$M_S || MM_S=$M_S
    [ $M_E -le 9 ] && MM_E=0$M_E || MM_E=$M_E
    
    # Submit job:
    ([ $Y_S -eq $Y_E ] && [ $M_S -eq $M_E ])\
	&& log_message="Submitting job for $Y_S/$MM_S." \
	|| log_message="Submitting job for period $Y_S/$MM_S - $Y_E/$MM_E."
    
    if [ $Y_S -eq $RESTARTYEAR ] && [ $M_S -eq $RESTARTMONTH ]; then
	
	# Beginning of run:

	# Previous month/year:
	if [ ${M_S} -ne 1 ]; then
	    M_EB=$((M_S-1))
	    Y_EB=$Y_S
	else
	    M_EB=12
	    Y_EB=$((Y_S-1))
	fi
	[ ${M_EB} -le 9 ] && MM_EB=0${M_EB} || MM_EB=${M_EB}

        # Check-File:
	echo '0' > ${LOGDIR}/${EXPNO}_run_${Y_EB}${MM_EB}_ok

	echo ${log_message} | tee ${RUN_JOB_NUMBERS_LOG}
	
	bsub -J chain_${EXPNO} -n $NCPUS -W $BRUTUS_RUNTIME_MODEL -o ${LOGDIR}/${EXPNO}_run_socol_script.log < ${MODELSCRIPTSDIR}/${RUNSCRIPT} | tee -a ${RUN_JOB_NUMBERS_LOG}
	
    else

	# Continuation of run:

	echo ${log_message} | tee -a ${RUN_JOB_NUMBERS_LOG}

	bsub -J chain_${EXPNO} -n $NCPUS -W $BRUTUS_RUNTIME_MODEL -o ${LOGDIR}/${EXPNO}_run_socol_script.log -w chain_${EXPNO} < ${MODELSCRIPTSDIR}/${RUNSCRIPT} | tee -a ${RUN_JOB_NUMBERS_LOG}

    fi

    # Increase month/year:
    if [ $M_E -le 11 ]; then
	M=$((M_E+1))
	YYYY=$Y_E
    else
	M=1
	YYYY=$((Y_E+1))
    fi

done    
    

exit
