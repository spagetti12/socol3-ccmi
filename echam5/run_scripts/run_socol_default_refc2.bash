#!/bin/bash

# Performs SOCOL simulation and calls *post_socol_default.bash*.
# *run_socol_default.bash*  is called from *call_socol_default.bash*. 

# Martin Schraner, ETH Zuerich, February 2010

###############################################################################


#################################
# Additional tasks at BRUTUS (1):
#################################

[ -f ${RUN_LOG} ] \
    && echo "==========================================" >> ${RUN_LOG} \
    || echo "==========================================" > ${RUN_LOG}
date >> ${RUN_LOG}
echo "Run from $HOSTNAME" >> ${RUN_LOG}
echo "Job start year: ${Y_S}" >> ${RUN_LOG}
echo "Job start month: ${M_S}" >> ${RUN_LOG}
echo "Job end year: ${Y_E}" >> ${RUN_LOG}
echo "Job end month: ${M_E}" >> ${RUN_LOG}
echo "Simulation number: ${EXPNO}" >> ${RUN_LOG}
echo "Starting model..." >> ${RUN_LOG}

HH_START_RUN=$(date +%H)         # Hour when program starts
MM_START_RUN=$(date +%M)         # Minutes when program starts

# Previous month/year:
if [ ${M_S} -ne 1 ]; then
    M_EB=$((M_S-1))
    Y_EB=$Y_S
else
    M_EB=12
    Y_EB=$((Y_S-1))
fi
[ ${M_EB} -le 9 ] && MM_EB=0${M_EB} || MM_EB=${M_EB}
[ ${M_E} -le 9 ] && MM_E=0${M_E} || MM_E=${M_E}

# Read in Check-File to see whether previous run has finished properly:
ok=$(cat ${LOGDIR}/${EXPNO}_run_${Y_EB}${MM_EB}_ok)
/bin/rm -f ${LOGDIR}/${EXPNO}_run_${Y_EB}${MM_EB}_ok

# Exit, if run of previous year has not finished properly:
if [ $ok -ne 0 ]; then
    error_message="ERROR: Previous month has not finished properly."
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job ${Y_S} ${M_S}
fi

# Set new default value for ${EXPNO}_run_${Y_E}${MM_E}_ok (changed at the end  
# of the script, if model run was sucessful):
echo "1" > ${LOGDIR}/${EXPNO}_run_${Y_E}${MM_E}_ok

################################################################################

##############################
# "Normal" run script (BEGIN):
##############################
ulimit -s unlimited

export MPLARGS="verbose=1;check=0"
export VPP_MBX_SIZE=6400000
export VPP_STATS=10
#
# Job file to run echam model on portlandgroup fortran compiler
# =============================================================
#
F_RECLUNIT=BYTE ; export F_RECLUNIT
MPIPROGINF=detail ; export MPIPROGINF
F_SYSLEN=600 ; export F_SYSLEN

# Some calculations:
Y_SM1=$((Y_S-1))
Y_EP1=$((Y_E+1))
EXPERIMENTEYEARP1=$((EXPERIMENTEYEAR+1))

cd $RUNMODELDIR
#
################################################################################
#
INI=$WORKDIR/data/T${RES}
INIAMIP=$WORKDIR/data/T${RES}/amip2
INIHADLEY=$WORKDIR/data/T${RES}/hadley
INISOCOL=$WORKDIR/data/SOCOL/T${RES}
INICCMI=$WORKDIR3/emissions/GHG_ODS
#
################################################################################
#
# specification of file structures
#
# stop execution after the first run time error
F_ERRCNT=0
export F_ERRCNT
#
################################################################################
#
rm -f unit.?? sst* ice* rrtadata *.codes atmout chem_initial
#
ln -s ${INI}/T${RES}L${LEV}_jan_spec.nc      unit.23
ln -s ${INI}/T${RES}_jan_surf.nc             unit.24
#
################################################################################
#
#  Initial file for chemical species (MEZON):

ln -s ${CHEMINITIALDIR}/${CHEMINITIAL_FILE}  chem_initial
#
################################################################################
#
ln -s  ${INI}/T${RES}_O3clim2.nc    unit.21
ln -s  ${INI}/T${RES}_VLTCLIM.nc    unit.90
ln -s  ${INI}/T${RES}_VGRATCLIM.nc  unit.91
ln -s  ${INI}/T${RES}_TSLCLIM2.nc   unit.92
ln -s  ${INI}/surrta_data           rrtadata
#
################################################################################
#
# Boundary conditions for SOCOL / SST/SIC:
source ${MODELSCRIPTSDIR}/${BCOND_LINKS}
#
################################################################################
#
# Prepare rerun files (if necessary):

M_S=${M_S#0}           # remove leading zero (if present)

if [ ${Y_S} -gt ${EXPERIMENTSYEAR} ] || [ ${M_S} -gt 1 ]; then
    if [ ${M_S} -eq 1 ]; then
	MM1=12
	YYYY1=${Y_SM1}
    else
	MM1=$((M_S-1))
	YYYY1=${Y_S}
    fi
    [ ${MM1} -le 9 ] && MM1=0${MM1}

    # a) ECHAM5:
    if [ ! -s ${RESTARTDIR}/${EXPNO}_rerun_${YYYY1}${MM1}_echam.nc.gz ]; then
	error_message="${EXPNO}_rerun_${YYYY1}${MM1}_echam.nc.gz does not exist!"
	echo ${error_message} >> ${RUN_LOG}
	echo ${error_message} > ${ERR_MES}
	kill_job ${Y_S} ${M_S}
    fi
    cp -f ${RESTARTDIR}/${EXPNO}_rerun_${YYYY1}${MM1}_echam.nc.gz rerun_${EXPNO}_echam.gz
    gunzip -f rerun_${EXPNO}_echam.gz

    if ([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ]); then 
    # b) chem1 (MEZON):
	if [ ! -s ${RESTARTDIR}/${EXPNO}_rerun_${YYYY1}${MM1}_chem1.nc.gz ]; then
	    error_message="${EXPNO}_rerun_${YYYY1}${MM1}_chem1.nc.gz does not exist!"
	    echo ${error_message} >> ${RUN_LOG}
	    echo ${error_message} > ${ERR_MES}
	    kill_job ${Y_S} ${M_S}
	fi
	cp -f ${RESTARTDIR}/${EXPNO}_rerun_${YYYY1}${MM1}_chem1.nc.gz rerun_${EXPNO}_tracer.gz
	gunzip -f rerun_${EXPNO}_tracer.gz

    # c) chem2 (MEZON): 
	if [ ! -s ${RESTARTDIR}/${EXPNO}_rerun_${YYYY1}${MM1}_chem2.nc.gz ]; then
	    error_message="${EXPNO}_rerun_${YYYY1}${MM1}_chem2.nc.gz does not exist!"
	    echo ${error_message} >> ${RUN_LOG}
	    echo ${error_message} > ${ERR_MES}
	    kill_job ${Y_S} ${M_S}
	fi
	cp -f ${RESTARTDIR}/${EXPNO}_rerun_${YYYY1}${MM1}_chem2.nc.gz rerun_${EXPNO}_chem2.gz
	gunzip -f rerun_${EXPNO}_chem2.gz

    # d) chem_m (MEZON):
    #    (no rerun fields included, but to be provided nevertheless):
#	if [ ! -s ${RESTARTDIR}/rerun_${EXPNO}_chem_m ]; then
#	    if [ -s $WORKDIR/rerun_chem/runnr_rerun_chem_m.nc ]; then
#		cp -f $WORKDIR/rerun_chem/runnr_rerun_chem_m.nc rerun_${EXPNO}_chem_m
#	    else
#		error_message="Neither ${EXPNO}_rerun_chem_m.nc nor $WORKDIR/rerun_chem/runnr_rerun_chem_m.nc exists!"
#		echo ${error_message} >> ${RUN_LOG}
#		echo ${error_message} > ${ERR_MES}
#		kill_job ${Y_S} ${M_S}
#	    fi
#	fi
    fi

    RERUN=.TRUE.     # Rerun switch; .false. for initial run, .true. for rerun
else
    RERUN=.FALSE.
fi

rm -f fortran_error_messages
#
################################################################################
#
# Set namelist and call model:

YYYY=$Y_S
M=$M_S

while ([ $YYYY -lt $Y_E ] || ([ $YYYY -eq $Y_E ] && [ $M -le $M_E ])); do

    [ $M -le 9 ] && MM=0$M || MM=$M

    export YYYY
    export MM

    # Check remaining run time and determine $runtime_ok
    source ${RUNSCRIPTDIR}/${CHECKRUNTIMESCRIPT}

    # Kill run jobs and restart if $runtime_ok=1
    [ $runtime_ok -eq 1 ] && source ${RUNSCRIPTDIR}/${KILLANDRESTARTSCRIPT}

    # Log files:
    SOCOLOUT_LOG=${LOGDIR}/${EXPNO}_socolout_${YYYY}${MM}.log
    SOCOLOUT_Y_LOG=${LOGDIR}/${EXPNO}_socolout_${YYYY}.log

    # Reset CO2FAC to 1 after RESTARTMONTH of RESTARTYEAR:
    ([ ${YYYY} -eq ${RESTARTYEAR} ] && [ $M -eq ${RESTARTMONTH} ]) \
	|| CO2FAC=1.0000 
		
    # ECHAM namelist:
    #  namelist control variables and output control for grid variables
    #  spectral variables are written out by default except liq. water
    #  for production runs set LABORT=.FALSE.
    [ -s namelist.echam ] && rm -f namelist.echam
    source ${MODELSCRIPTSDIR}/${ECHAM_NAMELIST}

    # Call model:
    echo "------------------------------------------" >> ${RUN_LOG}
    date >> ${RUN_LOG}
    echo "Calling SOCOL for Year $YYYY, month $MM..." >> ${RUN_LOG}
    ompirun ${MODELBINDIR}/${EXECUTABLENAME} >& ${SOCOLOUT_LOG}

    # Concatenate log files:
    if [ -s ${SOCOLOUT_Y_LOG} ]; then
	cat ${SOCOLOUT_LOG} >> ${SOCOLOUT_Y_LOG}
	rm -f ${SOCOLOUT_LOG}
    else
	mv -f ${SOCOLOUT_LOG} ${SOCOLOUT_Y_LOG}
    fi	    

    # Rename output files:
    #echam:
    mv -f ${EXPNO}_${YYYY}${MM}.01 ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}_echam
    mv -f ${EXPNO}_${YYYY}${MM}.01.codes ${MODELOUTPUTDIR}/${EXPNO}_echam.codes

    #mezon:
    mv -f ${EXPNO}_${YYYY}${MM}.01_tracer ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}_chem1
    mv -f ${EXPNO}_${YYYY}${MM}.01_chem1.nc ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}.01_chem1.nc
    mv -f ${EXPNO}_${YYYY}${MM}.01_chem2.nc ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}.01_chem2.nc
    mv -f ${EXPNO}_${YYYY}${MM}.01_tracer.codes ${MODELOUTPUTDIR}/${EXPNO}_chem1.codes
    # a) ECHAM5:
#    [ $RAW_FILETYPE_ECHAM -eq 1 ] && suffecham= || suffecham=.nc
#    mv -f ${EXPNO}_${YYYY}${MM}.01${suffecham} \
#	${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}_echam${suffecham}
#    [ $RAW_FILETYPE_ECHAM -eq 1 ] && \
#	mv -f ${EXPNO}_${YYYY}${MM}.01.codes ${MODELOUTPUTDIR}/${EXPNO}_echam.codes
#
#    # b) MEZON:
#    if ([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ]); then
#	[ $RAW_FILETYPE_CHEM -eq 1 ] && suffchem= || suffchem=.nc
#	mv -f ${EXPNO}_${YYYY}${MM}.01_tracer${suffchem} \
#	    ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}_chem1${suffchem}
#	mv -f ${EXPNO}_${YYYY}${MM}.01_chem1.nc \
#	    ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}.01_chem1.nc
#	mv -f ${EXPNO}_${YYYY}${MM}.01_chem2.nc \
#	    ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}.01_chem2.nc
#	#mv -f ${EXPNO}_${YYYY}${MM}.01_chem_m.nc \
#	#    ${MODELOUTPUTDIR}/${EXPNO}_chem_m_L${LEV}_${YYYY}${MM}.nc
#	if [ $RAW_FILETYPE_ECHAM -eq 1 ]; then
#	    mv -f ${EXPNO}_${YYYY}${MM}.01_tracer.codes \
#		${MODELOUTPUTDIR}/${EXPNO}_chem1.codes
#	 #   mv -f ${EXPNO}_${YYYY}${MM}.01_chem2.codes \
#	#	${MODELOUTPUTDIR}/${EXPNO}_chem2.codes
#	fi
#    fi

    echo "Renamed output files" >> ${RUN_LOG}

    # c) O3ORIG:
#    [ $RAW_FILETYPE_ECHAM -eq 1 ] && suffecham= || suffecham=.nc
#    mv -f ${EXPNO}_${YYYY}${MM}.01_o3orig${suffecham} \
#	${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}_o3orig${suffecham}
#    [ $RAW_FILETYPE_ECHAM -eq 1 ] && \
#	mv -f ${EXPNO}_${YYYY}${MM}.01_o3orig.codes ${MODELOUTPUTDIR}/${EXPNO}_o3orig.codes

    # Save rerun files:
    # a) ECHAM5:
    cp -f rerun_${EXPNO}_echam ${RESTARTDIR}/${EXPNO}_rerun_${YYYY}${MM}_echam.nc
    zipecham_ok=$(gzip -f ${RESTARTDIR}/${EXPNO}_rerun_${YYYY}${MM}_echam.nc; echo $?)
    
    if ([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ]); then
    # b) chem1 (MEZON):
	cp -f rerun_${EXPNO}_tracer ${RESTARTDIR}/${EXPNO}_rerun_${YYYY}${MM}_chem1.nc
	zipchem1_ok=$(gzip -f ${RESTARTDIR}/${EXPNO}_rerun_${YYYY}${MM}_chem1.nc; echo $?)
	
    # c) chem2 (MEZON):
	cp -f rerun_${EXPNO}_chem2 ${RESTARTDIR}/${EXPNO}_rerun_${YYYY}${MM}_chem2.nc
	zipchem2_ok=$(gzip -f ${RESTARTDIR}/${EXPNO}_rerun_${YYYY}${MM}_chem2.nc; echo $?)
    else
	zipchem1_ok=0
	zipchem2_ok=0
    fi
    
    RERUN=.TRUE.

    echo "Saved rerun files" >> ${RUN_LOG}
    

    # "Normal" run script (END)

################################################################################

    #################################
    # Additional tasks at BRUTUS (2):
    #################################


    # Check success of run:
    # ---------------------

    model_ok=0
    
    # Monthly mean chemistry output (nc):
#    if ([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ]); then
#	[ -e ${MODELOUTPUTDIR}/${EXPNO}_chem_m_L${LEV}_${YYYY}${MM}.nc ] \
#	    || model_ok=1
#    fi
    
    # 12h output data:
#    [ -e ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}_echam${suffecham} ] || model_ok=2
#    if ([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ]); then
#	[ -e ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}.01_chem1${suffchem} ] \
#	    || model_ok=3
#	[ -e ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}.01_chem2${suffchem} ] \
#	    || model_ok=4
#    fi

    echo "12h output data ok" >> ${RUN_LOG}
    
    # Rerun files:
    ([ $zipecham_ok -eq 0 ] && [ $zipchem1_ok -eq 0 ] \
	&& [ $zipchem2_ok -eq 0 ]) || model_ok=5
    [ -e ${RESTARTDIR}/${EXPNO}_rerun_${YYYY}${MM}_echam.nc.gz ] \
	|| model_ok=6
    if ([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ]); then
	[ -e ${RESTARTDIR}/${EXPNO}_rerun_${YYYY}${MM}_chem1.nc.gz ] \
	    || model_ok=7
	[ -e ${RESTARTDIR}/${EXPNO}_rerun_${YYYY}${MM}_chem2.nc.gz ] \
	    || model_ok=8
    fi
    
    echo "------------------------------------------" >> ${RUN_LOG}
    date >> ${RUN_LOG}
    echo "Status of model output: $model_ok" >> ${RUN_LOG}
    
    # Exit if something went wrong with the model output:
    if [ ${model_ok} -ne 0 ]; then 
	error_message="ERROR: Something went wrong with the run of Year ${YYYY} and Month ${MM}. Not all Output-Files are stored on ${MODELOUTPUTDIR}." >> ${RUN_LOG}
	echo ${error_message} >> ${RUN_LOG}
	echo ${error_message} > ${ERR_MES}
	kill_job ${YYYY} ${MM}
    fi
    

    # Call postprocessing:
    # --------------------

    if ([ $RAW_FILETYPE_ECHAM -eq 1 ] && [ $RAW_FILETYPE_CHEM -eq 1 ]); then 
	
	echo "------------------------------------------" >> ${RUN_LOG}
	date  >> ${RUN_LOG}
	echo "Calling postprocessing..." >> ${RUN_LOG}
	
	bsub -J chain_${EXPNO}_${YYYY}${MM}_post -n 1 -R"select[ib]" \
	    -W ${BRUTUS_RUNTIME_POST} \
	    -o ${LOGDIR}/${EXPNO}_post_socol_script.log \
	    < ${MODELSCRIPTSDIR}/${POSTSCRIPT}	
	
	echo "------------------------------------------" >> ${RUN_LOG}
	date  >> ${RUN_LOG}
	echo "Model month ${MM} finished " >> ${RUN_LOG}
	echo "==========================================" >> ${RUN_LOG}
	
    fi

    # Increase month/year:
    if [ $M -le 11 ]; then
	M=$((M+1))
    else
	echo "Model year ${YYYY} finished " >> ${RUN_LOG}
	echo "==========================================" >> ${RUN_LOG}
	M=1
	YYYY=$((YYYY+1))
    fi
    
done

# Model run sucessful -> set ${EXPNO}_run_${Y_E}${MM_E}_ok to 0:
echo "0" > ${LOGDIR}/${EXPNO}_run_${Y_E}${MM_E}_ok

exit
