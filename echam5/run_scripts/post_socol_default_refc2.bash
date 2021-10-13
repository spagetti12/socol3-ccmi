#!/bin/bash

# Performs postprocessing of SOCOL output.
# *post_socol_default.bash*  is called from *run_socol_default.bash*.

# Martin Schraner, ETH Zuerich, February 2010

###############################################################################


# Log files:

POST_LOG=${LOGDIR}/${EXPNO}_post_${YYYY}.log
    
[ ! -e $POST_LOG ] \
    && echo "==========================================" > $POST_LOG \
    || echo "==========================================" >> $POST_LOG
date >> $POST_LOG
echo "Simulation number: ${EXPNO}" >> $POST_LOG
echo "Current year: ${YYYY}" >> $POST_LOG
echo "Current month: ${MM}" >> $POST_LOG
echo "Starting postprocessing and storing output..." \
    >> $POST_LOG

M=${MM#0}           # remove leading zero (if present)

cd ${MODELOUTPUTDIR}


#1.1 Postprocessing of grib files:
###################################

echo "------------------------------------------" >> $POST_LOG
date >> ${POST_LOG}
echo "Postprocessing of echam, chem1, and chem2 grib files..." >> ${POST_LOG}

${BURNSCRIPTDIR}/afterd_afterm_echam5_default.bash >> $POST_LOG

# Check status:
post_grib_ok=0

if [ ${LAFTER_ECHAM_H} = .TRUE. ]; then
    [ -s ${AFTERPATH}/${EXPNO}_echam_h_${YYYY}${MM}.nc ] || post_grib_ok=1
fi
if [ ${LAFTER_ECHAM_C} = .TRUE. ]; then
    [ -s ${AFTERPATH}/${EXPNO}_echam_c_${YYYY}${MM}.nc ] || post_grib_ok=2
fi
if ([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ] \
    && [ ${LAFTER_CHEM1_H} = .TRUE. ]); then
    [ -s ${AFTERPATH}/${EXPNO}_chem1_h_${YYYY}${MM}.nc ] || post_grib_ok=3
fi
if ([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ] \
    && [ ${LAFTER_CHEM1_C} = .TRUE. ]); then
    [ -s ${AFTERPATH}/${EXPNO}_chem1_c_${YYYY}${MM}.nc ] || post_grib_ok=4
fi
if ([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ] \
    && [ ${LAFTER_CHEM2_H} = .TRUE. ]); then
    [ -s ${AFTERPATH}/${EXPNO}_chem2_h_${YYYY}${MM}.nc ] || post_grib_ok=5
fi
if ([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ] \
    && [ ${LAFTER_CHEM2_C} = .TRUE. ]); then
    [ -s ${AFTERPATH}/${EXPNO}_chem2_c_${YYYY}${MM}.nc ] || post_grib_ok=6
fi
# Check status:
echo "Status of postprocessing grib files: $post_grib_ok" >> ${POST_LOG}

# Warnings:
if [ ${post_grib_ok} -ne 0 ]; then
    error_message="WARNING: Something went wrong with postprocessing of grib files of Job ${EXPNO} for Year ${YYYY} and Month ${MM}"
    echo ${error_message} >> ${POST_LOG}
    echo ${error_message} > ${ERR_MES}
    /bin/mail -s "Postprocessing problems at BRUTUS" ${MAIL_NOTIFICATION} \
	< ${ERR_MES}
    exit
fi



# 2. Uploading files to HSM system and local PC:
################################################
   
# 2.1 Uploading grib-output files to local PC:
################################################
#
echo "------------------------------------------" >> ${POST_LOG}
date >> ${POST_LOG}
echo "Uploading GRIB-output files to local PC..." >> ${POST_LOG}

list=${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}_echam
([ ${LSOCOL} = .TRUE. ] && [ ${LCHEM} = .TRUE. ]) \
    && list="$list ${MODELOUTPUTDIR}/${EXPNO}_${YYYY}${MM}_chem?"

for file in $list; do

    # Zip file:
    status_ok=$(/bin/gzip ${file}; echo $?)
    echo "Status of zipping ${file}: $status_ok" >> ${POST_LOG}
    file=${file}.gz

    # Upload file:
    [ target=${LOCALPC_DIR} || target=${LOCALPC_DIR} 
    status_ok=$(scp -p ${file} \
        ${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME}:${target} \
        > /dev/null; echo $?)
    echo "Status of uploading ${file} to local PC: $status_ok" >> ${POST_LOG}

    # Warnings and removing file:
    if [ $status_ok -ne 0 ]; then
	error_message="WARNING: Something went wrong with uploading ${file} to local PC"
	echo ${error_message} >> ${POST_LOG}
	echo ${error_message} > ${ERR_MES}
	/bin/mail -s "Uploading problems at BRUTUS" ${MAIL_NOTIFICATION} \
	    < ${ERR_MES}
    else
	/bin/rm -f ${file}
    fi
    
done


# 2.2 Uploading restart files to LOCALPC:
#########################################

cd ${RESTARTDIR}

if ([ ${YYYY} -ne ${EXPERIMENTSYEAR} ] && [ ${M} -eq 1 ]); then

    echo "------------------------------------------" >> ${POST_LOG}
    date >> ${POST_LOG}
    echo "Uploading restart files (December) to local PC..." >> ${POST_LOG}


    # Make Archive-Files of December containing Restart-Files of Chemistry +
    # Dynamics:
    YYYYb=$((YYYY-1))
    MMb=12
    status_ok=$(/bin/tar -cvf \
        ${EXPNO}_rerun_${YYYYb}${MMb}.tar \
        ${EXPNO}_rerun_${YYYYb}${MMb}* > /dev/null; echo $?)
        echo "Status of restart-tar-file: $status_ok" >> ${POST_LOG}

    # Upload file:
    [ target=${LOCALPC_DIR}/restart || target=${LOCALPC_DIR}
    status_ok=$(scp -p ${EXPNO}_rerun_${YYYYb}${MMb}.tar \
        ${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME}:${target} \
        > /dev/null; echo $?)
    echo "Status of uploading ${EXPNO}_rerun_${YYYYb}${MMb}.tar to local PC: $status_ok" >> ${POST_LOG}

    # Warnings and removing file:
    if [ $status_ok -ne 0 ]; then
        error_message="WARNING: Something went wrong with uploading ${EXPNO}_rerun_${YYYYb}${MMb}.tar to local PC."
        echo ${error_message} >> ${POST_LOG}
        echo ${error_message} > ${ERR_MES}
        /bin/mail -s "Uploading problems at BRUTUS" ${MAIL_NOTIFICATION} \
            < ${ERR_MES}
    else
        /bin/rm -f ${file}
    fi    

elif [ ${M} -ne 1 ]; then

   # Remove restart files for January-November:
    MMb=$((M-1))
    [ ${MMb} -le 9 ] && MMb=0${MMb}
    /bin/rm -f ${EXPNO}_rerun_${YYYY}${MMb}*

fi


# 2.3 Uploading netCDF-files to local PC:
#########################################

cd ${MODELOUTPUTDIR}

echo "------------------------------------------" >> ${POST_LOG}
date >> ${POST_LOG}
echo "Uploading netCDF-files to local PC..." >> ${POST_LOG}

list=${AFTERPATH}/${EXPNO}*_[mdch]_${YYYY}${MM}.nc

for file in $list; do

    # Zip, if daily file:
    if [ "${file}" != "${file/_d_}" ]; then
	status_ok=$(/bin/gzip ${file}; echo $?)
	echo "Status of zipping ${file}: $status_ok" >> ${POST_LOG}
	file=${file}.gz
    fi

    # Upload file:
    [ "${file}" != "${file/_d_}" ] \
	&& target=${LOCALPC_DIR}/daily || target=${LOCALPC_DIR}   
    status_ok=$(scp -p ${file} \
	${LOCALPC_LOGINNAME}@${LOCALPC_PCNAME}:${target} \
	> /dev/null; echo $?)
    echo "Status of uploading ${file} to local PC: $status_ok" >> ${POST_LOG}

    # Warnings and removing file:
    if [ $status_ok -ne 0 ]; then
	error_message="WARNING: Something went wrong with uploading ${file} to local PC."
	echo ${error_message} >> ${POST_LOG}
	echo ${error_message} > ${ERR_MES}
	/bin/mail -s "Uploading problems at BRUTUS" ${MAIL_NOTIFICATION} \
	    < ${ERR_MES}
    else
	/bin/rm -f ${file}
    fi

done

# End of postprocessing and storing output data:
################################################   

echo "------------------------------------------" >> $POST_LOG
date >> ${POST_LOG}
echo "End of postprocessing and storing output data of year ${YYYY}, month ${MM}" >> ${POST_LOG}

exit
