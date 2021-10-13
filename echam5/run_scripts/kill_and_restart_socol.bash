#!/bin/bash

# Kills and restarts BRUTUS run jobs.

# Martin Schraner, ETH Zuerich, February 2010

###############################################################################


echo "==========================================" >> ${RUN_LOG}
date >> ${RUN_LOG}
echo "****************************************" >> ${RUN_LOG}
echo "* Killing and restarting run job(s)... *" >> ${RUN_LOG}
echo "****************************************" >> ${RUN_LOG}


# 1. Determine run job number(s) to be killed (after re-starting):

JOBNO_ALL=$(grep 'is submitted to queue' ${RUN_JOB_NUMBERS_LOG} | cut -f2 "-d ")
JOBNO_ALL=${JOBNO_ALL//</}
JOBNO_ALL=${JOBNO_ALL//>/}

echo "Job(s) to be killed: ${JOBNO_ALL}" >> ${RUN_LOG}


# 2. Write restart script:

# File name of new restart script:
RESTARTFILE=restart_run_${EXPNO}_${YYYY}${MM}.bash

echo "------------------------------------------" >> ${RUN_LOG}
date >> ${RUN_LOG}
echo "Creating new run script $RESTARTFILE..." >> ${RUN_LOG}

# Path and file names of available (re-)start scripts (newest first):
RESTARTFILE_OLD=$(/bin/ls -t $MODELSCRIPTSDIR/*start*$EXPNO*.bash 2> /dev/null)

# Error message if no script was found:
if [ $? -ne 0 ]; then
    error_message="ERROR: No (re-)start scripts found in $MODELSCRIPTSDIR. Could not restart BRUTUS run job(s)."
    echo ${error_message} >> ${RUN_LOG}
    echo ${error_message} > ${ERR_MES}
    kill_job ${YYYY} ${MM}
fi

# Path and file names of newest (re-)start script:
RESTARTFILE_OLD=${RESTARTFILE_OLD%%.bash*}.bash

# Copy script:
/bin/cp -p $RESTARTFILE_OLD $RUNSCRIPTDIR/$RESTARTFILE

# Replace RESTARTYEAR and RESTARTMONTH by current year / month:
[ -f out.txt ] && /bin/rm -f out.txt
sed -e "s/^RESTARTYEAR=..../RESTARTYEAR=${YYYY}/" $RUNSCRIPTDIR/$RESTARTFILE \
    > out.txt
/bin/mv -f out.txt $RUNSCRIPTDIR/$RESTARTFILE
sed -e "s/^RESTARTMONTH=../RESTARTMONTH=${MM}/" $RUNSCRIPTDIR/$RESTARTFILE \
    > out.txt
/bin/mv -f out.txt $RUNSCRIPTDIR/$RESTARTFILE

# Make file executable:
chmod +x $RUNSCRIPTDIR/$RESTARTFILE



# 3. Restart BRUTUS:

echo "------------------------------------------" >> ${RUN_LOG}
date >> ${RUN_LOG}
echo "Starting $RUNSCRIPTDIR/$RESTARTFILE..." >> ${RUN_LOG}

cd $RUNSCRIPTDIR
./$RESTARTFILE



# 4. Kill old jobs:

echo "------------------------------------------" >> ${RUN_LOG}
date >> ${RUN_LOG}
echo "Killing old job(s)..." >> ${RUN_LOG}

for JOBNO in $JOBNO_ALL; do
    bkill $JOBNO >> ${RUN_LOG} 2> /dev/null
done



# Do not continue old script:

exit


