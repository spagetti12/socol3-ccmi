#!/bin/bash

# Checks remaining runtime.
# If time is (probably) too short to finish the job within 
# $BRUTUS_RUNNINGTIME_MODEL the job is killed and re-started.

# Martin Schraner, ETH Zuerich, February 2010

###############################################################################


# 1. Calculate run time used up to now:

# Actual time:
HH_NOW=$(date +%H)
MM_NOW=$(date +%M)

# Remove leading zero if available:
HOUR_NOW=${HH_NOW#0}
MIN_NOW=${MM_NOW#0}
HOUR_START_RUN=${HH_START_RUN#0}
MIN_START_RUN=${MM_START_RUN#0}

HOUR_USED=$((HOUR_NOW-HOUR_START_RUN))
[ $HOUR_USED -lt 0 ] && HOUR_USED=$((HOUR_USED+24))

MIN_USED=$((MIN_NOW-MIN_START_RUN))
if [ $MIN_USED -lt 0 ]; then
    MIN_USED=$((MIN_USED+60))
    HOUR_USED=$((HOUR_USED-1))
fi

# Run time used up to now [min]:
MIN_USED=$((MIN_USED+HOUR_USED*60))



# 2. Calculate mean run time per month:

# Number of months calculated since run script has been started:
M_CALC=$((M-M_S+12*(YYYY-Y_S)))

# Mean run time per month [min]:
[ $M_CALC -ne 0 ] && MIN_PER_M=$(echo "scale =1; $MIN_USED/$M_CALC" | bc) \
    || MIN_PER_M=0



# 3. Calculate remaining runtime:

HH_RUNTIME=$(echo ${BRUTUS_RUNTIME_MODEL} | cut -c1-2)
MM_RUNTIME=$(echo ${BRUTUS_RUNTIME_MODEL} | cut -c4-5)

# Remove leading zero if available:
HOUR_RUNTIME=${HH_RUNTIME#0}
MIN_RUNTIME=${MM_RUNTIME#0}

# Total runtime in minutes:
MIN_RUNTIME=$((MIN_RUNTIME+HOUR_RUNTIME*60))

# Remaining runtime [min]:
MIN_REMAIN=$((MIN_RUNTIME-MIN_USED))



# 4. Determine critical remaining runtime:

MIN_REMAIN_CRIT=$(echo "scale=2; $BRUTUS_RUNTIME_KILL_FACT*$MIN_PER_M" | bc)
COND=$(echo "$MIN_REMAIN >= $MIN_REMAIN_CRIT" | bc)
[ $COND -eq 1 ] && runtime_ok=0 || runtime_ok=1



# 5. Write log message:

echo "------------------------------------------" >> ${RUN_LOG}
date >> ${RUN_LOG}
echo "Runtime used up to now [min]: ${MIN_USED}" >> ${RUN_LOG}
echo "Number of months finished within current run job: ${M_CALC}" \
    >> ${RUN_LOG}
echo "Mean runtime per month [min]: ${MIN_PER_M}" >> ${RUN_LOG}
echo "Remaining runtime [min]: ${MIN_REMAIN}" >> ${RUN_LOG}
echo "Status of runtime_ok: $runtime_ok" >> ${RUN_LOG}







