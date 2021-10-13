#!/bin/bash
#

cd ${MODELOUTPUTDIR}

# Ouput file names:
NAME_ECHAM_H=_echam_h_
NAME_ECHAM_C=_echam_c_
NAME_CHEM1_H=_chem1_h_
NAME_CHEM1_C=_chem1_c_
NAME_CHEM2_H=_chem2_h_
NAME_CHEM2_C=_chem2_c_

###########################
# 1. Hybrid (model) levels:
###########################

if [ ${LAFTER_ECHAM_H} = .TRUE. ]; then
    ECHAM_H_FILE=${EXPNO}_${YYYY}${MM}_echam
    ECHAM_H_OUT=${EXPNO}_echam_h_out_${YYYY}${MM}.nc
    if [ ! -f $ECHAM_H_FILE ]; then
        echo "File $MODELOUTPUTDIR/$ECHAM_H_FILE not found"
        exit
    fi
    
    ${BURNSCRIPTDIR}/after ${ECHAM_H_FILE} ${ECHAM_H_OUT} < \
        ${MODELSCRIPTSDIR}/${AFTER_ECHAM_H_NAMELIST}

    #Adapt time and levels for Hiphop/Ferret
    ${BURNSCRIPTDIR}/timesocol2hiphop_dv ${ECHAM_H_OUT}

    # Rename and move to ${AFTERPATH}:
    mv ${ECHAM_H_OUT} ${AFTERPATH}/${EXPNO}${NAME_ECHAM_H}${YYYY}${MM}.nc

fi


# Chemistry:
# chem1 fields:
if [ ${LAFTER_CHEM1_H} = .TRUE. ]; then
    CHEM1_H_FILE=${EXPNO}_${YYYY}${MM}.01_chem1.nc
    CHEM1_H_OUT=${EXPNO}_chem1_h_out_${YYYY}${MM}.nc

    if [ ! -f $CHEM1_H_FILE ]; then
        echo "File $MODELOUTPUTDIR/$CHEM1_H_FILE not found"
        exit
    fi

    cp ${EXPNO}_${YYYY}${MM}.01_chem1.nc ${EXPNO}_chem1_h_out_${YYYY}${MM}.nc
    CHEM1_H_OUT=${EXPNO}_chem1_h_out_${YYYY}${MM}.nc

# Rename and move to ${AFTERPATH}:
    mv ${CHEM1_H_OUT} ${AFTERPATH}/${EXPNO}${NAME_CHEM1_H}${YYYY}${MM}.nc

fi

# chem2 fields:
if [ ${LAFTER_CHEM2_H} = .TRUE. ]; then
    CHEM2_H_FILE=${EXPNO}_${YYYY}${MM}.01_chem2.nc
    CHEM2_H_OUT=${EXPNO}_chem2_h_out_${YYYY}${MM}.nc
    if [ ! -f $CHEM2_H_FILE ] ; then
        echo "File $MODELOUTPUTDIR/$CHEM2_H_FILE not found"
        exit
    fi
    cp ${EXPNO}_${YYYY}${MM}.01_chem2.nc ${EXPNO}_chem2_h_out_${YYYY}${MM}.nc

    # Rename and move to ${AFTERPATH}:
    mv ${CHEM2_H_OUT} ${AFTERPATH}/${EXPNO}${NAME_CHEM2_H}${YYYY}${MM}.nc
fi
    
#####################
# 2. Interpolated to CCMI levels:
#####################

# ECHAM fields:	
if [ ${LAFTER_ECHAM_C} = .TRUE. ]; then
    ECHAM_C_FILE=${EXPNO}_${YYYY}${MM}_echam
    ECHAM_C_OUT=${EXPNO}_echam_c_out_${YYYY}${MM}.nc
    if [ ! -f $ECHAM_C_FILE ]; then
	echo "File $MODELOUTPUTDIR/$ECHAM_C_FILE not found"
	exit
    fi
    ${BURNSCRIPTDIR}/after ${ECHAM_C_FILE} ${ECHAM_C_OUT} < \
	${MODELSCRIPTSDIR}/${AFTER_ECHAM_C_NAMELIST}

    #Adapt time and levels for Hiphop/Ferret
    ${BURNSCRIPTDIR}/timesocol2hiphop_dv ${ECHAM_C_OUT}

    #Rename and move to ${AFTERPATH}:
     mv ${ECHAM_C_OUT} ${AFTERPATH}/${EXPNO}${NAME_ECHAM_C}${YYYY}${MM}.nc
fi   

# Chemistry:
# chem1 fields:
if [ ${LAFTER_CHEM1_C} = .TRUE. ]; then
    CHEM1_C_FILE=${EXPNO}_${YYYY}${MM}.01_chem1.nc
    CHEM1_C_OUT=${EXPNO}_chem1_c_out_${YYYY}${MM}.nc
    if [ ! -f $CHEM1_C_FILE ]; then
	echo "File $MODELOUTPUTDIR/$CHEM1_C_FILE not found"
	exit
    fi
    
    # Extract variables: HARDWIRED
    cdo -f nc ml2plx,10,20,30,50,100,150,200,300,500,700,1000,1500,2000,3000,5000,7000,8000,9000,10000,11500,13000,15000,17000,20000,25000,30000,40000,50000,70000,85000,100000, -selvar,aps,O3,O,NO,NO2,HNO3,NO3,N2O5,HNO4,HO2,ClNO3,ClO,N,OH,Cl,HOCl,N2O,CH4,CO,HCl,H2,H2O2,Cl2,Cl2O2,CH3O2H,CH2O,BrO,BrNO3,BrCl,HBr,HOBr,Br,OD,CFC11,CFC12,CBRF3,CCL4,CH3CCL3,HCFC22,HCFC141b,HCFC142b,H1211,CH3BR,CH3CL,H2402,CHBR3,CH2BR2,O3ONPBL_1,O3ONMBL_2,O3OTRBL_3,O3OSMBL_4,O3OSPBL_5,O3ONPFT_6,O3ONMFT_7,O3OTRFT_8,O3OSMFT_9,O3OSPFT_10,O3ONPLS_11,O3ONMLS_12,O3OTRLS_13,O3OTRMS_14,O3OSMLS_15,O3OSPLS_16,O3ONPUS_17,O3ONMUS_18,O3OTRUS_19,O3OSMUS_20,O3OSPUS_21, ${CHEM1_C_FILE} ${CHEM1_C_OUT}

    # Rename and move to ${AFTERPATH}:
    mv ${CHEM1_C_OUT} ${AFTERPATH}/${EXPNO}${NAME_CHEM1_C}${YYYY}${MM}.nc
fi

# chem2 fields:
if [ ${LAFTER_CHEM2_C} = .TRUE. ]; then
    CHEM2_C_FILE=${EXPNO}_${YYYY}${MM}.01_chem2.nc
    CHEM2_C_OUT=${EXPNO}_chem2_c_out_${YYYY}${MM}.nc
    if [ ! -f $CHEM2_C_FILE ] ; then
	echo "File $MODELOUTPUTDIR/$CHEM2_C_FILE not found"
	exit
    fi

    # Extract variables: HARDWIRED
    cdo -f nc ml2plx,10,20,30,50,100,150,200,300,500,700,1000,1500,2000,3000,5000,7000,8000,9000,10000,11500,13000,15000,17000,20000,25000,30000,40000,50000,70000,85000,100000, -selvar,aps,RXN_1,RXN_2,RXN_3,RXN_4,RXN_5,RXN_6,RXN_7,RXN_8,RXN_9,RXN_10,RXN_12,RXN_13,RXN_14,RXN_15,RXN_17,RXN_18,RXN_19,RXN_20,RXN_21,RXN_22,RXN_23,RXN_24,RXN_25,RXN_26,RXN_27,RXN_28,RXN_29,RXN_30,RXN_31,RXN_32,RXN_33,RXN_11,RXN_34,RXN_16, ${CHEM2_C_FILE} ${CHEM2_C_OUT}

    # Rename and move to ${AFTERPATH}:
    mv ${CHEM2_C_OUT} ${AFTERPATH}/${EXPNO}${NAME_CHEM2_C}${YYYY}${MM}.nc
fi

    #Move grib files to ${AFTERPATH}:
    mv ${EXPNO}_${YYYY}${MM}_echam ${AFTERPATH}
    mv ${EXPNO}_${YYYY}${MM}_chem1 ${AFTERPATH}

exit