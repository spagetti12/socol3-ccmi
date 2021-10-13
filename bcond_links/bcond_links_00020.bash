#!/bin/bash
#
# Contains the links to the boundary condtions used in SOCOL.


# Martin Schraner, ETH Zuerich, February 2010

###############################################################################


# Greenhouse gases and ODSs:
############################
rm -f co2_y ch4_y n2o_y odsmem_y odscll_y odscls_y odsbr_y co2_trac

ln -s  ${INICCMI}/ccmi_co2_RCP6_socol_1950-2100.txt         co2_y
ln -s  ${INICCMI}/ccmi_ch4_RCP6_socol_1950-2100.txt         ch4_y
ln -s  ${INICCMI}/ccmi_n2o_RCP6_socol_1950-2100.txt         n2o_y

#CCMI
ln -s  ${INISOCOL}/../GHG_ODS/ccmi_odsmem_socol_1950-2100.txt      odsmem_y
ln -s  ${INISOCOL}/../GHG_ODS/ccmi_odscls_socol_1950-2100.txt      odscls_y
ln -s  ${INISOCOL}/../GHG_ODS/ccmi_odscll_socol_1950-2100.txt      odscll_y
ln -s  ${INISOCOL}/../GHG_ODS/ccmi_odsbr_socol_1950-2100.txt       odsbr_y

# Solar irradiation, Schumann-Runge bands/Lyman alpha line heating 
# parameterisation, lookup-tables for photolysis rates:
#######################################################
#rm -f sun_irrad_y sun_par_y photolysis* 
rm -f sun_irrad_y sun_par_y photolysis???? photolysis_mean photolysis_delta_e_corr

#CCMI
#ln -s ${INICCMI}/../../solar_data/sun_irrad_C1_1950_2011.txt       sun_irrad_y
#ln -s ${INICCMI}/../../solar_data/sun_par_C1_1950-2011.txt           sun_par_y

ln -s ${INISOCOL}/../CCMI_solar_data/irrad_echam5/sun_irrad_C1_1950_2011.txt       sun_irrad_y
ln -s ${INISOCOL}/../CCMI_solar_data/extra_heating/sun_par_C1_1950-2011.txt           sun_par_y

for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
    ln -s ${INISOCOL}/../CCMI_solar_data/photolysis/photo_C1_${YYYY}.nc photolysis${YYYY}
done

ln -s ${INISOCOL}/../sun/photolysis_socol_mean_1977_1998.nc photolysis_mean # mean 1977-1998 (two solar cycles)

# Delta-E photolysis corrections:
ln -s ${INISOCOL}/../sun/L${LEV}_photolysis_delta_e_corr_socol.nc photolysis_delta_e_corr

#CCMI
# Ionization data
#######################################################
rm -f geomag_data_1600_2100 ????_ionrate-nitrate.txt monthly_mean_ap gcr_data

ln -s ${INISOCOL}/../ionization/geomag_data_1600_2100.txt geomag_data_1600_2100
ln -s ${INISOCOL}/../ionization/gcr_data_1600_2100_constant.nc gcr_data
ln -s ${INISOCOL}/../CCMI_solar_data/ap_index/ap_annual_REF-C1.txt  monthly_mean_ap

echo metka
for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
    ln -s ${INISOCOL}/../CCMI_solar_data/ionrate/${YYYY}_ionrate_spe_REF_C1.txt  ${YYYY}_ionrate-nitrate.txt
   echo ${YYYY}_ionrate-nitrate.txt
done


# Stratospheric aerosols:
#########################
rm -f strataer*

for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
    ln -s ${INISOCOL}/strat_aerosols/CCMI_SOCOL3_4l_${YYYY}.nc strataer${YYYY}
done

#background aerosol climatology (1995-2002)
ln -s ${INISOCOL}/strat_aerosols/T${RES}L${LEV}_sad_nd_ext_omega_g_socol_mean_1995_2002.nc strataerbg

# Tropospheric aerosol climatology:
###################################
#[ $LSOCOL = .TRUE. ] && IAERO=10 || IAERO=2   # (2: Tanre, 10: SOCOL) 
rm -f tropoaer*

#CCMI
for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
    ln -s ${INISOCOL}/tropo_aerosols/T${RES}L${LEV}_GISS_AEROSOLS_${YYYY}.nc tropoaer${YYYY}
done

# aerosol climatology
ln -s ${INISOCOL}/tropo_aerosols/T${RES}L${LEV}_tropoaero_socol_clim.nc tropoaerclim

# CO and NOx emissions:
#######################
rm -f co_nox_emiss_surf_aircr???? nox_lightning lightn_fac

for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
    ln -s ${INISOCOL}/co_nox/T${RES}L${LEV}_co_nox_emiss_surf_aircr_socol_${YYYY}_ccmi.nc co_nox_emiss_surf_aircr${YYYY}
done

ln -s  ${INISOCOL}/co_nox/T${RES}L${LEV}_nox_lightning_socol.nc nox_lightning

ln -s  ${INICCMI}/../tropchem/CCMI_CO_tracer_T${RES}.nc CO_syn_emiss

#CCMI
ln -s  ${INISOCOL}/co_nox/T${RES}_fac_lightning_socol.nc lightn_fac

#  Input files for tropospheric chemistry:
##########################################
rm -f c5h8_emissions ch2o_emissions hcooh_emissions ch3cooh_emissions  #photolysis_rates_messy

# CH3COOH ####################
#For CH3COOH use 1960 emissions for 1950-1959, and 2008 emissions for 2009-2011. LR 04/07/13.
for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
    ln -s ${INISOCOL}/tropchem/ACCMIP_emissions_CH3COOH_T${RES}_${YYYY}.nc ch3cooh_emissions${YYYY}
done
    
# Formaldehyde and isoprene ###
#Use RCP 6.0 emissions from 2001-2010. LR 02/07/13.
for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
  if(($YYYY<2001)); then
     ln -s ${INISOCOL}/tropchem/ACCMIP_emissions_Isoprene_T42_${YYYY}.nc c5h8_emissions${YYYY}
     ln -s ${INISOCOL}/tropchem/ACCMIP_emissions_Formaldehyde_T42_${YYYY}.nc ch2o_emissions${YYYY}  
  else
    ln -s ${INISOCOL}/tropchem/RCP60/ACCMIP_emissions_RCP60_Isoprene_T42_${YYYY}.nc c5h8_emissions${YYYY}
    ln -s ${INISOCOL}/tropchem/RCP60/ACCMIP_emissions_RCP60_Formaldehyde_T42_${YYYY}.nc ch2o_emissions${YYYY}
  fi
done

# Others #########################
ln -s ${INISOCOL}/tropchem/IPCC_emissions_HCOOH_2000_T${RES}.nc hcooh_emissions
ln -s ${INISOCOL}/tropchem/Photolysis_tropchem_T${RES}.nc photolysis_rates_messy

# CH4 emissions:
#######################

rm -f ch4_emiss_surf_anthrop???? ch4_emiss_surf_bb???? ch4_emiss_surf_wet????

for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
    ln -s  ${INISOCOL}/ch4/EDGAR_v41_CH4_anthrop_T${RES}_${YYYY}.nc  ch4_emiss_surf_anthrop${YYYY}
    ln -s  ${INISOCOL}/ch4/GFED3.1_CH4_BB_T${RES}_${YYYY}.nc  ch4_emiss_surf_bb${YYYY}
    ln -s  ${INISOCOL}/ch4/MAIOLICA_CH4_wetlands_T${RES}_${YYYY}.nc  ch4_emiss_surf_wet${YYYY}
done

# Wetland parameterisation:

ln -s ${INISOCOL}/misc/global_wetland_frac_T${RES}.nc unit.40
ln -s ${INISOCOL}/misc/HR_T${RES}.nc unit.41
if  [ ${RES} -eq 31 ]; then
  ln -s ${INISOCOL}/misc/00001_Tsurf_annualmean_1981-88.nc unit.43
else
  ln -s ${INISOCOL}/misc/30008_1995_tsurf_annualmean.nc unit.43
fi

# CMDL CH4 climatology

ln -s ${INISOCOL}/ch4/CH4_clim_T${RES}.nc unit.44

# ECHAM5

ln -s  ${INISOCOL}/T${RES}_O3clim2.nc    unit.21
ln -s  ${INISOCOL}/misc/OH_ACCENT_T${RES}.nc unit.22
ln -s  ${INISOCOL}/T${RES}_VLTCLIM.nc    unit.90
ln -s  ${INISOCOL}/T${RES}_VGRATCLIM.nc  unit.91
ln -s  ${INISOCOL}/T${RES}_TSLCLIM2.nc   unit.92
ln -s  ${INISOCOL}/surrta_data           rrtadata

# NMVOC emissions:
########################
rm -f nmvoc_emiss_biogenic???? nmvoc_emiss_bb???? nmvoc_emiss_anthrop????

# Biomass burning and anthropogenic NMHCs  ###
#Use RCP 6.0 emissions from 2001-2010. LR 02/07/13.
for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
  if(($YYYY<2001)); then
     ln -s ${INISOCOL}/nmhc/ACCMIP_emissions_anthropogenic_NMHC_T42_${YYYY}.nc nmvoc_emiss_anthrop${YYYY}
     ln -s ${INISOCOL}/nmhc/ACCMIP_emissions_biomassburning_NMHC_T42_${YYYY}.nc  nmvoc_emiss_bb${YYYY}
  else
    ln -s ${INISOCOL}/nmhc/RCP60/ACCMIP_emissions_anthropogenic_RCP60_NMHC_T42_${YYYY}.nc nmvoc_emiss_anthrop${YYYY}
    ln -s ${INISOCOL}/nmhc/RCP60/ACCMIP_emissions_biomassburning_RCP60_NMHC_T42_${YYYY}.nc nmvoc_emiss_bb${YYYY}
  fi
done

# Biogenic emissions: Only have up to 2000. Between 2000-2010 NMHC biogenic emissions from the MEGAN model are constant, so use year 2000 emissions. LR 04/07/13.
for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
    ln -s ${INISOCOL}/nmhc/ACCMIP_emissions_biogenic_NMHC_T42_${YYYY}.nc  nmvoc_emiss_biogenic${YYYY}
done

# QBO:
######
rm -f qbo????

for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
    ln -s ${INISOCOL}/../QBO/QBO_socol_${YYYY}.nc qbo${YYYY}
done


# SST/SIC:
##########

## for Hadley (variable) sst and ice (LAMIP=T) use:
for ((YYYY=${Y_SM1}; YYYY<=${Y_EP1}; YYYY++)); do
    ln -s  ${INIHADLEY}/T${RES}_HadISST_sst_${YYYY}.nc sst${YYYY}
    ln -s  ${INIHADLEY}/T${RES}_HadISST_sic_${YYYY}.nc ice${YYYY}
done

# latitude bands for ozone origin diagnostics:
################################################################
ln -s ~/socolv3_ccmi/O3orig/regions4ozone_T42_7regions.nc unit.45
