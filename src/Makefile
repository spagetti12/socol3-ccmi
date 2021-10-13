# Generated automatically by createMakefiles.pl
 
PROG =	../bin/echam5

SRCS =	auxhyb.f90 call_submodels.f90 cloud.f90 clsst.f90 clveg.f90 \
	collect.f90 conteq.f90 control.f90 cover.f90 cuadjtq.f90 cuasc.f90 \
	cuasct.f90 cubase.f90 cubasmc.f90 cucall.f90 cuddraf.f90 cudlfs.f90 \
	cudtdq.f90 cududv.f90 cuentr.f90 cuentrt.f90 cuflx.f90 cuini.f90 \
	cumastr.f90 cumastrh.f90 cumastrt.f90 datim.f90 dyn.f90 ewd.f90 \
	fftd.f90 ffti.f90 geopot.f90 gpc.f90 hdiff.f90 helmo.f90 inhysi.f90 \
	inictl.f90 inidoc.f90 iniphy.f90 inipost.f90 init_decomposition.f90 \
	initemp.f90 initialize.f90 intaero.f90 ioinitial.f90 iorestart.f90 \
	labrun.f90 lake.f90 legtri.f90 licetemp.f90 ltd.f90 lti.f90 \
	m_solang.f90 master.f90 maxwind.f90 ml_flux.f90 ml_ocean.f90 \
	pgrad.f90 physc.f90 prerad.f90 pres.f90 presf.f90 rad_int.f90 \
	radheat.f90 radiation.f90 readfld.f90 rrtm_gasabs1a_140gp.f90 \
	rrtm_rrtm_140gp.f90 rrtm_rtrn1a_140gp.f90 rrtm_setcoef_140gp.f90 \
	rrtm_taumol_140gp.f90 scan1.f90 scan2.f90 sccd.f90 \
	scctp.f90 setdyn.f90 setgws.f90 setphys.f90 setrad.f90 si1.f90 \
	si2.f90 sicetemp.f90 soiltemp.f90 stepon.f90 subjob.f90 sudif.f90 \
	surf.f90 surftemp.f90 surrta.f90 surrtab.f90 surrtbg2.f90 \
	surrtftr.f90 surrtm.f90 surrtpk.f90 surrtrf.f90 susw.f90 sw.f90 \
	sw1s.f90 swclr.f90 swde.f90 swde_wr.f90 swni.f90 swr.f90 swtt.f90 swtt1.f90 \
	swu.f90 swuvo3.f90 sym1.f90 sym2.f90 tf1.f90 tf2.f90 vdiff.f90 \
	socol_boundary_conditions.f90 socol_chemini.f90 socol_col_o3.f90 \
	socol_dis.f90 socol_families_chem.f90 socol_famcorr_adv.f90 \
	socol_mezon.f90 socol_nrd.f90 socol_sediment.f90 

OBJS := $(SRCS:.f90=.o)

MOBJS =	$(MODULES)/m_alloc_mods.o $(MODULES)/mo_advection.o \
	$(MODULES)/mo_aero_gads.o $(MODULES)/mo_aero_tanre.o \
	$(MODULES)/mo_albedo.o $(MODULES)/mo_base_types.o \
	$(MODULES)/mo_buffer_fft.o $(MODULES)/mo_call_trans.o \
	$(MODULES)/mo_clim.o $(MODULES)/mo_cloud.o $(MODULES)/mo_column.o \
	$(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_convect_tables.o $(MODULES)/mo_couple.o \
	$(MODULES)/mo_cumulus_flux.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_diag_amip2.o $(MODULES)/mo_diag_dynamics.o \
	$(MODULES)/mo_diag_radiation.o $(MODULES)/mo_diag_tendency.o \
	$(MODULES)/mo_diff.o $(MODULES)/mo_doctor.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_fft991.o $(MODULES)/mo_fft992.o $(MODULES)/mo_field.o \
	$(MODULES)/mo_filename.o $(MODULES)/mo_forecast_switches.o \
	$(MODULES)/mo_gaussgrid.o $(MODULES)/mo_geoloc.o $(MODULES)/mo_gl1.o \
	$(MODULES)/mo_global_op.o $(MODULES)/mo_greenhouse_gases.o \
	$(MODULES)/mo_grib.o $(MODULES)/mo_gwspectrum.o $(MODULES)/mo_hdiff.o \
	$(MODULES)/mo_hyb.o $(MODULES)/mo_hydrology.o $(MODULES)/mo_interpo.o \
	$(MODULES)/mo_io.o $(MODULES)/mo_kind.o $(MODULES)/mo_legendre.o \
	$(MODULES)/mo_linked_list.o $(MODULES)/mo_lw_clop.o \
	$(MODULES)/mo_machine.o $(MODULES)/mo_memory_base.o \
	$(MODULES)/mo_memory_f.o $(MODULES)/mo_memory_g1a.o \
	$(MODULES)/mo_memory_g1b.o $(MODULES)/mo_memory_g2a.o \
	$(MODULES)/mo_memory_g2b.o $(MODULES)/mo_memory_g3a.o \
	$(MODULES)/mo_memory_g3b.o $(MODULES)/mo_memory_gl.o \
	$(MODULES)/mo_memory_ls.o $(MODULES)/mo_memory_sp.o \
	$(MODULES)/mo_memory_streams.o $(MODULES)/mo_midatm.o \
	$(MODULES)/mo_mpi.o $(MODULES)/mo_namelist.o $(MODULES)/mo_netcdf.o \
	$(MODULES)/mo_netcdfstream.o $(MODULES)/mo_nmi.o \
	$(MODULES)/mo_nudging.o $(MODULES)/mo_nudging_buffer.o \
	$(MODULES)/mo_nudging_constants.o $(MODULES)/mo_nudging_init.o \
	$(MODULES)/mo_nudging_io.o $(MODULES)/mo_nudging_pattern.o \
	$(MODULES)/mo_nudging_sst.o $(MODULES)/mo_nudging_utils.o \
	$(MODULES)/mo_o3_lwb.o $(MODULES)/mo_o3clim.o $(MODULES)/mo_orbit.o \
	$(MODULES)/mo_param_switches.o $(MODULES)/mo_parameters.o \
	$(MODULES)/mo_parrrtm.o $(MODULES)/mo_physc1.o $(MODULES)/mo_physc2.o \
	$(MODULES)/mo_port_test.o $(MODULES)/mo_post.o \
	$(MODULES)/mo_radiation.o $(MODULES)/mo_real_timer.o \
	$(MODULES)/mo_rrta_140gp.o $(MODULES)/mo_rrtab.o \
	$(MODULES)/mo_rrtbg2.o $(MODULES)/mo_rrtftr.o $(MODULES)/mo_rrtrf.o \
	$(MODULES)/mo_rrtwn.o $(MODULES)/mo_scan_buffer.o \
	$(MODULES)/mo_semi_impl.o $(MODULES)/mo_semi_lagrangian.o \
	$(MODULES)/mo_socol_ionization.o \
	$(MODULES)/mo_so4.o $(MODULES)/mo_solmon.o \
	$(MODULES)/mo_spectral.o $(MODULES)/mo_spitfire.o \
	$(MODULES)/mo_ssodrag.o $(MODULES)/mo_ssortns.o $(MODULES)/mo_sst.o \
	$(MODULES)/mo_sub_echam.o $(MODULES)/mo_sub_nml.o $(MODULES)/mo_sw.o \
	$(MODULES)/mo_sw_clop.o $(MODULES)/mo_test_trans.o \
	$(MODULES)/mo_time_base.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_time_conversion.o $(MODULES)/mo_time_event.o \
	$(MODULES)/mo_time_manager.o $(MODULES)/mo_timer.o \
	$(MODULES)/mo_timestamp.o $(MODULES)/mo_tmp_buffer.o \
	$(MODULES)/mo_tpcore.o $(MODULES)/mo_tracdef.o $(MODULES)/mo_tracer.o \
	$(MODULES)/mo_transpose.o $(MODULES)/mo_tropopause.o \
	$(MODULES)/mo_truncation.o $(MODULES)/mo_upper_sponge.o \
	$(MODULES)/mo_util_string.o $(MODULES)/mo_vegetation.o \
	$(MODULES)/mo_vsop87.o $(MODULES)/mo_socol_chem.o \
	$(MODULES)/mo_socol_constants.o $(MODULES)/mo_socol_co_nox.o \
	$(MODULES)/mo_socol_deposition.o $(MODULES)/mo_socol_dimensions.o \
	$(MODULES)/mo_socol_gcmfields.o $(MODULES)/mo_socol_ghg_ods.o \
	$(MODULES)/mo_socol_grid_calculations.o $(MODULES)/mo_socol_hetero.o \
	$(MODULES)/mo_socol_interpo.o $(MODULES)/mo_socol_lnox_main.o \
	$(MODULES)/mo_socol_namelist.o  $(MODULES)/mo_socol_qbo.o \
	$(MODULES)/mo_socol_readfile.o $(MODULES)/mo_socol_synth_tracers.o \
	$(MODULES)/mo_socol_strataerosols.o $(MODULES)/mo_socol_streams.o \
	$(MODULES)/mo_socol_sun.o $(MODULES)/mo_socol_time_control.o \
	$(MODULES)/mo_socol_tracers.o $(MODULES)/mo_socol_tropoaerosols.o \
	$(MODULES)/mo_socol_ch4.o $(MODULES)/mo_socol_isotope.o \
	$(MODULES)/mo_socol_ch4_streams.o $(MODULES)/mo_socol_no_gcm_condens_pscs.o \
	$(MODULES)/mo_socol_photo.o $(MODULES)/mo_socol_scav.o $(MODULES)/mo_socol_nmvoc.o \
        $(MODULES)/mo_ohclim.o $(MODULES)/mo_lpj_streams.o \
        $(MODULES)/mo_socol_o3orig.o

all: $(PROG)

$(PROG): $(OBJS) $(MOBJS) ../lib/libsupport.a
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(MOBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

auxhyb.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o
call_submodels.o: $(MODULES)/mo_kind.o $(MODULES)/mo_sub_echam.o \
	$(MODULES)/mo_sub_nml.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_tracer.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_socol_gcmfields.o \
	$(MODULES)/mo_socol_grid_calculations.o \
	$(MODULES)/mo_socol_namelist.o \
	$(MODULES)/mo_socol_no_gcm_condens_pscs.o \
	$(MODULES)/mo_socol_time_control.o \
        $(MODULES)/mo_socol_streams.o $(MODULES)/mo_socol_sun.o \
        $(MODULES)/mo_socol_tracers.o $(MODULES)/mo_socol_scav.o \
        $(MODULES)/mo_lpj_streams.o $(MODULES)/mo_socol_o3orig.o \
        $(MODULES)/mo_socol_ch4.o $(MODULES)/mo_socol_synth_tracers.o
cloud.o: $(MODULES)/mo_cloud.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_convect_tables.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_param_switches.o \
	$(MODULES)/mo_time_control.o
clsst.o: $(MODULES)/mo_column.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_geoloc.o $(MODULES)/mo_interpo.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_g3b.o $(MODULES)/mo_physc2.o \
	$(MODULES)/mo_radiation.o $(MODULES)/mo_sst.o
clveg.o: $(MODULES)/mo_clim.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_interpo.o $(MODULES)/mo_memory_g3b.o \
	$(MODULES)/mo_radiation.o
collect.o: $(MODULES)/mo_constants.o $(MODULES)/mo_couple.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_time_control.o
conteq.o: $(MODULES)/mo_control.o $(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o
control.o: $(MODULES)/m_alloc_mods.o $(MODULES)/mo_advection.o \
	$(MODULES)/mo_aero_tanre.o $(MODULES)/mo_clim.o \
	$(MODULES)/mo_column.o $(MODULES)/mo_control.o $(MODULES)/mo_couple.o \
	$(MODULES)/mo_decomposition.o $(MODULES)/mo_diag_amip2.o \
	$(MODULES)/mo_diag_tendency.o $(MODULES)/mo_doctor.o \
	$(MODULES)/mo_exception.o $(MODULES)/mo_fft991.o \
	$(MODULES)/mo_fft992.o $(MODULES)/mo_field.o \
	$(MODULES)/mo_gaussgrid.o $(MODULES)/mo_geoloc.o $(MODULES)/mo_grib.o \
	$(MODULES)/mo_hydrology.o $(MODULES)/mo_io.o $(MODULES)/mo_legendre.o \
	$(MODULES)/mo_machine.o $(MODULES)/mo_memory_f.o \
	$(MODULES)/mo_memory_g3b.o $(MODULES)/mo_memory_streams.o \
	$(MODULES)/mo_mpi.o $(MODULES)/mo_netcdf.o $(MODULES)/mo_nmi.o \
	$(MODULES)/mo_nudging_init.o $(MODULES)/mo_o3clim.o \
	$(MODULES)/mo_radiation.o $(MODULES)/mo_scan_buffer.o \
	$(MODULES)/mo_semi_lagrangian.o $(MODULES)/mo_so4.o \
	$(MODULES)/mo_spitfire.o $(MODULES)/mo_sst.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_timer.o \
	$(MODULES)/mo_tmp_buffer.o $(MODULES)/mo_tpcore.o \
	$(MODULES)/mo_tracer.o $(MODULES)/mo_ohclim.o
cover.o: $(MODULES)/mo_cloud.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_convect_tables.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_param_switches.o
cuadjtq.o: $(MODULES)/mo_constants.o $(MODULES)/mo_convect_tables.o \
	$(MODULES)/mo_kind.o
cuasc.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_cumulus_flux.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_time_control.o
cuasct.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_cumulus_flux.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_time_control.o
cubase.o: $(MODULES)/mo_constants.o $(MODULES)/mo_cumulus_flux.o \
	$(MODULES)/mo_kind.o
cubasmc.o: $(MODULES)/mo_constants.o $(MODULES)/mo_cumulus_flux.o \
	$(MODULES)/mo_kind.o
cucall.o: $(MODULES)/mo_constants.o $(MODULES)/mo_convect_tables.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_param_switches.o \
	$(MODULES)/mo_time_control.o
cuddraf.o: $(MODULES)/mo_constants.o $(MODULES)/mo_cumulus_flux.o \
	$(MODULES)/mo_kind.o
cudlfs.o: $(MODULES)/mo_constants.o $(MODULES)/mo_cumulus_flux.o \
	$(MODULES)/mo_kind.o
cudtdq.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_time_control.o \
	$(MODULES)/mo_tracer.o $(MODULES)/mo_control.o 
cududv.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o
cuentr.o: $(MODULES)/mo_constants.o $(MODULES)/mo_cumulus_flux.o \
	$(MODULES)/mo_kind.o
cuentrt.o: $(MODULES)/mo_constants.o $(MODULES)/mo_cumulus_flux.o \
	$(MODULES)/mo_kind.o
cuflx.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_physc2.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_socol_streams.o
cuini.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o
cumastr.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_convect_tables.o $(MODULES)/mo_cumulus_flux.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_time_control.o \
       $(MODULES)/mo_socol_lnox_main.o
cumastrh.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_convect_tables.o $(MODULES)/mo_cumulus_flux.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_time_control.o
cumastrt.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_convect_tables.o $(MODULES)/mo_cumulus_flux.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_time_control.o
dyn.o: $(MODULES)/mo_advection.o $(MODULES)/mo_column.o \
	$(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_decomposition.o $(MODULES)/mo_diag_tendency.o \
	$(MODULES)/mo_gaussgrid.o $(MODULES)/mo_geoloc.o \
	$(MODULES)/mo_global_op.o $(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_g3a.o $(MODULES)/mo_memory_g3b.o \
	$(MODULES)/mo_memory_gl.o $(MODULES)/mo_scan_buffer.o \
	$(MODULES)/mo_transpose.o
ewd.o: $(MODULES)/mo_buffer_fft.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_gaussgrid.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_truncation.o
fftd.o: $(MODULES)/mo_buffer_fft.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_fft991.o $(MODULES)/mo_fft992.o $(MODULES)/mo_kind.o
ffti.o: $(MODULES)/mo_buffer_fft.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_fft991.o $(MODULES)/mo_fft992.o $(MODULES)/mo_kind.o
geopot.o: $(MODULES)/mo_control.o $(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o
gpc.o: $(MODULES)/mo_control.o $(MODULES)/mo_nudging_sst.o
hdiff.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_decomposition.o $(MODULES)/mo_diag_tendency.o \
	$(MODULES)/mo_diff.o $(MODULES)/mo_hdiff.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_sp.o $(MODULES)/mo_semi_impl.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_truncation.o
helmo.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_decomposition.o $(MODULES)/mo_doctor.o \
	$(MODULES)/mo_exception.o $(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_tmp_buffer.o
inhysi.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o $(MODULES)/mo_semi_impl.o \
	$(MODULES)/mo_time_control.o
inictl.o: $(MODULES)/mo_advection.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_doctor.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_filename.o $(MODULES)/mo_mpi.o $(MODULES)/mo_namelist.o \
	$(MODULES)/mo_parameters.o $(MODULES)/mo_port_test.o \
	$(MODULES)/mo_time_base.o $(MODULES)/mo_time_control.o \
	$(INCLUDE)/runctl.inc
inidoc.o: $(MODULES)/mo_doctor.o $(MODULES)/mo_kind.o $(MODULES)/mo_mpi.o \
	$(MODULES)/mo_netcdf.o
iniphy.o: $(MODULES)/mo_cloud.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_cumulus_flux.o \
	$(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o $(MODULES)/mo_physc2.o \
	$(MODULES)/mo_vegetation.o
inipost.o: $(MODULES)/mo_doctor.o $(MODULES)/mo_memory_base.o \
	$(MODULES)/mo_memory_gl.o $(MODULES)/mo_memory_sp.o \
	$(MODULES)/mo_mpi.o $(MODULES)/mo_namelist.o $(MODULES)/mo_post.o \
	$(INCLUDE)/postctl.inc
init_decomposition.o: $(MODULES)/mo_advection.o $(MODULES)/mo_column.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_doctor.o $(MODULES)/mo_exception.o $(MODULES)/mo_mpi.o \
	$(MODULES)/mo_transpose.o $(MODULES)/mo_util_string.o
initemp.o: $(MODULES)/mo_clim.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_interpo.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_g3b.o $(MODULES)/mo_parameters.o \
	$(MODULES)/mo_physc2.o $(MODULES)/mo_radiation.o $(MODULES)/mo_sst.o \
	$(MODULES)/mo_time_control.o
initialize.o: $(MODULES)/m_alloc_mods.o $(MODULES)/mo_column.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_greenhouse_gases.o \
	$(MODULES)/mo_io.o $(MODULES)/mo_nudging_init.o \
	$(MODULES)/mo_radiation.o $(MODULES)/mo_so4.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_tracer.o
intaero.o: $(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_interpo.o $(MODULES)/mo_memory_g3b.o \
	$(MODULES)/mo_so4.o
ioinitial.o: $(MODULES)/mo_cloud.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_couple.o \
	$(MODULES)/mo_decomposition.o $(MODULES)/mo_doctor.o \
	$(MODULES)/mo_filename.o $(MODULES)/mo_gaussgrid.o \
	$(MODULES)/mo_hyb.o $(MODULES)/mo_io.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_base.o $(MODULES)/mo_memory_g3a.o \
	$(MODULES)/mo_memory_g3b.o $(MODULES)/mo_memory_gl.o \
	$(MODULES)/mo_memory_sp.o $(MODULES)/mo_mpi.o $(MODULES)/mo_tracer.o \
	$(MODULES)/mo_transpose.o $(MODULES)/mo_util_string.o
iorestart.o: $(MODULES)/mo_control.o $(MODULES)/mo_doctor.o \
	$(MODULES)/mo_io.o $(MODULES)/mo_memory_g1a.o \
	$(MODULES)/mo_memory_gl.o $(MODULES)/mo_mpi.o $(MODULES)/mo_timer.o \
	$(MODULES)/mo_tracer.o
labrun.o: $(MODULES)/mo_advection.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_doctor.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_filename.o $(MODULES)/mo_mpi.o \
	$(MODULES)/mo_param_switches.o $(MODULES)/mo_semi_impl.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_tracer.o \
	$(MODULES)/mo_truncation.o
lake.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_time_control.o
legtri.o: $(MODULES)/mo_kind.o
licetemp.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_param_switches.o $(MODULES)/mo_time_control.o
ltd.o: $(MODULES)/mo_decomposition.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_legendre.o $(MODULES)/mo_memory_f.o \
	$(MODULES)/mo_memory_ls.o
lti.o: $(MODULES)/mo_decomposition.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_legendre.o $(MODULES)/mo_memory_f.o \
	$(MODULES)/mo_memory_ls.o
m_solang.o: $(MODULES)/mo_decomposition.o $(MODULES)/mo_gaussgrid.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_physc1.o $(MODULES)/mo_radiation.o \
	$(MODULES)/mo_transpose.o
master.o: $(MODULES)/mo_couple.o $(MODULES)/mo_doctor.o \
	$(MODULES)/mo_exception.o $(MODULES)/mo_kind.o $(MODULES)/mo_mpi.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_util_string.o
maxwind.o: $(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_doctor.o $(MODULES)/mo_global_op.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_semi_impl.o $(MODULES)/mo_time_control.o
ml_flux.o: $(MODULES)/mo_decomposition.o $(MODULES)/mo_interpo.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_memory_g3b.o $(MODULES)/mo_sst.o
ml_ocean.o: $(MODULES)/mo_exception.o $(MODULES)/mo_kind.o
pgrad.o: $(MODULES)/mo_control.o $(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o
physc.o: $(MODULES)/mo_advection.o $(MODULES)/mo_cloud.o \
	$(MODULES)/mo_column.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_diag_amip2.o $(MODULES)/mo_diag_radiation.o \
	$(MODULES)/mo_diag_tendency.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_geoloc.o $(MODULES)/mo_hyb.o $(MODULES)/mo_hydrology.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_memory_g1a.o \
	$(MODULES)/mo_memory_g2a.o $(MODULES)/mo_memory_g3a.o \
	$(MODULES)/mo_memory_g3b.o $(MODULES)/mo_midatm.o $(MODULES)/mo_nmi.o \
	$(MODULES)/mo_param_switches.o $(MODULES)/mo_physc1.o \
	$(MODULES)/mo_physc2.o $(MODULES)/mo_radiation.o \
	$(MODULES)/mo_scan_buffer.o $(MODULES)/mo_spitfire.o \
	$(MODULES)/mo_ssortns.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_timer.o $(MODULES)/mo_tracer.o \
	$(MODULES)/mo_tropopause.o $(MODULES)/mo_socol_namelist.o \
	$(MODULES)/mo_socol_no_gcm_condens_pscs.o $(MODULES)/mo_socol_qbo.o  \
	$(MODULES)/mo_socol_sun.o $(MODULES)/mo_lpj_streams.o
prerad.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_diag_radiation.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_geoloc.o $(MODULES)/mo_kind.o $(MODULES)/mo_o3_lwb.o \
	$(MODULES)/mo_o3clim.o $(MODULES)/mo_orbit.o \
	$(MODULES)/mo_param_switches.o $(MODULES)/mo_physc1.o \
	$(MODULES)/mo_radiation.o $(MODULES)/mo_time_base.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_time_conversion.o \
	m_solang.o
pres.o: $(MODULES)/mo_control.o $(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o
presf.o: $(MODULES)/mo_control.o $(MODULES)/mo_kind.o
rad_int.o: $(MODULES)/mo_aero_gads.o $(MODULES)/mo_aero_tanre.o \
	$(MODULES)/mo_cloud.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_lw_clop.o $(MODULES)/mo_parrrtm.o \
	$(MODULES)/mo_radiation.o $(MODULES)/mo_sw.o $(MODULES)/mo_sw_clop.o \
        $(MODULES)/mo_socol_strataerosols.o
radheat.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_diag_tendency.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_radiation.o $(MODULES)/mo_time_control.o
radiation.o: $(MODULES)/mo_aero_gads.o $(MODULES)/mo_aero_tanre.o \
	$(MODULES)/mo_albedo.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_convect_tables.o $(MODULES)/mo_doctor.o \
	$(MODULES)/mo_exception.o $(MODULES)/mo_geoloc.o \
	$(MODULES)/mo_greenhouse_gases.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_o3_lwb.o $(MODULES)/mo_o3clim.o $(MODULES)/mo_physc1.o \
	$(MODULES)/mo_radiation.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_socol_ghg_ods.o $(MODULES)/mo_socol_tropoaerosols.o
readfld.o: $(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_doctor.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_field.o $(MODULES)/mo_io.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_mpi.o $(MODULES)/mo_test_trans.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_transpose.o
rrtm_gasabs1a_140gp.o: $(MODULES)/mo_kind.o $(MODULES)/mo_parrrtm.o \
	$(MODULES)/mo_rrtab.o
rrtm_rrtm_140gp.o: $(MODULES)/mo_kind.o $(MODULES)/mo_parrrtm.o
rrtm_rtrn1a_140gp.o: $(MODULES)/mo_kind.o $(MODULES)/mo_parrrtm.o \
	$(MODULES)/mo_rrtab.o $(MODULES)/mo_rrtftr.o $(MODULES)/mo_rrtwn.o
rrtm_setcoef_140gp.o: $(MODULES)/mo_kind.o $(MODULES)/mo_parrrtm.o \
	$(MODULES)/mo_rrtrf.o
rrtm_taumol_140gp.o: $(MODULES)/mo_kind.o $(MODULES)/mo_parrrtm.o \
	$(MODULES)/mo_rrta_140gp.o $(MODULES)/mo_rrtbg2.o \
	$(MODULES)/mo_rrtwn.o
scan1.o: $(MODULES)/mo_advection.o $(MODULES)/mo_call_trans.o \
	$(MODULES)/mo_column.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_diag_amip2.o $(MODULES)/mo_diag_dynamics.o \
	$(MODULES)/mo_diag_radiation.o $(MODULES)/mo_diag_tendency.o \
	$(MODULES)/mo_doctor.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_geoloc.o $(MODULES)/mo_global_op.o \
	$(MODULES)/mo_greenhouse_gases.o $(MODULES)/mo_hyb.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_legendre.o \
	$(MODULES)/mo_memory_g1a.o $(MODULES)/mo_memory_g1b.o \
	$(MODULES)/mo_memory_g2a.o $(MODULES)/mo_memory_g2b.o \
	$(MODULES)/mo_memory_g3b.o $(MODULES)/mo_memory_gl.o \
	$(MODULES)/mo_memory_ls.o $(MODULES)/mo_mpi.o \
	$(MODULES)/mo_port_test.o $(MODULES)/mo_radiation.o \
	$(MODULES)/mo_scan_buffer.o $(MODULES)/mo_semi_lagrangian.o \
	$(MODULES)/mo_solmon.o $(MODULES)/mo_spitfire.o \
	$(MODULES)/mo_test_trans.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_timer.o $(MODULES)/mo_tpcore.o $(MODULES)/mo_tracer.o \
	$(MODULES)/mo_socol_namelist.o $(MODULES)/mo_ohclim.o
scan2.o: $(MODULES)/mo_call_trans.o
sccd.o: $(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_diag_tendency.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_sp.o $(MODULES)/mo_tmp_buffer.o
scctp.o: $(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_diag_tendency.o $(MODULES)/mo_hyb.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_memory_sp.o
setdyn.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_doctor.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_fft991.o $(MODULES)/mo_fft992.o \
	$(MODULES)/mo_forecast_switches.o $(MODULES)/mo_gaussgrid.o \
	$(MODULES)/mo_hdiff.o $(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_mpi.o $(MODULES)/mo_namelist.o \
	$(MODULES)/mo_semi_impl.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_tmp_buffer.o $(MODULES)/mo_truncation.o \
	$(MODULES)/mo_upper_sponge.o $(INCLUDE)/dynctl.inc
setgws.o: $(MODULES)/mo_control.o $(MODULES)/mo_doctor.o \
	$(MODULES)/mo_exception.o $(MODULES)/mo_gwspectrum.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_mpi.o $(MODULES)/mo_namelist.o
setphys.o: $(MODULES)/mo_cloud.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_convect_tables.o \
	$(MODULES)/mo_cumulus_flux.o $(MODULES)/mo_doctor.o \
	$(MODULES)/mo_exception.o $(MODULES)/mo_kind.o $(MODULES)/mo_mpi.o \
	$(MODULES)/mo_namelist.o $(MODULES)/mo_param_switches.o \
	$(MODULES)/mo_physc1.o $(INCLUDE)/physctl.inc
setrad.o: $(MODULES)/mo_aero_gads.o $(MODULES)/mo_aero_tanre.o \
	$(MODULES)/mo_albedo.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_doctor.o \
	$(MODULES)/mo_exception.o $(MODULES)/mo_hyb.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_mpi.o $(MODULES)/mo_namelist.o $(MODULES)/mo_o3clim.o \
	$(MODULES)/mo_param_switches.o $(MODULES)/mo_radiation.o \
	$(MODULES)/mo_solmon.o $(MODULES)/mo_time_base.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_socol_tropoaerosols.o \
        $(MODULES)/mo_ohclim.o
si1.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_decomposition.o $(MODULES)/mo_diag_tendency.o \
	$(MODULES)/mo_gaussgrid.o $(MODULES)/mo_geoloc.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_g1a.o $(MODULES)/mo_scan_buffer.o \
	$(MODULES)/mo_semi_impl.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_tracer.o $(MODULES)/mo_transpose.o
si2.o: $(MODULES)/mo_buffer_fft.o $(MODULES)/mo_constants.o \
	$(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_diag_tendency.o $(MODULES)/mo_gaussgrid.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_semi_impl.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_truncation.o
sicetemp.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_param_switches.o \
	$(MODULES)/mo_physc2.o $(MODULES)/mo_time_control.o
socol_boundary_conditions.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_g3b.o $(MODULES)/mo_radiation.o \
	$(MODULES)/mo_socol_constants.o $(MODULES)/mo_socol_co_nox.o \
	$(MODULES)/mo_socol_deposition.o $(MODULES)/mo_socol_ghg_ods.o \
	$(MODULES)/mo_socol_grid_calculations.o \
	$(MODULES)/mo_socol_ionization.o $(MODULES)/mo_socol_streams.o \
	$(MODULES)/mo_socol_namelist.o $(MODULES)/mo_socol_qbo.o \
	$(MODULES)/mo_socol_strataerosols.o $(MODULES)/mo_socol_sun.o \
	$(MODULES)/mo_socol_time_control.o $(MODULES)/mo_socol_tracers.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_socol_ch4.o \
        $(MODULES)/mo_exception.o $(MODULES)/mo_socol_nmvoc.o \
	$(MODULES)/mo_socol_synth_tracers.o \
	$(MODULES)/mo_socol_lnox_main.o
socol_chemini.o: $(MODULES)/mo_exception.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_socol_chem.o $(MODULES)/mo_socol_constants.o \
	$(MODULES)/mo_socol_dimensions.o $(MODULES)/mo_socol_gcmfields.o \
        $(MODULES)/mo_socol_ghg_ods.o $(MODULES)/mo_socol_grid_calculations.o \
	$(MODULES)/mo_socol_hetero.o $(MODULES)/mo_socol_namelist.o \
	$(MODULES)/mo_socol_strataerosols.o $(MODULES)/mo_socol_streams.o \
	$(MODULES)/mo_socol_sun.o $(MODULES)/mo_socol_time_control.o \
        $(MODULES)/mo_socol_tracers.o $(MODULES)/mo_socol_ch4_streams.o \
	$(MODULES)/mo_socol_photo.o $(MODULES)/mo_memory_g3b.o \
	$(MODULES)/mo_geoloc.o $(MODULES)/mo_memory_g1a.o \
        $(MODULES)/mo_socol_nmvoc.o $(MODULES)/mo_socol_o3orig.f90
socol_col_o3.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_socol_constants.o $(MODULES)/mo_socol_dimensions.o \
	$(MODULES)/mo_socol_gcmfields.o $(MODULES)/mo_socol_streams.o \
	$(MODULES)/mo_socol_tracers.o
socol_dis.o: $(MODULES)/mo_constants.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_geoloc.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_socol_constants.o $(MODULES)/mo_socol_dimensions.o \
	$(MODULES)/mo_socol_grid_calculations.o \
	$(MODULES)/mo_socol_namelist.o $(MODULES)/mo_socol_sun.o
socol_families_chem.o: $(MODULES)/mo_kind.o $(MODULES)/mo_socol_streams.o \
	$(MODULES)/mo_socol_tracers.o $(MODULES)/mo_time_control.o
socol_famcorr_adv.o: $(MODULES)/mo_decomposition.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_g1a.o $(MODULES)/mo_scan_buffer.o \
	$(MODULES)/mo_socol_constants.o $(MODULES)/mo_socol_namelist.o \
	$(MODULES)/mo_socol_streams.o $(MODULES)/mo_socol_tracers.o \
	$(MODULES)/mo_time_control.o
socol_mezon.o: $(MODULES)/mo_constants.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_socol_constants.o \
	$(MODULES)/mo_socol_dimensions.o $(MODULES)/mo_socol_gcmfields.o \
	$(MODULES)/mo_socol_grid_calculations.o $(MODULES)/mo_socol_namelist.o \
	$(MODULES)/mo_socol_ionization.o $(MODULES)/mo_socol_ch4.o \
	$(MODULES)/mo_socol_time_control.o $(MODULES)/mo_socol_tracers.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_socol_o3orig.o \
	$(MODULES)/mo_socol_synth_tracers.o
socol_nrd.o: $(MODULES)/mo_kind.o $(MODULES)/mo_socol_chem.o \
	$(MODULES)/mo_socol_dimensions.o
socol_sediment.o: $(MODULES)/mo_constants.o $(MODULES)/mo_geoloc.o \
 	$(MODULES)/mo_kind.o $(MODULES)/mo_socol_constants.o \
	$(MODULES)/mo_socol_dimensions.o $(MODULES)/mo_socol_gcmfields.o \
        $(MODULES)/mo_socol_grid_calculations.o $(MODULES)/mo_socol_namelist.o \
	$(MODULES)/mo_socol_streams.o $(MODULES)/mo_socol_time_control.o \
     	$(MODULES)/mo_socol_tracers.o $(MODULES)/mo_socol_synth_tracers.o
soiltemp.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_parameters.o $(MODULES)/mo_physc2.o \
	$(MODULES)/mo_time_control.o
stepon.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_couple.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_diag_tendency.o $(MODULES)/mo_doctor.o \
	$(MODULES)/mo_exception.o $(MODULES)/mo_grib.o $(MODULES)/mo_hdiff.o \
	$(MODULES)/mo_hyb.o $(MODULES)/mo_hydrology.o $(MODULES)/mo_io.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_memory_f.o \
	$(MODULES)/mo_memory_sp.o $(MODULES)/mo_mpi.o $(MODULES)/mo_nmi.o \
	$(MODULES)/mo_nudging.o $(MODULES)/mo_nudging_buffer.o \
	$(MODULES)/mo_nudging_init.o $(MODULES)/mo_nudging_sst.o \
	$(MODULES)/mo_nudging_utils.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_timer.o $(MODULES)/mo_upper_sponge.o \
	$(MODULES)/mo_socol_time_control.o $(MODULES)/mo_socol_lnox_main.o \
	$(MODULES)/mo_socol_ch4.o $(MODULES)/mo_socol_o3orig.o
subjob.o: $(MODULES)/mo_control.o $(MODULES)/mo_exception.o \
	$(MODULES)/mo_filename.o $(MODULES)/mo_memory_base.o \
	$(MODULES)/mo_mpi.o $(MODULES)/mo_time_control.o
sudif.o: $(MODULES)/mo_control.o $(MODULES)/mo_diff.o \
	$(MODULES)/mo_exception.o
surf.o: $(MODULES)/mo_constants.o $(MODULES)/mo_control.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_param_switches.o \
	$(MODULES)/mo_physc2.o $(MODULES)/mo_semi_impl.o \
	$(MODULES)/mo_time_control.o $(MODULES)/mo_vegetation.o
surftemp.o: $(MODULES)/mo_kind.o
surrta.o: $(MODULES)/mo_mpi.o $(MODULES)/mo_rrta_140gp.o
surrtab.o: $(MODULES)/mo_kind.o $(MODULES)/mo_rrtab.o
surrtbg2.o: $(MODULES)/mo_kind.o $(MODULES)/mo_rrtbg2.o
surrtftr.o: $(MODULES)/mo_kind.o $(MODULES)/mo_rrtftr.o
surrtpk.o: $(MODULES)/mo_kind.o $(MODULES)/mo_rrtwn.o
surrtrf.o: $(MODULES)/mo_kind.o $(MODULES)/mo_rrtrf.o
susw.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o $(MODULES)/mo_sw.o
sw.o: $(MODULES)/mo_constants.o $(MODULES)/mo_kind.o $(MODULES)/mo_sw.o \
	$(MODULES)/mo_time_control.o
sw1s.o: $(MODULES)/mo_kind.o $(MODULES)/mo_sw.o
swclr.o: $(MODULES)/mo_aero_gads.o $(MODULES)/mo_aero_tanre.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_radiation.o $(MODULES)/mo_sw.o \
	$(MODULES)/mo_socol_strataerosols.o
swde.o: $(MODULES)/mo_kind.o
swde_wr.o: $(MODULES)/mo_kind.o
swni.o: $(MODULES)/mo_kind.o $(MODULES)/mo_sw.o
swr.o: $(MODULES)/mo_kind.o $(MODULES)/mo_sw.o
swtt.o: $(MODULES)/mo_kind.o $(MODULES)/mo_sw.o
swtt1.o: $(MODULES)/mo_kind.o $(MODULES)/mo_sw.o
swu.o: $(MODULES)/mo_kind.o $(MODULES)/mo_sw.o
swuvo3.o: $(MODULES)/mo_kind.o $(MODULES)/mo_sw.o
sym1.o: $(MODULES)/mo_buffer_fft.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_memory_f.o
sym2.o: $(MODULES)/mo_buffer_fft.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_kind.o $(MODULES)/mo_memory_f.o
tf1.o: $(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_diag_tendency.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_g1a.o $(MODULES)/mo_memory_g1b.o \
	$(MODULES)/mo_memory_g2a.o $(MODULES)/mo_memory_g2b.o \
	$(MODULES)/mo_memory_gl.o $(MODULES)/mo_scan_buffer.o \
	$(MODULES)/mo_semi_impl.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_tracer.o
tf2.o: $(MODULES)/mo_control.o $(MODULES)/mo_decomposition.o \
	$(MODULES)/mo_diag_tendency.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_memory_g1a.o $(MODULES)/mo_memory_g2a.o \
	$(MODULES)/mo_memory_gl.o $(MODULES)/mo_scan_buffer.o \
	$(MODULES)/mo_semi_impl.o $(MODULES)/mo_time_control.o \
	$(MODULES)/mo_tracer.o
vdiff.o: $(MODULES)/mo_constants.o $(MODULES)/mo_convect_tables.o \
	$(MODULES)/mo_exception.o $(MODULES)/mo_geoloc.o $(MODULES)/mo_kind.o \
	$(MODULES)/mo_param_switches.o $(MODULES)/mo_physc2.o \
	$(MODULES)/mo_radiation.o $(MODULES)/mo_semi_impl.o \
	$(MODULES)/mo_time_control.o \
	$(MODULES)/mo_tracer.o $(MODULES)/mo_vegetation.o

ifeq ($(strip $(ARCH)), SX)
lti.o: lti.f90
	$(F90) $(F90FLAGS) -Npi -c lti.f90
endif


ifeq ($(strip $(ARCH)), ES)
lti.o: lti.f90
	$(F90) $(F90FLAGS) -Npi -c lti.f90
endif

