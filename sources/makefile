OPSYS = g77
PROGRAM = mceep
DIR_PREFIX = /home/saw/MCEEP/JJ/mceep
#
ifeq ($(OPSYS),hp)
  F77 = f77
  FFLAGS = +E6 +ppu -c -I$(mceep_cmn)
  LINKFLAGS = +E6 +ppu -o
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lmathlib
endif
#
ifeq ($(OPSYS),ultrix)
  F77 = f77
  FFLAGS = -vms -c -I$(mceep_cmn)
  LINKFLAGS = -vms -o
  CERN_ROOT = /usr/local/cern/94b
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lmathlib
endif
#
ifeq ($(OPSYS),osf1)
  F77 = f77
  FFLAGS = -vms -c -I$(mceep_cmn)
  LINKFLAGS = -vms -o
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lmathlib
endif
#
ifeq ($(OPSYS),sun)
  F77 = f77
  FFLAGS = -ext_names=underscores -c -I$(mceep_cmn)
  LINKFLAGS = -o
  LIBS = -L$(CERN_ROOT)/lib \
         -lmathlib
endif
#
ifeq ($(OPSYS),absoft)
  F77 = f77
  FFLAGS = -c -s -N22 -N90 -B108 -f -I$(mceep_cmn)
  LINKFLAGS = -X -defsym -X MAIN__=main -o
  LIBS = -L$(CERN_ROOT)/lib \
         -lmathlib
endif
#
ifeq ($(OPSYS),g77)
  F77 = g77
  FFLAGS = -c -I$(mceep_cmn)
  LINKFLAGS = -o
  LIBS = -L$(CERN_ROOT)/lib -lmathlib
endif
#
ifeq ($(OPSYS),sgi)
  F77 = f77
  FFLAGS = -vms_cc -c -I$(mceep_cmn)
  LINKFLAGS = -vms_cc -o
  CERN_ROOT = /cern/pro
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lmathlib
endif
#
#
SRC_DIR1 = $(DIR_PREFIX)/sources
OBJ_DIR1 = $(DIR_PREFIX)/sources
SRC_DIR2 = $(DIR_PREFIX)/gross
OBJ_DIR2 = $(DIR_PREFIX)/gross
SRC_DIR3 = $(DIR_PREFIX)/cosy
OBJ_DIR3 = $(DIR_PREFIX)/cosy
fort_objs = \
          $(OBJ_DIR1)/arenhoevel.o \
          $(OBJ_DIR1)/borie_drechsel.o \
          $(OBJ_DIR1)/bremsgen.o \
          $(OBJ_DIR1)/c12_elastic.o \
          $(OBJ_DIR1)/calc_b.o \
          $(OBJ_DIR1)/compute_cosy_spin.o \
          $(OBJ_DIR1)/compute_transp_to_cosy.o \
          $(OBJ_DIR1)/decay.o \
          $(OBJ_DIR1)/denscorr.o \
          $(OBJ_DIR1)/deut_elastic.o \
          $(OBJ_DIR1)/deut_sabjes.o \
          $(OBJ_DIR1)/electro_prod.o \
          $(OBJ_DIR1)/eloss_e.o \
          $(OBJ_DIR1)/eloss_proton.o \
          $(OBJ_DIR1)/elossw.o \
          $(OBJ_DIR1)/epc.o \
          $(OBJ_DIR1)/ext_bremss.o \
          $(OBJ_DIR1)/fint.o \
          $(OBJ_DIR1)/formfact.o \
          $(OBJ_DIR1)/get_vect.o \
          $(OBJ_DIR1)/glande.o \
          $(OBJ_DIR1)/glandg.o \
          $(OBJ_DIR1)/glands.o \
          $(OBJ_DIR1)/gvaviv.o \
          $(OBJ_DIR1)/he3_elastic.o \
          $(OBJ_DIR1)/he4_elastic.o \
          $(OBJ_DIR1)/hrs.o \
          $(OBJ_DIR1)/hrs_inv.o \
          $(OBJ_DIR1)/hyd_elastic.o \
          $(OBJ_DIR1)/in_aperture.o \
          $(OBJ_DIR1)/index_accept_fcn.o \
          $(OBJ_DIR1)/input.o \
          $(OBJ_DIR1)/interp.o \
          $(OBJ_DIR1)/kinem_elas_ef.o \
          $(OBJ_DIR1)/kinem_elas_pf.o \
          $(OBJ_DIR1)/kinem.o \
          $(OBJ_DIR1)/kinfac_2hgross.o \
          $(OBJ_DIR1)/kinfac_pwvo.o \
          $(OBJ_DIR1)/Laget_Xsec.o \
          $(OBJ_DIR1)/lomon_ff.o \
          $(OBJ_DIR1)/mad12dfwd_rev.o \
          $(OBJ_DIR1)/mad12.o \
          $(OBJ_DIR1)/mad12_inv.o \
          $(OBJ_DIR1)/mad12nq.o \
          $(OBJ_DIR1)/mad12nq_inv.o \
          $(OBJ_DIR1)/mad35.o \
          $(OBJ_DIR1)/mad35_funcs.o \
          $(OBJ_DIR1)/mad.o \
          $(OBJ_DIR1)/mad_inv.o \
          $(OBJ_DIR1)/mad_nq12_f.o \
          $(OBJ_DIR1)/matrix.o \
          $(OBJ_DIR1)/mceep.o \
          $(OBJ_DIR1)/mceep_misc.o \
          $(OBJ_DIR1)/monte_inv.o \
          $(OBJ_DIR1)/monte_trans.o \
          $(OBJ_DIR1)/mo_tsai.o \
          $(OBJ_DIR1)/mscatt.o \
          $(OBJ_DIR1)/peaking.o \
          $(OBJ_DIR1)/peek.o \
          $(OBJ_DIR1)/peepi.o \
          $(OBJ_DIR1)/physics.o \
          $(OBJ_DIR1)/plot.o \
          $(OBJ_DIR1)/ps_beam.o \
          $(OBJ_DIR1)/pwia_vo.o \
          $(OBJ_DIR1)/qfsv.o \
          $(OBJ_DIR1)/random.o \
          $(OBJ_DIR1)/ranecu.o \
          $(OBJ_DIR1)/ranges.o \
          $(OBJ_DIR1)/read_cosy_input_file.o \
          $(OBJ_DIR1)/r-function.o \
          $(OBJ_DIR1)/rotate_pol.o \
          $(OBJ_DIR1)/schwinger.o \
          $(OBJ_DIR1)/seed_getnsav.o \
          $(OBJ_DIR1)/sep.o \
          $(OBJ_DIR1)/sep_inv.o \
          $(OBJ_DIR1)/septum_inv.o \
          $(OBJ_DIR1)/septum_trans.o \
          $(OBJ_DIR1)/setup_pion.o \
          $(OBJ_DIR1)/sigmacc1.o \
          $(OBJ_DIR1)/sigmacc2.o \
          $(OBJ_DIR1)/sig_mscat.o \
          $(OBJ_DIR1)/smott.o \
          $(OBJ_DIR1)/spectral.o \
          $(OBJ_DIR1)/spectrometer.o \
          $(OBJ_DIR1)/spin_1st_ord.o \
          $(OBJ_DIR1)/spin_dat.o \
          $(OBJ_DIR1)/spin_precess.o \
          $(OBJ_DIR1)/stats.o \
          $(OBJ_DIR1)/summary.o \
          $(OBJ_DIR1)/targ_geom.o \
          $(OBJ_DIR1)/targ_setup.o \
          $(OBJ_DIR1)/tg_path.o \
          $(OBJ_DIR1)/track.o \
          $(OBJ_DIR1)/tr_drift.o \
          $(OBJ_DIR1)/trit_elastic.o \
          $(OBJ_DIR1)/trptvslab.o \
          $(OBJ_DIR1)/tr_rot.o \
          $(OBJ_DIR2)/angfac.o \
          $(OBJ_DIR2)/boost.o \
          $(OBJ_DIR2)/currnt.o \
          $(OBJ_DIR2)/deep.o \
          $(OBJ_DIR2)/forfac.o \
          $(OBJ_DIR2)/gauleg.o \
          $(OBJ_DIR2)/pwia.o \
          $(OBJ_DIR2)/response.o \
          $(OBJ_DIR2)/rotation.o \
          $(OBJ_DIR2)/wave.o \
          $(OBJ_DIR3)/hrs_cosy.o \
          $(OBJ_DIR3)/hrs_hut.o \
          $(OBJ_DIR3)/hrs_recon.o \
          $(OBJ_DIR3)/locforunt.o \
          $(OBJ_DIR3)/musc_ext.o \
          $(OBJ_DIR3)/musc.o \
          $(OBJ_DIR3)/project.o \
          $(OBJ_DIR3)/rotate_haxis.o \
          $(OBJ_DIR3)/transp.o \
 
all_objs = $(fort_objs)
#
$(PROGRAM): $(all_objs)
	rm -f $@
	$(F77) $(LINKFLAGS) $@ $(all_objs) $(LIBS)
 
#
# Get the object files from the sources
#
$(OBJ_DIR1)/arenhoevel.o: $(SRC_DIR1)/arenhoevel.f \
               $(mceep_cmn)/aren_respon.cmn \
               $(mceep_cmn)/aren_ranges.cmn \
               $(mceep_cmn)/var.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/arenhoevel.f
$(OBJ_DIR1)/borie_drechsel.o: $(SRC_DIR1)/borie_drechsel.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/borie_drechsel.f
$(OBJ_DIR1)/bremsgen.o: $(SRC_DIR1)/bremsgen.f \
               $(mceep_cmn)/rada.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/bremsgen.f
$(OBJ_DIR1)/c12_elastic.o: $(SRC_DIR1)/c12_elastic.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/c12_elastic.f
$(OBJ_DIR1)/calc_b.o: $(SRC_DIR1)/calc_b.f \
               $(mceep_cmn)/eloss.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/calc_b.f
$(OBJ_DIR1)/compute_cosy_spin.o: $(SRC_DIR1)/compute_cosy_spin.f \
               $(mceep_cmn)/spin.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/compute_cosy_spin.f
$(OBJ_DIR1)/compute_transp_to_cosy.o: $(SRC_DIR1)/compute_transp_to_cosy.f \
               $(mceep_cmn)/masses.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/compute_transp_to_cosy.f
$(OBJ_DIR1)/decay.o: $(SRC_DIR1)/decay.f \
               $(mceep_cmn)/masses.cmn \
               $(mceep_cmn)/var.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/decay.f
$(OBJ_DIR1)/denscorr.o: $(SRC_DIR1)/denscorr.f \
               $(mceep_cmn)/eloss.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/denscorr.f
$(OBJ_DIR1)/deut_elastic.o: $(SRC_DIR1)/deut_elastic.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/deut_elastic.f
$(OBJ_DIR1)/deut_sabjes.o: $(SRC_DIR1)/deut_sabjes.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/deut_sabjes.f
$(OBJ_DIR1)/electro_prod.o: $(SRC_DIR1)/electro_prod.f \
               $(mceep_cmn)/var.cmn \
               $(mceep_cmn)/masses.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/electro_prod.f
$(OBJ_DIR1)/eloss_e.o: $(SRC_DIR1)/eloss_e.f \
               $(mceep_cmn)/eloss.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/eloss_e.f
$(OBJ_DIR1)/eloss_proton.o: $(SRC_DIR1)/eloss_proton.f \
               $(mceep_cmn)/masses.cmn \
               $(mceep_cmn)/eloss.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/eloss_proton.f
$(OBJ_DIR1)/elossw.o: $(SRC_DIR1)/elossw.f \
               $(mceep_cmn)/eloss.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/elossw.f
$(OBJ_DIR1)/epc.o: $(SRC_DIR1)/epc.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/epc.f
$(OBJ_DIR1)/ext_bremss.o: $(SRC_DIR1)/ext_bremss.f \
               $(mceep_cmn)/eloss.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/ext_bremss.f
$(OBJ_DIR1)/fint.o: $(SRC_DIR1)/fint.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/fint.f
$(OBJ_DIR1)/formfact.o: $(SRC_DIR1)/formfact.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/formfact.f
$(OBJ_DIR1)/get_vect.o: $(SRC_DIR1)/get_vect.f \
               $(mceep_cmn)/wc.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/get_vect.f
$(OBJ_DIR1)/glande.o: $(SRC_DIR1)/glande.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/glande.f
$(OBJ_DIR1)/glandg.o: $(SRC_DIR1)/glandg.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/glandg.f
$(OBJ_DIR1)/glands.o: $(SRC_DIR1)/glands.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/glands.f
$(OBJ_DIR1)/gvaviv.o: $(SRC_DIR1)/gvaviv.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/gvaviv.f
$(OBJ_DIR1)/he3_elastic.o: $(SRC_DIR1)/he3_elastic.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/he3_elastic.f
$(OBJ_DIR1)/he4_elastic.o: $(SRC_DIR1)/he4_elastic.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/he4_elastic.f
$(OBJ_DIR1)/hrs.o: $(SRC_DIR1)/hrs.f \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/hrs.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/hrs.f
$(OBJ_DIR1)/hrs_inv.o: $(SRC_DIR1)/hrs_inv.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/hrs_inv.f
$(OBJ_DIR1)/hyd_elastic.o: $(SRC_DIR1)/hyd_elastic.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/hyd_elastic.f
$(OBJ_DIR1)/in_aperture.o: $(SRC_DIR1)/in_aperture.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/in_aperture.f
$(OBJ_DIR1)/index_accept_fcn.o: $(SRC_DIR1)/index_accept_fcn.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/index_accept_fcn.f
$(OBJ_DIR1)/input.o: $(SRC_DIR1)/input.f \
               $(mceep_cmn)/input.cmn \
               $(mceep_cmn)/masses.cmn \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/wc.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/input.f
$(OBJ_DIR1)/interp.o: $(SRC_DIR1)/interp.f \
               $(mceep_cmn)/summary.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/interp.f
$(OBJ_DIR1)/kinem_elas_ef.o: $(SRC_DIR1)/kinem_elas_ef.f \
               $(mceep_cmn)/masses.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/kinem_elas_ef.f
$(OBJ_DIR1)/kinem_elas_pf.o: $(SRC_DIR1)/kinem_elas_pf.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/kinem_elas_pf.f
$(OBJ_DIR1)/kinem.o: $(SRC_DIR1)/kinem.f \
               $(mceep_cmn)/var.cmn \
               $(mceep_cmn)/masses.cmn \
               $(mceep_cmn)/lifetimes.cmn \
               $(mceep_cmn)/lifetimes_dat.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/kinem.f
$(OBJ_DIR1)/kinfac_2hgross.o: $(SRC_DIR1)/kinfac_2hgross.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/kinfac_2hgross.f
$(OBJ_DIR1)/kinfac_pwvo.o: $(SRC_DIR1)/kinfac_pwvo.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/kinfac_pwvo.f
$(OBJ_DIR1)/Laget_Xsec.o: $(SRC_DIR1)/Laget_Xsec.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/Laget_Xsec.f
$(OBJ_DIR1)/lomon_ff.o: $(SRC_DIR1)/lomon_ff.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/lomon_ff.f
$(OBJ_DIR1)/mad12dfwd_rev.o: $(SRC_DIR1)/mad12dfwd_rev.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mad12dfwd_rev.f
$(OBJ_DIR1)/mad12.o: $(SRC_DIR1)/mad12.f \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/mad.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mad12.f
$(OBJ_DIR1)/mad12_inv.o: $(SRC_DIR1)/mad12_inv.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mad12_inv.f
$(OBJ_DIR1)/mad12nq.o: $(SRC_DIR1)/mad12nq.f \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/mad.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mad12nq.f
$(OBJ_DIR1)/mad12nq_inv.o: $(SRC_DIR1)/mad12nq_inv.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mad12nq_inv.f
$(OBJ_DIR1)/mad35.o: $(SRC_DIR1)/mad35.f \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/mad.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mad35.f
$(OBJ_DIR1)/mad35_funcs.o: $(SRC_DIR1)/mad35_funcs.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mad35_funcs.f
$(OBJ_DIR1)/mad.o: $(SRC_DIR1)/mad.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mad.f
$(OBJ_DIR1)/mad_inv.o: $(SRC_DIR1)/mad_inv.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mad_inv.f
$(OBJ_DIR1)/mad_nq12_f.o: $(SRC_DIR1)/mad_nq12_f.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mad_nq12_f.f
$(OBJ_DIR1)/matrix.o: $(SRC_DIR1)/matrix.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/matrix.f
$(OBJ_DIR1)/mceep.o: $(SRC_DIR1)/mceep.f \
               $(mceep_cmn)/input.cmn \
               $(mceep_cmn)/var.cmn \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/summary.cmn \
               $(mceep_cmn)/labcoord.cmn \
               $(mceep_cmn)/masses.cmn \
               $(mceep_cmn)/wc.cmn \
               $(mceep_cmn)/eloss.cmn \
               $(mceep_cmn)/hrs.cmn \
               $(mceep_cmn)/sep.cmn \
               $(mceep_cmn)/mad.cmn \
               $(mceep_cmn)/rada.cmn \
               $(mceep_cmn)/var_dat.cmn \
               $(mceep_cmn)/spectrometer_dat.cmn \
               $(mceep_cmn)/masses_dat.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mceep.f
$(OBJ_DIR1)/mceep_misc.o: $(SRC_DIR1)/mceep_misc.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mceep_misc.f
$(OBJ_DIR1)/monte_inv.o: $(SRC_DIR1)/monte_inv.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/monte_inv.f
$(OBJ_DIR1)/monte_trans.o: $(SRC_DIR1)/monte_trans.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/monte_trans.f
$(OBJ_DIR1)/mo_tsai.o: $(SRC_DIR1)/mo_tsai.f \
               $(mceep_cmn)/masses.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mo_tsai.f
$(OBJ_DIR1)/mscatt.o: $(SRC_DIR1)/mscatt.f \
               $(mceep_cmn)/eloss.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mscatt.f
$(OBJ_DIR1)/peaking.o: $(SRC_DIR1)/peaking.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/peaking.f
$(OBJ_DIR1)/peek.o: $(SRC_DIR1)/peek.f \
               $(mceep_cmn)/masses.cmn \
               $(mceep_cmn)/var.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/peek.f
$(OBJ_DIR1)/peepi.o: $(SRC_DIR1)/peepi.f \
               $(mceep_cmn)/masses.cmn \
               $(mceep_cmn)/var.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/peepi.f
$(OBJ_DIR1)/physics.o: $(SRC_DIR1)/physics.f \
               $(mceep_cmn)/var.cmn \
               $(mceep_cmn)/lifetimes.cmn \
               $(mceep_cmn)/lifetimes_dat.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/physics.f
$(OBJ_DIR1)/plot.o: $(SRC_DIR1)/plot.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/plot.f
$(OBJ_DIR1)/ps_beam.o: $(SRC_DIR1)/ps_beam.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/ps_beam.f
$(OBJ_DIR1)/pwia_vo.o: $(SRC_DIR1)/pwia_vo.f \
               $(mceep_cmn)/masses.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/pwia_vo.f
$(OBJ_DIR1)/qfsv.o: $(SRC_DIR1)/qfsv.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/qfsv.f
$(OBJ_DIR1)/random.o: $(SRC_DIR1)/random.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/random.f
$(OBJ_DIR1)/ranecu.o: $(SRC_DIR1)/ranecu.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/ranecu.f
$(OBJ_DIR1)/ranges.o: $(SRC_DIR1)/ranges.f \
               $(mceep_cmn)/input.cmn \
               $(mceep_cmn)/var.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/ranges.f
$(OBJ_DIR1)/read_cosy_input_file.o: $(SRC_DIR1)/read_cosy_input_file.f \
               $(mceep_cmn)/spin.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/read_cosy_input_file.f
$(OBJ_DIR1)/r-function.o: $(SRC_DIR1)/r-function.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/r-function.f
$(OBJ_DIR1)/rotate_pol.o: $(SRC_DIR1)/rotate_pol.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/rotate_pol.f
$(OBJ_DIR1)/schwinger.o: $(SRC_DIR1)/schwinger.f \
               $(mceep_cmn)/var.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/schwinger.f
$(OBJ_DIR1)/seed_getnsav.o: $(SRC_DIR1)/seed_getnsav.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/seed_getnsav.f
$(OBJ_DIR1)/sep.o: $(SRC_DIR1)/sep.f \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/hrs.cmn \
               $(mceep_cmn)/sep.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/sep.f
$(OBJ_DIR1)/sep_inv.o: $(SRC_DIR1)/sep_inv.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/sep_inv.f
$(OBJ_DIR1)/septum_inv.o: $(SRC_DIR1)/septum_inv.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/septum_inv.f
$(OBJ_DIR1)/septum_trans.o: $(SRC_DIR1)/septum_trans.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/septum_trans.f
$(OBJ_DIR1)/setup_pion.o: $(SRC_DIR1)/setup_pion.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/setup_pion.f
$(OBJ_DIR1)/sigmacc1.o: $(SRC_DIR1)/sigmacc1.f \
               $(mceep_cmn)/masses.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/sigmacc1.f
$(OBJ_DIR1)/sigmacc2.o: $(SRC_DIR1)/sigmacc2.f \
               $(mceep_cmn)/masses.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/sigmacc2.f
$(OBJ_DIR1)/sig_mscat.o: $(SRC_DIR1)/sig_mscat.f \
               $(mceep_cmn)/masses.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/sig_mscat.f
$(OBJ_DIR1)/smott.o: $(SRC_DIR1)/smott.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/smott.f
$(OBJ_DIR1)/spectral.o: $(SRC_DIR1)/spectral.f \
               $(mceep_cmn)/summary.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/spectral.f
$(OBJ_DIR1)/spectrometer.o: $(SRC_DIR1)/spectrometer.f \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/wc.cmn \
               $(mceep_cmn)/masses.cmn \
               $(mceep_cmn)/input.cmn \
               $(mceep_cmn)/struct_hrs.cmn \
               $(mceep_cmn)/apertures.cmn \
               $(mceep_cmn)/track.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/spectrometer.f
$(OBJ_DIR1)/spin_1st_ord.o: $(SRC_DIR1)/spin_1st_ord.f \
               $(mceep_cmn)/input.cmn \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/masses.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/spin_1st_ord.f
$(OBJ_DIR1)/spin_dat.o: $(SRC_DIR1)/spin_dat.f \
               $(mceep_cmn)/masses.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/spin_dat.f
$(OBJ_DIR1)/spin_precess.o: $(SRC_DIR1)/spin_precess.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/spin_precess.f
$(OBJ_DIR1)/stats.o: $(SRC_DIR1)/stats.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/stats.f
$(OBJ_DIR1)/summary.o: $(SRC_DIR1)/summary.f \
               $(mceep_cmn)/summary.cmn \
               $(mceep_cmn)/var.cmn \
               $(mceep_cmn)/input.cmn \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/masses.cmn \
               $(mceep_cmn)/wc.cmn \
               $(mceep_cmn)/eloss.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/summary.f
$(OBJ_DIR1)/targ_geom.o: $(SRC_DIR1)/targ_geom.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/targ_geom.f
$(OBJ_DIR1)/targ_setup.o: $(SRC_DIR1)/targ_setup.f \
               $(mceep_cmn)/eloss.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/targ_setup.f
$(OBJ_DIR1)/tg_path.o: $(SRC_DIR1)/tg_path.f \
               $(mceep_cmn)/eloss.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/tg_path.f
$(OBJ_DIR1)/track.o: $(SRC_DIR1)/track.f \
               $(mceep_cmn)/wc.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/track.f
$(OBJ_DIR1)/tr_drift.o: $(SRC_DIR1)/tr_drift.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/tr_drift.f
$(OBJ_DIR1)/trit_elastic.o: $(SRC_DIR1)/trit_elastic.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/trit_elastic.f
$(OBJ_DIR1)/trptvslab.o: $(SRC_DIR1)/trptvslab.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/trptvslab.f
$(OBJ_DIR1)/tr_rot.o: $(SRC_DIR1)/tr_rot.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/tr_rot.f
$(OBJ_DIR2)/angfac.o: $(SRC_DIR2)/angfac.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR2)/angfac.f
$(OBJ_DIR2)/boost.o: $(SRC_DIR2)/boost.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR2)/boost.f
$(OBJ_DIR2)/currnt.o: $(SRC_DIR2)/currnt.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR2)/currnt.f
$(OBJ_DIR2)/deep.o: $(SRC_DIR2)/deep.f \
               $(mceep_cmn)/var.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR2)/deep.f
$(OBJ_DIR2)/forfac.o: $(SRC_DIR2)/forfac.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR2)/forfac.f
$(OBJ_DIR2)/gauleg.o: $(SRC_DIR2)/gauleg.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR2)/gauleg.f
$(OBJ_DIR2)/pwia.o: $(SRC_DIR2)/pwia.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR2)/pwia.f
$(OBJ_DIR2)/response.o: $(SRC_DIR2)/response.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR2)/response.f
$(OBJ_DIR2)/rotation.o: $(SRC_DIR2)/rotation.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR2)/rotation.f
$(OBJ_DIR2)/wave.o: $(SRC_DIR2)/wave.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR2)/wave.f
$(OBJ_DIR3)/hrs_cosy.o: $(SRC_DIR3)/hrs_cosy.f \
               $(mceep_cmn)/struct_hrs.cmn \
               $(mceep_cmn)/apertures.cmn \
               $(mceep_cmn)/track.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR3)/hrs_cosy.f
$(OBJ_DIR3)/hrs_hut.o: $(SRC_DIR3)/hrs_hut.f \
               $(mceep_cmn)/struct_hrs.cmn \
               $(mceep_cmn)/track.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR3)/hrs_hut.f
$(OBJ_DIR3)/hrs_recon.o: $(SRC_DIR3)/hrs_recon.f \
               $(mceep_cmn)/track.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR3)/hrs_recon.f
$(OBJ_DIR3)/locforunt.o: $(SRC_DIR3)/locforunt.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR3)/locforunt.f
$(OBJ_DIR3)/musc_ext.o: $(SRC_DIR3)/musc_ext.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR3)/musc_ext.f
$(OBJ_DIR3)/musc.o: $(SRC_DIR3)/musc.f \
               $(mceep_cmn)/ran.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR3)/musc.f
$(OBJ_DIR3)/project.o: $(SRC_DIR3)/project.f \
               $(mceep_cmn)/track.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR3)/project.f
$(OBJ_DIR3)/rotate_haxis.o: $(SRC_DIR3)/rotate_haxis.f \
               $(mceep_cmn)/track.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR3)/rotate_haxis.f
$(OBJ_DIR3)/transp.o: $(SRC_DIR3)/transp.f \
               $(mceep_cmn)/track.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR3)/transp.f
