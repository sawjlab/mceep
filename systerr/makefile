DIR_PREFIX = /home/saw/MCEEP/JJ/mceep
OPSYS = g77
#OPSYS = g77
PROGRAM = all
#DIR_PREFIX = /Users/ulmer/mceep
#
ifeq ($(OPSYS),hp)
  F77 = f77
  FFLAGS = +E6 +ppu -c -I$(mceep_cmn)
  LINKFLAGS = +E6 +ppu -o
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lgraflib -lgrafX11 -lpacklib -lmathlib \
         -Wl,-L/usr/lib/X11R5 -lX11 -lm
endif
#
ifeq ($(OPSYS),ultrix)
  F77 = f77
  FFLAGS = -vms -c -I$(mceep_cmn)
  LINKFLAGS = -vms -o
  CERN_ROOT = /usr/local/cern/94b
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lgraflib -lgrafX11 -lpacklib -lmathlib \
         -Wl,-L/usr/lib/X11R5 -lX11 -lm
endif
#
ifeq ($(OPSYS),osf1)
  F77 = f77
  FFLAGS = -vms -c -I$(mceep_cmn)
  LINKFLAGS = -vms -o
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lgraflib -lgrafX11 -lpacklib -lmathlib \
         -Wl,-L/usr/lib/X11R5 -lX11 -lm
endif
#
ifeq ($(OPSYS),sun)
  F77 = f77
  FFLAGS = -ext_names=underscores -c -I$(mceep_cmn)
  LINKFLAGS = -o
  LIBS = -L$(CERN_ROOT)/lib \
         -lgraflib -lgrafX11 -lpacklib -lmathlib \
         -L/usr/lib -lsocket -lnsl -lm
endif
#
ifeq ($(OPSYS),absoft)
  F77 = f77
  FFLAGS = -c -s -N22 -N90 -B108 -f -I$(mceep_cmn)
  LINKFLAGS = -X -defsym -X MAIN__=main -o
  LIBS = -L$(CERN_ROOT)/lib \
         -lpacklib -lmathlib \
         -L/usr/lib -lf2c -lm
endif
#
ifeq ($(OPSYS),g77)
  F77 = g77
  FFLAGS = -c -I$(mceep_cmn)
  LINKFLAGS = -o
  LIBS = -L$(CERN_ROOT)/lib -lpacklib -lmathlib \
         -L/usr/lib -lm
endif
#
ifeq ($(OPSYS),sgi)
  F77 = f77
  FFLAGS = -vms_cc -c -I$(mceep_cmn)
  LINKFLAGS = -vms_cc -o
  CERN_ROOT = /cern/pro
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lgraflib -lgrafX11 -lpacklib -lmathlib \
         -Wl,-L/usr/lib/X11R5 -lX11 -lm
endif
#
#
OBJ_DIR1 = $(DIR_PREFIX)/sources
OBJ_DIR2 = $(DIR_PREFIX)/gross
OBJ_DIR3 = $(DIR_PREFIX)/cosy
#
SRC_DIR4 = $(DIR_PREFIX)/systerr
OBJ_DIR4 = $(DIR_PREFIX)/systerr

fort_objs1 = \
          $(OBJ_DIR1)/Laget_Xsec.o \
          $(OBJ_DIR1)/arenhoevel.o \
          $(OBJ_DIR1)/borie_drechsel.o \
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
          $(OBJ_DIR1)/kinem.o \
          $(OBJ_DIR1)/kinem_elas_ef.o \
          $(OBJ_DIR1)/kinem_elas_pf.o \
          $(OBJ_DIR1)/kinfac_2hgross.o \
          $(OBJ_DIR1)/kinfac_pwvo.o \
          $(OBJ_DIR1)/lomon_ff.o \
          $(OBJ_DIR1)/mad.o \
          $(OBJ_DIR1)/mad12.o \
          $(OBJ_DIR1)/mad12_inv.o \
          $(OBJ_DIR1)/mad12dfwd_rev.o \
          $(OBJ_DIR1)/mad12nq.o \
          $(OBJ_DIR1)/mad12nq_inv.o \
          $(OBJ_DIR1)/mad35.o \
          $(OBJ_DIR1)/mad35_funcs.o \
          $(OBJ_DIR1)/mad_inv.o \
          $(OBJ_DIR1)/mad_nq12_f.o \
          $(OBJ_DIR1)/matrix.o \
          $(OBJ_DIR1)/mceep_misc.o \
          $(OBJ_DIR1)/mo_tsai.o \
          $(OBJ_DIR1)/monte_inv.o \
          $(OBJ_DIR1)/monte_trans.o \
          $(OBJ_DIR1)/mscatt.o \
          $(OBJ_DIR1)/peaking.o \
          $(OBJ_DIR1)/peek.o \
          $(OBJ_DIR1)/peepi.o \
          $(OBJ_DIR1)/physics.o \
          $(OBJ_DIR1)/plot.o \
          $(OBJ_DIR1)/ps_beam.o \
          $(OBJ_DIR1)/pwia_vo.o \
          $(OBJ_DIR1)/qfsv.o \
          $(OBJ_DIR1)/r-function.o \
          $(OBJ_DIR1)/random.o \
          $(OBJ_DIR1)/ranecu.o \
          $(OBJ_DIR1)/ranges.o \
          $(OBJ_DIR1)/read_cosy_input_file.o \
          $(OBJ_DIR1)/rotate_pol.o \
          $(OBJ_DIR1)/schwinger.o \
          $(OBJ_DIR1)/seed_getnsav.o \
          $(OBJ_DIR1)/sep.o \
          $(OBJ_DIR1)/sep_inv.o \
          $(OBJ_DIR1)/septum_inv.o \
          $(OBJ_DIR1)/septum_trans.o \
          $(OBJ_DIR1)/setup_pion.o \
          $(OBJ_DIR1)/sig_mscat.o \
          $(OBJ_DIR1)/sigmacc1.o \
          $(OBJ_DIR1)/sigmacc2.o \
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
          $(OBJ_DIR1)/tr_drift.o \
          $(OBJ_DIR1)/tr_rot.o \
          $(OBJ_DIR1)/track.o \
          $(OBJ_DIR1)/trit_elastic.o \
          $(OBJ_DIR1)/trptvslab.o \
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
          $(OBJ_DIR4)/systerr.o \

fort_objs2 = \
          $(OBJ_DIR4)/toterr.o \

all_objs = $(fort_objs)
#
$(PROGRAM): systerr toterr

systerr: $(fort_objs1)
	rm -f $@
	$(F77) $(LINKFLAGS) $@ $(fort_objs1) $(LIBS)

toterr: $(fort_objs2)
	rm -f $@
	$(F77) $(LINKFLAGS) $@ $(fort_objs2)

#
# Get the object files from the sources
#
$(OBJ_DIR4)/systerr.o: $(SRC_DIR4)/systerr.f \
               $(mceep_cmn)/var.cmn \
               $(mceep_cmn)/labcoord.cmn \
               $(mceep_cmn)/masses.cmn \
               $(mceep_cmn)/masses_dat.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR4)/systerr.f
$(OBJ_DIR4)/toterr.o: $(SRC_DIR4)/toterr.f
	rm -f $@
	$(F77) -o toterr.o $(FFLAGS) $(SRC_DIR4)/toterr.f


