OPSYS = g77
PROGRAM = mceep_util
DIR_PREFIX = /home/saw/MCEEP/JJ/mceep
#
ifeq ($(OPSYS),hp)
  F77 = f77
  FFLAGS = +E6 +ppu -c -I$(mceep_cmn)
  LINKFLAGS = +E6 +ppu -o
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lgraflib -lgrafX11 -lpacklib \
         -Wl,-L/usr/lib/X11R5 -lX11 -lm
endif
#
ifeq ($(OPSYS),ultrix)
  F77 = f77
  FFLAGS = -vms -c -I$(mceep_cmn)
  LINKFLAGS = -vms -o
  CERN_ROOT = /usr/local/cern/94b
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lgraflib -lgrafX11 -lpacklib \
         -Wl,-L/usr/lib/X11R5 -lX11 -lm
endif
#
ifeq ($(OPSYS),osf1)
  F77 = f77
  FFLAGS = -vms -c -I$(mceep_cmn)
  LINKFLAGS = -vms -o
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lgraflib -lgrafX11 -lpacklib \
         -Wl,-L/usr/lib/X11R5 -lX11 -lm
endif
#
ifeq ($(OPSYS),sun)
  F77 = f77
  FFLAGS = -ext_names=underscores -c -I$(mceep_cmn)
  LINKFLAGS = -o
  LIBS = -L$(CERN_ROOT)/lib \
         -lgraflib -lgrafX11 -lpacklib \
         -L/usr/lib -lsocket -lnsl -lm
endif
#
ifeq ($(OPSYS),absoft)
  F77 = f77
  FFLAGS = -c -s -N22 -N90 -B108 -f -I$(mceep_cmn)
  LINKFLAGS = -X -defsym -X MAIN__=main -o
  LIBS = -L$(CERN_ROOT)/lib \
         -lpacklib \
         -L/usr/lib -lf2c -lm
#        -L/usr/local/farm1/cernlib/99/lib -lf2c -L/usr/lib -lm
endif
#
ifeq ($(OPSYS),g77)
  F77 = g77
  FFLAGS = -c -I$(mceep_cmn)
  LINKFLAGS = -o
  LIBS = -L$(CERN_ROOT)/lib -lpacklib \
         -L/usr/lib/gcc-lib/i486-suse-linux/2.95.2 -lg2c \
         -L/usr/lib -lm
#        -L/usr/local/farm1/cernlib/99/lib -lf2c -L/usr/lib -lm
endif
#
ifeq ($(OPSYS),sgi)
  F77 = f77
  FFLAGS = -vms_cc -c -I$(mceep_cmn)
  LINKFLAGS = -vms_cc -o
  CERN_ROOT = /cern/pro
  LIBS = -Wl,-L$(CERN_ROOT)/lib \
         -lgraflib -lgrafX11 -lpacklib \
         -Wl,-L/usr/lib/X11R5 -lX11 -lm
endif
#
#
SRC_DIR1 = $(DIR_PREFIX)/utilities
OBJ_DIR1 = $(DIR_PREFIX)/utilities
fort_objs = \
          $(OBJ_DIR1)/gauss.o \
          $(OBJ_DIR1)/get_include_cwn.o \
          $(OBJ_DIR1)/make_new_cwn.o \
          $(OBJ_DIR1)/make_new_rwn.o \
          $(OBJ_DIR1)/mceep_to_hbook.o \
          $(OBJ_DIR1)/mceep_util.o \
          $(OBJ_DIR1)/ratio.o \
          $(OBJ_DIR1)/r-function.o \
 
all_objs = $(fort_objs)
#
$(PROGRAM): $(all_objs)
	rm -f $@
	$(F77) $(LINKFLAGS) $@ $(all_objs) $(LIBS)
 
#
# Get the object files from the sources
#
$(OBJ_DIR1)/gauss.o: $(SRC_DIR1)/gauss.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/gauss.f
$(OBJ_DIR1)/get_include_cwn.o: $(SRC_DIR1)/get_include_cwn.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/get_include_cwn.f
$(OBJ_DIR1)/make_new_cwn.o: $(SRC_DIR1)/make_new_cwn.f \
               $(mceep_cmn)/ntu.cmn \
               $(mceep_cmn)/ntu_user.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/make_new_cwn.f
$(OBJ_DIR1)/make_new_rwn.o: $(SRC_DIR1)/make_new_rwn.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/make_new_rwn.f
$(OBJ_DIR1)/mceep_to_hbook.o: $(SRC_DIR1)/mceep_to_hbook.f \
               $(mceep_cmn)/var.cmn \
               $(mceep_cmn)/spectrometer.cmn \
               $(mceep_cmn)/var_dat.cmn \
               $(mceep_cmn)/spectrometer_dat.cmn
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mceep_to_hbook.f
$(OBJ_DIR1)/mceep_util.o: $(SRC_DIR1)/mceep_util.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/mceep_util.f
$(OBJ_DIR1)/ratio.o: $(SRC_DIR1)/ratio.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/ratio.f
$(OBJ_DIR1)/r-function.o: $(SRC_DIR1)/r-function.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $(SRC_DIR1)/r-function.f
