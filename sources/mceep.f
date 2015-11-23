C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C  ***  PLEASE DO NOT DISTRIBUTE MODIFIED VERSIONS OF THIS PROGRAM  ***
C  ***             WITHOUT PERMISSION OF THE AUTHOR.                ***
C------------------------------------------------------------------------------
C
C       PROGRAM:  MCEEP
C       VERSION:  3.9
C
C       AUTHOR:      Paul E. Ulmer
C                    Department of Physics
C                    Old Dominion University
C                    Norfolk, VA 23529
C
C                    Phone:   (757) 683-5851
C                    Email:   ulmer@jlab.org
C
C       DATE:   V1.0    19-NOV-1991
C       MODIFICATIONS:
C               Version 1.01    06-DEC-1991
C                  New Options:
C                      Elastic scattering from deuterium - RWL/SVV
C                  Modules added:
C                      deut_elastic
C                  Modules changed:
C                      mceep,physics,spectral,spin_1st_ord,summary
C               Version 1.1     02-AUG-1992
C                  New Options:
C                      Track reconstruction for wire chambers - MN
C                      Relative ToF option - MN
C                      Second order Transport - MN
C                      Accidentals rates - MN
C                  Modules added:
C                      get_vect,track
C                  Modules changed: (in addition to .cmn files)
C                      input,mceep,spectrometer,summary,tr_rot,tr_drift,
C                      mceep_misc
C               Version 1.2     28-OCT-1992
C                  Changes:
C                      Fixed portability problem by replacing variable
C                         format descriptor by subroutine which reads
C                         character data from input line       - PEU
C                      Dipole nucleon form factors used throughout
C                         with GEn having the "correct" asymptotic
C                         form in routine formfact             - PEU
C                      Uses RANECU uniform random number generator,
C                         rather than VAX generator.  Now MCEEP uses
C                         only a single initial seed (actually
C                         pair of seeds).                      - PEU
C                  Modules changed:
C                      formfact,get_vect,hyd_elastic,kinem,input,
C                      mceep,plot,pwia_vo,random,ranges,setup_pion,
C                      sigmacc1,spectrometer,summary
C                  Modules added:
C                      ranecu,read_input_line,seed_getnsav
C               Version 2.0     07-OCT-1994
C                  Changes:
C                      Added N-Tuple option to allow replaying MCEEP
C                         data with PAW.  Currently, a file is written
C                         which then needs to be reformatted (outside
C                         MCEEP) for PAW.
C                  Modules added:
C
C               Version 2.1     09-SEP-1997
C                  Changes:  (Some of these changes have already been
C                             included in an unofficial release - these
C                             are marked with an *)
C                      ELECTRO_PROD has been modified by Glen Warren
C                             since previously a factor of i was missing
C                             from the born terms and the polarization
C                             observables were never boosted from the COM
C                             frame to the LAB frame.
C                      * Physics option 300 (d(e,e'p)n Gross/Van Orden)
C                             has been updated.
C                      * Option to bypass spectrometer traceback included
C                             so that cross sections and derived kinematics
C                             are calculated from the original coordinates.
C                      PEEPI has been modified to include boost to LAB frame.
C                      Now longer uses READ_INPUT_LINE as all reads from
C                             input file are now free format.
C                             (READ_INPUT_LINE was the routine which
C                             got around the use of variable format
C                             descriptors.)  This change to free format
C                             was required for portability to SUN.
C
C               Version 2.2     24-APR-1998
C                  Changes:
C                      DEC-VMS no longer supported
C                             This affects the method of distribution
C                             and not the sources (although logical
C                             names have now been replaced by explicit
C                             Unix directory names in the supplied sources
C                             rather than having this replacement done
C                             at "build-time").
C                      Linux (Absoft f77) now supported
C                             Only minor changes were required.
C
C               Version 2.30    30-MAR-1999
C                  Changes:
C                      Paris deuteron wave function added in spectral.f
C                      JLAB-HRS spectrometers (apertures and tgt-fp maps)
C                             have been added.
C                             New files:  monte_trans.f, hrs.f
C                             Modified:   input, spectrometer, summary
C
C               Version 3.0     11-JUL-1999
C                  Changes:
C                      Radiative effects added (internal and external rad.)
C                      Ionization energy loss and multiple scattering added.
C                      Multiple target foils now supported as well as
C                         standard JLAB Hall A cryotarget.
C                      p(e,e'K)Lambda physics option included.
C                          New files:
C                                    (Paul Ulmer:)
C                                         peaking.f, ps_beam.f
C                                         kinem_elas_ef.f
C                                         kinem_elas_pf.f
C                                         borie_drechsel.f
C                                         mo_tsai.f
C                                         targ_setup.f
C                                         sig_mscat.f
C                                         (kinem_elas.f is replaced by
C                                          kinem_elas_ef.f and kinem_elas_pf.f)
C                                    (Luminita Todor:)
C                                         calc_b.f, denscorr.f
C                                         eloss_e.f, eloss_proton.f, elossw.f
C                                         ext_bremss.f, mscatt.f, tg_path.f
C                                         (The next four routines are based
C                                         on CERNLIB routines.)
C                                         glande.f, glandg.f, glands.f
C                                         gvaviv.f
C                                    (Mark Jones:)
C                                         schwinger.f
C                                    (Pete Markowitz:)
C                                         peek.f
C                                         (From CERNLIB:) fint.f
C                                         peepi.f now split into
C                                           peepi.f and decay.f (PEU)
C                          Modified:      MANY routines
C                             
C               Version 3.1     18-NOV-1999
C                  Changes:
C                      JLAB-Hall A tuna can target model (3,4He) added.
C                      3He elastic scattering option added (A. Deur).
C                      Beam raster option added (F. Sabatie).
C                      LeRose electron arm inverse functions updated.
C                      "Bowtie" effect now corrected after applying
C                         LeRose inverse functions.
C                      Now allows option of including the multi-photon
C                         correction to the radiative tail (J. Templon
C                         et al. method).
C                      Photon cutoff energy is now read in instead of
C                         hardwired to 1 MeV (1 MeV is the default though).
C                      Bug fixes in energy loss routines for either
C                         non-standard target nuclei or neutral ejectiles
C                         (now neutral ejectiles have no mult. scatt. or
C                          energy loss - a "cheap" fix).
C                      Bug fix for missing mass, since for VCS half the
C                         events are thrown out as non-physical due to
C                         roundoff/resolution effects.
C                      Bug fix for spectrometer options.  Now, TRK,
C                         MAT, ROT, DFT and TOF are only called for
C                         kinematics (ICALL=2).  This is necessary for
C                         consistency among all these options.  This
C                         bug has been known to lead to serious errors
C                         when the TRK and HRS options were combined.
C                      Bug fix for file names.  Now accomodates longer
C                         paths/names.
C                      Bug fixes for p(e,e'K) - physics option 700.
C                      Bug fix:  avg. cross section in summary file
C                         used to be calculated by taking phase space
C                         weighted sum of sigma divided by sum of unity
C                         for each event.  Now I correctly divide by
C                         sum of phase space (with radiation this can
C                         make a BIG difference).
C                      Bug fix:  for radiative tail, the unradiated
C                         cross section in Mo/Tsai and Borie/Drechsel
C                         should be calculated using vertex kinematics.
C                         This was done, except that the beam energy
C                         and scattered electron energy were the
C                         asymptotic values.  Now everything is vertex.
C                      Bug fix:  for LH2/LD2 target and spectrometer angles
C                         between -90 and -180 deg, the target path (for
C                         energy loss calculation) was incorrect.
C                         This fix was incorporated after initial release.
C                      Bug fix:  for Meier-Hadjuk 3He spectral function
C                         events with recoil momentum between 0 and 5 MeV/c
C                         were thrown out.  Now the spectral function at
C                         5 MeV/c is used for all momenta in this range.
C                         This fix was incorporated after initial release.
C                      Bug fix:  schwinger routine converted Q^2 from
C                         (GeV/c)^2 to (MeV/c)^2 by mult. by 10^3 instead
C                         of 10^6.
C                         This fix was incorporated after initial release.
C                      Now calculate focal plane vector even if event fails
C                         the aperture cuts (hrs.f).  This makes the before
C                         cuts counters physically meaningful.
C                         This fix was incorporated after initial release.
C                      New files:
C                            he3_elastic.f  (A. Deur)
C                      Modified routines:
C                            calc_b.f, denscorr.f, eloss_e.f,
C                            eloss_proton.f, elossw.f, ext_bremss.f,
C                            forfac.f, get_vect.f, glandg.f,
C                            hrs.f, hrs_inv.f, input.f,
C                            kinem.f, mceep.f, mceep_misc.f, monte_inv.f,
C                            mscatt.f, peek.f, physics.f, plot.f,
C                            schwinger.f, seed_getnsav.f, setup_pion.f,
C                            spectral.f, spectrometer.f,
C                            summary.f, targ_setup.f, tg_path.f
C
C               Version 3.2     11-JAN-2000
C                  Changes:
C                      Bug fix for eloss_proton.f:  emax was calculated
C                            from variables with incompatible energy scales.
C                      Improvement of denscorr.f:  C0 now corrected for
C                            difference in A as well as density.  Also
C                            corrected two bug fixes in denscorr.
C                      Added Salme 3He(e,e'p) spectral function.
C                      Added kinematics variables:
C                            components of PREC in LAB system
C                            energy loss variables
C                      Now specify object point of each spectrometer
C                            allowing for possible mispointings.
C                            The BEAM_E1-BEAM_P3 coordinates now refer
C                            to the vertex coordinates in a modified LAB
C                            system with new origin at the object point.
C                            VERTEX_X-Z refer to the actual vertex
C                            coordinates in the unmodified LAB system.
C                      Incorporate ejectile spin precession via COSY
C                            spin matrix (Mark Jones).  Also produce
C                            output for helicity dependent and helicity
C                            independent pieces of each polarization
C                            component (Mark Jones).
C                      New files:
C                            compute_cosy_spin.f,
C                            compute_transp_to_cosy.f,
C                            read_cosy_input_file.f
C                      Modified routines:
C                            denscorr.f, deut_elastic.f,
C                            electro_prod.f, eloss_proton.f, 
C                            he3_elastic.f, hyd_elastic.f, input.f,
C                            kinem.f, mceep.f, mceep_misc.f,
C                            pwia_vo.f, physics.f,
C                            spectral.f, spectrometer.f, summary.f
C
C               Version 3.3     22-FEB-2000
C                  Changes:
C                      Allow for beam centroid offsets.
C                            Also account for these offsets in the
C                            case of target foils at non-zero angle.
C                      Allow for smearing of beam position with gaussian.
C                      Allow Kinematics Ntuples to include initial and
C                            final Transport coordinates.  This allows
C                            regular kinematic quantities to be cut
C                            on these variables via Ntuples.
C                      Increase maximum # variables per Ntuple
C                            from 20 to 100.
C                      Increase maximum # of elements per spectrometer
C                            from 30 to 50.
C                      Fixed bug in beam dE/dx correction (previously
C                            there was NO CORRECTION due to bug).
C                      Set focal plane vector to zero for events
C                            which fail HRS apertures.  This avoids
C                            numerical errors on some platforms.  The
C                            eloss/mscatt routines occasionally give
C                            large effects for ray coordinates and 
C                            the LeRose inverse functions are not
C                            designed to deal with events way outside
C                            the acceptance.  Note that these events
C                            fail the cuts anyway. 
C                      Modified 3He elastic routine (D. Higinbotham).
C                      Added Tritium, 4He and 12C elastic
C                            routines (D. Higinbotham).
C                      Modified nuclear masses in mceep.f (D. Higinbotham).
C                      Modified routines:
C                            he3_elastic.f, hrs.f, input.f, mceep.f,
C                            physics.f, plot.f, summary.f
C                      Added routines:
C                            c12_elastic.f, he4_elastic.f, trit_elastic.f
C
C               Version 3.4     02-NOV-2000
C                  Changes:
C                      Fixed "bowtie" correction and also correct for
C                            explicit X0 due to vertical raster,
C                            spectrometer mispointing, etc.
C                            Also, now one can use the "electron"
C                            spectrometer or "hadron" spectrometer
C                            for detection of either particle.
C                      GEP (for FFTYPE='DIPOLE') now employs fit to
C                            JLAB-Hall A GEP/GMP ratio in E93027
C                            (GMP is taken to be standard dipole
C                            form factor times MUP)
C                            (S. Dieterich/S. Strauch).
C                      Several new deuteron momentum distributions
C                            from Sabine Jeschonnek: Bonn, Argonne V18,
C                            CD Bonn, Gross/Van Orden.
C                      Now allow histogramming of focal plane coords.
C                            in addition to coords. before and after
C                            call to SPECTROMETER.
C                      MCEEP now beeps when it finishes (just like
C                            ESPACE - amazing).
C                      Routines SPLINE and SPLINT in wave.f renamed
C                            to GROSS_SPLINE and GROSS_SPLINT to avoid
C                            conflict with routines of EPIPROD (J.J. Kelly)
C                            with same names.
C                      Add histograms/ntuples of various kinematical
C                            quantities defined with respect to the
C                            virtual photon vertex (as opposed to the
C                            normally measured kinematical quantities
C                            which reflect radiation).
C                      Incorporate LeRose's new forward functions which
C                            have improved x0 dependence.
C                      Increased precision for Ntuple write statement.
C                      Fixed bug in elossw where the energy loss in the
C                            endcap (for target model=2) was negative
C                            for abs(proton_angle)>90.
C                      Added components of polarizations in scattering
C                            plane (total, hel. dep. and hel. indep.)
C                            (S. Strauch); changed tag for electron
C                            analyzing power from 8 to 11.
C                      Modified routines:
C                            elossw.f, formfact.f, hrs_inv.f,
C                            input.f (comments only), mceep.f,
C                            mceep_misc.f, monte_trans.f, physics.f,
C                            rotate_pol.f, spectral.f, spectrometer.f,
C                            summary.f, wave.f
C
C               Version 3.5     7-FEB-2001
C                  Changes:
C                      Incorporated Arenhoevel interpolated response
C                            functions for d(e,e'p)n.
C                      Added multiple weight Ntuple ('NTM').
C                  Changes below made after the official release of
C                  Version 3.5:
C                      Changed REAL-->DREAL in deut_sabjes.f to fix
C                            problems with g77 compiler (K. Fissum, 
C                            S. Jeschonnek).  Also fixed bug in
C                            get_deut_bonn (removed the argument list).
C                      Added angles of recoil momentum vector in LAB.
C                      Fixed bug in targ_geom.f.  The solid angle
C                            Jacobian had an error (pointed out by
C                            Marat Rvachev).  This error is very small
C                            in most cases (much less than 1%).
C                            (13-MAR-2001). 
C                      Added spectral functions for 16O from Udias and
C                            for 208Pb from Lapikas.  These were all
C                            coded by K. Fissum.
C                            (16-APR-2001). 
C                      Modified routines:
C                            deut_sabjes.f, input.f, kinem.f, mceep.f,
C                            mceep_misc.f, 
C                            physics.f, plot.f,
C                            ranges.f, spectral.f, summary.f, targ_geom.f
C                      Added routines:
C                            arenhoevel.f
C
C               Version 3.6     28-JUN-2001
C                  Changes:
C                      Allow histogramming of X and Y at HRS apertures.
C                      Changed default Q3 radius for HRS from 0.30 m
C                            to 0.28 m.
C                      Get arenhoevel responses from directory defined
C                            by environment variable mceep_aren
C                      Added JLAB-HRS R-functions of J. LeRose
C                      Fixed bug:  for energy loss with dE/dx correction
C                            the acceptances were enforced after making
C                            this correction, whereas they are now
C                            enforced prior to dE/dx correction.
C                      Include updated p(e,e'K) routine from P. Markowitz
C                      Modified routines:
C                            arenhoevel.f, hrs.f, input.f, mceep.f, 
C                            peek.f, spectrometer.f, summary.f
C                      Added routines:
C                            hrs.cmn, r-function.f
C
C               Version 3.7     3-MAY-2002
C                  Changes:
C                      Nucleon form factor routine, formfact.f, now
C                            refers to standard dipole as FFTYPE='DIPOLE',
C                            whereas the Hall A Fit to GE/GM is given by
C                            FFTYPE='HALLA1' (uses dipole GMp).  Another
C                            fit was also added called 'HALLA2' which
C                            uses alternate extraction of GMp and also the 
C                            Hall A GE/GM ratio to get GEp.  A third
C                            fit called 'HALLA3' (now the default)
C                            uses an alternate MMD form factor for GMp
C                            and the Hall A GE/GM ratio to get GEp.
C                            Other parametrizations added too.
C                      New LeRose transfer functions (forward and
C                            reverse) and new r-functions (data files
C                            only).
C                      Added Hall A MAD spectrometer (K. McCormick).
C                      Added Hall A polarized 3He target (K. McCormick).
C                      Modified routines:
C                            elossw.f, formfact.f, hyd_elastic.f, 
C                            input.f, mceep.f,
C                            monte_trans.f, monte_inv.f, pwia_vo.f,
C                            sigmacc1.f, spectrometer.f, summary.f,
C                            targ_setup.f, tg_path.f,
C                            ~/mceep/dat/r-function,
C                            ~/cmn/spectrometer.cmn
C                            ~/cmn/var_dat.cmn
C                      Added routines:  mad.f, mad35_funcs.f
C                            ~/cmn/mad.cmn
C
C               Version 3.8     24-APR-2003
C                  Changes:
C                      COSY model of JLab Hall A HRS (adapted from
C                            Hall C program, SIMC) added (W. Hinton).
C                      Small change made in r-function.f 
C                            (both in sources and in utilities directories) 
C                            allowing the program to compile on Sun
C                            systems.
C                      Bug fix in tg_path
C                      Modified routines:
C                            input.f, r-function.f, spectrometer.f, 
C                            summary.f,tg_path.f
C                            ~/mceep/utilities/r-function.f
C                            ~/mceep/cmn/spectrometer.cmn 
C                      Added routines:  
C                            all routines in directory ~/mceep/cosy,
C                            ~/cmn/apertures.cmn,
C                            ~/cmn/struct_hrs.cmn,
C                            ~/cmn/ran.cmn,
C                            ~/cmn/track.cmn,
C                            ~/dat/hrs_forward_cosy.dat,
C                            ~/dat/hrs_recon_cosy.dat
C                      06-AUG-2003:  P.E. Ulmer
C                            Fixed bug:  variable "fry" was never set
C                            This is needed by hrs_recon and should
C                            be the vertical position at the target.
C                            It's now set in spectrometer.f.  Previously
C                            the code would give NaN on occasion and
C                            obviously a bogus reconstruction.
C
C               Version 3.9     14-JUN-2006
C                  Changes:
C                      Added calculation of x_tg(5) in hrs_inv.f.
C                            Also, now use this information to get the
C                            interaction point in the target as seen
C                            by the spectrometers. (i.e. the BEAM_E1,
C                            BEAM_E2 and BEAM_E3 variables now refer to
C                            the interaction point, rather than the point
C                            where the particle intersects the target
C                            Transport (Z=0) plane. Same for proton variables).
C                      For the case of foils, the various walls were not
C                            handled properly in energy loss routines.
C                            There were several bugs, now fixed.
C                      The formatting of the cosy source files has been
C                            standardized for compatibility with the
C                            latest ABSOFT compiler (W. Boeglin).
C                      Added bremsstrahlung routine for elastic ep by
C                            Florian Weissbach and Kai Hencken.
C                      Added deForest sigma_CC2 offshell 
C                            cross section (P. Monaghan).
C                      Added cigar tube LH2/LD2 target (Model 5)
C                            (Hassan Ibrahim, December 7, 2004).
C                      Added MAD 12 degree configs. (std. and no quad tunes).
C                      Now allow for calculation of spectrometer
C                            acceptance functions (must bypass kinem
C                            in case of failure there).
C                      Fixed bug:  Window energy loss was doubly included!
C                      Fixed bug in c12_elastic routine (R. Feuerbach).
C                      Fixed bug in harmonic oscillator momentum distributions
C                            as the model numbers were incorrect.
C                      Added option to correct for most probable energy
C                            loss (mean energy loss correction was only
C                            option before)
C                      Added Laget grid for d(e,e'p) (E. Voutier/P. Ulmer)
C                      Fixed some portability issues with Sun-OS.
C                      Added the septum optics model ala' J. LeRose.
C                            Should be analogous to HRS model (P. Markowitz)
C                      Modified routines:
C                            calc_b.f, c12_elastic.f,
C                            eloss_e.f, eloss_proton.f,
C                            eloss_w.f, hrs_inv.f, input.f,
C                            kinem_elas_ef.f, mad.f, mceep.f, 
C                            peaking.f, physics.f,
C                            random.f, r-function.f, spectrometer.f,
C                            summary.f, 
C                            targ_setup.f, tg_path.f, trptvslab.f,
C                            ALL sources in mceep/cosy
C                            ~/cmn/eloss.cmn
C                            ~/cmn/input.cmn
C                            ~/cmn/mad.cmn
C                            ~/cmn/spectrometer.cmn
C                            ~/cmn/var.cmn
C                            ~/cmn/var_dat.cmn
C                      Added routines: bremsgen.f, index_accept_fcn.f, 
C                            Laget_Xsec.f,
C                            mad12.f, mad35.f,
C                            mad12dfwd_rev.f, mad_inv.f, mad12_inv.f,
C                            sep.f, sep_inv.f, septum_inv.f, septum_trans.f,
C                            sigmacc2.f,
C                            ~/cmn/rada.cmn
C                            ~/cmn/sep.cmn
C                            
C       PURPOSE:
C               Simulate yields, cross sections, polarizations and
C               asymmetries for (e,e'X) by Monte Carlo sampling.
C
C               Bound state and continuum scattering for the A-1 system
C               are handled as well as elastic scattering in A(e,e'A).
C               For elastic scattering one has the option to ignore
C               the hadron arm acceptance to simulate single-arm elastic
C               scattering.  Pion electroproduction from the proton
C               in p(e,e'p)pi0 is also handled.
C
C               The sampling variables are (for the general 6-dimensional
C               problem):
C                       PF_E            Electron momentum
C                       TH_E            Electron vertical angle
C                       PH_E            Electron horiz. angle
C
C                       PF_P            Ejectile momentum
C                       TH_P            Ejectile vertical angle
C                       PH_P            Ejectile horiz. angle
C
C               For bound states or elastic scattering the sampled
C               variables are a subset of the above; the rest are
C               determined from the kinematic constraints.
C
C               The code allows for an arbitrary geometry so that
C               out-of-plane experiments with either elevated detectors
C               or beam swingers can be simulated.  In addition,
C               extended target geometries are provided.  The
C               program also provides for a description of spectrometer
C               lines including effects from misalignments and resolution
C               smearing.  The spectrometer orientation is arbitrary
C               so that vertically, horizontally or otherwise bending
C               spectrometers can be incorporated.
C
C               Finally, the program allows for various cuts which can
C               either be applied globally or can be attached to specific
C               histograms.
C
C       INPUT:  See SUBROUTINE INPUT and comments therein for a complete
C               description of all input variables.
C
C
C       OUTPUT:  Topdrawer plottable files of user defined histograms.
C                Summary file (prefix = input file prefix, ext. = .sum).
C
C       ROUTINES CALLED:  definitely.
C
C------------------------------------------------------------------------------
C

      PROGRAM MCEEP
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
C
      COMMON /KINVAR/ KI,Q,Q2,QMU2,EP,PF,THETA_EP,CTH_PQ,PR
      COMMON /ELIMS/ EF_MIN,EF_MAX
      COMMON /PLIMS/ PF_MIN,PF_MAX
      COMMON /PFROOTS_D/ PF_ROOT_WT
      COMMON /PFROOTS_L/ FIRST_TIME
      COMMON /PFROOTS_I/ I_PF_ROOT
      COMMON /PARTICLE/ PROTON
      COMMON /L_TO_S/ ROT_P
      COMMON /TOF/ TOF_REL_TMP,TOFE,TOFP
      COMMON /RADPKAPPR/ cutoff,eg_max
      COMMON /RADFULL_L/ RADFULL
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
      COMMON /XTGT/ X_BEAM_E,X_BEAM_P,Y_BEAM_E,Y_BEAM_P
      COMMON /FPVEC/ FP_VEC_E,FP_VEC_P
      COMMON /ARENHOVEL/ FAIL_CELL
      COMMON /ARENFAIL/ I_AREN_FAIL
      COMMON /RFN/ E_RFN,P_RFN
      COMMON / LAGET_CNT/ LAGET_PS_FAIL,LAGET_GRID_FAIL

      DOUBLE PRECISION KI(3),Q(3),PF(3),PR(3)
      DOUBLE PRECISION SIG_DEN(50,500),SIG_DEN_2D(50,50,50)
      DOUBLE PRECISION SIGMA(0:200),SIGMA_MULT_WT(10)
      DOUBLE PRECISION NTARG,QI_E_SAV(6),QI_P_SAV(6)
      DOUBLE PRECISION N_A,E_RFN,P_RFN
      DOUBLE PRECISION KI_VEC(3),KF_VEC(3),KI_PS(3)
      DOUBLE PRECISION KG(3),EG_PK(2),KG_PK(2,3)
      DOUBLE PRECISION KG_UNIT(3),KG_PK_UNIT(2,3)
      DOUBLE PRECISION PF_E_I_TMP(2),PF_P_I_TMP(2)
      DOUBLE PRECISION PH_P_I_TMP(2),TH_P_I_TMP(2),DEF_DK_TMP(2)
      DOUBLE PRECISION E0_ASM,E0_VTX,PH_B_VTX,TH_B_VTX,KI_VTX(3)
      DOUBLE PRECISION EF_ASM,EF_VTX,PH_E_VTX,TH_E_VTX,KF_VTX(3)
      DOUBLE PRECISION PF_ASM,PF_VTX,PH_P_VTX,TH_P_VTX,KP_VTX(3)
      DOUBLE PRECISION RAD_NORM(2)
      DOUBLE PRECISION PF_ROOT_WT_TMP(2)
      DOUBLE PRECISION POL(3),POL_SPEC(3)
C
C     Helicity dependent and independent polarizations
C
      DOUBLE PRECISION POL_SHD(3),POL_SHI(3)
      DOUBLE PRECISION POL_SPEC_HD(3),POL_SPEC_HI(3)
      DOUBLE PRECISION POL_N_RHI,POL_T_RHI,POL_L_RHI
      DOUBLE PRECISION POL_N_RHD,POL_T_RHD,POL_L_RHD
      DOUBLE PRECISION POL_N_SHI,POL_T_SHI,POL_L_SHI
      DOUBLE PRECISION POL_N_SHD,POL_T_SHD,POL_L_SHD
      DOUBLE PRECISION POL_X_HI, POL_Y_HI, POL_Z_HI
      DOUBLE PRECISION POL_X_HD, POL_Y_HD, POL_Z_HD
      DOUBLE PRECISION POL_X, POL_Y, POL_Z
C
      DOUBLE PRECISION ROTX(3,3),ROTY(3,3),ROTZ(3,3),ROTXY(3,3)
      DOUBLE PRECISION YY(100),WW(100)
      DOUBLE PRECISION ROT_E(3,3),ROT_P(3,3),ROTI_E(3,3),ROTI_P(3,3)
      DOUBLE PRECISION TRHIST_DEN(2,50,500),TRHIST_DEN_2D(2,50,50,50)
      DOUBLE PRECISION TMAT(6,6),TMATINV(6,6)
      DOUBLE PRECISION ARR(500),ARR_2D(50,50),FWHM_SIG
      DOUBLE PRECISION TARG_CUM(10)
      DOUBLE PRECISION QQ_VEC_VTX(3),QQ_VTX,THETAE_VTX,THETAP_VTX
      DOUBLE PRECISION ENEW,EPRIMENEW,POPNEW,EP_VTX,EP_ASM,BETA_H_ASM
      DOUBLE PRECISION JACOB_SA
      DOUBLE PRECISION DELTA_X,DELTA_Y,FP_VEC_E(6),FP_VEC_P(6)
      DOUBLE PRECISION TRPT_T_E_SAV(6),TRPT_T_P_SAV(6)
C
      REAL RAN_UNIF,RAN_UNIF2
C
      INTEGER LUN(50),I_WC,IMAT
      INTEGER IDUM1(6),IDUM2(6),I_VAR_TMP(100),I_WT_TMP(100)
      INTEGER LUN_SCT(2,50)
      INTEGER I_VAR_TR(6),IPID_TR_DEN(2,50)
      INTEGER LUN_NTU(20),NUM_VAR_NTU(20)
      INTEGER LUN_NTM(20),NUM_VAR_NTM(20),NUM_WT_NTM(20)
      INTEGER IARM,IEL,NCHAR,I_PF_ROOT_TMP(2)
      INTEGER EOPT,HOPT,N_DAT_DIR
      INTEGER ICNT,IND_ACCEPT_FCN(26)
      INTEGER LAGET_PS_FAIL,LAGET_GRID_FAIL
      LOGICAL SINGLES,FAIL_KIN,PROTON,ANSWER,IN_APERTURE,FIRST_TIME
      LOGICAL ELASTIC,ACCEPT_CHECK,BOUND
      LOGICAL RADPEAKING,MULTIPHOTON,RADPEAKING_SAV,RADFULL,ACCEPT_FCN
      LOGICAL FAIL_RAD(2),FAIL_CELL
      LOGICAL WRITE_TMP_NTU(20)
C
      CHARACTER*512 DAT_DIR_TMP
      CHARACTER*100 DAT_DIR
      CHARACTER*300  TMP_FILE
      CHARACTER*1   CDUM
C
C ---------------------------------------------------------------------
C       Include variable declarations and common blocks for:
C             - Input variables
C             - Histogram variables
C             - Spectrometer analysis variables
C             - Summary file statistics, etc.
C             - Laboratory variables (equivalences)
C             - Masses
C             - Wire chambers
C             - Energy loss routines
C             - HRS aperture histogramming
C ---------------------------------------------------------------------
C
      INCLUDE 'input.cmn'
      INCLUDE 'var.cmn'
      INCLUDE 'spectrometer.cmn'
      INCLUDE 'summary.cmn'
      INCLUDE 'labcoord.cmn'
      INCLUDE 'masses.cmn'
      INCLUDE 'wc.cmn'
      INCLUDE 'eloss.cmn'
      INCLUDE 'hrs.cmn'
      INCLUDE 'sep.cmn'
      INCLUDE 'mad.cmn'
      INCLUDE 'rada.cmn'
C
C ---------------------------------------------------------------------
C       Particle ID equivalences for cross sections.
C ---------------------------------------------------------------------
C
      EQUIVALENCE (SIGMA(0),PHASE_SP)
      EQUIVALENCE (SIGMA(1),SIGMA_EEP)
      EQUIVALENCE (SIGMA(2),POL_N)
      EQUIVALENCE (SIGMA(3),POL_T)
      EQUIVALENCE (SIGMA(4),POL_L)              
      EQUIVALENCE (SIGMA(5),POL_SPEC(1)) !N comp. of spec. pol. (i.e.
      EQUIVALENCE (SIGMA(6),POL_SPEC(2)) !T comp. of spec. pol.   after
      EQUIVALENCE (SIGMA(7),POL_SPEC(3)) !L comp. of spec. pol.    precess) 
      EQUIVALENCE (SIGMA(8), POL_X)
      EQUIVALENCE (SIGMA(9), POL_Y)
      EQUIVALENCE (SIGMA(10),POL_Z)
      EQUIVALENCE (SIGMA(11),ASYMMETRY)   ! e- analyzing power
C
C     Helicity independent polarization components
C
      EQUIVALENCE (SIGMA(12),POL_N_RHI)
      EQUIVALENCE (SIGMA(13),POL_T_RHI)
      EQUIVALENCE (SIGMA(14),POL_L_RHI)              
      EQUIVALENCE (SIGMA(15),POL_SPEC_HI(1)) !N comp. of spec. pol. (i.e.
      EQUIVALENCE (SIGMA(16),POL_SPEC_HI(2)) !T comp. of spec. pol.   after
      EQUIVALENCE (SIGMA(17),POL_SPEC_HI(3)) !L comp. of spec. pol.  precess) 
      EQUIVALENCE (SIGMA(18),POL_X_HI)
      EQUIVALENCE (SIGMA(19),POL_Y_HI)
      EQUIVALENCE (SIGMA(20),POL_Z_HI) 
C
C     Helicity dependent polarization components
C
      EQUIVALENCE (SIGMA(22),POL_N_RHD)
      EQUIVALENCE (SIGMA(23),POL_T_RHD)
      EQUIVALENCE (SIGMA(24),POL_L_RHD)              
      EQUIVALENCE (SIGMA(25),POL_SPEC_HD(1)) !N comp. of spec. pol. (i.e.
      EQUIVALENCE (SIGMA(26),POL_SPEC_HD(2)) !T comp. of spec. pol.   after
      EQUIVALENCE (SIGMA(27),POL_SPEC_HD(3)) !L comp. of spec. pol.  precess) 
      EQUIVALENCE (SIGMA(28),POL_X_HD)
      EQUIVALENCE (SIGMA(29),POL_Y_HD)
      EQUIVALENCE (SIGMA(30),POL_Z_HD)
C
      EQUIVALENCE (SIGMA(50),ACC_E_P)    !accidentals for E,P 
      EQUIVALENCE (SIGMA(51),ACC_E_PIP)
      EQUIVALENCE (SIGMA(52),ACC_PIM_P)
      EQUIVALENCE (SIGMA(53),ACC_PIM_PIP)
      EQUIVALENCE (SIGMA(100),SIGMA_EE)
      EQUIVALENCE (SIGMA(101),SIGMA_EPIM)
      EQUIVALENCE (SIGMA(102),SIGMA_EP)
      EQUIVALENCE (SIGMA(103),SIGMA_EPIP)
C
      EQUIVALENCE (SIGMA(121),SIGMA_MULT_WT(1))  !For multiple weight
      EQUIVALENCE (SIGMA(122),SIGMA_MULT_WT(2))   !N-tuples ('NTM').
      EQUIVALENCE (SIGMA(123),SIGMA_MULT_WT(3))
      EQUIVALENCE (SIGMA(124),SIGMA_MULT_WT(4))   !Multiple calcs for
      EQUIVALENCE (SIGMA(125),SIGMA_MULT_WT(5))   !(e,e'p) can be stuffed
      EQUIVALENCE (SIGMA(126),SIGMA_MULT_WT(6))   !into these indices.
      EQUIVALENCE (SIGMA(127),SIGMA_MULT_WT(7))
      EQUIVALENCE (SIGMA(128),SIGMA_MULT_WT(8))
      EQUIVALENCE (SIGMA(129),SIGMA_MULT_WT(9))
      EQUIVALENCE (SIGMA(130),SIGMA_MULT_WT(10))
C
      EQUIVALENCE (POL(1),POL_N_S)     !Spect. pol. vector at target.
      EQUIVALENCE (POL(2),POL_T_S)      !(i.e. before precession)
      EQUIVALENCE (POL(3),POL_L_S)
C
C     Helicity independent polarization components
C
      EQUIVALENCE (POL_SHI(1),POL_N_SHI)     !Spect. pol. vector at target.
      EQUIVALENCE (POL_SHI(2),POL_T_SHI)      !(i.e. before precession)
      EQUIVALENCE (POL_SHI(3),POL_L_SHI)
C
C     Helicity dependent polarization components
C
      EQUIVALENCE (POL_SHD(1),POL_N_SHD)     !Spect. pol. vector at target.
      EQUIVALENCE (POL_SHD(2),POL_T_SHD)      !(i.e. before precession)
      EQUIVALENCE (POL_SHD(3),POL_L_SHD)
C
      PARAMETER (MP = 938.2796D0)             !Proton  mass in MeV
      PARAMETER (MN = 939.5731D0)             !Neutron mass in MeV
      PARAMETER (Q_E = 1.602D0)               !e- charge in uCoul (*1E+13)
      PARAMETER (N_A = 6.022D0)               !Avogadro's number (*1E-23)
C
      INCLUDE 'var_dat.cmn'
      INCLUDE 'spectrometer_dat.cmn'
      INCLUDE 'masses_dat.cmn'
C
C     Cutoff energy for radiation is now read in.
C
C     DATA CUTOFF /1.0D0/                     !Photon cutoff energy (MeV)
C
      PI = ASIN(1.D0)*2.D0
      FWHM_SIG = 2.*SQRT(2.*LOG(2.))          !FWHM --> gaussian sigma
C
C ---------------------------------------------------------------------
C       Decode value of $mceep_dat.
C ---------------------------------------------------------------------
C
      CALL GETENV('mceep_dat',DAT_DIR_TMP)
      CALL SQUEEZE(DAT_DIR_TMP,DAT_DIR_TMP,N_DAT_DIR)
      DAT_DIR = DAT_DIR_TMP(1:N_DAT_DIR)
CXXX  DAT_DIR = '../dat/'
C
C ---------------------------------------------------------------------
C       Get name of input file.
C ---------------------------------------------------------------------
C
      WRITE(6,10)
   10 FORMAT(/' Enter input file prefix (.inp) >')
      READ(5,15)FILEPRE
   15 FORMAT(A20)
      TMP_FILE = FILEPRE//'.inp'
      CALL SQUEEZE(TMP_FILE,TMP_FILE,NCHAR)
      INFILE = TMP_FILE(1:NCHAR)
C
C ---------------------------------------------------------------------
C       Read kinematics, acceptances, etc. from INFILE.
C ---------------------------------------------------------------------
C
      CALL INPUT
C
C ---------------------------------------------------------------------
C       Convert quantities.
C ---------------------------------------------------------------------
C
      DTH_E = DTH_E/1000.D0           !Convert mrad -> radians.
      DTH_P = DTH_P/1000.D0
      DPH_E = DPH_E/1000.D0
      DPH_P = DPH_P/1000.D0
C
      DTH_E_N = DTH_E_N/1000.D0       !Convert mrad -> radians.
      DTH_P_N = DTH_P_N/1000.D0
      DPH_E_N = DPH_E_N/1000.D0
      DPH_P_N = DPH_P_N/1000.D0
C
      TH_P_MIN = -DTH_P/2.D0          !Assume symmetric.
      TH_E_MIN = -DTH_E/2.D0
      PH_P_MIN = -DPH_P/2.D0
      PH_E_MIN = -DPH_E/2.D0
C
      TH_E = TH_E*PI/180.D0           !Convert deg. -> radians.
      TH_P = TH_P*PI/180.D0
      TH_B = TH_B*PI/180.D0
C
      PH_E = PH_E*PI/180.D0           !Convert deg. -> radians.
      PH_P = PH_P*PI/180.D0
      PH_B = PH_B*PI/180.D0
C
      PHI_TARG = PHI_TARG*PI/180.D0
C
      BEAMV = 2.0*BEAMV/1000.         !convert to meters (for dispersed beam)
      BEAMV = BEAMV/COS(TH_B)         !Eff. vert. extent for bent beam
      BEAMD = BEAMD/100.              !from % to fraction
      E0_MIN = E0*(1.-BEAMD)          !for
      E0_MAX = E0*(1.+BEAMD)            !dispersed
      DE0 = E0_MAX-E0_MIN                 !beam
C
      FWHM_X_B = 0.01*FWHM_X_B        !Convert cm -> meters
      FWHM_Y_B = 0.01*FWHM_Y_B
      DEL_X_B  = 0.01*DEL_X_B
      DEL_Y_B  = 0.01*DEL_Y_B
C
      FWHM_PH_B = 0.001*FWHM_PH_B     !Convert mr -> radians
      FWHM_TH_B = 0.001*FWHM_TH_B
      DEL_PH_B  = 0.001*DEL_PH_B
      DEL_TH_B  = 0.001*DEL_TH_B
C
      GSIG_X_B  = FWHM_X_B /FWHM_SIG  !Convert FWHM -> gaussian sigma
      GSIG_Y_B  = FWHM_Y_B /FWHM_SIG
      GSIG_PH_B = FWHM_PH_B/FWHM_SIG
      GSIG_TH_B = FWHM_TH_B/FWHM_SIG
      GSIG_E0   = FWHM_E0  /FWHM_SIG
C
C ---------------------------------------------------------------------
C       Get sines and cosines.
C ---------------------------------------------------------------------
C
      CTH_B = COS(TH_B)
      STH_B = SIN(TH_B)
      CPH_B = COS(PH_B)
      SPH_B = SIN(PH_B)
C
      CTH_E = COS(TH_E)
      STH_E = SIN(TH_E)
      CPH_E = COS(PH_E)
      SPH_E = SIN(PH_E)
C
      CTH_P = COS(TH_P)
      STH_P = SIN(TH_P)
      CPH_P = COS(PH_P)
      SPH_P = SIN(PH_P)
C
      CPHI_TARG = COS(PHI_TARG)
      TPHI_TARG = TAN(PHI_TARG)
      CPHIE_DIFF = COS(PH_E-PHI_TARG)
      CPHIP_DIFF = COS(PH_P-PHI_TARG)
C
C ---------------------------------------------------------------------
C       Ask questions and set logical flags.
C ---------------------------------------------------------------------
C
      CALL PHYS_CHOICE(ELASTIC,ACCEPT_CHECK,BOUND,PROTON,RADPEAKING,
     #                 MULTIPHOTON,RADFULL,ACCEPT_FCN)
      IF(ACCEPT_FCN) ELOSS_CALC = .FALSE.
      IF(ELASTIC) THEN
        IF(RADPEAKING) THEN
           CHAR_SIG_UNITS = 'fm^2/MeV/sr'
        ELSE
           CHAR_SIG_UNITS = 'fm^2/sr'
        ENDIF
      ELSEIF(BOUND) THEN
        IF(RADPEAKING) THEN
           CHAR_SIG_UNITS = 'fm^2/MeV^2/sr^2'
        ELSE
           CHAR_SIG_UNITS = 'fm^2/MeV/sr^2'
        ENDIF
      ELSE
        IF(RADPEAKING) THEN
           CHAR_SIG_UNITS = 'fm^2/MeV-(MeV/c-sr)^2'
        ELSE
           CHAR_SIG_UNITS = 'fm^2/(MeV/c-sr)^2'
        ENDIF
      ENDIF
C
      IF(RADFULL) THEN
         countall = 0    ! initialize
         countbad = 0    ! initialize
      ENDIF
C
      IF(.NOT. ELASTIC .AND. .NOT. ACCEPT_FCN) THEN
         SINGLES = ANSWER(' Compute singles and accidentals? (Y/N)>',
     #                    'Y','y')
         IF(SINGLES .AND. BOUND) WRITE(6,20)
   20    FORMAT(' Warning:  Singles being calculated for fixed Em ')
      ELSE
         SINGLES = .FALSE.
      ENDIF
C
C ---------------------------------------------------------------------
C       Retrieve initial seeds.  All randomly selected variables are
C       now based on the two starting seeds, rather than having
C       separate initial seeds for each variable.
C ---------------------------------------------------------------------
C
      CALL GETSEED
C
C ---------------------------------------------------------------------
C       Spectrometer analysis:  Read in Transport matrices and find
C       inverses where required.  Set up 1-D and 2-D histos and
C       scatter plots.  In addition, make logical array which defines
C       whether or not a given spectrometer element is to be applied.  This
C       array depends both on the value of OB (whether action is to affect
C       cross section ('S') or kinematics ('K') or both ('B')) and on
C       which call to SPECTROMETER (1st call - cross section;
C       2nd call - kinematic bins) is being made.
C ---------------------------------------------------------------------
C
      ICOUNT_SCT = 0      !Initialize # of TRANSPORT scatter plots
      ICOUNT_NTU = 0      !Initialize # of N-Tuple files (this counter
                            !runs over both TRANSPORT and reg. N-Tuples)
      ICOUNT_NTM = 0      !Initialize # of N-Tuple files for 'NTM' N-Tuples
      DO I=1,6
        I_VAR_TR(I) = I   !N-Tuple variables 1-6 written
      ENDDO                !(i.e. write full Transport vector)
      DO IARM=1,2
        IMAT  = 0          !reset matrix counter for each arm
        IF(NEL(IARM) .GT. 0) THEN
          DO IEL = 1,NEL(IARM)
            IF(OP(IARM,IEL).EQ.'RFN') THEN
              CALL INIT_R_FUNCTION(HRS_ID_INT(IARM,IEL),
     #                             HRS_COLL(HRS_ID_INT(IARM,IEL)))
            ELSEIF(OP(IARM,IEL).EQ.'MAT') THEN
              IMAT = IMAT + 1
              TMP_FILE=DAT_DIR(1:N_DAT_DIR)//'/'//TRNSPT_FIL(IARM,IEL)
              OPEN(UNIT=20,FILE=TMP_FILE,
     #               STATUS='OLD',FORM='FORMATTED')
              DO I=1,6
                READ(20,*)(TMATRIX(IARM,IMAT,I,J),J=1,6)
              ENDDO
              IF(INVERT(IARM,IEL)) THEN
                DO I=1,6
                  DO J=1,6
                    TMAT(I,J) = TMATRIX(IARM,IMAT,I,J)
                  ENDDO
                ENDDO
                CALL MATINV(TMAT,TMATINV,6)     !Get inverse matrix
                DO I=1,6
                  DO J=1,6
                    TMATRIX(IARM,IMAT,I,J) = TMATINV(I,J)
                  ENDDO
                 ENDDO
              ENDIF
              IF(ORDER(IARM,IEL) .EQ. 2) THEN !Read in the 2nd order mat.
                DO I=1,5
                  DO K=1,6
                    READ(20,*)(IDUM1(J),IDUM2(J),
     #                         TMATRIX2(IARM,IMAT,I,J,K),J=1,K)
                  ENDDO
                  READ(20,30) CDUM
30                FORMAT(A1)
                ENDDO
                IF(INVERT(IARM,IEL)) THEN
                  READ(20,30) CDUM
                  DO I=1,5
                    DO K=1,6
                      READ(20,*)(IDUM1(J),IDUM2(J),
     #                           TMATRIX2(IARM,IMAT,I,J,K),J=1,K)
                    ENDDO
                    READ(20,30) CDUM
                  ENDDO
                ENDIF
              ENDIF
              CLOSE(UNIT=20)
C
            ELSEIF(OP(IARM,IEL).EQ.'TRK') THEN
C
              TMP_FILE = DAT_DIR(1:N_DAT_DIR)//'/'//WC_FIL(IARM)
                            !File with W.C. prmtrs 
C
C ---------------------------------------------------------------------
C             The parameters read from WC_FIL are described in
C             get_vect.f.
C ---------------------------------------------------------------------
C
              OPEN(UNIT=20,FILE=TMP_FILE,
     #               STATUS='OLD',FORM='FORMATTED')
              READ(20,*) X_DELTA(IARM)
              READ(20,*) WC_OFF(IARM),WC_OFF_COORD(IARM),
     #                WC_OFF_CHMB(IARM)
              READ(20,*) WC_RES(IARM)
              READ(20,*) N_WC(IARM)
              WC_RES(IARM)  = WC_RES(IARM) /FWHM_SIG !convert to sigma
C
              DO I_WC = 1,N_WC(IARM)
                READ(20,*) WC_C_ANGLE(IARM,I_WC),WC_W_ANGLE(IARM,I_WC),
     #                     WC_LOC(IARM,I_WC),WC_RAD(IARM,I_WC)
C
                WC_C_CANGLE(IARM,I_WC) = COS(WC_C_ANGLE(IARM,I_WC)
     #                               *PI/180.D0)
                WC_C_SANGLE(IARM,I_WC) = SIN(WC_C_ANGLE(IARM,I_WC)
     #                               *PI/180.D0)
                WC_W_CANGLE(IARM,I_WC) = COS(WC_W_ANGLE(IARM,I_WC)
     #                               *PI/180.D0)
                WC_W_SANGLE(IARM,I_WC) = SIN(WC_W_ANGLE(IARM,I_WC)
     #                               *PI/180.D0)
                WC_MSCT(IARM,I_WC) = SIG_MSCAT(IARM,PF_E,PF_P,
     #                      WC_RAD(IARM,I_WC),WC_C_CANGLE(IARM,I_WC))
C
              ENDDO
              CLOSE(UNIT=20)
C
            ELSEIF(OP(IARM,IEL).EQ.'POL') THEN
              BEND_ANGLE(IARM,IEL) = BEND_ANGLE(IARM,IEL)*PI/180.D0
C
            ELSEIF(OP(IARM,IEL).EQ.'COS') THEN
              call read_cosy_input_file(cosy_spin_fil(iarm,iel),
     #              pf_p,dat_dir,n_dat_dir)
C
            ELSEIF(OP(IARM,IEL).EQ.'ROT') THEN
              ROT_ANGLE(IARM,IEL) = ROT_ANGLE(IARM,IEL)*PI/180.D0
              C_ROT_ANGLE(IARM,IEL) = COS(ROT_ANGLE(IARM,IEL))
              S_ROT_ANGLE(IARM,IEL) = SIN(ROT_ANGLE(IARM,IEL))
C
            ELSEIF(OP(IARM,IEL).EQ.'H1D') THEN
              CALL GET_PID_DEN(IPID_TR(IARM,IEL),IPID_TR_DEN(IARM,IEL))
C
            ELSEIF(OP(IARM,IEL).EQ.'H2D') THEN
              CALL GET_PID_DEN(IPID_TR(IARM,IEL),IPID_TR_DEN(IARM,IEL))
C
            ELSEIF(OP(IARM,IEL).EQ.'SCT' .AND. .NOT. ACCEPT_FCN) THEN
              ICOUNT_SCT = ICOUNT_SCT + 1
              LUN_SCT(IARM,IEL) = ICOUNT_SCT + 30
              CALL TOP_SCAT(LUN_SCT(IARM,IEL),TR_FILE(IARM,IEL),
     #          TR_AXIS(NCOORD(IARM,IEL,1)),TR_AXIS(NCOORD(IARM,IEL,2)),
     #          .FALSE.,TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL),
     #          TR_MINY(IARM,IEL),TR_MAXY(IARM,IEL))
C
            ELSEIF(OP(IARM,IEL).EQ.'NTU' .AND. .NOT. ACCEPT_FCN) THEN
              BEND_ANGLE_TOT(IARM,IEL) = BEND_ANGLE_TOT(IARM,IEL)
     #            *PI/180.D0
              ICOUNT_NTU = ICOUNT_NTU + 1
              LUN_NTU(ICOUNT_NTU) = ICOUNT_NTU + 50
              NUM_VAR_NTU(ICOUNT_NTU) = 6
              CALL NTU_SETUP(LUN_NTU(ICOUNT_NTU),TR_FILE(IARM,IEL),
     #            IPID_TR(IARM,IEL),.TRUE.,6,I_VAR_TR)
              IF(IPID_TR(IARM,IEL) .GE. 0) THEN
                  WRITE_TMP_NTU(ICOUNT_NTU) = .TRUE.
              ELSE
                  WRITE_TMP_NTU(ICOUNT_NTU) = .FALSE.
              ENDIF
              CALL GET_PID_DEN(IPID_TR(IARM,IEL),IPID_TR_DEN(IARM,IEL))
C
            ENDIF
C
            ACTIVE_EL(IARM,IEL,1) = .FALSE.      !Initialize
            ACTIVE_EL(IARM,IEL,2) = .FALSE.
            IF((OB(IARM,IEL).EQ.'S').OR.(OB(IARM,IEL).EQ.'s')
     #          .OR.(OB(IARM,IEL).EQ.'B').OR.(OB(IARM,IEL).EQ.'b'))
     #               ACTIVE_EL(IARM,IEL,1) = .TRUE.
            IF((OB(IARM,IEL).EQ.'K').OR.(OB(IARM,IEL).EQ.'k')
     #          .OR.(OB(IARM,IEL).EQ.'B').OR.(OB(IARM,IEL).EQ.'b'))
     #               ACTIVE_EL(IARM,IEL,2) = .TRUE.
            IF(ACCEPT_FCN) THEN
               ACTIVE_EL(IARM,IEL,1) = .TRUE.
            ENDIF
C
          ENDDO
        ENDIF
      ENDDO
C
      DO I=1,3          !Initialize spin rotation matrix
      DO J=1,3           !Only needed if POL or COS options are not used
         IF(I.EQ.J) THEN
            PMATRIX(I,J) = 1.D0
         ELSE
            PMATRIX(I,J) = 0.D0
         ENDIF
      ENDDO
      ENDDO
C
C ---------------------------------------------------------------------
C       Calculate rotation matrices from laboratory floor coordinates to
C       spectrometer (Transport) coordinates at target.  Conventions
C       are:
C   F
C    L          X - YxZ
C     O         Y - normal to laboratory floor (upwards = +)
C      O        Z - along nominal beam direction (i.e. for PH_B=TH_B=0)
C       R
C
C   S
C    P          X - along momentum dispersion (=-Y_floor for HRS2-CEBAF)
C     E         Y - ZxX
C      C        Z - along spectrometer central ray
C       T
C
C
C   Note:  THETA_BP(1) and THETA_BP(2) are the spectrometer bend plane
C          angles for the electron and hadron arm respectively.
C          For a spectrometer which bends vertically upward, the
C          bend plane angle is -90 degrees.  The corresponding rotation
C          below then insures that the Transport X axis points
C          vertically downward at the target.  MCEEP can accomodate a
C          spectrometer whose bend plane is at any arbitrary
C          orientation.  The common orientations are (as one looks
C          down the central ray into the spectrometer from the target):
C
C          BENDING ORIENTATION                THETA_BP (deg)
C
C           vertically upward                    -90
C           vertically downward                   90
C           leftward                             180
C           rightward                              0
C
C
C     The Y and X (passive) rotations below rotate the FLOOR coordinate
C     system so that its Z axis is aligned with the spectrometer
C     central ray.  The Z rotation then accounts for the orientation
C     of the spectrometer bend plane.
C
C ---------------------------------------------------------------------
C
      THETA_BP(1) = THETA_BP(1)*PI/180.D0       !convert to radians
      CALL ROTATE_Y(CPH_E, SPH_E,ROTY)
      CALL ROTATE_X(CTH_E,-STH_E,ROTX)
      CALL ROTATE_Z(COS(THETA_BP(1)),SIN(THETA_BP(1)),ROTZ)
      CALL MAT_MULT(3,3,3,ROTX,ROTY,ROTXY)
      CALL MAT_MULT(3,3,3,ROTZ,ROTXY,ROT_E)
C
      THETA_BP(2) = THETA_BP(2)*PI/180.D0       !convert to radians
      CALL ROTATE_Y(CPH_P, SPH_P,ROTY)
      CALL ROTATE_X(CTH_P,-STH_P,ROTX)
      CALL ROTATE_Z(COS(THETA_BP(2)),SIN(THETA_BP(2)),ROTZ)
      CALL MAT_MULT(3,3,3,ROTX,ROTY,ROTXY)
      CALL MAT_MULT(3,3,3,ROTZ,ROTXY,ROT_P)
C
C ---------------------------------------------------------------------
C       Calculate inverse matrices.
C ---------------------------------------------------------------------
C
      DO I=1,3
        DO J=1,3
          ROTI_E(I,J) = ROT_E(J,I)
          ROTI_P(I,J) = ROT_P(J,I)
        ENDDO
      ENDDO
C
C ---------------------------------------------------------------------
C       Target Setup.
C ---------------------------------------------------------------------
C
      TARG_TOT = 0.D0   !Initialize total target length
      DO I=1,NFOILS
         TARG_LO(I) = TARG_LO(I)/CTH_B   !Eff. targ. positions for bent beam
         TARG_HI(I) = TARG_HI(I)/CTH_B
         TARG_TOT = TARG_TOT + (TARG_HI(I)-TARG_LO(I))   ! in meters
         TARG_CUM(I) = TARG_TOT          !Cumulative thickness up to foil I
      ENDDO
C
      NTARG = ATARG - ZTARG           !Number of neutrons for target nucleus
C
      IF(ATARG .EQ. 1.) THEN
        MNUC = 938.279D0        !1H
      ELSEIF(ATARG .EQ. 2.) THEN
        MNUC = 1875.628D0       !2H
      ELSEIF(ATARG .EQ. 3. .AND. ZTARG .EQ. 1.) THEN
        MNUC = 2808.943D0       !3H
      ELSEIF(ATARG .EQ. 3. .AND. ZTARG .EQ. 2.) THEN
        MNUC = 2808.413D0       !3He
      ELSEIF(ATARG .EQ. 4. .AND. ZTARG .EQ. 2.) THEN
        MNUC = 3727.403D0       !4He
      ELSEIF(ATARG .EQ. 12 .AND. ZTARG .EQ. 6.) THEN
        MNUC = 11177.95D0       !12C
      ELSE
        MNUC = ATARG*(931.494D0-0.511D0)
      ENDIF
C
C ---------------------------------------------------------------------
C       Pack the common block /TG1/.
C
C          For the LH2/LD2 or He cryotargets, targ_lo and targ_hi
C          specify the locations of the cell entrance and exit walls.
C          This information is used to calculate the path length
C          traversed by each particle in the routine TG_PATH.
C
C          For multiple foils, we just pass in the starting
C          and ending positions of the first foil as "dummies".
C          For foils the path lengths are calculated differently,
C          so this poses no problem.
C ---------------------------------------------------------------------
C
      CALL TARG_SETUP(ZTARG,ATARG,DENS_TARG,
     #                TARG_LO(1),TARG_HI(1),ITARG_MOD)
C
C ---------------------------------------------------------------------
C       Energy loss calculation for target walls/windows.
C
C            For the LH2/LD2 cryotarget beer can (ITARG_MOD=2),
C            for the He  cryotarget (ITARG_MOD=3),
C            for the Polarized He3 target (ITARG_MOD=4),
C            and for the LH2/LD2 cryotarget cigar tube (ITARG_MOD=5),
C            the cell walls are included.
C
C            For foil(s) (ITARG_MOD=1), these contributions
C            are set to zero.
C ---------------------------------------------------------------------
C
      IF(ELOSS_CALC) THEN
        IF(ITARG_MOD .EQ. 1) THEN
           EOPT = 0     ! These are relevant only for the cryotarget
           HOPT = 0      ! and decide whether the particle passes
                          ! through an endcap or a side wall.
        ENDIF
        CALL ELOSSW(E0,PF_E,PH_E,PF_P,PH_P,ITARG_MOD) ! walls/windows
C
C       As of Version 3.9:  for foils, include air/mylar only in calc_b 
C
        CALL CALC_B(PH_E,ITARG_MOD)               ! 1/rad-length (cm^2/g)
        CALL DENSCORR(E0,M_ELEC,1)                ! beam
        CALL DENSCORR(PF_E,M_ELEC,2)              ! scatt. electron
        IF(EJECT_CHG .NE. 0) CALL DENSCORR(PF_P,EJECT_MASS,3) ! ejectile
      ENDIF
C
C ---------------------------------------------------------------------
C       Calculate solid angle for each arm (in sr).
C ---------------------------------------------------------------------
C
      IF(SA_SHAPE_E .EQ. 'E') THEN
        DOMEG_E = DTH_E*DPH_E*PI/4.D0
      ELSE
        DOMEG_E = DTH_E*DPH_E
      ENDIF
C
      IF(SA_SHAPE_P .EQ. 'E') THEN
        DOMEG_P = DTH_P*DPH_P*PI/4.D0
      ELSE
        DOMEG_P = DTH_P*DPH_P
      ENDIF
C
C ---------------------------------------------------------------------
C       Calculate bounds on electron and hadron momenta.
C ---------------------------------------------------------------------
C
      PF_E_MIN = PF_E*(1.+0.01D0*ACC_EM)
      PF_E_MAX = PF_E*(1.+0.01D0*ACC_EP)
      DPF_E = PF_E_MAX - PF_E_MIN
C
      PF_P_MIN = PF_P*(1.+0.01D0*ACC_PM)
      PF_P_MAX = PF_P*(1.+0.01D0*ACC_PP)
      DPF_P = PF_P_MAX - PF_P_MIN
C
C ---------------------------------------------------------------------
C     For radiative tail:
C
C       Determine maximum photon energy.  For now it is chosen to be
C       the minimum scattered electron energy, since this is
C       guaranteed to be consistent with event kinematics.
C       This excludes the very high energy part of the tail.
C
C       Also normalize the integral.  Since the photon energies
C       are sampled according to 1/k, the proper normalization
C       is log(eg_max/cutoff).
C
C       Finally, perform other radiation specific setup.
C ---------------------------------------------------------------------
C
      IF(RADPEAKING) THEN
         RADPEAKING_SAV = .TRUE.
         NCALC = 2    ! Do event loop twice (for tail and for peak)
         EG_MAX = PF_E_MIN
         RAD_NORM(1) = LOG(EG_MAX/CUTOFF)
         RAD_NORM(2) = 1.D0
      ELSE
         RADPEAKING_SAV = .FALSE.
         NCALC = 1    ! Do event loop only once since no tail calc.
         RAD_NORM(1) = 1.D0
      ENDIF
C
C ---------------------------------------------------------------------
C       Determine factor to convert cross sections to yield for each
C       particle ID.  Single arm cross sections are assumed to be in
C       nb/(MeV/c-sr).  Coincidence cross sections (SIGMA_EEP) returned
C       to MCEEP by the routine PHYSICS are in:
C                fm^2 sr^-1                   (elastic scattering)
C                fm^2 sr^-2 MeV^-1            (bound state)
C                fm^2 sr^-2 MeV^-1 (MeV/c)^-1 (continuum)
C       Note that the continuum cross section should be differential
C       in the hadron momentum (NOT kinetic energy) since it is the
C       momentum which is randomly sampled.
C
C       A NOTE about normalizations:
C         Y = Luminosity x (Acceptance/Number_into_accept) x Sum_of_sigma
C
C       where Y = Yield and Acceptance is the product of solid
C       angles and momentum/energy acceptances.  This expression is
C       valid since the second term in the product corresponds to the
C       average phase space volume per event.  Thus the above gives
C       the luminosity times the sum over events of [the cross section times
C       the phase space volume per event].
C
C       Number_into_accept is factored in later since, for elliptical
C       apertures, this number is not known ahead of time.
C ---------------------------------------------------------------------
C
      BEAMTIME = BEAMTIME * 3600.D0           !convert hours to seconds
      DUTY_FAC = DUTY_FAC * 0.01D0            !convert % to decimal
      DEL_TOF  = DEL_TOF  * 1.D-9             !convert nsec to sec
C
      SIGFACT = ALUM*(N_A/(Q_E*ATARG))*BEAMTIME !per picobarn
      SIGFACTOR(0) = SIGFACT*DOMEG_E*DOMEG_P*DPF_E*DPF_P
     #               * 1.D10
      IF(ELASTIC) THEN
        SIGFACTOR(0) = SIGFACTOR(0)/(DOMEG_P*DPF_E*DPF_P)
        SIGFACTOR(1) = SIGFACTOR(0)
      ELSEIF(BOUND) THEN
        SIGFACTOR(0) = SIGFACTOR(0)/DPF_P
        SIGFACTOR(1) = SIGFACTOR(0)
      ELSE
        SIGFACTOR(1) = SIGFACTOR(0)
      ENDIF
      DO I=2,11        !For polarizations in (e,e'p) or e- anal. power
        SIGFACTOR(I) = SIGFACTOR(1)
      ENDDO
      DO I=12,20       !For polarizations in (e,e'p)
        SIGFACTOR(I) = SIGFACTOR(1)
      ENDDO
      DO I=22,30       !For polarizations in (e,e'p)
        SIGFACTOR(I) = SIGFACTOR(1)
      ENDDO
      SIGFACTOR(100) = SIGFACT*DOMEG_E*DPF_E*1.D3
      SIGFACTOR(101) = SIGFACTOR(100)
      SIGFACTOR(102) = SIGFACT*DOMEG_P*DPF_P*1.D3
      SIGFACTOR(103) = SIGFACTOR(102)
      DO I=121,130
         SIGFACTOR(I) = SIGFACTOR(1)    ! (e,e'p) multiple calcs
      ENDDO
C
      IF(SINGLES) THEN
        IF(DUTY_FAC .EQ. 0.D0) THEN
           WRITE(6,500)
500        FORMAT(' Warning:  Duty Factor = 0, Setting to 100% ')
           DUTY_FAC = 1.D0
        ENDIF
        SIGFACTOR(50) = SIGFACTOR(100)*SIGFACTOR(102)*
     #                  (DEL_TOF/(DUTY_FAC*BEAMTIME))
        SIGFACTOR(51) = SIGFACTOR(100)*SIGFACTOR(103)*
     #                  (DEL_TOF/(DUTY_FAC*BEAMTIME))
        SIGFACTOR(52) = SIGFACTOR(101)*SIGFACTOR(102)*
     #                  (DEL_TOF/(DUTY_FAC*BEAMTIME))
        SIGFACTOR(53) = SIGFACTOR(101)*SIGFACTOR(103)*
     #                  (DEL_TOF/(DUTY_FAC*BEAMTIME))
      ENDIF
C 
C ---------------------------------------------------------------------
C       Initialize histograms and statistics.
C ---------------------------------------------------------------------
C
      DO IARM=1,2                     !TRANSPORT histograms
        IF(NEL(IARM) .GT. 0) THEN
          DO IEL = 1,NEL(IARM)
            IF(OP(IARM,IEL).EQ.'H1D') THEN
              I_OOB_TR(IARM,IEL) = 0
              DO J=1,NXCHAN_TR(IARM,IEL)
                TRHIST(IARM,IEL,J) = 0.
                TRHIST_DEN(IARM,IEL,J) = 0.
              ENDDO
            ELSEIF(OP(IARM,IEL).EQ.'H2D') THEN
              I_OOB_TR(IARM,IEL) = 0
              DO J=1,NXCHAN_TR(IARM,IEL)
              DO K=1,NYCHAN_TR(IARM,IEL)
                TRHIST_2D(IARM,IEL,J,K) = 0.
                TRHIST_DEN_2D(IARM,IEL,J,K) = 0.
              ENDDO
              ENDDO
            ELSEIF(OP(IARM,IEL).EQ.'SCT') THEN
              I_OOB_TR(IARM,IEL) = 0
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
      IF(NPLOTS .GT. 0) THEN
        DO I=1,NPLOTS
          IF(PLOT_TYPE(I).EQ.'P1D') THEN
            I_OOB(I) = 0
            DO J=1,NX_CHAN(I)
              SIG_ARR(I,J) = 0.
              SIG_DEN(I,J) = 0.
            ENDDO
          ELSEIF(PLOT_TYPE(I).EQ.'P2D') THEN
            I_OOB(I) = 0
            DO J=1,NX_CHAN(I)
            DO K=1,NY_CHAN(I)
              SIG_ARR_2D(I,J,K) = 0.
              SIG_DEN_2D(I,J,K) = 0.
            ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
      DO I=1,NUM_VAR
        VAR_CENT(I) = 0.
        VAR_AVG(I)  = 0.
        VAR_MIN(I)  = 1.D32
        VAR_MAX(I)  = -1.D32
      ENDDO
C
      I_AREN_FAIL     = 0   !Failed Arenhoevel interpolation
      I_EXTRAP        = 0   !Extrapolation statistic for interpolator
      I_INTERPEXP_NP  = 0   !Exponential interpolator failed due to 
                            ! non-positive y-value
      I_SPECT_OOR     = 0   !Spectral function out-of-range statistic
      LAGET_PS_FAIL   = 0   !Failed Laget phase-space
      LAGET_GRID_FAIL = 0   !Failed Laget grid boundary
C
      SUM_ACC   = 0.D0  !Number of tries into modified acceptance
      SUM_EV    = 0.D0  !  into actual accept. with physical solutions
      SUM_EV_C  = 0.D0  !  also passing global cuts and HRS apertures
      SUM_PS    = 0.D0  !Sum of phase space
      SUM_PS_C  = 0.D0  !Sum of phase space after global cuts
      SUM_EEP   = 0.D0  !Sum of (e,e'p) cross section
      SUM_EEP_C = 0.D0  !Sum of (e,e'p) cross section after global cuts
C
C ---------------------------------------------------------------------
C       Get default plot ranges, scales and offsets.
C ---------------------------------------------------------------------
C
      IF(BOUND) MISS_M = EM_BOUND             !Set missing mass
      IF(.NOT. ACCEPT_FCN) CALL RANGES(ELASTIC,BOUND)
C
C ---------------------------------------------------------------------
C       Open scatter plot and N-Tuple files.
C ---------------------------------------------------------------------
C
      IF(NPLOTS .GT. 0 .AND. .NOT. ACCEPT_FCN) THEN
        DO I=1,NPLOTS
          IF(PLOT_TYPE(I) .EQ. 'SCA') THEN
            LUN(I) = 40 + I
            CALL TOP_SCAT(LUN(I),PLOT_FIL(I),VAR_NAME(I_VAR(I,1)),
     #         VAR_NAME(I_VAR(I,2)),.FALSE.,X_MIN(I),X_MAX(I),
     #         Y_MIN(I),Y_MAX(I))
          ELSEIF(PLOT_TYPE(I) .EQ. 'NTU') THEN
            ICOUNT_NTU = ICOUNT_NTU + 1
            DO J=1,NVAR_NTU(I)
               I_VAR_TMP(J) = I_VAR(I,J)   !Repack into 1-D array
            ENDDO
            LUN_NTU(ICOUNT_NTU) = ICOUNT_NTU + 50
            NUM_VAR_NTU(ICOUNT_NTU) = NVAR_NTU(I)
            CALL NTU_SETUP(LUN_NTU(ICOUNT_NTU),PLOT_FIL(I),IPID(I),
     #         .FALSE.,NVAR_NTU(I),I_VAR_TMP)
            IF(IPID(I) .GE. 0) THEN
                WRITE_TMP_NTU(ICOUNT_NTU) = .TRUE.
            ELSE
                WRITE_TMP_NTU(ICOUNT_NTU) = .FALSE.
            ENDIF
          ELSEIF(PLOT_TYPE(I) .EQ. 'NTM') THEN
            ICOUNT_NTM = ICOUNT_NTM + 1
            DO J=1,NWT_NTM(I)
               I_WT_TMP(J)  = I_WT(I,J)    !Repack into 1-D array
            ENDDO
            DO J=1,NVAR_NTU(I)
               I_VAR_TMP(J) = I_VAR(I,J)   !Repack into 1-D array
            ENDDO
            LUN_NTM(ICOUNT_NTM)     = ICOUNT_NTM + 70
            NUM_WT_NTM(ICOUNT_NTM)  = NWT_NTM(I)
            NUM_VAR_NTM(ICOUNT_NTM) = NVAR_NTU(I)
            CALL NTM_SETUP(LUN_NTM(ICOUNT_NTM),PLOT_FIL(I),NWT_NTM(I),
     #         I_WT_TMP,NVAR_NTU(I),I_VAR_TMP)
          ENDIF
        ENDDO
      ENDIF
C
      NCOUNT_NTU = ICOUNT_NTU
      NCOUNT_NTM = ICOUNT_NTM
C
      IF(NCOUNT_NTU .GT. 10) THEN
         WRITE(6,*) ' Error: Total number of NTU Ntuples must be < 11 '
         STOP
      ENDIF
C
      IF(NCOUNT_NTM .GT. 10) THEN
         WRITE(6,*) ' Error: Total number of NTM tuples must be < 11 '
         STOP
      ENDIF
C
C ---------------------------------------------------------------------
C       For the calculation of the acceptance function, open
C       the sole Ntuple here.
C ---------------------------------------------------------------------
C
      IF(ACCEPT_FCN) THEN
         OPEN(UNIT=41,FILE='accept_fcn.ntu',STATUS='UNKNOWN',
     #        FORM='FORMATTED')
         WRITE(41,106) ' KINEM '
 106     FORMAT(A)
         IWEIGHT = 9999
         NVAR    = 26
         WRITE(41,107) IWEIGHT,NVAR
 107     FORMAT(1X,I5,1X,I5)
         CALL INDEX_ACCEPT_FCN(NVAR,IND_ACCEPT_FCN)
         WRITE(41,108) (IND_ACCEPT_FCN(I),I=1,NVAR)
 108     FORMAT(1X,100(I3,1X))
      ENDIF
C
      IF(BOUND) MISS_M = EM_BOUND             !Reset missing mass
C
C _____________________________________________________________________
C ---------------------------------------------------------------------
C       Start the event loop.
C
C       For RADPEAKING=.TRUE. the event loop is executed twice, first
C       for the tail calculation and second for the peak calculation.
C
C       Angles are chosen to be within collimator angular ranges.
C       This gives the nominal angles (labeled below by _N).  For
C       extended targets this defines the (X,Y,Z) position within
C       the collimator.  Given the positional coordinates of the
C       beam-target interaction point, one can then determine the
C       actual angles relative to the laboratory floor coordinate
C       system.  The interaction point thus defines the origin of
C       the coordinate system for subsequent kinematical calcuations.
C
C       All vector quantities are expressed in the so-called
C       "laboratory" coordinate system.  This system is the nominal
C       beam system (i.e. for beam bend angle = 0 ; this corresponds
C       to the Hall floor).  Axes are:
C
C               X - given by YxZ
C               Y - perpendicular to floor (upwards = +)
C               Z - nominal beam direction.
C
C       Positive beam bend angles imply an upward direction for the beam.
C
C       A NOTE ABOUT THE SAMPLING METHOD:
C
C       The so-called "collimator" angles are sampled.  These are
C       related to the "initial electron TRANSPORT coordinate system"
C       (see above), also refered to as LAB system, angles as:
C              ThetaC =  ThetaL - ThetaL_0
C              PhiC   = (PhiL   - PhiL_0) x cos(ThetaL)
C       where ThetaL_0 and PhiL_0 are the central spectrometer angles
C       in the lab.  In MCEEP's LAB system the y-axis is the polar axis.
C       The polar angle is 90-ThetaL.  Thus, the uniform sampling of
C       the collimator angles (dPhiC, dThetaC) results in a uniform
C       sampling of the quantity sin(90-ThetaL)dThetaL dPhiL,
C       which is just the solid angle in spherical coordinates.
C       Since the cross sections are assumed to be differential in
C       spherical solid angles this gives the required uniform
C       sampling of the solid angles.
C _____________________________________________________________________
C
      DO ICALC=1,NCALC    ! 1=rad-tail / 2="peak"
      IF(RADPEAKING .AND. ICALC .EQ. 1) THEN
         WRITE(6,'(a)') ' -------------------------------- '
         WRITE(6,'(a)') ' Doing Radiation Tail Calculation '
         WRITE(6,'(a)') ' -------------------------------- '
      ELSEIF(ICALC .EQ. 2) THEN
         WRITE(6,'(a)') ' ---------------------- '
         WRITE(6,'(a)') ' Doing Peak Calculation '
         WRITE(6,'(a)') ' ---------------------- '
      ENDIF
C
      DO IEVNT=1,N_EVENT
C
        IF(DFLOAT(IEVNT/10000)-DFLOAT(IEVNT)/10000. .EQ. 0.)
     #          WRITE(6,101)IEVNT
  101   FORMAT(' Event number: ',I9)
        FAIL_KIN   = .FALSE.
        FIRST_TIME = .TRUE.        !TRUE for first call to KINEM
C
C ---------------------------------------------------------------------
C          Determine beam configuration first.  This allows the
C          determination of the beam-target interaction point which
C          is needed for extracting the "actual" angles from the
C          "nominal" angles and vice versa.  "Nominal" angles are defined
C          by the intersection of the particle ray with the aperture
C          relative to the nominal target center (this is actually the
C          "Hall origin").  "Actual" angles are defined by the same
C          point relative to the beam-target interaction point
C          (actual vertex).
C
C          Get nominal Z position along target first.  For extended
C          targets, this is the Z position of the vertex.  However,
C          for foils with non-zero angle coupled with a beam X offset,
C          this will NOT be exactly the vertex Z.
C
C          Note that the beam X offset affect on the target length
C          for ITARG_MOD=3 (tuna can) is neglected here - it's a
C          cosine effect and therefore small.
C ---------------------------------------------------------------------
C
        CALL RANECU(RAN_UNIF,1)
        TARG0 = TARG_CUM(NFOILS)*DBLE(RAN_UNIF) ! targ. material
                                                 ! before interaction
        DO J=1,NFOILS
           IF(TARG0 .LE. TARG_CUM(J)) THEN
              IFOIL = J        ! Got the foil with the vertex
              GOTO 111
           ENDIF
        ENDDO
 111    IF(IFOIL .EQ. 1) THEN
           TARG_I = TARG_LO(IFOIL) + TARG0
        ELSE
           TARG_I = TARG_LO(IFOIL) + (TARG0-TARG_CUM(IFOIL-1))
        ENDIF
C
C ---------------------------------------------------------------------
C          For dispersed beam choose vertical position.  The energy
C          and vertical position are assumed to be correlated.
C ---------------------------------------------------------------------
C
        CALL RANECU(RAN_UNIF,1)
        IF(BEAMV .NE. 0.D0)THEN
          BEAMV_I = (DBLE(RAN_UNIF)-0.5D0)*BEAMV
          E0_N = (BEAMV_I/BEAMV+0.5D0)*DE0 + E0_MIN
        ELSE
          BEAMV_I = 0.D0
          E0_N = E0
        ENDIF
C
C ---------------------------------------------------------------------
C          Smear beam energy and angles by gaussian resolution
C          function and add offsets.  Fold in beam position offsets
C          later (after the raster kicks are determined).
C ---------------------------------------------------------------------
C
        E0_OFF   = GASDEVV(GSIG_E0)   + DEL_E0
        PH_B_OFF = GASDEVV(GSIG_PH_B) + DEL_PH_B
        TH_B_OFF = GASDEVV(GSIG_TH_B) + DEL_TH_B
C
        PH_B_N = PH_B   ! Save as nominal values
        TH_B_N = TH_B
C
        E0_I   = E0_N     + E0_OFF
        PH_B_I = PH_B_N   + PH_B_OFF
        TH_B_I = TH_B_N   + TH_B_OFF
C
C ---------------------------------------------------------------------
C          Get sines, cosines and tangets of beam angles.
C ---------------------------------------------------------------------
C
        CTH_B_I = COS(TH_B_I)
        CPH_B_I = COS(PH_B_I)
        STH_B_I = SIN(TH_B_I)
        SPH_B_I = SIN(PH_B_I)
        TTH_B_I = STH_B_I/CTH_B_I
        TPH_B_I = SPH_B_I/CPH_B_I
C
C ---------------------------------------------------------------------
C          Determine the DELTA_X and DELTA_Y kicks due to the raster
C          No raster implies DELTA_X=0 and DELTA_Y=0.
C          Note that the raster kicks are defined in the Z=0 LAB
C          system (i.e. at a vertical plane through the Hall center
C          normal to the nominal beam direction).
C ---------------------------------------------------------------------
C
        IF( RAST_SHAPE.EQ.'R' ) THEN
           CALL RANECU(RAN_UNIF,1)
           DELTA_X=0.5D0*RAST_AMP_X*SIN(DBLE(RAN_UNIF)*2.D0*PI)
           CALL RANECU(RAN_UNIF,1)
           DELTA_Y=0.5D0*RAST_AMP_Y*SIN(DBLE(RAN_UNIF)*2.D0*PI)
        ELSEIF( RAST_SHAPE.EQ.'E' ) THEN
           CALL RANECU(RAN_UNIF,1)
           CALL RANECU(RAN_UNIF2,1)
           DELTA_X=0.5D0*RAST_AMP_X*SQRT(DBLE(RAN_UNIF))*
     #              COS(DBLE(RAN_UNIF2)*2.D0*PI)
           DELTA_Y=0.5D0*RAST_AMP_Y*SQRT(DBLE(RAN_UNIF))*
     #              SIN(DBLE(RAN_UNIF2)*2.D0*PI)
        ELSE           
           DELTA_X=0.D0
           DELTA_Y=0.D0
        ENDIF
C
C ---------------------------------------------------------------------
C          Add beam offsets to above values.  These offsets (like
C          the raster kicks) are defined in the Z=0 LAB system.
C ---------------------------------------------------------------------
C
        DELTA_X = DELTA_X + GASDEVV(GSIG_X_B) + DEL_X_B
        DELTA_Y = DELTA_Y + GASDEVV(GSIG_Y_B) + DEL_Y_B + BEAMV_I
C
C ---------------------------------------------------------------------
C          Now get the vertex.
C
C          To do this, I take the intersection of the beam (defined
C          with the above X and Y offsets (wrt Z=0 plane) and angles
C          PH_B_I and TH_B_I) with a possibly tilted (for nonzero
C          PHI_TARG) vertical plane passing through the point
C          (X=0,Y=0,TARG_I).  This tilted plane is parallel to the
C          surface of the foil.  For cryotargets, PHI_TARG has been
C          set to zero, so this poses no problem.
C ---------------------------------------------------------------------
C
        VERTEX_X = (DELTA_X+TARG_I*TPH_B_I)/(1.D0+TPHI_TARG*TPH_B_I)
        VERTEX_Z = -VERTEX_X*TPHI_TARG + TARG_I
        VERTEX_Y = DELTA_Y + VERTEX_Z*TTH_B_I/CPH_B_I
C
C ---------------------------------------------------------------------
C       VERSION 3.2 and later:
C
C          For spectrometer mispointing effects:
C
C          Determine position of interaction point in modified LAB
C          systems.  A modified LAB system is defined for each
C          spectrometer by shifting the origin by the vector
C          OBJECT_PT.  This vector gives the spectrometer object point
C          relative to the nominal (i.e. unmodified) LAB system
C          and allows for possible spectrometer mispointing.  The
C          sampling of the angles is done relative to these interaction
C          points.  Thus, the spectrometer mispointing is affected
C          by moving the interaction point in the opposite direction.
C ---------------------------------------------------------------------
C
        BEAM_E1 = VERTEX_X - OBJECT_PT(1,1)   ! electron arm
        BEAM_E2 = VERTEX_Y - OBJECT_PT(1,2)
        BEAM_E3 = VERTEX_Z - OBJECT_PT(1,3)
C
        BEAM_P1 = VERTEX_X - OBJECT_PT(2,1)   ! hadron arm
        BEAM_P2 = VERTEX_Y - OBJECT_PT(2,2)
        BEAM_P3 = VERTEX_Z - OBJECT_PT(2,3)
C
        X_BEAM_E = -BEAM_E2   ! These are needed to make X0 correction
        X_BEAM_P = -BEAM_P2   ! in HRS_INV (Note: it's assumed that
                              ! the spectrometer X is vertically down).
        Y_BEAM_E =  BEAM_E1
        Y_BEAM_P =  BEAM_P1
C
C ---------------------------------------------------------------------
C          Perform energy loss calculation for beam.  This requires
C          determining the path length through the target for the beam.
C          Also determine the path lengths for the outgoing particles
C          here.  This necessarily involves an approximation if the
C          path lengths and energy are to be calculated in the same
C          routine:
C             Ideally, we should calculate the path length from the
C             actual angles for the outgoing particles.  However,
C             for elastic scattering the adjusted beam energy should
C             be used to get the hadron angles (i.e. the hadron path
C             length isn't known until the beam energy is determined).
C             Therefore, use the central angles for the spectrometers
C             to calculate the path lengths.
C ---------------------------------------------------------------------
C
        IF(ELOSS_CALC) THEN
C
           IF(ITARG_MOD .EQ. 1) THEN          ! Target foil(s)
C
              I_WALL_B = 1    ! Change this for version 3.9.
              I_WALL_E = 2     ! We still want the walls which 
              I_WALL_H = 4      ! are not target cell related.
                                 ! Here, 2 and 3 are equivalent
                                  ! as are 4 and 5.
C
              THIN = TARG0
C
              IF(ABS(PH_E-PHI_TARG) .LE. PI/2.D0) THEN
                 PATH_TMP = TARG_CUM(IFOIL) - TARG0
                 ETHO = PATH_TMP*CPHI_TARG/CPHIE_DIFF
              ELSE
                 IF(IFOIL .EQ. 1) THEN
                    PATH_TMP = TARG0
                 ELSE
                    PATH_TMP = TARG0 - TARG_CUM(IFOIL-1)
                 ENDIF
                 ETHO = -PATH_TMP*CPHI_TARG/CPHIE_DIFF
              ENDIF
C
              IF(EJECT_CHG .NE. 0) THEN     ! Don't calculate for neutrals
                IF(ABS(PH_P-PHI_TARG) .LE. PI/2.D0) THEN
                   PATH_TMP = TARG_CUM(IFOIL) - TARG0
                   HTHO = PATH_TMP*CPHI_TARG/CPHIP_DIFF
                ELSE
                   IF(IFOIL .EQ. 1) THEN
                      PATH_TMP = TARG0
                   ELSE
                      PATH_TMP = TARG0 - TARG_CUM(IFOIL-1)
                   ENDIF
                   HTHO = -PATH_TMP*CPHI_TARG/CPHIP_DIFF
                ENDIF
              ENDIF
C
           ELSE                               ! LH2/LD2 or He cryotarget
              I_WALL_B = 1   ! Entrance wall selected for beam energy loss
              THIN = TARG0
              CALL TG_PATH(VERTEX_X,VERTEX_Y,VERTEX_Z,
     #                     PH_E,ETHO,EOPT,ITARG_MOD)
              IF(EJECT_CHG .NE. 0) CALL TG_PATH(VERTEX_X,VERTEX_Y,
     #                           VERTEX_Z,PH_P,HTHO,HOPT,ITARG_MOD)
C
              IF(EOPT .EQ. 0) THEN
                 I_WALL_E = 2       ! sidewall for scattered electron
              ELSE
                 I_WALL_E = 3       ! endcap for scattered electron
              ENDIF
C
              IF(EJECT_CHG .NE. 0) THEN     ! Don't calculate for neutrals
                IF(HOPT .EQ. 0) THEN
                   I_WALL_H = 4       ! sidewall for hadron
                ELSE
                   I_WALL_H = 5       ! endcap for hadron
                ENDIF
              ENDIF
C
           ENDIF
C
           CALL EXT_BREMSS(THIN,E0_I,ERAD,I_WALL_B)
           CALL ELOSS_E(THIN,E0_I,1,ELOSS,1,I_WALL_B)
           CALL MULTI_SCANG(THIN,E0_I,BVARPHI,BVARTH,M_ELEC,1)
C
Cpeu           ELOSS_B_I = ELOSS + ERAD + ELWALK(1)  ! Beam energy loss
Cpeu           Walls are already included in eloss_e!
C
           ELOSS_B_I = ELOSS + ERAD  ! Beam energy loss
C
           E0_I   = E0_I   - ELOSS_B_I
C
           IF(E0_I .LE. 0.511D0) THEN   ! Protect against rare events
              FAIL_KIN = .TRUE.         ! with huge energy loss (ERAD)
              GOTO 999
           ENDIF
C
           PH_B_I = PH_B_I + BVARPHI
           TH_B_I = TH_B_I + BVARTH
C
        ENDIF
C
C ---------------------------------------------------------------------
C          Choose electron angles first.  For elastic scattering,
C          the electron angles (coupled with the beam) define
C          everything.  For knockout to bound states, the hadron
C          momentum is determined from the electron kinematics,
C          the hadron emission angles and the known missing mass.
C          Otherwise (for continuum case) all six variables (e- and
C          hadron momentum, in-plane angles and out-of-plane angles)
C          are randomly selected.
C ---------------------------------------------------------------------
C
        CALL RANECU(RAN_UNIF,1)
        TH_E_N =  TH_E_MIN + DBLE(RAN_UNIF)*DTH_E
C
        CALL RANECU(RAN_UNIF,1)
        PH_E_N =  PH_E_MIN + DBLE(RAN_UNIF)*DPH_E
C
C ---------------------------------------------------------------------
C          Check to see whether the chosen angles lie within the
C          solid-angle defining aperture.  For rectangular apertures
C          and randomly selected "nominal" angles this test is
C          unnecessary since every ray will fall within the aperture
C          by definition.  However, for elastic scattering (where the
C          nominal hadron angles are COMPUTED from actual electron
C          angles) this need not be the case.  For elliptical
C          apertures this method, although slightly inefficient,
C          guarantees a uniform sampling over the aperture.
C ---------------------------------------------------------------------
C
        IF(.NOT. IN_APERTURE(SA_SHAPE_E,PH_E_MIN,TH_E_MIN,
     #                       PH_E_N,TH_E_N) ) GOTO 999
C
C ---------------------------------------------------------------------
C          Add central angles to get absolute nominal angles in the LAB 
C          (as opposed to angles relative to the spectrometer central
C          ray (i.e. "collimator" angles)).  The division by
C          COS(theta) accounts for the difference between the
C          collimator angular range and the initial beam "transport"
C          system PHI range.
C ---------------------------------------------------------------------
C
        TH_E_N =  TH_E_N + TH_E
        PH_E_N =  PH_E_N/COS(TH_E_N) + PH_E
C
C ---------------------------------------------------------------------
C          Get actual angles in the LAB system from nominal angles,
C          interaction point and nominal drift distance to aperture.
C          The "actual" and "nominal" angles are the same for a
C          point target.
C ---------------------------------------------------------------------
C
        CALL NOM_TO_ACT(PH_E_N,TH_E_N,DRIFTE_N,BEAM_E1,BEAM_E2,BEAM_E3,
     #                  PH_E_I,TH_E_I,JACOB_SA_E)       !electron
C
C ---------------------------------------------------------------------
C          A NOTE about the counter SUM_ACC (used below) is necessary.
C
C          This counts the number of events which fall within
C          the acceptance.  The total acceptance volume divided
C          by this number then gives the effective acceptance per
C          event.  Here, "acceptance" refers to the set of variables
C          in which the cross section is differential.  Thus, for
C          elastic scattering the acceptance refers only to the
C          electron solid angle.  For a bound state
C          the acceptance refers to everything except the hadron
C          momentum.  With this prescription, the yields will be
C          properly normalized, since the cross section times the
C          effective acceptance per event times the luminosity times the
C          counting time gives the yield for an "effective" bin.
C          The sum over events (or "effective" bins) then gives the
C          total yield.
C ---------------------------------------------------------------------
C
C      Now we do the SAMPLING which is specific to the various
C      kinematical cases:  ELASTIC, BOUND, CONTINUUM.
C
C ---------------------------------------------------------------------
C     ELASTIC SCATTERING:
C
C          For elastic scattering, calculate the electron final
C          energy and the hadron momentum vector (=Q).  Then test
C          to see if momenta fall within acceptances (for elastic
C          scattering the hadron-arm acceptance tests are
C          overriden if the logical flag ACCEPT_CHECK is false;
C          this allows simulation of single-arm elastic scattering).
C
C          Now, also allows for radiation of a photon.  Since the
C          photon energy is directly sampled, we must still determine
C          the asymptotic scattered electron energy from the
C          elastic kinematics with radiation of a "known" photon.
C ---------------------------------------------------------------------
C
        IF(ELASTIC) THEN
           SUM_ACC = SUM_ACC + 1.D0    !number of tries within accept.
           CALL V3MAKE(E0_I,PH_B_I,TH_B_I,KI_VEC)  !make beam 3-vector
C
           IF(RADPEAKING) THEN
              CALL PEAKING(PH_B_I,TH_B_I,PH_E_I,TH_E_I,
     #                      EG_PK,KG_PK,KG_PK_UNIT)
              NCASES = 2
              DO I=1,2     ! Loop over pre- and post-radiation
                EG = EG_PK(I)
                DO J=1,3
                   KG(J)      = KG_PK(I,J)
                   KG_UNIT(J) = KG_PK_UNIT(I,J)
                ENDDO
                CALL PSEUDO_BEAM(E0_I,KI_VEC,EG,KG,EB_PS,KI_PS)
                CALL KINEM_ELAS_EF(RADPEAKING,EB_PS,KI_PS,E0_I,KI_VEC,
     #                         KG_UNIT,PH_E_I,TH_E_I,PF_E_I,DEF_DK)
                CALL V3MAKE(PF_E_I,PH_E_I,TH_E_I,KF_VEC)
                CALL KINEM_ELAS_PF(KI_PS,KF_VEC,PF_P_I,PH_P_I,TH_P_I)
                DEF_DK_TMP(I) = DEF_DK
                PF_E_I_TMP(I) = PF_E_I
                PF_P_I_TMP(I) = PF_P_I
                PH_P_I_TMP(I) = PH_P_I
                TH_P_I_TMP(I) = TH_P_I
              ENDDO
C
           ELSEIF(RADFULL) THEN
C
C             Here, I choose the vertex kinematics, constrained to
C             elastic scattering, and then let bremsgen modify the
C             particle momenta according to the radiated photon.
C             As the sampling occurs in the vertex variables, no
C             extra Jacobian is needed for the elastic cross section.
C
              NCASES = 1
              E0_VTX   = E0_I
              PH_B_VTX = PH_B_I
              TH_B_VTX = TH_B_I
              PH_E_VTX = PH_E_I
              TH_E_VTX = TH_E_I
              DO J=1,3
                 KI_VTX(J) = KI_VEC(J)
              ENDDO
              CALL KINEM_ELAS_EF(.FALSE.,E0_VTX,KI_VTX,E0_VTX,KI_VTX,
     #                  KG_UNIT,PH_E_VTX,TH_E_VTX,EF_VTX,DEF_DK)
              CALL V3MAKE(EF_VTX,PH_E_VTX,TH_E_VTX,KF_VTX)
              CALL V3DIF(KI_VTX,KF_VTX,QQ_VEC_VTX)
              QQ_VTX = DOT(QQ_VEC_VTX,QQ_VEC_VTX) - (E0_VTX-EF_VTX)**2
              CALL KINEM_ELAS_PF(KI_VTX,KF_VTX,
     #                           PF_VTX,PH_P_VTX,TH_P_VTX)
              CALL V3MAKE(PF_VTX,PH_P_VTX,TH_P_VTX,KP_VTX)
              EP_VTX     = SQRT(PF_VTX**2 + EJECT_MASS**2)
              THETAE_VTX = ACOS(DOT(KI_VTX,KF_VTX)/(E0_VTX*EF_VTX))
              THETAP_VTX = ACOS(DOT(KI_VTX,KP_VTX)/(E0_VTX*PF_VTX))
              IF(PH_E_VTX .LT. 0.) THETAE_VTX = -THETAE_VTX
              IF(PH_P_VTX .LT. 0.) THETAP_VTX = -THETAP_VTX
C
C             I pass the angles to BREMSGEN with their actual signs.
C             Inside BREMSGEN, I form new variables which are the
C             absolute values of these angles, since that's how the
C             BREMSGEN "angles" routine likes them (i.e. "angles" 
C             explicitly makes the electron "y" component of momentum 
C             ("y" defined in the bremsgen system) negative).
C
              CALL BREMSGEN(QQ_VTX,E0_VTX,THETAE_VTX,THETAP_VTX,CUTOFF,
     #              ENEW,EPRIMENEW,POPNEW)
              E0_I   = E0_VTX   - ENEW
              PF_E_I = EF_VTX   - EPRIMENEW
              EP_ASM = EP_VTX   - POPNEW
              PF_P_I = SQRT(EP_ASM**2 - EJECT_MASS**2)
C
              PH_B_I = PH_B_VTX + WINKEL6
              TH_B_I = TH_B_VTX + WINKEL5
              PH_E_I = PH_E_VTX + WINKEL2
              TH_E_I = TH_E_VTX + WINKEL1
              PH_P_I = PH_P_VTX + WINKEL4
              TH_P_I = TH_P_VTX + WINKEL3
C
           ELSE
              NCASES = 1
              CALL KINEM_ELAS_EF(.FALSE.,E0_I,KI_VEC,E0_I,KI_VEC,
     #                  KG_UNIT,PH_E_I,TH_E_I,PF_E_I,DEF_DK)
              CALL V3MAKE(PF_E_I,PH_E_I,TH_E_I,KF_VEC)
              CALL KINEM_ELAS_PF(KI_VEC,KF_VEC,PF_P_I,PH_P_I,TH_P_I)
           ENDIF
C
           DO I=1,NCASES      !Test kinematics (for rad. test pre and post)
             FAIL_RAD(I) = .FALSE.    !Initialize
C
             IF(RADPEAKING) THEN
                PF_E_I = PF_E_I_TMP(I)
                PF_P_I = PF_P_I_TMP(I)
                PH_P_I = PH_P_I_TMP(I)
                TH_P_I = TH_P_I_TMP(I)
             ENDIF
C
             IF((PF_E_I .LT. EF_MIN) .OR. (PF_E_I .GT. EF_MAX)) GOTO 777
             IF(ACCEPT_CHECK .AND. ((PF_P_I .LT. PF_MIN)
     #                       .OR. (PF_P_I .GT. PF_MAX))) GOTO 777
             CALL ACT_TO_NOM(PH_P_I,TH_P_I,DRIFTP_N,BEAM_P1,
     #                      BEAM_P2,BEAM_P3,PH_P_N,TH_P_N)
             PH_P_N =  (PH_P_N-PH_P)*COS(TH_P_N) !rel. to aperture center
             TH_P_N =  TH_P_N - TH_P             !rel. to aperture center
             IF(ACCEPT_CHECK .AND. .NOT. IN_APERTURE(SA_SHAPE_P,
     #            PH_P_MIN,TH_P_MIN,PH_P_N,TH_P_N) ) GOTO 777
C
             IF(RADFULL) THEN    ! Must also test electron aperture
                CALL ACT_TO_NOM(PH_E_I,TH_E_I,DRIFTE_N,BEAM_E1,
     #                          BEAM_E2,BEAM_E3,PH_E_N,TH_E_N)
                PH_E_N =  (PH_E_N-PH_E)*COS(TH_E_N) !rel. to aperture center
                TH_E_N =  TH_E_N - TH_E             !rel. to aperture center
                IF(.NOT. IN_APERTURE(SA_SHAPE_E,
     #              PH_E_MIN,TH_E_MIN,PH_E_N,TH_E_N) ) GOTO 777
             ENDIF
C
             GOTO 888
 777         FAIL_RAD(I) = .TRUE.
 888      ENDDO
C
          IF(RADPEAKING) THEN    ! Select either pre- or post-rad kinem.
             IF(FAIL_RAD(1) .AND. FAIL_RAD(2)) GOTO 999  !Get next event
             IF(FAIL_RAD(1)) THEN          !Select post-radiation
                I_RAD = 2
                RAD_WT = 1.D0
             ELSEIF(FAIL_RAD(2)) THEN      !Select pre-radiation
                I_RAD = 1
                RAD_WT = 1.D0
             ELSE                          !Both passed:  toss coin
                CALL RANECU(RAN_UNIF,1)
                IF(RAN_UNIF .LE. 0.5) THEN
                   I_RAD = 1
                ELSE
                   I_RAD = 2
                ENDIF
                RAD_WT = 2.D0     !Account for photon not selected
             ENDIF
             DEF_DK = DEF_DK_TMP(I_RAD)   !Jacobian
             PF_E_I = PF_E_I_TMP(I_RAD)
             PF_P_I = PF_P_I_TMP(I_RAD)
             PH_P_I = PH_P_I_TMP(I_RAD)
             TH_P_I = TH_P_I_TMP(I_RAD)
             EG = EG_PK(I_RAD)
             RAD_WT = RAD_WT/(1.D0/EG)   !Account for 1/k sampling
             DO J=1,3
                KG(J) = KG_PK(I_RAD,J)
             ENDDO
          ELSE
             IF(FAIL_RAD(1)) GOTO 999       !Get next event
             RAD_WT = 1.D0
          ENDIF
C
C ---------------------------------------------------------------------
C     BOUND STATE:
C          For bound state sample electron momentum and hadron angles.
C          If within aperture, calculate hadron momentum from
C          the known missing mass.
C
C          If we include radiation, then we calculate the final ejectile
C          momentum, consistent with the given radiated photon and
C          the assumed missing mass.  To do this, the photon quantities
C          (energy and momentum) are subtracted from the asymptotic
C          beam quantities, resulting in a "pseudo-beam".  For pre-
C          radiation this corresponds to the (possibly offshell)
C          incident electron.  For post-radiation, this is an artificial
C          quantity.  In both cases, the ejectile momentum can be
C          calculated from the difference between the pseudo-beam and
C          the asymptotic momentum of the scattered electron (i.e. this
C          necessarily corresponds to the virtual photon for pre- or
C          post-radiation).
C ---------------------------------------------------------------------
C
        ELSEIF(BOUND) THEN
           CALL RANECU(RAN_UNIF,1)
           TH_P_N =  DBLE(RAN_UNIF)*DTH_P + TH_P_MIN
           CALL RANECU(RAN_UNIF,1)
           PH_P_N =  DBLE(RAN_UNIF)*DPH_P + PH_P_MIN
           IF(.NOT. IN_APERTURE(SA_SHAPE_P,PH_P_MIN,TH_P_MIN,
     #                      PH_P_N,TH_P_N) ) GOTO 999
           CALL RANECU(RAN_UNIF,1)
           PF_E_I =  DBLE(RAN_UNIF)*DPF_E + PF_E_MIN
           SUM_ACC = SUM_ACC + 1.D0    !number of tries within accept.
           TH_P_N =  TH_P_N + TH_P
           PH_P_N =  PH_P_N/COS(TH_P_N) + PH_P
           CALL NOM_TO_ACT(PH_P_N,TH_P_N,DRIFTP_N,BEAM_P1,
     #                     BEAM_P2,BEAM_P3,PH_P_I,TH_P_I,JACOB_SA_P)
           MISS_M = EM_BOUND  !Reset missing mass; get hadron momentum
           CALL V3MAKE(E0_I,PH_B_I,TH_B_I,KI_VEC)    !make beam 3-vector
           CALL V3MAKE(PF_E_I,PH_E_I,TH_E_I,KF_VEC)  !make scatt e- 3-vector
C
           IF(RADPEAKING) THEN
              CALL PEAKING(PH_B_I,TH_B_I,PH_E_I,TH_E_I,
     #                      EG_PK,KG_PK,KG_PK_UNIT)
              DO I=1,2     ! Loop over pre- and post-radiation
                EG = EG_PK(I)
                DO J=1,3
                   KG(J) = KG_PK(I,J)
                ENDDO
                CALL PSEUDO_BEAM(E0_I,KI_VEC,EG,KG,EB_PS,KI_PS)
                CALL KINEM(EB_PS,KI_PS,PF_E_I,KF_VEC,PF_P_I,PH_P_I,
     #                 TH_P_I,ELASTIC,BOUND,.TRUE.,FAIL_KIN)
                IF(I.EQ.1) FIRST_TIME = .TRUE.  !Force new Pf calc.
                I_PF_ROOT_TMP(I) = I_PF_ROOT
                PF_ROOT_WT_TMP(I) = PF_ROOT_WT
                PF_P_I_TMP(I) = PF_P_I
                FAIL_RAD(I) = FAIL_KIN
              ENDDO
           ELSE
              CALL KINEM(E0_I,KI_VEC,PF_E_I,KF_VEC,PF_P_I,PH_P_I,TH_P_I,
     #                 ELASTIC,BOUND,.TRUE.,FAIL_KIN)
              FAIL_RAD(1) = FAIL_KIN
           ENDIF
C
           IF(RADPEAKING) THEN    ! Select either pre- or post-rad kinem.
             IF(FAIL_RAD(1) .AND. FAIL_RAD(2)) GOTO 999  !Get next event
             IF(FAIL_RAD(1)) THEN          !Select post-radiation
                I_RAD = 2
                RAD_WT = 1.D0
                PF_ROOT_WT = PF_ROOT_WT_TMP(2)
             ELSEIF(FAIL_RAD(2)) THEN      !Select pre-radiation
                I_RAD = 1
                RAD_WT = 1.D0
                PF_ROOT_WT = PF_ROOT_WT_TMP(1)
             ELSE                          !Both passed:  toss coin
                CALL RANECU(RAN_UNIF,1)
                IF(RAN_UNIF .LE. 0.5) THEN
                   I_RAD = 1
                ELSE
                   I_RAD = 2
                ENDIF
                RAD_WT = 2.D0   ! Account for photon not selected
                PF_ROOT_WT = PF_ROOT_WT_TMP(I_RAD)
             ENDIF
             PF_P_I = PF_P_I_TMP(I_RAD)
             I_PF_ROOT = I_PF_ROOT_TMP(I_RAD)
             EG = EG_PK(I_RAD)
             RAD_WT = RAD_WT/(1.D0/EG)  ! Account for 1/k sampling
             DO J=1,3
                KG(J) = KG_PK(I_RAD,J)
             ENDDO
           ELSE
             IF(FAIL_RAD(1)) GOTO 999       !Get next event
             RAD_WT = 1.D0
           ENDIF
C
C ---------------------------------------------------------------------
C     CONTINUUM CASE:
C          For continuum case sample electron momentum, hadron momentum
C          and hadron angles.  Perform acceptance test.
C
C          For radiation, the photon sampling is done here.
C ---------------------------------------------------------------------
C
        ELSE
           CALL RANECU(RAN_UNIF,1)
           TH_P_N =  DBLE(RAN_UNIF)*DTH_P + TH_P_MIN
           CALL RANECU(RAN_UNIF,1)
           PH_P_N =  DBLE(RAN_UNIF)*DPH_P + PH_P_MIN
           IF(.NOT. IN_APERTURE(SA_SHAPE_P,PH_P_MIN,TH_P_MIN,
     #                      PH_P_N,TH_P_N) ) GOTO 999
           CALL RANECU(RAN_UNIF,1)
           PF_E_I =  DBLE(RAN_UNIF)*DPF_E + PF_E_MIN
           CALL RANECU(RAN_UNIF,1)
           PF_P_I =  DBLE(RAN_UNIF)*DPF_P + PF_P_MIN
           SUM_ACC = SUM_ACC + 1.D0    !number of tries within accept.
           TH_P_N =  TH_P_N + TH_P
           PH_P_N =  PH_P_N/COS(TH_P_N) + PH_P
           CALL NOM_TO_ACT(PH_P_N,TH_P_N,DRIFTP_N,BEAM_P1,
     #                     BEAM_P2,BEAM_P3,PH_P_I,TH_P_I,JACOB_SA_P)
C
           IF(RADPEAKING) THEN
              CALL PEAKING(PH_B_I,TH_B_I,PH_E_I,TH_E_I,
     #                      EG_PK,KG_PK,KG_PK_UNIT)
              CALL RANECU(RAN_UNIF,1)       !Select pre- or post-rad
              IF(RAN_UNIF .LE. 0.5) THEN
                 I_RAD = 1
              ELSE
                 I_RAD = 2
              ENDIF
              RAD_WT = 2.D0     !Account for photon not selected
              EG = EG_PK(I_RAD)
              RAD_WT = RAD_WT/(1.D0/EG) ! Account for 1/k sampling
              DO J=1,3
                 KG(J) = KG_PK(I_RAD,J)
              ENDDO
           ELSE
              RAD_WT = 1.D0
           ENDIF
        ENDIF
C
C ---------------------------------------------------------------------
C          Save Laboratory coordinates for later.
C ---------------------------------------------------------------------
C
        DO I=1,6
          QI_E_SAV(I) = QI_E(I)
          QI_P_SAV(I) = QI_P(I)
        ENDDO
C
C ---------------------------------------------------------------------
C          Calculate Transport coordinates for each arm (at target).
C ---------------------------------------------------------------------
C
        IF(CALLSPEC) THEN
          CALL LAB_TO_TRPT(QI_E,ROT_E,PF_E,TRPT_T_E)
          CALL LAB_TO_TRPT(QI_P,ROT_P,PF_P,TRPT_T_P)
          IF(ACCEPT_FCN) THEN         !Save target coordinates
             DO I=1,6
                TRPT_T_E_SAV(I) = TRPT_T_E(I)
                TRPT_T_P_SAV(I) = TRPT_T_P(I)
             ENDDO
          ENDIF
C
C ---------------------------------------------------------------------
C          Do spectrometer analysis.
C ---------------------------------------------------------------------
C
          CALL SPECTROMETER(1)
C
C ---------------------------------------------------------------------
C          Calculate Laboratory coordinates from (modified) target
C          coordinates for each arm.
C ---------------------------------------------------------------------
C
          IF(APPLY_TO_LAB(1) .AND. .NOT. ACCEPT_FCN)
     #          CALL TRPT_TO_LAB(TRPT_T_E,ROTI_E,PF_E,QI_E)
          IF(APPLY_TO_LAB(2) .AND. .NOT. ACCEPT_FCN)
     #          CALL TRPT_TO_LAB(TRPT_T_P,ROTI_P,PF_P,QI_P)
C
        ENDIF
C
C ---------------------------------------------------------------------
C          If we are calculating the acceptance function, then
C          write the Ntuple here.  This will be the only output
C          of MCEEP in this case.  The Ntuple contains the solid
C          angle Jacobians and the initial and final Transport vectors.
C ---------------------------------------------------------------------
C
        IF(ACCEPT_FCN) THEN
           IF(.NOT. FAIL_GLOBAL_CUT) THEN
              WRITE(41,1112) JACOB_SA_E,JACOB_SA_P,
     #                       (TRPT_T_E_SAV(I),I=1,6),
     #                       (TRPT_T_P_SAV(I),I=1,6),
     #                       (TRPT_T_E(I),I=1,6),
     #                       (TRPT_T_P(I),I=1,6)
           ENDIF
           GOTO 999               ! Done; get next event
        ENDIF
C
C ---------------------------------------------------------------------
C          Calculate derived kinematical quantities for calculation
C          of the cross sections.
C
C          Note:  "ASM" suffix stands for asymptotic kinematics
C                 "VTX" suffix stands for vertex     kinematics
C
C          For the cross sections (note that these are the unradiated
C          cross sections even if the radiative tail is being
C          calculated, since the current radiative tail routines
C          involve factors times the unradiated cross sections), we
C          must use the vertex kinematics.  The multiplicative factors
C          in the Mo-Tsai and Borie-Drechsel formulas involve the
C          asymptotic variables however.
C ---------------------------------------------------------------------
C
        CALL V3MAKE(E0_I,PH_B_I,TH_B_I,KI_VEC)    !make beam 3-vector
        CALL V3MAKE(PF_E_I,PH_E_I,TH_E_I,KF_VEC)  !make scatt e- 3-vector
C
        IF(RADPEAKING) THEN
           IF(I_RAD .EQ. 1) THEN     !Pre-radiation
              E0_ASM = E0_I
              E0_VTX = E0_I - EG
              CALL V3DIF(KI_VEC,KG,KI_VTX)
              EF_ASM = PF_E_I
              EF_VTX = PF_E_I
              DO I=1,3
                 KF_VTX(I) = KF_VEC(I)
              ENDDO
           ELSE                      !Post-radiation
              E0_ASM = E0_I
              E0_VTX = E0_I
              DO I=1,3
                 KI_VTX(I) = KI_VEC(I)
              ENDDO
              EF_ASM = PF_E_I
              EF_VTX = PF_E_I + EG
              CALL V3SUM(KF_VEC,KG,KF_VTX)
           ENDIF
           PF_ASM   = PF_P_I
           PF_VTX   = PF_P_I
           PH_P_VTX = PH_P_I
           TH_P_VTX = TH_P_I
C
        ELSEIF(RADFULL) THEN         !Vertex values were already set
           E0_ASM = E0_I
           EF_ASM = PF_E_I
           PF_ASM = PF_P_I
C
        ELSE
           E0_ASM = E0_I
           E0_VTX = E0_I
           EF_ASM = PF_E_I
           EF_VTX = PF_E_I
           DO I=1,3
              KI_VTX(I) = KI_VEC(I)
              KF_VTX(I) = KF_VEC(I)
           ENDDO
           PF_ASM   = PF_P_I
           PF_VTX   = PF_P_I
           PH_P_VTX = PH_P_I
           TH_P_VTX = TH_P_I
        ENDIF
C
C ---------------------------------------------------------------------
C          Set the particle momentum magnitudes to vertex values!
C
C          This is necessary since many physics routines involve these
C          variables explicitly and since KINEM does not recalculate
C          them (i.e. E0_I, PF_E_I and PF_P_I are part of var.cmn and
C          many physics routines use these versions of the variables
C          to calculate cross sections).
C
C          Care must be taken in what follows to make sure that
C          the correct variables (VTX or ASM) are used.  For ICALL=2
C          E0_I, PF_E_I and PF_P_I are reset to the asymptotic values,
C          so the histogram binning is not affected by setting them
C          to the vertex variables here.
C ---------------------------------------------------------------------
C
        E0_I   = E0_VTX
        PF_E_I = EF_VTX
        PF_P_I = PF_VTX
C
        IF( E0_I   .LE. 0.D0  .OR. 
     #      PF_E_I .LE. 0.D0  .OR. 
     #      PF_P_I .LE. 0.D0       ) THEN  ! Protect against rare
           FAIL_KIN = .TRUE.                ! events with huge eloss
           GOTO 999
        ENDIF
C
        CALL KINEM(E0_VTX,KI_VTX,EF_VTX,KF_VTX,PF_VTX,
     #            PH_P_VTX,TH_P_VTX,ELASTIC,BOUND,.FALSE.,FAIL_KIN)
        IF(FAIL_KIN) GOTO 999
C
C ---------------------------------------------------------------------
C          Histogram various vertex kinematical quantities.
C          In the absense of radiation, these are the same
C          as the "measured" (i.e. asymptotic) kinematics.
C ---------------------------------------------------------------------  
C
        E0_V    = E0_VTX
        EF_V    = EF_VTX
        OMEGA_V = OMEGA
        QMAG_V  = QMAG
        QMU2_V  = QMU2_G
        PHIX_V  = PHI_X
        EMISS_V = EMISS
        MISSM_V = MISS_M
        PRMAG_V = PRMAG
        W_V     = W
        THCM_V  = THETA_CM
        THPQ_V  = THETA_PQ
        TSCAT_V = TSCAT
C
C ---------------------------------------------------------------------
C          Determine the phase-space weighting factor.
C
C          In the case of a bound state the weight, PF_ROOT_WT,
C          is either 1 or 2 depending upon whether 1 or 2 physical
C          solutions are found for the ejectile momentum within
C          the hadron arm momentum acceptance.  If 2 solutions are
C          found then the the phase space is effectively double for
C          the event.  See next paragraph for more details.
C          (PF_ROOT_WT = 1 if not BOUND).
C
C          For bound state calculations the ejectile momentum is
C          determined from the electron kinematics, the ejectile
C          emission angles and the missing mass.  In some cases
C          (especially for light recoils) there are two physical
C          solutions for the hadron momentum corresponding to ejectile
C          motion along and against q in the center of mass of the
C          (final) ejectile-recoil system.  When both solutions
C          are within the hadron arm momentum acceptance, one of them
C          is randomly selected by a coin toss.  Since each solution
C          is chosen only half the time (on average) a weighting factor
C          of 2 is applied.  If only one solution is within
C          the acceptance, the weighting factor (PF_ROOT_WT)
C          is equal to 1.  Cross section (polarization) histograms,
C          which are obtained by dividing the yield (partial yield)
C          histograms by the phase space (unpolarized yield)
C          histograms channel-by-channel, are not affected by this
C          weighting factor since the phase space histograms AND
C          the yields both include the PF_ROOT_WT factor.
C ---------------------------------------------------------------------
C
        IF(.NOT. BOUND) PF_ROOT_WT = 1.D0
        IF(ELASTIC) THEN           ! Solid angle Jacobians (ACT-->NOM)
           JACOB_SA = JACOB_SA_E
        ELSE
           JACOB_SA = JACOB_SA_E * JACOB_SA_P
        ENDIF
        PHASE_SP = PF_ROOT_WT*JACOB_SA
     #               *RAD_WT*RAD_NORM(ICALC)*DBLE(NCALC)
                                     !Phase space weighting
C
C ---------------------------------------------------------------------
C          Determine radiative weighting factors for cross section.
C
C          For reasons similar to those for PF_ROOT_WT (see comments
C          above), multiply by RAD_WT which accounts for
C          possible multiple solutions for radiated photon, only one
C          of which is actually selected.  For RADPEAKING and BOUND,
C          there are two possible Pf roots for pre- and for post-
C          radiation.
C
C          Note that NCALC is also a weighting factor, since for
C          RADPEAKING (NCALC=2), we are effectively only sampling
C          the peak (and the tail) half the time.
C
C          The multiplicative factors calculated by MO_TSAI and
C          by BORIE_DRECHSEL involve the asymptotic kinematics.
C          Since these are for peaking approx. only, TSCAT is
C          the same for asymptotic or vertex kinematics.
C
C          First, determine the weighting factor.
C
C ---------------------------------------------------------------------
C
        IF(RADPEAKING) THEN
           IF(ELASTIC) THEN
              CALL MO_TSAI(I_RAD,EG,E0_ASM,EF_ASM,TSCAT,RAD_FAC)
              RAD_FAC = RAD_FAC * ABS(DEF_DK)
           ELSE
              CALL BORIE_DRECHSEL(I_RAD,EG,E0_ASM,EF_ASM,RAD_FAC)
           ENDIF
        ELSE
           RAD_FAC = 1.0D0
        ENDIF
C
        TOTAL_WT =  RAD_FAC * PHASE_SP
C
        IF(RADFULL) THEN
           CALL SCHWINGER(-1,ATARG,CUTOFF,SCHWING,COR_MULTIPHOTON)
           TOTAL_WT = TOTAL_WT * SCHWING
        ENDIF
C
        IF( (ICALC.EQ.2) .OR. ((ICALC.EQ.1).AND.MULTIPHOTON) ) THEN
           IF(ELASTIC) THEN
              CALL SCHWINGER(-1,ATARG,CUTOFF,SCHWING,COR_MULTIPHOTON)
           ELSE
              CALL SCHWINGER( 1,ATARG,CUTOFF,SCHWING,COR_MULTIPHOTON)
           ENDIF
        ENDIF
        IF(ICALC.EQ.2) THEN
           TOTAL_WT = TOTAL_WT * SCHWING    ! Schwinger correction
        ENDIF
        IF(ICALC.EQ.1 .AND. MULTIPHOTON) THEN
           TOTAL_WT = TOTAL_WT * COR_MULTIPHOTON
        ENDIF
C
C ---------------------------------------------------------------------
C          Get singles cross sections, if requested.
C
C             For (e,e'p):  get (e,e'), (e,pi-), (e,p) and (e,pi+)
C             For (e,e'n):  get (e,e'), (e,pi-), (e,n) and (e,pi0)
C
C          The singles cross sections are calculated using
C          independent models and should NOT involve the
C          vertex kinematics.  However, as singles only make
C          sense for non-radiative problem, I don't bother
C          to explicitly use the ASM variables here.
C
C          NOTE:
C              Singles calcs. only make sense for CONTINUUM sampling
C                and for RADPEAKING = .FALSE.
C              Singles cross section histograms will be slightly
C                wrong since the phase space is for double-arm
C                and thus involves both solid angle jacobians,
C                whereas the singles yields only involve a single
C                jacobian.  The yields should be correct though.
C                (PHASE_SP=JACOB_SA_E*JACOB_SA_P for continuum
C                without radiation.)  The accidentals should also
C                be correct since this involves products of
C                cross sections and the double-arm jacobian naturally
C                arises.
C ---------------------------------------------------------------------
C
        IF(SINGLES)THEN              !singles cross sections
           CALL QFSV(E0_I,PF_E_I,TSCAT,PFERMI,EPS,EPSD,ZTARG,ATARG,
     #                  SIGMA_EE)
           CALL EPC(E0_I,PF_E_I,TSCAT,ZTARG,NTARG,3,SIGMA_EPIM)
           IF(PROTON) THEN
             CALL EPC(E0_I,PF_P_I,THETA_EP,ZTARG,NTARG,4,SIGMA_EP)
             CALL EPC(E0_I,PF_P_I,THETA_EP,ZTARG,NTARG,5,SIGMA_EPIP)
           ELSE
             CALL EPC(E0_I,PF_P_I,THETA_EP,ZTARG,NTARG,6,SIGMA_EP)
             CALL EPC(E0_I,PF_P_I,THETA_EP,ZTARG,NTARG,7,SIGMA_EPIP)
           ENDIF
           SIGMA_EE   = SIGMA_EE   * JACOB_SA_E
           SIGMA_EPIM = SIGMA_EPIM * JACOB_SA_E
           SIGMA_EP   = SIGMA_EP   * JACOB_SA_P
           SIGMA_EPIP = SIGMA_EPIP * JACOB_SA_P
         ENDIF
C
C ---------------------------------------------------------------------
C          Get Coincidence cross sections.
C
C          For Version 3.5 and beyond:
C              Some of the coding above has been reordered.
C              It was originally to deal with multiple calls to KINEM
C              required by Arenhoevel interpolator.  However, a
C              separate d(e,e'p)n kinematics routine was written
C              which avoids this need.  Since the reordering resulted
C              in prettier and more straightforward code, I have
C              kept it.
C ---------------------------------------------------------------------
C
        CALL PHYSICS(ELASTIC,BOUND,SPEC_FAC,POL_BEAM,
     #          SIGMA_EEP,SIGMA_MULT_WT,ASYMMETRY,
     #          POL_N,POL_T,POL_L,POL_N_S,POL_T_S,POL_L_S,
     #          POL_X,POL_Y,POL_Z,
     #          POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #          POL_N_RHD,POL_T_RHD,POL_L_RHD,
     #          POL_N_SHI,POL_T_SHI,POL_L_SHI,
     #          POL_N_SHD,POL_T_SHD,POL_L_SHD,
     #          POL_X_HI, POL_Y_HI, POL_Z_HI,
     #          POL_X_HD, POL_Y_HD, POL_Z_HD)
C
C ---------------------------------------------------------------------
C          For computing accidentals rates:
C
C          Convert product of cross sections to yields.
C
C             For e',p   :  SIGFACTOR(50) * (e,e')  * (e,p)
C             For e',pi+ :  SIGFACTOR(51) * (e,e')  * (e,pi+)
C             For pi-,p  :  SIGFACTOR(52) * (e,pi-) * (e,p)
C             For pi-,pi+:  SIGFACTOR(53) * (e,pi-) * (e,pi+)
C
C          As for the other reaction types, SIGFACTOR is put in
C          after the event loop.
C
C ---------------------------------------------------------------------
C
        ACC_E_P     = SIGMA_EE   * SIGMA_EP
        ACC_E_PIP   = SIGMA_EE   * SIGMA_EPIP
        ACC_PIM_P   = SIGMA_EPIM * SIGMA_EP
        ACC_PIM_PIP = SIGMA_EPIM * SIGMA_EPIP
C
C ---------------------------------------------------------------------
C          Now, weight various physics quantities.
C ---------------------------------------------------------------------
C
        SIGMA_EEP = SIGMA_EEP * TOTAL_WT
        DO I=1,10
           SIGMA_MULT_WT(I) = SIGMA_MULT_WT(I) * TOTAL_WT
        ENDDO
        ASYMMETRY = ASYMMETRY * TOTAL_WT
        POL_N     = POL_N     * TOTAL_WT
        POL_T     = POL_T     * TOTAL_WT
        POL_L     = POL_L     * TOTAL_WT
        POL_N_RHI = POL_N_RHI * TOTAL_WT
        POL_T_RHI = POL_T_RHI * TOTAL_WT
        POL_L_RHI = POL_L_RHI * TOTAL_WT
        POL_N_RHD = POL_N_RHD * TOTAL_WT
        POL_T_RHD = POL_T_RHD * TOTAL_WT
        POL_L_RHD = POL_L_RHD * TOTAL_WT
C
        POL_X     = POL_X     * TOTAL_WT
        POL_Y     = POL_Y     * TOTAL_WT
        POL_Z     = POL_Z     * TOTAL_WT
        POL_X_HI  = POL_X_HI  * TOTAL_WT
        POL_Y_HI  = POL_Y_HI  * TOTAL_WT
        POL_Z_HI  = POL_Z_HI  * TOTAL_WT
        POL_X_HD  = POL_X_HD  * TOTAL_WT
        POL_Y_HD  = POL_Y_HD  * TOTAL_WT
        POL_Z_HD  = POL_Z_HD  * TOTAL_WT
C
        POL_N_S   = POL_N_S   * TOTAL_WT     ! These are intermediate
        POL_T_S   = POL_T_S   * TOTAL_WT     ! quantities; the actual
        POL_L_S   = POL_L_S   * TOTAL_WT     ! histo weights are
        POL_N_SHI = POL_N_SHI * TOTAL_WT     ! found by precessing
        POL_T_SHI = POL_T_SHI * TOTAL_WT     ! these.  This precession
        POL_L_SHI = POL_L_SHI * TOTAL_WT     ! can only be done after
        POL_N_SHD = POL_N_SHD * TOTAL_WT     ! the second call to
        POL_T_SHD = POL_T_SHD * TOTAL_WT     ! SPECTROMETER, since the
        POL_L_SHD = POL_L_SHD * TOTAL_WT     ! spin matrix is based on
                                             ! asymptotic quantities.
C
C ---------------------------------------------------------------------
C          Restore kinematics to begin calculation of variables for
C          histogramming.
C ---------------------------------------------------------------------
C
        E0_I   = E0_N                          !Restore kinematics
        PH_B_I = PH_B_N
        TH_B_I = TH_B_N
C
        DO I=1,6                 ! Necessary even if SPECTROMETER
          QI_E(I) = QI_E_SAV(I)  ! is not called, since radiation
          QI_P(I) = QI_P_SAV(I)  ! sets PF_E_I/PF_P_I to vertex values!
        ENDDO
C
C ---------------------------------------------------------------------
C          Perform energy loss calculation and multiple scattering
C          for outgoing particles.  For now, what follows assumes that
C          all the multiple scattering happens at the vertex.
C ---------------------------------------------------------------------
C
        IF(ELOSS_CALC) THEN
C
          CALL EXT_BREMSS(ETHO,PF_E_I,ERAD,I_WALL_E)    ! electron
          CALL ELOSS_E(ETHO,PF_E_I,2,ELOSS,1,I_WALL_E)
          CALL MULTI_SCANG(ETHO,PF_E_I,EVARPHI,EVARTH,M_ELEC,1)
          IF(EOPT .EQ. 0) THEN
Cpeu           ELOSS_E_I = ERAD + ELOSS + ELWALK(2)  ! energy loss
Cpeu           Walls are already included in eloss_e!
             ELOSS_E_I = ERAD + ELOSS  ! energy loss
          ELSE
Cpeu           ELOSS_E_I = ERAD + ELOSS + ELWALK(3)  ! energy loss
Cpeu           Walls are already included in eloss_e!
             ELOSS_E_I = ERAD + ELOSS  ! energy loss
         ENDIF
          PF_E_I = PF_E_I - ELOSS_E_I
          IF(PF_E_I .LE. 0.D0) THEN   ! Protect against rare events
              FAIL_KIN = .TRUE.        ! with huge energy loss (ERAD)
              GOTO 999
          ENDIF
          PH_E_I = PH_E_I + EVARPHI
          TH_E_I = TH_E_I + EVARTH
C
          IF(EJECT_CHG .NE. 0) THEN      ! Don't calculate for neutrals
            P_ENER = SQRT(PF_P_I**2 + EJECT_MASS**2)  ! hadron total energy
            BETA_H_ASM = PF_P_I/P_ENER
            CALL ELOSS_H(HTHO,BETA_H_ASM,ELOSS,EJECT_CHG,1,I_WALL_H)
            CALL MULTI_SCANG(HTHO,PF_P_I,HVARPHI,HVARTH,EJECT_MASS,
     #                       EJECT_CHG)
            IF(HOPT .EQ. 0) THEN
Cpeu              ELOSS_P_I = ELOSS + ELWALK(4)  ! energy loss
Cpeu              Walls are already included in eloss_h!
              ELOSS_P_I = ELOSS  ! energy loss
            ELSE
Cpeu              ELOSS_P_I = ELOSS + ELWALK(5)  ! energy loss
Cpeu              Walls are already included in eloss_h!
              ELOSS_P_I = ELOSS  ! energy loss
            ENDIF
            P_ENER = P_ENER - ELOSS_P_I
            TEST_RADICAL = P_ENER**2 - EJECT_MASS**2
            IF(TEST_RADICAL .LT. 0.D0) THEN  !Protect against rare
               FAIL_KIN = .TRUE.             ! events with huge eloss
               GOTO 999
            ELSE
               PF_P_I = SQRT(TEST_RADICAL)
            ENDIF
            PH_P_I = PH_P_I + HVARPHI
            TH_P_I = TH_P_I + HVARTH
          ENDIF
C
        ENDIF
C
C ---------------------------------------------------------------------
C          Calculate Transport coordinates for each arm (at target).
C ---------------------------------------------------------------------
C
        IF(CALLSPEC) THEN
          CALL LAB_TO_TRPT(QI_E,ROT_E,PF_E,TRPT_T_E)
          CALL LAB_TO_TRPT(QI_P,ROT_P,PF_P,TRPT_T_P)
C
C ---------------------------------------------------------------------
C          Make histogram variables of initial Transport coordinates.
C ---------------------------------------------------------------------
C
          ICNT = 61      !First index of array VAR
          DO I=1,6
             VAR(ICNT) = TRPT_T_E(I)    !Electron arm
             ICNT = ICNT + 1
          ENDDO
          ICNT = 71
          DO I=1,6
             VAR(ICNT) = TRPT_T_P(I)    !Hadron arm
             ICNT = ICNT + 1
          ENDDO
C
C ---------------------------------------------------------------------
C          Do spectrometer analysis.
C ---------------------------------------------------------------------
C
          CALL SPECTROMETER(2)
          TOF_REL = TOF_REL_TMP   !Relative Time-of-Flight (nsec)
C
C ---------------------------------------------------------------------
C          Make histogram variables of intermediate Hall A-HRS
C          coords.
C ---------------------------------------------------------------------
C
          var(130) = eq1_ext_x
          var(131) = eq1_ext_y
          var(132) = edip_ent_x
          var(133) = edip_ent_y
          var(134) = edip_ext_x
          var(135) = edip_ext_y
          var(136) = eq3_ent_x
          var(137) = eq3_ent_y
          var(138) = eq3_ext_x
          var(139) = eq3_ext_y
c
          var(140) = hq1_ext_x
          var(141) = hq1_ext_y
          var(142) = hdip_ent_x
          var(143) = hdip_ent_y
          var(144) = hdip_ext_x
          var(145) = hdip_ext_y
          var(146) = hq3_ent_x
          var(147) = hq3_ent_y
          var(148) = hq3_ext_x
          var(149) = hq3_ext_y
C
C ---------------------------------------------------------------------
C          Make histogram variables of R-function variables
C          for Hall A HRSE and HRSH. 
C ---------------------------------------------------------------------
C
          var(150) = e_rfn
          var(151) = p_rfn
C
C ---------------------------------------------------------------------
C          Make histogram variables of coordinates at intermediate
C          locations in the MAD Hall A spectrometer.
C          Note:  not all these are filled for all MAD configurations.
C ---------------------------------------------------------------------
C
          var(152) = mag1_ent_x
          var(153) = mag1_ent_y
          var(154) = mag1_mid_x
          var(155) = mag1_mid_y
          var(156) = mag1_ext_x
          var(157) = mag1_ext_y
          var(158) = mag2_ent_x
          var(159) = mag2_ent_y
          var(160) = mag2_mid_x
          var(161) = mag2_mid_y
          var(162) = mag2_ext_x
          var(163) = mag2_ext_y
C
C ---------------------------------------------------------------------
C          Make histogram variables of intermediate Septum-HRS
C          coords.
C ---------------------------------------------------------------------
C 
          var(166) = sl_ep3_x
          var(167) = sl_ep3_y
          var(168) = sl_ep4_x
          var(169) = sl_ep4_y
          var(170) = sl_ep5_x
          var(171) = sl_ep5_y
          var(172) = sl_ep6_x
          var(173) = sl_ep6_y
          var(174) = sl_ep7_x
          var(175) = sl_ep7_y
C
          var(176) = sr_ep3_x
          var(177) = sr_ep3_y
          var(178) = sr_ep4_x
          var(179) = sr_ep4_y
          var(180) = sr_ep5_x
          var(181) = sr_ep5_y
          var(182) = sr_ep6_x
          var(183) = sr_ep6_y
          var(184) = sr_ep7_x
          var(185) = sr_ep7_y
C
C ---------------------------------------------------------------------
C          Make histogram variables of focal plane Transport coords.
C ---------------------------------------------------------------------
C
          ICNT = 81
          DO I=1,6
             VAR(ICNT) = FP_VEC_E(I)    !Electron arm
             ICNT = ICNT + 1
          ENDDO
          ICNT = 91
          DO I=1,6
             VAR(ICNT) = FP_VEC_P(I)    !Hadron arm
             ICNT = ICNT + 1
          ENDDO
C
C ---------------------------------------------------------------------
C          Make histogram variables of final Transport coordinates.
C          These coordinates are not necessarily at the FP, but
C          are at the final location defined by the entire spectrometer
C          element stack (usually these are therefore the traced-back
C          target coordinates).
C ---------------------------------------------------------------------
C
          ICNT = 101
          DO I=1,6
             VAR(ICNT) = TRPT_T_E(I)    !Electron arm
             ICNT = ICNT + 1
          ENDDO
          ICNT = 111
          DO I=1,6
             VAR(ICNT) = TRPT_T_P(I)    !Hadron arm
             ICNT = ICNT + 1
          ENDDO
C
C ---------------------------------------------------------------------
C          Calculate Laboratory coordinates from (modified) target
C          coordinates for each arm.
C ---------------------------------------------------------------------
C
          IF(APPLY_TO_LAB(1))
     #          CALL TRPT_TO_LAB(TRPT_T_E,ROTI_E,PF_E,QI_E)
          IF(APPLY_TO_LAB(2))
     #          CALL TRPT_TO_LAB(TRPT_T_P,ROTI_P,PF_P,QI_P)
C
        ENDIF
C       
C ---------------------------------------------------------------------
C          Perform spin precession for hadron arm.
C
C          The spin precession matrix PMATRIX is a 3x3 rotation
C          matrix; if spin precession is not requested in the
C          input file then this matrix is the unit matrix.  Note
C          that actually POL, POL_SHI and POL_SHD are 3-vectors,
C          the components of which represent the partial
C          cross sections for each of the three polarization
C          components.  Since these partial cross sections are simply
C          the polarizations times the unpolarized cross section
C          these vectors have the same orientation as the polarization
C          vectors; hence, it is valid to perform the rotation on
C          each vector of partial cross sections.
C ---------------------------------------------------------------------
C
        CALL MAT_MULT(3,3,1,PMATRIX,POL,POL_SPEC)
        CALL MAT_MULT(3,3,1,PMATRIX,POL_SHI,POL_SPEC_HI)
        CALL MAT_MULT(3,3,1,PMATRIX,POL_SHD,POL_SPEC_HD)
C
C ---------------------------------------------------------------------
C          Assign values for spin matrix histograms.
C ---------------------------------------------------------------------
C
        ICNT = 51
        DO I=1,3
           DO J=1,3
              VAR(ICNT) = PMATRIX(I,J)
              ICNT = ICNT + 1
           ENDDO
        ENDDO
C
C ---------------------------------------------------------------------
C          If energy loss/multiple scattering calculation was performed
C          enforce actual acceptances here.  Note that the "before cuts"
C          counters are only incremented for events within the
C          actual acceptances.
C ---------------------------------------------------------------------
C
        IF(ELOSS_CALC) THEN
C
           EF_MIN_N = PF_E*(1.+0.01*ACC_EM_N)
           EF_MAX_N = PF_E*(1.+0.01*ACC_EP_N)
           PH_E_MIN_N = -DPH_E_N/2.D0
           TH_E_MIN_N = -DTH_E_N/2.D0
C
           PF_MIN_N = PF_P*(1.+0.01*ACC_PM_N)
           PF_MAX_N = PF_P*(1.+0.01*ACC_PP_N)
           PH_P_MIN_N = -DPH_P_N/2.D0
           TH_P_MIN_N = -DTH_P_N/2.D0
C
           IF((PF_E_I .LT. EF_MIN_N)
     #             .OR. (PF_E_I .GT. EF_MAX_N))  GOTO 999
           IF(ACCEPT_CHECK .AND. ((PF_P_I .LT. PF_MIN_N)
     #             .OR. (PF_P_I .GT. PF_MAX_N))) GOTO 999
C
           CALL ACT_TO_NOM(PH_E_I,TH_E_I,DRIFTE_N,BEAM_E1,
     #                      BEAM_E2,BEAM_E3,PH_E_N,TH_E_N)
           PH_E_N =  (PH_E_N-PH_E)*COS(TH_E_N) !rel. to aperture center
           TH_E_N =  TH_E_N - TH_E             !rel. to aperture center
           IF(.NOT. IN_APERTURE(SA_SHAPE_E,
     #            PH_E_MIN_N,TH_E_MIN_N,PH_E_N,TH_E_N) ) GOTO 999
C
           CALL ACT_TO_NOM(PH_P_I,TH_P_I,DRIFTP_N,BEAM_P1,
     #                      BEAM_P2,BEAM_P3,PH_P_N,TH_P_N)
           PH_P_N =  (PH_P_N-PH_P)*COS(TH_P_N) !rel. to aperture center
           TH_P_N =  TH_P_N - TH_P             !rel. to aperture center
           IF(ACCEPT_CHECK .AND. .NOT. IN_APERTURE(SA_SHAPE_P,
     #            PH_P_MIN_N,TH_P_MIN_N,PH_P_N,TH_P_N) ) GOTO 999
C
        ENDIF
C
C ---------------------------------------------------------------------
C          Calculate derived kinematical quantities to determine values
C          of variables for incrementation.
C ---------------------------------------------------------------------
C
        IF(BEAMVAR.OR.CALLSPEC.OR.RADPEAKING.OR.ELOSS_CALC
     #                        .OR.RADFULL) THEN
C
C            Perform energy loss correction, if selected.  This is for
C            comparison to data analyzed with this correction.
C                     L. Todor - October 7, 1999
C
          IF(DE_DX_CORR) THEN
             DE_DX_CORR_B = CELOSS(1)      ! For histograms
             DE_DX_CORR_E = CELOSS(2)
             DE_DX_CORR_P = CELOSS(3)
          ELSEIF(ELOSS_MP_CORR) THEN
             DE_DX_CORR_B = ELOSS_MP(1)    ! For histograms
             DE_DX_CORR_E = ELOSS_MP(2)
             DE_DX_CORR_P = ELOSS_MP(3)
          ELSE
             DE_DX_CORR_B = 0.             ! For histograms
             DE_DX_CORR_E = 0.
             DE_DX_CORR_P = 0.
          ENDIF
C
          IF (DE_DX_CORR .OR. ELOSS_MP_CORR) THEN
             E0_I   = E0_I   - DE_DX_CORR_B
             PF_E_I = PF_E_I + DE_DX_CORR_E
             IF(EJECT_CHG .NE. 0) THEN   ! Not needed for neutrals
               PF_P_I = sqrt((sqrt(PF_P_I**2+EJECT_MASS**2)+
     #            DE_DX_CORR_P)**2 - EJECT_MASS**2)
             ENDIF
          ENDIF
                                                    !Only redo if necessary
          CALL V3MAKE(E0_I,PH_B_I,TH_B_I,KI_VEC)    !make beam 3-vector
          CALL V3MAKE(PF_E_I,PH_E_I,TH_E_I,KF_VEC)  !make scatt e- 3-vector
C
          CALL KINEM(E0_I,KI_VEC,PF_E_I,KF_VEC,PF_P_I,PH_P_I,TH_P_I,
     #              .FALSE.,.FALSE.,.FALSE.,FAIL_KIN)
          IF(FAIL_KIN) GOTO 999
        ENDIF
C
C ---------------------------------------------------------------------
C          Increment "before cuts" counters.
C ---------------------------------------------------------------------
C
        SUM_EV  = SUM_EV  + 1.0D0     !Sum of tries before cuts
        SUM_PS  = SUM_PS  + PHASE_SP  !Sum of phase space before cuts
        SUM_EEP = SUM_EEP + SIGMA_EEP !Sum of sigma_eep before cuts
C
C ---------------------------------------------------------------------
C          Check global cuts.
C ---------------------------------------------------------------------
C
        IF(FAIL_GLOBAL_CUT) GOTO 999            !Passed Transport cuts?
        IF(NCUTS .GT. 0) THEN
          DO I=1,NCUTS                            !Check regular histo cuts
            IF((VAR(ICUT_VAR(I)).LT.CUT_MIN(I))
     #             .OR. (VAR(ICUT_VAR(I)).GT.CUT_MAX(I))) GOTO 999
          ENDDO
        ENDIF
C
C ---------------------------------------------------------------------
C          Get min., max., avg. (unweighted) and centroid for each
C          variable (for summary file).
C ---------------------------------------------------------------------
C
        DO I=1,NUM_VAR
          VAR_CENT(I) = VAR_CENT(I) + VAR(I)*SIGMA_EEP
          VAR_AVG(I)  = VAR_AVG(I)  + VAR(I)
          IF(VAR(I) .GT. VAR_MAX(I)) VAR_MAX(I) = VAR(I)
          IF(VAR(I) .LT. VAR_MIN(I)) VAR_MIN(I) = VAR(I)
        ENDDO
C
C ---------------------------------------------------------------------
C          Increment sums for summary file.
C ---------------------------------------------------------------------
C
        SUM_EV_C  = SUM_EV_C  + 1.0D0     !Sum of tries after cuts
        SUM_PS_C  = SUM_PS_C  + PHASE_SP  !Sum of phase space after cuts
        SUM_EEP_C = SUM_EEP_C + SIGMA_EEP !Sum of sigma_eep after cuts
        IF(FAIL_CELL) THEN
           I_AREN_FAIL = I_AREN_FAIL + 1  !Only increment if passed cuts
        ENDIF
C
C ---------------------------------------------------------------------
C          Test for out-of-bounds and histogram specific cuts
C          and increment TRANSPORT histograms.
C ---------------------------------------------------------------------
C
        ICOUNT_NTU = 0               !Reset
        DO IARM=1,2                  !TRANSPORT histograms
          IF(NEL(IARM) .GT. 0) THEN
            DO IEL=1,NEL(IARM)
              IF(.NOT.TR_INC_FLAG(IARM,IEL)) GOTO 997   !Check cuts
C
              IF(OP(IARM,IEL).EQ.'H1D') THEN
                IXCHAN_TR = NINT((TR_VAR(IARM,IEL,1) -
     #            TR_MINX(IARM,IEL))/(TR_MAXX(IARM,IEL) -
     #            TR_MINX(IARM,IEL))*DFLOAT(NXCHAN_TR(IARM,IEL))+0.5D0)
                IF(IXCHAN_TR .LT. 1 .OR.
     #           IXCHAN_TR .GT. NXCHAN_TR(IARM,IEL)) THEN
                  I_OOB_TR(IARM,IEL) = I_OOB_TR(IARM,IEL) + 1
                ELSE
                  TRHIST(IARM,IEL,IXCHAN_TR)=TRHIST(IARM,IEL,IXCHAN_TR)
     #                  + SIGMA(ABS(IPID_TR(IARM,IEL)))
                  IF(IPID_TR(IARM,IEL) .LT. 0) THEN
                    TRHIST_DEN(IARM,IEL,IXCHAN_TR) =
     #                         TRHIST_DEN(IARM,IEL,IXCHAN_TR)
     #                         + SIGMA(IPID_TR_DEN(IARM,IEL))
                  ENDIF
                ENDIF
C
              ELSEIF(OP(IARM,IEL).EQ.'H2D')THEN
                IXCHAN_TR = NINT((TR_VAR(IARM,IEL,1) -
     #             TR_MINX(IARM,IEL))/(TR_MAXX(IARM,IEL) -
     #             TR_MINX(IARM,IEL))*DFLOAT(NXCHAN_TR(IARM,IEL))+0.5D0)
                IYCHAN_TR = NINT((TR_VAR(IARM,IEL,2) -
     #             TR_MINY(IARM,IEL))/(TR_MAXY(IARM,IEL) -
     #             TR_MINY(IARM,IEL))*DFLOAT(NYCHAN_TR(IARM,IEL))+0.5D0)
                IF(IXCHAN_TR .LT. 1 .OR.
     #             IXCHAN_TR .GT. NXCHAN_TR(IARM,IEL) .OR.
     #             IYCHAN_TR .LT. 1 .OR.
     #             IYCHAN_TR .GT. NYCHAN_TR(IARM,IEL)) THEN
                  I_OOB_TR(IARM,IEL) = I_OOB_TR(IARM,IEL) + 1
                ELSE
                  TRHIST_2D(IARM,IEL,IXCHAN_TR,IYCHAN_TR)
     #                  = TRHIST_2D(IARM,IEL,IXCHAN_TR,IYCHAN_TR)
     #                  + SIGMA(ABS(IPID_TR(IARM,IEL)))
                  IF(IPID_TR(IARM,IEL) .LT. 0) THEN
                    TRHIST_DEN_2D(IARM,IEL,IXCHAN_TR,IYCHAN_TR)
     #                  = TRHIST_DEN_2D(IARM,IEL,IXCHAN_TR,IYCHAN_TR)
     #                       + SIGMA(IPID_TR_DEN(IARM,IEL))
                  ENDIF
                ENDIF
C
              ELSEIF(OP(IARM,IEL).EQ.'SCT')THEN
                IF(TR_VAR(IARM,IEL,1).LT.TR_MINX(IARM,IEL) .OR.
     #             TR_VAR(IARM,IEL,1).GT.TR_MAXX(IARM,IEL) .OR.
     #             TR_VAR(IARM,IEL,2).LT.TR_MINY(IARM,IEL) .OR.
     #             TR_VAR(IARM,IEL,2).GT.TR_MAXY(IARM,IEL)) THEN
                     I_OOB_TR(IARM,IEL) = I_OOB_TR(IARM,IEL) + 1
                ELSE
                   WRITE(LUN_SCT(IARM,IEL),1111) TR_VAR(IARM,IEL,1),
     #                          TR_VAR(IARM,IEL,2)
                ENDIF
                IF(DFLOAT(IEVNT)/100. .EQ. DFLOAT(IEVNT/100))
     #              WRITE(LUN_SCT(IARM,IEL),'(''   plot '')')
C
              ELSEIF(OP(IARM,IEL).EQ.'NTU')THEN
                ICOUNT_NTU = ICOUNT_NTU + 1
                IF(IPID_TR(IARM,IEL) .LT. 0) THEN
                   WRITE(LUN_NTU(ICOUNT_NTU),1112)
     #                 SIGMA(ABS(IPID_TR(IARM,IEL))),
     #                 SIGMA(IPID_TR_DEN(IARM,IEL))
                   WRITE(LUN_NTU(ICOUNT_NTU),1112)
     #                      (TR_VAR(IARM,IEL,J),J=1,6)
1112               FORMAT(100(1X,1PE12.5))
                ELSE
                   WRITE(LUN_NTU(ICOUNT_NTU)+10,1112)
     #             SIGMA(IPID_TR(IARM,IEL))*SIGFACTOR(IPID_TR(IARM,IEL))
                   WRITE(LUN_NTU(ICOUNT_NTU)+10,1112)
     #                    (TR_VAR(IARM,IEL,J),J=1,6)
                ENDIF
C
              ENDIF
  997       ENDDO
          ENDIF
        ENDDO
C
C ---------------------------------------------------------------------
C          Test for out-of-bounds and histogram specific cuts
C          and increment regular 1-D histograms.
C ---------------------------------------------------------------------
C
        ICOUNT_NTM = 0               !Reset
        IF(NPLOTS .GT. 0) THEN
          DO I=1,NPLOTS
            IF(PLOT_TYPE(I).NE.'NTU'.AND.PLOT_TYPE(I).NE.'NTM') THEN
               IF(X_MIN(I).EQ.X_MAX(I)) GOTO 998
               IF(Y_MIN(I).EQ.Y_MAX(I) .AND. TWO_D(I)) GOTO 998
               IF(NCUTS_PL(I) .GT. 0) THEN    ! Check cuts
                 DO J=1,NCUTS_PL(I)
                   IF(VAR(ICUT_VAR_I(ICUT_IND_PL(I,J))) .LT.
     #                       CUT_MIN_I(ICUT_IND_PL(I,J))
     #                  .OR.VAR(ICUT_VAR_I(ICUT_IND_PL(I,J))) .GT.
     #                       CUT_MAX_I(ICUT_IND_PL(I,J)))
     #                          GOTO 998
                 ENDDO
               ENDIF
               XVALUE = VAR(I_VAR(I,1))*X_SC(I)+X_OFF(I)
               IF(TWO_D(I)) YVALUE = VAR(I_VAR(I,2))*Y_SC(I)+Y_OFF(I)
            ENDIF
C
            IF(PLOT_TYPE(I) .EQ. 'P1D') THEN
              IVARX = NINT((XVALUE-X_MIN(I))
     #            /(X_MAX(I)-X_MIN(I))*DFLOAT(NX_CHAN(I))+0.5D0)
              IF(IVARX .LT. 1 .OR. IVARX .GT. NX_CHAN(I)) THEN
                I_OOB(I) = I_OOB(I) + 1
              ELSE
                SIG_ARR(I,IVARX)=SIG_ARR(I,IVARX) + SIGMA(ABS(IPID(I)))
                IF(IPID(I) .LT. 0) THEN
                  SIG_DEN(I,IVARX)=SIG_DEN(I,IVARX) + SIGMA(IPID_DEN(I))
                ENDIF
              ENDIF
C
            ELSEIF(PLOT_TYPE(I) .EQ. 'P2D') THEN
              IVARX = NINT((XVALUE-X_MIN(I))
     #            /(X_MAX(I)-X_MIN(I))*DFLOAT(NX_CHAN(I))+0.5D0)
              IVARY = NINT((YVALUE-Y_MIN(I))
     #              /(Y_MAX(I)-Y_MIN(I))*DFLOAT(NY_CHAN(I))+0.5D0)
              IF(IVARX .LT. 1 .OR. IVARX .GT. NX_CHAN(I) .OR.
     #           IVARY .LT. 1 .OR. IVARY .GT. NY_CHAN(I)) THEN
                I_OOB(I) = I_OOB(I) + 1
              ELSE
                SIG_ARR_2D(I,IVARX,IVARY)=SIG_ARR_2D(I,IVARX,IVARY)
     #                  + SIGMA(ABS(IPID(I)))
                IF(IPID(I) .LT. 0) THEN
                  SIG_DEN_2D(I,IVARX,IVARY)=SIG_DEN_2D(I,IVARX,IVARY)
     #                  + SIGMA(IPID_DEN(I))
                ENDIF
              ENDIF
C
            ELSEIF(PLOT_TYPE(I) .EQ. 'SCA') THEN
              IF(XVALUE .LT. X_MIN(I) .OR. XVALUE .GT. X_MAX(I) .OR.
     #           YVALUE .LT. Y_MIN(I) .OR. YVALUE .GT. Y_MAX(I)) THEN
                I_OOB(I) = I_OOB(I) + 1
              ELSE
                WRITE(LUN(I),1111) XVALUE,YVALUE
 1111           FORMAT(1X,E14.5,3X,E14.5)
              ENDIF
              IF(DFLOAT(IEVNT)/100. .EQ. DFLOAT(IEVNT/100))
     #             WRITE(LUN(I),'(''   plot '')')
C
            ELSEIF(PLOT_TYPE(I) .EQ. 'NTU') THEN
              ICOUNT_NTU = ICOUNT_NTU + 1
              IF(IPID(I) .LT. 0) THEN
                 WRITE(LUN_NTU(ICOUNT_NTU),1112)
     #               SIGMA(ABS(IPID(I))),SIGMA(IPID_DEN(I))
                 WRITE(LUN_NTU(ICOUNT_NTU),1112)
     #               (VAR(I_VAR(I,J)),J=1,NVAR_NTU(I))
              ELSE
                 WRITE(LUN_NTU(ICOUNT_NTU)+10,1112)
     #               SIGMA(IPID(I))*SIGFACTOR(IPID(I))
                 WRITE(LUN_NTU(ICOUNT_NTU)+10,1112)
     #               (VAR(I_VAR(I,J)),J=1,NVAR_NTU(I))
              ENDIF
C
            ELSEIF(PLOT_TYPE(I) .EQ. 'NTM') THEN
              ICOUNT_NTM = ICOUNT_NTM + 1
              WRITE(LUN_NTM(ICOUNT_NTM)+10,1112)
     #            (SIGMA(I_WT(I,J))
     #                *SIGFACTOR(I_WT(I,J)),J=1,NWT_NTM(I))
              WRITE(LUN_NTM(ICOUNT_NTM)+10,1112)
     #            (VAR(I_VAR(I,J)),J=1,NVAR_NTU(I))
C
            ENDIF
  998     ENDDO
        ENDIF
C
  999 ENDDO         ! End of the event loop
C
      RADPEAKING = .FALSE.   ! Next calculation is for "peak"
      ENDDO         ! End of the ICALC loop
C
C _____________________________________________________________________
C ---------------------------------------------------------------------
C       End of the event loop.
C _____________________________________________________________________
C
C ---------------------------------------------------------------------
C       Close Acceptance Function Ntuple
C ---------------------------------------------------------------------
C
      IF(ACCEPT_FCN) THEN
         CLOSE(UNIT=41)
         WRITE(6,9900)
 9900    FORMAT(' Output:  accept_fcn.ntu ')
         GOTO 9999           ! All done
      ENDIF
C
C ---------------------------------------------------------------------
C       Close 2-dimensional scatter plot files and N-Tuple files.
C ---------------------------------------------------------------------
C
      IF(ICOUNT_SCT .GT. 0) THEN
        DO I=1,ICOUNT_SCT                       !TRANSPORT scatter plots
          WRITE(I+30,'(''   plot '')')
          CLOSE(UNIT=I+30)
        ENDDO
      ENDIF
C
      IF(NCOUNT_NTU .GT. 0) THEN
        DO I=1,NCOUNT_NTU                       !N-Tuple files ('NTU')
          IF(WRITE_TMP_NTU(I)) THEN
            REWIND(UNIT=I+60)
            DO J=1,3
              READ(I+60,30) CDUM
            ENDDO
            DO J=1,N_EVENT*NCALC
              READ(I+60,*,END=9990) XX
              WRITE(I+50,1112) XX/SUM_ACC
              READ(I+60,*) (YY(K),K=1,NUM_VAR_NTU(I))
              WRITE(I+50,1112) (YY(K),K=1,NUM_VAR_NTU(I))
            ENDDO
9990        CLOSE(UNIT=I+50)
            CLOSE(UNIT=I+60)
          ELSE
            CLOSE(UNIT=I+50)
          ENDIF
        ENDDO
      ENDIF
C
      IF(NCOUNT_NTM .GT. 0) THEN
        DO I=1,NCOUNT_NTM                       !N-Tuple files ('NTM')
           REWIND(UNIT=I+80)
           DO J=1,3
             READ(I+80,30) CDUM
           ENDDO
           DO J=1,N_EVENT*NCALC
             READ(I+80,*,END=9991) (WW(K),K=1,NUM_WT_NTM(I))
             WRITE(I+70,1112) (WW(K)/SUM_ACC,K=1,NUM_WT_NTM(I))
             READ(I+80,*) (YY(K),K=1,NUM_VAR_NTM(I))
             WRITE(I+70,1112) (YY(K),K=1,NUM_VAR_NTM(I))
           ENDDO
 9991      CLOSE(UNIT=I+70)
           CLOSE(UNIT=I+80)
        ENDDO
      ENDIF
C
      IF(NPLOTS .GT. 0) THEN
        DO I=1,NPLOTS
          IF(PLOT_TYPE(I) .EQ. 'SCA') THEN
            WRITE(LUN(I),'(''   plot '')')
            CLOSE(UNIT=LUN(I))
          ENDIF
        ENDDO
      ENDIF
C
C ---------------------------------------------------------------------
C       Create normalized yield or form ratio arrays (for average
C       cross section and polarization histograms).
C       Coincidence cross section histograms are in fm^2 per whatever.
C       Singles cross section histograms are in nb/MeV/c/sr.
C
C       The cross sections are averages for a given bin.  That is,
C       each time an event falls within the overall experimental
C       acceptance, an event counter is incremented and a running sum
C       of cross sections is made.  The average cross sections are
C       computed by dividing the sum for a given bin by the value of the
C       event counter (as given by the phase space array).  This
C       is, in contrast to the calculation of yields.
C       Here, the event counter is incremented for any event which falls
C       within the acceptance of those variables in which the cross
C       section is differential (SUM_ACC).  For elastic scattering,
C       for example, where some events which contribute to the electron
C       solid angle result in kinematics outside of the other
C       acceptances, the yield is not simply approximately equal to the
C       average cross section times the electron solid angle (times
C       luminosity times the beam time).
C ---------------------------------------------------------------------
C
      DO IARM=1,2                     !TRANSPORT histograms
        IF(NEL(IARM) .GT. 0) THEN
          DO IEL=1,NEL(IARM)
            IF(OP(IARM,IEL).EQ.'H1D')THEN
              IF(IPID_TR(IARM,IEL) .GE. 0) THEN !Only normalize YIELD histos
                DO J=1,NXCHAN_TR(IARM,IEL)
                  TRHIST(IARM,IEL,J)=TRHIST(IARM,IEL,J)
     #                    *SIGFACTOR(IPID_TR(IARM,IEL))/SUM_ACC
                ENDDO
              ELSE
                DO J=1,NXCHAN_TR(IARM,IEL)
                  IF(TRHIST_DEN(IARM,IEL,J) .NE. 0.) THEN
                    TRHIST(IARM,IEL,J) = TRHIST(IARM,IEL,J)/
     #                            TRHIST_DEN(IARM,IEL,J)
                  ELSE
                    TRHIST(IARM,IEL,J) = 0.D0
                  ENDIF
                ENDDO
              ENDIF
C
            ELSEIF(OP(IARM,IEL).EQ.'H2D')THEN
              IF(IPID_TR(IARM,IEL) .GE. 0) THEN !Only normalize YIELD histos
                DO J=1,NXCHAN_TR(IARM,IEL)
                DO K=1,NYCHAN_TR(IARM,IEL)
                  TRHIST_2D(IARM,IEL,J,K)=TRHIST_2D(IARM,IEL,J,K)
     #                    *SIGFACTOR(IPID_TR(IARM,IEL))/SUM_ACC
                ENDDO
                ENDDO
              ELSE
                DO J=1,NXCHAN_TR(IARM,IEL)
                DO K=1,NYCHAN_TR(IARM,IEL)
                  IF(TRHIST_DEN_2D(IARM,IEL,J,K) .NE. 0.) THEN
                    TRHIST_2D(IARM,IEL,J,K) = TRHIST_2D(IARM,IEL,J,K)/
     #                            TRHIST_DEN_2D(IARM,IEL,J,K)
                  ELSE
                    TRHIST_2D(IARM,IEL,J,K) = 0.D0
                  ENDIF
                ENDDO
                ENDDO
              ENDIF
C
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
      IF(NPLOTS .GT. 0) THEN
        DO I=1,NPLOTS                    !Regular kinematics histograms
          IF(PLOT_TYPE(I) .EQ. 'P1D') THEN
            IF(IPID(I) .GE. 0) THEN      !Only normalize YIELD histos
              DO J=1,NX_CHAN(I)
                SIG_ARR(I,J) = SIG_ARR(I,J)*SIGFACTOR(IPID(I))/SUM_ACC
              ENDDO
            ELSE
              DO J=1,NX_CHAN(I)
                IF(SIG_DEN(I,J) .NE. 0.) THEN
                  SIG_ARR(I,J) = SIG_ARR(I,J)/SIG_DEN(I,J)
                ELSE
                  SIG_ARR(I,J) = 0.D0
                ENDIF
              ENDDO
            ENDIF
C
          ELSEIF(PLOT_TYPE(I) .EQ. 'P2D') THEN
            IF(IPID(I) .GE. 0) THEN      !Only normalize YIELD histos
              DO J=1,NX_CHAN(I)
              DO K=1,NY_CHAN(I)
                SIG_ARR_2D(I,J,K) = SIG_ARR_2D(I,J,K)*SIGFACTOR(IPID(I))
     #                            /SUM_ACC
              ENDDO
              ENDDO
            ELSE
              DO J=1,NX_CHAN(I)
              DO K=1,NY_CHAN(I)
                IF(SIG_DEN_2D(I,J,K) .NE. 0.) THEN
                 SIG_ARR_2D(I,J,K) = SIG_ARR_2D(I,J,K)/SIG_DEN_2D(I,J,K)
                ELSE
                  SIG_ARR_2D(I,J,K) = 0.D0
                ENDIF
              ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDDO
      ENDIF
C
C ---------------------------------------------------------------------
C       Write the arrays.
C ---------------------------------------------------------------------
C
      DO IARM=1,2                     !TRANSPORT histograms
        IF(NEL(IARM) .GT. 0) THEN
          DO IEL=1,NEL(IARM)
            IF(OP(IARM,IEL).EQ.'H1D') THEN
              DO J=1,NXCHAN_TR(IARM,IEL)
                ARR(J) = TRHIST(IARM,IEL,J)
              ENDDO
              CALL TOP_1D(TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL),
     #              NXCHAN_TR(IARM,IEL),I_OOB_TR(IARM,IEL),
     #              TR_FILE(IARM,IEL),REACT_NAME(IPID_TR(IARM,IEL)),
     #              TR_AXIS(NCOORD(IARM,IEL,1)),ARR,0.,.FALSE.)
C
            ELSEIF(OP(IARM,IEL).EQ.'H2D') THEN
              DO J=1,NXCHAN_TR(IARM,IEL)
              DO K=1,NYCHAN_TR(IARM,IEL)
                ARR_2D(J,K) = TRHIST_2D(IARM,IEL,J,K)
              ENDDO
              ENDDO
              CALL TOP_2D(TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL),
     #                    TR_MINY(IARM,IEL),TR_MAXY(IARM,IEL),
     #                    NXCHAN_TR(IARM,IEL),NYCHAN_TR(IARM,IEL),
     #                    TR_FILE(IARM,IEL),ARR_2D)
C
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
      IF(NPLOTS .GT. 0) THEN
        DO I=1,NPLOTS
          IF(PLOT_TYPE(I) .EQ. 'P1D') THEN
            DO J=1,NX_CHAN(I)
              ARR(J) = SIG_ARR(I,J)
            ENDDO
            CALL TOP_1D(X_MIN(I),X_MAX(I),NX_CHAN(I),I_OOB(I),
     #                  PLOT_FIL(I),REACT_NAME(IPID(I)),
     #                  VAR_NAME(I_VAR(I,1)),ARR,0.,.FALSE.)
C
          ELSEIF(PLOT_TYPE(I) .EQ. 'P2D') THEN
            DO J=1,NX_CHAN(I)
            DO K=1,NY_CHAN(I)
              ARR_2D(J,K) = SIG_ARR_2D(I,J,K)
            ENDDO
            ENDDO
            CALL TOP_2D(X_MIN(I),X_MAX(I),Y_MIN(I),Y_MAX(I),
     #                  NX_CHAN(I),NY_CHAN(I),PLOT_FIL(I),ARR_2D)
C
          ENDIF
        ENDDO
      ENDIF
C
C ---------------------------------------------------------------------
C       Write the final seeds to a file.
C ---------------------------------------------------------------------
C
      CALL SAVESEED
C
C ---------------------------------------------------------------------
C       Write the summary file.
C ---------------------------------------------------------------------
C
      CALL SUMMARY(ELASTIC,BOUND,RADPEAKING_SAV,MULTIPHOTON,RADFULL,
     #             ACCEPT_CHECK)
C
 9999 WRITE(6,*) CHAR(7)
C
      STOP
      END
















