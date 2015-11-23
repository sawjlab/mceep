C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE INPUT
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    13-OCT-1990
C       PURPOSE: Open and read the input file
C
C       INPUT VARIABLES:
C
C               N_EVENT:     Number of events (i.e. tries).
C
C               NITER1:      RANGES preprocessor iterations for Theta_p
C               NITER2:      RANGES preprocessor iterations for Phi_p
C               NITER3:      RANGES preprocessor iterations for Theta_e
C               NITER4:      RANGES preprocessor iterations for Phi_e
C               NITER5:      RANGES preprocessor iterations for e- momentum
C               NITER6:      RANGES preprocessor iterations for hadron momentum
C
C               EJECT_MASS:  Mass of ejectile (MeV).
C               EJECT_CHG:   Charge of ejectile (integer) - for mult. scatt.
C               EM_BOUND:    Missing mass for bound state (MeV).
C
C               E0:          Incident Energy (MeV).
C               PH_B:        Beam in-plane angle (deg.).
C               TH_B:        Beam out-of-plane angle (deg.).
C               PF_E:        Electron central momentum (MeV/c).
C               PH_E:        Electron scattering angle (deg.).
C               TH_E:        Electron out-of-plane angle (deg.).
C               PF_P:        Hadron central momentum (MeV/c).
C               PH_P:        Hadron angle w.r.t. Incident e- (deg.).
C               TH_P:        Hadron out-of-plane angle (deg.).
C
C               ACC_EP(M):   Electron spect. momentum accept. in +/- %
C               ACC_PP(M):   Hadron          momentum accept. in +/- %
C
C               SA_SHAPE_E:  Shape of e- aperture ('E'=ell., 'R'=rect.)
C               SA_SHAPE_P:  Shape of p  aperture ('E'=ell., 'R'=rect.)
C                            (Note: for elliptical apertures, the
C                                   acceptances below indicate the
C                                   lengths of each side of the
C                                   superscribed rectangle (i.e. TWICE
C                                   the semi-minor/major axis))
C               DPH_E:       Full Phi   angular acceptance for electron (mrad).
C               DTH_E:       Full Theta angular acceptance for electron (mrad).
C               DPH_P:       Full Phi   angular acceptance for hadron   (mrad).
C               DTH_P:       Full Theta angular acceptance for hadron   (mrad).
C
C               ALUM:        Average luminosity (uA-g/cm^2)
C               BEAMTIME:    Acquisition time (hours).
C               SPEC_FAC:    Spectroscopic factor (mult. factor for sigma_eeN)
C
C               PFERMI:      Fermi momentum (MeV/c).
C               EPS:         Nucleon separation energy (MeV).
C               EPSD:        Delta separation energy (MeV).
C
C               ATARG:       Atomic weight of target nucleus.
C               ZTARG:       Atomic number of target nucleus.
C               DENS_TARG:   Density of target material (g/cm^3)
C               ITARG_MOD:   Integer specifying target model.
C                            =1 --> foil(s)
C                            =2 --> JLAB Hall A H/D cryotarget (beer can)
C                            =3 --> JLAB Hall A He target (tuna can)
C                            =4 --> JLAB HALL A Polarized He3 target
C                            =5 --> JLAB HALL A H/D cryotarget (cigar tube)
C               ELOSS_MOD:   Integer determining energy loss model:
C                            =0 --> don't do energy loss calc.
C                            =1 --> do full energy loss calc.
C                            =2 --> do calculation but remove mean dE/dx
C                                      for comparison to analyzers which
C                                      perform a this correction.
C                            =3 --> do calculation but remove the most
C                                      probable energy loss
C                                      for comparison to analyzers which
C                                      perform this correction.
C                                   
C               (For ITARG_MOD = 1:)
C                        PHI_TARG:    Angle of all foils (0=normal to beam).
C                        NFOILS:      Number of foils (normally 1)
C                        (For J=1, NFOILS:)
C                        TARG_LO(J):  Foil starting position along beam (m).
C                        TARG_HI(J):  Foil ending   position along beam (m).
C               (For ITARG_MOD = 2,3,4,5:)
C                        TARG_LO(1):  Cell starting position along beam (m).
C                        TARG_HI(1):  Cell ending   position along beam (m).
C
C               DRIFTE_N:    Nominal drift to S.A. aperture for e- arm (m).
C               DRIFTP_N:    Nominal drift to S.A. aperture for p  arm (m).
C
C               POL_BEAM:    Beam polarization (range is -1 to +1).
C               BEAMV:       Vertical extent of dispersed beam (+/- mm).
C               BEAMD:       Beam momentum spread for dispersed beam (+/-%).
C               DUTY_FAC:    Duty factor in % (for calculation of accidentals)
C               DEL_TOF:     ToF window in nsec (for calc. of accidentals)
C
C               FWHM_var:    FWHM resolution width for beam quantities.
C                            (units: positions (cm), angles (mr) energy(%))
C
C               DEL_var:     Offsets for beam quantities.
C                            (units: positions (cm), angles (mr) energy(%))
C
C               RAST_SHAPE:  Tag specifying Rectangular ('R') or 
C                            Elliptical ('E') shape or simply NO raster ('N').
C                            IMPORTANT NOTE:
C                            The elliptical raster in Hall A is a sqrt()
C                            function modulated by a sine/cosine wave. It is
C                            right now hardwired in the code. In order to
C                            get a _really_ uniform elliptical shape, one
C                            has to take a saw-wave function and modulate
C                            it with a sine/cosine wave. If one wants to
C                            simulate this, one has to create a new shape
C                            and simply replace SQRT(DBLE(RAN_UNIF)) by
C                            DBLE(RAN_UNIF) in mceep.f.
C               RAST_AMP_X:  Raster size amplitude in the horizontal plane.
C                            In the rectangular case, the peak-to-peak size.
C                            In the elliptical case, it is the "diameter".
C                            (units:  meters)
C               RAST_AMP_Y:  Raster size amplitude in the vertical plane.
C                            In the rectangular case, the peak-to-peak size.
C                            In the elliptical case, it is the "diameter".
C                            (units:  meters)
C
C                  ---------  SPECTROMETER ANALYSIS INPUT SECTION ----------
C
C               C_PID:       Tag specifying electron ('E') or hadron ('P') arm.
C               APPLY_TO_LAB: Logical variable which determines whether the
C                              final Transport vector is used to compute
C                              a new set of laboratory coordinates.  This
C                              allows the user to look at focal plane
C                              coordinates (for example) without having to
C                              apply an inverse map to get back to the
C                              target.
C                                      T = apply to get new lab coords
C                                      F = don't apply
C               NEL:         Number operations to be performed for given arm.
C               THETA_BP:    Bend plane orientation of spectrometer (degrees).
C               OBJECT_PT(1): X-coord of object pt of spectr in LAB system (cm)
C               OBJECT_PT(2): Y-coord of object pt of spectr in LAB system (cm)
C               OBJECT_PT(3): Z-coord of object pt of spectr in LAB system (cm)
C
C                 For I=1,NEL:
C                     OP:  Operation to be performed:
C                          = 'HRS' --> JLAB HRS TGT-->FP including apertures.
C                          = 'SEP' --> JLAB SEP+HRS TGT-->FP w/ apertures.
C                          = 'SIM' --> JLAB HRS TGT-->FP-->TGT incl apertures.
C                          = 'MAD' --> JLAB MAD TGT-->FP including apertures.
C                          = 'HRI' --> JLAB HRS FP-->TGT.
C                          = 'SPI' --> JLAB SEP+HRS FP-->TGT.
C                          = 'MDI' --> JLAB MAD FP-->TGT.
C                          = 'RFN' --> JLAB HRS R-function aperture test.
C                          = 'MAT' --> Multiplication by Transport matrix.
C                          = 'POL' --> Mult. by Spin precession matrix.
C                          = 'COS' --> Mult. by COSY Spin precession matrix.
C                          = 'DFT' --> Obtain Trans. vector after drift elm.
C                          = 'ROT' --> Rotate Trans. vector about given axis.
C                          = 'TRK' --> Track Reconstruction
C                          = 'TOF' --> Compute relative time of flight. 
C                          = 'RES' --> Fold in resolution function.
C                          = 'OFF' --> Add offset.
C                          = 'H1D' --> 1-dimensional histogram.
C                          = 'H2D' --> 2-dimensional histogram.
C                          = 'SCT' --> 2-dimensional scatter plot.
C                          = 'NTU' --> N-Tuple file (write vector).
C                          = 'CUT' --> cut on Tranport coordinate.
C                     (For OP='SEP')
C                                    HRS_ID:  Which HRS spectrometer?
C                                        = 'E':  electron or left arm
C                                        = 'H':  hadron or right arm
C                                    APER_SEP_TEST(J):  test aperture?
C                                                   T = perform test
C                                                   F = don't perform test
C                                        J=1:  septum entrance
C                                        J=2:  1/4 through septum
C                                        J=3:  1/2 through septum
C                                        J=4:  3/4 through septum
C                                        J=5:  septum exit
C                                        J=6:  Q1 exit
C                                        J=7:  Dipole entrance
C                                        J=8:  Dipole exit
C                                        J=9:  Q3 entrance
C                                        J=10: Q3 exit
C                     (For OP='HRS')
C                                    HRS_ID:  Which HRS spectrometer?
C                                        = 'E':  electron arm
C                                        = 'H':  hadron arm
C                                    APERTURE_TEST(J):  test aperture?
C                                                   T = perform test
C                                                   F = don't perform test
C                                        J=1:  Q1 exit
C                                        J=2:  Dipole entrance
C                                        J=3:  Dipole exit
C                                        J=4:  Q3 entrance
C                                        J=5:  Q3 exit
C                     (For OP='SIM')
C                                    HRS_ID:  Which HRS spectrometer?
C                                        = 'E':  electron arm
C                                        = 'H':  hadron arm
C                                    HUT_TEST(J):  perform which hut task?
C                                                   T = perform task
C                                                   F = don't perform task
C                                        J=1:  Multiple scattering
C                                        J=2:  Wire smearing
C                                        J=3:  Collimator test
C
C                     (For OP='MAD')
C                                    MAD_CONFIG: Which configuration?
C                                        = 1: 35 degree
C                                        = 2: 12 degree
C                                        = 3: 12 degree, no quads
C                                    For MAD_CONFIG = 1:
C                                       APER_MAD_TEST(J):  test aperture?
C                                                   T = perform test
C                                                   F = don't perform test
C                                           J=1:  magnet 1 entrance
C                                           J=2:  magnet 1 exit
C                                           J=3:  magnet 2 entrance
C                                           J=4:  magnet 2 exit
C                                    For MAD_CONFIG = 2 or 3
C                                       APER_MAD_TEST(J):  test aperture?
C                                                   T = perform test
C                                                   F = don't perform test
C                                           J=1:  magnet 1 entrance
C                                           J=2:  magnet 1 middle
C                                           J=3:  magnet 1 exit
C                                           J=4:  magnet 2 entrance
C                                           J=5:  magnet 2 middle
C                                           J=6:  magnet 2 exit
C                     (For OP='HRI')
C                                    HRS_ID:  Which HRS spectrometer?
C                                        = 'E':  electron arm
C                                        = 'H':  hadron arm
C                     (For OP='SPI')
C                                    HRS_ID:  Which HRS spectrometer?
C                                        = 'E':  electron arm
C                                        = 'H':  hadron arm
C                     (For OP='MDI')
C                                    MAD_INV_CONFIG: Which configuration?
C                                        = 1: 35 degree
C                                        = 2: 12 degree
C                     (For OP='RFN')
C                                    HRS_ID:  Which HRS spectrometer?
C                                        = 'E':  electron arm
C                                        = 'H':  hadron arm
C                                    HRS_COLL:
C                                        = 'T':  6.0 msr collimator
C                                        = 'F':  no front-end collimator
C                                        = 'S':  septum + HRS
C                                    RFN_CUT:  R-function cut value
C                     (For OP='MAT')
C                                    ORDER:      Order of matrix (1 or 2).
C                                    INVERT:     Take inv. of matrix:
C                                                   T = invert
C                                                   F = don't invert
C                                    TRNSPT_FIL: Filename for Transp. matrix.
C                     (For OP='POL')
C                                    BEND_ANGLE: Central ray bend angle (in
C                                                degrees) for each element.
C                     (For OP='COS')
C                                    COSY_SPIN_FIL: File for COSY Spin matrix.
C                     (For OP='DFT') 
C                                    DRIFT_LENGTH: Drift distance in cm.
C                     (For OP='ROT')
C                                    ROT_AXIS:     Rotation axis (X,Y or Z).
C                                    ROT_ANGLE:    Angle of rotation. 
C                     (For OP='TRK')
C                                    WC_FIL: Filename for W.C. parameters.
C                     (For OP='TOF')
C                                    CRAY_PATH:  Central ray path length (m).
C                     (For OP='RES')
C                                    OB:    Object on which operation applies.
C                                    NCOORD:     Index of coordinate to affect.
C                                    RESMOD:     1 or 2
C                                    RESVAL:     For RESMOD=1: FWHM of gaussian
C                                                For RESMOD=2: X/X0
C                     (For OP='OFF')
C                                    OB:    Object on which operation applies.
C                                    NCOORD:     Index of coordinate to affect.
C                                    OFF:        Offset to add.
C                     (For OP='H1D')
C                                    IPID_TR:    Weights for histogram.
C                                    TR_MINX:    Minimum for histogram
C                                    TR_MAXX:    Maximum for histogram
C                                    NXCHAN_TR:  Number of channels
C                                    NCUTS_TR:   Number of cuts to apply
C                                    (For J=1, NCUTS_TR:)
C                                       ICUT_IND_TR:   Index for cut
C                                    NCOORD(1):  Number of Trnspt. coord. (1-6)
C                                    TR_FILE:    Filename for histogram
C                     (For OP='H2D')
C                                    IPID_TR:    Weights for histogram.
C                                    TR_MINX:    Minimum X for histogram
C                                    TR_MAXX:    Maximum X for histogram
C                                    TR_MINY:    Minimum Y for histogram
C                                    TR_MAXY:    Maximum Y for histogram
C                                    NXCHAN_TR:  Number of X channels
C                                    NYCHAN_TR:  Number of Y channels
C                                    NCUTS_TR:   Number of cuts to apply
C                                    (For J=1, NCUTS_TR:)
C                                       ICUT_IND_TR:   Index for cut
C                                    NCOORD(1):  X Trnspt. coord. # (1-6)
C                                    NCOORD(2):  Y Trnspt. coord. # (1-6)
C                                    TR_FILE:    Filename for histogram
C                     (For OP='SCT')
C                                    TR_MINX:    Minimum X for histogram
C                                    TR_MAXX:    Maximum X for histogram
C                                    TR_MINY:    Minimum Y for histogram
C                                    TR_MAXY:    Maximum Y for histogram
C                                    NCUTS_TR:   Number of cuts to apply
C                                    (For J=1, NCUTS_TR:)
C                                       ICUT_IND_TR:   Index for cut
C                                    NCOORD(1):  X Trnspt. coord. # (1-6)
C                                    NCOORD(2):  Y Trnspt. coord. # (1-6)
C                                    TR_FILE:    Filename for histogram
C                     (For OP='NTU')
C                                    IPID_TR:         Weight for N-Tuple
C                                    BEND_ANGLE_TOT:  Total bend angle
C                                                     from target (deg)
C                                    TR_FILE:         Filename for N-Tuple
C                     (For OP='CUT')
C                                    CUT_DOMAIN: (G)lobal or (S)pecific
C                                    CUT_TYPE:   (R)ectangular or (E)lliptical
C                                    (For CUT_DOMAIN='G' and CUT_TYPE='R')
C                                         TR_CUT_MIN: Minimum for cut
C                                         TR_CUT_MAX: Maximum for cut
C                                         NCOORD(1):  Coordinate to cut
C                                    (For CUT_DOMAIN='G' and CUT_TYPE='E')
C                                         TR_CUT_X:   "Semi-axis" for X
C                                         TR_CUT_Y:   "Semi-axis" for Y
C                                         NCOORD(1):  X Coordinate for cut
C                                         NCOORD(2):  Y Coordinate for cut
C                                    (For CUT_DOMAIN='S' and CUT_TYPE='R')
C                                         TR_CUT_MIN: Minimum for cut
C                                         TR_CUT_MAX: Maximum for cut
C                                         NCOORD(1):  Coordinate to cut
C                                         TR_CUT_IND: Cut index
C                                    (For CUT_DOMAIN='S' and CUT_TYPE='E')
C                                         TR_CUT_X:   "Semi-axis" for X
C                                         TR_CUT_Y:   "Semi-axis" for Y
C                                         NCOORD(1):  X Coordinate for cut
C                                         NCOORD(2):  Y Coordinate for cut
C                                         TR_CUT_IND: Cut index
C
C               The possibilities for OB are:
C                            = 'K'  --> operation affects kinematics only.
C                            = 'S'  --> operation affects sigma only.
C                            = 'B'  --> operation affects both.
C                            = 'D'  --> default.
C                                           --> affects kin. for OP='RES'
C                                           --> affects sig. for OP='OFF'
C
C               (FWHM and OFF have units of percent for delta, mr for angles
C               and cm for positions (standard Transport units.)
C
C               (NCOORD = 1,2,...6 --> x,theta,y,phi,l,delta (Transport vector))
C
C               (Above spectrometer analysis section must be repeated in the
C               input file for the other arm.)
C
C
C               NCUTS:       Number of global cuts
C                  For I=1,NCUTS:
C                       ICUT_VAR(I): Number of variable to cut on
C                       CUT_MIN(I):  Lower limit for cut
C                       CUT_MAX(I):  Upper limit for cut
C
C               NCUTS_I:     Number of cuts to be applied to specific histos
C                  For I=1,NCUTS_I:
C                       ICUT_IND(I):   Index specifying number of cut
C                       ICUT_VAR_I(I): Number of variable to cut on
C                       CUT_MIN_I(I):  Lower limit for cut
C                       CUT_MAX_I(I):  Upper limit for cut
C
C               NPLOTS:      Number of plots
C                  For I=1,NPLOTS:
C                       PLOT_TYPE(I): 3 character ID specifying type of plot
C                            = 'P1D'  --> 1-dimensional histogram.
C                            = 'P2D'  --> 2-dimensional histogram.
C                            = 'SCA'  --> 2-D scatter plot.
C                            = 'NTU'  --> N-Tuple file.
C                       (For PLOT_TYPE='P1D')
C                            IPID(I):     Particle ID for histogram
C                            X_MIN(I):    X Low lim (after scaling)
C                            X_MAX(I):    X Hi  lim (after scaling)
C                            NX_CHAN(I):  Number of X channels
C                            X_SC(I):     X scale factor
C                            X_OFF(I):    X offset (in scaled units)
C                            NCUTS_PL(I): Number of cuts to apply
C                             (For J=1,NCUTS_PL)
C                               ICUT_IND_PL(I,J):  Ind of Jth cut for plot I
C                            I_VAR(I,1):  X variable index
C                            PLOT_FIL(I): Plot filename
C                       (For PLOT_TYPE='P2D')
C                            IPID(I):     Particle ID for histogram
C                            X_MIN(I):    X Low lim (after scaling)
C                            X_MAX(I):    X Hi  lim (after scaling)
C                            Y_MIN(I):    Y Low lim (after scaling)
C                            Y_MAX(I):    Y Hi  lim (after scaling)
C                            NX_CHAN(I):  Number of X channels
C                            NY_CHAN(I):  Number of Y channels
C                            X_SC(I):     X scale factor
C                            Y_SC(I):     Y scale factor
C                            X_OFF(I):    X offset (in scaled units)
C                            Y_OFF(I):    Y offset (in scaled units)
C                            NCUTS_PL(I): Number of cuts to apply
C                             (For J=1,NCUTS_PL)
C                               ICUT_IND_PL(I,J):  Ind of Jth cut for plot I
C                            I_VAR(I,1):  X variable index
C                            I_VAR(I,2):  Y variable index
C                            PLOT_FIL(I): Plot filename
C                       (For PLOT_TYPE='SCA')
C                            X_MIN(I):    X Low lim (after scaling)
C                            X_MAX(I):    X Hi  lim (after scaling)
C                            Y_MIN(I):    Y Low lim (after scaling)
C                            Y_MAX(I):    Y Hi  lim (after scaling)
C                            X_SC(I):     X scale factor
C                            Y_SC(I):     Y scale factor
C                            X_OFF(I):    X offset (in scaled units)
C                            Y_OFF(I):    Y offset (in scaled units)
C                            NCUTS_PL(I): Number of cuts to apply
C                             (For J=1,NCUTS_PL)
C                               ICUT_IND_PL(I,J):  Ind of Jth cut for plot I
C                            I_VAR(I,1):  X variable index
C                            I_VAR(I,2):  Y variable index
C                            PLOT_FIL(I): Plot filename
C                       (For PLOT_TYPE='NTU')
C                            IPID(I):     Particle ID for N-Tuple
C                            NVAR_NTU(I): Number of coordinates in N-Tuple
C                             (For J=1,NVAR_NTU)
C                               I_VAR(I,J):  Variable index
C                            PLOT_FIL(I): N-Tuple filename
C                       (For PLOT_TYPE='NTM')
C                            NWT_NTM(I):  Number of weights in N-Tuple
C                             (For J=1,NWT_NTM)
C                               I_WT(I,J):   Weight index (i.e. PID)
C                            NVAR_NTU(I): Number of coordinates in N-Tuple
C                             (For J=1,NVAR_NTU)
C                               I_VAR(I,J):  Variable index
C                            PLOT_FIL(I): N-Tuple filename
C
C
C       PID (particle ID) refers to reaction type.  The following
C       reactions are handled, though the (e,p) singles routine
C       is currently too slow to be practical. (Note:  negative
C       particle ID tags imply ratios of histograms which result
C       in average cross sections and polarizations.)
C
C            SIGMA(0)   = Phase space (i.e. unit weighting)
C            SIGMA(1)   = (e,e'N)
C            SIGMA(2)   = (e,e'N) for N-type pol. in R frame (total)
C            SIGMA(3)   = (e,e'N) for T-type pol. in R frame (total)
C            SIGMA(4)   = (e,e'N) for L-type pol. in R frame (total)
C            SIGMA(5)   = (e,e'N) for N-type pol. in S frame (total)
C            SIGMA(6)   = (e,e'N) for T-type pol. in S frame (total)
C            SIGMA(7)   = (e,e'N) for L-type pol. in S frame (total)
C            SIGMA(8)   = (e,e'N) for X pol. in scatt. plane (total)
C            SIGMA(9)   = (e,e'N) for Y pol. in scatt. plane (total)
C            SIGMA(10)  = (e,e'N) for Z pol. in scatt. plane (total)
C            SIGMA(11)  = Electron analyzing power
C            SIGMA(12)  = (e,e'N) for N-type pol. in R frame (hel. indep.)
C            SIGMA(13)  = (e,e'N) for T-type pol. in R frame (hel. indep.)
C            SIGMA(14)  = (e,e'N) for L-type pol. in R frame (hel. indep.)
C            SIGMA(15)  = (e,e'N) for N-type pol. in S frame (hel. indep.)
C            SIGMA(16)  = (e,e'N) for T-type pol. in S frame (hel. indep.)
C            SIGMA(17)  = (e,e'N) for L-type pol. in S frame (hel. indep.)
C            SIGMA(18)  = (e,e'N) for X pol. in scatt. plane (hel. indep.)
C            SIGMA(19)  = (e,e'N) for Y pol. in scatt. plane (hel. indep.)
C            SIGMA(20)  = (e,e'N) for Z pol. in scatt. plane (hel. indep.)
C            SIGMA(22)  = (e,e'N) for N-type pol. in R frame (hel. dep.)
C            SIGMA(23)  = (e,e'N) for T-type pol. in R frame (hel. dep.) 
C            SIGMA(24)  = (e,e'N) for L-type pol. in R frame (hel. dep.) 
C            SIGMA(25)  = (e,e'N) for N-type pol. in S frame (hel. dep.) 
C            SIGMA(26)  = (e,e'N) for T-type pol. in S frame (hel. dep.) 
C            SIGMA(27)  = (e,e'N) for L-type pol. in S frame (hel. dep.) 
C            SIGMA(28)  = (e,e'N) for X pol. in scatt. plane (hel. dep.)
C            SIGMA(29)  = (e,e'N) for Y pol. in scatt. plane (hel. dep.)
C            SIGMA(30)  = (e,e'N) for Z pol. in scatt. plane (hel. dep.)
C            SIGMA(50)  = Accidentals for  e',p     or  e'n
C            SIGMA(51)  = Accidentals for  e',pi+   or  e',pi0
C            SIGMA(52)  = Accidentals for  pi-,p    or  pi-,n
C            SIGMA(53)  = Accidentals for  pi-,pi+  or  pi-,pi0
C            SIGMA(100) = (e,e')
C            SIGMA(101) = (e,pi-)
C            SIGMA(102) = (e,p)   for (e,e'p),  or (e,n)   for (e,e'n)
C            SIGMA(103) = (e,pi+) for (e,e'p),  or (e,pi0) for (e,e'n)
C            SIGMA(121-130) = (e,e'p) for multiple weight N-tuples ('NTM')
C
C       Here, the R frame is the Reaction frame used for the actual
C       calculation of the polarization observables.  This frame is
C       event dependent.  The S frame is the Spectrometer frame used
C       to refer the polarizations to a fixed frame in the laboratory,
C       useful for the incorporation of the polarimeter as well as
C       spin transport in the spectrometer.  For the definition of
C       these reference frames see SUBROUTINE ROTATE_POL.
C
C -----------------------------------------------------------------------
C
C     Modified to include the cigar tube target (Model 5).
C     Hassan Ibrahim, December 7, 2004
C

      SUBROUTINE INPUT
      IMPLICIT NONE
C
      DOUBLE PRECISION FWHM_SIG,THETABP,RESMOD,RESVAL,SIG_MSCAT
      DOUBLE PRECISION ACCEPT_ENLARGE,FACT_TMP,OBJECTPT(3),DDUMMY
      INTEGER          NDUMMY,IEL,NELMTS,N_PID,I,J
      LOGICAL          APPLYTOLAB,DUMMYL(10),RADFULL
      CHARACTER*3      DUMMY1
      CHARACTER*1      C1DUMMY(2),C_PID
      CHARACTER*7      FFTYPE
C
      INCLUDE 'input.cmn'
      INCLUDE 'masses.cmn'
      INCLUDE 'spectrometer.cmn'
      INCLUDE 'wc.cmn'
C
      COMMON /FF_MOD/ FFTYPE
      COMMON /RADFULL_L/ RADFULL
C
      DATA FFTYPE /'HALLA3 '/                    !Nucleon form factor model
C
      FWHM_SIG = 2.D0*SQRT(2.D0*LOG(2.D0))       !FWHM --> gaussian sigma
C
      OPEN(UNIT=1,FILE=INFILE,STATUS='OLD',FORM='FORMATTED')
C
C -----------------------------------------------------------------------
C       Read Kinematics and acceptances.
C -----------------------------------------------------------------------
C 
      READ(1,*)N_EVENT
C
      READ(1,*)NITER1,NITER2,NITER3,NITER4,NITER5,NITER6
      IF(NITER1 .EQ. 0) NITER1 = 5      !Default to 5(+1) iters. (from 0 to 5)
      IF(NITER2 .EQ. 0) NITER2 = 5
      IF(NITER3 .EQ. 0) NITER3 = 5
      IF(NITER4 .EQ. 0) NITER4 = 5
      IF(NITER5 .EQ. 0) NITER5 = 5
      IF(NITER6 .EQ. 0) NITER6 = 5
C
      READ(1,*)EJECT_MASS,EJECT_CHG,EM_BOUND
      READ(1,*)E0,PH_B,TH_B,PF_E,PH_E,TH_E,PF_P,PH_P,TH_P
      READ(1,*)ACC_EP_N,ACC_EM_N,ACC_PP_N,ACC_PM_N
      READ(1,*)SA_SHAPE_E,SA_SHAPE_P,DPH_E_N,DTH_E_N,DPH_P_N,DTH_P_N
      IF(SA_SHAPE_E .EQ. 'E' .OR. SA_SHAPE_E .EQ. 'e') THEN
        SA_SHAPE_E = 'E'      !Force capitals
      ELSE
        SA_SHAPE_E = 'R'      !Default to rectangular aperture
      ENDIF
      IF(SA_SHAPE_P .EQ. 'E' .OR. SA_SHAPE_P .EQ. 'e') THEN
        SA_SHAPE_P = 'E'      !Force capitals
      ELSE
        SA_SHAPE_P = 'R'      !Default to rectangular aperture
      ENDIF
C
C ----------------------------------------------------------------------
C       Read counting rate parameters, singles cross section input
C       and target parameters.
C ----------------------------------------------------------------------
C
      READ(1,*)ALUM,BEAMTIME,SPEC_FAC         !rate parameters
      READ(1,*)PFERMI,EPS,EPSD                !for singles X-sections
      READ(1,*)ATARG,ZTARG,DENS_TARG,ITARG_MOD,ELOSS_MOD
                            !target A,Z,density,targ_model,eloss_model
      IF(ITARG_MOD .EQ. 1) THEN
         READ(1,*) PHI_TARG,NFOILS,(TARG_LO(J),TARG_HI(J),J=1,NFOILS)
      ELSEIF(ITARG_MOD .EQ. 2 .OR. ITARG_MOD .EQ. 3 .OR. ITARG_MOD 
     #       .EQ. 4 .OR. ITARG_MOD .EQ. 5) THEN
         PHI_TARG = 0.D0
         NFOILS = 1
         READ(1,*) TARG_LO(1),TARG_HI(1)
      ELSE
         WRITE(6,'(a)') ' Invalid target model selected '
         STOP
      ENDIF
C
      DE_DX_CORR    = .FALSE.   ! Default
      ELOSS_MP_CORR = .FALSE.   ! Default
      IF(ELOSS_MOD .EQ. 0) THEN
         ELOSS_CALC = .FALSE.
      ELSE
         ELOSS_CALC = .TRUE.
         IF(ELOSS_MOD .EQ. 2) DE_DX_CORR    = .TRUE.
         IF(ELOSS_MOD .EQ. 3) ELOSS_MP_CORR = .TRUE.
      ENDIF
C
      READ(1,*)DRIFTE_N,DRIFTP_N              !for extended targets
C
C ----------------------------------------------------------------------
C       If energy loss/multiple scattering calculation is done, enlarge
C       acceptances to account for in-scattering.  In mceep.f, cuts
C       will ensure that the actual acceptances are finally used.
C
C       RADFULL also requires this since the particle momenta
C       sampled are vertex quantities and then changed by the
C       radiated photon.
C ----------------------------------------------------------------------
C
      IF(ELOSS_CALC .OR. RADFULL) THEN
         ACCEPT_ENLARGE = 1.1D0    ! 10% increase in each acceptance
      ELSE
         ACCEPT_ENLARGE = 1.0D0    ! Use nominal values
      ENDIF
      FACT_TMP = (ACCEPT_ENLARGE-1.D0)/2.D0
C
      ACC_EP = ACC_EP_N + FACT_TMP*(ACC_EP_N-ACC_EM_N)
      ACC_EM = ACC_EM_N - FACT_TMP*(ACC_EP_N-ACC_EM_N)
      ACC_PP = ACC_PP_N + FACT_TMP*(ACC_PP_N-ACC_PM_N)
      ACC_PM = ACC_PM_N - FACT_TMP*(ACC_PP_N-ACC_PM_N)
C
      DPH_E  = DPH_E_N * ACCEPT_ENLARGE
      DTH_E  = DTH_E_N * ACCEPT_ENLARGE
      DPH_P  = DPH_P_N * ACCEPT_ENLARGE
      DTH_P  = DTH_P_N * ACCEPT_ENLARGE
C
C ----------------------------------------------------------------------
C       Read beam polarization, dispersion, variations and raster param.
C ----------------------------------------------------------------------
C
      READ(1,*)POL_BEAM,BEAMV,BEAMD,DUTY_FAC,DEL_TOF
      READ(1,*)FWHM_X_B,FWHM_Y_B,FWHM_PH_B,FWHM_TH_B,FWHM_E0
      READ(1,*) DEL_X_B, DEL_Y_B, DEL_PH_B, DEL_TH_B, DEL_E0
      READ(1,*)RAST_SHAPE,RAST_AMP_X,RAST_AMP_Y

      BEAMVAR = .TRUE.
      IF(   FWHM_E0.EQ.0. .AND. FWHM_PH_B.EQ.0. .AND. FWHM_TH_B.EQ.0.
     # .AND. DEL_E0.EQ.0. .AND.  DEL_PH_B.EQ.0. .AND.  DEL_TH_B.EQ.0.)
     #      BEAMVAR = .FALSE.
      FWHM_E0 = E0*(FWHM_E0/100.D0)     !from % to MeV
      DEL_E0  = E0*(DEL_E0 /100.D0)     !from % to MeV

      IF ( RAST_SHAPE.EQ.'R' .OR. RAST_SHAPE.EQ.'r' ) THEN
         RAST_SHAPE='R'
      ELSEIF ( RAST_SHAPE.EQ.'E' .OR. RAST_SHAPE.EQ.'e' ) THEN
         RAST_SHAPE='E'
      ELSE
         RAST_SHAPE='N'         ! This is to enforce capitals
      ENDIF                     ! By default, NO RASTER
C
C ----------------------------------------------------------------------
C       Read spectrometer analysis input.
C ----------------------------------------------------------------------
C
      DO I=1,2
        READ(1,*)C_PID,APPLYTOLAB,NELMTS,THETABP,(OBJECTPT(J),J=1,3)
        IF((C_PID .EQ. 'E') .OR. (C_PID .EQ. 'e')) THEN
          N_PID = 1
        ELSE
          N_PID = 2
        ENDIF
        NEL(N_PID) = NELMTS
        THETA_BP(N_PID) = THETABP
        APPLY_TO_LAB(N_PID) = APPLYTOLAB
        DO J=1,3
           OBJECT_PT(N_PID,J) = OBJECTPT(J)*0.01D0  ! cm --> meters
        ENDDO
        IF(NEL(N_PID) .GT. 0) THEN
          DO IEL=1,NEL(N_PID)
            READ(1,*)DUMMY1  !Read once to determine format
            BACKSPACE 1        !Reread previous line
C
            IF((DUMMY1 .EQ. 'HRS') .OR. (DUMMY1 .EQ. 'hrs')) THEN
              READ(1,*)OP(N_PID,IEL),HRS_ID(N_PID,IEL),
     #                           (DUMMYL(J),J=1,5)
              IF(HRS_ID(N_PID,IEL).EQ.'H'
     #                .OR. HRS_ID(N_PID,IEL).EQ.'h') THEN
                 HRS_ID_INT(N_PID,IEL) = 2
              ELSE
                 HRS_ID_INT(N_PID,IEL) = 1    !Default to electron arm
              ENDIF
              DO J=1,5
                 APERTURE_TEST(HRS_ID_INT(N_PID,IEL),J) = DUMMYL(J)
              ENDDO
              OP(N_PID,IEL) = 'HRS'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
            ELSEIF((DUMMY1 .EQ. 'SEP') .OR. (DUMMY1 .EQ. 'sep')) THEN
              READ(1,*)OP(N_PID,IEL),HRS_ID(N_PID,IEL),
     #                           (DUMMYL(J),J=1,10)
              IF(HRS_ID(N_PID,IEL).EQ.'H'
     #                .OR. HRS_ID(N_PID,IEL).EQ.'h') THEN
                 HRS_ID_INT(N_PID,IEL) = 2
              ELSE
                 HRS_ID_INT(N_PID,IEL) = 1    !Default to electron arm
              ENDIF
              DO J=1,10
                 APER_SEP_TEST(HRS_ID_INT(N_PID,IEL),J) = DUMMYL(J)
              ENDDO
              OP(N_PID,IEL) = 'SEP'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
            ELSEIF((DUMMY1 .EQ. 'MAD') .OR. (DUMMY1 .EQ. 'mad')) THEN
              READ(1,*)OP(N_PID,IEL),NDUMMY
              BACKSPACE 1       !Reread previous line
              IF(NDUMMY .EQ. 1) THEN
                 READ(1,*)OP(N_PID,IEL),MAD_CONFIG,
     #                    (APER_MAD_TEST(J),J=1,4)
              ELSEIF(NDUMMY .EQ. 2 .OR. NDUMMY .EQ. 3) THEN
                 READ(1,*)OP(N_PID,IEL),MAD_CONFIG,
     #                    (APER_MAD_TEST(J),J=1,6)
              ENDIF
              OP(N_PID,IEL) = 'MAD'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C JLab Hall A HRS COSY choice
            ELSEIF((DUMMY1 .EQ. 'SIM') .OR. (DUMMY1 .EQ. 'sim')) THEN
              READ(1,*)OP(N_PID,IEL),HRS_ID(N_PID,IEL),
     #                           (DUMMYL(J),J=1,3)
              IF(HRS_ID(N_PID,IEL).EQ.'H'
     #                .OR. HRS_ID(N_PID,IEL).EQ.'h') THEN
                 HRS_ID_INT(N_PID,IEL) = 2
              ELSE
                 HRS_ID_INT(N_PID,IEL) = 1    !Default to electron arm
              ENDIF
              DO J=1,3
                 HUT_TEST(HRS_ID_INT(N_PID,IEL),J) = DUMMYL(J)
              ENDDO
              OP(N_PID,IEL) = 'SIM'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'HRI') .OR. (DUMMY1 .EQ. 'hri')) THEN
              READ(1,*)OP(N_PID,IEL),HRS_ID(N_PID,IEL)
              IF(HRS_ID(N_PID,IEL).EQ.'H'
     #                .OR. HRS_ID(N_PID,IEL).EQ.'h') THEN
                 HRS_ID_INT(N_PID,IEL) = 2
              ELSE
                 HRS_ID_INT(N_PID,IEL) = 1    !Default to electron arm
              ENDIF
              OP(N_PID,IEL) = 'HRI'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
            ELSEIF((DUMMY1 .EQ. 'SPI') .OR. (DUMMY1 .EQ. 'spi')) THEN
              READ(1,*)OP(N_PID,IEL),HRS_ID(N_PID,IEL)
              IF(HRS_ID(N_PID,IEL).EQ.'H'
     #                .OR. HRS_ID(N_PID,IEL).EQ.'h') THEN
                 HRS_ID_INT(N_PID,IEL) = 2
              ELSE
                 HRS_ID_INT(N_PID,IEL) = 1    !Default to electron arm
              ENDIF
              OP(N_PID,IEL) = 'SPI'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'MDI') .OR. (DUMMY1 .EQ. 'mdi')) THEN
              READ(1,*)OP(N_PID,IEL),MAD_INV_CONFIG
              OP(N_PID,IEL) = 'MDI'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'RFN') .OR. (DUMMY1 .EQ. 'rfn')) THEN
              READ(1,*)OP(N_PID,IEL),HRS_ID(N_PID,IEL),
     #                     C1DUMMY(1),DDUMMY
              IF(HRS_ID(N_PID,IEL).EQ.'H'
     #                .OR. HRS_ID(N_PID,IEL).EQ.'h') THEN
                 HRS_ID_INT(N_PID,IEL) = 2
              ELSE
                 HRS_ID_INT(N_PID,IEL) = 1    !Default to electron arm
              ENDIF
              HRS_COLL(HRS_ID_INT(N_PID,IEL)) = C1DUMMY(1)
              RFN_CUT(HRS_ID_INT(N_PID,IEL))  = DDUMMY
              OP(N_PID,IEL) = 'RFN'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'MAT') .OR. (DUMMY1 .EQ. 'mat')) THEN
              READ(1,*)OP(N_PID,IEL),ORDER(N_PID,IEL),INVERT(N_PID,IEL),
     #             TRNSPT_FIL(N_PID,IEL)
              OP(N_PID,IEL) = 'MAT'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'POL') .OR. (DUMMY1 .EQ. 'pol')) THEN
              READ(1,*)OP(N_PID,IEL),BEND_ANGLE(N_PID,IEL)
              OP(N_PID,IEL) = 'POL'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'COS') .OR. (DUMMY1 .EQ. 'cos')) THEN
              READ(1,*)OP(N_PID,IEL),COSY_SPIN_FIL(N_PID,IEL)
              OP(N_PID,IEL) = 'COS'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'DFT') .OR. (DUMMY1 .EQ. 'dft')) THEN
              READ(1,*)OP(N_PID,IEL),DRIFT_LENGTH(N_PID,IEL)
              OP(N_PID,IEL) = 'DFT'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'ROT') .OR. (DUMMY1 .EQ. 'rot')) THEN
              READ(1,*)OP(N_PID,IEL),ROT_AXIS(N_PID,IEL),
     #           ROT_ANGLE(N_PID,IEL)
              IF (ROT_AXIS(N_PID,IEL) .EQ. 'x' .OR.
     #               ROT_AXIS(N_PID,IEL) .EQ. 'X') THEN
                         ROT_AXIS(N_PID,IEL) = 'X' 
              ELSEIF (ROT_AXIS(N_PID,IEL) .EQ. 'y' .OR.
     #               ROT_AXIS(N_PID,IEL) .EQ. 'Y') THEN
                         ROT_AXIS(N_PID,IEL) = 'Y' 
              ELSEIF (ROT_AXIS(N_PID,IEL) .EQ. 'z' .OR.
     #              ROT_AXIS(N_PID,IEL) .EQ. 'Z') THEN
                         ROT_AXIS(N_PID,IEL) = 'Z' 
              ELSE
                  WRITE(6,211)
211               FORMAT(/' Illegal rotation axis entry ')
                  STOP
              ENDIF
              OP(N_PID,IEL) = 'ROT'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'TRK') .OR. (DUMMY1 .EQ. 'trk')) THEN
              READ(1,*)OP(N_PID,IEL),WC_FIL(N_PID)
              OP(N_PID,IEL) = 'TRK'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'TOF') .OR. (DUMMY1 .EQ. 'tof')) THEN
              READ(1,*)OP(N_PID,IEL),CRAY_PATH(N_PID)
              OP(N_PID,IEL) = 'TOF'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'RES') .OR. (DUMMY1 .EQ. 'res')) THEN
              READ(1,*)OP(N_PID,IEL),OB(N_PID,IEL),
     #               NCOORD(N_PID,IEL,1),RESMOD,RESVAL
              OP(N_PID,IEL) = 'RES'           !Make all CAPS
              IF(RESMOD .EQ. 1) THEN
                 GSIG(N_PID,IEL) = RESVAL/FWHM_SIG     !FWHM --> sigma
              ELSE
                 GSIG(N_PID,IEL) =
     #                   SIG_MSCAT(N_PID,PF_E,PF_P,RESVAL,1.D0)
              ENDIF
              IF(OB(N_PID,IEL).EQ.'D'.OR.OB(N_PID,IEL).EQ.'d')
     #               OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'OFF') .OR. (DUMMY1 .EQ. 'off')) THEN
              READ(1,*)OP(N_PID,IEL),OB(N_PID,IEL),
     #               NCOORD(N_PID,IEL,1),OFF(N_PID,IEL)
              OP(N_PID,IEL) = 'OFF'           !Make all CAPS
              IF(OB(N_PID,IEL).EQ.'D'.OR.OB(N_PID,IEL).EQ.'d')
     #               OB(N_PID,IEL) = 'S'
C
            ELSEIF((DUMMY1 .EQ. 'H1D') .OR. (DUMMY1 .EQ. 'h1d')) THEN
              READ(1,*)OP(N_PID,IEL),IPID_TR(N_PID,IEL),
     #              TR_MINX(N_PID,IEL),TR_MAXX(N_PID,IEL),
     #              NXCHAN_TR(N_PID,IEL),NDUMMY
              BACKSPACE 1        !Reread previous line
              IF(NDUMMY .GT. 0) THEN
                READ(1,*)OP(N_PID,IEL),IPID_TR(N_PID,IEL),
     #              TR_MINX(N_PID,IEL),TR_MAXX(N_PID,IEL),
     #              NXCHAN_TR(N_PID,IEL),NCUTS_TR(N_PID,IEL),
     #              (ICUT_IND_TR(N_PID,IEL,J),J=1,NCUTS_TR(N_PID,IEL)),
     #              NCOORD(N_PID,IEL,1),TR_FILE(N_PID,IEL)
              ELSE
                READ(1,*)OP(N_PID,IEL),IPID_TR(N_PID,IEL),
     #              TR_MINX(N_PID,IEL),TR_MAXX(N_PID,IEL),
     #              NXCHAN_TR(N_PID,IEL),NCUTS_TR(N_PID,IEL),
     #              NCOORD(N_PID,IEL,1),TR_FILE(N_PID,IEL)
              ENDIF
              OP(N_PID,IEL) = 'H1D'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'H2D') .OR. (DUMMY1 .EQ. 'h2d')) THEN
              READ(1,*)OP(N_PID,IEL),IPID_TR(N_PID,IEL),
     #              TR_MINX(N_PID,IEL),TR_MAXX(N_PID,IEL),
     #              TR_MINY(N_PID,IEL),TR_MAXY(N_PID,IEL),
     #              NXCHAN_TR(N_PID,IEL),NYCHAN_TR(N_PID,IEL),
     #              NDUMMY
              BACKSPACE 1        !Reread previous line
              IF(NDUMMY .GT. 0) THEN
                READ(1,*)OP(N_PID,IEL),IPID_TR(N_PID,IEL),
     #              TR_MINX(N_PID,IEL),TR_MAXX(N_PID,IEL),
     #              TR_MINY(N_PID,IEL),TR_MAXY(N_PID,IEL),
     #              NXCHAN_TR(N_PID,IEL),NYCHAN_TR(N_PID,IEL),
     #              NCUTS_TR(N_PID,IEL),
     #              (ICUT_IND_TR(N_PID,IEL,J),J=1,NCUTS_TR(N_PID,IEL)),
     #              NCOORD(N_PID,IEL,1),NCOORD(N_PID,IEL,2),
     #              TR_FILE(N_PID,IEL)
              ELSE
                READ(1,*)OP(N_PID,IEL),IPID_TR(N_PID,IEL),
     #              TR_MINX(N_PID,IEL),TR_MAXX(N_PID,IEL),
     #              TR_MINY(N_PID,IEL),TR_MAXY(N_PID,IEL),
     #              NXCHAN_TR(N_PID,IEL),NYCHAN_TR(N_PID,IEL),
     #              NCUTS_TR(N_PID,IEL),
     #              NCOORD(N_PID,IEL,1),NCOORD(N_PID,IEL,2),
     #              TR_FILE(N_PID,IEL)
              ENDIF
              OP(N_PID,IEL) = 'H2D'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'SCT') .OR. (DUMMY1 .EQ. 'sct')) THEN
              READ(1,*)OP(N_PID,IEL),
     #              TR_MINX(N_PID,IEL),TR_MAXX(N_PID,IEL),
     #              TR_MINY(N_PID,IEL),TR_MAXY(N_PID,IEL),
     #              NDUMMY
              BACKSPACE 1        !Reread previous line
              IF(NDUMMY .GT. 0) THEN
                READ(1,*)OP(N_PID,IEL),
     #              TR_MINX(N_PID,IEL),TR_MAXX(N_PID,IEL),
     #              TR_MINY(N_PID,IEL),TR_MAXY(N_PID,IEL),
     #              NCUTS_TR(N_PID,IEL),
     #              (ICUT_IND_TR(N_PID,IEL,J),J=1,NCUTS_TR(N_PID,IEL)),
     #              NCOORD(N_PID,IEL,1),NCOORD(N_PID,IEL,2),
     #              TR_FILE(N_PID,IEL)
              ELSE
                READ(1,*)OP(N_PID,IEL),
     #              TR_MINX(N_PID,IEL),TR_MAXX(N_PID,IEL),
     #              TR_MINY(N_PID,IEL),TR_MAXY(N_PID,IEL),
     #              NCUTS_TR(N_PID,IEL),
     #              NCOORD(N_PID,IEL,1),NCOORD(N_PID,IEL,2),
     #              TR_FILE(N_PID,IEL)
              ENDIF
              OP(N_PID,IEL) = 'SCT'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'NTU') .OR. (DUMMY1 .EQ. 'ntu')) THEN
              READ(1,*)OP(N_PID,IEL),IPID_TR(N_PID,IEL),
     #              BEND_ANGLE_TOT(N_PID,IEL),TR_FILE(N_PID,IEL)
              OP(N_PID,IEL) = 'NTU'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
C
            ELSEIF((DUMMY1 .EQ. 'CUT') .OR. (DUMMY1 .EQ. 'cut')) THEN
              READ(1,*)OP(N_PID,IEL),(C1DUMMY(J),J=1,2)
              BACKSPACE 1        !Reread previous line
              IF( (C1DUMMY(1).EQ.'G' .OR. C1DUMMY(1).EQ.'g') .AND.
     #            (C1DUMMY(2).EQ.'R' .OR. C1DUMMY(2).EQ.'r') ) THEN
                 READ(1,*)OP(N_PID,IEL),CUT_DOMAIN(N_PID,IEL),
     #               CUT_TYPE(N_PID,IEL),TR_CUT_MIN(N_PID,IEL),
     #               TR_CUT_MAX(N_PID,IEL),NCOORD(N_PID,IEL,1)
                     CUT_DOMAIN(N_PID,IEL) = 'G'
                     CUT_TYPE(N_PID,IEL)   = 'R'
              ELSEIF( (C1DUMMY(1).EQ.'G' .OR. C1DUMMY(1).EQ.'g') .AND.
     #                (C1DUMMY(2).EQ.'E' .OR. C1DUMMY(2).EQ.'e') ) THEN
                 READ(1,*)OP(N_PID,IEL),CUT_DOMAIN(N_PID,IEL),
     #               CUT_TYPE(N_PID,IEL),TR_CUT_X(N_PID,IEL),
     #               TR_CUT_Y(N_PID,IEL),NCOORD(N_PID,IEL,1),
     #               NCOORD(N_PID,IEL,2)
                     CUT_DOMAIN(N_PID,IEL) = 'G'
                     CUT_TYPE(N_PID,IEL)   = 'E'
              ELSEIF( (C1DUMMY(1).EQ.'S' .OR. C1DUMMY(1).EQ.'s') .AND.
     #                (C1DUMMY(2).EQ.'R' .OR. C1DUMMY(2).EQ.'r') ) THEN
                 READ(1,*)OP(N_PID,IEL),CUT_DOMAIN(N_PID,IEL),
     #               CUT_TYPE(N_PID,IEL),TR_CUT_MIN(N_PID,IEL),
     #               TR_CUT_MAX(N_PID,IEL),NCOORD(N_PID,IEL,1),
     #               TR_CUT_IND(N_PID,IEL)
                     CUT_DOMAIN(N_PID,IEL) = 'S'
                     CUT_TYPE(N_PID,IEL)   = 'R'
              ELSEIF( (C1DUMMY(1).EQ.'S' .OR. C1DUMMY(1).EQ.'s') .AND.
     #                (C1DUMMY(2).EQ.'E' .OR. C1DUMMY(2).EQ.'e') ) THEN
                 READ(1,*)OP(N_PID,IEL),CUT_DOMAIN(N_PID,IEL),
     #               CUT_TYPE(N_PID,IEL),TR_CUT_X(N_PID,IEL),
     #               TR_CUT_Y(N_PID,IEL),NCOORD(N_PID,IEL,1),
     #               NCOORD(N_PID,IEL,2),TR_CUT_IND(N_PID,IEL)
                     CUT_DOMAIN(N_PID,IEL) = 'S'
                     CUT_TYPE(N_PID,IEL)   = 'E'
              ENDIF
              OP(N_PID,IEL) = 'CUT'           !Make all CAPS
              OB(N_PID,IEL) = 'K'
            ENDIF
C
            IF(OB(N_PID,IEL).EQ.'K'.OR.OB(N_PID,IEL).EQ.'k')THEN
              OB_NAM(N_PID,IEL) = 'Kin. bin'
            ELSEIF(OB(N_PID,IEL).EQ.'S'.OR.OB(N_PID,IEL).EQ.'s')THEN
              OB_NAM(N_PID,IEL) = 'X-sect. '
            ELSE
              OB_NAM(N_PID,IEL) = 'Both    '
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IF(NEL(1) .EQ. 0 .AND. NEL(2) .EQ. 0) THEN
        CALLSPEC = .FALSE.       !need not do spectrometer analysis
      ELSE
        CALLSPEC = .TRUE.
      ENDIF
C
C ----------------------------------------------------------------------
C       Read global cuts.
C ----------------------------------------------------------------------
C
      READ(1,*)NCUTS
      IF(NCUTS .GT. 0) THEN
        DO I=1,NCUTS
          READ(1,*)ICUT_VAR(I),CUT_MIN(I),CUT_MAX(I)
        ENDDO
      ENDIF
C
C ----------------------------------------------------------------------
C       Read specific cuts.
C ----------------------------------------------------------------------
C
      READ(1,*)NCUTS_I
      IF(NCUTS_I .GT. 0) THEN
        DO I=1,NCUTS_I
          READ(1,*)ICUT_IND(I),ICUT_VAR_I(ICUT_IND(I)),
     #            CUT_MIN_I(ICUT_IND(I)),CUT_MAX_I(ICUT_IND(I))
        ENDDO
      ENDIF
C
C ----------------------------------------------------------------------
C       Read input for various PLOTS.
C       NOTE:
C            For negative PID (i.e. ratio histograms) the PID for
C            spectrum which divides yield histos must be found.  The
C            ID for the "denominator" spectrum is found by the routine
C            GET_PID_DEN.  The ID returned is either
C            0 (for cross sections) or 1 (for polarizations or 
C            electron analyzing power).
C ----------------------------------------------------------------------
C
      READ(1,*)NPLOTS
      IF(NPLOTS .GT. 0) THEN
        DO I=1,NPLOTS
          TWO_D(I) = .FALSE. !Initialize flag
          READ(1,*)DUMMY1  !Read once to determine format
          BACKSPACE 1        !Reread previous line
C                                                1-D HISTOS
          IF((DUMMY1 .EQ. 'P1D') .OR. (DUMMY1 .EQ. 'p1d')) THEN
            READ(1,*) PLOT_TYPE(I),IPID(I),X_MIN(I),X_MAX(I),
     #                NX_CHAN(I),X_SC(I),X_OFF(I),NDUMMY
            BACKSPACE 1        !Reread previous line
            IF(NDUMMY .GT. 0) THEN
              READ(1,*) PLOT_TYPE(I),IPID(I),X_MIN(I),X_MAX(I),
     #                  NX_CHAN(I),X_SC(I),X_OFF(I),NCUTS_PL(I),
     #                  (ICUT_IND_PL(I,J),J=1,NCUTS_PL(I)),
     #                  I_VAR(I,1),PLOT_FIL(I)
            ELSE
              READ(1,*) PLOT_TYPE(I),IPID(I),X_MIN(I),X_MAX(I),
     #                  NX_CHAN(I),X_SC(I),X_OFF(I),NCUTS_PL(I),
     #                  I_VAR(I,1),PLOT_FIL(I)
            ENDIF
            CALL GET_PID_DEN(IPID(I),IPID_DEN(I))
            PLOT_TYPE(I) = 'P1D'          !Make all CAPS
C                                                2-D HISTOS
          ELSEIF((DUMMY1 .EQ. 'P2D') .OR. (DUMMY1 .EQ. 'p2d')) THEN
            READ(1,*) PLOT_TYPE(I),IPID(I),X_MIN(I),X_MAX(I),
     #                Y_MIN(I),Y_MAX(I),NX_CHAN(I),NY_CHAN(I),
     #                X_SC(I),Y_SC(I),X_OFF(I),Y_OFF(I),NDUMMY
            BACKSPACE 1        !Reread previous line
            IF(NDUMMY .GT. 0) THEN
              READ(1,*) PLOT_TYPE(I),IPID(I),X_MIN(I),X_MAX(I),
     #                  Y_MIN(I),Y_MAX(I),NX_CHAN(I),NY_CHAN(I),
     #                  X_SC(I),Y_SC(I),X_OFF(I),Y_OFF(I),NCUTS_PL(I),
     #                  (ICUT_IND_PL(I,J),J=1,NCUTS_PL(I)),
     #                  I_VAR(I,1),I_VAR(I,2),PLOT_FIL(I)
            ELSE
              READ(1,*) PLOT_TYPE(I),IPID(I),X_MIN(I),X_MAX(I),
     #                  Y_MIN(I),Y_MAX(I),NX_CHAN(I),NY_CHAN(I),
     #                  X_SC(I),Y_SC(I),X_OFF(I),Y_OFF(I),NCUTS_PL(I),
     #                  I_VAR(I,1),I_VAR(I,2),PLOT_FIL(I)
            ENDIF
            CALL GET_PID_DEN(IPID(I),IPID_DEN(I))
            TWO_D(I) = .TRUE.
            PLOT_TYPE(I) = 'P2D'          !Make all CAPS
C                                                2-D SCATTER PLOTS
          ELSEIF((DUMMY1 .EQ. 'SCA') .OR. (DUMMY1 .EQ. 'sca')) THEN
            READ(1,*) PLOT_TYPE(I),X_MIN(I),X_MAX(I),
     #                Y_MIN(I),Y_MAX(I),
     #                X_SC(I),Y_SC(I),X_OFF(I),Y_OFF(I),NDUMMY
            BACKSPACE 1        !Reread previous line
            IF(NDUMMY .GT. 0) THEN
              READ(1,*) PLOT_TYPE(I),X_MIN(I),X_MAX(I),
     #                  Y_MIN(I),Y_MAX(I),
     #                  X_SC(I),Y_SC(I),X_OFF(I),Y_OFF(I),NCUTS_PL(I),
     #                  (ICUT_IND_PL(I,J),J=1,NCUTS_PL(I)),
     #                  I_VAR(I,1),I_VAR(I,2),PLOT_FIL(I)
            ELSE
              READ(1,*) PLOT_TYPE(I),X_MIN(I),X_MAX(I),
     #                  Y_MIN(I),Y_MAX(I),
     #                  X_SC(I),Y_SC(I),X_OFF(I),Y_OFF(I),NCUTS_PL(I),
     #                  I_VAR(I,1),I_VAR(I,2),PLOT_FIL(I)
            ENDIF
            TWO_D(I) = .TRUE.
            PLOT_TYPE(I) = 'SCA'          !Make all CAPS
C                                                N-TUPLE FILES
          ELSEIF((DUMMY1 .EQ. 'NTU') .OR. (DUMMY1 .EQ. 'ntu')) THEN
            READ(1,*) PLOT_TYPE(I),IPID(I),NVAR_NTU(I),
     #                (I_VAR(I,J),J=1,NVAR_NTU(I)),
     #                PLOT_FIL(I)
            CALL GET_PID_DEN(IPID(I),IPID_DEN(I))
            PLOT_TYPE(I) = 'NTU'          !Make all CAPS
C
          ELSEIF((DUMMY1 .EQ. 'NTM') .OR. (DUMMY1 .EQ. 'ntm ')) THEN
            READ(1,*) PLOT_TYPE(I),NWT_NTM(I),
     #                (I_WT(I,J),J=1,NWT_NTM(I)),
     #                NVAR_NTU(I),(I_VAR(I,J),J=1,NVAR_NTU(I)),
     #                PLOT_FIL(I)
            PLOT_TYPE(I) = 'NTM'          !Make all CAPS
C
          ENDIF
        ENDDO
      ENDIF
C
C ----------------------------------------------------------------------
C       Read comment lines (up to 10 lines will appear in the summary 
C       file).
C ----------------------------------------------------------------------
C
      NCOMMENTS = 0     !Initialize number of comment lines.
      DO I=1,10
        READ(1,304,END=999) COMMENT(I)
304     FORMAT(A80)
        NCOMMENTS = I
      ENDDO
C
C
999   CLOSE(UNIT=1)
C
      RETURN
      END
