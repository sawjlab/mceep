C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE PHYS_CHOICE
C
C       AUTHOR:  P.E. ULMER
C       DATE:    AUG-30-1990
C
C       PURPOSE: Choose physics option from menu.
C
C       References for various options:
C
C         (100) Deforest - CC1:  Taber de Forest, Jr.,
C                                Nucl. Phys. A392, 232 (1983).
C         (101) Deforest - CC2:  Taber de Forest, Jr.,
C                                Nucl. Phys. A392, 232 (1983).
C         (200) Van Orden - PW:  A. Picklesimer and J.W. Van Orden,
C                                Phys. Rev. C 40, 290 (1989).
C         (300) Deuterium - IA:  V. Dmitrasinovic and F. Gross,
C                                Phys. Rev. C 40, 2479 (1989).
C         (350) Deuterium:       S. Jeschonnek
C         (500) Pion production: R.W. Lourie,
C                                Nucl. Phys. A509, 653 (1990).
C         (500) Pion production: EPIPROD , J. J. Kelly.
C               (Currently, Lourie option used.  Kelly option not
C                yet available in standard release of MCEEP.)
C         (600) p(e,e'pi+)n      Brauel et al.,
C                                Z. Phys. C. 3, 101 (1979).
C         (700) p(e,e'K+)Y       G.Niculescu et.al., 
C                                Phys. Rev. Lett. 81, 4576 (1998) and
C                                J.Cha (PhD theses, HALL C web page, 2001),
C                                in turn similar to P.Brauel et.al., 
C                                Z.Phys. C, 3 p. 101 (1979).
C         (801) A(e,e'p) from J.J.Kelly (LEA) response functions
C         (802) A(e,e'p) from J.M.Laget response functions
C         (803) A(e,e'p) from J.M.Udias response functions
C         (804) d(e,e'p)n Arenhoevel unpolarized responses
C         (814) d(e,e'p)n Laget      unpolarized responses
C
C       Includes options for elastic scattering off:
C                proton
C                deuteron
C                triton
C                3He
C                4He
C                12C
C
C------------------------------------------------------------------------------
C
      SUBROUTINE PHYS_CHOICE(ELASTIC,ACCEPT_CHECK,BOUND,PROTON,
     #                       RADPEAKING,MULTIPHOTON,RADFULL,ACCEPT_FCN)
C
      COMMON /AMPL/    ZRFF1,ZRFF2,ZRFF3,ZRFF4,ZRFF5,ZRFF6,
     2                 ZIFF1,ZIFF2,ZIFF3,ZIFF4,ZIFF5,ZIFF6
      COMMON /DEUT_FF/ FC,FM,FQ,Q_2
      COMMON /N_DEUT_FF/ N_D_PTS
      COMMON /PHYSVAR/ IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
      COMMON /SURVIVE_L/ DECAY_CHECK
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
      COMMON /RADPKAPPR/ cutoff,eg_max

      COMMON / GRID_DEU / Q2_STEP,OM_STEP,TA_STEP
      COMMON / GRID_PAR / NQ2,NOM,NTA
      COMMON / GRID_SIG / SIG_L,SIG_T,SIG_LT,SIG_TT
      COMMON / LAGET_MOD/ LAGET_INTP,LAGET_PWIA,LAGET_FSI,LAGET_MEC
C
      LOGICAL          ELASTIC,ACCEPT_CHECK,BOUND,PROTON,RADPEAKING
      LOGICAL          MULTIPHOTON,ANSWER,DECAY_CHECK,RADFULL
      LOGICAL          ACCEPT_FCN
      CHARACTER*100    DAT_DIR
      CHARACTER*200    FILENAME
      CHARACTER*200    RES_OPT_FILE
      CHARACTER*200    FC_FILE,FM_FILE,FQ_FILE
      CHARACTER*1      ANS
      DOUBLE PRECISION FC(100),FQ(100),FM(100),Q_2(100)
      DOUBLE PRECISION DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7
      DOUBLE PRECISION CUTOFF,EG_MAX
      DOUBLE PRECISION SIG_L(100,500,91),SIG_T(100,500,91)
      DOUBLE PRECISION SIG_LT(100,500,91), SIG_TT(100,500,91)
      DOUBLE PRECISION Q2_STEP,OM_STEP,TA_STEP
      REAL             ZRFF1(10,11,19),ZRFF2(10,11,19),ZRFF3(10,11,19),
     3                 ZRFF4(10,11,19),ZRFF5(10,11,19),ZRFF6(10,11,19),
     4                 ZIFF1(10,11,19),ZIFF2(10,11,19),ZIFF3(10,11,19),
     5                 ZIFF4(10,11,19),ZIFF5(10,11,19),ZIFF6(10,11,19)
      INTEGER          IREAD,IQ2,IANG
      INTEGER          N_D_PTS
      INTEGER          NQ2,NOM,NTA
      INTEGER          LAGET_INTP,LAGET_PWIA,LAGET_FSI,LAGET_MEC
C
C----------------------------------------------------------------------
C     Choose between elastic scattering A(e,e'A), (e,e'X) to bound
C     state or (e,e'X) to continuum.  As of Version 3.9, can also
C     choose "Acceptance".
C
C     If elastic scattering is chosen, inquire whether hadron arm
C     acceptances are to be enforced.  If not, then the results
C     of MCEEP correspond to single-arm elastic scattering, A(e,e');
C     if so then results correspond to A(e,e'A) where detection of
C     the recoil nucleus is required.
C----------------------------------------------------------------------
C
      RADPEAKING   = .FALSE.       !initialize
      MULTIPHOTON  = .FALSE.       !initialize
      RADFULL      = .FALSE.       !initialize
      ELASTIC      = .FALSE.       !initialize
      BOUND        = .FALSE.       !initialize
      PROTON       = .FALSE.       !initialize
      ACCEPT_CHECK = .TRUE.        !default to A(e,e'A) not A(e,e')
      DECAY_CHECK  = .TRUE.        !default pion decay in p(e,e'pi+)n
                                    !or kaon decay in p(e,e'K)Lambda
      ACCEPT_FCN   = .FALSE.       !initialize
C
      WRITE(6,100)
  100 FORMAT(' Elast (E), Bnd State (B), Contin (C), Acceptance (A)? ')
      READ(5,101) ANS
  101 FORMAT(A1)
      IF((ANS .EQ. 'E') .OR. (ANS .EQ. 'e')) ELASTIC = .TRUE.
      IF((ANS .EQ. 'B') .OR. (ANS .EQ. 'b')) BOUND   = .TRUE.
C
      IF((ANS .EQ. 'A') .OR. (ANS .EQ. 'a')) THEN
         ACCEPT_FCN = .TRUE.
         WRITE(6,105)
 105     FORMAT(' You have chosen to compute the acceptance profiles ',
     #       /,'  No kinematics or cross sections will be computed ',
     #       /,'  Energy loss/multiple scatt will be turned OFF ')
         GOTO 999
      ENDIF
C
      IF(.NOT. ELASTIC) THEN
         WRITE(6,102)
  102    FORMAT(' Radiate? (0=No, 1=Yes, 2=Yes & multiphoton corr.) ')
         READ(5,*) IRAD
      ELSE
         WRITE(6,107)
  107    FORMAT(' Radiate? ',/,
     #          ' 0=No  ',/,
     #          ' 1=Yes ',/,
     #          ' 2=Yes & multiphoton correction ',/,
     #          ' 3=full angular dist. (ep elastic ONLY) ')
         READ(5,*) IRAD
      ENDIF
C
      IF(IRAD.EQ.1 .OR. IRAD.EQ.2) RADPEAKING  = .TRUE.
      IF(IRAD.EQ.2) MULTIPHOTON = .TRUE.
      IF(IRAD.EQ.3) RADFULL     = .TRUE.
C
      IF(RADPEAKING .OR. RADFULL) THEN
         WRITE(6,103)
  103    FORMAT(' NOTE:  Cross sections will be additionally ',/,
     #          '        differential in the photon energy. ',/,
     #          ' Cutoff energy in MeV (RET for 1 MeV default) > ')
         READ(5,104) CUTOFF
  104    FORMAT(D9.3)
         IF(CUTOFF .EQ. 0.D0) CUTOFF = 1.D0
      ENDIF
C
      IF(ELASTIC) THEN
C
    4   WRITE(6,1001)
 1001   FORMAT(
     #  /' Proton Elastic --------------- 10'
     #  /' Deuteron Elastic ------------- 20'
     #  /' Tritium Elastic -------------- 30'
     #  /' 3He Elastic ------------------ 40'
     #  /' 4He Elastic ------------------ 50'
     #  /' 12C Elastic ------------------ 60'
     #  //' Enter Option>')
C
        READ(5,1011,ERR=4)IELAST_OPT
 1011   FORMAT(I2)
C
        IF(IELAST_OPT .EQ. 20)THEN
C
           DO J=1,100
              FC(J)  = 0.
              FM(J)  = 0.
              FQ(J)  = 0.
              Q_2(J) = 0.
           ENDDO
C
C ------------------------------------------------------------------------
C       Ask for deuteron form factor model (from Tjon), then read them
C       in from files in DAT$D
C ------------------------------------------------------------------------
C
           WRITE(6,'(a)')' Model: ia,iamec,rsc,rscmec (1,2,3,4) >'
           READ(5,*)IMODEL
           IF(IMODEL .EQ. 1)THEN
              FC_FILE =  DAT_DIR(1:N_DAT_DIR)//'/ia.fc'
              FQ_FILE =  DAT_DIR(1:N_DAT_DIR)//'/ia.fq'
              FM_FILE =  DAT_DIR(1:N_DAT_DIR)//'/ia.fm'
           ELSEIF(IMODEL .EQ. 2) THEN
              FC_FILE =  DAT_DIR(1:N_DAT_DIR)//'/iamec.fc'
              FQ_FILE =  DAT_DIR(1:N_DAT_DIR)//'/iamec.fq'
              FM_FILE =  DAT_DIR(1:N_DAT_DIR)//'/iamec.fm'
           ELSEIF(IMODEL .EQ. 3) THEN
              FC_FILE =  DAT_DIR(1:N_DAT_DIR)//'/rsc.fc'
              FQ_FILE =  DAT_DIR(1:N_DAT_DIR)//'/rsc.fq'
              FM_FILE =  DAT_DIR(1:N_DAT_DIR)//'/rsc.fm'
           ELSE
              FC_FILE =  DAT_DIR(1:N_DAT_DIR)//'/rscmec.fc'
              FQ_FILE =  DAT_DIR(1:N_DAT_DIR)//'/rscmec.fq'
              FM_FILE =  DAT_DIR(1:N_DAT_DIR)//'/rscmec.fm'
           ENDIF
C
           OPEN(UNIT=11,FILE=FM_FILE,STATUS='OLD')
           OPEN(UNIT=12,FILE=FC_FILE,STATUS='OLD')
           OPEN(UNIT=13,FILE=FQ_FILE,STATUS='OLD')
C
           DO J=1,100
              READ(11,*,END=9) Q_2(J),FM(J)
              READ(12,*) Q_2(J),FC(J)
              READ(13,*) Q_2(J),FQ(J)
              N_D_PTS = J
           ENDDO
C
9          DO ILUN = 11,13
              CLOSE(UNIT=ILUN)
           ENDDO
C
        ENDIF
C
        ACCEPT_CHECK =
     #      ANSWER(' Enforce hadron-arm accept. tests? (Y/N)>','Y','y')
        GOTO 999    !no need for further selections
C
      ENDIF
C
    5 WRITE(6,1000)
 1000 FORMAT(
     #  /' S(p,E) X Sigma_CC1 (unpolarized) --------------- 100'
     #  /' S(p,E) X Sigma_CC2 (unpolarized) --------------- 101'
     #  /' S(p,E) X Sigma_CC1 (polarized - Van Orden) ----- 200'
     #  /' Gross/Dmitrasinovic Deuterium IA --------------- 300'
     #  /' S. Jeschonnek Deuterium PWBA ------------------- 350'
cxx  #  /' Enter response functions from files (not ready)  400'
     #  /' p(e,e''p)pi0 ------------------------------------ 500'
     #  /' p(e,e''pi+)n ------------------------------------ 600'
     #  /' p(e,e''K+)Y ------------------------------------- 700'
c     #  /' A(e,e''p) J.J. KELLY (LEA) response functions --- 801'
c     #  /' A(e,e''p) J.M. LAGET response functions --------- 802'
c     #  /' A(e,e''p) J.M. Udias response functions --------- 803'
     #  /' d(e,e''p)n Arenhoevel unpolarized responses ----- 804'
     #  /' d(e,e''p)n Laget      unpolarized responses ----- 814'
     # //' Enter Option>')
C
      READ(5,1010,ERR=5)IPHYS_OPT
 1010 FORMAT(I3)
C
      IF(IPHYS_OPT .EQ. 100 .OR. IPHYS_OPT .EQ. 101
     #      .OR. IPHYS_OPT .EQ. 200) THEN
         PROTON =ANSWER(' (e,e''p) or (e,e''n)? (P/N)>','P','p')
      ELSEIF(IPHYS_OPT .EQ. 300 .OR. IPHYS_OPT .EQ. 350
     #      .OR. IPHYS_OPT .EQ. 500 .OR. IPHYS_OPT .EQ. 804
     #      .OR. IPHYS_OPT .EQ. 814) THEN
         PROTON = .TRUE.
      ELSEIF(IPHYS_OPT .EQ. 600 .OR. IPHYS_OPT .EQ. 700) THEN
         PROTON = .FALSE.
      ENDIF

C ------------------------------------------------------------------------
C       Sigma_CC1 (unpolarized) X S(p,E)
C ------------------------------------------------------------------------
C
      IF(IPHYS_OPT .EQ. 100) THEN
        CALL SPECT_SETUP(BOUND)
C
C ------------------------------------------------------------------------
C       Sigma_CC2 (unpolarized) X S(p,E)
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 101) THEN
        CALL SPECT_SETUP(BOUND)
C
C ------------------------------------------------------------------------
C       Sigma_CC1 (polarized - Van Orden) X S(p,E)
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 200) THEN
        CALL SPECT_SETUP(BOUND)
C
C ------------------------------------------------------------------------
C       Gross/Dmitrasinovic Deuterium IA
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 300) THEN
        IF(.NOT. BOUND) WRITE(6,3000)
        FILENAME = DAT_DIR(1:N_DAT_DIR)//'/wave_fourx.dat'
        OPEN(UNIT=15,FILE=FILENAME,STATUS='OLD')
        CALL WAVEIN(15)
        CALL WAVSET
        CLOSE(UNIT=15)
C
C ------------------------------------------------------------------------
C       Sabine Jeschonnek Deuterium PWBA
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 350) THEN
        IF(.NOT. BOUND) WRITE(6,3000)
        WRITE(6,2500) 
 2500   FORMAT(' Reading parameters from file deut_sabjes.dat ')
        FILENAME = DAT_DIR(1:N_DAT_DIR)//'/deut_sabjes.dat'
        OPEN(UNIT=15,FILE=FILENAME,STATUS='OLD')
        CALL GET_SABJES_PAR
        CLOSE(UNIT=15)
C
C ------------------------------------------------------------------------
C       Enter response functions from files (not implemented yet).
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 400) THEN
        CALL GET_RESPONSE
C
C ------------------------------------------------------------------------
C       p(e,e'p)pi0
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 500) THEN
CXXX        CALL INIT_MCEEP_MODE   ! Kelly option
        IF(.NOT. BOUND) WRITE(6,3000)
        WRITE(6,'(A)')'$ Resonance option file> '
        READ(5,'(A)')RES_OPT_FILE
        RES_OPT_FILE = DAT_DIR(1:N_DAT_DIR)//'/'//RES_OPT_FILE
        OPEN(UNIT=15,FILE=RES_OPT_FILE,STATUS='OLD')
        CALL SETUP_PION(15)
        CLOSE(UNIT=15)
C
C ------------------------------------------------------------------------
C       p(e,e'pi+)n
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 600) THEN
        IF(.NOT. BOUND) WRITE(6,3000)
        DECAY_CHECK =
     #       ANSWER(' Pion Decay Global Cut On? (Y/N)>','Y','y')
C
 3000   FORMAT(' Warning:  This option is explicitly bound state ')
C
C ------------------------------------------------------------------------
C       p(e,e'K+)Y
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 700) THEN
        IF(.NOT. BOUND) WRITE(6,3100)
        DECAY_CHECK =
     #       ANSWER(' Kaon Decay Global Cut On? (Y/N)>','Y','y')
C
 3100 FORMAT(' Warning: This option is only Lambda final state ')
C
C ------------------------------------------------------------------------
C       A(e,e'p) KELLY response functions
C ------------------------------------------------------------------------
C
C      ELSEIF(IPHYS_OPT .EQ. 801) THEN
C         CALL GET_LEA_RESPONSE
C
C ------------------------------------------------------------------------
C       A(e,e'p) LAGET response functions
C ------------------------------------------------------------------------
C
C      ELSEIF(IPHYS_OPT .EQ. 802) THEN
C         CALL GET_LAGET_RESPONSE
C
C ------------------------------------------------------------------------
C       A(e,e'p) UDIAS response functions
C ------------------------------------------------------------------------
C
C      ELSEIF(IPHYS_OPT .EQ. 803) THEN
C         CALL GET_UDIAS_RESPONSE
C
C ------------------------------------------------------------------------
C       d(e,e'p)n Arenhoevel unpolarized response functions
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 804) THEN
         WRITE(6,5001)
 5001    FORMAT(' Reading Arenhoevel files ')
         CALL AREN_READ
         WRITE(6,5002)
 5002    FORMAT(' Finished reading Arenhoevel files ')
C
C ------------------------------------------------------------------------
C       d(e,e'p)n Laget unpolarized response functions
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 814) THEN
         WRITE(6,'(A)') ' Linear (1) or Log (2) Pmiss Interpolation? '
         READ(5,*) LAGET_INTP
         WRITE(6,'(A)') 
     #       ' Scatter from Neutron (0), Proton (1) or Both (2)? '
         READ(5,*) LAGET_PWIA
         WRITE(6,'(A)') ' no-FSI (0) or FSI (1)? '
         READ(5,*) LAGET_FSI
         WRITE(6,'(A)') ' no-MEC (0) or MEC (1)? '
         READ(5,*) LAGET_MEC

         WRITE(6,'(A)') ' Start loading the data grid...'  
         CALL GRID_LOAD(SIG_L,SIG_T,SIG_LT,SIG_TT)            
         WRITE(6,'(A)') ' End of data grid load !'  

      ELSE
        GOTO 5          !nonexistent option selected
      ENDIF
C
  999 RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       Subroutine PHYSICS
C
C       AUTHOR:  P.E. ULMER
C       DATE:    AUG-30-1990
C
C       PURPOSE:
C               Get (e,e'N) cross sections and polarizations
C               according to choices made in subroutine PHYS_CHOICE.
C
C               Coincidence cross sections (SIGMA_EEP) are returned
C               with the following units:
C                      fm^2 sr^-2 MeV^-1            (bound state)
C                      fm^2 sr^-2 MeV^-1 (MeV/c)^-1 (continuum)
C               Note that the continuum cross section should be
C               differential in the hadron momentum (NOT kinetic energy)
C               since it is the momentum which is randomly sampled.
C               The recoil factor for the continuum case is equal to
C               the ejectile total energy divided by the momentum.
C               Division by this factor converts the cross section
C               from being differential in energy (as, for example
C               is the case for the deForest CC1 cross section) to
C               differential in momentum.
C ---------------------------------------------------------------------
C
      SUBROUTINE PHYSICS(ELASTIC,BOUND,SPEC_FAC,POL_BEAM,
     #     SIGMA_EEP,SIGMA_MULT_WT,ASYMMETRY,
     #     POL_N,POL_T,POL_L,POL_N_S,POL_T_S,POL_L_S,
     #           POL_X,POL_Y,POL_Z,
     #           POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #           POL_N_RHD,POL_T_RHD,POL_L_RHD,
     #           POL_N_HI,POL_T_HI,POL_L_HI,
     #           POL_N_HD,POL_T_HD,POL_L_HD,
     #           POL_X_HI,POL_Y_HI,POL_Z_HI,
     #           POL_X_HD,POL_Y_HD,POL_Z_HD)
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      COMMON /PHYSVAR/ IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
      COMMON /KINVAR/  KI,Q,Q2,QMU2,EP,PF,THETA_EP,CTH_PQ,PR
      COMMON /SURVIVE_L/ DECAY_CHECK
      COMMON /SURVIVE_D/ SURVIVE_PROB
      COMMON /ARENHOVEL/ FAIL_CELL
      COMMON / GRID_SIG / SIG_L,SIG_T,SIG_LT,SIG_TT
C
      DOUBLE PRECISION KI(3),Q(3),PF(3),PR(3),SURVIVE_PROB
      DOUBLE PRECISION SIGMA_AREN(7),SIGMA_MULT_WT(10)
      DOUBLE PRECISION SIG_L(100,500,91),SIG_T(100,500,91)
      DOUBLE PRECISION SIG_LT(100,500,91), SIG_TT(100,500,91)

      LOGICAL ELASTIC,BOUND,DECAY_CHECK,FAIL_CELL
C
      INCLUDE 'var.cmn'
      INCLUDE 'lifetimes.cmn'
C
      PARAMETER (HBARC = 197.3286D0)
      PARAMETER (PI=3.14159265359D0) 
C
      INCLUDE 'lifetimes_dat.cmn'
C
C ------------------------------------------------------------------------
C       Elastic scattering option
C ------------------------------------------------------------------------
C
      IF(ELASTIC) THEN
         IF(IELAST_OPT .EQ. 10)THEN
            CALL HYD_ELASTIC(E0_I,PF_E_I,TSCAT,POL_BEAM,SIGMA_EEP,
     #                       ASYMMETRY,POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                       POL_N_RHD,POL_T_RHD,POL_L_RHD)
         ELSEIF(IELAST_OPT .EQ. 20)THEN
            CALL DEUT_ELASTIC(E0_I,PF_E_I,TSCAT,POL_BEAM,SIGMA_EEP,
     #                       ASYMMETRY,POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                       POL_N_RHD,POL_T_RHD,POL_L_RHD)
         ELSEIF(IELAST_OPT .EQ. 30)THEN
            CALL TRIT_ELASTIC(E0_I,PF_E_I,TSCAT,POL_BEAM,SIGMA_EEP,
     #                       ASYMMETRY,POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                       POL_N_RHD,POL_T_RHD,POL_L_RHD)
         ELSEIF(IELAST_OPT .EQ. 40)THEN
            CALL HE3_ELASTIC(E0_I,PF_E_I,TSCAT,POL_BEAM,SIGMA_EEP,
     #                       ASYMMETRY,POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                       POL_N_RHD,POL_T_RHD,POL_L_RHD)
         ELSEIF(IELAST_OPT .EQ. 50)THEN
            CALL HE4_ELASTIC(E0_I,PF_E_I,TSCAT,POL_BEAM,SIGMA_EEP,
     #                       ASYMMETRY,POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                       POL_N_RHD,POL_T_RHD,POL_L_RHD)
         ELSEIF(IELAST_OPT .EQ. 60)THEN
            CALL C12_ELASTIC(E0_I,PF_E_I,TSCAT,POL_BEAM,SIGMA_EEP,
     #                       ASYMMETRY,POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                       POL_N_RHD,POL_T_RHD,POL_L_RHD)
         ENDIF
         CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
     #                   POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                   POL_N_HI,POL_T_HI,POL_L_HI,
     #                   POL_X_HI,POL_Y_HI,POL_Z_HI)
         CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
     #                   POL_N_RHD,POL_T_RHD,POL_L_RHD,
     #                   POL_N_HD,POL_T_HD,POL_L_HD,
     #                   POL_X_HD,POL_Y_HD,POL_Z_HD)
         POL_N = POL_N_RHI + POL_N_RHD
         POL_T = POL_T_RHI + POL_T_RHD
         POL_L = POL_L_RHI + POL_L_RHD
         POL_N_S = POL_N_HI + POL_N_HD
         POL_T_S = POL_T_HI + POL_T_HD
         POL_L_S = POL_L_HI + POL_L_HD
         POL_X   = POL_X_HI  + POL_X_HD
         POL_Y   = POL_Y_HI  + POL_Y_HD
         POL_Z   = POL_Z_HI  + POL_Z_HD
         RETURN
      ENDIF
C
C ------------------------------------------------------------------------
C       Sigma_CC1 (unpolarized) X S(p,E)
C ------------------------------------------------------------------------
C
      IF(IPHYS_OPT .EQ. 100) THEN
        CALL OFF_SHELL_D(PF_E_I,Q2,QMU2,EP,PREC,THETA_PQ,
     #             TSCAT,SIGLPT,SIGLT,SIGTT)
        SIGMA_EEP = (SIGLPT+SIGLT*COS(PHI_X)+SIGTT*COS(2.D0*PHI_X))
     #                  *SPECTRAL(BOUND,SPEC_FAC,PRMAG,MISS_M)/RECFAC
C
C ------------------------------------------------------------------------
C       Sigma_CC2 (unpolarized) X S(p,E)
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 101) THEN
        CALL OFF_SHELL_D_CC2(E0_I,OMEGA,PF_P_I,PH_E_I,PH_P_I,
     #             TH_E_I,TH_P_I,PREC,SIGLPT,SIGLT,SIGTT)
        SIGMA_EEP = (SIGLPT+SIGLT*COS(PHI_X)+SIGTT*COS(2.D0*PHI_X))
     #                  *SPECTRAL(BOUND,SPEC_FAC,PRMAG,MISS_M)/RECFAC
C
C ------------------------------------------------------------------------
C       Sigma_CC1 (polarized - Van Orden) X S(p,E)
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 200) THEN
        CALL PWIA_VO(PF_E_I,OMEGA,Q2,QMU2,PF_P_I,PREC,
     #     TSCAT,THETA_PQ,PHI_X,POL_BEAM,SIGMA_EEP,ASYMMETRY,
     #     POL_N_RHI,POL_T_RHI,POL_L_RHI,POL_N_RHD,POL_T_RHD,POL_L_RHD)
        SPECTRAL_FCN = SPECTRAL(BOUND,SPEC_FAC,PRMAG,MISS_M)/RECFAC
        SIGMA_EEP = SIGMA_EEP*SPECTRAL_FCN
        ASYMMETRY = ASYMMETRY*SPECTRAL_FCN
        POL_N_RHI = POL_N_RHI*SPECTRAL_FCN
        POL_T_RHI = POL_T_RHI*SPECTRAL_FCN
        POL_L_RHI = POL_L_RHI*SPECTRAL_FCN
        POL_N_RHD = POL_N_RHD*SPECTRAL_FCN
        POL_T_RHD = POL_T_RHD*SPECTRAL_FCN
        POL_L_RHD = POL_L_RHD*SPECTRAL_FCN
        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
     #                  POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                  POL_N_HI,POL_T_HI,POL_L_HI,
     #                  POL_X_HI,POL_Y_HI,POL_Z_HI)
        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
     #                  POL_N_RHD,POL_T_RHD,POL_L_RHD,
     #                  POL_N_HD,POL_T_HD,POL_L_HD,
     #                  POL_X_HD,POL_Y_HD,POL_Z_HD)
        POL_N   = POL_N_RHI + POL_N_RHD
        POL_T   = POL_T_RHI + POL_T_RHD
        POL_L   = POL_L_RHI + POL_L_RHD
        POL_N_S = POL_N_HI  + POL_N_HD
        POL_T_S = POL_T_HI  + POL_T_HD
        POL_L_S = POL_L_HI  + POL_L_HD
        POL_X   = POL_X_HI  + POL_X_HD
        POL_Y   = POL_Y_HI  + POL_Y_HD
        POL_Z   = POL_Z_HI  + POL_Z_HD
C
C ------------------------------------------------------------------------
C       Gross/Dmitrasinovic Deuterium IA
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 300) THEN
        CALL DEEP_GROSS(QMU2,CTH_PQ,POL_BEAM,
     #     SIGMA_EEP,ASYMMETRY,POL_N_RHI,POL_T_RHI,
     #     POL_L_RHI,POL_N_RHD,POL_T_RHD,POL_L_RHD)
        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
     #                  POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                  POL_N_HI,POL_T_HI,POL_L_HI,
     #                  POL_X_HI,POL_Y_HI,POL_Z_HI)
        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
     #                  POL_N_RHD,POL_T_RHD,POL_L_RHD,
     #                  POL_N_HD,POL_T_HD,POL_L_HD,
     #                  POL_X_HD,POL_Y_HD,POL_Z_HD)
        POL_N   = POL_N_RHI + POL_N_RHD
        POL_T   = POL_T_RHI + POL_T_RHD
        POL_L   = POL_L_RHI + POL_L_RHD
        POL_N_S = POL_N_HI  + POL_N_HD
        POL_T_S = POL_T_HI  + POL_T_HD
        POL_L_S = POL_L_HI  + POL_L_HD
        POL_X   = POL_X_HI  + POL_X_HD
        POL_Y   = POL_Y_HI  + POL_Y_HD
        POL_Z   = POL_Z_HI  + POL_Z_HD
C
C ------------------------------------------------------------------------
C       Sabine Jeschonnek Deuterium PWBA
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 350) THEN
        CALL WQPWBA(E0_I/1000.D0,QMAG/1000.D0,OMEGA/1000.D0,
     #              PRMAG/1000.D0,PF_P_I/1000.D0,
     #              THETA_PQ,PHI_X,THETA_PRQ,PHI_X+PI,
     #              POL_BEAM,SIGMA_EEP)
C
        ASYMMETRY = 0.D0  ! not yet calculated
        POL_N   = 0.d0    ! not yet calculated
        POL_T   = 0.d0    ! not yet calculated
        POL_L   = 0.d0    ! not yet calculated
        POL_N_S = 0.d0    ! not yet calculated
        POL_T_S = 0.d0    ! not yet calculated
        POL_L_S = 0.d0    ! not yet calculated
        POL_X   = 0.d0    ! not yet calculated
        POL_Y   = 0.d0    ! not yet calculated
        POL_Z   = 0.d0    ! not yet calculated
C
C ------------------------------------------------------------------------
C       Enter response functions from files (not yet implemented).
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 400) THEN
        RETURN
C
C ------------------------------------------------------------------------
C       p(e,e'p)pi0
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 500) THEN
        CALL ELECTRO_PROD(POL_BEAM,SIGMA_EEP,ASYMMETRY,
     #                    POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                    POL_N_RHD,POL_T_RHD,POL_L_RHD)
cxxx        CALL epiprod_xsec(W,QMU2_G,THETA_CM,PHI_X,EPSILON,POL_BEAM,
cxxx     #      SIGMA_EEP,ASYMMETRY,POL_N_RHI,POL_T_RHI,
cxxx     #      POL_L_RHI,POL_N_RHD,POL_T_RHD,POL_L_RHD)
        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
     #                  POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #                  POL_N_HI,POL_T_HI,POL_L_HI,
     #                  POL_X_HI,POL_Y_HI,POL_Z_HI)
        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
     #                  POL_N_RHD,POL_T_RHD,POL_L_RHD,
     #                  POL_N_HD,POL_T_HD,POL_L_HD,
     #                  POL_X_HD,POL_Y_HD,POL_Z_HD)
        POL_N   = POL_N_RHI + POL_N_RHD
        POL_T   = POL_T_RHI + POL_T_RHD
        POL_L   = POL_L_RHI + POL_L_RHD
        POL_N_S = POL_N_HI  + POL_N_HD
        POL_T_S = POL_T_HI  + POL_T_HD
        POL_L_S = POL_L_HI  + POL_L_HD
        POL_X   = POL_X_HI  + POL_X_HD
        POL_Y   = POL_Y_HI  + POL_Y_HD
        POL_Z   = POL_Z_HI  + POL_Z_HD
C
C ------------------------------------------------------------------------
C       p(e,e'pi+)n
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 600) THEN
        CALL PEEPI(SIGMA_EEP)
        IF(DECAY_CHECK) THEN
          CALL DECAY(LIFE_PION)
          SIGMA_EEP = SIGMA_EEP*SURVIVE_PROB  !fold in pion survival prob.
        ENDIF
C
C ------------------------------------------------------------------------
C       p(e,e'K+)Y
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 700) THEN
        CALL PEEK(SIGMA_EEP)
        IF(DECAY_CHECK) THEN
          CALL DECAY(LIFE_KAON)
          SIGMA_EEP = SIGMA_EEP*SURVIVE_PROB  !fold in kaon survival prob.
        ENDIF
C
C ------------------------------------------------------------------------
C       LEA Response Functions (J.J. Kelly)
C ------------------------------------------------------------------------
C
C      ELSEIF(IPHYS_OPT .EQ. 801) THEN
C        CALL LEA_PHYSICS(PF_E_I,OMEGA,Q2,QMU2,PF_P_I,PREC,
C     #     TSCAT,THETA_PQ,PHI_X,RECFAC,POL_BEAM,SIGMA_EEP,ASYMMETRY,
C     #     POL_N_RHI,POL_T_RHI,POL_L_RHI,POL_N_RHD,POL_T_RHD,POL_L_RHD,
C     #     POL_X_HI, POL_Y_HI, POL_Z_HI, POL_X_HD, POL_Y_HD, POL_Z_HD)
C        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
C     #                  POL_N_RHI,POL_T_RHI,POL_L_RHI,
C     #                  POL_N_HI,POL_T_HI,POL_L_HI,
C     #                  POL_X_HI,POL_Y_HI,POL_Z_HI)
C        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
C     #                  POL_N_RHD,POL_T_RHD,POL_L_RHD,
C     #                  POL_N_HD,POL_T_HD,POL_L_HD,
C     #                  POL_X_HD,POL_Y_HD,POL_Z_HD)
C        POL_N   = POL_N_RHI + POL_N_RHD
C        POL_T   = POL_T_RHI + POL_T_RHD
C        POL_L   = POL_L_RHI + POL_L_RHD
C        POL_N_S = POL_N_HI  + POL_N_HD
C        POL_T_S = POL_T_HI  + POL_T_HD
C        POL_L_S = POL_L_HI  + POL_L_HD
C        POL_X   = POL_X_HI  + POL_X_HD
C        POL_Y   = POL_Y_HI  + POL_Y_HD
C        POL_Z   = POL_Z_HI  + POL_Z_HD
C
C ------------------------------------------------------------------------
C       LAGET Response Functions (J.M.Laget)
C ------------------------------------------------------------------------
C
C      ELSEIF(IPHYS_OPT .EQ. 802) THEN
C        CALL LAGET_PHYSICS(PF_E_I,OMEGA,Q2,QMU2,PF_P_I,PREC,
C     #     TSCAT,THETA_PQ,PHI_X,W,RECFAC,POL_BEAM,SIGMA_EEP,ASYMMETRY,
C     #     POL_N_RHI,POL_T_RHI,POL_L_RHI,POL_N_RHD,POL_T_RHD,POL_L_RHD)
C
C        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
C     #                  POL_N_RHI,POL_T_RHI,POL_L_RHI,
C     #                  POL_N_HI,POL_T_HI,POL_L_HI,
C     #                  POL_X_HI,POL_Y_HI,POL_Z_HI)
C        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
C     #                  POL_N_RHD,POL_T_RHD,POL_L_RHD,
C     #                  POL_N_HD,POL_T_HD,POL_L_HD,
C     #                  POL_X_HD,POL_Y_HD,POL_Z_HD)
C
C        POL_N   = POL_N_RHI + POL_N_RHD
C        POL_T   = POL_T_RHI + POL_T_RHD
C        POL_L   = POL_L_RHI + POL_L_RHD
C        POL_N_S = POL_N_HI  + POL_N_HD
C        POL_T_S = POL_T_HI  + POL_T_HD
C        POL_L_S = POL_L_HI  + POL_L_HD
C        POL_X   = POL_X_HI  + POL_X_HD
C        POL_Y   = POL_Y_HI  + POL_Y_HD
C        POL_Z   = POL_Z_HI  + POL_Z_HD
C
C ------------------------------------------------------------------------
C       UDIAS Response Functions (J.M. Udias)
C ------------------------------------------------------------------------
C
C      ELSEIF(IPHYS_OPT .EQ. 803) THEN
C        CALL UDIAS_PHYSICS(PF_E_I,OMEGA,Q2,QMU2,PF_P_I,PREC,
C     #     TSCAT,THETA_PQ,PHI_X,W,RECFAC,POL_BEAM,SIGMA_EEP,ASYMMETRY,
C     #     POL_N_RHI,POL_T_RHI,POL_L_RHI,POL_N_RHD,POL_T_RHD,POL_L_RHD)
C
C        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
C     #                  POL_N_RHI,POL_T_RHI,POL_L_RHI,
C     #                  POL_N_HI,POL_T_HI,POL_L_HI,
C     #                  POL_X_HI,POL_Y_HI,POL_Z_HI)
C        CALL ROTATE_POL(ELASTIC,PHI_X,THETA_PQ,KI,Q,PF,
C     #                  POL_N_RHD,POL_T_RHD,POL_L_RHD,
C     #                  POL_N_HD,POL_T_HD,POL_L_HD,
C     #                  POL_X_HD,POL_Y_HD,POL_Z_HD)
C
C        POL_N   = POL_N_RHI + POL_N_RHD
C        POL_T   = POL_T_RHI + POL_T_RHD
C        POL_L   = POL_L_RHI + POL_L_RHD
C        POL_N_S = POL_N_HI  + POL_N_HD
C        POL_T_S = POL_T_HI  + POL_T_HD
C        POL_L_S = POL_L_HI  + POL_L_HD
C        POL_X   = POL_X_HI  + POL_X_HD
C        POL_Y   = POL_Y_HI  + POL_Y_HD
C        POL_Z   = POL_Z_HI  + POL_Z_HD
C
C ------------------------------------------------------------------------
C       ARENHOEVEL Unpolarized Response Functions (H. Arenhoevel)
C
C       Interpolation is done on ratio AREN/PWIA which should be
C       much flatter than the AREN cross section. 
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 804) THEN
         ISPEC_OPT = 10     ! This is for PWIA
         CALL AREN_INTERP(fail_cell,sigma_aren)
         IF(.NOT. FAIL_CELL) THEN
            SIGMA_EEP = SIGMA_AREN(7)   ! Default is "FULL" calculation
            DO I=1,7
               SIGMA_MULT_WT(I) = SIGMA_AREN(I)
            ENDDO
         ELSE
            SIGMA_EEP = 0.D0
            DO I=1,7
               SIGMA_MULT_WT(I) = 0.D0
            ENDDO
         ENDIF
C
C ------------------------------------------------------------------------
C       Laget Unpolarized Response Functions (J.-M. Laget)
C ------------------------------------------------------------------------
C
      ELSEIF(IPHYS_OPT .EQ. 814) THEN
        CALL LAGET_XSEC(SIGMA_EEP,E0_I,QMU2,OMEGA,THETA_CM,PHI_X,
     +                   SIG_L,SIG_T,SIG_LT,SIG_TT)
        SIGMA_EEP = SIGMA_EEP * 1.D-4   ! ub -> fm^2

      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE GET_RESPONSE
C
C       AUTHOR:  P.E. ULMER
C       DATE:    AUG-30-1990
C
C       PURPOSE: Get response functions for (e,e'N).
C
C       The user needs to specify weighting factors so that
C       cross sections and polarizations can be determined from
C       the response functions.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_RESPONSE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      STOP
      END










