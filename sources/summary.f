C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE SUMMARY
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    15-OCT-1990
C       PURPOSE: Write summary of MCEEP results to a file.
C -----------------------------------------------------------------------
C
      SUBROUTINE SUMMARY(ELASTIC,BOUND,RADPEAKING,MULTIPHOTON,RADFULL,
     #                   ACCEPT_CHECK)
      IMPLICIT NONE
C
      COMMON /PHYSVAR/ IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
      COMMON /PARTICLE/ PROTON
      COMMON /RADPKAPPR/ CUTOFF,EG_MAX
      COMMON /ARENFAIL/ I_AREN_FAIL
      COMMON / LAGET_CNT/ LAGET_PS_FAIL,LAGET_GRID_FAIL
      COMMON / LAGET_MOD/ LAGET_INTP,LAGET_PWIA,LAGET_FSI,LAGET_MEC
C
      REAL VERSION_NUMBER
C
      DOUBLE PRECISION PI,Q_E,N_A,CUTOFF,EG_MAX
      DOUBLE PRECISION SUM,CENT_X,CENT_Y,FWHM_SIG
      DOUBLE PRECISION ARR(500),ARR_2D(50,50),KI_VEC(3),KF_VEC(3)
C
      INTEGER IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL,I_WC,NCALC
      INTEGER I_TR,IA,IEL,NCHAR,IARM,I,J,K,I_AREN_FAIL
      INTEGER LAGET_PS_FAIL,LAGET_GRID_FAIL
      INTEGER LAGET_INTP,LAGET_PWIA,LAGET_FSI,LAGET_MEC
C
      LOGICAL ELASTIC,BOUND,RADPEAKING,MULTIPHOTON,FAIL_KIN,PROTON
      LOGICAL ACCEPT_CHECK,LG,RADFULL
C
      CHARACTER*300 TMP_FILE
      CHARACTER*300 PL_F,TR_F,WC_F
C
C -----------------------------------------------------------------------
C       Include various common blocks
C -----------------------------------------------------------------------
C
      INCLUDE 'summary.cmn'
      INCLUDE 'var.cmn'
      INCLUDE 'input.cmn'
      INCLUDE 'spectrometer.cmn'
      INCLUDE 'masses.cmn'
      INCLUDE 'wc.cmn'
      INCLUDE 'eloss.cmn'
C
      PARAMETER (Q_E = 1.602D0)               !e- charge (uCoul*1E-13)
      PARAMETER (N_A = 6.022D0)               !Avogadro's number (*1E-23)
      DATA VERSION_NUMBER /3.90/              !Program version number
      PI = ASIN(1.D0)*2.D0
      FWHM_SIG = 2.D0*SQRT(2.D0*LOG(2.D0))    !FWHM --> gaussian sigma
C
C ---------------------------------------------------------------------
C       Form centroids and averages for each variable in VAR.CMN
C ---------------------------------------------------------------------
C
      DO I=1,NUM_VAR
        IF(SUM_EEP_C.NE.0.) VAR_CENT(I) = VAR_CENT(I)/SUM_EEP_C
        IF(SUM_EV_C .NE.0.) VAR_AVG(I)  = VAR_AVG(I)/SUM_EV_C
      ENDDO
C
C ---------------------------------------------------------------------
C       Open summary file.
C ---------------------------------------------------------------------
C
      TMP_FILE = FILEPRE//'.sum'
      CALL SQUEEZE(TMP_FILE,SUMFILE,NCHAR)
      OPEN(UNIT=1,FILE=SUMFILE,STATUS='UNKNOWN',FORM='FORMATTED')
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C       WRITE SUMMARY FILE (INPUT SECTION).
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      WRITE(1,7000) VERSION_NUMBER
 7000 FORMAT(//,T22,' MCEEP PROGRAM VERSION NUMBER: ',2X,F5.2)
      WRITE(1,8001)
 8001 FORMAT(///,T22,' I N P U T     S U M M A R Y '//)
C
      IF(NCOMMENTS .GT. 0) THEN
        DO I=1,NCOMMENTS
          WRITE(1,8060) COMMENT(I)
 8060     FORMAT(1X,A79)
        ENDDO
      ENDIF
C
C      CALL SQUEEZE(INFILE,INFILE,NCHAR)       !Extract non-blank characters
      WRITE(1,8002)INFILE(1:NCHAR)
 8002 FORMAT(//' Input file: ',A)
C
      IF(RADPEAKING) THEN
         NCALC = 2
      ELSE
         NCALC = 1
      ENDIF
C
      WRITE(1,8100) NITER1,NITER2,NITER3,NITER4,NITER5,NITER6
 8100 FORMAT(/' Number of iterations for RANGES preprocessor: ',/,
     #         T10,' Theta_p:         ',T35,I3,/,
     #         T10,' Phi_p:           ',T35,I3,/,
     #         T10,' Theta_e:         ',T35,I3,/,
     #         T10,' Phi_e:           ',T35,I3,/,
     #         T10,' e- momentum:     ',T35,I3,/,
     #         T10,' hadron momentum: ',T35,I3)
C
C ---------------------------------------------------------------------
C       Radiative effects included?
C       Elastic scattering, bound state or continuum?
C       Write appropriate information.
C ---------------------------------------------------------------------
C
      IF(RADPEAKING) THEN
        WRITE(1,2008) CUTOFF
 2008   FORMAT(//' INTERNAL BREMSSTRAHLUNG INCLUDED, CUTOFF = ',
     #           F9.3,' MeV ')
        IF(MULTIPHOTON) THEN
           WRITE(1,*) '    Multiphoton correction included '
        ELSE
           WRITE(1,*) '    Multiphoton correction NOT included '
        ENDIF
      ELSEIF(RADFULL) THEN
        WRITE(1,2015)
 2015   FORMAT(//' INTERNAL BREM. INCLUDED WITH FULL ANG. DIST. ')
      ELSE
        WRITE(1,2009)
 2009   FORMAT(//' INTERNAL BREMSSTRAHLUNG NOT INCLUDED ')
      ENDIF
C
      IF(ELASTIC) THEN
        WRITE(1,2010)
 2010   FORMAT(//' ELASTIC SCATTERING: ')
        IF(IELAST_OPT .EQ. 10) THEN
           WRITE(1,2011)
 2011      FORMAT(T15,' from proton ')
        ELSEIF(IELAST_OPT .EQ. 20) THEN
           WRITE(1,2012)IMODEL
 2012      FORMAT(T15,' from deuteron, model # = ',I2)
        ENDIF
        IF(ACCEPT_CHECK) THEN
           WRITE(1,2013)
 2013      FORMAT(T15,' Hadron arm acceptances enforced ')
        ELSE
           WRITE(1,2014)
 2014      FORMAT(T15,' Hadron arm acceptances NOT enforced ')
        ENDIF
      ELSEIF(BOUND) THEN
        WRITE(1,8004)EM_BOUND
 8004   FORMAT(//' BOUND STATE ',2X,' E_m = ',F10.3, ' MeV ')
      ELSE
        WRITE(1,8005)
 8005   FORMAT(//' CONTINUUM ')
      ENDIF
C
      IF(.NOT. ELASTIC) THEN
         WRITE(1,8080)IPHYS_OPT
 8080    FORMAT(/' Physics Option selected: ',T30,I4)
         IF(ISPEC_OPT .NE. 0) WRITE(1,8081)ISPEC_OPT
 8081    FORMAT(' Spectral Function: ',T30,I4)
C
         IF(IPHYS_OPT .NE. 600) THEN
           IF(PROTON) THEN
             WRITE(1,8082)
 8082        FORMAT(' (e,e''p) selected ')
           ELSE
             WRITE(1,8083)
 8083        FORMAT(' (e,e''n) selected ')
           ENDIF
         ELSE
           WRITE(1,8084)
 8084      FORMAT(' (e,e''pi+) selected ')
         ENDIF
      ENDIF
C
      IF(IPHYS_OPT .EQ. 814) THEN   ! Laget d(e,ep)n grid
         WRITE(1,8085) LAGET_INTP,LAGET_PWIA,LAGET_FSI,LAGET_MEC
 8085    FORMAT(' Laget Choices: ',
     #        /,'      LAGET_INTP = ',I3,
     #        /,'      LAGET_PWIA = ',I3,
     #        /,'      LAGET_FSI  = ',I3,
     #        /,'      LAGET_MEC  = ',I3,//)
      ENDIF
C
C ---------------------------------------------------------------------
C       Write masses/charges
C ---------------------------------------------------------------------
C
      WRITE(1,8090)MNUC,EJECT_MASS,EJECT_CHG
 8090 FORMAT(//' MASSES: ',/,
     #   T10,' Target mass      = ',T30,F12.3,' MeV ',/,
     #   T10,' Ejectile mass    = ',T30,F12.3,' MeV ',/,
     #   T10,' Ejectile charge  = ',T30,I3)
C
C ---------------------------------------------------------------------
C       Write target parameters.
C ---------------------------------------------------------------------
C
      WRITE(1,8092)INT(ATARG),INT(ZTARG),DENS_TARG,ITARG_MOD,
     #     PHI_TARG*180./PI,
     #     NFOILS,(J,TARG_LO(J),TARG_HI(J),J=1,NFOILS)
 8092 FORMAT(//' TARGET PARAMETERS: ',/,
     #   T10,' Target A       = ',T30,I4,/,
     #   T10,' Target Z       = ',T30,I4,/,
     #   T10,' Target Density = ',T30,F8.4,' g/cm^3 ',/,
     #   T10,' Target Model   = ',T30,I4,/,
     #   T10,' Target Angle   = ',T30,F8.2,' deg ',/,
     #   T10,' # Foils        = ',T30,I4,//,
     #   T10,' Foil ',T20,' Start (m) ',T35,' End (m) ',/,
     #   T10,' ---- ',T20,' --------- ',T35,' ------- ',/,
     #   10(T12,I2,T20,F8.4,T35,F8.4,/))
C
      IF(ELOSS_MOD .EQ. 0) THEN
         WRITE(1,'(A)') ' Energy loss not included '
      ELSEIF(ELOSS_MOD .EQ. 1) THEN
         WRITE(1,'(A)') ' Energy loss included '
      ELSEIF(ELOSS_MOD .EQ. 2) THEN
         WRITE(1,'(A)') ' Energy loss included: mean dE removed '
      ELSEIF(ELOSS_MOD .EQ. 3) THEN
         WRITE(1,'(A)') ' Energy loss included: most prob dE removed '
      ELSE
         WRITE(1,'(A)') ' Invalid energy loss model selected '
      ENDIF
C
C ---------------------------------------------------------------------
C       Write kinematics.
C ---------------------------------------------------------------------
C
      WRITE(1,8006)E0,PH_B*180./PI,TH_B*180./PI,
     #               PF_E,PH_E*180./PI,TH_E*180./PI,
     #               PF_P,PH_P*180./PI,TH_P*180./PI
 8006 FORMAT(//' KINEMATICS ',
     #   T17,' Momentum ',T35,' Phi ',T50,' Theta ',/,
     #   T17,'   MeV/c  ',T35,' deg ',T50,'  deg  ',/,
     #   T17,' -------- ',T35,' --- ',T50,' ----- ',/,
     #   ' Beam:     ',T17,F9.3,T30,F10.3,T46,F10.3,/,
     #   ' Electron: ',T17,F9.3,T30,F10.3,T46,F10.3,/,
     #   ' Hadron:   ',T17,F9.3,T30,F10.3,T46,F10.3,/)
C
C ---------------------------------------------------------------------
C       Write acceptances.
C ---------------------------------------------------------------------
C
      WRITE(1,8007)ACC_EP_N,ACC_EM_N,DPH_E_N*1000.,DTH_E_N*1000.,
     #             SA_SHAPE_E,
     #             ACC_PP_N,ACC_PM_N,DPH_P_N*1000.,DTH_P_N*1000.,
     #             SA_SHAPE_P
 8007 FORMAT(///' ACCEPTANCES ',
     #   T17,' Momentum ',T36,' Phi ',T47,' Theta ',T58,' Aperture ',/,
     #   T17,'   +/- %  ',T36,'  mr ',T47,'   mr  ',T58,'  Shape   ',/,
     #   T17,' -------- ',T36,' --- ',T47,' ----- ',T58,' -------- ',/,
     #   ' Electron: ',T14,F6.1,1X,F6.1,T33,F8.2,T45,F8.2,T62,A1,/,
     #   ' Hadron:   ',T14,F6.1,1X,F6.1,T33,F8.2,T45,F8.2,T62,A1,/)
C
C ---------------------------------------------------------------------
C       Write rate parameters, etc.
C ---------------------------------------------------------------------
C
      WRITE(1,8008)ALUM*(N_A/(Q_E*ATARG)),BEAMTIME/3600.,SPEC_FAC
 8008 FORMAT(//,  T10,' Luminosity     = ',T30,E12.3,
     #                                   ' [10^(36) cm^-2 sec^-1] ',/,
     #            T10,' Beam time      = ',T30,F12.3,' hrs ',/,
     #            T10,' Spectr. factor = ',T30,F12.3,/)
C
C ---------------------------------------------------------------------
C       Write drift parameters.
C ---------------------------------------------------------------------
C
      WRITE(1,8061) DRIFTE_N,DRIFTP_N
 8061 FORMAT(//,T15,' DRIFTS TO APERTURES ',/,
     #          T15,' ------ -- --------- ',/,
     #          T10,' Drift to E aperture = ',T40,F12.3,'  m ',/,
     #          T10,' Drift to P aperture = ',T40,F12.3,'  m ',/)
C
C ---------------------------------------------------------------------
C       Write beam parameters.
C ---------------------------------------------------------------------
C
      WRITE(1,8062)POL_BEAM,BEAMV*1000.,BEAMD*100.,
     #             DUTY_FAC*100.,DEL_TOF*1.D9,
     #             FWHM_X_B*100.,FWHM_Y_B*100.,FWHM_PH_B*1000.,
     #             FWHM_TH_B*1000.,100.D0*FWHM_E0/E0,
     #             DEL_X_B*100.,DEL_Y_B*100.,DEL_PH_B*1000.,
     #             DEL_TH_B*1000.,100.D0*DEL_E0/E0
 8062 FORMAT(//,T30,' BEAM PARAMETERS ',/,
     #          T30,' ---- ---------- ',/,
     #     T10,' Polarization                 = ',T50,F12.3,/,
     #     T10,' Eff. vertical spread at tgt. = ',T50,F12.3,' mm ',/,
     #     T10,' Total dispersion at tgt.     = ',T50,F12.3,' %  ',/,
     #     T10,' Duty factor                  = ',T50,F12.3,' %  ',/,
     #     T10,' ToF window                   = ',T50,F12.3,' nsec ',//,
     #     T10,' FWHM for               X_B   = ',T50,F12.3,' cm ',/,
     #     T10,'                        Y_B   = ',T50,F12.3,' cm ',/,
     #     T10,'                        PH_B  = ',T50,F12.3,' mr ',/,
     #     T10,'                        TH_B  = ',T50,F12.3,' mr ',/,
     #     T10,'                        E0    = ',T50,F12.3,' %  ',//,
     #     T10,' Offset for             X_B   = ',T50,F12.3,' cm ',/,
     #     T10,'                        Y_B   = ',T50,F12.3,' cm ',/,
     #     T10,'                        PH_B  = ',T50,F12.3,' mr ',/,
     #     T10,'                        TH_B  = ',T50,F12.3,' mr ',/,
     #     T10,'                        E0    = ',T50,F12.3,' %  ',/)
      WRITE(1,8063)RAST_SHAPE,RAST_AMP_X,RAST_AMP_Y
 8063 FORMAT( T10,' RASTER          Shape        = ',T56,A1,/,
     #     T10,'                 Full-width X = ',T50,F12.5,' m  ',/,
     #     T10,'                 Full-width Y = ',T50,F12.5,' m  ',/)
C
C ---------------------------------------------------------------------
C       Write spectrometer analysis input.
C ---------------------------------------------------------------------
C
      WRITE(1,8009)
 8009 FORMAT(//,72('-'),/,
     #         '  ELECTRON ARM SPECTROMETER ELEMENTS: ',/)
      WRITE(1,8150) THETA_BP(1)*180.D0/PI,
     #              (OBJECT_PT(1,J)*100.D0,J=1,3)
 8150 FORMAT(T16,' Bend-plane angle = ',F8.2,' deg ',/,
     #       T16,' Coordinates of Object point in LAB system (cm): '
     #      ,/,T25,' X:  ',F8.3
     #      ,/,T25,' Y:  ',F8.3
     #      ,/,T25,' Z:  ',F8.3)
      IF(APPLY_TO_LAB(1)) THEN
         WRITE(1,8070)
 8070    FORMAT(T16,' Apply to LAB coordinates?  Yes ')
      ELSE
         WRITE(1,8071)
 8071    FORMAT(T16,' Apply to LAB coordinates?  No ')
      ENDIF
      WRITE(1,8010)
 8010 FORMAT(/,T1 ,' Type/ID',T12,'Affects ',T22,'X-Var',T28,'Y-Var',
     #           T35,' Cut ranges/Axes ',T54,'Value',
     #           T63,'Filename',/,
     #           T1 ,' -------',T12,'------- ',T22,'-----',T28,'-----',
     #           T35,' --------------- ',T54,'-----',
     #           T63,'--------')
      I_TR = 0
      DO IA=1,2
        IF(IA .EQ. 2)THEN
          WRITE(1,8011)
 8011     FORMAT(//,72('-'),/,
     #             '  HADRON ARM SPECTROMETER ELEMENTS: ',/)
          WRITE(1,8150) THETA_BP(2)*180.D0/PI,
     #              (OBJECT_PT(2,J)*100.D0,J=1,3)
          IF(APPLY_TO_LAB(2)) THEN
             WRITE(1,8070)
          ELSE
             WRITE(1,8071)
          ENDIF
          WRITE(1,8010)
        ENDIF
        IF(NEL(IA) .GT. 0) THEN
          DO IEL=1,NEL(IA)
            IF(OP(IA,IEL).EQ.'HRS') THEN
              WRITE(1,8091)OP(IA,IEL),HRS_ID(IA,IEL),OB_NAM(IA,IEL),
     #                    (APERTURE_TEST(HRS_ID_INT(IA,IEL),J),J=1,5)
 8091         FORMAT(T3,A3,A1,T12,A8,T54,5L1)
C HRS + Septum choice 
            ELSEIF(OP(IA,IEL).EQ.'SEP') THEN
              WRITE(1,8094)OP(IA,IEL),HRS_ID(IA,IEL),OB_NAM(IA,IEL),
     #                    (APER_SEP_TEST(HRS_ID_INT(IA,IEL),J),J=1,10)
 8094         FORMAT(T3,A3,A1,T12,A8,T54,10L1)
C COSY choice 
            ELSEIF(OP(IA,IEL).EQ.'SIM') THEN
              WRITE(1,8091)OP(IA,IEL),HRS_ID(IA,IEL),OB_NAM(IA,IEL),
     #                    (HUT_TEST(HRS_ID_INT(IA,IEL),J),J=1,3)
            ELSEIF(OP(IA,IEL).EQ.'HRI') THEN
              WRITE(1,8091)OP(IA,IEL),HRS_ID(IA,IEL),OB_NAM(IA,IEL)
            ELSEIF(OP(IA,IEL).EQ.'SPI') THEN
              WRITE(1,8094)OP(IA,IEL),HRS_ID(IA,IEL),OB_NAM(IA,IEL)
            ELSEIF(OP(IA,IEL).EQ.'MAD') THEN
              IF(MAD_CONFIG .EQ. 1) THEN
                WRITE(1,8095)OP(IA,IEL),MAD_CONFIG,OB_NAM(IA,IEL),
     #                    (APER_MAD_TEST(J),J=1,4)
              ELSEIF(MAD_CONFIG .EQ. 2 .OR. MAD_CONFIG .EQ. 3) THEN
                WRITE(1,8095)OP(IA,IEL),MAD_CONFIG,OB_NAM(IA,IEL),
     #                    (APER_MAD_TEST(J),J=1,6)
              ENDIF
 8095         FORMAT(T3,A3,I1,T12,A8,T54,6L1)
            ELSEIF(OP(IA,IEL).EQ.'MDI') THEN
              WRITE(1,8095)OP(IA,IEL),MAD_INV_CONFIG,OB_NAM(IA,IEL)
            ELSEIF(OP(IA,IEL).EQ.'RFN') THEN
              WRITE(1,8093)OP(IA,IEL),HRS_ID(IA,IEL),OB_NAM(IA,IEL),
     #                    HRS_COLL(HRS_ID_INT(IA,IEL)),
     #                    RFN_CUT(HRS_ID_INT(IA,IEL))
 8093         FORMAT(T3,A3,A1,T12,A8,T36,A1,T38,F12.5)
            ELSEIF(OP(IA,IEL).EQ.'OFF') THEN
              WRITE(1,8012)OP(IA,IEL),OB_NAM(IA,IEL),
     #                    TR_LABEL(NCOORD(IA,IEL,1)),OFF(IA,IEL)
 8012         FORMAT(T3,A3,T12,A8,T22,A5,T53,F7.3)
            ELSEIF(OP(IA,IEL).EQ.'RES') THEN
              WRITE(1,8013)OP(IA,IEL),OB_NAM(IA,IEL),
     #                    TR_LABEL(NCOORD(IA,IEL,1)),GSIG(IA,IEL)
 8013         FORMAT(T3,A3,T12,A8,T22,A5,T53,F7.3)
            ELSEIF(OP(IA,IEL).EQ.'MAT') THEN
              CALL SQUEEZE(TRNSPT_FIL(IA,IEL),TR_F,NCHAR)
              IF(INVERT(IA,IEL))THEN
                TR_F = '(Inv) '//TR_F(1:NCHAR)
                NCHAR = NCHAR + 6
              ELSE
                TR_F = TR_F(1:NCHAR)
              ENDIF
              WRITE(1,8030)OP(IA,IEL),ORDER(IA,IEL),TR_F(1:NCHAR)
 8030         FORMAT(T3,A3,1X,I1,T63,A)
            ELSEIF(OP(IA,IEL).EQ.'ROT') THEN
              WRITE(1,8031)OP(IA,IEL),ROT_AXIS(IA,IEL),
     #                     ROT_ANGLE(IA,IEL)*180.D0/PI
 8031         FORMAT(T3,A3,2X,A1,T53,F7.3)
            ELSEIF(OP(IA,IEL).EQ.'DFT') THEN
              WRITE(1,8032)OP(IA,IEL),DRIFT_LENGTH(IA,IEL)
 8032         FORMAT(T3,A3,T53,F7.2)
            ELSEIF(OP(IA,IEL).EQ.'TRK') THEN
              CALL SQUEEZE(WC_FIL(IA),WC_F,NCHAR)
              WRITE(1,8014)OP(IA,IEL),WC_F(1:NCHAR)
 8014         FORMAT(T3,A3,T63,A)
            ELSEIF(OP(IA,IEL).EQ.'POL') THEN
              WRITE(1,8032)OP(IA,IEL),BEND_ANGLE(IA,IEL)*180.D0/PI
            ELSEIF(OP(IA,IEL).EQ.'TOF') THEN
              WRITE(1,8032)OP(IA,IEL),CRAY_PATH(IA)
            ELSEIF(OP(IA,IEL).EQ.'CUT') THEN
              IF(CUT_DOMAIN(IA,IEL).EQ.'G' .AND.
     #             CUT_TYPE(IA,IEL).EQ.'R') THEN
                 WRITE(1,8015)OP(IA,IEL),'All hist',
     #                TR_LABEL(NCOORD(IA,IEL,1)),TR_CUT_MIN(IA,IEL),
     #                TR_CUT_MAX(IA,IEL)
 8015            FORMAT(T3,A3,T12,A8,T22,A5,T33,F7.2,T41,' to',T45,F7.2)
              ELSEIF(CUT_DOMAIN(IA,IEL).EQ.'G' .AND.
     #                 CUT_TYPE(IA,IEL).EQ.'E') THEN
                 WRITE(1,8016)OP(IA,IEL),'All hist',
     #                TR_LABEL(NCOORD(IA,IEL,1)),
     #                TR_LABEL(NCOORD(IA,IEL,2)),
     #                TR_CUT_X(IA,IEL),TR_CUT_Y(IA,IEL)
 8016            FORMAT(T3,A3,T12,A8,T22,A5,T28,A5,T33,F7.2,T41,' and',
     #                  T45,F7.2)
              ELSEIF(CUT_DOMAIN(IA,IEL).EQ.'S' .AND.
     #                 CUT_TYPE(IA,IEL).EQ.'R') THEN
                 WRITE(1,8017)OP(IA,IEL),TR_CUT_IND(IA,IEL),
     #                TR_LABEL(NCOORD(IA,IEL,1)),TR_CUT_MIN(IA,IEL),
     #                TR_CUT_MAX(IA,IEL)
 8017            FORMAT(T3,A3,1X,I2,T22,A5,T33,F7.2,T41,' to',T45,F7.2)
              ELSEIF(CUT_DOMAIN(IA,IEL).EQ.'S' .AND.
     #                 CUT_TYPE(IA,IEL).EQ.'E') THEN
                 WRITE(1,8018)OP(IA,IEL),TR_CUT_IND(IA,IEL),
     #                TR_LABEL(NCOORD(IA,IEL,1)),
     #                TR_LABEL(NCOORD(IA,IEL,2)),
     #                TR_CUT_X(IA,IEL),TR_CUT_Y(IA,IEL)
 8018            FORMAT(T3,A3,1X,I2,T22,A5,T28,A5,T33,F7.2,T41,' and',
     #                  T45,F7.2)
              ENDIF
            ELSEIF(OP(IA,IEL).EQ.'H1D') THEN
              I_TR = I_TR + 1
              CALL SQUEEZE(TR_FILE(IA,IEL),TR_F,NCHAR)
              WRITE(1,8019)OP(IA,IEL),I_TR,TR_LABEL(NCOORD(IA,IEL,1)),
     #              TR_F(1:NCHAR)
 8019         FORMAT(T3,A3,1X,I2,T22,A5,T63,A)
            ELSEIF(OP(IA,IEL).EQ.'H2D') THEN
              I_TR = I_TR + 1
              CALL SQUEEZE(TR_FILE(IA,IEL),TR_F,NCHAR)
              WRITE(1,8020)OP(IA,IEL),I_TR,TR_LABEL(NCOORD(IA,IEL,1)),
     #                TR_LABEL(NCOORD(IA,IEL,2)),TR_F(1:NCHAR)
 8020         FORMAT(T3,A3,1X,I2,T22,A5,T28,A5,T63,A)
            ELSEIF(OP(IA,IEL).EQ.'SCT') THEN
              I_TR = I_TR + 1
              CALL SQUEEZE(TR_FILE(IA,IEL),TR_F,NCHAR)
              WRITE(1,8021)OP(IA,IEL),I_TR,TR_LABEL(NCOORD(IA,IEL,1)),
     #                TR_LABEL(NCOORD(IA,IEL,2)),TR_F(1:NCHAR)
 8021         FORMAT(T3,A3,1X,I2,T22,A5,T28,A5,T63,A)
            ELSEIF(OP(IA,IEL).EQ.'NTU') THEN
              I_TR = I_TR + 1
              CALL SQUEEZE(TR_FILE(IA,IEL),TR_F,NCHAR)
              WRITE(1,8088)OP(IA,IEL),I_TR,BEND_ANGLE_TOT(IA,IEL)
     #                *180.D0/PI,TR_F(1:NCHAR)
 8088         FORMAT(T3,A3,1X,I2,T53,F7.3,T63,A)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      WRITE(1,9022)
 9022 FORMAT(//,72('-'),/,
     #      T18,' WIRE CHAMBER PARAMETERS FOR TRK OPTION ',//,
     #    ' IARM',T8,' <x|delta>',T20,' Offset',T29,' Coord_off',T40,
     #    ' Chmbr_off',T51,' Resol.',T60,' # Chmbr'
     #    ,/,T9,' (cm/%)',T19,' (mr,cm,%)',T50,' FWHM-cm'
     #    ,/,' ----',T8,' ---------',T20,' ------',T29,' ---------',
     #    T40,' ---------',T51,' -------',T60,' -------')
      DO IARM=1,2
        WRITE(1,9023)IARM,X_DELTA(IARM),WC_OFF(IARM),WC_OFF_COORD(IARM),
     #         WC_OFF_CHMB(IARM),WC_RES(IARM)*FWHM_SIG,N_WC(IARM)
      ENDDO
 9023 FORMAT(2X,I2,3X,F8.2,4X,F6.2,8X,I2,9X,I2,2X,F8.2,7X,I2)
      WRITE(1,9024)
 9024 FORMAT(//,T11,' IARM',T19,' Chamber #',T30,' Location',T41,
     #    ' Tilt angle',T53,' Wire angle',T65,' Sigma Mscat'
     #    ,/,T32,' (cm)',T43,' (deg)',T56,' (deg)',T68,' (mr)'
     #    ,/,T11,' ----', T19,' ---------'
     #    ,T30,' ---------',T41,' ----------',T53,' ----------'
     #    ,T65,' -----------')
      DO IARM=1,2
         IF (N_WC(IARM) .EQ. 0) THEN
           WRITE(1,9025)IARM,N_WC(IARM)
         ELSE
           LG = .TRUE.
           DO I_WC = 1,N_WC(IARM)
              IF (LG) THEN
                WRITE(1,9025)IARM,I_WC,WC_LOC(IARM,I_WC),
     #                WC_C_ANGLE(IARM,I_WC),WC_W_ANGLE(IARM,I_WC),
     #                WC_MSCT(IARM,I_WC)
              ELSE
                WRITE(1,9026)I_WC,WC_LOC(IARM,I_WC),
     #                WC_C_ANGLE(IARM,I_WC),WC_W_ANGLE(IARM,I_WC),
     #                WC_MSCT(IARM,I_WC)
              ENDIF
              LG = .FALSE.
            ENDDO
          ENDIF
      ENDDO
9025  FORMAT(12X,I2,8X,I2,6X,F6.2,6X,F6.2,7X,F6.2,7X,F6.2)
9026  FORMAT(22X,I2,6X,F6.2,6X,F6.2,7X,F6.2,7X,F6.2)
      WRITE(1,8022)
 8022 FORMAT(//,T17,' TRANSPORT HISTOGRAM LIMITS, ETC. ',//,
     #    '  ID',T6,' Reaction',T19,' X or Y ',T28,' Minimum ',T37,
     #    ' Maximum ',T46,' # channels ',T58,' TRANSPORT Cut Index ',/,
     #    '  --',T6,' --------',T19,' ------ ',T28,' ------- ',T37,
     #    ' ------- ',T46,' ---------- ',T58,' ------------------- ')
      I_TR = 0
      DO IARM=1,2
        IF(NEL(IARM) .GT. 0) THEN
          DO IEL=1,NEL(IARM)
            IF(OP(IARM,IEL).EQ.'H1D') THEN
              I_TR = I_TR + 1
              IF(NCUTS_TR(IARM,IEL) .GT. 0) THEN
                WRITE(1,8023)I_TR,REACT_NAME(IPID_TR(IARM,IEL)),
     #          TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL),NXCHAN_TR(IARM,IEL),
     #           (ICUT_IND_TR(IARM,IEL,J),J=1,NCUTS_TR(IARM,IEL))
              ELSE
                WRITE(1,8023)I_TR,REACT_NAME(IPID_TR(IARM,IEL)),
     #           TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL),NXCHAN_TR(IARM,IEL)
              ENDIF
 8023         FORMAT(1X,I3,T6,A8,T22,' X',T29,F7.2,T38,F7.2,T49,I5,
     #             T58,20(1X,I3))
            ELSEIF(OP(IARM,IEL).EQ.'H2D') THEN
              I_TR = I_TR + 1
              IF(NCUTS_TR(IARM,IEL) .GT. 0) THEN
                WRITE(1,8023)I_TR,REACT_NAME(IPID_TR(IARM,IEL)),
     #          TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL),NXCHAN_TR(IARM,IEL),
     #           (ICUT_IND_TR(IARM,IEL,J),J=1,NCUTS_TR(IARM,IEL))
              ELSE
                WRITE(1,8023)I_TR,REACT_NAME(IPID_TR(IARM,IEL)),
     #           TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL),NXCHAN_TR(IARM,IEL)
              ENDIF
              WRITE(1,8024)TR_MINY(IARM,IEL),TR_MAXY(IARM,IEL),
     #          NYCHAN_TR(IARM,IEL)
 8024         FORMAT(T22,' Y',T29,F7.2,T38,F7.2,T49,I5)
            ELSEIF(OP(IARM,IEL).EQ.'SCT') THEN
              I_TR = I_TR + 1
              IF(NCUTS_TR(IARM,IEL) .GT. 0) THEN
                WRITE(1,8025)I_TR,TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL),
     #           (ICUT_IND_TR(IARM,IEL,J),J=1,NCUTS_TR(IARM,IEL))
              ELSE
                WRITE(1,8025)I_TR,TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL)
              ENDIF
 8025         FORMAT(1X,I3,T22,' X',T29,F7.2,T38,F7.2,T58,20(1X,I3))
              WRITE(1,8026)TR_MINY(IARM,IEL),TR_MAXY(IARM,IEL)
 8026         FORMAT(T22,' Y',T29,F7.2,T38,F7.2)
            ELSEIF(OP(IARM,IEL).EQ.'NTU') THEN
              I_TR = I_TR + 1
              WRITE(1,8027)I_TR,REACT_NAME(IPID_TR(IARM,IEL))
 8027         FORMAT(1X,I3,T6,A8)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C ---------------------------------------------------------------------
C       Write CUTS information.
C ---------------------------------------------------------------------
C
      WRITE(1,8040)
 8040 FORMAT(T32,' GLOBAL CUTS ',//,
     #          T15,' Variable',T37,' Xmin',T52,' Xmax',/,
     #          T15,' --------',T37,' ----',T52,' ----')
      IF(NCUTS .GT. 0) THEN
        DO I=1,NCUTS
          WRITE(1,8041)VAR_NAME(ICUT_VAR(I)),CUT_MIN(I),CUT_MAX(I)
 8041     FORMAT(T16,A8,T33,E12.4,T48,E12.4)
        ENDDO
      ENDIF
C
      WRITE(1,8042)
 8042 FORMAT(//,T25,' HISTOGRAM SPECIFIC CUTS ',//,
     #    T10,' Cut Index',T25,' Variable',T45,' Min',T60,' Max',/,
     #    T10,' ---------',T25,' --------',T45,' ---',T60,' ---')
      IF(NCUTS_I .GT. 0) THEN
        DO I=1,NCUTS_I
          WRITE(1,8043)ICUT_IND(I),VAR_NAME(ICUT_VAR_I(ICUT_IND(I))),
     #          CUT_MIN_I(ICUT_IND(I)), CUT_MAX_I(ICUT_IND(I))
 8043     FORMAT(T13,I3,T26,A8,T42,E11.4,T57,E11.4)
        ENDDO
      ENDIF
C
C ---------------------------------------------------------------------
C       Write kinematics histograms input.
C ---------------------------------------------------------------------
C
      WRITE(1,8050)
 8050 FORMAT(//,T10,
     #      ' K I N E M A T I C S      H I S T O G R A M S ',/,
     #  T10,' -------------------      ------------------- ')
      WRITE(1,8051)
 8051 FORMAT(//'  ID',T6,' Type',T13,' Reaction   ',T27,' X Variable',
     #     T39,' Y Variable',T51,' File name',T67,' Cuts',/,
     #         '  --',T6,' ----',T13,' --------   ',T27,' ----------',
     #     T39,' ----------',T51,' ---------',T67,' ----')
C
      IF(NPLOTS .GT. 0) THEN
        DO I=1,NPLOTS
          CALL SQUEEZE(PLOT_FIL(I),PL_F,NCHAR)
C
          IF(PLOT_TYPE(I) .EQ. 'P1D') THEN
            IF(NCUTS_PL(I) .GT. 0) THEN
              WRITE(1,8052)I,PLOT_TYPE(I),REACT_NAME(IPID(I)),
     #          VAR_NAME(I_VAR(I,1)),
     #          PL_F(1:NCHAR),(ICUT_IND_PL(I,J),J=1,NCUTS_PL(I))
            ELSE
              WRITE(1,8052)I,PLOT_TYPE(I),REACT_NAME(IPID(I)),
     #          VAR_NAME(I_VAR(I,1)),
     #          PL_F(1:NCHAR)
            ENDIF
 8052       FORMAT(1X,I3,T8,A3,T14,A8,T30,A8,T52,A,T66,
     #          30(1X,I3))
C
          ELSEIF(PLOT_TYPE(I) .EQ. 'P2D') THEN
            IF(NCUTS_PL(I) .GT. 0) THEN
              WRITE(1,8053)I,PLOT_TYPE(I),REACT_NAME(IPID(I)),
     #          VAR_NAME(I_VAR(I,1)),VAR_NAME(I_VAR(I,2)),
     #          PL_F(1:NCHAR),(ICUT_IND_PL(I,J),J=1,NCUTS_PL(I))
            ELSE
              WRITE(1,8053)I,PLOT_TYPE(I),REACT_NAME(IPID(I)),
     #          VAR_NAME(I_VAR(I,1)),VAR_NAME(I_VAR(I,2)),
     #          PL_F(1:NCHAR)
            ENDIF
 8053       FORMAT(1X,I3,T8,A3,T14,A8,T30,A8,T41,A8,T52,A,T66,
     #          30(1X,I3))
C
          ELSEIF(PLOT_TYPE(I) .EQ. 'SCA') THEN
            IF(NCUTS_PL(I) .GT. 0) THEN
              WRITE(1,8054)I,PLOT_TYPE(I),
     #          VAR_NAME(I_VAR(I,1)),VAR_NAME(I_VAR(I,2)),
     #          PL_F(1:NCHAR),(ICUT_IND_PL(I,J),J=1,NCUTS_PL(I))
            ELSE
              WRITE(1,8054)I,PLOT_TYPE(I),
     #          VAR_NAME(I_VAR(I,1)),VAR_NAME(I_VAR(I,2)),
     #          PL_F(1:NCHAR)
            ENDIF
 8054       FORMAT(1X,I3,T8,A3,T30,A8,T41,A8,T52,A,T66,
     #          30(1X,I3))
C
          ELSEIF(PLOT_TYPE(I) .EQ. 'NTU') THEN
            WRITE(1,8099)I,PLOT_TYPE(I),REACT_NAME(IPID(I)),
     #          PL_F(1:NCHAR)
 8099       FORMAT(1X,I3,T8,A3,T14,A8,T52,A)
C
          ELSEIF(PLOT_TYPE(I) .EQ. 'NTM') THEN
            WRITE(1,8200)I,PLOT_TYPE(I),PL_F(1:NCHAR)
 8200       FORMAT(1X,I3,T8,A3,T52,A)
C
          ENDIF
        ENDDO
      ENDIF
C
      WRITE(1,8055)
 8055 FORMAT(//,T21,' HISTOGRAM LIMITS, ETC. ',//,
     #    '  ID',T6,' X or Y ',T17,' Minimum ',T28,' Maximum ',T37,
     #    ' # channels ',T51,' Scale ',T61,' Offset ',/,
     #    '  --',T6,' ------ ',T17,' ------- ',T28,' ------- ',T37,
     #    ' ---------- ',T51,' ----- ',T61,' ------ ')
C
      IF(NPLOTS .GT. 0) THEN
        DO I=1,NPLOTS
          IF(PLOT_TYPE(I) .EQ. 'P1D') THEN
            WRITE(1,8056)I,X_MIN(I),X_MAX(I),NX_CHAN(I),X_SC(I),X_OFF(I)
 8056       FORMAT(1X,I3,T9,' X',T16,E10.3,T27,E10.3,T40,I5,T49,E10.3,
     #           T60,E10.3)
C
          ELSEIF(PLOT_TYPE(I) .EQ. 'P2D') THEN
            WRITE(1,8056)I,X_MIN(I),X_MAX(I),NX_CHAN(I),X_SC(I),X_OFF(I)
            WRITE(1,8057)Y_MIN(I),Y_MAX(I),NY_CHAN(I),Y_SC(I),Y_OFF(I)
 8057       FORMAT(T9,' Y',T16,E10.3,T27,E10.3,T40,I5,T49,E10.3,
     #           T60,E10.3)
C
          ELSEIF(PLOT_TYPE(I) .EQ. 'SCA') THEN
            WRITE(1,8058)I,X_MIN(I),X_MAX(I),X_SC(I),X_OFF(I)
 8058       FORMAT(1X,I3,T9,' X',T16,E10.3,T27,E10.3,T49,E10.3,
     #           T60,E10.3)
            WRITE(1,8059)Y_MIN(I),Y_MAX(I),Y_SC(I),Y_OFF(I)
 8059       FORMAT(T9,' Y',T16,E10.3,T27,E10.3,T49,E10.3,T60,E10.3)
C
          ENDIF
        ENDDO
      ENDIF
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C       WRITE SUMMARY FILE (OUTPUT SECTION).
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
      WRITE(1,9001)
 9001 FORMAT(//,72('-'),/,3x,66('*'),/,72('-'),
     #       //,T22,' O U T P U T     S U M M A R Y ',//)
C
C ---------------------------------------------------------------------
C       Get derived kinematical quantities for central kinematics.
C       Note that the central kinematics are computed from the
C       hadron arm central momentum (and not the missing mass)
C       whether or not the bound state case is selected.
C ---------------------------------------------------------------------
C
      CALL V3MAKE(E0,PH_B,TH_B,KI_VEC)    !make beam 3-vector
      CALL V3MAKE(PF_E,PH_E,TH_E,KF_VEC)  !make scatt e- 3-vector
C
      CALL KINEM(E0,KI_VEC,PF_E,KF_VEC,PF_P,PH_P,TH_P,
     #                        .FALSE.,.FALSE.,.FALSE.,FAIL_KIN)
C
      WRITE(1,9002)
     #          TSCAT*180./PI,
     #             VAR_AVG(14)*180./PI,VAR_CENT(14)*180./PI,
     #             VAR_MIN(14)*180./PI,VAR_MAX(14)*180./PI,
     #          OMEGA,
     #             VAR_AVG(15),VAR_CENT(15),VAR_MIN(15),VAR_MAX(15),
     #          QMAG,
     #             VAR_AVG(16),VAR_CENT(16),VAR_MIN(16),VAR_MAX(16),
     #          QMU2_G,
     #             VAR_AVG(19),VAR_CENT(19),
     #             VAR_MIN(19),VAR_MAX(19),
     #          PHIQ*180./PI,
     #             VAR_AVG(18)*180./PI,VAR_CENT(18)*180./PI,
     #             VAR_MIN(18)*180./PI,VAR_MAX(18)*180./PI,
     #          THEQ*180./PI,
     #             VAR_AVG(17)*180./PI,VAR_CENT(17)*180./PI,
     #             VAR_MIN(17)*180./PI,VAR_MAX(17)*180./PI,
     #          X,
     #             VAR_AVG(20),VAR_CENT(20),VAR_MIN(20),VAR_MAX(20),
     #          EPSILON,
     #             VAR_AVG(21),VAR_CENT(21),VAR_MIN(21),VAR_MAX(21),
     #          PF_P,
     #             VAR_AVG(12),VAR_CENT(12),VAR_MIN(12),VAR_MAX(12),
     #          PRMAG,
     #             VAR_AVG(26),VAR_CENT(26),VAR_MIN(26),VAR_MAX(26),
     #          THETA_PQ*180./PI,
     #             VAR_AVG(22)*180./PI,VAR_CENT(22)*180./PI,
     #             VAR_MIN(22)*180./PI,VAR_MAX(22)*180./PI,
     #          PHI_X*180./PI,
     #             VAR_AVG(23)*180./PI,VAR_CENT(23)*180./PI,
     #             VAR_MIN(23)*180./PI,VAR_MAX(23)*180./PI,
     #          MISS_M,
     #             VAR_AVG(24),VAR_CENT(24),VAR_MIN(24),VAR_MAX(24),
     #          W*0.001,
     #             VAR_AVG(28)*0.001,VAR_CENT(28)*0.001,
     #             VAR_MIN(28)*0.001,VAR_MAX(28)*0.001,
     #          THETA_CM*180./PI,
     #             VAR_AVG(29)*180./PI,VAR_CENT(29)*180./PI,
     #             VAR_MIN(29)*180./PI,VAR_MAX(29)*180./PI
C
 9002 FORMAT(/T15,' Kinematical quantities (after global cuts)',//,
     #          '   Variable     ',T21,' Central',2X,' Average',1X,
     #            ' Centroid',2X,' Minimum',2X,' Maximum',1X,'Units',/,
     #          '   --------     ',T21,' -------',2X,' -------',1X,
     #            ' --------',2X,' -------',2X,' -------',1X,'-----',/,
     #          ' Scatt. angle ',T20,5(F9.3,1X),' deg',/,
     #          ' Omega        ',T20,5(F9.3,1X),' MeV',/,
     #          ' Q-vector     ',T20,5(F9.3,1X),' MeV/c',/,
     #          ' Q_mu^2       ',T20,5(F9.3,1X),' GeV/c^2',/,
     #          ' Phi_Q        ',T20,5(F9.3,1X),' deg',/,
     #          ' Theta_Q      ',T20,5(F9.3,1X),' deg',/,
     #          ' X            ',T20,5(F9.3,1X),/,
     #          ' Long. Pol.   ',T20,5(F9.3,1X),//,
     #          ' P_f          ',T20,5(F9.3,1X),' MeV/c',/,
     #          ' P_recoil     ',T20,5(F9.3,1X),' MeV/c',/,
     #          ' Theta_pq     ',T20,5(F9.3,1X),' deg',/,
     #          ' Phi_x        ',T20,5(F9.3,1X),' deg',/,
     #          ' M_miss       ',T20,5(F9.3,1X),' MeV',/,
     #          ' W            ',T20,5(F9.3,1X),' GeV',/,
     #          ' Theta_CM     ',T20,5(F9.3,1X),' deg')
C
C ---------------------------------------------------------------------
C        Write summary of statistics (average cross sections, etc.).
C ---------------------------------------------------------------------
C
      WRITE(1,9003)
 9003 FORMAT(//,T30,' STATISTICS ',/,
     #          T30,' ---------- ',/,
     #          T15,' NOTE: if energy loss calc. performed ',/,
     #          T15,'       then acceptances initially enlarged. ',//)
C
      WRITE(1,9004)N_EVENT*NCALC,
     #          NINT(SUM_EV_C)
C
 9004 FORMAT(/' # tries: ',T50,I8,/,
     #        ' # passed global cuts + HRS apertures: ',T50,I8)
C
      WRITE(1,9014) SUM_EEP_C*SIGFACTOR(1)/(SUM_ACC*BEAMTIME),
     #              SUM_EEP_C/SUM_PS_C,CHAR_SIG_UNITS,
     #              I_EXTRAP,I_INTERPEXP_NP,I_SPECT_OOR
C
 9014 FORMAT(/' (e,e''N) rate after cuts: ',T50,E13.5,' sec^-1',/,
     #        ' sigma after global cuts: ',T40,E13.5,1X,A21,//,
     #        ' Bad interpolation (extrapolation): ',T50,I8,/,
     #        ' Expon. interpolation failures: ',T50,I8,/,
     #        ' Spectral function out-of-range: ',T50,I8,//)
C
      IF(IPHYS_OPT .EQ. 804) WRITE(1,9015) I_AREN_FAIL
 9015 FORMAT(/' Failed Arenhoevel interpolation: ',T50,I8,//)
      IF(IPHYS_OPT .EQ. 814) WRITE(1,9016) LAGET_PS_FAIL,LAGET_GRID_FAIL
 9016 FORMAT(/' Failed Laget interpolation: (before cuts) ',
     #      /,'        Phase Space:      ',I8,
     #      /,'        Grid  Boundary:   ',I8,//)
 
C
C ---------------------------------------------------------------------
C        Write TRANSPORT histogram summary.
C ---------------------------------------------------------------------
C
      WRITE(1,9005)
 9005 FORMAT(T14,' TRANSPORT HISTOGRAM SUMMARY ',//,
     # '  ID',T10,'  OOB ',T22,' Sum',T36,' X-cent. ',T51,' Y-cent. ',/,
     # '  --',T10,'  --- ',T22,' ---',T36,' ------- ',T51,' ------- ')
      I_TR = 0
      DO IARM=1,2
        IF(NEL(IARM) .GT. 0) THEN
          DO IEL=1,NEL(IARM)
            IF(OP(IARM,IEL).EQ.'H1D') THEN
              I_TR = I_TR + 1
              DO J=1,NXCHAN_TR(IARM,IEL)
                ARR(J) = TRHIST(IARM,IEL,J)
              ENDDO
              CALL SPECSTATS_1D(TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL),
     #                  NXCHAN_TR(IARM,IEL),ARR,CENT_X,SUM)
              WRITE(1,9006)I_TR,I_OOB_TR(IARM,IEL),SUM,CENT_X
 9006         FORMAT(1X,I3,T8,I7,T20,E10.3,T35,E10.3)
            ELSEIF(OP(IARM,IEL).EQ.'H2D') THEN
              I_TR = I_TR + 1
              DO J=1,NXCHAN_TR(IARM,IEL)
              DO K=1,NYCHAN_TR(IARM,IEL)
                ARR_2D(J,K) = TRHIST_2D(IARM,IEL,J,K)
              ENDDO
              ENDDO
              CALL SPECSTATS_2D(TR_MINX(IARM,IEL),TR_MAXX(IARM,IEL),
     #                TR_MINY(IARM,IEL),TR_MAXY(IARM,IEL),
     #                NXCHAN_TR(IARM,IEL),NYCHAN_TR(IARM,IEL),
     #                ARR_2D,CENT_X,CENT_Y,SUM)
              WRITE(1,9007)I_TR,I_OOB_TR(IARM,IEL),SUM,CENT_X,CENT_Y
 9007         FORMAT(1X,I3,T8,I7,T20,E10.3,T35,E10.3,T50,E10.3)
            ELSEIF(OP(IARM,IEL).EQ.'SCT') THEN
              I_TR = I_TR + 1
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C ---------------------------------------------------------------------
C        Write Kinematics histogram summary.
C ---------------------------------------------------------------------
C
      WRITE(1,9008)
 9008 FORMAT(//,T14,' KINEMATICS HISTOGRAM SUMMARY ',//,
     # '  ID',T10,'  OOB ',T22,' Sum',T36,' X-cent. ',T51,' Y-cent. ',/,
     # '  --',T10,'  --- ',T22,' ---',T36,' ------- ',T51,' ------- ')
      IF(NPLOTS .GT. 0) THEN
        DO I=1,NPLOTS
          IF(PLOT_TYPE(I) .EQ. 'P1D') THEN
            DO J=1,NX_CHAN(I)
              ARR(J) = SIG_ARR(I,J)
            ENDDO
            CALL SPECSTATS_1D(X_MIN(I),X_MAX(I),NX_CHAN(I),ARR,CENT_X,
     #           SUM)
            WRITE(1,9006)I,I_OOB(I),SUM,CENT_X
          ELSEIF(PLOT_TYPE(I) .EQ. 'P2D') THEN
            DO J=1,NX_CHAN(I)
            DO K=1,NY_CHAN(I)
              ARR_2D(J,K) = SIG_ARR_2D(I,J,K)
            ENDDO
            ENDDO
            CALL SPECSTATS_2D(X_MIN(I),X_MAX(I),Y_MIN(I),Y_MAX(I),
     #           NX_CHAN(I),NY_CHAN(I),ARR_2D,CENT_X,CENT_Y,SUM)
            WRITE(1,9007)I,I_OOB(I),SUM,CENT_X,CENT_Y
          ENDIF
        ENDDO
      ENDIF
C
      WRITE(1,9080)
 9080 FORMAT(//)
C
C ---------------------------------------------------------------------
C        Write Energy loss information, if included.
C ---------------------------------------------------------------------
C
      IF(ELOSS_MOD .NE. 0) THEN
       WRITE(1,9090) ELWALK(1),ELWALK(2),ELWALK(4),ELWALK(3),ELWALK(5)
      ENDIF
 9090 FORMAT(T10,' ENERGY LOSS FROM OTHER THAN TARGET MATERIAL ',/,
     #       T10,' ------------------------------------------- ',/,
     #       T10,' Beam:                       ',T40,F10.4,' MeV',/,
     #       T10,' Scatt. Electron (sidewall): ',T40,F10.4,' MeV'/,
     #       T10,' Hadron          (sidewall): ',T40,F10.4,' MeV'/,
     #       T10,' Scatt. Electron   (endcap): ',T40,F10.4,' MeV'/,
     #       T10,' Hadron            (endcap): ',T40,F10.4,' MeV'/,
     #       //)
C
      CLOSE(UNIT=1)
C
      RETURN
      END
