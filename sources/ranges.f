C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C       Subroutine RANGES:
C
C               Calculates default ranges for plotting variables
C               given kinematics and spectrometer acceptances.
C
C               For elastic scattering, A(e,e'A), the histogram
C               ranges are calculated as though the continuum
C               (e,e'X) case were selected.  This will overestimate
C               the limits since elastic scattering is more
C               kinematically restricted than continuum.
C
C ------------------------------------------------------------------
C
      SUBROUTINE RANGES(ELASTIC,BOUND)
C
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
C
      COMMON /KINVAR/ KI,Q,Q2,QMU2,EP,PF,THETA_EP,CTH_PQ,PR
      COMMON /ELIMS/ EF_MIN,EF_MAX
      COMMON /PLIMS/ PF_MIN,PF_MAX
      COMMON /PFROOTS_D/ PF_ROOT_WT
      COMMON /PFROOTS_L/ FIRST_TIME
C
      DOUBLE PRECISION PL_MINX(50),PL_MAXX(50),PL_MINY(50),PL_MAXY(50)
      DOUBLE PRECISION KI(3),Q(3),PF(3),PR(3),KI_VEC(3),KF_VEC(3)
C
      LOGICAL FAIL_KIN,ELASTIC,BOUND,FAIL_RANGES,FIRST_TIME
C
C ---------------------------------------------------------------------
C       Get input information and arrays for plots and cuts.
C ---------------------------------------------------------------------
C
      INCLUDE 'input.cmn'
      INCLUDE 'var.cmn'
C
      PI = ASIN(1.D0)*2.D0
C
C ---------------------------------------------------------------------
C       Initialize counter which determines whether limits
C       could be found (if no or only one physical solution was found
C       then limits are meaningless - this can happen especially for
C       a bound state where all "solutions" give rise to hadron
C       momenta outside the acceptance).
C ---------------------------------------------------------------------
C
      NPASS = 0
C
C ---------------------------------------------------------------------
C       Limits for angular variables.
C ---------------------------------------------------------------------
C
      THEP_MIN = -DTH_P/2.            !Assume symmetric.
      THEE_MIN = -DTH_E/2.
      PHIP_MIN = -DPH_P/2.
      PHIE_MIN = -DPH_E/2.
C
C ---------------------------------------------------------------------
C       Calculate bounds on electron and hadron momenta.
C ---------------------------------------------------------------------
C
      EF_MIN = PF_E*(1.+0.01*ACC_EM)
      EF_MAX = PF_E*(1.+0.01*ACC_EP)
      DEF = EF_MAX - EF_MIN
C
      PF_MIN = PF_P*(1.+0.01*ACC_PM)
      PF_MAX = PF_P*(1.+0.01*ACC_PP)
      DPF = PF_MAX - PF_MIN
C
C ---------------------------------------------------------------------
C       Initialize limits and set default scale factors.
C       (Default offsets = 0)
C ---------------------------------------------------------------------
C
      DO I=1,NPLOTS
        IF(PLOT_TYPE(I).NE.'NTU'.AND.PLOT_TYPE(I).NE.'NTM') THEN
           PL_MINX(I) =  1.D16
           PL_MAXX(I) = -1.D16
           IF(X_SC(I)  .EQ. 0.) X_SC(I)  = 1.D0
           IF(TWO_D(I)) THEN
             PL_MINY(I) =  1.D16
             PL_MAXY(I) = -1.D16
             IF(Y_SC(I)  .EQ. 0.) Y_SC(I)  = 1.D0
           ENDIF
        ENDIF
      ENDDO
C
C ---------------------------------------------------------------------
C       Start the loops to determine ranges of variables.
C ---------------------------------------------------------------------
C
      DO J1=0,NITER1                  !Theta_p
        DO J2=0,NITER2                  !Phi_p
          DO J3=0,NITER3                  !Theta_e
            DO J4=0,NITER4                  !Phi_e
              DO J5=0,NITER5                  !electron momentum
                DO J6=0,NITER6                  !hadron momentum
C
                  FIRST_TIME = .TRUE.    !TRUE for first call to KINEM
C
                  TH_P_I = DTH_P*DFLOAT(J1)/DFLOAT(NITER1) + THEP_MIN +
     #              TH_P
                  PH_P_I = (DPH_P*DFLOAT(J2)/DFLOAT(NITER2)
     #               + PHIP_MIN)/COS(TH_P_I) + PH_P
                  TH_E_I = DTH_E*DFLOAT(J3)/DFLOAT(NITER3) + THEE_MIN +
     #              TH_E
                  PH_E_I = (DPH_E*DFLOAT(J4)/DFLOAT(NITER4)
     #               + PHIE_MIN)/COS(TH_E_I) + PH_E
                  PF_E_I = DEF*DFLOAT(J5)/DFLOAT(NITER5) + EF_MIN
                  PF_P_I = DPF*DFLOAT(J6)/DFLOAT(NITER6) + PF_MIN
C
C ---------------------------------------------------------------------
C          The COS(TH_P_I/TH_E_I) above accounts for difference between
C          collimator angular range and the initial beam "transport"
C          system PHI range.
C
C
C          Determine derived kinematical quantities.
C ---------------------------------------------------------------------
C
                  CALL V3MAKE(E0,PH_B,TH_B,KI_VEC)    !make beam 3-vector
                  CALL V3MAKE(PF_E_I,PH_E_I,TH_E_I,KF_VEC)  !make scatt
                                                             !e- 3-vector
C
                  CALL KINEM(E0,KI_VEC,PF_E_I,KF_VEC,PF_P_I,PH_P_I,
     #                        TH_P_I,ELASTIC,BOUND,.FALSE.,FAIL_KIN)
                  IF(FAIL_KIN) GOTO 999
                  NPASS = NPASS + 1            !Found a physical solution
C
C ---------------------------------------------------------------------
C          Determine minimum and maximum values for the variables.
C ---------------------------------------------------------------------
C
                  DO I=1,NPLOTS
                    IF(PLOT_TYPE(I).NE.'NTU'
     #                    .AND.PLOT_TYPE(I).NE.'NTM') THEN
                       IF(VAR(I_VAR(I,1)) .LT. PL_MINX(I))
     #                     PL_MINX(I) = VAR(I_VAR(I,1))
                       IF(VAR(I_VAR(I,1)) .GT. PL_MAXX(I))
     #                     PL_MAXX(I) = VAR(I_VAR(I,1))
                       IF(TWO_D(I)) THEN
                           IF(VAR(I_VAR(I,2)) .LT. PL_MINY(I))
     #                        PL_MINY(I) = VAR(I_VAR(I,2))
                           IF(VAR(I_VAR(I,2)) .GT. PL_MAXY(I))
     #                        PL_MAXY(I) = VAR(I_VAR(I,2))
                       ENDIF
                    ENDIF
                  ENDDO
C
C ---------------------------------------------------------------------
C
  999           ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C ---------------------------------------------------------------------
C       Are limits meaningful?
C ---------------------------------------------------------------------
C
      IF((NPASS .EQ. 0) .OR. (NPASS .EQ. 1)) THEN
        FAIL_RANGES = .TRUE.
      ELSE
        FAIL_RANGES = .FALSE.
      ENDIF
C
C ---------------------------------------------------------------------
C       Widen limits to be safe.
C ---------------------------------------------------------------------
C
      IF(.NOT. FAIL_RANGES) THEN
        DO I=1,NPLOTS
          IF(PLOT_TYPE(I).NE.'NTU'.AND.PLOT_TYPE(I).NE.'NTM') THEN
             DELTAX = PL_MAXX(I) - PL_MINX(I)
             PL_MAXX(I) = PL_MAXX(I) + DELTAX*0.1
             PL_MINX(I) = PL_MINX(I) - DELTAX*0.1
             IF(TWO_D(I)) THEN
               DELTAY = PL_MAXY(I) - PL_MINY(I)
               PL_MAXY(I) = PL_MAXY(I) + DELTAY*0.1
               PL_MINY(I) = PL_MINY(I) - DELTAY*0.1
             ENDIF
          ENDIF
        ENDDO
      ENDIF
C
      IFAIL = 0       !initialize
C
      DO I=1,NPLOTS
        IF(PLOT_TYPE(I).NE.'NTU'.AND.PLOT_TYPE(I).NE.'NTM') THEN
           IF(X_MIN(I) .EQ. 0.) THEN    !limit not given in input file
             IF(FAIL_RANGES) THEN
               IFAIL = IFAIL + 1
               IF(IFAIL .EQ. 1) WRITE(6,1259)
               WRITE(6,1250)I,VAR_NAME(I_VAR(I,1))
               READ(5,*)X_MIN(I)
             ELSE
               X_MIN(I) = PL_MINX(I)*X_SC(I) + X_OFF(I)
             ENDIF
           ENDIF
           IF(X_MAX(I) .EQ. 0.) THEN    !limit not given in input file
             IF(FAIL_RANGES) THEN
               IFAIL = IFAIL + 1
               IF(IFAIL .EQ. 1) WRITE(6,1259)
               WRITE(6,1251)I,VAR_NAME(I_VAR(I,1))
               READ(5,*)X_MAX(I)
             ELSE
               X_MAX(I) = PL_MAXX(I)*X_SC(I) + X_OFF(I)
             ENDIF
           ENDIF
C
           IF(TWO_D(I)) THEN
             IF(Y_MIN(I) .EQ. 0.) THEN    !limit not given in input file
               IF(FAIL_RANGES) THEN
                 IFAIL = IFAIL + 1
                 IF(IFAIL .EQ. 1) WRITE(6,1259)
                 WRITE(6,1252)I,VAR_NAME(I_VAR(I,2))
                 READ(5,*)Y_MIN(I)
               ELSE
                 Y_MIN(I) = PL_MINY(I)*Y_SC(I) + Y_OFF(I)
               ENDIF
             ENDIF
             IF(Y_MAX(I) .EQ. 0.) THEN    !limit not given in input file
               IF(FAIL_RANGES) THEN
                 IFAIL = IFAIL + 1
                 IF(IFAIL .EQ. 1) WRITE(6,1259)
                 WRITE(6,1253)I,VAR_NAME(I_VAR(I,2))
                 READ(5,*)Y_MAX(I)
               ELSE
                 Y_MAX(I) = PL_MAXY(I)*Y_SC(I) + Y_OFF(I)
               ENDIF
             ENDIF
           ENDIF
        ENDIF
      ENDDO
C
 1259 FORMAT(/' Histogram range routine failed ',/,
     #          ' -----------------------------')
 1250 FORMAT(' Histogram #',I3,3X,' Variable ',A8,4X,' Enter Xmin >')
 1251 FORMAT(' Histogram #',I3,3X,' Variable ',A8,4X,' Enter Xmax >')
 1252 FORMAT(' Histogram #',I3,3X,' Variable ',A8,4X,' Enter Ymin >')
 1253 FORMAT(' Histogram #',I3,3X,' Variable ',A8,4X,' Enter Ymax >')
C
      RETURN
      END
