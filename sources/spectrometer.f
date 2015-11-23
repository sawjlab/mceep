C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C       Subroutine SPECTROMETER:
C
C       AUTHOR: P.E. Ulmer
C       DATE:   V1.0    8-OCT-1990
C       PURPOSE:
C                   Performs spectrometer simulation:
C
C                   From Transport vectors at the target this routine
C                   gives resulting Transport vectors after all
C                   "spectrometer" elements have been incorporated.
C                   The input routine, INPUT, lists options of this
C                   routine.  Briefly, this routine allows for
C                   matrix multiplication (to simulate spectrometer 
C                   map), multiplication by a spin rotation matrix
C                   (currently, an approximate solution is incorporated
C                   which assumes the spin precession angle is
C                   proportional to the net bend angle for a given
C                   element; the user can include more sophistocated
C                   prescriptions by modifying the user subroutines
C                   in ../sources/spin_dat.for),
C                   offsets of a particular Transport ray coordinate
C                   (to simulate misalignments, etc.) and folding of
C                   a particular ray coordinate with a gaussian 
C                   resolution function (to simulate multiple 
C                   scattering, etc.).  This routine provides 1-D and 
C                   2-D histograms as well as scatter plots and cuts on
C                   any ray coordinate which can either be attached to a
C                   specific Transport histogram or all histograms
C                   (global cut).  Global cuts affect the ordinary
C                   kinematics histograms as well.
C
C                   NOTE:  ICALL specifies whether the cross section 
C                          is being evaluated (ICALL=1) or the 
C                          kinematic variables(ICALL=2).
C
C ------------------------------------------------------------------
C
      SUBROUTINE SPECTROMETER(ICALL)
      IMPLICIT NONE
C
      DOUBLE PRECISION TVEC_TMP(6),TVEC_TMP2(6),GASDEVV,TEST_ELL
      DOUBLE PRECISION PMAT_EL(3,3),PMAT(3,3),C_RAY(3),C_RAY_ROT(3)
      DOUBLE PRECISION C,BETA,TOFE,TOFP,PMOM
      DOUBLE PRECISION TOF_REL_TMP,ANGLE,TEST_RFN,E_RFN,P_RFN
      DOUBLE PRECISION COSY_VECT(5)
      DOUBLE PRECISION X_BEAM_E,X_BEAM_P,Y_BEAM_E,Y_BEAM_P
      DOUBLE PRECISION X_BEAM,Y_BEAM,FP_VEC_E(6),FP_VEC_P(6)
C
      REAL             RFUNCTION,Y_TMP,DEL_TMP,THE_TMP,PHI_TMP
C
      INTEGER          ICALL,IARM,IMAT,IJ,IEL,I,J,K,IMAD
      LOGICAL          TR_CUT_FLAG(60)

C ---------------------------------------------------------------------
C The arguments for Hall A HRS COSY calculation
C ---------------------------------------------------------------------

      DOUBLE PRECISION x_tg,y_tg,z_tg        !target values (cm)
      DOUBLE PRECISION dpp                   !delta p/p (%)
      DOUBLE PRECISION dxdz,dydz             !X,Y slope in spectrometer
      DOUBLE PRECISION x_fp,y_fp,dx_fp,dy_fp !Focal plane values to return
      DOUBLE PRECISION p_spec,fry            !spectrometer setting 
      DOUBLE PRECISION m2                    !particle mass
      LOGICAL*4        ms_flag,wcs_flag,col_flag !mult scatt, 
                                                 ! VDC smearing
                                                 ! and collimator flag
      LOGICAL*4	       ok_hrs                !true if particle makes it 
                                             ! through spectrometer
C
C ---------------------------------------------------------------------
C       Spectrometer analysis common area.
C ---------------------------------------------------------------------
C
      COMMON /XTGT/ X_BEAM_E,X_BEAM_P,Y_BEAM_E,Y_BEAM_P
      COMMON /FPVEC/ FP_VEC_E,FP_VEC_P
      COMMON /TOF/ TOF_REL_TMP,TOFE,TOFP
      COMMON /RFN/ E_RFN,P_RFN
C
      INCLUDE 'spectrometer.cmn'
      INCLUDE 'wc.cmn'
      INCLUDE 'masses.cmn'
      INCLUDE 'input.cmn'
C
C ---------------------------------------------------------------------
C       New COSY related commons
C ---------------------------------------------------------------------
C
      INCLUDE 'struct_hrs.cmn'
      INCLUDE 'apertures.cmn'
      INCLUDE 'track.cmn'
C
      PARAMETER (C=30.D0)       !speed of light in cm/nsec
      FAIL_GLOBAL_CUT = .FALSE. !Initialize flag
C
C ---------------------------------------------------------------------
C     Perform initializations.
C ---------------------------------------------------------------------
C
      DO I=1,3                   !Init. spin rotation matrix to identity
        DO J=1,3
          IF (I.EQ.J) THEN
             PMAT(I,J) = 1.D0
          ELSE
             PMAT(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
C ---------------------------------------------------------------------
C     Loop over each arm (electron, then proton (IARM = 1,2)) and
C     over individual elements (Matrix, CUT, Histogram, ...).
C
C     Determine whether a given element is to be applied in this
C     particular call to SPECTROMETER by checking the logical flag
C     ACTIVE_EL.
C ---------------------------------------------------------------------
C
      DO IARM=1,2
        C_RAY(1) = 0.D0   !Init. central ray vector's direction
        C_RAY(2) = 0.D0
        C_RAY(3) = 1.D0
C
        IMAT = 0        !Init. Transport matrix number for each arm
        NPOL_EL = 0     !Init. spin rotation element #
        IF(NEL(IARM) .GT. 0) THEN
          DO IEL=1,NEL(IARM)
            IF( ACTIVE_EL(IARM,IEL,ICALL)) THEN
C
C ---------------------------------------------------------------------
C            Multiply by Transport matrix.
C ---------------------------------------------------------------------
C
              IF(OP(IARM,IEL).EQ.'MAT') THEN
                IMAT = IMAT + 1
                DO I=1,6
                  TVEC_TMP(I) = 0.D0
                  DO J=1,6
                    TVEC_TMP(I) = TVEC_TMP(I)
     #                          + TVEC(J,IARM)*TMATRIX(IARM,IMAT,I,J)
                  ENDDO
                ENDDO
                IF (ORDER(IARM,IEL) .EQ. 2) THEN   !2nd order matrix
                  DO I=1,5
                    DO K=1,6
                      DO IJ=1,K
                        TVEC_TMP(I) = TVEC_TMP(I) + 
     #                                 TMATRIX2(IARM,IMAT,I,IJ,K)*
     #                                 TVEC(IJ,IARM)*TVEC(K,IARM)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDIF
                DO I=1,6
                    TVEC(I,IARM) = TVEC_TMP(I)
                ENDDO
C
C ----------------------------------------------------------------------
C            Get spin rotation matrix for this element and form a
C            running product.
C ----------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'POL')THEN
                IF(IARM .EQ. 1) THEN
                  WRITE(6,100)
100               FORMAT(/' Illegal: spin trans. for hadron arm only ')
                  STOP
                ENDIF
                NPOL_EL = NPOL_EL + 1
                CALL SPIN_PRECESS(NPOL_EL,BEND_ANGLE(IARM,IEL),PMAT_EL)
                CALL MAT_MULT(3,3,3, PMAT_EL, PMAT, PMATRIX)
                DO I = 1,3     !set PMAT to PMATRIX to allow running prod.
                  DO J = 1,3
                    PMAT(I,J) = PMATRIX(I,J)
                  ENDDO
                ENDDO
C
C ----------------------------------------------------------------------
C            Get spin rotation matrix from COSY.
C            Use for target polarizations to focal plane polarizations
C            in transport coordinates.  COSY matrix parametrized in
C            quantities at the target.
C ----------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'COS')THEN
                IF(IARM .EQ. 1) THEN
                  WRITE(6,100)
                  STOP
                ENDIF
                CALL COMPUTE_TRANSP_TO_COSY(TVEC,PF_P,COSY_VECT)
                CALL COMPUTE_COSY_SPIN(COSY_VECT,PMATRIX)
C
C ----------------------------------------------------------------------
C            Take into account the effect of a drift element.
C ----------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'DFT') THEN
                DO I=1,6
                  TVEC_TMP(I) = TVEC(I,IARM)
                ENDDO
                CALL DRIFT(TVEC_TMP,TVEC_TMP2,DRIFT_LENGTH(IARM,IEL))
                DO I=1,6
                  TVEC(I,IARM) = TVEC_TMP2(I)
                ENDDO
C
C ---------------------------------------------------------------------
C            Perform rotation about a given axis by the given angle.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'ROT') THEN
                DO I=1,6
                  TVEC_TMP(I) = TVEC(I,IARM)
                ENDDO
                CALL ROT_TRPT(TVEC_TMP,TVEC_TMP2,ROT_AXIS(IARM,IEL),
     #                C_ROT_ANGLE(IARM,IEL),S_ROT_ANGLE(IARM,IEL),
     #                C_RAY,C_RAY_ROT,.TRUE.)
                DO I=1,6
                  TVEC(I,IARM) = TVEC_TMP2(I)
                ENDDO
                DO I=1,3
                  C_RAY(I) = C_RAY_ROT(I)
                ENDDO
C
C ---------------------------------------------------------------------
C            Get JLAB Hall A HRS focal plane vector from target vector
C            and check that all apertures were passed.
C            Also, store focal plane vector for histogramming
C            (Version 3.4 and later).
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'HRS') THEN
                DO I=1,6
                  TVEC_TMP(I) = TVEC(I,IARM)
                ENDDO
                CALL HRS(HRS_ID_INT(IARM,IEL),TVEC_TMP,TVEC_TMP2,
     #                FAIL_GLOBAL_CUT)
                DO I=1,6
                  TVEC(I,IARM) = TVEC_TMP2(I)
                ENDDO
                IF(IARM .EQ. 1) THEN
                   DO I=1,6
                     FP_VEC_E(I) = TVEC(I,IARM)
                   ENDDO
                ELSE
                   DO I=1,6
                     FP_VEC_P(I) = TVEC(I,IARM)
                   ENDDO
                ENDIF
C ---------------------------------------------------------------------
C            Get JLAB Hall A HRS+septum focal plane vector from target vector
C            and check that all apertures were passed.
C            Also, store focal plane vector for histogramming
C            (Version 3.4 and later).
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'SEP') THEN
                DO I=1,6
                  TVEC_TMP(I) = TVEC(I,IARM)
                ENDDO
                CALL SEPTUM(HRS_ID_INT(IARM,IEL),TVEC_TMP,TVEC_TMP2,
     #                FAIL_GLOBAL_CUT)
                DO I=1,6
                  TVEC(I,IARM) = TVEC_TMP2(I)
                ENDDO
                IF(IARM .EQ. 1) THEN
                   DO I=1,6
                     FP_VEC_E(I) = TVEC(I,IARM)
                   ENDDO
                ELSE
                   DO I=1,6
                     FP_VEC_P(I) = TVEC(I,IARM)
                   ENDDO
                ENDIF
C
C---------------------------------------------------------------------
C            Get JLAB Hall A HRS focal plane vector from target vector
C            and check that all apertures were passed.
C            Also, store focal plane vector for histogramming
C            Unpack array to change to COSY system,
C            corrected for the difference in coordinate systems.
C--------------------------------------------------------------------
C
C NEW COSY TRANSPORT
              ELSEIF(OP(IARM,IEL).EQ.'SIM') THEN

C Set the multiple scattering and 
C VDC resolution smearing and collimator flags

                 IF(HUT_TEST(HRS_ID_INT(IARM,IEL),1)) THEN
                    MS_FLAG =.TRUE.
                 ELSE
                    MS_FLAG = .FALSE.
                 ENDIF
                 IF(HUT_TEST(HRS_ID_INT(IARM,IEL),2)) THEN
                    WCS_FLAG =.TRUE.
                 ELSE
                    WCS_FLAG = .FALSE.
                 ENDIF
                 IF(HUT_TEST(HRS_ID_INT(IARM,IEL),3)) THEN
                    COL_FLAG =.TRUE.
                 ELSE
                    COL_FLAG = .FALSE.
                 ENDIF

C Fill the vectors and change units (cm,mrad)

                X_TG = TVEC(1,IARM)
                Y_TG = TVEC(3,IARM) 
                Z_TG = TVEC(5,IARM)
                DYDZ = TVEC(4,IARM)
                DYDZ = DYDZ/1000.D0
                DXDZ = TVEC(2,IARM)
                DXDZ = DXDZ/1000.D0
                DPP  = TVEC(6,IARM)

C Set the momentum units for 'COSY'

                IF(IARM.EQ.1) THEN
                  P_SPEC = PF_E/1000.d0
                  FRY    = X_BEAM_E/100.d0   ! convert to cm (+ = down)
                ELSE
                  P_SPEC = PF_P/1000.d0
                  FRY    = X_BEAM_P/100.d0   ! convert to cm (+ = down)
                ENDIF

                CALL HRS_COSY(p_spec,dpp,x_tg,y_tg,z_tg,dxdz,
     >             dydz,x_fp,dx_fp,y_fp,dy_fp,
     >             m2,ms_flag,wcs_flag,col_flag,fry,iarm,ok_hrs)

C Set fail global cut value

                IF(.NOT.OK_HRS) THEN
                   FAIL_GLOBAL_CUT = .TRUE.
                ENDIF

                IF(.NOT.FAIL_GLOBAL_CUT) THEN
                   TVEC(1,IARM) = X_TG
                   TVEC(3,IARM) = Y_TG
                   TVEC(5,IARM) = Z_TG
                   TVEC(4,IARM) = DYDZ*1000.d0
                   TVEC(2,IARM) = DXDZ*1000.d0
                   TVEC(6,IARM) = DPP

C Fill focal plane vectors

                   IF(IARM .EQ. 1) THEN
                      FP_VEC_E(1) = X_FP
                      FP_VEC_E(2) = DX_FP*1000.d0
                      FP_VEC_E(3) = Y_FP
                      FP_VEC_E(4) = DY_FP*1000.d0
                      FP_VEC_E(5) = 0.D0
                      FP_VEC_E(6) = DPP
                   ELSE
                      FP_VEC_P(1) = X_FP
                      FP_VEC_P(2) = DX_FP*1000.d0
                      FP_VEC_P(3) = Y_FP
                      FP_VEC_P(4) = DY_FP*1000.d0
                      FP_VEC_P(5) = 0.D0
                      FP_VEC_P(6) = DPP
                   ENDIF
                ENDIF
C
C ---------------------------------------------------------------------
C            Get JLAB Hall A MAD focal plane vector from target vector
C            and check that all apertures were passed.
C            Also, store focal plane vector for histogramming
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'MAD') THEN
                DO I=1,6
                  TVEC_TMP(I) = TVEC(I,IARM)
                ENDDO
                IMAD = MAD_CONFIG
                CALL MAD(IMAD,TVEC_TMP,TVEC_TMP2,
     #                FAIL_GLOBAL_CUT)
                DO I=1,6
                  TVEC(I,IARM) = TVEC_TMP2(I)
                ENDDO
                IF(IARM .EQ. 1) THEN
                   DO I=1,6
                     FP_VEC_E(I) = TVEC(I,IARM)
                   ENDDO
                ELSE
                   DO I=1,6
                     FP_VEC_P(I) = TVEC(I,IARM)
                   ENDDO
                ENDIF
C
C ---------------------------------------------------------------------
C            Get target vector from JLAB Hall A HRS focal plane vector.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'HRI') THEN
                DO I=1,6
                  TVEC_TMP(I) = TVEC(I,IARM)
                ENDDO
                IF(IARM .EQ. 2) THEN
                   ANGLE = PH_P      ! Needed for "bowtie" correction
                   X_BEAM = X_BEAM_P
                   Y_BEAM = Y_BEAM_P
                ELSE
                   ANGLE = PH_E      ! Needed for "bowtie" correction
                   X_BEAM = X_BEAM_E
                   Y_BEAM = Y_BEAM_E
                ENDIF
                CALL HRS_INV(HRS_ID_INT(IARM,IEL),ANGLE,X_BEAM,Y_BEAM,
     #                       TVEC_TMP,TVEC_TMP2)
                DO I=1,6
                  TVEC(I,IARM) = TVEC_TMP2(I)
                ENDDO
C
C ---------------------------------------------------------------------
C            Get target vector from Hall A HRS+septum focal plane vector.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'SPI') THEN
                DO I=1,6
                  TVEC_TMP(I) = TVEC(I,IARM)
                ENDDO
                IF(IARM .EQ. 2) THEN
                   ANGLE = PH_P      ! Needed for "bowtie" correction
                   X_BEAM = X_BEAM_P
                   Y_BEAM = Y_BEAM_P
                ELSE
                   ANGLE = PH_E      ! Needed for "bowtie" correction
                   X_BEAM = X_BEAM_E
                   Y_BEAM = Y_BEAM_E
                ENDIF
                CALL SEP_INV(HRS_ID_INT(IARM,IEL),ANGLE,X_BEAM,Y_BEAM,
     #                       TVEC_TMP,TVEC_TMP2)
                DO I=1,6
                  TVEC(I,IARM) = TVEC_TMP2(I)
                ENDDO
C
C ---------------------------------------------------------------------
C            Get target vector from JLAB Hall A MAD focal plane vector.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'MDI') THEN
                DO I=1,6
                  TVEC_TMP(I) = TVEC(I,IARM)
                ENDDO
                IF(IARM .EQ. 2) THEN
                   ANGLE = PH_P      ! Needed for "bowtie" correction
                   X_BEAM = X_BEAM_P
                   Y_BEAM = Y_BEAM_P
                ELSE
                   ANGLE = PH_E      ! Needed for "bowtie" correction
                   X_BEAM = X_BEAM_E
                   Y_BEAM = Y_BEAM_E
                ENDIF
                IMAD = MAD_INV_CONFIG
                CALL MAD_INV(IMAD,ANGLE,X_BEAM,Y_BEAM,
     #                       TVEC_TMP,TVEC_TMP2)
                DO I=1,6
                  TVEC(I,IARM) = TVEC_TMP2(I)
                ENDDO
C
C ---------------------------------------------------------------------
C            Apply JLAB Hall A HRS apertures via R-functions.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'RFN') THEN
                Y_TMP   = REAL(TVEC(3,IARM)/100.D0)
                DEL_TMP = REAL(TVEC(6,IARM)/100.D0)
                THE_TMP = REAL(TAN(TVEC(2,IARM)/1000.D0))
                PHI_TMP = REAL(TAN(TVEC(4,IARM)/1000.D0))
                TEST_RFN = DBLE(RFUNCTION(HRS_ID_INT(IARM,IEL),
     #                          Y_TMP,DEL_TMP,THE_TMP,PHI_TMP))
                IF(TEST_RFN .LT. RFN_CUT(HRS_ID_INT(IARM,IEL)))
     #                FAIL_GLOBAL_CUT = .TRUE.
                IF(IARM .EQ. 1) THEN
                   E_RFN = TEST_RFN
                ELSE
                   P_RFN = TEST_RFN
                ENDIF
C
C ---------------------------------------------------------------------
C             Reconstruct particle trajectory.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'TRK') THEN
                DO I=1,6
                  TVEC_TMP(I) = TVEC(I,IARM)
                ENDDO
                CALL GET_VECT(TVEC_TMP,ICALL,IARM)
                CALL TRACK(IARM,TVEC_TMP2)
                DO I=1,6
                  TVEC(I,IARM) = TVEC_TMP2(I)
                ENDDO
C 
C ---------------------------------------------------------------------
C            Compute relative time of flight (in nsec).
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'TOF') THEN
                IF(IARM .EQ. 1) THEN      !Beta = 1 (electron)
                  TOFE = (TVEC(5,IARM) + CRAY_PATH(IARM)*100.D0)/C
                  TOF_REL_TMP = TOFE
                ELSE
                  PMOM = PF_P*(1.D0 + TVEC(6,2)/100.D0)
                  BETA = PMOM/SQRT(PMOM**2 + EJECT_MASS**2)
                  TOFP = (TVEC(5,IARM)+CRAY_PATH(IARM)*100.D0)/(BETA*C)
                  TOF_REL_TMP = TOFE - TOFP
                ENDIF
C 
C ---------------------------------------------------------------------
C            Offset a given coordinate.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'OFF') THEN
                TVEC(NCOORD(IARM,IEL,1),IARM) =
     #                  TVEC(NCOORD(IARM,IEL,1),IARM) + OFF(IARM,IEL)
C
C ---------------------------------------------------------------------
C            Fold in resolution function to a given coordinate.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'RES') THEN
                TVEC(NCOORD(IARM,IEL,1),IARM) =
     #           TVEC(NCOORD(IARM,IEL,1),IARM) + GASDEVV(GSIG(IARM,IEL))
C
C ---------------------------------------------------------------------
C            Check various cuts.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'CUT') THEN
                IF(CUT_DOMAIN(IARM,IEL).EQ.'G' .AND.
     #               CUT_TYPE(IARM,IEL).EQ.'R') THEN
                   IF(TVEC(NCOORD(IARM,IEL,1),IARM)
     #                .LT. TR_CUT_MIN(IARM,IEL) .OR.
     #                TVEC(NCOORD(IARM,IEL,1),IARM)
     #                .GT. TR_CUT_MAX(IARM,IEL))
     #                FAIL_GLOBAL_CUT = .TRUE.
                ELSEIF(CUT_DOMAIN(IARM,IEL).EQ.'G' .AND.
     #                   CUT_TYPE(IARM,IEL).EQ.'E') THEN
                   TEST_ELL = (TVEC(NCOORD(IARM,IEL,1),IARM)/
     #                        TR_CUT_X(IARM,IEL))**2 +
     #                        (TVEC(NCOORD(IARM,IEL,2),IARM)/
     #                      TR_CUT_Y(IARM,IEL))**2
                   IF(TEST_ELL .GT. 1.D0) FAIL_GLOBAL_CUT = .TRUE.
                ELSEIF(CUT_DOMAIN(IARM,IEL).EQ.'S' .AND.
     #                   CUT_TYPE(IARM,IEL).EQ.'R') THEN
                   IF(TVEC(NCOORD(IARM,IEL,1),IARM)
     #                .GE. TR_CUT_MIN(IARM,IEL) .AND.
     #                TVEC(NCOORD(IARM,IEL,1),IARM)
     #                .LE. TR_CUT_MAX(IARM,IEL)) THEN
                         TR_CUT_FLAG(TR_CUT_IND(IARM,IEL)) = .TRUE.
                   ELSE
                         TR_CUT_FLAG(TR_CUT_IND(IARM,IEL)) = .FALSE.
                   ENDIF
                ELSEIF(CUT_DOMAIN(IARM,IEL).EQ.'S' .AND.
     #                   CUT_TYPE(IARM,IEL).EQ.'E') THEN
                   TEST_ELL = (TVEC(NCOORD(IARM,IEL,1),IARM)/
     #                        TR_CUT_X(IARM,IEL))**2 +
     #                        (TVEC(NCOORD(IARM,IEL,2),IARM)/
     #                        TR_CUT_Y(IARM,IEL))**2
                   IF(TEST_ELL .LE. 1.D0) THEN
                       TR_CUT_FLAG(TR_CUT_IND(IARM,IEL)) = .TRUE.
                   ELSE
                         TR_CUT_FLAG(TR_CUT_IND(IARM,IEL)) = .FALSE.
                   ENDIF
                ENDIF
C
C ---------------------------------------------------------------------
C            Get value of X variable for 1-D histogram.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'H1D') THEN
                TR_VAR(IARM,IEL,1) = TVEC(NCOORD(IARM,IEL,1),IARM)
C
C ---------------------------------------------------------------------
C            Get value of X and Y variables for 2-D histogram.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'H2D') THEN
                TR_VAR(IARM,IEL,1) = TVEC(NCOORD(IARM,IEL,1),IARM)
                TR_VAR(IARM,IEL,2) = TVEC(NCOORD(IARM,IEL,2),IARM)
C
C ---------------------------------------------------------------------
C            Get value of X and Y variables for scatter plot.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'SCT') THEN
                TR_VAR(IARM,IEL,1) = TVEC(NCOORD(IARM,IEL,1),IARM)
                TR_VAR(IARM,IEL,2) = TVEC(NCOORD(IARM,IEL,2),IARM)
C
C ---------------------------------------------------------------------
C            Get values of all 6 TRANSPORT coordinates for N-Tuple.
C ---------------------------------------------------------------------
C
              ELSEIF(OP(IARM,IEL).EQ.'NTU') THEN
                DO I=1,6
                  TR_VAR(IARM,IEL,I) = TVEC(I,IARM)
                ENDDO
C
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C ----------------------------------------------------------------------
C       Decide whether a histogram or scatter plot should be inremented 
C       based on the status of cuts determined above.
C ----------------------------------------------------------------------
C
      IF(ICALL .EQ. 2) THEN
        DO IARM=1,2
          IF(NEL(IARM) .GT. 0) THEN
            DO IEL=1,NEL(IARM)
              TR_INC_FLAG(IARM,IEL) = .FALSE.
              IF(OP(IARM,IEL).EQ.'H1D' .OR.
     #           OP(IARM,IEL).EQ.'H2D' .OR.
     #           OP(IARM,IEL).EQ.'SCT') THEN
                    TR_INC_FLAG(IARM,IEL) = .TRUE.  !initialize
                    IF(NCUTS_TR(IARM,IEL) .GT. 0) THEN
                      DO J=1,NCUTS_TR(IARM,IEL)
                        IF(.NOT. TR_CUT_FLAG(ICUT_IND_TR(IARM,IEL,J)))
     #                      TR_INC_FLAG(IARM,IEL) = .FALSE.
                      ENDDO
                    ENDIF
              ELSEIF(OP(IARM,IEL).EQ.'NTU') THEN
                    TR_INC_FLAG(IARM,IEL) = .TRUE.  !Don't apply specific
                                                   ! cuts to N-Tuples
              ENDIF
             ENDDO
          ENDIF
        ENDDO
      ENDIF
C
      RETURN
      END

