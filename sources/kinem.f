C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C     Subroutine KINEM:
C
C               Calculates various kinematical quantities.
C
C     Modifications:
C         10-JUL-1999 (PEU)
C         Now allows for internal radiation from the electron.
C         The photon can be emitted in any direction
C         (i.e. not necessarily peaking approximation).
C         This implies that the virtual electron
C         is no longer on its mass shell:  for radiated photons of
C         high energy, the electron energy and magnitude of
C         3-momentum are no longer equal.  So now the 3-momenta and
C         energies of the beam and scattered electron must both
C         be passed into this routine.
C
C     INPUT:
C            EB:           Energy of beam
C            KI_VEC(3):    3-momentum of beam
C            EF:           Energy of scattered electron
C            KF_VEC(3):    3-momentum of scattered electron
C            PHP,THP:      Angles of ejectile
C
C            PFP:          Magnitude of 3-momentum of ejectile
C                            (calculated for bound state case)
C
C ------------------------------------------------------------------
C
      SUBROUTINE KINEM(EB,KI_VEC,EF,KF_VEC,PFP,PHP,THP,
     #     ELASTIC,BOUND,JUSTPF,FAIL_KIN)
C
      IMPLICIT NONE
C
      COMMON /KINVAR/  KI,Q,Q2,QMU2,EP,PF,THETA_EP,CTH_PQ,PR
      COMMON /PHYSVAR/ IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
      COMMON /SURVIVE_L/ DECAY_CHECK
      COMMON /SURVIVE_D/ SURVIVE_PROB
C
      DOUBLE PRECISION EB,KI_VEC(3),EF,KF_VEC(3),PFP,PHP,THP
      DOUBLE PRECISION PB,PFE
      DOUBLE PRECISION KI(3),Q(3),Q2,QMU2,EP,PF(3),THETA_EP,CTH_PQ
      DOUBLE PRECISION PR(3),Q_X_KP(3),KPERP(3),SURVIVE_PROB
      DOUBLE PRECISION PI,MD,MP
      DOUBLE PRECISION DOT,ARG_W,TEST_CTH_PQ,BETA,GAMMA
      DOUBLE PRECISION E_CM,P_CM,P_CM_PARA,P_CM_PERP,QDOTPR
      DOUBLE PRECISION CTH_PRQ,TEST_CTH_PRQ,KIDOTQ,ARGNUM,ARGDEN
cxxx      DOUBLE PRECISION ROT(3,3),PR_ROT(3),CANG,SANG
C
      INTEGER IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL,I
      LOGICAL FAIL_KIN,ELASTIC,BOUND,JUSTPF,DECAY_CHECK
C
      PARAMETER (PI=3.14159265359D0)
      PARAMETER (MD = 1875.6D0)               !Deuteron mass
      PARAMETER (MP = 938.28D0)               !Proton mass
C
C ---------------------------------------------------------------------
C       Arrays for plots and cuts.
C ---------------------------------------------------------------------
C
      INCLUDE 'var.cmn'
      INCLUDE 'masses.cmn'
      INCLUDE 'lifetimes.cmn'
C
      INCLUDE 'lifetimes_dat.cmn'
C
      FAIL_KIN = .FALSE.                      !initialize
C
      DO I=1,3
         KI(I) = KI_VEC(I)
      ENDDO
C
C ---------------------------------------------------------------------
C       Compute magnitudes of 3-momenta of initial and scattered
C       electrons and scattering angle.
C ---------------------------------------------------------------------
C
      PB  = SQRT(DOT(KI_VEC,KI_VEC))
      PFE = SQRT(DOT(KF_VEC,KF_VEC))
      TSCAT = ACOS(DOT(KI_VEC,KF_VEC)/(PB*PFE))      !Scattering angle
C
C ---------------------------------------------------------------------
C       Calculate Q, Omega and W.
C ---------------------------------------------------------------------
C
      OMEGA = EB - EF                         !energy transfer
      CALL V3DIF(KI_VEC,KF_VEC,Q)             !3-momentum transfer vector
      Q2 = DOT(Q,Q)
      QMAG = SQRT(Q2)
      QMU2 = Q2-OMEGA**2                      !4-mom. transf.^2
      QMU2_G = QMU2/1.D6                      !4-mom. transf.^2 in GeV^2/c^2
      EPSILON = 1.D0/(1.D0+2.D0*QMAG**2*TAN(TSCAT/2.D0)**2/QMU2)
      X = QMU2/(2.D0*MP*OMEGA)
      THEQ = ASIN(Q(2)/QMAG)                  !Q angle wrt floor
      PHIQ = ATAN(Q(1)/Q(3))                  !Q angle wrt beam (in-plane)
      ARG_W = MNUC**2+2.D0*OMEGA*MNUC-QMU2
      IF(ARG_W .GE. 0.D0) THEN
        W = SQRT(ARG_W)            !invariant mass of hadronic final state
      ELSE
        FAIL_KIN = .TRUE.
        GOTO 999
      ENDIF
C
C ---------------------------------------------------------------------
C       Get ejectile momentum vector and total energy in lab.
C ---------------------------------------------------------------------
C
      CALL V3MAKE(1.0D0,PHP,THP,PF)           !unit vector along PF
      CTH_PQ = DOT(PF,Q)/(QMAG)
      IF(BOUND) THEN
        IF(W .LT. MNUC+MISS_M) THEN
           FAIL_KIN = .TRUE.
           GOTO 999
        ENDIF
        CALL PFIN(MISS_M,OMEGA,QMAG,CTH_PQ,PFP,FAIL_KIN)
        IF(FAIL_KIN .OR. JUSTPF) GOTO 999
      ENDIF
      DO I=1,3
        PF(I) = PF(I)*PFP
      ENDDO
      EP = SQRT(PFP**2 + EJECT_MASS**2)
      BETA_H = PFP/EP     ! v/c for hadron
C
C ---------------------------------------------------------------------
C       Get missing energy.
C ---------------------------------------------------------------------
C
      EMISS = OMEGA - (EP-EJECT_MASS)
C
C ---------------------------------------------------------------------
C       Get angle between ejectile and Q and between ejectile and beam.
C ---------------------------------------------------------------------
C
      TEST_CTH_PQ = 1.D0 - CTH_PQ
      IF(TEST_CTH_PQ .LT. 0.D0 .AND. TEST_CTH_PQ .GT. -1.D-12)
     #     CTH_PQ = 1.D0   !protect against roundoff (esp. for elastic)
      THETA_PQ = ACOS(CTH_PQ)             !angle between ejectile and Q
      THETA_EP = ACOS(DOT(KI,PF)/(PB*PFP)) !angle between ejectile and beam
C
      T_GEV2 = (QMAG-PFP*COS(THETA_PQ))**2 
     #            + (PFP*SIN(THETA_PQ))**2 - (OMEGA-EP)**2
      T_GEV2 = T_GEV2/1.E6                !convert to (GeV//c)**2
C
C ---------------------------------------------------------------------
C       Calculate quantities in Center Of Mass system.
C ---------------------------------------------------------------------
C
      BETA  = QMAG/(OMEGA+MNUC)              !center of mass velocity
      GAMMA = (OMEGA+MNUC)/W                 !gamma = sqrt(1/(1-beta^2))
C
      E_CM = GAMMA*(EP-BETA*PFP*CTH_PQ)      !ejectile COM energy
      P_CM_PARA = GAMMA*(PFP*CTH_PQ-BETA*EP) !component of P_CM along Q
      P_CM_PERP = PFP*SIN(THETA_PQ)          !component of P_CM perp. to Q
      P_CM = SQRT(P_CM_PARA**2+P_CM_PERP**2) !ejectile COM momentum
      IF(P_CM .NE. 0.D0) THEN
        THETA_CM = ACOS(P_CM_PARA/P_CM)      !angle of ejectile wrt Q in COM
      ELSE
        THETA_CM = 0.D0   !choose arbitrarily (may be necessary for elastic)
      ENDIF
C
C ---------------------------------------------------------------------
C       Get recoil momentum.
C ---------------------------------------------------------------------
C
      CALL V3DIF(Q,PF,PR)
      PRMAG = SQRT(DOT(PR,PR))
      QDOTPR = DOT(Q,PR)
      PREC = -PRMAG*DSIGN(1.D0,QDOTPR)        !neg. along Q
C
C ---------------------------------------------------------------------
C       Get recoil momentum components in the Laboratory system.
C
C       The "commented out" lines below roughly give the components
C       of PREC in a basis with z-axis along Q.  However, since Q
C       may not be in the plane of the Laboratory floor, a rotation
C       about Y by the angle between Q and the beam does not exactly
C       give the correct components.  It's easier to express PREC
C       in the Lab basis, so do that instead.
C ---------------------------------------------------------------------
C
      PREC_X = PR(1)
      PREC_Y = PR(2)
      PREC_Z = PR(3)
C
      PREC_TH = ASIN(PR(2)/PRMAG)              !Pr angle wrt floor
C
      IF(PR(3).EQ. 0.D0) THEN
         PREC_PH = PI/2.D0
      ELSE
         PREC_PH = ABS(ATAN(PR(1)/PR(3)))      !Pr angle wrt beam (in-plane)
      ENDIF
      IF(PR(1).GE.0.D0 .AND. PR(3).GE. 0.D0) THEN
         PREC_PH = PREC_PH
      ELSEIF(PR(1).GE.0.D0 .AND. PR(3).LT. 0.D0) THEN
         PREC_PH = PI - PREC_PH
      ELSEIF(PR(1).LT.0.D0 .AND. PR(3).LT. 0.D0) THEN
         PREC_PH = PREC_PH - PI
      ELSEIF(PR(1).LT.0.D0 .AND. PR(3).GE. 0.D0) THEN
         PREC_PH = -PREC_PH
      ENDIF
C
cxxx      CANG = DOT(Q,KI_VEC)/(QMAG*PB)
cxxx      SANG = SIN(ACOS(CANG))
cxxx      CALL ROTATE_Y(CANG,-SANG,ROT)
cxxx      CALL MAT_MULT(3,3,1,ROT,PR,PR_ROT)  ! PR_ROT wrt Q axis now
cxxx      PREC_X = PR_ROT(1)
cxxx      PREC_Y = PR_ROT(2)
cxxx      PREC_Z = PR_ROT(3)
C
C ---------------------------------------------------------------------
C       Get angle between Q and recoil momentum.
C ---------------------------------------------------------------------
C
      IF(PRMAG .EQ. 0.D0) THEN
         THETA_PRQ = 0.D0
      ELSE
         CTH_PRQ = QDOTPR/(QMAG*PRMAG)
         TEST_CTH_PRQ = 1.D0 - CTH_PRQ
         IF(TEST_CTH_PRQ .LT. 0.D0 .AND. TEST_CTH_PRQ .GT. -1.D-12)
     #     CTH_PRQ =  1.D0   !protect against roundoff
         TEST_CTH_PRQ = 1.D0 + CTH_PRQ
         IF(TEST_CTH_PRQ .LT. 0.D0 .AND. TEST_CTH_PRQ .GT. -1.D-12)
     #     CTH_PRQ = -1.D0   !protect against roundoff
         THETA_PRQ = ACOS(CTH_PRQ)     !angle between recoil and Q
      ENDIF
C
C ---------------------------------------------------------------------
C       Get recoil factor (see SUBROUTINE FREC for an explanation).
C       For continuum get missing mass as well.
C ---------------------------------------------------------------------
C
      CALL FREC(BOUND,PR,PF,OMEGA,EP,PFP,RECFAC)
      IF(.NOT. BOUND) THEN
        CALL MISSM(ELASTIC,OMEGA,PFP,PREC,MISS_M,MMSQ,FAIL_KIN)
        IF(FAIL_KIN) GOTO 999
      ELSE
        MMSQ=MISS_M*MISS_M
      ENDIF
C
C ---------------------------------------------------------------------
C       Calculate azimuthal angle, PHI_X.
C ---------------------------------------------------------------------
C
      KIDOTQ = DOT(KI,Q)/QMAG
      DO I=1,3
        KPERP(I) = KI(I) - KIDOTQ*Q(I)/QMAG
      ENDDO
      Q_X_KP(1) = (Q(2)*KPERP(3) - Q(3)*KPERP(2))/QMAG
      Q_X_KP(2) = (Q(3)*KPERP(1) - Q(1)*KPERP(3))/QMAG
      Q_X_KP(3) = (Q(1)*KPERP(2) - Q(2)*KPERP(1))/QMAG
      ARGNUM = DOT(PF,Q_X_KP)
      ARGDEN = DOT(PF,KPERP)
      IF(ARGDEN .NE. 0.)PHI_X = ATAN(ARGNUM/ARGDEN)
      IF(ARGDEN .LT. 0.)THEN
        PHI_X = PHI_X + PI
      ELSEIF((ARGNUM .LT. 0.) .AND. (ARGDEN .GT. 0.))THEN
        PHI_X = PHI_X + 2.*PI
      ELSEIF((ARGNUM .LT. 0.) .AND. (ARGDEN .EQ. 0.))THEN
        PHI_X = 1.5D0*PI
      ELSEIF((ARGNUM .GT. 0.) .AND. (ARGDEN .EQ. 0.))THEN
        PHI_X = 0.5D0*PI
      ELSEIF((ARGNUM .EQ. 0.) .AND. (ARGDEN .EQ. 0.))THEN
        PHI_X = 0.
      ENDIF
C
C ---------------------------------------------------------------------
C     For p(e,e'pi+)n or p(e,e'K)Lambda calculate the ejectile
C     survival probability.
C ---------------------------------------------------------------------
C
      IF(DECAY_CHECK) THEN
         IF(IPHYS_OPT .EQ. 600) THEN
            CALL DECAY(LIFE_PION)
         ELSEIF(IPHYS_OPT .EQ. 700) THEN
            CALL DECAY(LIFE_KAON)
         ENDIF
         SURVIVAL_PROB = SURVIVE_PROB
      ENDIF
C
  999 RETURN
      END
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C       Subroutine PFIN:
C
C               Calculates ejectile final momentum for bound state.
C
C               In some cases there are two physical solutions for
C               the momentum, one corresponding to the ejectile
C               moving along q in the Center-of-Mass and the other
C               against q.  In this case a coin is tossed to decide
C               which root to take and the event carries a weight
C               of 2 in calculating the cross section.
C ------------------------------------------------------------------
C
      SUBROUTINE PFIN(MISS_M,OMEGA,QMAG,CTH_PQ,PF_P_I,FAIL_KIN)
C
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
      COMMON /PLIMS/ PF_MIN,PF_MAX
      COMMON /PFROOTS_D/ PF_ROOT_WT
      COMMON /PFROOTS_L/ FIRST_TIME
      COMMON /PFROOTS_I/ I_PF_ROOT
      REAL RAN_UNIF
      LOGICAL FAIL_KIN,PLUS,MINUS,FIRST_TIME
C
      INCLUDE 'masses.cmn'
C
      SAVE PLUS,MINUS
C
      FAIL_KIN = .FALSE.         !initialize
C
      Q2 = QMAG*QMAG      !q-3vector squared
C
      CONST1 = (MISS_M + MNUC - EJECT_MASS)**2   !constants to determine
      CONST2 = (OMEGA + MNUC)**2                 !   ejectile momentum
      CONST3 = Q2-CONST2+CONST1-EJECT_MASS**2
C
      A1 = 4.*(CONST2-Q2*CTH_PQ**2)           !coeff's. of quadratic
      A2 = 4.*QMAG*CONST3*CTH_PQ              !   to determine
      A3 = 4.*CONST2*EJECT_MASS**2-CONST3**2  !   ejectile momentum
C
C ------------------------------------------------------------------
C     Find all real roots of quadratic equation for momentum.
C ------------------------------------------------------------------
C
      DISCRIM = A2**2-4.D0*A1*A3
      IF(DISCRIM .LT. 0.) THEN                 !no physical solutions
        FAIL_KIN = .TRUE.
        GOTO 999
      ELSE
        FACT = 1.D0/(2.D0*A1)
        SQRT_DISC = SQRT(DISCRIM)
        PF_PLUS  = FACT*(-A2 + SQRT_DISC)
        PF_MINUS = FACT*(-A2 - SQRT_DISC)
      ENDIF
C
C ------------------------------------------------------------------
C     If this is the first call to this routine for this event
C     then decide which roots lie within hadron-arm momentum
C     acceptance (ejectile momentum is also constrained to be >0).
c     If both roots are within acceptance then toss a coin to decide
C     which one to choose and set logical flags accordingly.
C     In this case the event must carry an additional
C     weighting factor of 2 since each root is selected half the
C     time.
C
C     If this is not the first call to this routine for this event
C     select same root as before (to prevent switching solutions
C     in mid-event).
C ------------------------------------------------------------------
C
      IF(FIRST_TIME) THEN          !first call for this event
        PLUS     = .FALSE.         !   initialize
        MINUS    = .FALSE.         !   flags
        IF((PF_PLUS .GE.PF_MIN) .AND. (PF_PLUS .LE.PF_MAX)
     #                          .AND. (PF_PLUS .GE.0.D0)) PLUS  = .TRUE.
        IF((PF_MINUS.GE.PF_MIN) .AND. (PF_MINUS.LE.PF_MAX)
     #                          .AND. (PF_MINUS.GE.0.D0)) MINUS = .TRUE.
        PF_ROOT_WT = 1.D0
        IF(PLUS .AND. MINUS) THEN        !both roots pass - toss coin
            CALL RANECU(RAN_UNIF,1)
            ROOT_TEST = RAN_UNIF
            IF(ROOT_TEST .GE. 0.5D0) THEN
              MINUS = .FALSE.
            ELSE
              PLUS  = .FALSE.
            ENDIF
            PF_ROOT_WT = 2.D0
        ENDIF
      ELSE
         IF(I_PF_ROOT .EQ. 1) THEN
            PLUS  = .TRUE.
            MINUS = .FALSE.
         ELSE
            PLUS  = .FALSE.
            MINUS = .TRUE.
         ENDIF
      ENDIF
C
C ------------------------------------------------------------------
C     Select ejectile momentum based on logical flags determined
C     in the first call to this routine for this event.
C ------------------------------------------------------------------
C
      IF(PLUS .AND. .NOT. MINUS) THEN       !only + root passes
        PF_P_I = PF_PLUS
        I_PF_ROOT = 1
      ELSEIF(MINUS .AND. .NOT. PLUS) THEN   !only - root passes
        PF_P_I = PF_MINUS
        I_PF_ROOT = 2
      ELSEIF((.NOT. MINUS) .AND. (.NOT. PLUS)) THEN
        FAIL_KIN = .TRUE.                   !no acceptable solution
        I_PF_ROOT = 0
        GOTO 999
      ENDIF
C
C ------------------------------------------------------------------
C     Test to see whether PF_P_I is within hadron arm momentum
C     acceptance.  This test must be performed for ALL calls to
C     this routine since the kinematics can change slightly from
C     call to call for a given event due to the effect of various
C     errors (offsets, etc.) introduced via SPECTROMETER.
C ------------------------------------------------------------------
C
      IF(.NOT. FIRST_TIME) THEN  !perform test for subsequent calls
        IF((PF_P_I.LT.PF_MIN).OR.(PF_P_I.GT.PF_MAX))FAIL_KIN = .TRUE.
      ENDIF
C
  999 FIRST_TIME = .FALSE.
      RETURN
      END
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C       Subroutine FREC:
C
C           Calculate recoil factor for bound state or continuum.
C
C       The recoil factors assume that the cross section starts out
C       being differential in the hadron energy.  For a bound state
C       the recoil factor is required to do the integration over
C       missing mass to arrive at a 5-fold differential cross section.
C       For the continuum the "recoil factor" converts the 6-fold
C       differential cross section from being differential in hadron
C       energy to being differential in hadron momentum.  This is
C       required since the momentum, rather than the energy, is
C       sampled within the main Monte Carlo loop.
C ------------------------------------------------------------------
C
      SUBROUTINE FREC(BOUND,PR,PF,OMEGA,EP,PF_P_I,RECFAC)
C
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
      DOUBLE PRECISION PF(3),PR(3)
      LOGICAL BOUND
C
      INCLUDE 'masses.cmn'
C
      IF(BOUND) THEN
        RECFAC = ABS(1.D0-EP*DOT(PR,PF)/((OMEGA-EP+MNUC)*PF_P_I**2))
      ELSE
        RECFAC = EP/PF_P_I
      ENDIF
C
      RETURN
      END
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C       Subroutine MISSM:
C
C               Calculates missing mass for continuum.
C ------------------------------------------------------------------
C
      SUBROUTINE MISSM(ELASTIC,OMEGA,PF_P_I,PREC,MISS_M,MMSQ,FAIL_KIN)
C
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
      LOGICAL ELASTIC,FAIL_KIN
C
      INCLUDE 'masses.cmn'
C
      FAIL_KIN = .FALSE.                      !initialize
C
C ------------------------------------------------------------------
C     For elastic scattering:
C
C         The missing mass is simply Omega-T (where T is the kinetic
C         energy of the recoil nucleus.
C
C     For bound state:
C
C         If RADICAL < 0 there is no physical solution for the
C         missing mass.  However, to avoid having events being
C         thrown out due to round-off resulting in negative
C         RADICALs, take absolute value if RADICAL is negative
C         but close to zero.
C    
C     Note: Changed from original 3.0 version by F. Sabatie to take
C           into account that the missing mass can be negative due
C           to resolution effects. Cut is set to MM^2 at 5000 MeV^2
C           for cases where the target and ejectile masses are equal.
C ------------------------------------------------------------------
C
      IF(ELASTIC) THEN
        MISS_M = OMEGA-SQRT(PF_P_I**2+EJECT_MASS**2)+MNUC
        MMSQ   = MISS_M * ABS(MISS_M)  ! make signed quantity
      ELSE
        RADICAL = (OMEGA-SQRT(PF_P_I**2+EJECT_MASS**2)+MNUC)**2-PREC**2
        IF(ABS(MNUC-EJECT_MASS) .LT. 0.01 D0) THEN
          IF(RADICAL .LT. -5000.D0) THEN         ! resolution effects
            FAIL_KIN = .TRUE.                    ! can drive the MM negative
            GOTO 999                             ! in VCS
          ELSE
            MMSQ    = RADICAL     ! make signed quantity
            RADICAL = ABS(RADICAL)
            MISS_M  = SQRT(RADICAL)
          ENDIF
        ELSE
          IF(RADICAL .LT. -1.D-20) THEN
            FAIL_KIN = .TRUE.
            GOTO 999
          ELSE
            RADICAL = ABS(RADICAL)
            MISS_M  = SQRT(RADICAL)-MNUC+EJECT_MASS
            MMSQ    = MISS_M * MISS_M 
          ENDIF
        ENDIF
      ENDIF
C
  999 RETURN
      END
