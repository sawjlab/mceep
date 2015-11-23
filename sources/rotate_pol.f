C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       SUBROUTINE:  ROTATE_POL
C       AUTHOR:      P.E. Ulmer
C       DATE:        4-DEC-1990
C
C       PURPOSE:
C               Determine polarization components in the (fixed) spectrometer
C               reference frame from the components in the reaction frame
C               (used in the calculation of these observables).
C               As of version 3.4:  also calculate components in
C               electron scattering plane (S. Strauch).
C
C       Four right-handed reference frames are defined (SP, L, R and S)
C       as follows.
C
C          The SCATTERING PLANE frame, SP, is given by:
C
C          S          X - YxZ
C           C         Y - along k_i x k_f, where k_i and k_f are the
C            A            the incident and scattered electron momenta
C             T       Z - along Q direction
C
C          The LABORATORY frame, L, (the same frame used by MCEEP for all
C          kinematical vectors) is given by:
C
C          L          X - YxZ
C           A         Y - normal to laboratory floor (upwards = +)
C            B        Z - along nominal beam direction (i.e. for PH_B=TH_B=0)
C
C          The SPECTROMETER frame, S, (also used by MCEEP for the Transport
C          vectors) is given by:
C
C          S
C           P         X - along momentum dispersion
C            E        Y - ZxX
C             C       Z - along spectrometer central ray
C              T
C
C          The REACTION frame, R, is defined as follows:
C
C          R
C           E
C            A        X ---> normal to the reaction plane
C             C       Y ---> Z x X
C              T      Z ---> along the hadron momentum
C               I    
C                O
C                 N
C
C                     Note that by convention the X-axis points upward
C                     for coplanar kinematics with phi_x = 0 (i.e. hadrons
C                     forward of q).
C
C------------------------------------------------------------------------------
C
      SUBROUTINE ROTATE_POL(ELASTIC,PHI_X,THE_NQ,KI,Q,PF,
     &                      PN,PT,PL,PN_S,PT_S,PL_S,PX,PY,PZ)
      IMPLICIT NONE
C
      COMMON /L_TO_S/ ROT_P
C
      DOUBLE PRECISION PN,PT,PL,PN_S,PT_S,PL_S,PX,PY,PZ
      DOUBLE PRECISION PHI_X, THE_NQ
      DOUBLE PRECISION POL_REACT(3),POL_SPECT(3),POL_LAB(3)
      DOUBLE PRECISION KI(3),Q(3),PF(3),PF_MAG,DOT
      DOUBLE PRECISION Q_X_KI(3),Q_X_PF(3),Q_X_KI_MAG,Q_X_PF_MAG
      DOUBLE PRECISION R(3,3),ROT_P(3,3)
      DOUBLE PRECISION R1(3),R2(3),R3(3)
      DOUBLE PRECISION CST, CSP, SNT, SNP
      INTEGER I,J
      LOGICAL ELASTIC
C
      EQUIVALENCE (R(1,1),R1(1))
      EQUIVALENCE (R(1,2),R2(1))
      EQUIVALENCE (R(1,3),R3(1))
C
C------------------------------------------------------------------------------
C     Rotate from the REACTION to the SCATTERING PLANE frame
C     (S. Strauch).
C------------------------------------------------------------------------------
C
      IF (ELASTIC) THEN     ! No rotation necessary
         PX = PT
         PY = PN
         PZ = PL
      ELSE
         CSP  = COS(PHI_X)
         SNP  = SIN(PHI_X)
         CST  = COS(THE_NQ)
         SNT  = SIN(THE_NQ)
         PX = PL * SNT * CSP + PT * CST * CSP - PN * SNP
         PY = PL * SNT * SNP + PT * CST * SNP + PN * CSP
         PZ = PL * CST       - PT * SNT
      ENDIF
C
C------------------------------------------------------------------------------
C     Now do other rotations successively starting at REACTION frame
C     and ending at SPECTROMETER frame.
C------------------------------------------------------------------------------
C
      POL_REACT(1) = PN
      POL_REACT(2) = PT
      POL_REACT(3) = PL
C
C------------------------------------------------------------------------------
C     First rotate from the REACTION to the LABORATORY frame.
C
C     R1, R2 and R3 are three unit vectors pointing along the
C     X,  -Y  and Z axes of the REACTION frame respectively.
C     R1, R2 and R3 correspond to the usual N, T and L polarizations
C     respectively.  The minus sign for R2 is incorporated to
C     make N,L,T a right-handed system (rather than N,T,L) which
C     is the usual convention (i.e. N = L x T).
C     The (i,j) element of the rotation matrix R corresponds
C     to the ith component of the jth R unit vector (this is
C     guaranteed by the EQUIVALENCE statements above).
C
C     NOTE:  For elastic scattering there is only one plane,
C            the electron scattering plane.  Therefore, for this
C            case, the "reaction" vectors refer to the scattering
C            plane.
C------------------------------------------------------------------------------
C
      PF_MAG = SQRT(DOT(PF,PF))             !Mag. of hadron momentum
      DO I=1,3
        R3(I) = PF(I)/PF_MAG
      ENDDO
C
      IF(ELASTIC) THEN
        CALL VEC_CROSS(Q,KI,Q_X_KI)           !Get q x ki
        Q_X_KI_MAG = SQRT(DOT(Q_X_KI,Q_X_KI)) !Mag. of the cross product
        DO I=1,3
          R1(I) = Q_X_KI(I)/Q_X_KI_MAG
        ENDDO
      ELSE
        CALL VEC_CROSS(Q,PF,Q_X_PF)           !Get q x p
        Q_X_PF_MAG = SQRT(DOT(Q_X_PF,Q_X_PF)) !Mag. of the cross product
        DO I=1,3
          R1(I) = Q_X_PF(I)/Q_X_PF_MAG
        ENDDO
      ENDIF
C
      CALL VEC_CROSS(R3,R1,R2)
C
C------------------------------------------------------------------------------
C     Reverse the direction of R2 to make (N,L,T) right-handed (see
C     comments above).
C------------------------------------------------------------------------------
C
      DO I=1,3
        R2(I) = -R2(I)
      ENDDO
      DO I=1,3
        POL_LAB(I) = 0.D0
        DO J=1,3
          POL_LAB(I) = POL_LAB(I) + R(I,J)*POL_REACT(J)
        ENDDO
      ENDDO
C
C------------------------------------------------------------------------------
C     Now, rotate from the LABORATORY to the SPECTROMETER frame.
C------------------------------------------------------------------------------
C
      DO I=1,3
        POL_SPECT(I) = 0.D0
        DO J=1,3
          POL_SPECT(I) = POL_SPECT(I) + ROT_P(I,J)*POL_LAB(J)
        ENDDO
      ENDDO
C
      PN_S = POL_SPECT(1)
      PT_S = POL_SPECT(2)
      PL_S = POL_SPECT(3)
C
      RETURN
      END
