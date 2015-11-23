C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       SUBROUTINE:  KINEM_ELAS_PF
C
C       AUTHOR: P.E. Ulmer
C       DATE:   10-JUL-1999
C
C       PURPOSE:
C               Calculate momentum vector for recoiling nucleus for
C               elastic scattering in A(e,e'A).
C               (This routine and KINEM_ELAS_EF replace the former
C               KINEM_ELAS.)
C
C               This routine requires that the electron final energy
C               has already been calculated (from KINEM_ELAS_EF),
C               from the elastic kinematics (with or without radiation).
C
C
C       INPUT:
C               KI_PS(3):  3-momentum of pseudo-beam
C               KF(3):     Asymptotic scattered e- 3-momentum
C
C       OUTPUT:
C               PF:        Momentum of recoiling nucleus ("ejectile")
C               PHIP,THEP  Angles of recoiling nucleus 
C
C------------------------------------------------------------------------------
C
      SUBROUTINE KINEM_ELAS_PF(KI_PS,KF,
     #                      PF,PHIP,THEP)
      IMPLICIT NONE
C
      DOUBLE PRECISION KI_PS(3),KF(3),PF,PHIP,THEP
      DOUBLE PRECISION PF_VEC(3),DOT
      INTEGER          I
C
C------------------------------------------------------------------------------
C     Compute the momentum of recoiling nucleus, which is equal
C     to the momentum transfer, Q.
C------------------------------------------------------------------------------
C
      DO I=1,3
         PF_VEC(I) = KI_PS(I) - KF(I)  ! 3-momentum of recoiling nucleus
      ENDDO
C
C------------------------------------------------------------------------------
C     Calculate magnitude of momentum for recoiling nucleus as well
C     as angles in the LAB system.
C------------------------------------------------------------------------------
C
      PF = SQRT(DOT(PF_VEC,PF_VEC))
      THEP = ASIN(PF_VEC(2)/PF)
      PHIP = ACOS(PF_VEC(3)/(PF*COS(THEP)))
      IF(PF_VEC(1) .LT. 0.D0) PHIP = - PHIP
C
      RETURN
      END



