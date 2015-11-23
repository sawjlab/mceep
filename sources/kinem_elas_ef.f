C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       SUBROUTINE:  KINEM_ELAS_EF
C
C       AUTHOR: P.E. Ulmer
C       DATE:   10-JUL-1999
C
C       PURPOSE:
C               Calculate electron final energy for elastic scattering
C               in A(e,e'A). Also, for RADPEAKING=.TRUE., calculate the
C               Jacobian, def_dk, which may be needed to convert the
C               radiative tail cross section to the proper differential.
C
C               (This routine and KINEM_ELAS_PF replace the former
C               KINEM_ELAS.)
C
C               Allows for radiation on the electron side.  This
C               implies that the input beam vector (i.e. "pseudo-beam")
C               describes a possibly off-mass-shell electron (except
C               in peaking approximation).  Therefore, the energy and
C               3-momentum must both be passed to this routine.  For
C               the case of no radiation, the regular elastic kinematics
C               are recovered.
C
C       INPUT:
C               RADPEAKING: Calculate Jacobian?
C               EB_PS:      Energy of the pseudo-beam
C               KI_PS(3):   3-momentum of pseudo-beam
C               EB:         Beam energy     (only for Jacobian)
C               KI(3):      Beam 3-momentum (only for Jacobian)
C               KG_UNIT(3): Unit vector along photon. (only for Jacobian)
C               PHIE,THEE:  Asymptotic scattered e- angles
C
C       OUTPUT:
C               EF:      Asymptotic  scattered e- energy
C               DEF_DK:  Jacobian of scattered e- energy wrt photon energy
C
C------------------------------------------------------------------------------
C
      SUBROUTINE KINEM_ELAS_EF(RADPEAKING,EB_PS,KI_PS,EB,KI,
     #                         KG_UNIT,PHIE,THEE,EF,DEF_DK)
      IMPLICIT NONE
C
      DOUBLE PRECISION EB_PS,KI_PS(3),EB,KI(3),KG_UNIT(3),PHIE,THEE
      DOUBLE PRECISION EF,DEF_DK
      DOUBLE PRECISION PB_PS2,PS_DOT_KF,PS_MU2,KI_DOT_KG,KF_DOT_KG,DOT
      DOUBLE PRECISION KF_UNIT_VEC(3),NUM,DEN
      LOGICAL          RADPEAKING
C
      INCLUDE 'masses.cmn'
C
C------------------------------------------------------------------------------
C     Construct unit vector along asymptotic scattered electron momentum.
C------------------------------------------------------------------------------
C
      CALL V3MAKE(1.D0,PHIE,THEE,KF_UNIT_VEC)   !scatt. electron
C
C------------------------------------------------------------------------------
C     Calculate electron final energy.  Terms of order m_e^2 are
C     ignored (where m_e is the electron mass).
C
C     Quantities:
C        PB_PS2:     square of magnitude of 3-momentum for pseudo-beam.
C        PS_DOT_KF:  reduces to EB*COS(scattering_angle) for no radiation.
C        PS_MU2:     square of 4-momentum of pseudo-beam.
C        KI_DOT_KG:  (inc.   e- momentum) DOT (Photon unit vector).
C        KF_DOT_KG:  (scatt. e- momentum) DOT (Photon unit vector).
C------------------------------------------------------------------------------
C
      PB_PS2 = DOT(KI_PS,KI_PS)
      PS_DOT_KF = DOT(KI_PS,KF_UNIT_VEC)
      PS_MU2 = EB_PS**2 - PB_PS2
      NUM = EJECT_MASS*EB_PS + PS_MU2/2.D0
      DEN = EB_PS - PS_DOT_KF + EJECT_MASS
      EF = NUM/DEN                 !asymptotic e- final energy
C
      IF(RADPEAKING) THEN             !calculate Jacobian
         KI_DOT_KG = DOT(KI,KG_UNIT)
         KF_DOT_KG = EF*DOT(KF_UNIT_VEC,KG_UNIT)
         DEF_DK = (EF/NUM) * ( -EJECT_MASS - (EB-EF)
     #              + (KI_DOT_KG-KF_DOT_KG) )
      ENDIF
C
      RETURN
      END
