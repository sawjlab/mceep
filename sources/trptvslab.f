C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C       Subroutine LAB_TO_TRPT:
C
C       AUTHOR:      P.E. Ulmer
C       DATE:        9-JUL-1991
C       MODIFICATIONS:
C           Formerly calculated the "Transport" angles as
C
C                TRPT_T(2) = (QI(2)-TH)*1000.          (Theta_Transport)
C                TRPT_T(4) = (QI(4)-PH)*CTH*1000.      (Phi_Transport)
C
C           where CTH is the cosine of the spectrometer out-of-plane
C           angle and TH and PH are the spectrometer central angles.
C           These approximations, although accurate for small
C           acceptances, are unreliable for moderately large
C           acceptances.  The current method gives the exact Transport
C           angles from the Laboratory angles.  Thus, all the
C           approximations are made in the spectrometer field map
C           whereas the coordinate transformations are exact.
C           Note that to obtain the Transport coordinates one must
C           perform a traceback of the ray to the Transport XY plane.
C
C       PURPOSE:
C          Calculates Transport coordinates from Laboratory coordinates.
C
C          Conventions for the coordinate systems are:
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
C ---------------------------------------------------------------------
C
      SUBROUTINE LAB_TO_TRPT(QI,ROT,PF,TRPT_T)
C
      IMPLICIT NONE
      DOUBLE PRECISION  QI(6),TRPT_T(6),ROT(3,3),PVEC_L(3),POS_T(3)
      DOUBLE PRECISION  PF,CTHL,STHL,CPHL,SPHL,TTHET,TPHIT
      DOUBLE PRECISION  PVEC_T(3)
      INTEGER           I,J,INDJ
C
C ---------------------------------------------------------------------
C       Get sines and cosines of LAB angles.
C ---------------------------------------------------------------------
C
      CTHL = COS(QI(2))
      STHL = SIN(QI(2))
      CPHL = COS(QI(4))
      SPHL = SIN(QI(4))
C
C ---------------------------------------------------------------------
C       Construct unit momentum vector in LAB system.
C ---------------------------------------------------------------------
C
      PVEC_L(1) = CTHL*SPHL
      PVEC_L(2) = STHL
      PVEC_L(3) = CTHL*CPHL
C
C ---------------------------------------------------------------------
C       Rotate to get unit momentum vector in Transport system.
C ---------------------------------------------------------------------
C
      DO I=1,3
        PVEC_T(I) = 0.D0           !initialize
        DO J=1,3
          PVEC_T(I) = PVEC_T(I) + ROT(I,J)*PVEC_L(J)
        ENDDO
      ENDDO
C
C ---------------------------------------------------------------------
C       Calculate particle angles in spectrometer (Transport)
C       coordinate system in mr.
C ---------------------------------------------------------------------
C
      TTHET = PVEC_T(1)/PVEC_T(3)       !tan(theta_T)
      TPHIT = PVEC_T(2)/PVEC_T(3)       !tan(phi_T)
      TRPT_T(2) =  1000.D0*ATAN(TTHET)  !theta_T in mr
      TRPT_T(4) =  1000.D0*ATAN(TPHIT)  !phi_T   in mr
C
C ---------------------------------------------------------------------
C       Calculate dispersion in percent.
C ---------------------------------------------------------------------
C
      TRPT_T(6) = (QI(6)-PF)*100.D0/PF
C
C ---------------------------------------------------------------------
C       Calculate position of interaction point in spectrometer
C       (Transport) coordinate systems in cm.
C ---------------------------------------------------------------------
C
      DO I=1,3
        POS_T(I) = 0.D0            !Initialize
        DO J=1,3
          INDJ = 2*J-1
          POS_T(I) = POS_T(I) + 100.D0*ROT(I,J)*QI(INDJ)
        ENDDO
      ENDDO
C
C ---------------------------------------------------------------------
C      Obtain positional (X,Y) components of Transport vector after 
C      traceback to the XY plane of the rotated system.
C
C      NOTE:  POS_T(1), POS_T(2) and POS_T(3) are the X, Y and Z
C             positional coordinates in the Transport system
C             before traceback, whereas TRPT_T(1) and TRPT_T(3)
C             are the X and Y Transport coordinates (i.e. after
C             traceback).  Also, TRPT_T(5) is the path length and
C             NOT the Z-coordinate.
C ---------------------------------------------------------------------
C
      TRPT_T(1) = POS_T(1)-POS_T(3)*TTHET
      TRPT_T(3) = POS_T(2)-POS_T(3)*TPHIT
C
C ---------------------------------------------------------------------
C     Determine differential path length corresponding to the traceback.
C ---------------------------------------------------------------------

      TRPT_T(5) = - POS_T(3)*SQRT(1.D0 + TTHET*TTHET + TPHIT*TPHIT)
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C       Subroutine TRPT_TO_LAB:
C
C       AUTHOR:      P.E. Ulmer
C       DATE:        9-JUL-1991
C       MODIFICATIONS:
C           Formerly calculated the "Laboratory" angles as
C
C                QI(2) = 0.001*TRPT_T(2) + TH       (Theta_Lab)
C                QI(4) = 0.001*TRPT_T(4)/CTH + PH   (Phi_Lab)
C
C           where CTH is the cosine of the spectrometer out-of-plane
C           angle and TH and PH are the spectrometer central angles.
C           These approximations, although accurate for small
C           acceptances, are unreliable for moderately large
C           acceptances.  The current method gives the exact Laboratory
C           angles from the Transport angles.  Thus, all the
C           approximations are made in the spectrometer field map
C           whereas the coordinate transformations are exact.
C
C           The positional coordinates returned by this routine
C           are determined by performing a rotation of the
C           Transport X and Y (in the Transport plane Z=0 by
C           definition) coordinates without any traceback to the
C           actual interaction point.  This is because with only
C           information from the spectrometer focal plane, one
C           can only determine the intersection of the ray with
C           the Transport XY plane and NOT the interaction point.
C
C       19-FEB-2004    PEU
C           Now calculate the actual interaction point, as opposed
C           to the intersection of the particle ray with the
C           target Transport system (Z=0).  This change made for
C           version 3.9.
C
C
C       PURPOSE:
C          Calculates Laboratory coordinates from Transport coordinates.
C
C          Conventions for the coordinate systems are:
C   F
C    L          X - YxZ
C     O         Y - normal to laboratory floor (upwards = +)
C      O        Z - along nominal beam direction (i.e. for PH_B=TH_B=0)
C       R
C   S
C    P          X - along momentum dispersion (=-Y_floor for HRS2-CEBAF)
C     E         Y - ZxX
C      C        Z - along spectrometer central ray
C       T
C
C ---------------------------------------------------------------------
C
      SUBROUTINE TRPT_TO_LAB(TRPT_T,ROTI,PF,QI)
C
      IMPLICIT NONE
      DOUBLE PRECISION  QI(6),TRPT_T(6),ROTI(3,3),PVEC_L(3),PVEC_T(3)
      DOUBLE PRECISION  PF,TTHET,TPHIT,TNORM,CTHEL,STHEL,CPHIL,SPHIL
      DOUBLE PRECISION  POS_T(3)
      INTEGER           I,J,INDI
C
      TTHET = TAN(TRPT_T(2)*0.001D0)       !tan(Theta_Transport)
      TPHIT = TAN(TRPT_T(4)*0.001D0)       !tan(Phi_Transport)
      TNORM = SQRT(1.D0+TTHET**2+TPHIT**2)
C
C ----------------------------------------------------------------------
C      Define the position vector based on the interaction point.
C      Note that this requires a traceback (v3.9).  In version 3.8
C      and before, the position vector was taken to be the
C      intersection of the particle ray with the target
C      Transport (Z=0) plane.
C
C      NOTE:  Here, for simplicity, it is assumed that the drift
C      for the traceback is equal to TRPT_T(5), which is actually
C      the path length difference relative to the target Transport
C      (Z=0) plane.  Note also the sign of the drift!
C ----------------------------------------------------------------------
C
       POS_T(1) =  TRPT_T(1) - TRPT_T(5)*TTHET
       POS_T(2) =  TRPT_T(3) - TRPT_T(5)*TPHIT
       POS_T(3) = -TRPT_T(5)
C
C ---------------------------------------------------------------------
C       Calculate position of interaction point in Laboratory
C       coordinate system in meters.
C ---------------------------------------------------------------------
C
      DO I=1,3                             !Initialize 
        INDI = 2*I-1
        QI(INDI) = 0.D0    
        DO J=1,3
          QI(INDI) = QI(INDI) + 0.01D0*ROTI(I,J)*POS_T(J)
        ENDDO
      ENDDO
C
C ---------------------------------------------------------------------
C       Construct unit momentum vector in Transport system.
C ---------------------------------------------------------------------
C
      PVEC_T(1) = TTHET/TNORM
      PVEC_T(2) = TPHIT/TNORM
      PVEC_T(3) = 1.D0/TNORM
C
C ---------------------------------------------------------------------
C       Rotate to get unit momentum vector in Laboratory system.
C ---------------------------------------------------------------------
C
      DO I=1,3
        PVEC_L(I) = 0.D0                   !initialize
        DO J=1,3
          PVEC_L(I) = PVEC_L(I) + ROTI(I,J)*PVEC_T(J)
        ENDDO
      ENDDO
C
C ---------------------------------------------------------------------
C       Calculate particle angles in Laboratory coordinate system
C       in radians.
C ---------------------------------------------------------------------
C
      STHEL = PVEC_L(2)
      CTHEL = SQRT(1.D0-STHEL**2)
      CPHIL = PVEC_L(3)/CTHEL
      SPHIL = PVEC_L(1)/CTHEL
C
      QI(2) = ASIN(STHEL)
      QI(4) = ACOS(CPHIL)
      IF(SPHIL .LT. 0.D0) QI(4) = -QI(4)
C
C ---------------------------------------------------------------------
C       Calculate momentum in MeV/c.
C ---------------------------------------------------------------------
C
      QI(6) = (1.D0 + 0.01D0*TRPT_T(6))*PF
C
      RETURN
      END
