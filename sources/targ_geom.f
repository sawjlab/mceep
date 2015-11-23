C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       SUBROUTINE ACT_TO_NOM
C
C       AUTHOR: P.E. Ulmer
C       DATE:   7-MAY-1991
C
C       PURPOSE:
C               Calculate "nominal" angles from "actual" angles in the
C               "LAB" coordinate system.  These angles also depend on the
C               interaction point within the extended target and the
C               "nominal" drift distance to solid-angle defining aperture.
C
C           The origin (Pt. O) of the "LAB" coordinate system is taken
C           to be the nominal center of the target.  The interaction
C           point (Pt. I) then defines a vector from Pt. O.  The
C           solid-angle defining aperture is assumed to be a section of
C           a sphere centered at Pt. O of radius DRIFTN (nominal drift
C           from Pt. O to aperture).  The intersection of the surface of
C           this sphere with a line drawn from Pt. I and oriented at
C           the "actual" laboratory angles (THE_ACT, PHI_ACT) determines
C           the entry point of the particle at the solid-angle defining
C           aperture (Pt. P).  The angles of the vector from Pt. O
C           to Pt. P in the "LAB" system are the "nominal" angles.
C           Tests can be applied to these angles to determine whether
C           they reside within the solid-angle defining aperture.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE ACT_TO_NOM(PHI_ACT,THE_ACT,DRIFTN,X_I,Y_I,Z_I,
     #                      PHI_NOM,THE_NOM)
      IMPLICIT NONE
C
      DOUBLE PRECISION PHI_ACT,THE_ACT,PHI_NOM,THE_NOM
      DOUBLE PRECISION CPHI_ACT,CTHE_ACT,SPHI_ACT,STHE_ACT
      DOUBLE PRECISION X_I,Y_I,Z_I
      DOUBLE PRECISION DRIFTN,DRIFTA,R,D
C
C------------------------------------------------------------------------------
C     Get sines and cosines of "actual" angles.
C------------------------------------------------------------------------------
C
      CPHI_ACT = COS(PHI_ACT)
      CTHE_ACT = COS(THE_ACT)
      SPHI_ACT = SIN(PHI_ACT)
      STHE_ACT = SIN(THE_ACT)
C
C------------------------------------------------------------------------------
C     Compute actual drift to collimator from the interaction point
C     (X_I, Y_I, Z_I), the "actual" angles and the nominal drift.
C     The "actual drift" is defined as the distance from the interaction
C     point to Pt. P (see introductory comments).
C------------------------------------------------------------------------------
C
      R = X_I*CTHE_ACT*SPHI_ACT + Y_I*STHE_ACT + Z_I*CTHE_ACT*CPHI_ACT
      D = DRIFTN**2 - (X_I**2 + Y_I**2 + Z_I**2)
      DRIFTA = SQRT(D + R**2) - R            !Actual drift
C
C------------------------------------------------------------------------------
C     Compute nominal angles (see introductory comments).
C------------------------------------------------------------------------------
C
      THE_NOM = ASIN((DRIFTA*STHE_ACT + Y_I)/DRIFTN)
      PHI_NOM = ACOS((DRIFTA*CTHE_ACT*CPHI_ACT + Z_I)/
     #               (DRIFTN*COS(THE_NOM)))
      IF(PHI_ACT .LT. 0.D0) PHI_NOM = - PHI_NOM
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       SUBROUTINE NOM_TO_ACT
C
C       AUTHOR: P.E. Ulmer
C       DATE:   7-MAY-1991
C
C       MODIFICATIONS:
C
C               22-AUG-1999  P.E. Ulmer
C               Calculate solid angle Jacobian.  This was omitted
C               previously and this can lead to significant error
C               for long targets.  (Thanks to Jeff Templon for
C               pointing out the omission.)
C
C       PURPOSE:
C               Calculate "actual" angles from "nominal" angles,
C               interaction point within extended target and "nominal"
C               drift distance to solid-angle defining aperture.
C               For a point target, the nominal angles and actual angles
C               are the same.
C
C               Also, get the solid angle Jacobian between the
C               actual angles and the nominal angles.  This is necessary
C               since the nominal solid angle is uniformly sampled,
C               but we need uniform sampling in the actual solid angle
C               (since the cross section is differential in the latter).
C
C           The origin (Pt. O) of the "LAB" coordinate system is taken
C           to be the nominal center of the target.  The interaction
C           point (Pt. I) then defines a vector from Pt. O.  The
C           solid-angle defining aperture is assumed to be a section of
C           a sphere centered at Pt. O of radius DRIFTN (nominal drift
C           from Pt. O to aperture).  "Nominal" angles (which have been
C           chosen somewhere within the aperture angular ranges) are
C           defined relative to Pt. O.  The intersection of the ray which
C           starts at Pt. O and is oriented at the nominal angles with
C           the section of the sphere subtended by the aperture defines
C           Pt. P.  The angles of the vector from Pt. I to Pt. P represent
C           the "actual" laboratory angles (THE_ACT, PHI_ACT).  The
C           sampling of so-called "nominal" angles facilitates simulation
C           of experiments employing an extended interaction region.
C           This has the advantage that the sampled points are guaranteed
C           to lie within the aperture.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE NOM_TO_ACT(PHI_NOM,THE_NOM,DRIFTN,X_I,Y_I,Z_I,
     #                      PHI_ACT,THE_ACT,DOMEGA_DOMEGN)
      IMPLICIT NONE
C
      DOUBLE PRECISION PHI_ACT,THE_ACT,PHI_NOM,THE_NOM
      DOUBLE PRECISION CPHI_NOM,CTHE_NOM,SPHI_NOM,STHE_NOM
      DOUBLE PRECISION CPHI_ACT,CTHE_ACT,SPHI_ACT,STHE_ACT
      DOUBLE PRECISION X_I,Y_I,Z_I,INT_TO_COL(3),DOT
      DOUBLE PRECISION DRIFTN,DRIFTA,TEST_CPHI_ACT
      DOUBLE PRECISION DRX_DTHEN,DRY_DTHEN,DRZ_DTHEN
      DOUBLE PRECISION DRX_DPHIN,DRY_DPHIN,DRZ_DPHIN
      DOUBLE PRECISION DDRIFTA_DTHEN,DDRIFTA_DPHIN
      DOUBLE PRECISION DTHEA_DTHEN,DPHIA_DTHEN
      DOUBLE PRECISION DTHEA_DPHIN,DPHIA_DPHIN
      DOUBLE PRECISION DOMEGA_DOMEGN
C
C------------------------------------------------------------------------------
C     Get sines and cosines of nominal angles.
C------------------------------------------------------------------------------
C
      CPHI_NOM = COS(PHI_NOM)
      CTHE_NOM = COS(THE_NOM)
      SPHI_NOM = SIN(PHI_NOM)
      STHE_NOM = SIN(THE_NOM)
C
C ---------------------------------------------------------------------
C     Determine (X,Y,Z) position of ray at entrance to
C     collimator relative to interaction point.
C ---------------------------------------------------------------------
C
      INT_TO_COL(1) = DRIFTN*CTHE_NOM*SPHI_NOM - X_I
      INT_TO_COL(2) = DRIFTN*STHE_NOM          - Y_I
      INT_TO_COL(3) = DRIFTN*CTHE_NOM*CPHI_NOM - Z_I
C
C ---------------------------------------------------------------------
C     Determine actual drift distance, defined as the distance from
C     the interaction point to the solid-angle defining aperture.
C ---------------------------------------------------------------------
C
      DRIFTA = SQRT(DOT(INT_TO_COL,INT_TO_COL))  !Actual drift
C
C ---------------------------------------------------------------------
C     Determine "actual" angles of particle in LAB system.
C ---------------------------------------------------------------------
C
      STHE_ACT = INT_TO_COL(2)/DRIFTA
      THE_ACT = ASIN(STHE_ACT)
      CTHE_ACT = COS(THE_ACT)
      CPHI_ACT = INT_TO_COL(3)/(DRIFTA*CTHE_ACT)
      TEST_CPHI_ACT = 1.D0 - CPHI_ACT
      IF(TEST_CPHI_ACT .LT. 0.D0 .AND. TEST_CPHI_ACT .GT. -1.D-12)
     #     CPHI_ACT =  1.D0   !protect against roundoff
      TEST_CPHI_ACT = 1.D0 + CPHI_ACT
      IF(TEST_CPHI_ACT .LT. 0.D0 .AND. TEST_CPHI_ACT .GT. -1.D-12)
     #     CPHI_ACT = -1.D0   !protect against roundoff
      PHI_ACT = ACOS(CPHI_ACT)
      IF(INT_TO_COL(1) .LT. 0.) PHI_ACT = - PHI_ACT
      SPHI_ACT = SIN(PHI_ACT)
C
C ---------------------------------------------------------------------
C     Now determine the solid angle Jacobian between the actual angles
C     and the nominal angles.
C
C     First, get the derivatives of the INT_TO_COL vector with respect
C     to the nominal angles.  (As a shorthand, INT_TO_COL(1) is
C     referred to as RX here, and likewise for other components.)
C ---------------------------------------------------------------------
C
      DRX_DTHEN = -DRIFTN*STHE_NOM*SPHI_NOM
      DRX_DPHIN =  DRIFTN*CTHE_NOM*CPHI_NOM
C
      DRY_DTHEN =  DRIFTN*CTHE_NOM
      DRY_DPHIN =  0.D0
C
      DRZ_DTHEN = -DRIFTN*STHE_NOM*CPHI_NOM
      DRZ_DPHIN = -DRIFTN*CTHE_NOM*SPHI_NOM
C
      DDRIFTA_DTHEN = ( INT_TO_COL(1)*DRX_DTHEN
     #                + INT_TO_COL(2)*DRY_DTHEN
     #                + INT_TO_COL(3)*DRZ_DTHEN ) / DRIFTA
C
      DDRIFTA_DPHIN = ( INT_TO_COL(1)*DRX_DPHIN
     #                + INT_TO_COL(2)*DRY_DPHIN
     #                + INT_TO_COL(3)*DRZ_DPHIN ) / DRIFTA
C
      DTHEA_DTHEN = (DRY_DTHEN-INT_TO_COL(2)*DDRIFTA_DTHEN/DRIFTA)
     #             / (DRIFTA*CTHE_ACT)
      DTHEA_DPHIN = (DRY_DPHIN-INT_TO_COL(2)*DDRIFTA_DPHIN/DRIFTA)
     #             / (DRIFTA*CTHE_ACT)
C
      DPHIA_DTHEN = - ( STHE_ACT*CPHI_ACT*DTHEA_DTHEN
     #              +(DRZ_DTHEN-INT_TO_COL(3)*DDRIFTA_DTHEN/DRIFTA)
     #               /DRIFTA ) / ( CTHE_ACT*SPHI_ACT )
      DPHIA_DPHIN = - ( STHE_ACT*CPHI_ACT*DTHEA_DPHIN
     #              +(DRZ_DPHIN-INT_TO_COL(3)*DDRIFTA_DPHIN/DRIFTA)
     #               /DRIFTA ) / ( CTHE_ACT*SPHI_ACT )
C
C ---------------------------------------------------------------------
C     The Jacobian below had an error pointed out by Marat Rvachev.
C     The last factor was previously DPHIA_DPHIN but has been
C     changed to DPHIA_DTHEN.  Apparently, ordinarily the error
C     was small (less than 0.3% for the cases tested by Marat).  This
C     is not surprising since the error was in the "cross terms".  The
C     fix was made on 3/13/01.
C ---------------------------------------------------------------------
C
      DOMEGA_DOMEGN = ABS( (CTHE_ACT/CTHE_NOM) *
     #             (DTHEA_DTHEN*DPHIA_DPHIN-DTHEA_DPHIN*DPHIA_DTHEN) )
C
      RETURN
      END
