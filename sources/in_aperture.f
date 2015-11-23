C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       LOGICAL FUNCTION IN_APERTURE
C
C       AUTHOR: P.E. Ulmer
C       DATE:   7-MAY-1991
C
C       PURPOSE:
C               Test to see whether specified angles lie within
C               solid-angle defining aperture.  Currently, rectangular
C               and elliptical apertures are treated.
C
C------------------------------------------------------------------------------
C
      LOGICAL FUNCTION IN_APERTURE(APER_TYPE,PHI_MIN,THE_MIN,
     #                           PHI_NOM,THE_NOM)
      IMPLICIT NONE
C
      DOUBLE PRECISION PHI_MIN,THE_MIN,PHI_NOM,THE_NOM,TEST_ELL
      CHARACTER*1 APER_TYPE
C
      IN_APERTURE = .TRUE.                     !initialize
C
      IF(APER_TYPE .EQ. 'E') THEN       !elliptical aperture
        TEST_ELL = (PHI_NOM/PHI_MIN)**2 + (THE_NOM/THE_MIN)**2
        IF(TEST_ELL .GT. 1.D0) IN_APERTURE = .FALSE.  !outside of ellipse
C
      ELSE       !rectangular aperture (assumed symmetric)
        IF((PHI_NOM .LT. PHI_MIN) .OR. (PHI_NOM .GT. -PHI_MIN) .OR.
     #     (THE_NOM .LT. THE_MIN) .OR. (THE_NOM .GT. -THE_MIN))
     #            IN_APERTURE = .FALSE.
      ENDIF
C
      RETURN
      END
