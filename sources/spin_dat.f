C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C        PROGRAM:  SPIN_DAT
C        AUTHOR:   M. Nozar
C        VERSION:
C        DATE:     19-JUL-1991
C         
C        PURPOSE:  This is a sample program consisting of 10 subroutines
C                  each of which returns a 3x3 rotation matrix.  Every
C                  matrix corresponds to an element of a spectrometer.
C                  The elements of each  matrix are functions of the 
C                  transport coordinates.  Note if the Transport matrix
C                  precedes the spin matrix (in the input file) for a
C                  given element then the spin matrix is a function of
C                  the final (i.e. after the element) Transport vector.
C                  However, if the Transport matrix follows the spin matrix
C                  then the spin matrix is a function of the initial
C                  (i.e. before the element) Transport vector.
C                  
C                  (At this time, all the matrices are chosen to be unit
C                  matrices ==> independent of the transport vector.
C                  The user must modify these routines to give actual
C                  spin precession matrices.)
C
C ----------------------------------------------------------------------
C     The spin matrix for the first element
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_1(BEND_0,PMAT_EL)
      IMPLICIT NONE
C
      INCLUDE           'masses.cmn'
CXXX  INCLUDE           'input.cmn'
CXXX  INCLUDE           'spectrometer.cmn'
C
      DOUBLE PRECISION  PMAT_EL(3,3), BEND_0
      DOUBLE PRECISION  G_FACTOR
      INTEGER I,J
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
C
      DO I = 1,3
        DO J = 1,3
          IF(I .EQ. J) THEN
            PMAT_EL(I,J) = 1.D0
          ELSE
            PMAT_EL(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C     The spin matrix for the second element
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_2(BEND_0,PMAT_EL)
      IMPLICIT NONE
C
      INCLUDE           'masses.cmn'
CXXX  INCLUDE           'input.cmn'
CXXX  INCLUDE           'spectrometer.cmn'
C
      DOUBLE PRECISION  PMAT_EL(3,3), BEND_0
      DOUBLE PRECISION  G_FACTOR
      INTEGER I,J
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
C
      DO I = 1,3
        DO J = 1,3
          IF(I .EQ. J) THEN
            PMAT_EL(I,J) = 1.D0
          ELSE
            PMAT_EL(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C     The spin matrix for the third element
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_3(BEND_0,PMAT_EL)
      IMPLICIT NONE
C
      INCLUDE           'masses.cmn'
CXXX  INCLUDE           'input.cmn'
CXXX  INCLUDE           'spectrometer.cmn'
C
      DOUBLE PRECISION  PMAT_EL(3,3), BEND_0
      DOUBLE PRECISION  G_FACTOR
      INTEGER I,J
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
C
      DO I = 1,3
        DO J = 1,3
          IF(I .EQ. J) THEN
            PMAT_EL(I,J) = 1.D0
          ELSE
            PMAT_EL(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C     The spin matrix for the fourth element
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_4(BEND_0,PMAT_EL)
      IMPLICIT NONE
C
      INCLUDE           'masses.cmn'
CXXX  INCLUDE           'input.cmn'
CXXX  INCLUDE           'spectrometer.cmn'
C
      DOUBLE PRECISION  PMAT_EL(3,3), BEND_0
      DOUBLE PRECISION  G_FACTOR
      INTEGER I,J
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
C
      DO I = 1,3
        DO J = 1,3
          IF(I .EQ. J) THEN
            PMAT_EL(I,J) = 1.D0
          ELSE
            PMAT_EL(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C     The spin matrix for the fifth element
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_5(BEND_0,PMAT_EL)
      IMPLICIT NONE
C
      INCLUDE           'masses.cmn'
CXXX  INCLUDE           'input.cmn'
CXXX  INCLUDE           'spectrometer.cmn'
C
      DOUBLE PRECISION  PMAT_EL(3,3), BEND_0
      DOUBLE PRECISION  G_FACTOR
      INTEGER I,J
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
C
      DO I = 1,3
        DO J = 1,3
          IF(I .EQ. J) THEN
            PMAT_EL(I,J) = 1.D0
          ELSE
            PMAT_EL(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C     The spin matrix for the sixth element
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_6(BEND_0,PMAT_EL)
      IMPLICIT NONE
C
      INCLUDE           'masses.cmn'
CXXX  INCLUDE           'input.cmn'
CXXX  INCLUDE           'spectrometer.cmn'
C
      DOUBLE PRECISION  PMAT_EL(3,3), BEND_0
      DOUBLE PRECISION  G_FACTOR
      INTEGER I,J
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
C
      DO I = 1,3
        DO J = 1,3
          IF(I .EQ. J) THEN
            PMAT_EL(I,J) = 1.D0
          ELSE
            PMAT_EL(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C     The spin matrix for the seventh element
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_7(BEND_0,PMAT_EL)
      IMPLICIT NONE
C
      INCLUDE           'masses.cmn'
CXXX  INCLUDE           'input.cmn'
CXXX  INCLUDE           'spectrometer.cmn'
C
      DOUBLE PRECISION  PMAT_EL(3,3), BEND_0
      DOUBLE PRECISION  G_FACTOR
      INTEGER I,J
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
C
      DO I = 1,3
        DO J = 1,3
          IF(I .EQ. J) THEN
            PMAT_EL(I,J) = 1.D0
          ELSE
            PMAT_EL(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C     The spin matrix for the eighth element
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_8(BEND_0,PMAT_EL)
      IMPLICIT NONE
C
      INCLUDE           'masses.cmn'
CXXX  INCLUDE           'input.cmn'
CXXX  INCLUDE           'spectrometer.cmn'
C
      DOUBLE PRECISION  PMAT_EL(3,3), BEND_0
      DOUBLE PRECISION  G_FACTOR
      INTEGER I,J
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
C
      DO I = 1,3
        DO J = 1,3
          IF(I .EQ. J) THEN
            PMAT_EL(I,J) = 1.D0
          ELSE
            PMAT_EL(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C     The spin matrix for the ninth element
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_9(BEND_0,PMAT_EL)
      IMPLICIT NONE
C
      INCLUDE           'masses.cmn'
CXXX  INCLUDE           'input.cmn'
CXXX  INCLUDE           'spectrometer.cmn'
C
      DOUBLE PRECISION  PMAT_EL(3,3), BEND_0
      DOUBLE PRECISION  G_FACTOR
      INTEGER I,J
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
C
      DO I = 1,3
        DO J = 1,3
          IF(I .EQ. J) THEN
            PMAT_EL(I,J) = 1.D0
          ELSE
            PMAT_EL(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C     The spin matrix for the tenth element
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_10(BEND_0,PMAT_EL)
      IMPLICIT NONE
C
      INCLUDE           'masses.cmn'
CXXX  INCLUDE           'input.cmn'
CXXX  INCLUDE           'spectrometer.cmn'
C
      DOUBLE PRECISION  PMAT_EL(3,3), BEND_0
      DOUBLE PRECISION  G_FACTOR
      INTEGER I,J
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
C
      DO I = 1,3
        DO J = 1,3
          IF(I .EQ. J) THEN
            PMAT_EL(I,J) = 1.D0
          ELSE
            PMAT_EL(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
