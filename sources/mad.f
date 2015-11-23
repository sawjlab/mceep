C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE MAD
C
C       AUTHOR:  Paul Ulmer
C       DATE:    14-DEC-2004
C       PURPOSE: Calls appropriate MAD subroutine.
C
C                MAD35:    Version 3.7 and later
C                MAD12:    Version 3.9 and later
C                MAD12NQ:  Version 3.9 and later
C
C -----------------------------------------------------------------------
C
      SUBROUTINE MAD(IMAD,x_tg,x_fp,fail_apertures)
      IMPLICIT NONE
C      
      DOUBLE PRECISION x_tg(6),x_fp(6)
      INTEGER IMAD
      LOGICAL fail_apertures
C
      IF(IMAD .EQ. 1) THEN
         CALL MAD35(x_tg,x_fp,fail_apertures)
      ELSEIF(IMAD .EQ. 2) THEN
         CALL MAD12(x_tg,x_fp,fail_apertures)
      ELSEIF(IMAD .EQ. 3) THEN
         CALL MAD12NQ(x_tg,x_fp,fail_apertures)
      ENDIF
C
      RETURN
      END

