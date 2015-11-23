C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE MAD_INV
C
C       AUTHOR:  Paul Ulmer
C       DATE:    14-DEC-2004
C       PURPOSE: Calls appropriate MAD Inverse subroutine.
C
C                MAD35_INV:    Not yet implemented
C                MAD12_INV:    Version 3.9 and later
C                MAD12NQ_INV:  Version 3.9 and later
C
C -----------------------------------------------------------------------
C
      SUBROUTINE MAD_INV(IMAD,angle,x_beam,y_beam,x_fp,x_tg)
      IMPLICIT NONE
C      
      DOUBLE PRECISION x_tg(6),x_fp(6),angle,x_beam,y_beam
      INTEGER IMAD
C
ccc      IF(IMAD .EQ. 1) THEN
ccc         CALL MAD35_INV(angle,x_beam,y_beam,x_fp,x_tg)
      IF(IMAD .EQ. 2) THEN
         CALL MAD12_INV(angle,x_beam,y_beam,x_fp,x_tg)
      ELSEIF(IMAD .EQ. 3) THEN
         CALL MAD12NQ_INV(angle,x_beam,y_beam,x_fp,x_tg)
      ENDIF
C
      RETURN
      END

