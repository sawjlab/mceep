C
C       Make gaussian distributed weighting function
C       (normalized to unity).
C
C
      SUBROUTINE GAUSS
      IMPLICIT NONE
      DOUBLE PRECISION PI,X0,SIG,ANORM,X_LO,X_HI,DEL,X,GAUSSIAN
      INTEGER          I,NPTS
C
      CHARACTER*80 OUTFILE
      PARAMETER (PI = 3.1415926D0)
      DATA NPTS /100/
      WRITE(6,10) '$Enter mean and width (sigma) >'
   10 FORMAT(A)
      READ(5,*)X0,SIG
      WRITE(6,10) '$Enter filename > '
      READ(5,20) OUTFILE
   20 FORMAT(A20)
      OPEN(UNIT=1,FILE=OUTFILE,STATUS='NEW',FORM='FORMATTED')
C
      ANORM = 1.D0/(SQRT(2.D0*PI)*SIG)
      X_LO = X0 - 5.D0*SIG
      X_HI = X0 + 5.D0*SIG
      DEL = (X_HI-X_LO)/FLOAT(NPTS)
      DO I = 0,NPTS
        X = FLOAT(I)*DEL + X_LO
        GAUSSIAN = ANORM*EXP(-(X-X0)**2/(2.D0*SIG*SIG))
        WRITE(1,100) X,GAUSSIAN
 100    FORMAT(2X,E12.4,2X,E12.4)
      ENDDO
      CLOSE(UNIT=1)
      RETURN
      END
