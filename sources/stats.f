C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       Obtain statistics for a given 1-D spectrum passed in as array ARR.
C------------------------------------------------------------------------------
C
      SUBROUTINE SPECSTATS_1D(XMIN,XMAX,NBINSX,ARR,CENT_X,SUM)
      IMPLICIT NONE
C
      DOUBLE PRECISION ARR(500),X,XMIN,XMAX,CENT_X,SUM,SUM_X,DELTAX
      INTEGER NBINSX,I
C
      SUM = 0.                !Initialize sums
      SUM_X = 0.
C
      DELTAX = (XMAX-XMIN)/DFLOAT(NBINSX)     !bin width
      DO I=1,NBINSX
        X = (DFLOAT(I)-0.5D0)*DELTAX + XMIN
        SUM_X = SUM_X + X*ARR(I)
        SUM = SUM + ARR(I)
      ENDDO
C
      IF(SUM .NE. 0.) THEN
        CENT_X = SUM_X/SUM             !Centroid
      ELSE
        CENT_X = 0.
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       Obtain statistics for a given 2-D spectrum passed in as array ARR_2D.
C------------------------------------------------------------------------------
C
      SUBROUTINE SPECSTATS_2D(XMIN,XMAX,YMIN,YMAX,NBINSX,NBINSY,ARR_2D,
     #                        CENT_X,CENT_Y,SUM)
      IMPLICIT NONE
C
      DOUBLE PRECISION ARR_2D(50,50),X,XMIN,XMAX,Y,YMIN,YMAX
      DOUBLE PRECISION CENT_X,CENT_Y,SUM,SUM_X,SUM_Y,DELTAX,DELTAY
      INTEGER NBINSX,NBINSY,I,J
C
      SUM = 0.                !Initialize sums
      SUM_X = 0.
      SUM_Y = 0.
C
      DELTAX = (XMAX-XMIN)/DFLOAT(NBINSX)     !bin width in X
      DELTAY = (YMAX-YMIN)/DFLOAT(NBINSY)     !bin width in Y
      DO I=1,NBINSX
      DO J=1,NBINSY
        X = (DFLOAT(I)-0.5D0)*DELTAX + XMIN
        Y = (DFLOAT(J)-0.5D0)*DELTAY + YMIN
        SUM_X = SUM_X + X*ARR_2D(I,J)
        SUM_Y = SUM_Y + Y*ARR_2D(I,J)
        SUM = SUM + ARR_2D(I,J)
      ENDDO
      ENDDO
C
      IF(SUM .NE. 0.) THEN
        CENT_X = SUM_X/SUM             !Centroid in X
        CENT_Y = SUM_Y/SUM             !Centroid in Y
      ELSE
        CENT_X = 0.
        CENT_Y = 0.
      ENDIF
C
      RETURN
      END
