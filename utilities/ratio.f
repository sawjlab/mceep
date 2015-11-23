C
C ----------------------------------------------------------------------
C
C       Subroutine RATIO
C
C       AUTHOR:  P.E. Ulmer
C       DATE:    10-MAR-1993
C
C       Reads two 1-D MCEEP generated Topdrawer files and produces
C       either a ratio or difference histrogram.  For ratios, it
C       also produces an average for all but NELIM channels
C       on each side.
C
C ----------------------------------------------------------------------
C
      SUBROUTINE RATIO
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(1000),RAT(1000)
      CHARACTER*80 FILE1,FILE2,OUTFILE,DUMMY
      CHARACTER*1 ANS
      LOGICAL RATFLAG
C
      DATA NELIM /5/    !Number of channels on each side to eliminate
                        !for ratio histogram
C
      WRITE(6,10) '$Enter file1 >'
   10 FORMAT(A)
      READ(5,20) FILE1
   20 FORMAT(A20)
      WRITE(6,10) '$Enter file2 >'
      READ(5,20) FILE2
      WRITE(6,10) '$Enter output file >'
      READ(5,20) OUTFILE
      WRITE(6,10) '$Ratio (R) or Difference (D)? (R/D) >'
      READ(5,30) ANS
   30 FORMAT(A1)
      RATFLAG = .TRUE.
      IF((ANS .EQ. 'D') .OR. (ANS .EQ. 'd')) RATFLAG = .FALSE.
C
      OPEN(UNIT=1,FILE=FILE1,STATUS='OLD',FORM='FORMATTED')
      OPEN(UNIT=2,FILE=FILE2,STATUS='OLD',FORM='FORMATTED')
C
      OPEN(UNIT=3,FILE=OUTFILE,STATUS='NEW',FORM='FORMATTED')
      WRITE(3,'(''   set size 7.88 by 10.5 '')')
      WRITE(3,'(''   set intensity 1 '')')
      IF(RATFLAG) WRITE(3,
     #  '(''   title left '''' ratio-1 (percent) '' '' '''' '')')
C
      DO I=1,7
        READ(1,20) DUMMY
        READ(2,20) DUMMY
      ENDDO
C
      NTERMS = 0
C
      DO I=1,1000
        READ(1,*,ERR=999,END=999) X(I),Y1
        READ(2,*) X(I),Y2
        IF(RATFLAG) THEN
          IF(Y2 .NE. 0.D0) THEN
            NTERMS = NTERMS + 1
            RAT(NTERMS) = Y1/Y2
            WRITE(3,100) X(I),(RAT(NTERMS)-1.D0)*100.D0
  100       FORMAT(3X,2E12.4)
          ENDIF
        ELSE
          DIFF = Y1-Y2
          WRITE(3,100) X(I),DIFF
        ENDIF
      ENDDO
C
  999 RAT_SUM = 0.D0
      RAT_SUM_SQ = 0.D0
      IF(RATFLAG) THEN
        DO I=1+NELIM,NTERMS-NELIM  !Eliminate NELIM channels on each side
          RAT_SUM = RAT_SUM + RAT(I)
          RAT_SUM_SQ = RAT_SUM_SQ + RAT(I)**2
        ENDDO
      ENDIF
C
      CLOSE(UNIT=1)
      CLOSE(UNIT=2)
C
      WRITE(3,'(''   histogram solid '')')
      IF(RATFLAG)THEN
        WRITE(3,'(''   set order x y '')')
        WRITE(3,50) X(1),0.D0
   50   FORMAT(3X,F12.4,3X,F12.4)
        WRITE(3,50) X(NTERMS),0.D0
        WRITE(3,'(''   join dash '')')
        RAT_SUM = RAT_SUM/DFLOAT(NTERMS-2*NELIM)
        SIGMA = SQRT(RAT_SUM_SQ/DFLOAT(NTERMS-2*NELIM)-RAT_SUM**2)
        WRITE(6,110) RAT_SUM,SIGMA
  110   FORMAT(/' Avg. ratio  = ',E12.4,/,' Std. dev. = ',E12.4,//)
        WRITE(3,
     #  '(''   title 2.5 4.0 '''' Avg. ratio: '',T35,E12.4,'' '''' '')')
     #        RAT_SUM
        WRITE(3,
     #  '(''   title 2.5 3.5 '''' Std. dev.:  '',T35,E12.4,'' '''' '')')
     #        SIGMA
      ENDIF
C
      CLOSE(UNIT=3)
C
      RETURN
      END
