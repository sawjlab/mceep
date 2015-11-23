C
C ---------------------------------------------------------------------
C
C       PROGRAM MCEEP_UTIL
C
C       AUTHOR:  P.E. Ulmer
C       DATE:    10-MAR-1993
C
C       Companion programs for MCEEP.
C
C ---------------------------------------------------------------------
C
      WRITE(6,5)
    5 FORMAT(//'      MCEEP  UTILITIES  PROGRAM ')
   10 WRITE(6,1000)
 1000 FORMAT(/
     #   /'  MCEEP_TO_HBOOK ............................ 1 '
     #   /'  RATIO ..................................... 2 '
     #   /'  GAUSS FILE ................................ 3 '
     #   /'  ADD R-FUNCTIONS TO ROW-WISE N-TUPLE ....... 4 '
     #   /'  GET INCLUDE FILE FOR COLUMN-WISE N-TUPLE .. 5 '
     #   /'  ADD R-FUNCTIONS TO COLUMN-WISE N-TUPLE .... 6 '
     #  //'  Enter Option>')
      READ(5,1010,ERR=10,END=999)I
 1010 FORMAT(I1)
C
      GOTO(100,200,300,400,500,600)I
      GOTO 999
C
  100 CALL MCEEP_TO_HBOOK
      GOTO 10
  200 CALL RATIO
      GOTO 10
  300 CALL GAUSS
      GOTO 10
  400 CALL MAKE_NEW_RWN
      GOTO 10
  500 CALL GET_INCLUDE_CWN
      GOTO 10
  600 CALL MAKE_NEW_CWN
      GOTO 10
C
  999 STOP
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       Extracts non-blank characters of any character string.
C------------------------------------------------------------------------------
C
      SUBROUTINE SQUEEZE(CHARIN,CHAROUT,NCHAR)
      IMPLICIT NONE
C
      INTEGER NCHAR,N,J
      CHARACTER*80 CHARIN,CHAROUT
      N = 0
      DO J=1,80
        IF(CHARIN(J:J).NE.' ') THEN
          N = N+1
          CHAROUT(N:N) = CHARIN(J:J)
        ENDIF
      ENDDO
      NCHAR = N
      RETURN
      END
