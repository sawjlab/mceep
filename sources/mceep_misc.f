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
      CHARACTER*300 CHARIN,CHAROUT
      N = 0
      DO J=1,300
        IF(CHARIN(J:J).NE.' ') THEN
          N = N+1
          CHAROUT(N:N) = CHARIN(J:J)
        ENDIF
      ENDDO
C
      DO J=N+1,300
         CHAROUT(J:J) = ' '
      ENDDO
C
      NCHAR = N
      RETURN
      END
C
C------------------------------------------------------------------------------
C       LOGICAL FUNCTION ANSWER
C
C       AUTHOR: P.E. Ulmer
C       DATE:   31-AUG-1990
C       PURPOSE:
C               Sets logical flag depending upon the answer to a question.
C------------------------------------------------------------------------------
C
      LOGICAL FUNCTION ANSWER(QUESTION,ANS_1,ANS_2)
      IMPLICIT NONE
C
      CHARACTER*(*) QUESTION
      CHARACTER*1 ANS_1,ANS_2,ANS
C
      WRITE(6,10) QUESTION
   10 FORMAT(A)
      READ(5,20)ANS
   20 FORMAT(A1)
      IF(ANS .EQ. ANS_1 .OR. ANS .EQ. ANS_2) THEN
        ANSWER = .TRUE.
      ELSE
        ANSWER = .FALSE.
      ENDIF
      RETURN
      END
C
C------------------------------------------------------------------------------
C       Subroutine GET_PID_DEN
C
C       AUTHOR: P.E. Ulmer
C       DATE:   V1.0    29-AUG-1990
C       MODIFICATIONS:
C
C       PURPOSE:
C               Gets particle ID types for spectrum which divides
C               yield histos.  If the desired histogram is average
C               cross section then IPID_DEN1 = 0 (i.e. divide
C               yield histos by phase space).  If the desired histogram
C               is average polarization or e- analyzing power
C               then IPID_DEN1 = 1 (i.e. divide partial yield histos by
C               unpolarized yield array).
C
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_PID_DEN(IPID1,IPID_DEN1)
      IMPLICIT NONE
C
      INTEGER IPID1,IPID_DEN1,ISIG_AVG(19),IPOL_AVG(28),NSIG,NPOL,J
C
      DATA ISIG_AVG/-1,-50,-51,-52,-53,-100,-101,-102,-103,
     #              -121,-122,-123,-124,-125,-126,-127,-128,-129,-130/
                                              !PIDs for cross section histos
      DATA IPOL_AVG/-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,
     #             -15,-16,-17,-18,-19,-20,-22,-23,-24,-25,-26,-27,
     #             -28,-29,-30/
                                              !PIDs for pol/asymm  histos
      DATA NSIG,NPOL/19,28/
                                              !Number of each type of histo
C
      IF(IPID1 .LT. 0)THEN            !Do for avgd. sigma and pol/asym only
        DO J=1,NSIG
          IF(IPID1 .EQ. ISIG_AVG(J)) IPID_DEN1 = 0
        ENDDO
        DO J=1,NPOL
          IF(IPID1 .EQ. IPOL_AVG(J)) IPID_DEN1 = 1
        ENDDO
      ENDIF
C
      RETURN
      END
