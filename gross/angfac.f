C
C_________________________________________________________________________
C
C*****angfac
C
C     This subroutine calculates the angular dependent factors needed in
C     CURRNT.
C
C     INPUT
C          thetap..................angle of nucleon in final state
C          theta...................angle of nucleon in initial state
C
C     OUTPUT
C          del1(lp,l)..............angular factor 1
C          del2(lp,l)..............angular factor 2
C          del3(lp,l)..............angular factor 3
C
C                   lp=0...........positive helicity
C                     =1...........negative helicity
C                   l=0............positive helicity
C                    =1............negative helicity
C
      SUBROUTINE ANGFAC(THETAP,THETA,DEL1,DEL2,DEL3)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION DEL1(0:1,0:1),DEL2(0:1,0:1),DEL3(0:1,0:1),
     &          DHALFP(0:1,0:1),DHALF(0:1,0:1)
C
      DATA ROOT2/1.414213562373095/
C
      CALL ROTHLF(THETAP,DHALFP)
      CALL ROTHLF(THETA,DHALF)
C
      DO 300 L=0,1
        DO 200 LP=0,1
          DEL1(LP,L)=0.0
          DEL2(LP,L)=0.0
          DO 100 I=0,1
            DEL1(LP,L)=DEL1(LP,L)+DHALFP(I,LP)*DHALF(I,L)
            DEL2(LP,L)=DEL2(LP,L)+(1.-2.*I)*DHALFP(I,LP)*DHALF(I,L)
  100     CONTINUE
          DEL3(LP,L)=-ROOT2*DHALFP(0,LP)*DHALF(1,L)
  200   CONTINUE
  300 CONTINUE
C
      RETURN
      END
