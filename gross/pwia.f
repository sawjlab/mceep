C*****PWIA
C
C     This program returns the necessary helicity matrix elements as
C     calculated in the plane wave impluse approximation.
C
      SUBROUTINE PWIA(QL,OMEGAL,W,THETAC,PHIC,IFF,CURME)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMPLEX CURME(0:1,0:1,0:1,-1:1)
C
      DIMENSION SJA(0:1,0:1,0:1,0:1,0:1),SJB(0:1,0:1,0:1,0:1,0:1),
     &          ABA(0:1,0:1,0:1,-1:1),ABB(0:1,0:1,0:1,-1:1)
C
      DOUBLE PRECISION MN,MT
      COMMON/DPARA/MN,MT
C
      PARAMETER(PI=3.141592653589793)
      PARAMETER(ROOT2=1.414213562373095,ROOT2I=1./ROOT2)
C
C                       calculate boost parameters
C
      SHZTAZ=-QL/W
      CHZTAZ=SQRT(SHZTAZ**2+1.)
C
C                       calculate cm momenta and angles
C
      PC=SQRT(0.25*W**2-MN**2)
C      write(6,*)pc
      QC=CHZTAZ*QL+SHZTAZ*OMEGAL
      OMEGAC=CHZTAZ*OMEGAL+SHZTAZ*QL
C
C            calculate momentum and angles for struck nucleon
C
      CSTHC=COS(THETAC)
C
      P1=SQRT(PC**2+QC**2-2.*PC*QC*CSTHC)
      CSTH1=(PC*CSTHC-QC)/P1
      THETA1=ACOS(CSTH1)
C
      P2=SQRT(PC**2+QC**2+2.*PC*QC*CSTHC)
      CSTH2=(PC*CSTHC+QC)/P2
      THETA2=ACOS(CSTH2)
C
C              obtain current operators and boosted wave functions
C
      CALL CURRNT(QC,OMEGAC,PC,THETAC,P1,THETA1,0,0,IFF,1,SJA)
      CALL CURRNT(QC,OMEGAC,PC,THETAC,P2,THETA2,1,1,IFF,1,SJB)
      CALL BOOST(PC,THETAC,SHZTAZ,CHZTAZ,0,1,ABA)
      CALL BOOST(PC,THETAC,SHZTAZ,CHZTAZ,1,1,ABB)
C
C                         construct matrix elements
C
C                 Note: The factor root2i is an isospin factor
C
      DO 600 LDEU=-1,1
        DO 500 LGAM=0,1
          DO 400 LNUC2=0,1
            DO 300 LNUC1=0,1
              CURME(LNUC1,LNUC2,LGAM,LDEU)=(0.,0.)
              DO 200 NUCP=0,1
                DO 100 IENRGY=0,1
                  CURME(LNUC1,LNUC2,LGAM,LDEU)=
     &                 CURME(LNUC1,LNUC2,LGAM,LDEU)
     &                 +ROOT2I*(SJA(LNUC1,NUCP,0,IENRGY,LGAM)
     &                 *ABA(NUCP,LNUC2,IENRGY,-LDEU)
     &                 +SJB(LNUC2,NUCP,0,IENRGY,LGAM)
     &                 *ABB(LNUC1,NUCP,IENRGY,-LDEU))
  100           CONTINUE
  200         CONTINUE
C              write(10,1) lnuc1,lnuc2,lgam,ldeu,
C     &                    curme(lnuc1,lnuc2,lgam,ldeu)
C   1          format(4i3,2(1x,1pe12.5))
  300       CONTINUE
  400     CONTINUE
  500   CONTINUE
  600 CONTINUE
C
      RETURN
      END
