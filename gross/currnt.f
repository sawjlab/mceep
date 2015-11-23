C
C_________________________________________________________________________
C
C*****currnt
C
C     This subroutine returns the matrix elements of the relativistic single
C     nucleon current operator between positive and negative energy Dirac
C     spinors.
C
C     INPUT
C           qc...........momentum transfer in frame of calculation
C           omegac.......energy transfer in frame of calculation
C           pp...........momentum of final nucleon
C           thetap.......angle of final nucleon
C           p............momentum of initial nucleon
C           theta........angle of initial nucleon
C           itype=0......particle 2 on shell
C                 1......particle 1 on shell
C           isosp=0......proton
C                 1......neutron
C                 2......isoscalar
C                 3......isovector
C           iff=1,6......code to choose type of form factor in DIRCFF
C           irel=1.......relativistic current operator
C               =0.......nonrelativistic current operator
C
C     OUTPUT
C           sj(lnucp,lnuc,inrgyp,inrgy,lgam)
C
C               lnucp=0..............positive nucleon helicity
C                    =1..............negative nucleon helicity
C               lnuc=0...............positive nucleon helicity
C                   =1...............negative nucleon helicity
C               inrgyp=0.............positive energy projection
C                     =1.............negative energy projection
C               inrgy=0..............positive energy projection
C                    =1..............negative energy projection
C               lgam=0,1.............photon helicity
C
C
      SUBROUTINE CURRNT(QC,OMEGAC,PP,THETAP,P,THETA,ITYPE,
     &                 ISOSP,IFF,IREL,SJ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION SJ(0:1,0:1,0:1,0:1,0:1),F(0:1,2),DEL1(0:1,0:1),
     &          DEL2(0:1,0:1),DEL3(0:1,0:1)
C
      DOUBLE PRECISION MN,MT
      COMMON/DPARA/MN,MT
C
C               calculate square of four-momentum transfer
C
      Q2=QC**2-OMEGAC**2
      Q=SQRT(Q2)
      Q2GEV=Q2*1.E-6
C
C                     call form factor routine
C
      CALL DIRCFF(Q2GEV,F,IFF)
C
      IF(ISOSP.EQ.0.OR.ISOSP.EQ.1) THEN    ! proton or neutron
        F1=F(ISOSP,1)
        F2=F(ISOSP,2)
      ELSE IF(ISOSP.EQ.3) THEN             ! isoscalar
        F1=F(0,1)+F(1,1)
        F2=F(0,2)+F(1,2)
      ELSE IF(ISOSP.EQ.4) THEN             ! isovector
        F1=F(0,1)-F(1,1)
        F2=F(0,2)-F(1,2)
      END IF
C
      FAC1=F1*QC/Q
      FAC2=F1*OMEGAC/Q
      FAC3=0.5*F2*Q/MN
      FAC4=F1
      FAC5=0.5*F2*OMEGAC/MN
      FAC6=0.5*F2*QC/MN
C
      CALL ANGFAC(THETAP,THETA,DEL1,DEL2,DEL3)
C
      IF(IREL.EQ.1) THEN                ! relativistic current operator
C
        EPP=SQRT(PP**2+MN**2)
        EP=SQRT(P**2+MN**2)
        FAC=SQRT((EP+MN)*(EPP+MN))/(2.*MN)
        ETAPP=PP/(EPP+MN)
        ETAP=P/(EP+MN)
C
        DO LNUC=0,1
          GP=(1.-2.*LNUC)*ETAP
          IP=ABS(ITYPE-LNUC)
          DO LNUCP=0,1
            GPP=(1.-2.*LNUCP)*ETAPP
            IPP=ABS(ITYPE-LNUCP)
C
            GAM1=1.+GPP*GP
            GAM2=1.-GPP*GP
            GAM3=GPP+GP
            GAM4=GPP-GP
C
            SJ(LNUCP,LNUC,0,0,0)=FAC*(FAC1*GAM1*DEL1(IPP,IP)-(FAC2*GAM3
     &                         +FAC3*GAM4)*DEL2(IPP,IP))
            SJ(LNUCP,LNUC,0,1,0)=FAC*(FAC1*GAM3*DEL1(IPP,IP)-(FAC2*GAM1
     &                         -FAC3*GAM2)*DEL2(IPP,IP))
            SJ(LNUCP,LNUC,1,0,0)=FAC*(FAC1*GAM3*DEL1(IPP,IP)-(FAC2*GAM1
     &                         +FAC3*GAM2)*DEL2(IPP,IP))
            SJ(LNUCP,LNUC,1,1,0)=FAC*(FAC1*GAM1*DEL1(IPP,IP)-(FAC2*GAM3
     &                         -FAC3*GAM4)*DEL2(IPP,IP))
            SJ(LNUCP,LNUC,0,0,1)=-FAC*(FAC4*GAM3-FAC5*GAM4
     &                             +FAC6*GAM2)*DEL3(IPP,IP)
            SJ(LNUCP,LNUC,0,1,1)=-FAC*(FAC4*GAM1+FAC5*GAM2
     &                             -FAC6*GAM4)*DEL3(IPP,IP)
            SJ(LNUCP,LNUC,1,0,1)=-FAC*(FAC4*GAM1-FAC5*GAM2
     &                             +FAC6*GAM4)*DEL3(IPP,IP)
            SJ(LNUCP,LNUC,1,1,1)=-FAC*(FAC4*GAM3+FAC5*GAM4
     &                             -FAC6*GAM2)*DEL3(IPP,IP)
          END DO
        END DO
C
      ELSE IF(IREL.EQ.0) THEN         ! nonrelativistic current operator
C
        FAC=SQRT((EP+MN)*(EPP+MN))/(2.*MN)
        ETAPP=PP/(2.*MN)
        ETAP=P/(2.*MN)
C
        DO LNUC=0,1
          GP=(1.-2.*LNUC)*ETAP
          IP=ABS(ITYPE-LNUC)
          DO LNUCP=0,1
            GPP=(1.-2.*LNUCP)*ETAPP
            IPP=ABS(ITYPE-LNUCP)
C
            GAM3=GPP+GP
C
            SJ(LNUCP,LNUC,0,0,0)=FAC1*DEL1(IPP,IP)
            SJ(LNUCP,LNUC,0,1,0)=0.0
            SJ(LNUCP,LNUC,1,0,0)=0.0
            SJ(LNUCP,LNUC,1,1,0)=0.0
            SJ(LNUCP,LNUC,0,0,1)=-(FAC4*GAM3+FAC6)*DEL3(IPP,IP)
            SJ(LNUCP,LNUC,0,1,1)=0.0
            SJ(LNUCP,LNUC,1,0,1)=0.0
            SJ(LNUCP,LNUC,1,1,1)=0.0
          END DO
        END DO
C
      ENDIF
C
      RETURN
      END
