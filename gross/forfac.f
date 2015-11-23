C
C_______________________________________________________________________
C
C*****sachs_gross
C
      SUBROUTINE SACHS_GROSS(Q2,G,IFF)
C
C     This subroutine returns Sachs form factors for a variety of
C     parameterizations.
C
C     9-NOV-1999 (PEU)
C     Changed name of subroutine from SACHS to SACHS_GROSS to avoid
C     conflict with J.J. Kelly's EPIPROD routines, soon to be
C     incorporated into MCEEP.
C
C
C     INPUT
C             Q2.............negative of square of four-momentum transfer
C                            in Gev**2. Therefore Q2>0 for spacelike
C                            momentum transfers.
C             iff=1..........Hohler 8.2 parameterization with renormaliztion
C                            according to Friar.
C                =2..........Gari and Krumpelmann
C                =3..........Iachello, Jackson and Lande
C                =4..........dipole parameterization
C                =5..........dipole parameterization, F_1 neutron=0
C                =6..........dipole parameteriztion by Cioffi degli Atti
C                =7..........dipole G_E neutron=0.0
C
C     OUTPUT
C             g(iso,itype)...form factors
C
C             iso=0..........proton
C                 1..........neutron
C             itype=0........electric form factor
C                  =1........magnetic form factor
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION G(0:1,0:1),F(2,2)
C
      IF(IFF.EQ.1) THEN
        CALL HOHLER(Q2,F)
        CALL FTOG(Q2,F,G)
      ELSE IF(IFF.EQ.2) THEN
        CALL GARIK(Q2,F)
        CALL FTOG(Q2,F,G)
      ELSE IF(IFF.EQ.3) THEN
        CALL IJL(Q2,F)
        CALL FTOG(Q2,F,G)
      ELSE
        IFFM3=IFF-3
        CALL DIPOLE(Q2,G,IFFM3)
      ENDIF
C
      RETURN
      END
C
C_______________________________________________________________________
C
C*****dircff
C
      SUBROUTINE DIRCFF(Q2,FP,IFF)
C
C     This subroutine returns Dirac and Pauli form factors for a variety of
C     parameterizations.
C
C
C     INPUT
C             Q2.............negative of square of four-momentum transfer
C                            in Gev**2. Therefore Q2>0 for spacelike
C                            momentum transfers.
C             iff=1..........Hohler 8.2 parameterization with renormaliztion
C                            according to Friar.
C                =2..........Gari and Krumpelmann
C                =3..........Iachello, Jackson and Lande
C                =4..........dipole parameterization
C                =5..........dipole parameterization, F_1 neutron=0
C                =6..........dipole parameteriztion by Cioffi degli Atti
C                =7..........dipole G_E neutron=0.0
C                =0..........F_1p=1 all others =0
C
C     OUTPUT
C             fp(iso,itype)..form factors
C
C             iso=0..........proton
C                 1..........neutron
C             itype=1........Dirac form factor
C                  =2........Pauli form factor
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION G(0:1,0:1),F(2,2),FP(0:1,2)
C
      IF(IFF.EQ.0) THEN
        FP(0,1)=1.
        FP(0,2)=0.
        FP(1,1)=0.
        FP(1,2)=0.
      ELSE IF(IFF.EQ.1) THEN
        CALL HOHLER(Q2,F)
        CALL FTOFP(F,FP)
      ELSE IF(IFF.EQ.2) THEN
        CALL GARIK(Q2,F)
        CALL FTOFP(F,FP)
      ELSE IF(IFF.EQ.3) THEN
        CALL IJL(Q2,F)
        CALL FTOFP(F,FP)
      ELSE
        CALL DIPOLE(Q2,G,IFF-3)
        CALL GTOFP(Q2,G,FP)
      ENDIF
C
      RETURN
      END

C
C_______________________________________________________________________
C
C*****ftog
C
      SUBROUTINE FTOG(Q2,F,G)
C
C     This subroutine translates the representation of f to that of g.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION MN
C
      DIMENSION F(2,2),G(0:1,0:1),FP(0:1,2)
C
      DATA MN/.939/
C
      TAU=Q2/(4.*MN**2)
C
      CALL FTOFP(F,FP)
C
      DO 100 I=0,1
        G(I,0)=FP(I,1)-TAU*FP(I,2)
        G(I,1)=FP(I,1)+FP(I,2)
  100 CONTINUE
C
      RETURN
      END
C
C_______________________________________________________________________
C
C*****ftofp
C
      SUBROUTINE FTOFP(F,FP)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION F(2,2),FP(0:1,2)
C
      DO 100 I=1,2
        FP(0,I)=0.5*(F(2,I)+F(1,I))
        FP(1,I)=0.5*(F(2,I)-F(1,I))
  100 CONTINUE
C
      RETURN
      END
C
C_______________________________________________________________________
C
C*****gtofp
C
      SUBROUTINE GTOFP(Q2,G,FP)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION MN
C
      DIMENSION G(0:1,0:1),FP(0:1,2)
C
      DATA MN/.939/
C
      TAU=Q2/(4.*MN**2)
C
      DO 100 I=0,1
        FP(I,1)=(G(I,0)+TAU*G(I,1))/(1.+TAU)
        FP(I,2)=(G(I,1)-G(I,0))/(1.+TAU)
  100 CONTINUE
C
      RETURN
      END
C
C_______________________________________________________________________
C
C*****hohler
C
      SUBROUTINE HOHLER(Q2,F)
C
C     This subroutine calculates the electromagnetic form factors
C     using the Hohler parameterization 8.2.  Input is in GeV**2
C     and the square of the four-momentum transfer should be
C     positive.
C
C     REFERENCE: G. Hohler, et al., Nucl. Phys. B114, 505 (1976).
C
C     INPUT
C             Q2.............negative of square of four-momentum transfer
C                            in Gev**2. Therefore Q2>0 for spacelike
C                            momentum transfers.
C
C     OUTPUT
C             f(iso,itype)...form factors
C
C             iso=1..........isovector
C                 2..........isoscalar
C             itype=1........Dirac form factor F_1
C                  =2........Pauli form factor F_2
C
C     7-16-87
C            This version has been modified such that it takes into
C            account changes made by J. Friar to this parameterization
C            to give better behavior at the q2=0 point.
C
      IMPLICIT DOUBLE PRECISION(A-H,K-Z)
      DIMENSION SRT(3,2),A(3,2,2),F(2,2),FROV(2)
      DATA((SRT(I,J),I=1,3),J=1,2)/1.21,2.45,2.95,0.783,1.02,1.8/
      DATA(((A(I,J,IJ),I=1,3),J=1,2),IJ=1,2)/0.05,-0.52,0.28,0.71
     #  ,-0.64,-0.13,-1.99,0.1956,0.1856,-0.11,0.13,-0.02/
C
      FROV(1)=(0.955+0.09/(1.+Q2/0.355)**2)/(1.+Q2/0.536)/2.
      FROV(2)=(5.335+0.962/(1.+Q2/0.268))/(1.+Q2/.603)/2.
C
      DO 50 I1=1,2
        DO 50 I2=1,2
          F(I1,I2)=FROV(I2)*(1-(-1)**I1)/2.
          DO 50 I3=1,3
   50 F(I1,I2)=F(I1,I2)+A(I3,I1,I2)/(SRT(I3,I1)**2+Q2)
C
C                  Renormalization according to Friar
C
      F(1,1)=2.*F(1,1)*0.9954211
      F(2,1)=2.*F(2,1)*0.9930881
      F(1,2)=2.*F(1,2)*1.0052166
      F(2,2)=2.*F(2,2)*0.98891518
C
      RETURN
      END
C
C_______________________________________________________________________
C
C*****GariK
C
      SUBROUTINE GARIK(Q2,F)
C
C     This subroutine produces the form factors of Gari and Krumpelmann.
C
C
C     REFERENCE: M. Gari and W. Krumpelmann, Z. Phys. A322 (1985);
C                                            Phys. Lett. 173B, 10 (1986).
C
C     INPUT
C             Q2.............negative of square of four-momentum transfer
C                            in Gev**2. Therefore Q2>0 for spacelike
C                            momentum transfers.
C
C     OUTPUT
C             f(iso,itype)...form factors
C
C             iso=1..........isovector
C                 2..........isoscalar
C             itype=1........Dirac form factor F_1
C                  =2........Pauli form factor F_2
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION MRHO,MOMEGA,KAPV,KAPS,KAPRHO,KAPOMG,LAM1,
     &                 LAM2,LAMQCD
C
      DIMENSION F(2,2)
C
      DATA MRHO/0.7661/,MOMEGA/0.784/,KAPV/3.706/,KAPS/-0.12/,
     &     KAPRHO/6.62/,KAPOMG/0.163/,CRHO/0.377/,COMG/0.411/,
     &     LAM1/0.795/,LAM2/2.27/,LAMQCD/0.29/
C
      Q2TILD=Q2*LOG((LAM2**2+Q2)/LAMQCD**2)/LOG((LAM2/LAMQCD)**2)
C
      D1=1./(1.+Q2TILD/LAM1**2)
      D2=1./(1.+Q2TILD/LAM2**2)
C
      F1=D1*D2
      F2=F1*D2
C
      DRHO=1./(1.+Q2/MRHO**2)
      DOMG=1./(1.+Q2/MOMEGA**2)
C
      F(1,1)=(CRHO*DRHO+1.-CRHO)*F1
      F(2,1)=(COMG*DOMG+1.-COMG)*F1
      F(1,2)=(CRHO*KAPRHO*DRHO+KAPV-CRHO*KAPRHO)*F2
      F(2,2)=(COMG*KAPOMG*DOMG+KAPS-COMG*KAPOMG)*F2
C
      RETURN
      END
C
C_______________________________________________________________________
C
C*****ijl
C
      SUBROUTINE IJL(Q2,F)
C
C     This subroutine produces the form factors of Iachello, Jackson and
C     Lande.
C
C     REFERENCE: F. Iachello, A. D. Jackson and A. Lande, Phys. Lett. 43B,
C                191 (1973).
C
C     INPUT
C             Q2.............negative of square of four-momentum transfer
C                            in Gev**2. Therefore Q2>0 for spacelike
C                            momentum transfers.
C
C     OUTPUT
C             f(iso,itype)...form factors
C
C             iso=1..........isovector
C                 2..........isoscalar
C             itype=1........Dirac form factor F_1
C                  =2........Pauli form factor F_2
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION MPI,MRHO,MOMEGA,MPHI,KAPV,KAPS
C
      DIMENSION F(2,2)
C
      DATA MPI/.13957/,MRHO/0.765/,MOMEGA/0.784/,MPHI/1.019/,
     &     GAMRHO/.112/,GAM/0.25/,BTARHO/0.672/,BTAOMG/1.102/,
     &     BTAPHI/0.112/,ALPPHI/-0.052/,KAPS/-0.120/,KAPV/3.706/,
     &     PI/3.141592741012573/
C
      G=1./(1.+GAM*Q2)**2
C
      IF(Q2.GT.1.E-7) THEN
        ALPHA=2./PI*SQRT((Q2+4.*MPI**2)/Q2)*LOG((SQRT(Q2+4.*MPI**2)+
     &        SQRT(Q2))/(2.*MPI))
      ELSE
        ALPHA=0.
      ENDIF
C
      DRHO=(MRHO**2+8.*GAMRHO*MPI/PI)/(MRHO**2+Q2+(4.*MPI**2+Q2)*
     &GAMRHO*ALPHA/MPI)
      DOMG=1./(1.+Q2/MOMEGA**2)
      DPHI=1./(1.+Q2/MPHI**2)
C
      F(1,1)=G*(1.-BTARHO+BTARHO*DRHO)
      F(2,1)=G*(1-BTAOMG-BTAPHI+BTAOMG*DOMG+BTAPHI*DPHI)
      F(1,2)=G*KAPV*DRHO
      F(2,2)=G*((KAPS-ALPPHI)*DOMG+ALPPHI*DPHI)
C
      RETURN
      END
C
C_______________________________________________________________________
C
C*****dipole
C
      SUBROUTINE DIPOLE(Q2,G,IFF)
C
C     This subroutine produces several form factors which are expressed
C     as dipole fits to the Sachs form factors.
C
C
C     INPUT
C             Q2.............negative of square of four-momentum transfer
C                            in Gev**2. Therefore Q2>0 for spacelike
C                            momentum transfers.
C             iff=1..........Preston and Bhaduri
C                =2..........F_1 neutron=0
C                =3..........Cioffi degli Atti
C                =4..........G_E neutron=0
C
C     OUTPUT
C             g(iso,itype)...form factors
C
C             iso=0..........proton
C                 1..........neutron
C             itype=0........electric form factor
C                  =1........magnetic form factor
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION MN,MUP,MUN
C
      DIMENSION G(0:1,0:1)
C
      DATA MN/0.939/,MUP/2.793/,MUN/-1.913/
C
      TAU=Q2/(4.*MN**2)
C
C                     as in Preston and Bhaduri
C
      IF(IFF.EQ.1) THEN
        G(0,0)=1./(1.+Q2/0.71)**2
        G(1,0)=-MUN*TAU*G(0,0)/(1.+5.6*TAU)
        G(0,1)=MUP*G(0,0)
        G(1,1)=MUN*G(0,0)
C
C           dipole with neutron Dirac form factor set to zero
C
      ELSE IF(IFF.EQ.2) THEN
        G(0,0)=1./(1.+Q2/0.71)**2
        G(1,0)=-MUN*TAU*G(0,0)
        G(0,1)=MUP*G(0,0)
        G(1,1)=MUN*G(0,0)
C
C           dipole fit according to Cioffi degli Atti
C
      ELSE IF(IFF.EQ.3) THEN
        G(0,0)=1./(1.+Q2/0.695)**2
        G(1,0)=-MUN*TAU*G(0,0)/(1.+4.*TAU)
        G(0,1)=MUP*G(0,0)
        G(1,1)=MUN*G(0,0)
C
C                     G_E neutron = 0.
C
      ELSE IF(IFF.EQ.4) THEN
        G(0,0)=1./(1.+Q2/0.71)**2
        G(1,0)=0.0
        G(0,1)=MUP*G(0,0)
        G(1,1)=MUN*G(0,0)
      ENDIF
      RETURN
      END
