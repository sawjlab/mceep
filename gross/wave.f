C________________________________________________________________________
C
C*****wavein
C
C          This subroutine reads the tabulated deuteron wave function
C       components from nunit. The wave function is normalized such that
C       it is in units of (MeV)^{-3/2}. The wave function table is then
C       reorganized in a form more convenient for the expression of the
C       Gross wave functions in the helicity basis.
C
      SUBROUTINE WAVEIN(NUNIT)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER(NWAVE=4,MNDATA=32)
      PARAMETER(PI=3.141592653589793,TWOPI3=(2.*PI)**3)
C
      COMMON/WAVPAR/Q(MNDATA),Z(MNDATA,NWAVE),Z2(MNDATA,NWAVE),
     &              ZDATA(MNDATA,NWAVE),QLN(MNDATA),WT(MNDATA),
     &              ANRMNR,NDATA
C
      DIMENSION PROB(NWAVE)
C
      READ(NUNIT,*) NDATA
      DO I=1,NDATA
        READ(NUNIT,*) Q(I),(ZDATA(I,J),J=1,NWAVE),WT(I)
      END DO
C
C                        normalize wave function
C
      DO IWAVE=1,NWAVE
        PROB(IWAVE)=0.0
      END DO
      OFF=0.0
C
      DO IDATA=1,NDATA
        DO IWAVE=1,NWAVE
          PROB(IWAVE)=PROB(IWAVE)+WT(IDATA)*Q(IDATA)**2*
     &                ZDATA(IDATA,IWAVE)**2
        END DO
        OFF=OFF+WT(IDATA)*Q(IDATA)**2*ZDATA(IDATA,1)*ZDATA(IDATA,2)
      END DO
C
      ANORM=0.0
      DO IWAVE=1,NWAVE
        ANORM=ANORM+PROB(IWAVE)
      END DO
C
      OFF=OFF/(PROB(1)+PROB(2))
      DO IWAVE=1,NWAVE
        PROB(IWAVE)=PROB(IWAVE)/ANORM
      END DO
C
      ANRMNR=1./SQRT(PROB(1)+PROB(2))

C
C
C                Note: The wave function must be normalized
C                      with the factor of (2.*pi)**3 to be
C                      consistent with the fourier transform
C                      conventions used in constructing the
C                      cross section expressions.
C
      ANORM=SQRT(ANORM/TWOPI3)
C
C      write(6,*) anorm
C      write(6,*) prob
C      write(6,*) anrmnr,off
C
      DO I=1,NDATA
        DO J=1,NWAVE
          ZDATA(I,J)=ZDATA(I,J)/ANORM
        END DO
      END DO
C
      RETURN
      END
C
C_________________________________________________________________________
C
C*****wavset
C
C           This subroutine sets up a cubic spline fit to the wave
C        function components as a function of the natural log of the
C        magnitude of the momentum. The routine GROSS_SPLINE is called to
C        initialize the spline parameters for this fit.
C
      SUBROUTINE WAVSET
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER(NWAVE=4,MNDATA=32)
C
      COMMON/WAVPAR/Q(MNDATA),Z(MNDATA,NWAVE),Z2(MNDATA,NWAVE),
     &              ZDATA(MNDATA,NWAVE),QLN(MNDATA),WT(MNDATA),
     &              ANRMNR,NDATA
C
      COMMON/XTRAP/XPAR(NWAVE,2,2),IPUP,LWAVE(NWAVE)
C
      LWAVE(1)=0
      LWAVE(2)=2
      LWAVE(3)=1
      LWAVE(4)=1
C
      YP=1.E30
      IPUP=NDATA-3
C
      DO 200 I=1,NDATA
        QLN(I)=LOG(Q(I))
C
        DO 100 J=1,NWAVE
          Z(I,J)=ZDATA(I,J)
  100   CONTINUE
C
  200 CONTINUE
C
      CALL GROSS_SPLINE(QLN,Z,NDATA,NWAVE,YP,YP,Z2)
C
C           find parameters for fit to low p such that each
C           wave function goes as p**l
C
      ILOW=2
      IUP=3
      DO J=1,NWAVE
        XPAR(J,1,1)=(ZDATA(IUP,J)/Q(IUP)**LWAVE(J)
     &              -ZDATA(ILOW,J)/Q(ILOW)**LWAVE(J))
     &              /(Q(IUP)-Q(ILOW))
        XPAR(J,2,1)=(Q(IUP)*ZDATA(ILOW,J)/Q(ILOW)**LWAVE(J)
     &              -Q(ILOW)*ZDATA(IUP,J)/Q(IUP)**LWAVE(J))
     &              /(Q(IUP)-Q(ILOW))
      END DO
C
C           find parameters for exponential fit to tail
C
      ILOW=IPUP
      IUP=ILOW+1
      DO J=1,NWAVE
        XPAR(J,1,2)=(LOG(ABS(ZDATA(IUP,J)))-LOG(ABS(
     &                  ZDATA(ILOW,J))))/(QLN(IUP)-QLN(ILOW))
        XPAR(J,2,2)=(QLN(IUP)*LOG(ABS(ZDATA(ILOW,J)))-QLN(ILOW)*
     &                  LOG(ABS(ZDATA(IUP,J))))/(QLN(IUP)-QLN(ILOW))
      END DO
C
C      write(6,1) ((xpar(j,1,i),i=1,2),j=1,nwave)
C 1    format(2(2x,1pe12.5))
C
      RETURN
      END
C
C__________________________________________________________________________
C
C*****wav
C
C          This subroutine returns the for wave function components in
C       the array zp as interpolated for a given in input momentum p.
C       p is assumed to be in MeV. The cubic spline routine GROSS_SPLINT
C       performs the interpolation.
C
      SUBROUTINE WAV(P,ZP)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER(NWAVE=4,MNDATA=32)
C
      COMMON/WAVPAR/Q(MNDATA),Z(MNDATA,NWAVE),Z2(MNDATA,NWAVE),
     &              ZDATA(MNDATA,NWAVE),QLN(MNDATA),WT(MNDATA),
     &              ANRMNR,NDATA
C
      COMMON/XTRAP/XPAR(NWAVE,2,2),IPUP,LWAVE(NWAVE)
C
      DIMENSION ZP(NWAVE)
C

      IF(P.NE.0.0) PLN=LOG(P)
C
      IF(P.LT.Q(2)) THEN
C
        DO I=1,NWAVE
          ZP(I)=P**LWAVE(I)*(XPAR(I,1,1)*P+XPAR(I,2,1))
        END DO
C
      ELSE IF(PLN.GT.QLN(IPUP)) THEN
C
        DO I=1,NWAVE
          ZP(I)=SIGN(EXP(XPAR(I,1,2)*PLN+XPAR(I,2,2)),
     &               ZDATA(IPUP,I))
        END DO
C
      ELSE
        CALL GROSS_SPLINT(QLN,Z,Z2,NDATA,NWAVE,PLN,ZP)
      ENDIF
C
C
      RETURN
      END
C
C______________________________________________________________________
C
C*****deuwav
C
C         This subroutine constructs the Gross deuteron wave functions.
C
C       INPUT
C           p.............magnitude of momentum in MeV
C           irel=1........relativistic wave functions
C                0........nonrelativistic wave functions
C
C       OUTPUT
C           a(lnuc,lnucp,irho).........array containing deuteron wave
C                                      functions
C
C             lnuc=0.........positive nucleon helicity
C                 =1.........negative nucleon helicity
C             lnucp=0........positive nucleon helicity
C                  =1........negative nucleon helicity
C             ienrgy=0.......positive energy projections
C                    1.......negative energy projections
C
C
      SUBROUTINE DEUWAV(P,IREL,A)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER(NWAVE=4,MNDATA=32)
C
      COMMON/WAVPAR/Q(MNDATA),Z(MNDATA,NWAVE),Z2(MNDATA,NWAVE),
     &              ZDATA(MNDATA,NWAVE),QLN(MNDATA),WT(MNDATA),
     &              ANRMNR,NDATA
C
      DIMENSION ZP(NWAVE),A(0:1,0:1,0:1)
C
      DOUBLE PRECISION MN,MT
      COMMON/DPARA/MN,MT
C
      PARAMETER(ROOT2=1.414213562373095,ROOT3=1.732050807568877)
      PARAMETER(PI=3.141592653589793,FAC=0.199471140)
C
C           obtain the wave function components and the Dirac spinors
C
      CALL WAV(P,ZP)
C
C              construct new combinations of wave matrix
C
      IF(IREL.EQ.1) THEN               ! relativistic wave functions
        POMN=P/MN
        EPOMN=SQRT(POMN**2+1.)
        A(0,0,0)=FAC*(ZP(1)+ROOT2*ZP(2)+ROOT3*POMN*ZP(4))
        A(1,1,0)=A(0,0,0)
        A(1,0,0)=FAC*(ROOT2*ZP(1)-ZP(2)-ROOT3*POMN*ZP(3))
        A(0,1,0)=A(1,0,0)
        A(0,0,1)=-FAC*EPOMN*ROOT3*ZP(4)
        A(1,1,1)=-A(0,0,1)
        A(1,0,1)=FAC*EPOMN*ROOT3*ZP(3)
        A(0,1,1)=-A(1,0,1)
      ELSE IF(IREL.EQ.0) THEN          ! nonrelativistic wave functions
        A(0,0,0)=FAC*ANRMNR*(ZP(1)+ROOT2*ZP(2))
        A(1,1,0)=A(0,0,0)
        A(1,0,0)=FAC*ANRMNR*(ROOT2*ZP(1)-ZP(2))
        A(0,1,0)=A(1,0,0)
        A(0,0,1)=0.0
        A(1,1,1)=0.0
        A(1,0,1)=0.0
        A(0,1,1)=0.0
      END IF
C
      RETURN
      END
C
C______________________________________________________________________
C
      SUBROUTINE GROSS_SPLINE(X,Y,N,M,YP1,YPN,Y2)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N,M),Y2(N,M),U(NMAX)
C
      DO 100 J=1,M
C
        IF (YP1.GT..99E30) THEN
          Y2(1,J)=0.
          U(1)=0.
        ELSE
          Y2(1,J)=-0.5
          U(1)=(3./(X(2)-X(1)))*((Y(2,J)-Y(1,J))/(X(2)-X(1))-YP1)
        ENDIF
        DO 11 I=2,N-1
          SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
          P=SIG*Y2(I-1,J)+2.
          Y2(I,J)=(SIG-1.)/P
          U(I)=(6.*((Y(I+1,J)-Y(I,J))/(X(I+1)-X(I))-(Y(I,J)-Y(I-1,J))
     *        /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
   11   CONTINUE
        IF (YPN.GT..99E30) THEN
          QN=0.
          UN=0.
        ELSE
          QN=0.5
          UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N,J)-Y(N-1,J))/(X(N)-X(N-1)))
        ENDIF
        Y2(N,J)=(UN-QN*U(N-1))/(QN*Y2(N-1,J)+1.)
        DO 12 K=N-1,1,-1
          Y2(K,J)=Y2(K,J)*Y2(K+1,J)+U(K)
   12   CONTINUE
C
  100 CONTINUE
      RETURN
      END
C
C__________________________________________________________________________
C
      SUBROUTINE GROSS_SPLINT(XA,YA,Y2A,N,M,X,Y)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION XA(N),YA(N,M),Y2A(N,M),Y(M)
      KLO=1
      KHI=N
    1 IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
C
      DO 100 J=1,M
        Y(J)=A*YA(KLO,J)+B*YA(KHI,J)+
     *      ((A**3-A)*Y2A(KLO,J)+(B**3-B)*Y2A(KHI,J))*(H**2)/6.
  100 CONTINUE
C
      RETURN
      END
