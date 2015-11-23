      SUBROUTINE GAULEG(X1,X2,X,W,N)
C
C     This subroutine generates gauss legendre points and weights on the
C     interval x1 to x2.  The routine is "Numerical Recipes" by Press,
C     Flannery, Teukolsky, and Vettering.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      REAL*4 X1,X2,X(N),W(N)
      DIMENSION X(N),W(N)
C
      PARAMETER(PI=3.141592653589793)
C
      PARAMETER (EPS=3.D-16)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=COS(PI*(I-.25D0)/(N+.5D0))
    1   CONTINUE
        P1=1.D0
        P2=0.D0
        DO 11 J=1,N
          P3=P2
          P2=P1
          P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
   11   CONTINUE
        PP=N*(Z*P1-P2)/(Z*Z-1.D0)
        Z1=Z
        Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
   12 CONTINUE
      RETURN
      END
