C
C ----------------------------------------------------------------
C     Function RAN2_MCEEP
C
C       (Renamed from RAN2 as of Version 3.9, to avoid conflict
C        with RAN2 from bremsgen.f)
C
C       random number routines from Numerical Recipes book
C
C       I) uniform deviate on [0,1) with no sequential correlations
C ----------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RAN2_MCEEP(IDUM)
      IMPLICIT NONE
      DOUBLE PRECISION RM
      INTEGER IDUM,M,IA,IC,IR,IFF,J,IY
C
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.D0/M)
      DIMENSION IR(97)
      DATA IFF/0/
      IF(IDUM .LT. 0 .OR. IFF .EQ. 0)THEN
        IFF = 1
        IDUM = MOD(IC-IDUM,M)
        DO J=1,97
          IDUM = MOD(IA*IDUM+IC,M)
          IR(J) = IDUM
        ENDDO
        IDUM = MOD(IA*IDUM+IC,M)
        IY = IDUM
      ENDIF
      J = 1 + (97*IY)/M
      IF(J .GT. 97 .OR. J .LT. 1)PAUSE
      IY = IR(J)
      RAN2_MCEEP = IY*RM
      IDUM = MOD(IA*IDUM+IC,M)
      IR(J) = IDUM
      RETURN
      END
C
C ---------------------------------------------------------------------
C     Function GASDEV
C
C       II) normal deviate (zero mean and specified sigma), uses RANECU
C           normalized to unity on (-infty,infty)
C ---------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION GASDEV(SIGMA)
      IMPLICIT NONE
      DOUBLE PRECISION SIGMA,V1,V2,R,FAC,GSET
      REAL RAN_VEC(2)
      INTEGER ISET
C
      DATA ISET/0/
      IF(ISET .EQ. 0)THEN
    1   CALL RANECU(RAN_VEC,2)
        V1 = 2.D0*RAN_VEC(1) - 1.D0
        V2 = 2.D0*RAN_VEC(2) - 1.D0
        R = V1*V1 + V2*V2
        IF(R .GE. 1.D0)GOTO 1
        FAC = SIGMA*SQRT(-2.D0*LOG(R)/R)
        GSET = V1*FAC
        GASDEV = V2*FAC
        ISET = 1
      ELSE
        GASDEV = GSET
        ISET = 0
      ENDIF
      RETURN
      END
C
C ----------------------------------------------------------------
C     Function GASDEVV
C ----------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION GASDEVV(SIGMA)
      IMPLICIT NONE
      DOUBLE PRECISION SIGMA,V1,V2,R,FAC
      REAL RAN_VEC(2)
C
    1 CALL RANECU(RAN_VEC,2)
      V1 = 2.D0*RAN_VEC(1) - 1.D0
      V2 = 2.D0*RAN_VEC(2) - 1.D0
      R = V1*V1 + V2*V2
      IF(R .GE. 1.D0)GOTO 1
      FAC = SIGMA*SQRT(-2.D0*LOG(R)/R)
      GASDEVV = V2*FAC
      RETURN
      END
