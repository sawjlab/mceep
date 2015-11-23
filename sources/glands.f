*
* $Id: glands.F,v 1.1.1.1 1995/10/24 10:21:25 cernlib Exp $
*
* $Log: glands.F,v $
* Revision 1.1.1.1  1995/10/24 10:21:25  cernlib
* Geant
*
*
*#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      FUNCTION GLANDS(X)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Copy of the CERN library routine DSTLAN (G110)                *
C.    *                                                                *
C.    *    ==>Called by : GVACOE                                       *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION P1(0:4),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3)
      DIMENSION Q1(0:4),Q2(0:3),Q3(0:3),Q4(0:3),Q5(0:3),Q6(0:3)
      DIMENSION A1(1:3),A2(1:3)
C
      DATA (P1(I),I=0,4),(Q1(J),J=0,4)
     1/ 0.25140 91491E+0,-0.62505 80444E-1, 0.14583 81230E-1,
     2 -0.21088 17737E-2, 0.74112 47290E-3,
     3  1.0             ,-0.55711 75625E-2, 0.62253 10236E-1,
     4 -0.31373 78427E-2, 0.19314 96439E-2/
C
      DATA (P2(I),I=0,3),(Q2(J),J=0,3)
     1/ 0.28683 28584E+0, 0.35643 63231E+0, 0.15235 18695E+0,
     2  0.22513 04883E-1,
     3  1.0             , 0.61911 36137E+0, 0.17207 21448E+0,
     4  0.22785 94771E-1/
C
      DATA (P3(I),I=0,3),(Q3(J),J=0,3)
     1/ 0.28683 29066E+0, 0.30038 28436E+0, 0.99509 51941E-1,
     2  0.87338 27185E-2,
     3  1.0             , 0.42371 90502E+0, 0.10956 31512E+0,
     4  0.86938 51567E-2/
C
      DATA (P4(I),I=0,3),(Q4(J),J=0,3)
     1/ 0.10003 51630E+1, 0.45035 92498E+1, 0.10858 83880E+2,
     2  0.75360 52269E+1,
     3  1.0             , 0.55399 69678E+1, 0.19335 81111E+2,
     4  0.27213 21508E+2/
C
      DATA (P5(I),I=0,3),(Q5(J),J=0,3)
     1/ 0.10000 06517E+1, 0.49094 14111E+2, 0.85055 44753E+2,
     2  0.15321 53455E+3,
     3  1.0             , 0.50099 28881E+2, 0.13998 19104E+3,
     4  0.42000 02909E+3/
C
      DATA (P6(I),I=0,3),(Q6(J),J=0,3)
     1/ 0.10000 00983E+1, 0.13298 68456E+3, 0.91621 49244E+3,
     2 -0.96050 54274E+3,
     3  1.0             , 0.13398 87843E+3, 0.10559 90413E+4,
     4  0.55322 24619E+3/
C
      DATA (A1(I),I=1,3)
     1/-0.45833 33333E+0, 0.66753 47222E+0,-0.16417 41416E+1/
C
      DATA (A2(I),I=1,3)
     1/ 1.0             ,-0.42278 43351E+0,-0.20434 03138E+1/
C
      V=X
      IF(V .LT. -5.5) THEN
       U=EXP(V+1.0)
       GLANDS=0.3989422803*EXP(-1.0/U)*SQRT(U)*
     1        (1.0+(A1(1)+(A1(2)+A1(3)*U)*U)*U)
      ELSE IF(V .LT. -1.0) THEN
       U=EXP(-V-1.0)
       GLANDS=(EXP(-U)/SQRT(U))*
     1        (P1(0)+(P1(1)+(P1(2)+(P1(3)+P1(4)*V)*V)*V)*V)/
     2        (Q1(0)+(Q1(1)+(Q1(2)+(Q1(3)+Q1(4)*V)*V)*V)*V)
      ELSE IF(V .LT. 1.0) THEN
       GLANDS=(P2(0)+(P2(1)+(P2(2)+P2(3)*V)*V)*V)/
     1        (Q2(0)+(Q2(1)+(Q2(2)+Q2(3)*V)*V)*V)
      ELSE IF(V .LT. 4.0) THEN
       GLANDS=(P3(0)+(P3(1)+(P3(2)+P3(3)*V)*V)*V)/
     1        (Q3(0)+(Q3(1)+(Q3(2)+Q3(3)*V)*V)*V)
      ELSE IF(V .LT. 12.0) THEN
       U=1.0/V
       GLANDS=(P4(0)+(P4(1)+(P4(2)+P4(3)*U)*U)*U)/
     1        (Q4(0)+(Q4(1)+(Q4(2)+Q4(3)*U)*U)*U)
      ELSE IF(V .LT. 50.0) THEN
       U=1.0/V
       GLANDS=(P5(0)+(P5(1)+(P5(2)+P5(3)*U)*U)*U)/
     1        (Q5(0)+(Q5(1)+(Q5(2)+Q5(3)*U)*U)*U)
      ELSE IF(V .LT. 300.0) THEN
       U=1.0/V
       GLANDS=(P6(0)+(P6(1)+(P6(2)+P6(3)*U)*U)*U)/
     1        (Q6(0)+(Q6(1)+(Q6(2)+Q6(3)*U)*U)*U)
      ELSE
       U=1.0/(V-V*LOG(V)/(V+1.0))
       GLANDS=1.0-(A2(1)+(A2(2)+A2(3)*U)*U)*U
      END IF
c      write(6,*) 'glandsin',GLANDS
      RETURN
      END
 
