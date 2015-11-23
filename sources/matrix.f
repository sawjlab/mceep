C
C --------------------------------------------------------------------
C --------------------------------------------------------------------
C       SUBROUTINE MATINV
C
C       Purpose:
C               Invert a symmetric matrix
C
C       Usage:
C               CALL MATINV( ARRAY, ARRAY_INV, NORDER )
C
C       Description of parameters:
C               ARRAY     - input matrix
C               ARRAY_INV - inverse matrix
C               NORDER    - degree of matrix
C
C       Subroutines and function subprograms required:
C               none
C
C       Based on routine of Bevington.
C
C --------------------------------------------------------------------
C
      SUBROUTINE MATINV (ARRAY, ARRAY_INV, NORDER  )
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ARRAY(NORDER,NORDER),ARRAY_INV(NORDER,NORDER)
      INTEGER IK(20), JK(20)
C
      DO I=1,NORDER
        DO J=1,NORDER
          ARRAY_INV(I,J) = ARRAY(I,J)
        ENDDO
      ENDDO
C
      DO K=1,NORDER
C
C       FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX
C
        AMAX = 0.D0
   21   DO I=K,NORDER
          DO J=K, NORDER
            IF(  ABS(AMAX) -  ABS(ARRAY_INV(I,J) ) ) 24,24,30
   24       AMAX = ARRAY_INV( I,J )
            IK(K) = I
            JK(K) = J
   30     ENDDO
        ENDDO
C
C          INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAY_INV(K,K)
C
        IF( AMAX ) 41,140, 41
   41   I=IK(K)
        IF( I-K) 21,51,43
C
   43   DO J=1,NORDER
          SAVE = ARRAY_INV( K,J)
          ARRAY_INV(K,J) = ARRAY_INV(I,J)
          ARRAY_INV(I,J) = -SAVE
        ENDDO
C
   51   J = JK(K)
        IF(J-K) 21,61,53
C
   53   DO I=1,NORDER
          SAVE = ARRAY_INV(I,K)
          ARRAY_INV(I,K) = ARRAY_INV(I,J)
          ARRAY_INV(I,J) = -SAVE
        ENDDO
C
C          ACCUMULATE ELEMENTS OF INVERSE MATRIX
C
   61   DO I=1,NORDER
          IF(I-K) 63,70,63
   63     ARRAY_INV(I,K) = -ARRAY_INV(I,K)/AMAX
   70   ENDDO
C
   71   DO I=1,NORDER
          DO J=1,NORDER
            IF( I-K) 74,80,74
   74       IF( J-K) 75,80,75
   75       ARRAY_INV(I,J) = ARRAY_INV(I,J)
     #                  + ARRAY_INV(I,K)*ARRAY_INV(K,J)
   80     ENDDO
        ENDDO
C
   81   DO J=1,NORDER
          IF(J-K) 83,90,83
   83     ARRAY_INV(K,J) = ARRAY_INV(K,J)/AMAX
   90   ENDDO
C
        ARRAY_INV(K,K) = 1./AMAX
      ENDDO
C
C       RESTORE ORDERING OF MATRIX
C
  101 DO L=1, NORDER
        K = NORDER-L+1
        J = IK(K)
        IF(J-K) 111,111,105
C
  105   DO I=1,NORDER
          SAVE = ARRAY_INV(I,K)
          ARRAY_INV(I,K) = -ARRAY_INV(I,J)
          ARRAY_INV(I,J) = SAVE
        ENDDO
C
  111   I = JK(K)
        IF(I-K) 130,130,113
C
  113   DO J=1,NORDER
          SAVE = ARRAY_INV(K,J)
          ARRAY_INV(K,J) = -ARRAY_INV(I,J)
          ARRAY_INV(I,J) = SAVE
        ENDDO
  130 ENDDO
  140 RETURN
      END
C
C ----------------------------------------------------------------------
C     
C     SUBROUTINE MAT_MULT
C     AUTHOR:   M. Nozar
C     DATE:     19-JUL-1991
C     PURPOSE:
C               Computes product of any two matrices:
C               matprod(mxp) = matrix1(mxn) * matrix2(nxp)         
C -----------------------------------------------------------------------
C
      SUBROUTINE MAT_MULT(m,n,p,MATRIX1,MATRIX2,MATPROD)
      IMPLICIT NONE
C      
      INTEGER            m,n,p,I,J,K  
      DIMENSION          MATRIX1(m,n),MATRIX2(n,p),MATPROD(m,p)
      DOUBLE PRECISION   MATRIX1, MATRIX2, MATPROD
C
C ----------------------------------------------------------------------
C     Multiply two matrices
C ----------------------------------------------------------------------
C
      DO I = 1,m
        DO J = 1,p
          MATPROD(I,J) = 0.D0
          DO K = 1,n
            MATPROD(I,J) = MATPROD(I,J)+MATRIX1(I,K)*MATRIX2(K,J) 
          ENDDO         
        ENDDO
      ENDDO
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C      SUBROUTINE  ROTATE_X
C      AUTHOR:     M. Nozar and P.E. Ulmer  
C      DATE:       18-JUL-91
C      PURPOSE:    
C                  Perform PASSIVE (i.e. rotation of the coordinate 
C                  system) rotation by angle theta about the X axis
C                  where C_TH and S_TH are the cosine and sine of theta.
C
C                  NOTE: The following convention is used.
C                        Positive angles imply rotation of the Y axis 
C                        into the direction of the original Z axis
C                        (i.e. right hand rule: positive angles imply
C                        positive sense of rotation with respect to X).
C   
C ----------------------------------------------------------------------
C     
      SUBROUTINE   ROTATE_X(C_TH,S_TH,R_X)
      IMPLICIT          NONE
C         
      DOUBLE PRECISION  C_TH,S_TH,R_X(3,3)
      INTEGER           I,J
C
C ----------------------------------------------------------------------
C     Initialize to the unit matrix.
C ----------------------------------------------------------------------
C
      DO I = 1,3
        DO J = 1,3
           IF (I.EQ.J) THEN
              R_X(I,J) = 1.D0
           ELSE
              R_X(I,J) = 0.D0
           ENDIF
        ENDDO
      ENDDO
C 
      R_X(2,2) =  C_TH
      R_X(2,3) =  S_TH
      R_X(3,2) = -S_TH
      R_X(3,3) =  C_TH
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C     SUBROUTINE  ROTATE_Y
C     AUTHOR:     M. Nozar and P.E. Ulmer  
C     DATE:       18-JUL-91
C     PURPOSE:    
C                 Perform PASSIVE (i.e. rotation of the coordinate 
C                 system) rotation by angle theta about the Y axis
C                 where C_TH and S_TH are the cosine and sine of theta.
C
C                 NOTE: The following convention is used.
C                       Positive angles imply rotation of the Z axis 
C                       into the direction of the original X axis
C                       (i.e. right hand rule: positive angles imply
C                       positive sense of rotation with respect to Y).
C   
C ----------------------------------------------------------------------
C
      SUBROUTINE   ROTATE_Y(C_TH,S_TH,R_Y)
      IMPLICIT          NONE
C        
      DOUBLE PRECISION  C_TH,S_TH,R_Y(3,3)
      INTEGER           I,J
C
C ----------------------------------------------------------------------
C     Initialize to the unit matrix.
C ----------------------------------------------------------------------
C
      DO I = 1,3
        DO J = 1,3
          IF (I.EQ.J) THEN
             R_Y(I,J) = 1.D0
          ELSE
             R_Y(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C      
      R_Y(1,1) =  C_TH
      R_Y(1,3) = -S_TH
      R_Y(3,1) =  S_TH
      R_Y(3,3) =  C_TH
C 
      RETURN
      END
C
C ----------------------------------------------------------------------
C     SUBROUTINE  ROTATE_Z
C     AUTHOR:     M. Nozar and P.E. Ulmer  
C     DATE:       18-JUL-91
C     PURPOSE:    
C                 Perform PASSIVE (i.e. rotation of the coordinate 
C                 system) rotation by angle theta about the Z axis
C                 where C_TH and S_TH are the cosine and sine of theta.
C
C                 NOTE: The following convention is used.
C                       Positive angles imply rotation of the X axis 
C                       into the direction of the original Y axis
C                       (i.e. right hand rule: positive angles imply
C                       positive sense of rotation with respect to Z).
C   
C ----------------------------------------------------------------------
C
      SUBROUTINE   ROTATE_Z(C_TH,S_TH,R_Z)
      IMPLICIT          NONE
C         
      DOUBLE PRECISION  C_TH,S_TH,R_Z(3,3)
      INTEGER           I,J
C
C ----------------------------------------------------------------------
C     Initialize to the unit matrix.
C ----------------------------------------------------------------------
C
      DO I = 1,3
        DO J = 1,3
          IF (I.EQ.J) THEN
             R_Z(I,J) = 1.D0
          ELSE
             R_Z(I,J) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
      R_Z(1,1) =  C_TH
      R_Z(1,2) =  S_TH
      R_Z(2,1) = -S_TH
      R_Z(2,2) =  C_TH
C 
      RETURN
      END

C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     SUBROUTINE V3DIF
C     PURPOSE:
C                Takes difference of two 3-vectors.
C-----------------------------------------------------------------------
C
      SUBROUTINE V3DIF(V1,V2,VD)
      DOUBLE PRECISION V1(3),V2(3),VD(3)
      DO I=1,3
        VD(I) = V1(I) - V2(I)
      ENDDO
      RETURN
      END
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     SUBROUTINE V3SUM
C     PURPOSE:
C                Takes sum of two 3-vectors.
C-----------------------------------------------------------------------
C
      SUBROUTINE V3SUM(V1,V2,VS)
      DOUBLE PRECISION V1(3),V2(3),VS(3)
      DO I=1,3
        VS(I) = V1(I) + V2(I)
      ENDDO
      RETURN
      END
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     SUBROUTINE V3MAKE
C     PURPOSE:     
C                Makes a 3-vector from a magnitude and two angles 
C                defined relative to the "standard" coordinate system:
C
C       All vector quantities are expressed in the "initial electron
C       TRANSPORT coordinate system".   This system is the nominal
C       beam system (i.e. for beam bend angle = 0 ; this corresponds
C       to the Hall floor).  Axes are:
C
C               X - given by YxZ
C               Y - perpendicular to floor (upwards = +)
C               Z - nominal beam direction.
C
C       Positive beam bend angles imply an upward direction for the beam.
C       All angles must be in radians.
C-----------------------------------------------------------------------
C
      SUBROUTINE V3MAKE(VMAG,VPHI,VTHE,VEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION VEC(3)
C
      CTHE = COS(VTHE)
      STHE = SIN(VTHE)
      CPHI = COS(VPHI)
      SPHI = SIN(VPHI)
C
      VEC(1) = VMAG*CTHE*SPHI
      VEC(2) = VMAG*STHE
      VEC(3) = VMAG*CTHE*CPHI
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C    Function DOT(A,B)
C
C    Purpose:
C
C             Takes dot product of two 3-vectors.
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DOT(A,B)
      IMPLICIT DOUBLE PRECISION (A-H,K,O-Z)
      DOUBLE PRECISION A(3),B(3)
      DOT = 0.
      DO I = 1,3
        DOT = DOT + A(I)*B(I)
      ENDDO
      RETURN
      END
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C    SUBROUTINE VEC_CROSS(A,B,C)
C
C    Purpose:
C
C         Takes vector cross product of two 3-vectors:
C              C = A X B
C-----------------------------------------------------------------------
C
      SUBROUTINE VEC_CROSS(A,B,C)
      IMPLICIT NONE
      DOUBLE PRECISION A(3),B(3),C(3)
C
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
C
      RETURN
      END
