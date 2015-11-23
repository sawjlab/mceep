C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C     SUBROUTINE     ROT_TRPT
C     AUTHOR:        M.D. Nozar and P.E. Ulmer
C     DATE:          6-AUG-1991
C     PURPOSE:       
C                    This subroutine will perform a PASSIVE rotation 
C                    about a given axis by R_ANGLE where 
C                    C_ANGLE and S_ANGLE are the cosine and
C                    sine of R_ANGLE.
C
C                    The Transport coordinates are expressed relative to
C                    the rotated coordinate system.  Note that the
C                    positional coordinates (X and Y) reflect the
C                    intersection of the particle ray with the XY plane
C                    of the rotated system (i.e. Z=0 plane).  For this
C                    purpose a traceback is performed from the
C                    original position to the rotated XY coordinate
C                    plane.  The path length, L, also reflects this
C                    traceback.  Thus, the path length represents the
C                    distance to the rotated system for the particle
C                    relative to that for the central ray.
C
C                    NOTE: Conventions used are specified in subroutines
C                          ROTATE_X, ROTATE_Y, ROTATE_Z.
C
C ----------------------------------------------------------------------
C
       SUBROUTINE ROT_TRPT(VECT,VECT_ROT,R_AXIS,C_ANGLE,S_ANGLE,
     #              C_RAY,C_RAY_ROT,COMPUTE_PATH)
C
       IMPLICIT          NONE
C
       CHARACTER*1       R_AXIS
       DOUBLE PRECISION  R_MAT(3,3),POS_TMP(3),VEL_TMP(3),POS(3),VEL(3)
       DOUBLE PRECISION  VECT(6),VECT_ROT(6),C_ANGLE,S_ANGLE
       DOUBLE PRECISION  TRACE(3),TRACE_LEN,DOT
       DOUBLE PRECISION  TTHET,TPHIT,C_RAY(3),C_RAY_ROT(3)
       LOGICAL           COMPUTE_PATH
C
C ----------------------------------------------------------------------
C      Define the position vector.  Note that, by definition, the
C      Z coordinate is zero.  This is because the Transport coordinates
C      are defined in the XY plane (i.e. Z=0 plane).
C ----------------------------------------------------------------------
C
       POS_TMP(1) = VECT(1)
       POS_TMP(2) = VECT(3)
       POS_TMP(3) = 0.D0
C
C ----------------------------------------------------------------------
C      Define a vector in the direction of the particle's velocity.
C ----------------------------------------------------------------------
C
       VEL_TMP(1) = TAN(1.D-3*VECT(2))  !Transport angles are in mr
       VEL_TMP(2) = TAN(1.D-3*VECT(4))  !Transport angles are in mr
       VEL_TMP(3) = 1.D0
C
C ----------------------------------------------------------------------
C      Obtain the rotation matrix, R_MAT.
C ----------------------------------------------------------------------
C
       IF (R_AXIS .EQ. 'X') THEN
        CALL ROTATE_X(C_ANGLE,S_ANGLE,R_MAT)
       ELSEIF (R_AXIS .EQ. 'Y') THEN
        CALL ROTATE_Y(C_ANGLE,S_ANGLE,R_MAT)
       ELSE
        CALL ROTATE_Z(C_ANGLE,S_ANGLE,R_MAT)
       ENDIF   
C
C ----------------------------------------------------------------------
C      Calculate components of various vectors after rotation.
C
C            POS:    the position  vector of the particle
C            VEL:    the direction vector of the particle
C            C_RAY:  the direction vector of the central ray
C ----------------------------------------------------------------------
C
       CALL MAT_MULT(3,3,1,R_MAT,POS_TMP,POS)
       CALL MAT_MULT(3,3,1,R_MAT,VEL_TMP,VEL)
       IF(COMPUTE_PATH) CALL MAT_MULT(3,3,1,R_MAT,C_RAY,C_RAY_ROT)
C
C ----------------------------------------------------------------------
C      Determine the Transport angles in the rotated system
C      in terms of the rotated directional vector.
C ----------------------------------------------------------------------
C
       TTHET = VEL(1)/VEL(3)                 !tan(theta_T)
       TPHIT = VEL(2)/VEL(3)                 !tan(phi_T)
       VECT_ROT(2) = 1.D3 * ATAN(TTHET)      !in mr
       VECT_ROT(4) = 1.D3 * ATAN(TPHIT)      !in mr
C
C ----------------------------------------------------------------------
C      Obtain positional (X,Y) components of Transport vector after
C      traceback to the XY plane of the rotated system.
C ----------------------------------------------------------------------
C
       VECT_ROT(1) = POS(1) - POS(3) * TTHET
       VECT_ROT(3) = POS(2) - POS(3) * TPHIT
C
C ----------------------------------------------------------------------
C      Define a vector, TRACE, which points from the original position
C      to the traced-back position (i.e. a vector in the direction
C      of the traceback).  The angle between this vector and the
C      central ray directional vector determines whether the path
C      length, L, should be increased or decreased.
C ----------------------------------------------------------------------
C
       TRACE(1) = -POS(3) * TTHET
       TRACE(2) = -POS(3) * TPHIT
       TRACE(3) = -POS(3)
C
C ----------------------------------------------------------------------
C      Determine differential path length, TRACE_LEN, corresponding
C      to the traceback.
C ----------------------------------------------------------------------
C
       IF(COMPUTE_PATH) THEN
         TRACE_LEN = ABS(POS(3)) * SQRT(1.D0+TTHET**2+TPHIT**2)
C
         IF (DOT(TRACE,C_RAY_ROT) .GE. 0.D0) THEN      
            VECT_ROT(5) = VECT(5) + TRACE_LEN
         ELSE
            VECT_ROT(5) = VECT(5) - TRACE_LEN
         ENDIF
       ELSE
         VECT_ROT(5) = VECT(5)    !use original path length
       ENDIF
C
       VECT_ROT(6) = VECT(6)      !delta doesn't change under rotations
C
       RETURN     
       END
