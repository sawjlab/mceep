C
C ----------------------------------------------------------------------
C     SUBROUTINE       TRACK
C     AUTHOR:          M.D. Nozar and P.E. Ulmer      
C     DATE:            February 29, 1992
C     PURPOSE:         
C                      Reconstruct the trajectory from the wire chamber
C                      positional coordinates obtained with the routine
C                      GET_VECT.
C                      These coordinates correspond to the intersection
C                      of the ray with each of the wire chambers
C                      in the detector stack.
C                      
C                      The Transport vector of the particle is computed
C                      at the location just prior to the TRK element in
C                      the input file.
C
C      NOTE:           This routine is specifically for a two-chamber
C                      system (where each chamber consists of a pair
C                      of orthogonal wire planes).
C
C ----------------------------------------------------------------------
C
       SUBROUTINE TRACK(IARM,TVEC_RECONST)
C
       IMPLICIT          NONE
C
       INCLUDE 'wc.cmn'
C
       INTEGER  IARM
C
       DOUBLE PRECISION  VECT_CH_1(3),VECT_CH_2(3)
       DOUBLE PRECISION  VECT_TR_1(3),VECT_TR_2(3),TVEC_RECONST(6)
       DOUBLE PRECISION  DELTA_X,DELTA_Y,DELTA_Z
       DOUBLE PRECISION  ROT_MAT(3,3),T_THETA,T_PHI
C
C ----------------------------------------------------------------------
C      Form positon vectors for the intersection of the ray with
C      each of the two chambers.
C ----------------------------------------------------------------------
C
       VECT_CH_1(1) = TVEC_CHMBR(1,1)
       VECT_CH_1(2) = TVEC_CHMBR(1,3)
       VECT_CH_1(3) = 0.D0
C
       VECT_CH_2(1) = TVEC_CHMBR(2,1)
       VECT_CH_2(2) = TVEC_CHMBR(2,3)
       VECT_CH_2(3) = 0.D0
C
C ----------------------------------------------------------------------
C      Determine the coordinates of the "chamber" position vectors in
C      the Transport system.
C
C      These are simple passive rotations (with no traceback since the
C      orientation of the ray has not been determined yet).
C ----------------------------------------------------------------------
C
       CALL ROTATE_Y(WC_C_CANGLE(IARM,1),-WC_C_SANGLE(IARM,1),ROT_MAT) 
       CALL MAT_MULT(3,3,1,ROT_MAT,VECT_CH_1,VECT_TR_1)
C
       CALL ROTATE_Y(WC_C_CANGLE(IARM,2),-WC_C_SANGLE(IARM,2),ROT_MAT)
       CALL MAT_MULT(3,3,1,ROT_MAT,VECT_CH_2,VECT_TR_2)
C
C ----------------------------------------------------------------------
C      Add the wire chamber locations to the z-coordinates so that both
C      vectors refer to the "initial" (i.e. just before the TRK element)
C      Transport system.
C ----------------------------------------------------------------------
C
       VECT_TR_1(3) = VECT_TR_1(3) + WC_LOC(IARM,1)
       VECT_TR_2(3) = VECT_TR_2(3) + WC_LOC(IARM,2)
C
       DELTA_X = VECT_TR_2(1) - VECT_TR_1(1)
       DELTA_Y = VECT_TR_2(2) - VECT_TR_1(2)
       DELTA_Z = VECT_TR_2(3) - VECT_TR_1(3)
C
       T_THETA = DELTA_X/DELTA_Z              !Tangent THETA
       T_PHI   = DELTA_Y/DELTA_Z              !Tangent PHI
C
C ----------------------------------------------------------------------
C      Reconstruct the Transport vector.  For X and Y this involves
C      tracing back along the particle trajectory from the intersection
C      point with the first wire chamber to the "initial" Transport XY
C      plane.  Here, again, "initial" refers to the location just prior
C      to the TRK element in the input file.
C
C      The L coordinate (path length) should correspond to the "initial"
C      Transport coordinate system.  Since the first wire chamber is
C      at WC_LOC(1) relative to the "initial" Transport system, we
C      must subtract from L the path length traversed between these
C      two systems.
C ----------------------------------------------------------------------
C
       TVEC_RECONST(1) = VECT_TR_1(1)-VECT_TR_1(3)*T_THETA
       TVEC_RECONST(2) = 1.D3*ATAN(T_THETA)     !Transport THETA (mr)
       TVEC_RECONST(3) = VECT_TR_1(2)-VECT_TR_1(3)*T_PHI
       TVEC_RECONST(4) = 1.D3*ATAN(T_PHI)       !Transport PHI   (mr)
       TVEC_RECONST(5) = TVEC_CHMBR(1,5)
     #          - WC_LOC(IARM,1)*(SQRT(1.D0+T_THETA**2+T_PHI**2)-1.D0)
       TVEC_RECONST(6) = TVEC_RECONST(1)/X_DELTA(IARM)
C
       RETURN
       END
