C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C     SUBROUTINE     SPIN_1ST_ORDER
C     AUTHOR:        M.D. Nozar and P.E. Ulmer
C     DATE:          22-JUL-1991
C     PURPOSE:
C                    Produce a spin rotation matrix to obtain the 
C                    components of the rotated spin in the Transport
C                    system after the particle passes through a
C                    magnetic element.  This routine
C                    determines the spin precession angle from the
C                    angle between the momentum vectors before and
C                    after the magnetic element.  The angle between
C                    the entrance and exit momentum vectors is deduced
C                    via a 1st order Transport matrix.
C
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_1ST_ORDER(BEND_0,PMAT_EL)
C
      IMPLICIT           NONE
C
      COMMON /PHYSVAR/ IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
C
      INCLUDE           'input.cmn'
      INCLUDE           'spectrometer.cmn'
      INCLUDE           'masses.cmn'
C
      DOUBLE PRECISION   P1(3),P2(3),P2_TMP(3),P1_X_P2(3),P1_X_P2_MAG
      DOUBLE PRECISION   AXIS(3),TVEC_TMP(6),PMAT_EL(3,3)
      DOUBLE PRECISION   BEND_0,S_BEND_0,C_BEND_0
      DOUBLE PRECISION   TH_P1_P2,C_TH_P1_P2
      DOUBLE PRECISION   TH_AXIS,S_TH_AXIS,C_TH_AXIS
      DOUBLE PRECISION   PH_AXIS,S_PH_AXIS,C_PH_AXIS
      DOUBLE PRECISION   CHI,R_ANGLE,G_FACTOR,GAMMA
      DOUBLE PRECISION   RMAT_Y(3,3),RMAT_Z(3,3)
      DOUBLE PRECISION   M_MULT1(3,3),M_MULT2(3,3),M_MULT3(3,3)
      DOUBLE PRECISION   M_MULT4(3,3)
      DOUBLE PRECISION   T_TH_1,T_PH_1,P1_NORM,P_MOM
      DOUBLE PRECISION   T_TH_2,T_PH_2,P2_NORM
      DOUBLE PRECISION   DOT,PI
C
      INTEGER            IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL,I,J
C
      PARAMETER          (PI  = 3.141592653590D0)
C
      DATA G_FACTOR /5.58667D0/   !gyromagnetic ratio (proton)
      IF(IELAST_OPT .EQ. 20) G_FACTOR = 0.8574D0   !deuteron
C
C ----------------------------------------------------------------------
C     Calculate particle momentum vector before magnetic element
C     from the current Transport vector.
C ----------------------------------------------------------------------
C
      P_MOM   = PF_P * (1.D0 + TVEC(6,2)/100.D0)      !mag. of momentum
      T_TH_1  = TAN(TVEC(2,2) * 1.D-3)   !tangent of Transport Theta
      T_PH_1  = TAN(TVEC(4,2) * 1.D-3)   !tangent of Transport Phi
      P1_NORM = P_MOM/SQRT(1.D0 + T_TH_1**2 + T_PH_1**2) 
      P1(1)   = P1_NORM * T_TH_1         !form components
      P1(2)   = P1_NORM * T_PH_1          !of momentum vector
      P1(3)   = P1_NORM * 1.D0             !in Transport system
C
C ----------------------------------------------------------------------
C     Calculate particle momentum vector after magnetic element.
C     First, obtain the Transport vector after the element by
C     multiplication of the current vector by a 1st order Transport
C     matrix.
C ----------------------------------------------------------------------
C
      DO I=1,6                      !get Transport vector after element
         TVEC_TMP(I) = 0.
         DO J=1,6
            TVEC_TMP(I) = TVEC_TMP(I)
     #                   + TVEC(J,2)*TMATRIX(2,NPOL_EL,I,J)           
         ENDDO
      ENDDO
C
      T_TH_2    = TAN(TVEC_TMP(2)*1.D-3)    !tangent of Transport Theta
      T_PH_2    = TAN(TVEC_TMP(4)*1.D-3)    !tangent of Transport Phi
      P2_NORM   = P_MOM/SQRT(1.D0+T_TH_2**2+T_PH_2**2) 
      P2_TMP(1) = P2_NORM*T_TH_2     ! P2_TMP is the momentum vector
      P2_TMP(2) = P2_NORM*T_PH_2      ! after the element relative to
      P2_TMP(3) = P2_NORM*1.D0         ! the new Transport system.
C 
C ----------------------------------------------------------------------
C     Find the components of the exit momentum vector in the
C     entrance Transport system (where exit and entrance are defined
C     relative to the current magnetic element) by performing a
C     passive (i.e. coordinate) rotation about the Transport Y axis
C     by the bend angle for this element.  P2 is then the momentum
C     vector after the element expressed in the same coordinate system
C     as for P1.  Note that, no matter what the orientation of the
C     spectrometer bend plane, the sense of rotation to get back to
C     the entrance Transport system is always positive.
C ----------------------------------------------------------------------
C   
      C_BEND_0 = COS(BEND_0)
      S_BEND_0 = SIN(BEND_0)
      CALL ROTATE_Y(C_BEND_0,S_BEND_0,RMAT_Y)   
      CALL MAT_MULT(3,3,1,RMAT_Y,P2_TMP,P2)           
C 
C ----------------------------------------------------------------------
C     Find the components of unit vector AXIS (rotation axis for spin
C     and momentum).
C ----------------------------------------------------------------------
C   
      CALL VEC_CROSS(P1,P2,P1_X_P2)
      P1_X_P2_MAG = SQRT(DOT(P1_X_P2,P1_X_P2))
C
      IF(P1_X_P2_MAG .NE. 0.D0) THEN
        AXIS(1) = P1_X_P2(1)/P1_X_P2_MAG
        AXIS(2) = P1_X_P2(2)/P1_X_P2_MAG
        AXIS(3) = P1_X_P2(3)/P1_X_P2_MAG
      ELSE
        AXIS(1) = P1(1)/P_MOM    !In this case direction of AXIS
        AXIS(2) = P1(2)/P_MOM     !is arbitrarily chosen along P1
        AXIS(3) = P1(3)/P_MOM      !since there is no spin rotation.
      ENDIF
C
C ---------------------------------------------------------------------
C     Obtain TH_AXIS and PH_AXIS, the spherical coordinate polar
C     and azimuthal angles of AXIS in the entrance (i.e. at entrance of 
C     magnetic element) transport system. 
C ----------------------------------------------------------------------
C  
      C_TH_AXIS = AXIS(3)
      TH_AXIS   = ACOS(C_TH_AXIS)      !polar angle of rotation axis
      S_TH_AXIS = SIN(TH_AXIS)
      C_PH_AXIS   = AXIS(1)/S_TH_AXIS
      PH_AXIS   = ACOS(C_PH_AXIS) !azimuthal angle of rot. axis
      IF (AXIS(2).LT.0.D0)  PH_AXIS = 2.D0 * PI - PH_AXIS
      S_PH_AXIS = SIN(PH_AXIS)
C
C ----------------------------------------------------------------------
C     Calculate the momentum bend angle (TH_P1_P2), the spin precession
C     angle (Chi) and the total spin rotation angle 
C     (R_ANGLE = TH_P1_P2 + CHI).
C ----------------------------------------------------------------------
C
      C_TH_P1_P2 = DOT(P1,P2)/P_MOM**2
      IF (C_TH_P1_P2 .GT. 1.D0 .AND. C_TH_P1_P2-1.D0 .LT. 1.E-12) THEN
        C_TH_P1_P2 = 1.D0      !protect against roundoff
      ELSEIF(C_TH_P1_P2.LT.-1.D0 .AND. C_TH_P1_P2+1.D0.GT. -1.E-12) THEN
        C_TH_P1_P2 = -1.D0     !protect against roundoff
      ENDIF
      TH_P1_P2 = ACOS(C_TH_P1_P2)
      GAMMA    = SQRT((P_MOM/EJECT_MASS)**2 + 1.D0)   !Lorentz factor
      CHI      = (G_FACTOR - 2.D0) * GAMMA * TH_P1_P2/2.D0
      R_ANGLE  = TH_P1_P2 + CHI
C
C ----------------------------------------------------------------------
C     Obtain the rotation matrix PMAT_EL by performing appropriate 
C     passive and active rotations.
C
C     First perform passive rotations about Z and Y axes to bring
C     the Transport Z axis along the spin rotation axis, "AXIS".
C     Then rotation about "AXIS" is equivalent to rotation about
C     one of the Transport coordinate axes (Z in particular).
C ----------------------------------------------------------------------
C      
      CALL ROTATE_Z(C_PH_AXIS,S_PH_AXIS,RMAT_Z)     
      CALL ROTATE_Y(C_TH_AXIS,S_TH_AXIS,RMAT_Y)
      CALL MAT_MULT(3,3,3,RMAT_Y,RMAT_Z,M_MULT1)
C
C ----------------------------------------------------------------------
C     Now, perform (active) rotation of spin about Z axis (i.e.
C     about "AXIS").
C ----------------------------------------------------------------------
C
      CALL ROTATE_Z(COS(R_ANGLE),-SIN(R_ANGLE),RMAT_Z)
      CALL MAT_MULT(3,3,3,RMAT_Z,M_MULT1,M_MULT2)
C
C ----------------------------------------------------------------------
C     Now, invert previous passive rotations to restore coordinate
C     system to the entrance Transport system.  The resulting
C     product matrix gives the spin rotation matrix in terms of
C     the entrance Transport coordinate system.
C ----------------------------------------------------------------------
C   
      CALL ROTATE_Y(C_TH_AXIS,-S_TH_AXIS,RMAT_Y)
      CALL MAT_MULT(3,3,3,RMAT_Y,M_MULT2,M_MULT3)
      CALL ROTATE_Z(C_PH_AXIS,-S_PH_AXIS,RMAT_Z)
      CALL MAT_MULT(3,3,3,RMAT_Z,M_MULT3,M_MULT4)
C
C ----------------------------------------------------------------------
C     Finally, perform passive rotation about Y axis to express
C     rotated spin vector in terms of exit Transport coordinate
C     system.  The product of the resulting matrix with the spin vector
C     at the entrance to the element gives the spin components
C     after the element in the exit Transport system.
C ----------------------------------------------------------------------
C
      CALL ROTATE_Y(C_BEND_0,-S_BEND_0,RMAT_Y)   
      CALL MAT_MULT(3,3,3,RMAT_Y,M_MULT4,PMAT_EL)
C
      RETURN
      END     
