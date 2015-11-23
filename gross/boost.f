C
C_______________________________________________________________________
C
C*****boost
C
C            This subroutine produces Gross deuteron wave functions
C         which have been given an arbitrary boost along the z-axis.
C
C         INPUT
C              pp.........momentum in boosted frame in MeV
C              thetap.....polar angle in boosted frame in radians
C              shztaz.....hyperbolic sine of rapidity of boost
C              chztaz.....hyperbolic cosine of rapidity of boost
C              itype=0....particle 2 on shell
C                   =1....particle 1 on shell
C              irel=1.....relativistic boosted wave function
C                  =0.....nonrelativistic wave function
C
C         OUTPUT
C           aboost(lnuc1,lnuc2,ienrgy,ldeu)......array containing boosted
C                                                deuteron wave functions.
C
C             lnuc1=0.........positive nucleon helicity for particle 1
C                  =1.........negative nucleon helicity for particle 1
C             lnuc2=0.........positive nucleon helicity for particle 2
C                  =1.........negative nucleon helicity for particle 2
C             ienrgy=0........positive energy projection
C                   =1........negative energy projection
C             ldeu=-1,0,1....deuteron helicity
C
C
      SUBROUTINE BOOST(PP,THETAP,SHZTAZ,CHZTAZ,ITYPE,IREL,ABOOST)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION CHI(0:1),ABOOST(0:1,0:1,0:1,-1:1),D1(-1:1,-1:1),
     &          DHALFP(0:1,0:1),DHALFM(0:1,0:1),A(0:1,0:1,0:1)
C
      DOUBLE PRECISION MN,MT
      COMMON/DPARA/MN,MT
C
      DATA PI/3.141592653589793238/
C
C             calculate deuteron rest frame polar angle and momentum
C
      KAPPA=2*ITYPE-1
C
      SNTHP=SIN(THETAP)
      CSTHP=COS(THETAP)
C
      IF(IREL.EQ.1) THEN                ! relativistic kinematics
C
        SHZTAP=PP/MN
        CHZTAP=SQRT(SHZTAP**2+1.)
C
        TNTH=KAPPA*SHZTAP*SNTHP/(KAPPA*CHZTAZ*SHZTAP*CSTHP-SHZTAZ*
     &       CHZTAP)
        THETA=ATAN(TNTH)+(1.-SIGN(1.D0,TNTH))*PI/2.
C
        CHZTA=CHZTAZ*CHZTAP-KAPPA*SHZTAZ*SHZTAP*CSTHP
        SHZTA=SQRT(CHZTA**2-1.)
        P=MN*SHZTA
C              calculate the Wigner rotation angles
C
C
        DO I=0,1
          KP=1-2*I
          SGNFAC=KP*SHZTAZ
C          tnchi=shztaz*sin(theta)/(shztaz*chzta
C     &          +kp*cos(theta)*chztaz*shzta)
          TNCHI=SHZTAZ*SIN(THETA)/(KP*SHZTA*CHZTAZ
     &          +SHZTAZ*CHZTA*COS(THETA))
          CHI(I)=ATAN(TNCHI)
          IF(SGNFAC.GT.0.0.AND.TNCHI.LT.0.0) THEN
            CHI(I)=CHI(I)+PI
          ELSEIF(SGNFAC.LT.0.0.AND.TNCHI.GT.0.0) THEN
            CHI(I)=CHI(I)-PI
          ENDIF
        END DO
C
      ELSE IF(IREL.EQ.0) THEN            ! nonrelativistic kinematics
C
        PUC=MT*SHZTAZ
        P=SQRT(PP**2-KAPPA*PUC*PP*CSTHP+0.25*PUC**2)
        CSTH=(PP*CSTHP-0.5*KAPPA*PUC)/P
        THETA=ACOS(CSTH)
        CHI(0)=0.0
        CHI(1)=0.0
C
      END IF
C
C               get deuteron wave functions in rest frame
C
      CALL DEUWAV(P,IREL,A)
C
C                     get spin one rotation matrix
C
      CALL ROTONE(THETA,D1)
C
C           get spin one half rotation matrices for two Wigner angles
C
      CALL ROTHLF(CHI(0),DHALFP)
      CALL ROTHLF(CHI(1),DHALFM)
C
C                construct boosted wave functions
C
      IF(ITYPE.EQ.0) THEN
        DO LDEU=-1,1
          DO IENRGY=0,1
            DO LNUC1=0,1
              DO LNUC2=0,1
                ABOOST(LNUC1,LNUC2,IENRGY,LDEU)=0.0
                DO LP1=0,1
                  DO LP2=0,1
                    LP1M2=LP2-LP1
                    ABOOST(LNUC1,LNUC2,IENRGY,LDEU)=
     &                    ABOOST(LNUC1,LNUC2,IENRGY,LDEU)+
     &                    DHALFP(LNUC1,LP1)*D1(LDEU,LP1M2)*
     &                    DHALFM(LP2,LNUC2)*A(LP2,LP1,IENRGY)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      ELSE IF(ITYPE.EQ.1) THEN
        DO LDEU=-1,1
          DO IENRGY=0,1
            DO LNUC1=0,1
              DO LNUC2=0,1
                ABOOST(LNUC1,LNUC2,IENRGY,LDEU)=0.0
                DO LP1=0,1
                  DO LP2=0,1
                    LP1M2=LP2-LP1
                    ABOOST(LNUC1,LNUC2,IENRGY,LDEU)=
     &                    ABOOST(LNUC1,LNUC2,IENRGY,LDEU)+
     &                    DHALFP(LNUC1,LP1)*D1(LDEU,LP1M2)*
     &                    DHALFM(LP2,LNUC2)*A(LP1,LP2,IENRGY)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      ENDIF
C
C
      RETURN
      END
