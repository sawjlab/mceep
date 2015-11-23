C
C---------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE OFF_SHELL_D_CC2
C
C       AUTHOR: P.E. Ulmer
C       DATE:   1-APR-1993
C       MODIFICATIONS:
C               March 2004 - Peter Monaghan; make subroutine 
C               compatible with MCEEP.
C       PURPOSE:
C               Calculates elementary OFF-SHELL (e,p) or (e,n)
C               cross-section using the "CC2" prescription given by DeForest
C               in Nucl. Phys.  A392 (1983)
C
C               The quantities returned are those parts of the cross
C               section which do not depend on the ejectile azimuthal angle.
C               The cross section can be reconstructed from these
C               three quantities by multiplying each by
C               either 1, cos(phi_x) or cos(2phi_x) and summing.  (PEU)
C
C       VARIABLES:
C       INPUT:
C       OUTPUT:
C
C       ROUTINES CALLED:
C               FORM_FACT
C------------------------------------------------------------------------------
C
      SUBROUTINE OFF_SHELL_D_CC2(E0_I,OMEGA,PF_P_I,PH_E_I,PH_P_I,
     #     TH_E_I,TH_P_I,PREC,SIGLPT,SIGLT,SIGTT)
C
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
C
      COMMON /PARTICLE/PROTON
      COMMON /FF_MOD/ FFTYPE
C
      LOGICAL PROTON
      CHARACTER*7 FFTYPE
      DOUBLE PRECISION KI(4),KF(4),PF(4),Q(4),QBAR(4),PBAR(4),PSUM(4)
C
      PARAMETER (HBARC = 197.3286)
      PARAMETER (ALPHA = 7.297E-3)            !Fine structure constant
      PARAMETER (KAPPA_P = 1.793)             !Proton anomalous moment
      PARAMETER (KAPPA_N = -1.91)             !Neutron anomalous moment
C
      INCLUDE 'masses.cmn'
C
      DATA M_EFF /1./
C
C ----------------------------------------------------------------------
C       Set up kinematic factors etc.
C ----------------------------------------------------------------------
C
      EPROT  = SQRT(PF_P_I**2+EJECT_MASS**2)   !Ejectile energy
      E_BAR  = SQRT(PREC**2+EJECT_MASS**2)     !Bound nucleon energy
      OMEGA_BAR = EPROT - E_BAR
      PF_MAG = SQRT(EPROT**2-EJECT_MASS**2)   !redefine PF_MAG 
C
C ----------------------------------------------------------------------
C     Form components of electron initial 4-momentum vector.
C ----------------------------------------------------------------------
C
      KI(1) = 0.D0
      KI(2) = 0.D0
      KI(3) = E0_I
      KI(4) = E0_I
C
C ----------------------------------------------------------------------
C     Form components of electron final 4-momentum vector.
C ----------------------------------------------------------------------
C
      EF = E0_I - OMEGA
      KF(1) = EF*COS(TH_E_I)*SIN(PH_E_I)
      KF(2) = EF*SIN(TH_E_I)
      KF(3) = EF*COS(TH_E_I)*COS(PH_E_I)
      KF(4) = EF
C
C ----------------------------------------------------------------------
C     Form components of ejectile 4-momentum vector.
C ----------------------------------------------------------------------
C
      PF(1) =  PF_MAG*COS(TH_P_I)*SIN(PH_P_I)
      PF(2) =  PF_MAG*SIN(TH_P_I)
      PF(3) =  PF_MAG*COS(TH_P_I)*COS(PH_P_I)
      PF(4) =  EPROT
C
C ----------------------------------------------------------------------
C     Form components of momentum transfer 4-vector.
C     Form components of proton initial 4-momentum vector.
C ----------------------------------------------------------------------
C
      DO I=1,4
         Q(I)     = KI(I) - KF(I)
         QBAR(I) = Q(I)
         IF(I.EQ.4) QBAR(I) = OMEGA_BAR
         PBAR(I) = PF(I) - QBAR(I)
         PSUM(I) = PBAR(I) + PF(I)
      ENDDO
C
C ----------------------------------------------------------------------
C	Form 3-vector products.
C ----------------------------------------------------------------------
C
      KIKF  = 0.D0
      QVEC2 = 0.D0
      PFQ   = 0.D0
      DO I=1,3
         KIKF  = KIKF  + KI(I)*KF(I)
         QVEC2 = QVEC2 + Q(I)**2
         PFQ   = PFQ   + PF(I)*Q(I)
      ENDDO
      C_TH_SCAT = KIKF/(E0_I*EF)
      TH_SCAT   = ACOS(C_TH_SCAT)
      QVEC      = SQRT(QVEC2)
      C_TH_PQ   = PFQ/(PF_MAG*QVEC)
      TH_PQ     = ACOS(C_TH_PQ)
C
C ----------------------------------------------------------------------
C     Form dot products of 4-vectors.
C ----------------------------------------------------------------------
C
      PBAR_DOT_4_PF  = -DOT_4(PBAR,PF)
      PBAR_DOT_4_Q   = -DOT_4(PBAR,Q)
      PF_DOT_4_Q     = -DOT_4(PF,Q)
      Q2           = -DOT_4(Q,Q)
      QBAR_DOT_4_Q   = -DOT_4(QBAR,Q)
      PSUM_DOT_4_Q   = -DOT_4(PSUM,Q)
C
      EBP    = E_BAR*EPROT
      Q4_MAG = SQRT(Q2)
      QQ     = Q2/QVEC2
      TAN_E2 = TAN(TH_SCAT/2.D0)**2
      Q_MU = SQRT(Q2)
C
C ----------------------------------------------------------------------
C       Get form factors.
C ----------------------------------------------------------------------
C
      CALL FORM_FACT(FFTYPE,Q_MU,M_EFF,GEP,GMP,GEN,GMN,GAP,GAN,
     #                       F1P,F2P,F1N,F2N)
      IF(PROTON)THEN
        F1 = F1P
        F2 = F2P
        KAPPA = KAPPA_P
      ELSE
        F1 = F1N
        F2 = F2N
        KAPPA = KAPPA_N
      ENDIF
      F2_FAC   = (KAPPA*F2)**2/(4.D0*EJECT_MASS**2)
      F1F2_FAC = F1*KAPPA*F2
C
C ----------------------------------------------------------------------
C	Form response functions.
C ----------------------------------------------------------------------
C
      W_C = (1.D0/EBP)*( (EBP+0.5D0*(PBAR_DOT_4_PF+EJECT_MASS**2))*F1**2
     #      - 0.5D0*QVEC2*F1F2_FAC
     #      - ( (PBAR_DOT_4_Q*EPROT+PF_DOT_4_Q*E_BAR)*OMEGA
     #          -EBP*Q2+PBAR_DOT_4_Q*PF_DOT_4_Q
     #          -0.5D0*(PBAR_DOT_4_PF-EJECT_MASS**2)*QVEC2 ) * F2_FAC )
      W_T = (1.D0/EBP)*( -(PBAR_DOT_4_PF+EJECT_MASS**2)*F1**2
     #      + QBAR_DOT_4_Q*F1F2_FAC
     #      + (2.D0*PBAR_DOT_4_Q*PF_DOT_4_Q
     #      -  (PBAR_DOT_4_PF-EJECT_MASS**2)*Q2)*F2_FAC )
      W_S = ((PF_MAG*SIN(TH_PQ))**2/EBP)*(F1**2+Q2*F2_FAC)
      W_I = (PF_MAG*SIN(TH_PQ)/EBP)*( -(E_BAR+EPROT)*F1**2
     #      + (PSUM_DOT_4_Q*OMEGA-(E_BAR+EPROT)*Q2)*F2_FAC )
C
C ----------------------------------------------------------------------
C       Calculate Off-shell cross-section.
C       Here, the Mott cross section is in MeV^-2 sr^-1.
C ----------------------------------------------------------------------
C
      SIGMA_MOTT = 4.D0*((ALPHA*COS(TH_SCAT/2.D0)*EF)/Q2)**2
      K_FACT = EPROT*PF_MAG
      FACTOR = K_FACT*SIGMA_MOTT*HBARC**2
C
C ----------------------------------------------------------------------
C       Collect terms multiplying COS(nPHI_X), n=(0,1,2).
C
C         K*SIGMA_EP is
C               SIGLPT + SIGLT*COS(PHI_X) + SIGTT*COS(2PHI_X)
C
C       Units for cross section terms are (fm^2/sr)*MeV^2.
C ----------------------------------------------------------------------
C
      SIGLPT = FACTOR*(QQ**2*W_C + (QQ/2.D0+TAN_E2)*(W_T+W_S))
      SIGLT  = FACTOR*(QQ*SQRT(QQ+TAN_E2)*W_I)
      SIGTT  = FACTOR*(QQ*W_S/2.D0)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C     Function DOT_4(A,B)
C
C     Purpose:
C
C       Takes dot product of two 4-vectors.
C       Metric:  DOT_4(A,A) = A_time**2 - A_space .dot. A_space
C------------------------------------------------------------------------------
C
      FUNCTION DOT_4(A,B)
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
      DOUBLE PRECISION A(4),B(4)
      SPACE = 0.D0
      DO I = 1,3
         SPACE = SPACE + A(I)*B(I)
      ENDDO
      TIME = A(4)*B(4)
      DOT_4 = TIME - SPACE
      RETURN
      END


