C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE OFF_SHELL_D
C
C       AUTHOR: R.W. Lourie
C       DATE:   V1.0    18-AUG-1986
C       MODIFICATIONS:
C               none
C       PURPOSE:
C               Calculates elementary OFF-SHELL (e,p) or (e,n)
C               cross-section using the prescription given by DeForest in
C               Nucl. Phys.  A392 (1983)
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
      SUBROUTINE OFF_SHELL_D(EFIN,Q2,QMU2,EPROT,PREC,THEPQ,
     #     SCAT_ANG,SIGLPT,SIGLT,SIGTT)
C
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
C
      COMMON /PARTICLE/PROTON
      COMMON /FF_MOD/ FFTYPE
C
      LOGICAL PROTON
      CHARACTER*7 FFTYPE
C
      PARAMETER (PI = 3.14159)
      PARAMETER (HBARC = 197.3286)
      PARAMETER (ALPHA = 7.297E-3)            !Fine structure constant
      PARAMETER (KAPPA_P = 1.793)             !Proton anomalous moment
      PARAMETER (KAPPA_N = -1.91)             !Neutron anomalous moment
C
      INCLUDE 'masses.cmn'
C
      DATA M_EFF /1./
C
C       Set up kinematic factors etc.
C
      PF_MAG = SQRT(EPROT**2-EJECT_MASS**2)
      E_BAR = SQRT(PREC**2+EJECT_MASS**2)     !Initial ejectile energy
      OMEGA_BAR = EPROT - E_BAR
      Q_MU_BAR = SQRT(Q2 - OMEGA_BAR**2)
      EBP = E_BAR*EPROT
      X = Q_MU_BAR**2/4./EJECT_MASS**2
      Q_MU = SQRT(QMU2)
      QQ = QMU2/Q2
      TAN_E2 = TAN(SCAT_ANG/2.)**2
C
C       Put together various response functions
C
      CALL FORM_FACT(FFTYPE,Q_MU,M_EFF,GEP,GMP,GEN,GMN,GAP,GAN,
     #                       F1P,F2P,F1N,F2N)
C
      IF(PROTON)THEN
        F1 = F1P
        F2 = F2P
        KAPPA = KAPPA_P
      ELSE
        F1 = F1N
        F2 = F2N
        KAPPA = KAPPA_N
      ENDIF
      FSUMSQ = F1**2 + KAPPA**2*X*F2**2
      FSUM   = F1    + KAPPA*F2
C
      W_C = ((E_BAR+EPROT)**2*FSUMSQ-Q2*FSUM**2)/(4.*EBP)
      W_T = FSUM**2*Q_MU_BAR**2/2./EBP
      W_S = FSUMSQ*PF_MAG**2*(SIN(THEPQ))**2/EBP
      W_I = -(E_BAR+EPROT)*FSUMSQ*PF_MAG*SIN(THEPQ)/EBP
C
C       Calculate Off-shell cross-section.
C       Here, the Mott cross section is in MeV^-2 sr^-1.
C
      SIGMA_MOTT = 4.*((ALPHA*COS(SCAT_ANG/2.)*EFIN)/QMU2)**2
      K_FACT = EPROT*PF_MAG
      FACTOR = K_FACT*SIGMA_MOTT*HBARC**2
C
C       Collect terms multiplying COS(nPHI_X), n=(0,1,2).
C
C         K*SIGMA_EP is
C               SIGLPT + SIGLT*COS(PHI_X) + SIGTT*COS(2PHI_X)
C
C       Units for cross section terms are (fm^2/sr)*MeV^2.
C
      SIGLPT = FACTOR*(QQ**2*W_C + (QQ/2.+TAN_E2)*(W_T+W_S))
      SIGLT  = FACTOR*(QQ*SQRT(QQ+TAN_E2)*W_I)
      SIGTT  = FACTOR*(QQ*W_S/2.)
C
      RETURN
      END
