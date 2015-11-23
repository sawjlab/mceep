C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE PWIA_VO
C
C       AUTHOR: P.E. Ulmer
C       DATE:   V1.0    04-SEPT-1990
C       PURPOSE:
C               Calculates "semirelativistic plane-wave impulse
C               approximation" response functions using prescription
C               given by A. Picklesimer and J.W. Van Orden in
C               Phys. Rev. C 40, 290 (1989).  (See Equation 3.18 in that
C               reference).
C
C               The quantities returned need to be multiplied by a
C               momentum distribution to obtain the response functions.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE PWIA_VO(EFIN,OMEGA,Q2,QMU2,PNUC_F,PREC,
     #     SCAT_ANG,THE_NQ,PHI_X,POL_BEAM,SIGMA,ASYM,
     #     PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD)
C
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
C
      COMMON /PARTICLE/PROTON
      COMMON /FF_MOD/ FFTYPE
C
      LOGICAL PROTON
      CHARACTER*7 FFTYPE
      PARAMETER (PI = 3.14159D0)
      PARAMETER (KAPPA_P = 1.793)             !Proton anomalous moment
      PARAMETER (KAPPA_N = -1.91)             !Neutron anomalous moment
C
      INCLUDE 'masses.cmn'
C
      DATA M_EFF /1.D0/
C
C --------------------------------------------------------------------------
C       Set up kinematic factors etc.
C --------------------------------------------------------------------------
C
      QMU = SQRT(QMU2)
      QVEC = SQRT(Q2)
      ENUC_F = SQRT(PNUC_F**2+EJECT_MASS**2)    !Final nucleon energy
      ENUC_I = SQRT(PREC**2+EJECT_MASS**2)      !Initial nucleon energy
      OMEGA_BAR = ENUC_F - ENUC_I
      QMU_BAR2 = Q2 - OMEGA_BAR**2
      TAU_BAR = QMU_BAR2/(4.D0*EJECT_MASS*EJECT_MASS)
      STHE_NQ = SIN(THE_NQ)
      CTHE_NQ = COS(THE_NQ)
      SPHIX = SIN(PHI_X)
      CPHIX = COS(PHI_X)
      S2PHIX = SIN(2.D0*PHI_X)
      C2PHIX = COS(2.D0*PHI_X)
C
C --------------------------------------------------------------------------
C       Get nucleon form factors
C --------------------------------------------------------------------------
C
      CALL FORM_FACT(FFTYPE,QMU,M_EFF,GEP,GMP,GEN,GMN,GAP,GAN,
     #                       F1P,F2P,F1N,F2N)
C
      IF(PROTON)THEN
        GE2 = GEP**2
        GM2 = GMP**2
        F1 = F1P
        F2 = F2P
        KAPPA = KAPPA_P
      ELSE
        GE2 = GEN**2
        GM2 = GMN**2
        F1 = F1N
        F2 = F2N
        KAPPA = KAPPA_N
      ENDIF
      F2 = KAPPA*F2           !Van Orden F2(0) = Kappa
      GM = SQRT(GM2)
      F1_F2 = F1**2 + TAU_BAR*F2**2
C
C --------------------------------------------------------------------------
C       Construct response functions
C --------------------------------------------------------------------------
C
      RL     = F1_F2*((ENUC_I+ENUC_F)/(2.D0*EJECT_MASS))**2
     #          - GM2*Q2/(4.D0*EJECT_MASS**2)
      RL_N   = 0.D0
      RT     = F1_F2*(PNUC_F*STHE_NQ/EJECT_MASS)**2 + 2.D0*TAU_BAR*GM2
      RT_N   = 0.D0
      RTT    = F1_F2*(PNUC_F*STHE_NQ/EJECT_MASS)**2
      RTT_N  = 0.D0
      RTT_T  = 0.D0
      RTT_L  = 0.D0
      RLT    = -F1_F2*(ENUC_I+ENUC_F)*PNUC_F*STHE_NQ/EJECT_MASS**2
      RLT_N  = 0.D0
      RLT_T  = 0.D0
      RLT_L  = 0.D0
      RLTP   = 0.D0
      RLTP_N = (F1+F2*(ENUC_F*OMEGA_BAR - QVEC*PNUC_F*CTHE_NQ)
     #                  /(2.D0*EJECT_MASS**2))*GM*QVEC/EJECT_MASS
      RLTP_T = -(F1*CTHE_NQ + F2*(ENUC_F*OMEGA_BAR*CTHE_NQ-QVEC*
     #                  PNUC_F)/(2.D0*EJECT_MASS**2))*GM*QVEC/EJECT_MASS
      RLTP_L = -(F1*ENUC_F/EJECT_MASS + F2*OMEGA_BAR/(2.D0*EJECT_MASS))
     #                  *GM*QVEC*STHE_NQ/EJECT_MASS
      RTTP_T = (F1*OMEGA_BAR/EJECT_MASS
     #     - F2*ENUC_F*QMU_BAR2/(2.D0*EJECT_MASS**3))*GM*STHE_NQ
      RTTP_L = (F1*(PNUC_F*QVEC-OMEGA_BAR*ENUC_F*CTHE_NQ)/EJECT_MASS**2
     #                  + F2*QMU_BAR2*CTHE_NQ/(2.D0*EJECT_MASS**2))*GM
C
C --------------------------------------------------------------------------
C       Get kinematic factors and Mott cross section.
C --------------------------------------------------------------------------
C
      CALL KINFAC_PWVO(QVEC,OMEGA,SCAT_ANG,V_L,V_T,V_TT,V_LT,V_LTP,
     #          V_TTP)
      CALL SMOTT(SCAT_ANG,EFIN,QMU2,SIG_MOTT)
      SIGFAC = SIG_MOTT*EJECT_MASS*PNUC_F/(16.D0*PI**3)
C
C --------------------------------------------------------------------------
C       A word about normalizations:
C
C       Van Orden's scaler momentum distribution (see Equation 3.16
C       in the reference cited above) is normalized so that:
C              integral of 1/(2pi)^3 n_nlj(p) d^3p = 1
C       Also the momentum distribution appearing in the response
C       functions is n^S_nlj = m/E n_nlj
C       Since the MCEEP momentum distributions are normalized so that:
C              integral of phi(p)^2 d^3p = 1
C       we need to multiply the Van Orden response functions by
C       a factor of m/E (2pi)^3
C
C       But wait, there's more:
C
C       Van Orden normalizes his cross section so that when SUMMED
C       over the two polarization states (up and down) for any
C       orientation (N, L or T) it gives the unpolarized cross section.
C       Thus, if we want the unpolarized cross section we must sum over
C       the five unpolarized response functions and multiply by a
C       factor of TWO.
C
C       Furthermore, the polarization components should also be
C       multiplied by 2 since:
C               Sigma_VO = 1/2 Sigma_0 (1 + P dot S)
C       where Sigma_0 is the unpolarized cross section and P is the
C       polarization vector defined as:
C               (Sigma_N/Sigma_0 , Sigma_L/Sigma_0 , Sigma_T/Sigma_0).
C --------------------------------------------------------------------------
C
      SIGFAC = SIGFAC*2.D0*(2.D0*PI)**3*EJECT_MASS/ENUC_I
C
C --------------------------------------------------------------------------
C       Form cross section and polarizations.  Also form electron
C       analyzing power (ASYM) which is zero in this (PWIA) calculation.
C --------------------------------------------------------------------------
C
      SIGMA = SIGFAC*(V_L*RL + V_T*RT + V_TT*RTT*C2PHIX
     #          + V_LT*RLT*CPHIX)
      ASYM  = SIGFAC*(POL_BEAM*V_LTP*RLTP*SPHIX)
C
      PN_HI = SIGFAC*(V_L*RL_N + V_T*RT_N + V_TT*RTT_N*C2PHIX
     #          + V_LT*RLT_N*CPHIX )
      PT_HI = SIGFAC*(V_TT*RTT_T*S2PHIX + V_LT*RLT_T*SPHIX)
      PL_HI = SIGFAC*(V_TT*RTT_L*S2PHIX + V_LT*RLT_L*SPHIX)
C
      PN_HD = SIGFAC*POL_BEAM*V_LTP*RLTP_N*SPHIX
      PT_HD = SIGFAC*(POL_BEAM*V_LTP*RLTP_T*CPHIX 
     #          + POL_BEAM*V_TTP*RTTP_T)
      PL_HD = SIGFAC*(POL_BEAM*V_LTP*RLTP_L*CPHIX 
     #          + POL_BEAM*V_TTP*RTTP_L)
C
      RETURN
      END
