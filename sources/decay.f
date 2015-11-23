
C DECAY.FOR
C
C AUTHOR: G.M. Huber
C DATE:   93.02.15
C
C MODIFICATIONS:
C         11-JAN-1995   P.E. Ulmer
C         Implementation in MCEEP now weights the cross section by
C         the pion survival probability rather than applying a test
C         for survival (i.e. a cut).  Therefore, no longer set logical
C         variable LIVE.
C
C         07-Jul-1999  P. Markowitz
C         Adapted to kaon survival probability
C
C         12-AUG-1999  P.E. Ulmer
C         Made general for any unstable particle by passing lifetime
C         to routine.
C
C PURPOSE:
C calculates pion decay according to e^{-tof/tau}
C requires tof definition line in the hadron spectrometer stack
C
C VARIABLES:
C   INPUT:
C       lifetime      !lifetime of particle in rest frame (ns)
C       pf_p_i        !momentum of unstable particle (ejectile) (MeV/c)
C       tofp          !tof in hadron spectrometer               (ns)
C
C   OUTPUT:
C       survive_prob  !survival probability of ejectile
C
C ROUTINES CALLED:
C       none

      SUBROUTINE DECAY(LIFETIME)
C
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
C
      INCLUDE 'masses.cmn'
      INCLUDE 'var.cmn'
C
      DOUBLE PRECISION MEANTIME,SURVIVE_PROB,TOF_REL_TMP,TOFE,TOFP
      DOUBLE PRECISION EF_P_I,GAMMA,LIFETIME
C
      COMMON /SURVIVE_D/ SURVIVE_PROB
      COMMON /TOF/  TOF_REL_TMP,TOFE,TOFP
C
      IF (TOFP.EQ.0.) RETURN
C
C calculate time dilated mean lifetime
C
      EF_P_I   = SQRT(PF_P_I**2 + EJECT_MASS**2)
      GAMMA    = EF_P_I/EJECT_MASS
      MEANTIME = LIFETIME*GAMMA            ! mean lifetime of ejectile (ns)
C
C ejectile survival probability
C
      SURVIVE_PROB = EXP(-TOFP/MEANTIME)
C
      RETURN
      END
