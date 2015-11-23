C
C -----------------------------------------------------------------------
C     This subroutine produces cross sections for deuteron
C     electrodisintegration in Impulse Approximation using the
C     response functions of V. Dmitrasinovic and Franz Gross
C     from Phys. Rev. C 40, 2479 (1989).
C     (See Equation 95 in that reference).
C -----------------------------------------------------------------------
C
      SUBROUTINE DEEP_GROSS(QMU2,CTHEPQ,POL_BEAM,SIGMA,ASYM,
     #                      PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION RESP(0:8,0:3,9)
C
      DOUBLE PRECISION MN,MT
      COMMON/DPARA/MN,MT
C
      INCLUDE 'var.cmn'
C
      DATA MN/938.2796D0/
      DATA MT/1875.63D0/
      DATA PI/3.141592653589793238D0/
C
      W_RESPNS = W    !Make sure that W calc. in KINEM is not altered
C
C -----------------------------------------------------------------------
C       Convert beam polarization (from -1 to +1) to electron helicity
C       in the Dmitrasinovic/Gross convention (from -1/2 to + 1/2).
C -----------------------------------------------------------------------
C
      E_HEL_GROSS = POL_BEAM/2.D0
C
      EP = SQRT(PF_P_I**2+MN**2)                !proton energy in the lab
C
      IFF=4      !Choose dipole form factors (Galster GEn)
      ILAB=1     !get laboratory frame response functions
                 ! --> no need to perform Wigner rotations on pol.
                 !     components due to boost along Q.
C
C -----------------------------------------------------------------------
C       Get response functions
C -----------------------------------------------------------------------
C
      CALL RESPNS(PF_P_I,QMU2,THETA_PQ,W_RESPNS,IFF,ILAB,RESP)
C
C -----------------------------------------------------------------------
C       Get kinematical factors and Mott cross section
C -----------------------------------------------------------------------
C
      RECOIL = (W/MT)*(1.D0/(1.D0 +
     #          (OMEGA*PF_P_I-EP*QMAG*CTHEPQ)/(MT*PF_P_I)))
      CALL KINFAC_2HGROSS(QMAG,OMEGA,TSCAT,S_T,S_LT,S_TP,S_LTP)
      CALL SMOTT(TSCAT,PF_E_I,QMU2,SIG_MOTT)
C
      SIGFAC = SIG_MOTT*PF_P_I*RECOIL*QMU2/(4.D0*PI*MT*QMAG*QMAG)
C
C -----------------------------------------------------------------------
C       Build the cross section, electron analyzing power (ASYM)
C       and polarizations.  In this (PWIA) calculation ASYM=0.
C -----------------------------------------------------------------------
C
      SIGMA = SIGFAC*(       RESP(0,0,1) + S_T*RESP(0,0,2)
     &          - 0.5D0*     RESP(0,0,3)*COS(2.D0*PHI_X)
     &          + S_LT *     RESP(0,0,4)*COS(PHI_X))
      ASYM  = SIGFAC*(S_LTP* RESP(0,0,5)*SIN(PHI_X)*2.D0*E_HEL_GROSS)
C
      PN_HI    = SIGFAC*(    RESP(0,1,1) + S_T*RESP(0,1,2)
     &          - 0.5D0*     RESP(0,1,3)*COS(2.D0*PHI_X)
     &          + S_LT *     RESP(0,1,4)*COS(PHI_X) )
      PT_HI    = SIGFAC*(-0.5D0*RESP(0,2,6)*SIN(2.D0*PHI_X)
     &          + S_LT*         RESP(0,2,7)*SIN(PHI_X) )
      PL_HI    = SIGFAC*(-0.5D0*RESP(0,3,6)*SIN(2.D0*PHI_X)
     &          + S_LT*         RESP(0,3,7)*SIN(PHI_X) )
C
      PN_HD    = SIGFAC*(S_LTP*RESP(0,1,5)*SIN(PHI_X)*2.D0*E_HEL_GROSS)
      PT_HD    = SIGFAC*( S_TP*RESP(0,2,9) * 2.D0*E_HEL_GROSS
     &          + S_LTP*       RESP(0,2,8)*COS(PHI_X)*2.D0*E_HEL_GROSS)
      PL_HD    = SIGFAC*( S_TP*RESP(0,3,9) * 2.D0*E_HEL_GROSS
     &          + S_LTP*       RESP(0,3,8)*COS(PHI_X)*2.D0*E_HEL_GROSS)
C
      RETURN
      END












