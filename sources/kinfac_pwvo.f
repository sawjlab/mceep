C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       Form kinematical weighting factors for response functions using
C       conventions of A. Picklesimer and J.W. Van Orden,
C       Phys. Rev. C 40, 290 (1989).
C
C       Note:  the factor V_T as given in that reference contains a
C       mistake (it should be tan^2 instead of tan).
C------------------------------------------------------------------------------
C
      SUBROUTINE KINFAC_PWVO(Q,W,THETA,V_L,V_T,V_TT,V_LT,V_LTP,V_TTP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      TTHT2 = TAN(THETA/2.D0)
      QMU2 = Q**2-W**2                !4-momentum transf.^2
      QRAT = QMU2/Q**2
C
      V_L   = QRAT*QRAT
      V_T   = QRAT/2.D0 + TTHT2*TTHT2
      V_TT  = QRAT/2.D0
      V_LT  = QRAT*SQRT(QRAT + TTHT2*TTHT2)
      V_LTP = QRAT*TTHT2
      V_TTP = TTHT2*SQRT(QRAT + TTHT2*TTHT2)
C
      RETURN
      END
