C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       Form kinematical weighting factors for response functions
C       (conventions of V. Dmitrasinovic and Franz Gross:
C       Phys. Rev. C40, 2479 (1989) - see Equations 95 and 96).
C------------------------------------------------------------------------------
C
      SUBROUTINE KINFAC_2HGROSS(QVEC,OMEG,THETA_SCAT,
     #          S_T,S_LT,S_TP,S_LTP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ROOT2 = 1.414213562D0)
C
      XI    = QVEC*TAN(THETA_SCAT/2.D0)/SQRT(QVEC*QVEC-OMEG*OMEG)
      XI2   = XI*XI
      RADCL = SQRT(1.D0+XI2)
C
      S_T   =  0.5D0 + XI2
      S_LT  = -RADCL/ROOT2
      S_TP  =  XI*RADCL
      S_LTP = -XI/ROOT2
C
      RETURN
      END
