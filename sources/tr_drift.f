C
C ----------------------------------------------------------------------
C     SUBROUTINE       DRIFT
C     AUTHOR:          M.D. Nozar and P.E. Ulmer      
C     DATE:            2-AUG-1991
C     PURPOSE:         
C                      Obtain the transport vector of the particle
C                      after a drift through a field-free region.
C
C                      Note:
C                             1-  The transformation is exact.
C                             2-  Only the X,Y, and L components of   
C                                 the transport vector are affected  
C                                 upon passage through a drift element.
C ----------------------------------------------------------------------
C
       SUBROUTINE DRIFT(VECT,VECT_DR,DRIFT_DIST)
C
       IMPLICIT          NONE
C                                
       DOUBLE PRECISION  VECT(6),VECT_DR(6),DRIFT_DIST,TTHET,TPHIT
C
       TTHET = TAN(1.D-3 * VECT(2))   !tan(theta_T)
       TPHIT = TAN(1.D-3 * VECT(4))   !tan(phi_T)
C
       VECT_DR(1) = VECT(1) + DRIFT_DIST*TTHET
       VECT_DR(3) = VECT(3) + DRIFT_DIST*TPHIT
       VECT_DR(5) = VECT(5) + DRIFT_DIST*
     #                          (SQRT(1.D0+TTHET**2+TPHIT**2)-1.D0)
C
       VECT_DR(2) = VECT(2)
       VECT_DR(4) = VECT(4)
       VECT_DR(6) = VECT(6)
C
       RETURN
       END
