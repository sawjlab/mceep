      SUBROUTINE schwinger(IKIND,ATARG,CUTOFF,COR_TOT,COR_MULTI)
C------------------------------------------------------------------------------
C Radiative corrections:
C           only SCHWINGER CORRECTION
C------------------------------------------------------------------------------
C  The Spence function routines were taken from aeexb.
C This routine calculates the Schwinger radiative correction 
C assuming that the experimental spectrum has a
C specified resolution (real photon energy cutoff).
C  Also, S. Penner, Proc. of the 18th
C Scottish Univ. Summer School in Phys. (1977)284 for the kinematical recoil
C correction.  The routine returns only the virtual photon correction if the
C real photon energy cutoff is not positive.
C------------------------------------------------------------------------------
      IMPLICIT NONE
c
      INCLUDE 'var.cmn'
c 
c     from var.cmn need
c          OMEGA  = energy transfer (MeV)
c          QMU2_G = four-momentum transfer (GeV/c)^2
c          TSCAT  = electron scattering angle (radians)
c          E0_I   = incident electron energy (MeV)
c          PF_E_I = momentum of scattered electron (MeV)
c     input
c          ATARG = atomic weight of target
c          CUTOFF = cutoff energy 
c
c     output
c          cor_tot   = schwinger correction
c          cor_multi = multi-photon correction
c
      DOUBLE PRECISION cor_real,cor_virtual,cor_tot,cor_multi
      DOUBLE PRECISION atarg
      DOUBLE PRECISION cutoff
      DOUBLE PRECISION Q2       ! four momentum in (MeV/c)^2
      DOUBLE PRECISION E       ! scattered electron energy
      DOUBLE PRECISION LN,SPEN
      DOUBLE PRECISION DM      ! target mass in MeV/c2
      DOUBLE PRECISION PSI,KAPPA
      DOUBLE PRECISION DELTA_VIRTUAL,DELTA_REAL
      DOUBLE PRECISION ARG
      INTEGER          ikind

      
c          
C Parameters:
      DOUBLE PRECISION      M_E,M_E2,ALPHA,PI,A1,A2,A3,A4,AMU
      PARAMETER ( M_E   = 5.11D-1                      )
      PARAMETER ( M_E2  = M_E*M_E                      )
      PARAMETER ( ALPHA = 1.D0/.137036D3               )
      PARAMETER ( PI    = 3.141592653589793238462643D0 )
      PARAMETER ( A1    = ALPHA/PI                     )
      PARAMETER ( A2    = 2.D0*A1                      )
      PARAMETER ( A3    = 2.8D1/9.D0 + PI*PI/6.D0      )
      PARAMETER ( A4    = 1.3D1/6.D0                   )
      PARAMETER ( AMU    = 931.502D0                   )
C
C Calculate the virtual photon radiative correction.
      E             = DSQRT (PF_E_I**2.D0  + M_E2)
      Q2            = QMU2_G/1.0d-6 ! in units of (MeV/c)^2
      LN            = DLOG( Q2/M_E2 )
      ARG           = DCOS( TSCAT/2.D0 )
      ARG           = ARG * ARG
      IF ( IKIND.GE.0 ) THEN                              !inelastic scattering
         DELTA_VIRTUAL = A3 - A4 * LN + 5.D-1 * (DLOG( E0_I/E ))**2.D0
     &                 - SPEN(ARG)
      ELSE                                                !  elastic scattering
         DELTA_VIRTUAL = A3 - A4 * LN - SPEN(ARG)
      ENDIF
      DELTA_VIRTUAL = A1 * DELTA_VIRTUAL
      COR_VIRTUAL   = 1.D0 - DELTA_VIRTUAL
      IF ( CUTOFF .LE. 0. ) THEN      ! no real photon radiative correction
         COR_TOT = COR_VIRTUAL
         RETURN
      ENDIF
C Calculate the real photon radiative correction.
      DM = atarg*amu  ! target mass in MeV
      PSI           = 1.D0 + E0_I/DM * 2.D0*(1.D0 - ARG)
      IF ( IKIND.GE.0 ) THEN  ! inelastic
         KAPPA      = 1.D0 + DBLE(OMEGA)/DM * 2.D0*(1.D0 - ARG)
         DELTA_REAL = A2 * DLOG( DSQRT(KAPPA*E0_I*E)/PSI/CUTOFF )
     &                   * ( LN - 1.D0 )
      ELSE ! elastic
         DELTA_REAL = A2 * DLOG( E0_I/PSI**1.5D0/CUTOFF )
     &                   * ( LN - 1.D0 )
      ENDIF
C Exponentiate the real photon radiative correction and return the total
C correction.
        COR_REAL  = EXP(-DELTA_REAL)
        COR_TOT   = COR_REAL * COR_VIRTUAL
C
C Calculate multi-photon correction (ala J. Templon)  - PEU
C    I include the virtual photon correction to
C    the tail here as well (PEU).
C
        COR_MULTI = COR_VIRTUAL * ((1.D0-COR_REAL)/DELTA_REAL)
C
      RETURN
      END
c
C-----------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION SPEN( Y )
C-----------------------------------------------------------------------------
CEV This routine computes the Spence function for 0 .lt. y .le. 1,
C   by an interpolation on tabulated values.
C
      IMPLICIT NONE
C Passed variable:
      DOUBLE PRECISION    Y
C Local variables:
      DOUBLE PRECISION    Y1,P
      INTEGER N
      DOUBLE PRECISION    ORD(12)
      SAVE ORD
      DATA ORD/-9.76D-2, 0.D0,    1.026D-1,2.110D-1,3.261D-1,4.493D-1,
     &          5.882D-1,7.276D-1,8.894D-1,1.0748D0,1.2997D0,1.6449D0/
C
      IF (Y .EQ. 1.D0) THEN
         SPEN = ORD(12)
         RETURN
      ENDIF
      Y1   = 1.D1*Y + 2.D0
      N    = INT(Y1)
      P    = Y1 - DFLOAT(N)
      SPEN =  ORD(N-1)*P*(P  - 1.D0)/2.D0 + 
     &        ORD(N)  *(1.D0 - P**2)      +
     &        ORD(N+1)*P*(P  + 1.D0)/2.D0
C
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SPENCE( X )
CEV This routine calculates the Spence function for all real values of the
C   argument X.  The method applied is an advanced version of the method of
C   N.T. Meister and D.R. Yennie, Phys. Rev. 130, 1210(1963).
C
      IMPLICIT NONE
C Passed variable:
      DOUBLE PRECISION X
C Local variables:
      DOUBLE PRECISION PISQ(3),Y,Y1,SPCE
C The following numbers are (1/3)pi**2, (1/6)pi**2, and -(1/12)pi**2,
C in this order:
      DATA   PISQ/3.2898681336964D0,1.6449340668482D0,-.8224670334241D0/
      SAVE   PISQ
C
      IF ( X.GT.2.D0 ) THEN
         Y  = 1.D0 / X
         SPENCE = PISQ(1) - SPCE(Y) - 5.D-1*DLOG(X)**2.D0
      ELSE IF ( X.GT.1.D0 ) THEN
         Y  = (X - 1.D0) / X
         Y1 = (X - 1.D0) * Y
         SPENCE = PISQ(2) + SPCE(Y) - 5.D-1*DLOG(X)*DLOG(Y1)
      ELSE IF ( X.EQ.1.D0 ) THEN
         SPENCE = PISQ(2)
      ELSE IF ( X.GT.5.D-1 ) THEN
         Y  = 1.D0 - X
         SPENCE = PISQ(2) - SPCE(Y) - DLOG(X)*DLOG(Y)
      ELSE IF ( X.LT.-1.D0 ) THEN
         Y  = 1.D0 - X
         Y1 = 1.D0 / Y
         SPENCE = SPCE(Y1) - PISQ(2) - 5.D-1*DLOG(Y)*DLOG(X*X/Y)
      ELSE IF ( X.EQ.-1.D0 ) THEN
         SPENCE = PISQ(3)
      ELSE IF ( X.LT.0.D0 ) THEN
         Y  = X / (X - 1.D0)
         Y1 = 1.D0 - X
         SPENCE = - SPCE(Y) - 5.D-1*DLOG(Y1)**2.D0
      ELSE
         SPENCE = SPCE(X)
      ENDIF
C
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SPCE( X )
CEV Series expansion of the Spence function for the argument X in the range
C   from -1 to 1, inclusive.  See N.T. Meister and D.R. Yennie, Phys. Rev. 130,
C   1210(1963).
C
      IMPLICIT NONE
C Passed variable:
      DOUBLE PRECISION X
C Local variables:
      INTEGER I
      DOUBLE PRECISION TOP
C
      SPCE = 0.D0
      TOP  = 1.D0
      DO I=1,20
         TOP  =  TOP * DABS(X)
         SPCE = SPCE + TOP/DFLOAT(I)**2.D0
      ENDDO
C
      RETURN
      END
