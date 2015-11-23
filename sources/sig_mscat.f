C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       FUNCTION SIG_MSCAT
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    12-SEP-1999
C
C       PURPOSE: Calculates sigma of multiple scattering distribution
C                based on gaussian approximation (needed for 'TRK'
C                option).
C       INPUT:
C               ia:      IARM (explicitly assumes electron for IARM=1
C                              and ejectile for IARM=2)
C               emom:    particle momentum in MeV/c (electron)
C               pmom:    particle momentum in MeV/c (ejectile)
C               x_x0:    x/x0 for wire chamber
C                             x0 = radiation length of "chamber material"
C                             x  = actual thickness of wire chamber
C               cangle:  cos of wire chamber angle (to get effective
C                        thickness of wire chamber seen by particle).
C
C       OUTPUT:
C               sig_mscat:  sigma of gaussian mult. scatt. distribution
C                           in mrad
C
C -----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION SIG_MSCAT(ia,emom,pmom,x_x0,cangle)
      IMPLICIT NONE
C
      DOUBLE PRECISION emom,pmom,x_x0,cangle
      DOUBLE PRECISION x_x0_part,p,z,beta
      INTEGER          ia
C
      INCLUDE 'masses.cmn'
C
      x_x0_part = x_x0/cangle   ! x/x0 seen by particle
C
      if(ia .eq. 1) then        ! electron
         p    = emom
         z    = 1.d0
         beta = 1.d0
      else                      ! ejectile
         p    = pmom
         z    = eject_chg
         beta = p/sqrt(p**2+eject_mass**2)
      endif
C
      sig_mscat = 1000.d0*(13.6d0/(beta*p))*z*sqrt(x_x0_part)
     #            * (1.d0 + 0.038d0*log(x_x0_part))  ! in mrad
C
      return
      end
