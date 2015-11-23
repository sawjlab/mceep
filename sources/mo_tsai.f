C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE MO_TSAI
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    14-AUG-1999
C
C       PURPOSE: Uses the formalism of Mo and Tsai for the
C                radiative tail in (e,e') elastic scattering
C                in peaking approximation.
C
C                L.W. Mo and Y.S. Tsai, Rev. Mod. Phys. 41, 205 (1969).
C
C       INPUT:
C               irad:  1 for pre-radiation and 2 for post-radiation
C               eg:    photon energy
C               e0:    asymptotic beam energy
C               ef:    asymptotic scattered electron energy
C               theta: electron scattering angle
C
C       OUTPUT:
C               radfac:  factor which multiplies the unradiated cross
C                          section.
C
C -----------------------------------------------------------------------
C
      SUBROUTINE MO_TSAI(irad,eg,e0,ef,theta,radfac)
      IMPLICIT NONE
c
      DOUBLE PRECISION eg,e0,ef,theta,radfac
      DOUBLE PRECISION alpha,pi,me,fact0,fact_ang,fact,sp,x,t
      INTEGER          irad
c
      INCLUDE 'masses.cmn'
c
      PARAMETER (alpha  = 7.29735d-3) 
      PARAMETER (pi  = 3.141592653590d0) 
      PARAMETER (me  = 0.511d0)           ! electron mass (MeV)
c
      fact0    = alpha/pi
      fact_ang = 1.d0-cos(theta)
      sp       = e0*ef*fact_ang
c
      if(irad .eq. 1) then
         x    = (e0-eg)/e0
         fact =    (eject_mass + (e0-eg)*fact_ang)
     #            /(eject_mass - ef*fact_ang)
      else
         x    = ef/(ef+eg)
         fact = 1.d0
      endif
c
      t = fact0*(0.5d0*(1.0d0+x*x)*log(2.d0*sp/(me*me)) - x)
c
      radfac = t*fact/eg
c
      return
      end
