C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE BORIE_DRECHSEL
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    12-JUL-1999
C
C       PURPOSE: Uses the formalism of Borie and Drechsel for the
C                radiative cross section in (e,e'p) in peaking
C                approximation.
C
C                E. Borie and D. Drechsel, Nucl. Phys. A167, 369 (1971).
C
C       INPUT:
C               irad:  1 for pre-radiation and 2 for post-radiation
C               eg:    photon energy
C               e0:    asymptotic beam energy
C               ef:    asymptotic scattered electron energy
C
C       OUTPUT:
C               radfac:  factor which multiplies the unradiated cross
C                          section.
C
C -----------------------------------------------------------------------
C
      SUBROUTINE BORIE_DRECHSEL(irad,eg,e0,ef,radfac)
      IMPLICIT NONE
C
      DOUBLE PRECISION eg,e0,ef,radfac
      DOUBLE PRECISION alpha,pi,me,fact
      INTEGER          irad
C
      PARAMETER (alpha  = 7.29735d-3) 
      PARAMETER (pi  = 3.141592653590d0) 
      PARAMETER (me  = 0.511d0)           ! electron mass (MeV)
C
      fact = alpha/(pi*eg)
      if(irad .eq. 1) then
         radfac = fact * ((e0**2 + (e0-eg)**2)/(e0**2))
     #            * log(2.d0*e0/me)
      else
         radfac = fact * (((ef+eg)**2 + ef**2)/((ef+eg)**2))
     #            * log(2.d0*ef/me)
      endif
C
      return
      end
