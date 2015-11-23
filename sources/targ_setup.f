C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE TARG_SETUP
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    4-AUG-1999
C
C       PURPOSE: Set up some parameters for target to include
C                in common block:  /tg1/
C
C       VARIABLES:
C                   zin:       location of cell entrance window (m)
C                   zout:      location of cell exit     window (m)
C                   atg:       mass   number of target
C                   ztg:       atomic number of target
C                   radius:    radius of target cell (m)
C                   imep:      mean excitation energy of target (eV)
C                   rho:       density of target (g/cm^3)
C
C                   imodel:    1=foil(s), 2=LH2/LD2 cryo, 3=He cryo
C                              4=polarized He3 target
C -----------------------------------------------------------------------
C
C     Modified to include the cigar tube target (Model 5).
C     Hassan Ibrahim, December 7, 2004
C

      SUBROUTINE TARG_SETUP(ztarg,atarg,density,z1,z2,imodel)
      IMPLICIT NONE
C
      DOUBLE PRECISION ztarg,atarg,density,z1,z2
      INTEGER          imodel
C
      INCLUDE 'eloss.cmn'
C
      ztg  = nint(ztarg)
      atg  = atarg
      rho  = density
      zin  = z1
      zout = z2
C
      if(imodel .eq. 1) then
         radius = 0.d0          ! set to zero for foil(s)
      elseif(imodel .eq. 2) then
         radius = 0.0318d0      ! for cryotarget beer can (meters)
      elseif(imodel .eq. 3) then
         radius = 0.0516d0      ! for He tuna can (meters)
      elseif(imodel .eq. 4) then 
         radius = 0.00955d0        ! for Polarized He3 target (meters)
      elseif(imodel .eq. 5) then
         radius = 0.02033d0        ! for cryotarget cigar tube (meters)
      endif
C
      return
      end
