C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE PSEUDO_BEAM
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    08-JUL-1999
C
C       PURPOSE: The pseudo-beam is the difference between the
C                beam 4-momentum and the photon 4-momentum.
C                For pre-radiation, this is just the effective beam
C                4-momentum before the scattering.  For post-radiation,
C                this is an artificial quantity useful for determining
C                the modified kinematics.  In both cases, only in
C                peaking approximation does the pseudo-beam 4-vector
C                represent an on-shell electron.
C
C       INPUT:   Photon and asymptotic beam kinematics.
C       OUTPUT:  "Pseudo-beam" 4-vector:  eb_ps and ki_ps(3)
C
C -----------------------------------------------------------------------
C
      SUBROUTINE PSEUDO_BEAM(eb,ki,eg,kg,eb_ps,ki_ps)
      IMPLICIT NONE
C
      DOUBLE PRECISION eb,ki(3),eg,kg(3)
      DOUBLE PRECISION eb_ps,ki_ps(3)
      INTEGER          i
C
C -----------------------------------------------------------------------
C     Make a "pseudo-beam" 4-momentum by removing photon momentum/energy.
C     For pre-radiation, this actually describes the beam at the virtual
C     photon vertex.  For post-radiation, this quantity is artificial
C     but can still be used to determine the vertex kinematics.
C -----------------------------------------------------------------------
C
      do i=1,3
         ki_ps(i) = ki(i) - kg(i)
      enddo
      eb_ps = eb - eg       ! energy
C
      return
      end
