C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE PEAKING
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    09-JUL-1999
C
C       PURPOSE: Sample the photon 3-momentum for internal Bremsstrahlung.
C                The peaking approximation is used to sample the photon
C                angles which amounts to alternately choosing the
C                incident and scattered beam directions for the direction
C                of the radiated photon.
C
C                Notes: For the elastic and bound state cases,
C                       acceptance tests must be applied to the two
C                       kinematics for pre- and post-radiation.  If
C                       both kinematics pass, then one will be selected
C                       by a random coin toss and a weight of two
C                       will be assigned to the one selected.
C
C                       The maximum photon energy is chosen to be the
C                       minimum asymptotic scattered electron momentum.
C                       This ignores the part of the tail at very
C                       high energy.
C
C       INPUT:   Incident beam (before radiation) angles and
C                scattered electron (after radiation) angles.
C
C       OUTPUT:  Photon energy & 3-momentum & unit-vector:
C                     eg_pk(2), kg_pk(2,3), kg_pk_unit(2,3)
C
C                  eg_pk(1) and eg_pk(2) correspond to the photon energies
C                  assuming pre- and post-radiation, respectively.
C
C                  kg_pk(1,i) and kg_pk(2,i) are the photon 3-momenta
C                  defined similarly.
C
C                  kg_pk_unit(1,i) and kg_pk_unit(2,i) are unit vectors
C                  along the photon momentum for the two cases.
C
C -----------------------------------------------------------------------
C
      SUBROUTINE PEAKING(phb,thb,phe,the,eg_pk,kg_pk,kg_pk_unit)
      IMPLICIT NONE
C
      COMMON /RADPKAPPR/ cutoff,eg_max
C
      DOUBLE PRECISION phb,thb,phe,the
      DOUBLE PRECISION eg_pk(2),kg_pk(2,3),kg_pk_unit(2,3)
      DOUBLE PRECISION phg(2),thg(2),cutoff,eg_max,kg_tmp(3)
      REAL             ran_unif
      INTEGER          i,j
C
C -----------------------------------------------------------------------
C     Determine maximum photon energy and photon angles.
C -----------------------------------------------------------------------
C
      phg(1) = phb      ! pre-radiation
      thg(1) = thb
C
      phg(2) = phe      ! post-radiation
      thg(2) = the
C
C -----------------------------------------------------------------------
C     Determine photon energy based on 1/E weighting.
C     Build photon 3-momentum and unit vector along photon direction.
C
C     (This unit vector is useful for calculating other quantities,
C     such as the Jacobian of EF wrt photon energy for elastic scattering,
C     in a way which does not give large errors as the
C     photon energy --> 0.)
C -----------------------------------------------------------------------
C
      do i=1,2
        call RANECU(ran_unif,1)
        eg_pk(i) = cutoff*(eg_max/cutoff)**DBLE(ran_unif)
        call V3MAKE(1.d0,phg(i),thg(i),kg_tmp)  ! unit-vector for photon
        do j=1,3
           kg_pk_unit(i,j) = kg_tmp(j)
           kg_pk(i,j)      = eg_pk(i)*kg_pk_unit(i,j)
        enddo
      enddo
C
      return
      end

