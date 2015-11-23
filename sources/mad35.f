C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE MAD35
C
C       AUTHOR:  Kathy McCormick
C                (based on hrs.f by Paul E. Ulmer)
C       DATE:    09-APRIL-2002
C       PURPOSE: Produces focal plane Transport vector for the
C                JLAB Hall A MAD spectrometer (35 deg configuration)
C                given the target Transport
C                vector.  Also gives the value of the logical flag
C                fail_apertures consistent with the size of various
C                apertures.
C
C       USES functions of J. LeRose to map spectrometer coordinates.
C
C       INPUT:   x_tg:            Target Transport vector
C       OUTPUT:  x_fp:            Focal plane Transport vector
C                fail_apertures:  True if any aperture is not passed.
C
C       Modifications:
C
C       15-DEC-2004   P.E. Ulmer
C       Incorporate X offsets for aperture tests.  This is necessary, 
C       since the central ray is NOT centered on the magnet axis.
C
C -----------------------------------------------------------------------
C
      SUBROUTINE MAD35(x_tg,x_fp,fail_apertures)
      IMPLICIT NONE
C
      DOUBLE PRECISION x_tg(6),x_fp(6)
      DOUBLE PRECISION mag1_radius,mag2_radius
      DOUBLE PRECISION x1,t1,y1,p1
      DOUBLE PRECISION xoffset(13)  ! central ray location wrt magnet axis
      REAL x(5)
      REAL x_m35_0_2,y_m35_0_2,x_m35_0_4,y_m35_0_4 
                                          !entrance and exit of mag1
      REAL x_m35_0_5,y_m35_0_5,x_m35_0_11,y_m35_0_11 
                                          !entrance and exit of mag2
      REAL x_m35_0_13,t_m35_0_13,y_m35_0_13,p_m35_0_13 
                                          !x,th,y,phi at 1st chamber
      INTEGER i
      LOGICAL fail_apertures
C
      INCLUDE 'spectrometer.cmn'
      INCLUDE 'mad.cmn'
C
C -----------------------------------------------------------------------
C     Parameters describing location in the local (magnet) coordinate 
C     system of the central trajectory (+x is down).
C
C     The index into this array refers to the plane # in LeRose's
C     set of planes.
C -----------------------------------------------------------------------
C
      DATA xoffset /0.d-3,-70.027237d-3,42.759571d-3,-120.611023d-3,
     #              -120.390579d-3,-26.365952d-3,63.262596d-3,
     #              84.443748d-3,92.592392d-3,92.738144d-3,
     #              89.496307d-3,85.170967d-3,401.250183d-3/    ! in m
C
C -----------------------------------------------------------------------
C     Parameters describing apertures.
C        1st magnet:  circle (entrance and exit).
C        2nd magnet:  circle (entrance and exit).
C -----------------------------------------------------------------------
C
      DATA mag1_radius /0.60d0/   ! circle (meters)
      DATA mag2_radius /0.60d0/   ! circle (meters)
C
C -----------------------------------------------------------------------
C     Initialize the focal plane vector.  This is necessary
C     since if the aperture tests fail, some elements of the focal
C     plane vector are arbitrary.  These arbitrary values could lead to
C     huge numbers upon applying the inverse transforms which could,
C     in turn, give rise to floating overflows.
C -----------------------------------------------------------------------
C
      do i=1,6
         x_fp(i) = 0.d0
      enddo
C
C -----------------------------------------------------------------------
C     Repack the array and change units to match LeRose's conventions.
C -----------------------------------------------------------------------
C
      x(1) = REAL(x_tg(1)/100.d0)        ! cm --> meters
      x(2) = REAL(tan(x_tg(2)/1000.d0))  ! tan(theta)
      x(3) = REAL(x_tg(3)/100.d0)        ! cm --> meters
      x(4) = REAL(tan(x_tg(4)/1000.d0))  ! tan(phi)
      x(5) = REAL(x_tg(6)/100.d0)        ! % --> fraction
C
C -----------------------------------------------------------------------
C     Check all apertures.
C -----------------------------------------------------------------------
C
         x1 = DBLE(x_m35_0_2(x,5)) + xoffset(2)           ! mag1 entrance
         y1 = DBLE(y_m35_0_2(x,5))
      if(x1*x1+y1*y1 .gt. mag1_radius*mag1_radius
     &               .and. aper_mad_test(1)) then
         fail_apertures = .true.
         goto 123
      endif
C
C histograms of X and Y
         mag1_ent_x = x1
         mag1_ent_y = y1
C
         x1 = DBLE(x_m35_0_4(x,5)) + xoffset(4)           ! mag1 exit
         y1 = DBLE(y_m35_0_4(x,5))
      if(x1*x1+y1*y1 .gt. mag1_radius*mag1_radius
     &               .and. aper_mad_test(2)) then
         fail_apertures = .true.
         goto 123
      endif
C
C histograms of X and Y
         mag1_ext_x = x1
         mag1_ext_y = y1
C
C -----------------------------------------------------------------------
         x1 = DBLE(x_m35_0_5(x,5)) + xoffset(5)           ! mag2 entrance
         y1 = DBLE(y_m35_0_5(x,5))
      if(x1*x1+y1*y1 .gt. mag2_radius*mag2_radius
     &               .and. aper_mad_test(3)) then
         fail_apertures = .true.
         goto 123
      endif
C
C histograms of X and Y
         mag2_ent_x = x1
         mag2_ent_y = y1
C
         x1 = DBLE(x_m35_0_11(x,5)) + xoffset(11)         ! mag2 exit
         y1 = DBLE(y_m35_0_11(x,5))
      if(x1*x1+y1*y1 .gt. mag2_radius*mag2_radius
     &               .and. aper_mad_test(4)) then
         fail_apertures = .true.
         goto 123
      endif
C
C histograms of X and Y
         mag2_ext_x = x1
         mag2_ext_y = y1
C
C -----------------------------------------------------------------------
C     Get focal plane Transport vector.
C -----------------------------------------------------------------------
C
         x1 = DBLE(x_m35_0_13(x,5))             ! mad focal plane
         t1 = DBLE(t_m35_0_13(x,5))
         y1 = DBLE(y_m35_0_13(x,5))
         p1 = DBLE(p_m35_0_13(x,5))
C
C -----------------------------------------------------------------------
C     Repack array and change units back to MCEEP conventions.
C -----------------------------------------------------------------------
C
      x_fp(1) = x1*100.d0          ! m --> cm
      x_fp(2) = atan(t1)*1000.d0   ! mr
      x_fp(3) = y1*100.d0          ! m --> cm
      x_fp(4) = atan(p1)*1000.d0   ! mr
      x_fp(6) = x_tg(6)            ! delta not modified
C
C -----------------------------------------------------------------------
C     Calculate path length (in cm) here.
C     Overall path is L=L0+x_fp(5) where L0 is the central ray
C     path length.
C
C     For now, set to value at target.
C -----------------------------------------------------------------------
C
      x_fp(5) = 0.d0 + x_tg(5)
C
 123  RETURN
      END
      







