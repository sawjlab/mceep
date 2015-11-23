C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE MAD12NQ
C
C       AUTHOR:  Paul Ulmer
C                (based on mad35.f by Kathy McCormick)
C       DATE:    14-DEC-2004
C       PURPOSE: Produces focal plane Transport vector for the
C                JLAB Hall A MAD spectrometer (12 deg configuration)
C                with QUADS turned OFF given the target Transport
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
C -----------------------------------------------------------------------
C
      SUBROUTINE MAD12NQ(x_tg,x_fp,fail_apertures)
      IMPLICIT NONE
C
      DOUBLE PRECISION x_tg(6),x_fp(6)
      DOUBLE PRECISION mag1_radius,mag2_radius
      DOUBLE PRECISION x1,t1,y1,p1,l1
      DOUBLE PRECISION xoffset(15)  ! central ray location wrt magnet axis
      REAL x(5)
      REAL x_mnq_0_2, t_mnq_0_2, y_mnq_0_2, p_mnq_0_2, pl_mnq_0_2
                                          !entrance of mag1
      REAL x_mnq_0_3, t_mnq_0_3, y_mnq_0_3, p_mnq_0_3, pl_mnq_0_3
                                          !middle   of mag1
      REAL x_mnq_0_4, t_mnq_0_4, y_mnq_0_4, p_mnq_0_4, pl_mnq_0_4
                                          !exit     of mag1
      REAL x_mnq_0_5, t_mnq_0_5, y_mnq_0_5, p_mnq_0_5, pl_mnq_0_5
                                          !entrance of mag2
      REAL x_mnq_0_7, t_mnq_0_7, y_mnq_0_7, p_mnq_0_7, pl_mnq_0_7
                                          !middle   of mag2
      REAL x_mnq_0_9, t_mnq_0_9, y_mnq_0_9, p_mnq_0_9, pl_mnq_0_9
                                          !exit     of mag2
      REAL x_mnq_0_11,t_mnq_0_11,y_mnq_0_11,p_mnq_0_11,pl_mnq_0_11
                                          !x,th,y,phi,l at 1st chamber
      REAL x_mnq_0_15,t_mnq_0_15,y_mnq_0_15,p_mnq_0_15,pl_mnq_0_15
                                          !x,th,y,phi,l at calorimeter
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
C     set of planes (so some offsets are set to zero, since LeRose 
C     doesn't provide functions for every plane).
C -----------------------------------------------------------------------
C
      DATA xoffset /0.d-3,-0.671598d-3,112.408371d-3,-21.088808d-3,
     #       71.681625d-3,0.d-3,243.803665d-3,0.d-3,48.602173d-3,
     #       0.d-3,-88.020126d-3,0.d-3,0.d-3,0.d-3,-50.74712d-3/ ! in m
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
      x1 = DBLE(x_mnq_0_2(x,5)) + xoffset(2)           ! mag1 entrance
      y1 = DBLE(y_mnq_0_2(x,5))
      if(x1*x1+y1*y1 .gt. mag1_radius*mag1_radius
     &               .and. aper_mad_test(1)) then
         fail_apertures = .true.
         goto 234
      endif
      mag1_ent_x = x1      ! X histogram
      mag1_ent_y = y1      ! Y histogram
C -------------------------------------------------------------
      x1 = DBLE(x_mnq_0_3(x,5)) + xoffset(3)           ! mag1 middle
      y1 = DBLE(y_mnq_0_3(x,5))
      if(x1*x1+y1*y1 .gt. mag1_radius*mag1_radius
     &               .and. aper_mad_test(2)) then
         fail_apertures = .true.
         goto 234
      endif
      mag1_mid_x = x1      ! X histogram
      mag1_mid_y = y1      ! Y histogram
C -------------------------------------------------------------
      x1 = DBLE(x_mnq_0_4(x,5)) + xoffset(4)           ! mag1 exit
      y1 = DBLE(y_mnq_0_4(x,5))
      if(x1*x1+y1*y1 .gt. mag1_radius*mag1_radius
     &               .and. aper_mad_test(3)) then
         fail_apertures = .true.
         goto 234
      endif
      mag1_ext_x = x1      ! X histogram
      mag1_ext_y = y1      ! Y histogram
C -------------------------------------------------------------
      x1 = DBLE(x_mnq_0_5(x,5)) + xoffset(5)           ! mag2 entrance
      y1 = DBLE(y_mnq_0_5(x,5))
      if(x1*x1+y1*y1 .gt. mag2_radius*mag2_radius
     &               .and. aper_mad_test(4)) then
         fail_apertures = .true.
         goto 234
      endif
      mag2_ent_x = x1      ! X histogram
      mag2_ent_y = y1      ! Y histogram
C -------------------------------------------------------------
      x1 = DBLE(x_mnq_0_7(x,5)) + xoffset(7)           ! mag2 middle
      y1 = DBLE(y_mnq_0_7(x,5))
      if(x1*x1+y1*y1 .gt. mag2_radius*mag2_radius
     &               .and. aper_mad_test(5)) then
         fail_apertures = .true.
         goto 234
      endif
      mag2_mid_x = x1      ! X histogram
      mag2_mid_y = y1      ! Y histogram
C -------------------------------------------------------------
      x1 = DBLE(x_mnq_0_9(x,5)) + xoffset(9)           ! mag2 exit
      y1 = DBLE(y_mnq_0_9(x,5))
      if(x1*x1+y1*y1 .gt. mag2_radius*mag2_radius
     &               .and. aper_mad_test(6)) then
         fail_apertures = .true.
         goto 234
      endif
      mag2_ext_x = x1      ! X histogram
      mag2_ext_y = y1      ! Y histogram
C
C -----------------------------------------------------------------------
C     Get focal plane Transport vector.
C -----------------------------------------------------------------------
C
      x1 = DBLE(x_mnq_0_11(x,5))             ! mad focal plane
      t1 = DBLE(t_mnq_0_11(x,5))
      y1 = DBLE(y_mnq_0_11(x,5))
      p1 = DBLE(p_mnq_0_11(x,5))
      l1 = DBLE(pl_mnq_0_11(x,5))
C
C -----------------------------------------------------------------------
C     Repack array and change units back to MCEEP conventions.
C -----------------------------------------------------------------------
C
      x_fp(1) = x1*100.d0          ! m --> cm
      x_fp(2) = atan(t1)*1000.d0   ! mr
      x_fp(3) = y1*100.d0          ! m --> cm
      x_fp(4) = atan(p1)*1000.d0   ! mr
      x_fp(5) = l1*100.d0          ! m --> cm
      x_fp(6) = x_tg(6)            ! delta not modified
C
 234  RETURN
      END
      







