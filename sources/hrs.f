C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE HRS
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    25-MAR-1999
C       PURPOSE: Produces focal plane Transport vector for the
C                JLAB Hall A spectrometers given the target Transport
C                vector.  Also gives the value of the logical flag
C                fail_apertures consistent with the size of various
C                apertures.
C
C       MODS:    14-DEC-1999 PEU
C                If the aperture cuts fail, the subroutine
C                still calculates the focal plane vector, so that
C                the before cuts counters are physically meaningful.
C
C                28-JAN-2000 PEU
C                If the aperture cuts fail, the subroutine
C                now sets the focal plane vector to zero (as it
C                used to).  Otherwise floating errors can result
C                on some platforms since the eloss/mscatt routines
C                occasionally give large effects which the LeRose
C                polynomials are not meant to handle.  Note that
C                events outside the apertures fail the cuts anyway,
C                so all histograms and Ntuples are unaffected.
C
C                15-MAY-2001 PEU
C                Reduce radius of Q3 from 0.30 m to 0.28 m (v3.6).
C                This seems to give better agreement with white
C                spectrum delta scan connected with E94-004.
C                Also, added histograms of X and Y at intermediate
C                apertures
C
C                24-JAN-2003 PEU
C                Restore radius of Q3 to 0.30 m from 0.28 m (v3.8)
C                since new LeRose transfer functions model the
C                acceptance better, without need for this "fudge".
C
C       USES functions of J. LeRose to map spectrometer coordinates.
C
C       INPUT:   x_tg:            Target Transport vector
C       OUTPUT:  x_fp:            Focal plane Transport vector
C                fail_apertures:  True if any aperture is not passed.
C
C -----------------------------------------------------------------------
C
      SUBROUTINE HRS(ispect,x_tg,x_fp,fail_apertures)
      IMPLICIT NONE
C
      DOUBLE PRECISION x_tg(6),x_fp(6)
      DOUBLE PRECISION q1_radius(2),q3_radius(2)
      DOUBLE PRECISION dip1(2),dip2(2),dip3(2),dip4(2)
      DOUBLE PRECISION x1,t1,y1,p1,ylim
      REAL X(5)
      REAL X_E_Q1EX,Y_E_Q1EX,X_E_DENT,Y_E_DENT,X_E_DEXT,Y_E_DEXT
      REAL X_E_Q3EN,Y_E_Q3EN,X_E_Q3EX,Y_E_Q3EX
      REAL X_E_FP,T_E_FP,Y_E_FP,P_E_FP
      REAL X_H_Q1EX,Y_H_Q1EX,X_H_DENT,Y_H_DENT,X_H_DEXT,Y_H_DEXT
      REAL X_H_Q3EN,Y_H_Q3EN,X_H_Q3EX,Y_H_Q3EX
      REAL X_H_FP,T_H_FP,Y_H_FP,P_H_FP
      INTEGER ispect,i
      LOGICAL fail_apertures
C
      INCLUDE 'spectrometer.cmn'
      INCLUDE 'hrs.cmn'
C
C -----------------------------------------------------------------------
C     Parameters describing apertures.
C        1st quad:  circle (exit).
C        3rd quad:  circle (entrance and exit).
C        Dipole:    trapezoid (entrance and exit).
C
C     Each parameter is a two valued array where the first value
C     gives the parameter for the electron arm and the second for
C     the hadron arm.  The parameters are taken to be identical for
C     each arm however.
C
C     Index defining which spectrometer is being used:
C         ispect = 1  (electron arm)
C         ispect = 2  (hadron arm)
C     Note that this allows either arm to be used for detection of
C     either particle (i.e. "electron" arm could be used to detect
C     hadrons) since this index can be set to 1 or 2 regardless of
C     the value of IARM.
C -----------------------------------------------------------------------
C
      DATA q1_radius /0.1492d0,0.1492d0/ ! circle (meters)
Cv3.6      DATA q3_radius /0.280d0,0.280d0/   ! circle (meters)
      DATA q3_radius /0.300d0,0.300d0/   ! circle (meters)
C
      DATA dip1 /0.40d0,0.40d0/          ! trapezoid
      DATA dip2 /0.125d0,0.125d0/
      DATA dip3 /125.d0,125.d0/
      DATA dip4 /840.d0,840.d0/
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
      x(4) = REAL(tan(x_tg(4)/1000.d0))  ! tan(theta)
      x(5) = REAL(x_tg(6)/100.d0)        ! % --> fraction
C
C -----------------------------------------------------------------------
C     Check all apertures.
C -----------------------------------------------------------------------
C
C
      if(ispect .eq. 1) then
         x1 = DBLE(X_E_Q1EX(x,5))           ! E_Q1 exit
         y1 = DBLE(Y_E_Q1EX(x,5))
      else
         x1 = DBLE(X_H_Q1EX(x,5))           ! H_Q1 exit
         y1 = DBLE(Y_H_Q1EX(x,5))
      endif
      if(x1*x1+y1*y1 .gt. q1_radius(ispect)*q1_radius(ispect)
     &               .and. aperture_test(ispect,1)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y
         eq1_ext_x = x1
         eq1_ext_y = y1
      else
         hq1_ext_x = x1
         hq1_ext_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(X_E_DENT(x,5))           ! E_Dipole entrance
         y1 = DBLE(Y_E_DENT(x,5))
      else
         x1 = DBLE(X_H_DENT(x,5))           ! H_Dipole entrance
         y1 = DBLE(Y_H_DENT(x,5))
      endif
      ylim = dip2(ispect)*(1.d0-(dip3(ispect)*x1/dip4(ispect)))
      if((abs(x1) .gt. dip1(ispect) .or. abs(y1) .gt. ylim)
     &               .and. aperture_test(ispect,2)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y
         edip_ent_x = x1
         edip_ent_y = y1
      else
         hdip_ent_x = x1
         hdip_ent_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(X_E_DEXT(x,5))           ! E_Dipole exit
         y1 = DBLE(Y_E_DEXT(x,5))
      else
         x1 = DBLE(X_H_DEXT(x,5))           ! H_Dipole exit
         y1 = DBLE(Y_H_DEXT(x,5))
      endif
      ylim = dip2(ispect)*(1.d0-(dip3(ispect)*x1/dip4(ispect)))
      if((abs(x1) .gt. dip1(ispect) .or. abs(y1) .gt. ylim)
     &               .and. aperture_test(ispect,3)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y
         edip_ext_x = x1
         edip_ext_y = y1
      else
         hdip_ext_x = x1
         hdip_ext_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(X_E_Q3EN(x,5))           ! E_Q3 entrance
         y1 = DBLE(Y_E_Q3EN(x,5))
      else
         x1 = DBLE(X_H_Q3EN(x,5))           ! H_Q3 entrance
         y1 = DBLE(Y_H_Q3EN(x,5))
      endif
      if(x1*x1+y1*y1 .gt. q3_radius(ispect)*q3_radius(ispect)
     &               .and. aperture_test(ispect,4)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y
         eq3_ent_x = x1
         eq3_ent_y = y1
      else
         hq3_ent_x = x1
         hq3_ent_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(X_E_Q3EX(x,5))           ! E_Q3 exit
         y1 = DBLE(Y_E_Q3EX(x,5))
      else
         x1 = DBLE(X_H_Q3EX(x,5))           ! H_Q3 exit
         y1 = DBLE(Y_H_Q3EX(x,5))
      endif
      if(x1*x1+y1*y1 .gt. q3_radius(ispect)*q3_radius(ispect)
     &               .and. aperture_test(ispect,5)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! for histograms of X and Y
         eq3_ext_x = x1
         eq3_ext_y = y1
      else
         hq3_ext_x = x1
         hq3_ext_y = y1
      endif
C
C -----------------------------------------------------------------------
C     Get focal plane Transport vector.
C -----------------------------------------------------------------------
C
      if(ispect .eq. 1) then
         x1 = DBLE(X_E_FP(x,5))             ! Electron arm
         t1 = DBLE(T_E_FP(x,5))
         y1 = DBLE(Y_E_FP(x,5))
         p1 = DBLE(P_E_FP(x,5))
      else
         x1 = DBLE(X_H_FP(x,5))             ! Hadron arm
         t1 = DBLE(T_H_FP(x,5))
         y1 = DBLE(Y_H_FP(x,5))
         p1 = DBLE(P_H_FP(x,5))
      endif
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
C     Overall path is L=L0+x_fp(5) where L0=24.73 m
C     Assume same calculation for both arms for the moment.
C -----------------------------------------------------------------------
C
      x_fp(5) = 100.d0 * DBLE(3.216389*x(1)
     &                      - 5.150459*x(4)
     &                      - 0.580672*x(5)) + x_tg(5) 
C
 123  RETURN
      END
      







