C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE SEP
C
C       AUTHOR:  Pete Markowitz 
C       DATE:    08-MAY-2006
C       PURPOSE: Produces focal plane Transport vector for the
C                JLAB Hall A spectrometers + septa given the target Transport
C                vector.  Also gives the value of the logical flag
C                fail_apertures consistent with the size of various
C                apertures.
C
C       MODS:    27-APR-2006 PECM
C                Based on v3.9 hrs.f of 25-Mar-1999 by Paul Ulmer
C                
C                
C
C       USES functions of J. LeRose to map septum+spectrometer coordinates.
C
C       INPUT:   x_tg:            Target Transport vector
C       OUTPUT:  x_fp:            Focal plane Transport vector
C                fail_apertures:  True if any aperture is not passed.
C
C -----------------------------------------------------------------------
C
      SUBROUTINE SEPTUM(ispect,x_tg,x_fp,fail_apertures)
      IMPLICIT NONE
C
      DOUBLE PRECISION x_tg(6),x_fp(6)
      DOUBLE PRECISION yep(2),xep3l(2),xep3u(2)
      DOUBLE PRECISION xep4l(2),xep4u(2),xep5l(2),xep5u(2)
      DOUBLE PRECISION xep6l(2),xep6u(2),xep7l(2),xep7u(2)
      DOUBLE PRECISION q1_radius(2),q3_radius(2)
      DOUBLE PRECISION dip1u(2),dip1l(2),dip2(2),dip3(2),dip4(2),dip5(2)
      DOUBLE PRECISION x1,t1,y1,p1,ylim
      REAL X(5)
      REAL x_sl_ep3, y_sl_ep3, x_sr_ep3, y_sr_ep3
      REAL x_sl_ep4, y_sl_ep4, x_sr_ep4, y_sr_ep4
      REAL x_sl_ep5, y_sl_ep5, x_sr_ep5, y_sr_ep5
      REAL x_sl_ep6, y_sl_ep6, x_sr_ep6, y_sr_ep6
      REAL x_sl_ep7, y_sl_ep7, x_sr_ep7, y_sr_ep7
      REAL X_SL_Q1EX,Y_SL_Q1EX,X_SL_DENT,Y_SL_DENT,X_SL_DEXT,Y_SL_DEXT
      REAL X_SL_Q3EN,Y_SL_Q3EN,X_SL_Q3EX,Y_SL_Q3EX
      REAL X_SL_FP,T_SL_FP,Y_SL_FP,P_SL_FP
      REAL X_SR_Q1EX,Y_SR_Q1EX,X_SR_DENT,Y_SR_DENT,X_SR_DEXT,Y_SR_DEXT
      REAL X_SR_Q3EN,Y_SR_Q3EN,X_SR_Q3EX,Y_SR_Q3EX
      REAL X_SR_FP,T_SR_FP,Y_SR_FP,P_SR_FP
      INTEGER ispect,i
      LOGICAL fail_apertures
C
      INCLUDE 'spectrometer.cmn'
      INCLUDE 'hrs.cmn'
      INCLUDE 'sep.cmn'
C
C -----------------------------------------------------------------------
C     Parameters describing apertures.
C        Septum (ep3,4,5,6,7): rectangular apertures
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
C     ep3: -0.1486233 < x < -0.08869672
c          -0.110 < y < 0.110
c     ep4: -0.1792231 < x < -0.1089169
c          -0.110 < y < 0.110
c     ep5: -0.2209211 < x < -0.1353789
c          -0.110 < y < 0.110
c     ep6: -0.2763536 < x < -0.1697464
c          -0.110 < y < 0.110
c     ep7: -0.3485396 < x < -0.2156404
c          -0.110 < y < 0.110
c     q1ex is a circle of radius 0.1492 m
c     dent is a trapazoid:
c                                   -5.22008 < x < -4.98099
c             -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784  
c     dext is also a trapazoid: 
c                                   -0.46188 < x < 0.46188
c                   -(-0.01610808*x + 0.125) < y < -0.01610808*x + 0.125
c     q3en is a circle of radius 0.3 m
c     q3ex is a circle of radius 0.3 m
cC -----------------------------------------------------------------------
C
      DATA yep /0.110d0,0.110d0/ ! -.0110<y<0.0110 throughout septum 
      DATA xep3l /-0.1486233d0,-0.1486233d0/ ! xep3 lower bound septum entrance
      DATA xep4l /-0.1792231d0,-0.1792231d0/ ! ep4 lower 1/4 through septum
      DATA xep5l /-0.2209211d0,-0.2209211d0/ ! ep5 lower 1/2 through septum
      DATA xep6l /-0.2763536d0,-0.2763536d0/ ! ep6 lower 3/4 through septum
      DATA xep7l /-0.3485396d0,-0.3485396d0/ ! ep7 lower septum exit
      DATA xep3u /-0.08869672d0,-0.08869672d0/ ! xep3 upper 
      DATA xep4u /-0.1089169d0,-0.1089169d0/ ! ep4 upper
      DATA xep5u /-0.1353789d0,-0.1353789d0/ ! 
      DATA xep6u /-0.1697464d0,-0.1697464d0/ ! 
      DATA xep7u /-0.2156404d0,-0.2156404d0/ ! 

      DATA q1_radius /0.1492d0,0.1492d0/ ! circle (meters)
Cv3.6      DATA q3_radius /0.280d0,0.280d0/   ! circle (meters)
      DATA q3_radius /0.300d0,0.300d0/   ! circle (meters)
C
      DATA dip1u /-4.98099d0,-4.98099d0/         ! trapezoid x lower entrance
      DATA dip1l /-5.22008d0,-5.22008d0/         ! trapezoid x upper entrance
      DATA dip2 /-0.192436784d0,-0.192436784d0/   ! trapezoid y entrance
      DATA dip3 /-.01610808d0,-.01610808d0/
      DATA dip4 /125.d0,125.d0/
      DATA dip5 /0.46188d0,0.46188d0/
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
         x1 = DBLE(x_sl_ep3(x,5))           ! e=l left septum near entrance
         y1 = DBLE(y_sl_ep3(x,5))
      else
         x1 = DBLE(x_sr_ep3(x,5))           ! h=r right septum near entrance
         y1 = DBLE(y_sr_ep3(x,5))
      endif

      if((abs(y1) .gt. yep(ispect) .or. x1 .lt. xep3l(ispect) .or.
     &      x1 .gt. xep3u(ispect) ) .and. aper_sep_test(ispect,1)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y
         sl_ep3_x = x1
         sl_ep3_y = y1
      else
         sr_ep3_x = x1
         sr_ep3_y = y1
      endif
      if(ispect .eq. 1) then
         x1 = DBLE(x_sl_ep4(x,5))           ! e=l left septum 1/4 through
         y1 = DBLE(y_sl_ep4(x,5))
      else
         x1 = DBLE(x_sr_ep4(x,5))           ! h=r right septum 1/4 through
         y1 = DBLE(y_sr_ep4(x,5))
      endif

      if((abs(y1) .gt. yep(ispect) .or. x1 .lt. xep4l(ispect) .or.
     &      x1 .gt. xep4u(ispect) ) .and. aper_sep_test(ispect,2)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y
         sl_ep4_x = x1
         sl_ep4_y = y1
      else
         sr_ep4_x = x1
         sr_ep4_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(x_sl_ep5(x,5))           ! e=l left septum 1/2 through
         y1 = DBLE(y_sl_ep5(x,5))
      else
         x1 = DBLE(x_sr_ep5(x,5))           ! h=r right septum 1/2 through
         y1 = DBLE(y_sr_ep5(x,5))
      endif

      if((abs(y1) .gt. yep(ispect) .or. x1 .lt. xep5l(ispect) .or.
     &      x1 .gt. xep5u(ispect) ) .and. aper_sep_test(ispect,3)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y
         sl_ep5_x = x1
         sl_ep5_y = y1
      else
         sr_ep5_x = x1
         sr_ep5_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(x_sl_ep6(x,5))           ! e=l left septum 3/4 through
         y1 = DBLE(y_sl_ep6(x,5))
      else
         x1 = DBLE(x_sr_ep6(x,5))           ! h=r right septum 3/4 through
         y1 = DBLE(y_sr_ep6(x,5))
      endif

      if((abs(y1) .gt. yep(ispect) .or. x1 .lt. xep6l(ispect) .or.
     &      x1 .gt. xep6u(ispect) ) .and. aper_sep_test(ispect,4)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y
         sl_ep6_x = x1
         sl_ep6_y = y1
      else
         sr_ep6_x = x1
         sr_ep6_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(x_sl_ep7(x,5))           ! e=l left septum exit
         y1 = DBLE(y_sl_ep7(x,5))
      else
         x1 = DBLE(x_sr_ep7(x,5))           ! h=r right septum exit
         y1 = DBLE(y_sr_ep7(x,5))
      endif

      if((abs(y1) .gt. yep(ispect) .or. x1 .lt. xep7l(ispect) .or.
     &      x1 .gt. xep7u(ispect) ) .and. aper_sep_test(ispect,5)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y
         sl_ep7_x = x1
         sl_ep7_y = y1
      else
         sr_ep7_x = x1
         sr_ep7_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(x_sl_q1ex(x,5))           ! E_Q1 exit left==electron
         y1 = DBLE(y_sl_q1ex(x,5))
      else
         x1 = DBLE(x_sr_q1ex(x,5))           ! H_Q1 exit
         y1 = DBLE(y_sr_q1ex(x,5))
      endif
      if(x1*x1+y1*y1 .gt. q1_radius(ispect)*q1_radius(ispect)
     &               .and. aper_sep_test(ispect,6)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y -- kept same names
         eq1_ext_x = x1           ! as without septum
         eq1_ext_y = y1
      else
         hq1_ext_x = x1
         hq1_ext_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(X_SL_DENT(x,5))           ! E_Dipole entrance
         y1 = DBLE(Y_SL_DENT(x,5))
      else
         x1 = DBLE(X_SR_DENT(x,5))           ! H_Dipole entrance
         y1 = DBLE(Y_SR_DENT(x,5))
      endif
      ylim = dip2(ispect)*(1.d0 + x1)
      if((x1 .lt. dip1l(ispect) .or. x1 .gt. dip1u(ispect) .or. 
     &   abs(y1) .gt. ylim) .and. aper_sep_test(ispect,7)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y -- same as w/o septum
         edip_ent_x = x1
         edip_ent_y = y1
      else
         hdip_ent_x = x1
         hdip_ent_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(X_SL_DEXT(x,5))           ! E_Dipole exit
         y1 = DBLE(Y_SL_DEXT(x,5))
      else
         x1 = DBLE(X_SR_DEXT(x,5))           ! H_Dipole exit
         y1 = DBLE(Y_SR_DEXT(x,5))
      endif
      ylim = dip3(ispect)*x1 + dip4(ispect)
      if((abs(x1) .gt. dip5(ispect) .or. abs(y1) .gt. ylim)
     &               .and. aper_sep_test(ispect,8)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y -- same as w/o septum
         edip_ext_x = x1
         edip_ext_y = y1
      else
         hdip_ext_x = x1
         hdip_ext_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(X_SL_Q3EN(x,5))           ! E_Q3 entrance
         y1 = DBLE(Y_SL_Q3EN(x,5))
      else
         x1 = DBLE(X_SR_Q3EN(x,5))           ! H_Q3 entrance
         y1 = DBLE(Y_SR_Q3EN(x,5))
      endif
      if(x1*x1+y1*y1 .gt. q3_radius(ispect)*q3_radius(ispect)
     &               .and. aper_sep_test(ispect,9)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then      ! histograms of X and Y -- same as w/o septum
         eq3_ent_x = x1
         eq3_ent_y = y1
      else
         hq3_ent_x = x1
         hq3_ent_y = y1
      endif
C
      if(ispect .eq. 1) then
         x1 = DBLE(X_SL_Q3EX(x,5))           ! E_Q3 exit
         y1 = DBLE(Y_SL_Q3EX(x,5))
      else
         x1 = DBLE(X_SR_Q3EX(x,5))           ! H_Q3 exit
         y1 = DBLE(Y_SR_Q3EX(x,5))
      endif
      if(x1*x1+y1*y1 .gt. q3_radius(ispect)*q3_radius(ispect)
     &               .and. aper_sep_test(ispect,10)) then
         fail_apertures = .true.
         goto 123
      endif
C
      if(ispect .eq. 1) then     ! histograms of X and Y -- same as w/o septum
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
         x1 = DBLE(X_SL_FP(x,5))             ! Electron arm
         t1 = DBLE(T_SL_FP(x,5))
         y1 = DBLE(Y_SL_FP(x,5))
         p1 = DBLE(P_SL_FP(x,5))
      else
         x1 = DBLE(X_SR_FP(x,5))             ! Hadron arm
         t1 = DBLE(T_SR_FP(x,5))
         y1 = DBLE(Y_SR_FP(x,5))
         p1 = DBLE(P_SR_FP(x,5))
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
C     ??Can't figure out above -- just add 79 cm for now          pecm
C -----------------------------------------------------------------------
C
      x_fp(5) = 100.d0 * DBLE(3.216389*x(1)
     &                      - 5.150459*x(4)
     &                      - 0.580672*x(5)) + x_tg(5) + 79.d0 
C
 123  RETURN
      END
      







