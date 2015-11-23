C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE HRS_INV
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    20-APR-1999
C       PURPOSE: Produces target Transport vector for the
C                JLAB Hall A spectrometers given the focal plane
C                Transport vector.
C
C       05-NOV-1999    PEU
C            Modified to perform "bowtie" correction (as of MCEEP
C            Version 3.1).  This correction requires the spectrometer
C            central in-plane angle (now passed) and corrects for
C            the non-zero X0 "seen by the spectrometer" due to the
C            reference Transport system not being along the target.
C
C       07-MAR-2000    PEU
C            Fixed error in "bowtie" correction and now also correct
C            for explicit X0 due to raster, mispointing, etc.
C            (as of Version 3.4).
C
C       19-FEB-2004    PEU
C            Add calculation for x_tg(5), the path length.  Prior
C            to this, it was set to zero.
C            (new as of Version 3.9).
C
C
C       USES inverse functions of J. LeRose to map spectrometer
C            coordinates.
C
C       INPUT:   spec_angle:      Spect. (in-plane) central angle (rad)
C                x_fp:            Focal plane Transport vector
C       OUTPUT:  x_tg:            Target Transport vector
C
C -----------------------------------------------------------------------
C
      SUBROUTINE HRS_INV(ispect,spec_angle,x_beam,y_beam,x_fp,x_tg)
      IMPLICIT NONE
C
C
      DOUBLE PRECISION x_tg(6),x_fp(6),spec_angle
      DOUBLE PRECISION x_x0(2),x_delta(2),theta_x0(2)
      DOUBLE PRECISION theta_theta(2),theta_delta(2)
      DOUBLE PRECISION cangle,t_theta0,sphisum
      DOUBLE PRECISION delta_corr,theta_corr
      DOUBLE PRECISION x_beam,y_beam
      REAL x(4)
      REAL E_TXFIT,E_THETA,E_Y00,E_PHI,E_DELTA
      REAL H_TXFIT,H_THETA,H_Y00,H_PHI,H_DELTA
      INTEGER ispect
C
C -----------------------------------------------------------------------
C     The following DATA are needed to correct for the "bowtie" effect.
C -----------------------------------------------------------------------
C    
      DATA x_x0        /-2.181d0,-2.170d0/   ! <x|x0>
      DATA x_delta     /11.905d0,11.890d0/   ! <x|delta>     (cm/%)
cxxx  DATA theta_x0    /-1.000d0,-1.010d0/   ! <theta|x0>    (mr/cm)
      DATA theta_x0    /-0.503d0,-0.485d0/   ! <theta|x0>    (mr/cm)
      DATA theta_theta /-0.47d0,-0.47d0/     ! <theta|theta>
      DATA theta_delta /19.33d0,19.28d0/     ! <theta|delta> (mr/%)
C
C -----------------------------------------------------------------------
C     Index defining which spectrometer is being used:
C         ispect = 1  (electron arm)
C         ispect = 2  (hadron arm)
C     Note that this allows either arm to be used for detection of
C     either particle (i.e. "electron" arm could be used to detect
C     hadrons) since this index can be set to 1 or 2 regardless of
C     the value of IARM.
C
C
C     Repack the array and change units to match LeRose's conventions.
C -----------------------------------------------------------------------
C
      x(1) = REAL(x_fp(1)/100.d0)        ! cm --> meters
      x(2) = REAL(tan(x_fp(2)/1000.d0))  ! tan(theta)
      x(3) = REAL(x_fp(3)/100.d0)        ! cm --> meters
      x(4) = REAL(tan(x_fp(4)/1000.d0))  ! tan(phi)
C
C -----------------------------------------------------------------------
C     Get target Transport vector.
C -----------------------------------------------------------------------
C
      if(ispect .eq. 1) then
         x(2) = x(2) - e_txfit(x,1)  ! theta now wrt local central ray
c
         x_tg(2) = dble(e_theta(x,4))
         x_tg(3) = dble(e_y00(x,4))
         x_tg(4) = dble(e_phi(x,4))
         x_tg(6) = dble(e_delta(x,4))
      else
         x(2) = x(2) - h_txfit(x,1)  ! theta now wrt local central ray
c
         x_tg(2) = dble(h_theta(x,4))
         x_tg(3) = dble(h_y00(x,4))
         x_tg(4) = dble(h_phi(x,4))
         x_tg(6) = dble(h_delta(x,4))
      endif
C
C -----------------------------------------------------------------------
C     Change units back to MCEEP conventions.
C -----------------------------------------------------------------------
C
      t_theta0 = x_tg(2)    ! save this for bowtie correction
C
      x_tg(2) = atan(x_tg(2))*1000.d0   ! mr
      x_tg(3) = x_tg(3)*100.d0          ! m --> cm
      x_tg(4) = atan(x_tg(4))*1000.d0   ! mr
      x_tg(6) = x_tg(6)*100.d0          ! fraction --> %
C
C -----------------------------------------------------------------------
C     Perform "bowtie" (i.e. depth of field) correction.
C     This correction accounts for the fact that an extended target
C     can give rise to an effective X0 relative to the spectrometer
C     Transport reference plane.
C
C     Also correct for explicit X0 due to vertical beam raster, beam
C     vertical position offset/spread or spectrometer vertical
C     mispointing.
C
C     I add the corrections, where LeRose would subtract them.  This
C     difference arises from our opposite definition of the sign of the
C     spectrometer central angle.  My definition is that + angles are
C     on the left side of the beam.  This same difference in convention
C     means my sphisum reflects the sum of the spectrometer angle
C     and phi_tg, rather than the difference.
C -----------------------------------------------------------------------
C
      cangle   = cos(spec_angle)
      sphisum  = sin(spec_angle+x_tg(4)/1000.d0)
C
      x_tg(1) = x_beam*100.d0 + t_theta0
     #          * (x_tg(3)*cangle-y_beam*100.d0) / sphisum
C
      theta_corr = (((theta_delta(ispect)*x_x0(ispect)/x_delta(ispect))
     #      - theta_x0(ispect)) / theta_theta(ispect))*x_tg(1)
C
C
C -----------------------------------------------------------------------
C     Iterate with improved theta0
C -----------------------------------------------------------------------
C
      t_theta0 = tan((x_tg(2) + theta_corr)/1000.d0)
      x_tg(1) = x_beam*100.d0 + t_theta0
     #          * (x_tg(3)*cangle-y_beam*100.d0) / sphisum
C
      theta_corr = (((theta_delta(ispect)*x_x0(ispect)/x_delta(ispect))
     #      - theta_x0(ispect)) / theta_theta(ispect))*x_tg(1)
C
      delta_corr = -(    x_x0(ispect) /     x_delta(ispect))*x_tg(1)
C
      x_tg(6) = x_tg(6) + delta_corr
      x_tg(2) = x_tg(2) + theta_corr
C
C -----------------------------------------------------------------------
C     Calculate the path length by simply inverting
C     the formula found in hrs.f (i.e. the forward map).  Assume the
C     same functional form for both arms.  Prior to v3.9, x_tg(5) was
C     set to zero.
C -----------------------------------------------------------------------
C
      x_tg(5) = x_fp(5)
     &            - 100.d0 * DBLE( 3.216389 * (x_tg(1)/100.d0)
     &                           - 5.150459 * tan(x_tg(4)/1000.d0)
     &                           - 0.580672 * (x_tg(6)/100.d0) )
C
 999  RETURN
      END
