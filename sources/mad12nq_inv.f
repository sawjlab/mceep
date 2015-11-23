C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE MAD12NQ_INV
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    14-DEC-2004
C       PURPOSE: Produces target Transport vector for the
C                JLAB MAD spectrometer (12 deg configuration
C                with QUADS turned OFF) given the focal plane 
C                Transport vector.
C
C       USES inverse functions of J. LeRose to map spectrometer
C            coordinates.
C
C       INPUT:   spec_angle:      Spect. (in-plane) central angle (rad)
C                x_beam:          vertical beam position (m)
C                y_beam:          horizontal beam position (m)
C                x_fp:            Focal plane Transport vector
C       OUTPUT:  x_tg:            Target Transport vector
C
C       Currently, spec_angle and y_beam are not used, though
C       they are passed in case they are eventually needed for the
C       "bowtie" (or depth-of-field) correction.
C
C       Apparently, txfit is not needed in the no quads tune.
C
C -----------------------------------------------------------------------
C
      SUBROUTINE MAD12NQ_INV(spec_angle,x_beam,y_beam,x_fp,x_tg)
      IMPLICIT NONE
C
C
      DOUBLE PRECISION x_tg(6),x_fp(6),spec_angle,x_beam,y_beam
C
      REAL x(5)
      REAL delta_nq,theta_nq,y00_nq,phi_nq
C
      x(1) = REAL(x_fp(1)/100.d0)        ! cm --> meters
      x(2) = REAL(tan(x_fp(2)/1000.d0))  ! tan(theta)
      x(3) = REAL(x_fp(3)/100.d0)        ! cm --> meters
      x(4) = REAL(tan(x_fp(4)/1000.d0))  ! tan(phi)
      x(5) = REAL(x_beam)                ! already in meters
C
C -----------------------------------------------------------------------
C     Get target Transport vector.
C -----------------------------------------------------------------------
C
      x_tg(2) = dble(theta_nq(x,5))
      x_tg(3) = dble(y00_nq(x,5))
      x_tg(4) = dble(phi_nq(x,5))
      x_tg(6) = dble(delta_nq(x,5))
C
      x_tg(5) = 0.d0                ! set path length to zero for now
      x_tg(1) = 0.d0                ! set vertical pos. to zero for now
C
C -----------------------------------------------------------------------
C     Change units back to MCEEP conventions.
C -----------------------------------------------------------------------
C
      x_tg(1) = x_tg(1)*100.d0          ! m --> cm
      x_tg(2) = atan(x_tg(2))*1000.d0   ! mr
      x_tg(3) = x_tg(3)*100.d0          ! m --> cm
      x_tg(4) = atan(x_tg(4))*1000.d0   ! mr
      x_tg(5) = x_tg(5)*100.d0          ! m --> cm
      x_tg(6) = x_tg(6)*100.d0          ! fraction --> %
C
 999  RETURN
      END
