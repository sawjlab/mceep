	subroutine rotate_haxis(rotang,xp0,yp0)

!-------------------------------------------------------------------------------------
!
! ROTATE_HAXIS- Calculate new trajectory coordinates in reference frame rotated
!   about horizontal axis by angle ROTANG relative to central ray.
!
! *** Right-handed TRANSPORT coordinates are assumed! ***
!
!   ROTANG is an angle about the  negative Y-axis.
!   For the SOS BM01 entrance, it is a negative number.
!
!   Input trajectory is: X = XS + ALPHA*(Z-ZS)
!                        Y = YS + BETA *(Z-ZS)
!                        Z = ZS is current point.
!
!   Output traject is:  XP = XP0 + ALPHA_P*ZP
!                       YP = YP0 + BETA_P *ZP
!                       ZP = 0 gives intersection of track with rotated plane.
!
! ROTANG (R*4):	Rotation angle in degrees.
!
! D. Potterveld, 15-Mar-1993.
!-------------------------------------------------------------------------------------

	implicit none

	INCLUDE 'track.cmn'

	DOUBLE PRECISION rotang,xp0,yp0,xi
	DOUBLE PRECISION alpha,beta,alpha_p,beta_p,sin_th,cos_th,tan_th

	save

! ============================= Executable Code ================================

	tan_th = tan(rotang*3.14159/180.d0)
	sin_th = sin(rotang*3.14159/180.d0)
	cos_th = cos(rotang*3.14159/180.d0)

	alpha  = dxdzs
	beta   = dydzs

	alpha_p= (alpha + tan_th)/(1.d0 - alpha*tan_th)
	beta_p = beta/(cos_th - alpha*sin_th)

	xp0    = xs*(cos_th + alpha_p*sin_th)
	yp0    = ys + xs*beta_p*sin_th

	return
	end


