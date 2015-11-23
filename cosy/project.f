	subroutine project(X_NEW,Y_NEW,Z_DRIFT)

!-------------------------------------------------------------------------------------
!
! PROJECT - Calculate new track transverse coordinates after drifting in a field
!   free region for a distance z_drift from current track location.
!-------------------------------------------------------------------------------------

	implicit	none

	INCLUDE 'track.cmn'

	DOUBLE PRECISION X_NEW,Y_NEW,Z_DRIFT, DXDZS_PROJECT, 
     >                   DYDZS_PROJECT
!	DOUBLE PRECISION NEW
	save
	

	DXDZS_PROJECT = DXDZS
	DYDZS_PROJECT = DYDZS

! ============================= Executable Code ================================

	X_NEW = XS + Z_DRIFT * DXDZS
	Y_NEW = YS + Z_DRIFT * DYDZS

	return
	end
