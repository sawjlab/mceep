      SUBROUTINE HRS_COSY(p_spec, dpp, x_tg, y_tg, z_tg, dxdz, dydz,
     >		x_fp,dx_fp,y_fp,dy_fp,
     >	        m2, ms_flag, wcs_flag, col_flag,fry, hrs,ok_hrse)


! Monte-Carlo of HRS spectrometer.
!

! Author: David Meekins April 2000
! based on code by David Potterveld
!
! Modification History:
!
!  units for this file are percents, cm, and mrads.
!
C-______________________________________________________________________________

      IMPLICIT    NONE 

        INCLUDE 'struct_hrs.cmn'
	INCLUDE 'apertures.cmn'
	INCLUDE 'track.cmn'

	DOUBLE PRECISION x_offset_pipes/0.0/,y_offset_pipes/0.0/
! Spectrometer definitions - for double arm monte carlo compatability
	INTEGER*4 HRS	
        INTEGER*4 HRSE
! Collimator (rectangle) dimensions and offsets.

	DOUBLE PRECISION   Z_OFF		!offset in distance from target to front of slit
	PARAMETER (Z_OFF  = 0.d0)

! z-position of important apertures.
	DOUBLE PRECISION  z_entr,z_exit
	PARAMETER (z_entr = 110.9d0 + z_off)	!nominally 1.109 m
	PARAMETER (z_exit = z_entr + 8.0d0)	!8.0 cm thick

! Math constants

	DOUBLE PRECISION pi,d_r,r_d,root
	PARAMETER (pi = 3.141592654)
	PARAMETER (d_r = pi/180.)
	PARAMETER (r_d = 180./pi)
	PARAMETER (root = 0.707106781)		!square root of 1/2

! Masses.

	DOUBLE PRECISION me,mp,md,mpi,mka
	DOUBLE PRECISION me2,mp2,md2,mpi2,mka2
	PARAMETER (me   = 0.511E-03)
	PARAMETER (mp   = 0.9382796d0)
	PARAMETER (md   = 1.875628d0)
	PARAMETER (mpi  = 0.1395675d0)
	PARAMETER (mka  = 0.493646d0)

! The arguments
	DOUBLE PRECISION		x_tg,y_tg,z_tg		       !(cm)
	DOUBLE PRECISION		dpp     		       !delta p/p (%)
	DOUBLE PRECISION 		dxdz,dydz		       !X,Y slope in spectrometer
	DOUBLE PRECISION		x_fp,y_fp,dx_fp,dy_fp	       !Focal plane values to return
	DOUBLE PRECISION		p_spec			       !spectrometer setting
	DOUBLE PRECISION		tg_rad_len		       !target length in r.l.
	DOUBLE PRECISION                fry                            !vertical position@tgt (+y=down)
	LOGICAL*4	                ms_flag, wcs_flag,col_flag     !particle, m_scat, wc_smear
	LOGICAL*4	                ok_hrse,OK			       !true if particle makes it
 
! Local declarations.

	INTEGER*4	chan/1/,n_classes
	LOGICAL*4	first_time_hrse/.true./
	DOUBLE PRECISION dpp_recon,dth_recon,dph_recon	               !reconstructed quantities
	DOUBLE PRECISION y_recon
	DOUBLE PRECISION p				               !More kinematic variables.
	DOUBLE PRECISION xt,yt				               !temporaries
	DOUBLE PRECISION m2				               !particle mass
	DOUBLE PRECISION ZDRIFT,ztmp				

	LOGICAL*4 fit_reverse	/.false./	                       !Flag for fitting reverse coeff.

        DOUBLE PRECISION PROJECT_xs,PROJECT_YS		               !track x,y,z positions (cm)
        DOUBLE PRECISION x_new,y_new

	LOGICAL check_dipole
	EXTERNAL check_dipole  

	save		!Remember it all!

C ================================ Executable Code =============================


! Initialize some stuff
        HRSE = 1
	ME2 = ME*ME 
	MP2 = MP*MP
	MD2 = MD*MD
	MPI2 = MPI*MPI
	MKA2 = MKA*MKA
        IF (HRS.EQ.1) THEN 
           M2 = ME2
        ELSE
           M2 = MP2
        ENDIF   
! Initialize ok_hrse to zero
	OK_HRSE =.FALSE.

! particle momentum
        DPPS = DPP
	P = P_SPEC * (1.d0 + DPPS/100.d0)*1000.d0

! and the rest....
	XS = X_TG                      !cm 
	YS = Y_TG                      !cm 
	ZS = Z_TG                      !cm 
        DXDZS = DXDZ                   !mrad
	DYDZS = DYDZ                   !mrad 

! Read in transport coefficients.
	IF (FIRST_TIME_HRSE) THEN
	  CALL TRANSP_INIT(HRSE,N_CLASSES)
	  CLOSE (UNIT=CHAN)
	  if (n_classes.ne.12) 
     >         stop 'MC_HRSL, wrong number of transport classes'
	  FIRST_TIME_HRSE = .FALSE.
	ENDIF
C------------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

C Begin transporting particle.


C Do transformations, checking against apertures.
! Circular apertures before slitbox (only important for no collimator)
	ZDRIFT = 65.686
	ztmp = ZDRIFT

	CALL PROJECT(PROJECT_XS,PROJECT_YS,ZDRIFT)
        XS = PROJECT_XS
        YS = PROJECT_YS
 
	IF (sqrt(XS*XS+YS*YS).GT.7.3787) THEN
	  GOTO 500
	ENDIF

	ZDRIFT = 80.436 - ztmp
	ztmp = 80.436

	CALL PROJECT(PROJECT_XS,PROJECT_YS,ZDRIFT)
        XS = PROJECT_XS
        YS = PROJECT_YS 
	IF (sqrt(XS*XS+YS*YS).GT.7.4092) THEN
           GOTO 500
        ENDIF

!ENTER COLLIMATOR
        If(col_flag)then
! Check front of fixed slit, at about 1.109 meter
           ZDRIFT = Z_ENTR-ZTMP

           call PROJECT(PROJECT_XS,PROJECT_YS,ZDRIFT)
           XS =  PROJECT_XS
           YS =  PROJECT_YS

           IF (ABS(YS-Y_OFFSET).GT.H_ENTR) THEN
              GOTO 500
           ENDIF

           IF (ABS(XS-X_OFFSET).GT.V_ENTR) THEN
              GOTO 500
           ENDIF

         
! Check back of fixed slit, at about 1.189 meter

           ZDRIFT = Z_EXIT-Z_ENTR
           CALL PROJECT(PROJECT_XS,PROJECT_YS,ZDRIFT)
           XS = PROJECT_XS
           YS = PROJECT_YS

           IF (ABS(YS-Y_OFFSET).GT.H_EXIT) THEN
              GOTO 500
           ENDIF
           IF (ABS(XS-X_OFFSET).GT.V_EXIT) THEN
              GOTO 500
           ENDIF
 
        ENDIF
! Aperture before Q1 (can only check this if next transformation is DRIFT).

        ztmp = 135.064
        ZDRIFT = ztmp - z_exit
        CALL PROJECT(PROJECT_XS,PROJECT_YS,ZDRIFT) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS

        IF (sqrt(XS*XS+YS*YS).GT.12.5222) THEN
           GOTO 500
        ENDIF

! Go to Q1 IN  mag bound.
	IF (.not.adrift(HRSE,1)) 
     >       write(6,*) 'Transformation #1 is NOT a drift'
	ZDRIFT = driftdist(HRSE,1) - ztmp
	CALL PROJECT(PROJECT_XS,PROJECT_YS,ZDRIFT) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS
	IF ((XS*XS + YS*YS).GT.r_Q1*r_Q1) THEN
	  GOTO 500
	ENDIF

! Check aperture at 2/3 of Q1.
          
        CALL TRANSP(HRSE,2)
	IF ((XS*XS + YS*YS).GT.R_Q1*R_Q1) THEN
	   GOTO 500
	ENDIF
	


! Go to Q1 OUT mag boundary.

	   CALL TRANSP(HRSE,3)
	   IF ((XS*XS + YS*YS).GT.R_Q1*R_Q1) THEN
	      GOTO 500
	   ENDIF


! Apertures after Q1, before Q2 (can only check this if next trans. is DRIFT).
	  
	zdrift = 300.464 - 253.16		!Q1 exit is z=253.16
	ztmp = zdrift				!distance from Q1 exit
	CALL PROJECT(PROJECT_XS,PROJECT_YS,zdrift) !PROJECT
        XS = PROJECT_XS
        YS = PROJECT_YS

 	IF (sqrt(XS*XS+YS*YS).GT.14.9225) THEN
           GOTO 500
	ENDIF

	zdrift = 314.464 - 300.464
	ztmp = ztmp + zdrift			!distance from Q1 exit.
	CALL PROJECT(PROJECT_XS,PROJECT_YS,zdrift) !PROJECT
        XS = PROJECT_XS
        YS = PROJECT_YS 
	IF (sqrt(XS*XS+YS*YS).GT.20.9550) THEN
	  GOTO 500
	ENDIF

! Go to Q2 IN  mag bound.

	IF (.NOT.adrift(HRSE,4)) 
     >       write(6,*) 'Transformation #4 is NOT a drift'
	zdrift = driftdist(HRSE,4) - ztmp
	CALL PROJECT(PROJECT_XS,PROJECT_YS,zdrift) !PROJECT
        XS = PROJECT_XS
        YS = PROJECT_YS
	IF ((XS*XS + YS*YS).GT.r_Q2*r_Q2) THEN
	  GOTO 500
	ENDIF


! Check aperture at 2/3 of Q2.

	   CALL TRANSP(HRSE,5)
	   IF ((XS*XS + YS*YS).GT.R_Q2*R_Q2) THEN
	      GOTO 500
	   ENDIF




! Go to Q2 OUT mag boundary.

	   CALL TRANSP(HRSE,6)
	   IF ((XS*XS + YS*YS).GT.R_Q2*R_Q2) THEN
	      GOTO 500
	   ENDIF




! Apertures after Q2, before D1 (can only check this if next trans. is DRIFT).

        zdrift = 609.664 - 553.020 !Q2 exit is z=553.02
	ztmp = zdrift				!distance from Q2 exit
	call PROJECT(PROJECT_XS,PROJECT_YS,zdrift) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS 
	if (sqrt(XS*XS+YS*YS).gt.30.0073) then
           goto 500
	endif

	zdrift = 641.800 - 609.664
	ztmp = ztmp + zdrift			!distance from Q2 exit.
	call PROJECT(PROJECT_XS,PROJECT_YS,zdrift) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS 
	if (sqrt(XS*XS+YS*YS).gt.30.0073) then
	  goto 500
	endif

	zdrift = 819.489 - 641.800
	ztmp = ztmp + zdrift			!distance from Q2 exit.
	call PROJECT(PROJECT_XS,PROJECT_YS,zdrift) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS 
	if (abs(XS).gt.50.0 .or. abs(YS).gt.15.0) then
	  goto 500
      endif


! Go to D1 IN magnetic boundary.
! Find intersection with rotated aperture plane.
! Aperture has elliptical form.
! Go to D1 IN magnetic boundary, Find intersection with rotated aperture plane.

!	if (.not.adrift(hrse,7)) write(6,*) 'Transformation #7 is NOT a drift/I made it this far'
	zdrift = driftdist(HRSE,7) - ztmp
        call PROJECT(PROJECT_XS,PROJECT_YS,zdrift) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS 
	xt=XS
	yt=YS
	call rotate_haxis(-30.0,xt,yt)
	if (abs(xt-2.500).gt.52.5) then		! -50 < x < +55
	  goto 500
	endif
	if ( (abs(yt)+0.01861*xt) .gt. 12.5 ) then !tan(1.066) ~ 0.01861
	  goto 500
	endif



! Go to D1 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	   CALL TRANSP(HRSE,8)
	     xt=XS
             yt=YS
	   CALL ROTATE_HAXIS(30.0,XT,YT)
           IF (abs(xt-2.500).GT.52.5) THEN ! -50 < x < +55
	      GOTO 500
	   ENDIF

            IF ( (abs(yt)+0.01861*xt) .GT. 12.5 ) THEN !tan(1.066) ~ 0.01861
	     GOTO 500
	   ENDIF 




! Apertures after D1, before Q3 (can only check this if next trans. is DRIFT).

	zdrift = 1745.33546 - 1655.83446	!D1 exit is z=1655.83446
	ztmp = zdrift				!distance from D1 exit
	CALL PROJECT(PROJECT_XS,PROJECT_YS,zdrift) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS 

	IF (sqrt(XS*XS+YS*YS).GT.30.3276) THEN
	  GOTO 500
	ENDIF


	IF (abs(XS).GT.50.0 .OR. abs(YS).GT.15.0) THEN
	  GOTO 500
	ENDIF

	zdrift = 1759.00946 - 1745.33546
	ztmp = ztmp + zdrift			!distance from D1 exit

	CALL PROJECT(PROJECT_XS,PROJECT_YS,zdrift) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS 

	IF (sqrt(XS*XS+YS*YS).GT.30.3276) then
	  GOTO 500
	ENDIF



! Go to Q3 IN  mag bound.

	if (.not.adrift(hrse,9)) 
     >       write(6,*) 'Transformation #9 is NOT a drift'
	zdrift = driftdist(hrse,9) - ztmp
	call PROJECT(PROJECT_XS,PROJECT_YS,zdrift) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS 
	if ((XS*XS + YS*YS).gt.r_Q3*r_Q3) then
	  goto 500
	endif



! Check apeture at 2/3 of Q3
	call transp(HRSE,10)
	IF ((XS*XS + YS*YS).gt.r_Q3*r_Q3) then
	  goto 500
	endIF


! Go to Q3 OUT mag boundary.

	   CALL TRANSP(HRSE,11)
	   IF ((XS*XS + YS*YS).GT.R_Q3*R_Q3) THEN
	      GOTO 500
	   ENDIF

! Apertures after Q3 (can only check this if next trans. is DRIFT).

	ZDRIFT = 2080.38746 - 1997.76446	!Q3 exit is z=1997.76446
	ztmp = ZDRIFT				!distance from Q3 exit
	call PROJECT(PROJECT_XS,PROJECT_YS,ZDRIFT) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS 
	IF (abs(XS).gt.35.56 .or. abs(YS).gt.17.145) then

	  goto 500
	endIF

 

! Vacuum window is 15.522cm before FP (which is at VDC1)
	zdrift = 2327.47246 - 2080.38746
	ztmp = ztmp + zdrift			!distance from Q3 exit
	call PROJECT(PROJECT_XS,PROJECT_YS,ZDRIFT) !PROJECT 
        XS = PROJECT_XS
        YS = PROJECT_YS 
	if (abs(XS).gt.99.76635 .or. abs(YS).gt.17.145) then
	  goto 500
	endif

! If we get this far, the particle is in the hut.

!	if (.not.adrift(hrse,12)) write(6,*) 'Transformation #12 is NOT a drift'

	zdrift = driftdist(hrse,12) - ztmp    !distance left to go.

        CALL HRS_HUT(M2,P,X_FP,DX_FP,Y_FP,DY_FP,MS_FLAG,
     >         WCS_FLAG,OK,-zdrift)
        IF (.NOT.OK) GOTO 500

! replace XS,YS,... with 'tracked' quantities. 
          XS    = X_FP
          YS    = Y_FP
          DXDZS = DX_FP
          DYDZS = DY_FP

! Apply offset to y_fp (detectors offset w.r.t optical axis).
! In ESPACE, the offset is taken out for recon, but NOT for y_fp in ntuple,
! so we do not apply it to YS (which goes to recon), but do shift it for y_fp.
! IF THIS IS TRUE OFFSET, WE SHOULD SHIFT DETECTOR APERTURES - NEED TO CHECK!!!!
! But in general the dectectors don't limit the acceptance, so we should be OK

	y_fp = y_fp -0.48		!VDC center is +4.8mm.




C Reconstruct target quantities.

        CALL HRS_RECON(DPP_RECON,DTH_RECON,DPH_RECON,Y_RECON,fry)


C Fill output to return to main code 
        DXDZ = DPH_RECON
        DYDZ = DTH_RECON
        Y_TG = Y_RECON
        X_TG = 0.d0
	dpp = dpp_recon
        OK_HRSE = .TRUE.
C We are done with this event, whether GOOD or BAD.

 500	   CONTINUE

C ALL done!

	END




















