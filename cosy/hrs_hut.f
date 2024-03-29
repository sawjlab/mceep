	subroutine hrs_hut (m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		ok_hut,zinit)

C----------------------------------------------------------------------
C
C Monte-Carlo of HRSR detector hut.
C
C	The particle is stepped through the detector (using project), and
C	multiple scattering is applied for each detector or air gap.
C	If particle decay in enabled, then project.f also checks for
C	decay of particle.  The particle starts at z=zinit.  This
C	needs to be before the first mult. scattering point (the exit window)
C	or the decay distance is negative, and things are BAD.
C
C----------------------------------------------------------------------

	implicit 	none

	INCLUDE 'struct_hrs.cmn'
	INCLUDE 'track.cmn'

C Math constants

	real*8 pi,d_r,r_d

	parameter (pi = 3.141592654)
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)

c	real*8 gauss1			!external functions

C all parameters, later to take from .parm files
C----------------------------------------------------------------------
C HRSR_MATERIALS
C CTP parameter file containing the materials of all the HRSR detectors.
C For all materials AFTER the bend only have to do multiple scattering.
C     radlen = 1 radiation length (in cm)
C     thick  = thickness in cm
C In case a "+" is added to the comment, the thickness is a guess only.
C----------------------------------------------------------------------
C spectrometer exit window, .1 mm of titanium (X0=3.56cm)
	DOUBLE PRECISION hfoil_exit_radlen,hfoil_exit_thick
	parameter (hfoil_exit_radlen = 3.56d0)
	parameter (hfoil_exit_thick  = 0.00010d0*2.54d0)

C spectrometer air gaps
	DOUBLE PRECISION hair_radlen
	parameter (hair_radlen = 30420.d0)

C chamber entrance foil, .18 mm of Kapton Foil (or polyimide film)
	DOUBLE PRECISION hdc_entr_radlen,hdc_entr_thick
	parameter (hdc_entr_radlen = 34.4d0)
	parameter (hdc_entr_thick  = 0.00018d0*2.54d0)

C chamber gas, 35/65 ethane/argon
	DOUBLE PRECISION hdc_radlen,hdc_thick
	parameter (hdc_radlen = 16700.d0)
!	parameter (hdc_thick  = 1.8)
	parameter (hdc_thick  = 1.5d0)	!made up number, but want 1/2 thickness
					! < 5 cm to avoid negative drifts.

C effective wire plane, 25 micron W+Be/Cu gives <t>=pi*(.0025/2)**2
	DOUBLE PRECISION hdc_wire_radlen,hdc_wire_thick
	parameter (hdc_wire_radlen = 0.35)	!Assuming all Tungsten
	parameter (hdc_wire_thick  = 0.0000049)

C effective cathode plane, Be/Cu
	DOUBLE PRECISION hdc_cath_radlen,hdc_cath_thick
	parameter (hdc_cath_radlen = 7.2)	!'Ave' of Be/Cu
	parameter (hdc_cath_thick  = 0.000177)

C chamber exit foil, .18 mm of Kapton Foil (or polyimide film)
	DOUBLE PRECISION hdc_exit_radlen,hdc_exit_thick
	parameter (hdc_exit_radlen = 34.4)
	parameter (hdc_exit_thick  = 0.00018*2.54)

C hodoscopes
	DOUBLE PRECISION hscin_radlen
	parameter (hscin_radlen =  42.4)

C Cherenkov entrance foil, 40 mil of Tedlar Film
	DOUBLE PRECISION hcer_entr_radlen,hcer_entr_thick
	parameter (hcer_entr_radlen = 8.90)
	parameter (hcer_entr_thick  = 0.040*2.54)

C Cherenkov, 0.5 atm of CO2 for better e/pi
	DOUBLE PRECISION hcer_radlen
	parameter (hcer_radlen = 36620.0)

C Cherenkov, 2 atm of Freon for pi/p
C	DOUBLE PRECISION hcer_radlen
C	parameter (hcer_radlen = 2405.0)

C Cherenkov mirror, 75 mu plus 2 cm of Rotacell +
	DOUBLE PRECISION hcer_mir_radlen,hcer_mir_thick
	parameter (hcer_mir_radlen = 400.0)
	parameter (hcer_mir_thick  = 2.0)

C Cherenkov exit foil
	DOUBLE PRECISION hcer_exit_radlen,hcer_exit_thick
	parameter (hcer_exit_radlen = 8.90)
	parameter (hcer_exit_thick  = 0.040*2.54)

C shower counter
	DOUBLE PRECISION hcal_radlen
	parameter (hcal_radlen = 14.83)

C Wire chamber resolutions (sigma)

	DOUBLE PRECISION hdc_sigma(1:12)
     >    /0.0225,0.0225,0.0225,0.0225,0.0225,0.0225,
     >	   0.0225,0.0225,0.0225,0.0225,0.0225,0.0225/

C Wire plane positions, construct hdc_zpos array using these parameters

	integer*4 hdc_nr_cham,hdc_nr_plan
	parameter (hdc_nr_cham = 2)
	parameter (hdc_nr_plan = 6)

	DOUBLE PRECISION hdc_1_zpos,hdc_1_left,hdc_1_right,
     >    hdc_1_top,hdc_1_bot
	DOUBLE PRECISION hdc_1x_offset,hdc_1y_offset
	DOUBLE PRECISION hdc_2_zpos,hdc_2_left,hdc_2_right,
     >   hdc_2_top,hdc_2_bot
	DOUBLE PRECISION hdc_2x_offset,hdc_2y_offset

	DOUBLE PRECISION hdc_del_plane

C Drift chamber 1 is the focal plane, so shift all zpos values by 25cm
	parameter (hdc_1_zpos = -25.d0 + 25.d0)
	parameter (hdc_2_zpos =  25.0 + 25.0)
	parameter (hdc_del_plane = hdc_thick + 
     >             hdc_wire_thick + hdc_cath_thick)
	parameter (hdc_1_left  =  14.4)
	parameter (hdc_1_right = -14.4)
	parameter (hdc_1y_offset = 0.000)
	parameter (hdc_1_top   = -105.6)
	parameter (hdc_1_bot   =  105.6)
	parameter (hdc_1x_offset = 0.000)
	parameter (hdc_2_left  =  14.4)
	parameter (hdc_2_right = -14.4)
	parameter (hdc_2y_offset = 0.000)
	parameter (hdc_2_top   = -105.6)
	parameter (hdc_2_bot   =  105.6)
	parameter (hdc_2x_offset = 0.000)

C Scintillator positions and thicknesses

	DOUBLE PRECISION hscin_1x_zpos
	DOUBLE PRECISION hscin_1x_thick
	DOUBLE PRECISION hscin_1x_left,hscin_1x_right,hscin_1x_offset
	DOUBLE PRECISION hscin_2x_zpos
	DOUBLE PRECISION hscin_2x_thick
	DOUBLE PRECISION hscin_2x_left,hscin_2x_right,hscin_2x_offset

	parameter (hscin_1x_zpos =  95.0 + 25.0)
	parameter (hscin_2x_zpos = 288.3 + 25.0)
	parameter (hscin_1x_thick = 0.5*1.067) ! 1.067 for overlap
	parameter (hscin_2x_thick = 0.5*1.067)
	parameter (hscin_1x_left  =  18.0)
	parameter (hscin_1x_right = -18.0)
	parameter (hscin_1x_offset = 0.00)     !up-down offset to hscin_1y.
	parameter (hscin_2x_left  =  30.0)
	parameter (hscin_2x_right = -30.0)
	parameter (hscin_2x_offset = 0.00)


C Cherenkov position

	DOUBLE PRECISION hcer_zentrance,hcer_zmirror,hcer_zexit
	parameter (hcer_zentrance = 137.0 + 25.0)
	parameter (hcer_zmirror   = 197.0 + 25.0)
	parameter (hcer_zexit     = 237.0 + 25.0)

C Calorimeter position

	DOUBLE PRECISION hcal_1pr_zpos,hcal_2ta_zpos,hcal_3ta_zpos,
     >                   hcal_4ta_zpos
	DOUBLE PRECISION hcal_left,hcal_right,hcal_top,hcal_bottom
	parameter (hcal_1pr_zpos = 302.3 + 25.0)
	parameter (hcal_2ta_zpos = 337.3 + 25.0)
	parameter (hcal_3ta_zpos = 372.3 + 25.0)
	parameter (hcal_4ta_zpos = 407.3 + 25.0)
	parameter (hcal_left     =  7.5)
	parameter (hcal_right    = -7.5)
	parameter (hcal_top      = -7.5)
	parameter (hcal_bottom   =  7.5)

C The arguments

	DOUBLE PRECISION	p,m2			!momentum and mass of particle
	DOUBLE PRECISION	xt,yt			!temporary variables.
	DOUBLE PRECISION	x_fp,y_fp,dx_fp,dy_fp	!Focal plane values to return
	DOUBLE PRECISION	xcal,ycal		!Position of track at calorimeter.
	DOUBLE PRECISION	zinit			!Initial z-position (Not at F.P.)
	DOUBLE PRECISION	pathlen
	logical ms_flag			!mult. scattering flag.
	logical wcs_flag		!wire chamber smearing flag
	logical decay_flag		!check for decay
	logical ok_hut			!true if particle makes it

C Local declarations.

	integer*4 i,iplane,jchamber,npl_off

	logical dflag				!has particle decayed?

	DOUBLE PRECISION	resmult
	DOUBLE PRECISION	tmpran1,tmpran2			!temporary random numbers
	DOUBLE PRECISION	radw,drift
	DOUBLE PRECISION	random(4)			!space for 4 random numbers

	DOUBLE PRECISION nsig_max
	parameter(nsig_max=99.0d0)	!max #/sigma for gaussian ran #s.

C These have to be real*4 for the CERNLIB lfit routine

	real*4	badf				!temporaries
	real*4  xfp4,yfp4,dxfp4,dyfp4           !real*4 versions of fp track
	real*4	xdc(12),ydc(12),zdc(12)		!positions at d.c. planes
        DOUBLE PRECISION PROJECT_XS, PROJECT_YS
C ================================ Executable Code =============================

C Initialize ok_hut to zero

	ok_hut = .false.

C Initialize the xdc and ydc arrays to zero

	do i=1,12
	  xdc(i) = 0.
	  ydc(i) = 0.
	enddo

C------------------------------------------------------------------------------C
C                           Top of loop through hut                            C
C------------------------------------------------------------------------------C

C Scatter in spectrometer exit foil, which is located at zinit.
C As usual, neglect effect of nonzero dydzs and dxdzs on radw.

	radw = hfoil_exit_thick/hfoil_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C Go to first drift chamber set
C For simplicity, perform air MS (probably negligeable) at before drift
C instead of 1/2 way through.

	drift = (hdc_1_zpos - 0.5*hdc_nr_plan*hdc_del_plane) - zinit
	radw = drift/hair_radlen
	if (ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(PROJECT_xs,PROJECT_ys,DRIFT)
        XS=PROJECT_XS
        YS=PROJECT_YS

	jchamber = 1
	radw = hdc_entr_thick/hdc_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	npl_off = (jchamber-1)*hdc_nr_plan
	do iplane = 1,hdc_nr_plan
	  radw = hdc_cath_thick/hdc_cath_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_cath_thick
	  call project(PROJECT_xs,PROJECT_ys,drift)
          XS=PROJECT_XS
          YS=PROJECT_YS
	  radw = hdc_wire_thick/hdc_wire_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  if(wcs_flag) then
             call gaussians(random(1),random(2))
	  else
             random(1) = 0.
             random(2) = 0.
	  endif
	  xdc(npl_off+iplane) = xs + hdc_sigma(npl_off+iplane)*random(1)
	  ydc(npl_off+iplane) = ys + hdc_sigma(npl_off+iplane)*random(2)
	  if (iplane.eq.1 .or. iplane.eq.3 .or. iplane.eq.5) then
	    xdc(npl_off+iplane) = 0.   !y plane, no x information
	  else
	    ydc(npl_off+iplane) = 0.   !x-like plane, no y info
	  endif

	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_wire_thick
	  call project(PROJECT_xs,PROJECT_ys,drift)
          XS=PROJECT_XS
          YS=PROJECT_YS
	enddo
	radw = hdc_exit_thick/hdc_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

!rotate 45 degrees to compare to VDCs.  CHECK SIGN AND SIZE OF ROTATAION!!!
	xt=xs
	yt=ys
	call rotate_haxis(45.0d0,xt,yt)

	if (xt.gt.(hdc_1_bot-hdc_1x_offset) .or.
     >      xt.lt.(hdc_1_top-hdc_1x_offset) .or.
     >      yt.gt.(hdc_1_left-hdc_1y_offset) .or.
     >      yt.lt.(hdc_1_right-hdc_1y_offset) ) then
	  goto 500
	endif
	radw = hdc_cath_thick/hdc_cath_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C at last cathode foil of first drift chamber set, drift to next

	drift = hdc_2_zpos - hdc_1_zpos - hdc_nr_plan*hdc_del_plane
	call project(PROJECT_xs,PROJECT_ys,drift)
        XS=PROJECT_XS
        YS=PROJECT_YS
	radw = drift/hair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	jchamber = 2
	radw = hdc_entr_thick/hdc_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	npl_off = (jchamber-1)*hdc_nr_plan
	do iplane = 1,hdc_nr_plan
	  radw = hdc_cath_thick/hdc_cath_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_cath_thick
	  call project(PROJECT_xs,PROJECT_ys,drift)
          XS=PROJECT_XS
          YS=PROJECT_YS
	  radw = hdc_wire_thick/hdc_wire_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  if(wcs_flag) then
             call gaussians(random(1),random(2))
	  else
	    random(1) = 0.
	    random(2) = 0.
	  endif
	  xdc(npl_off+iplane) = XS + hdc_sigma(npl_off+iplane)*random(1)
	  ydc(npl_off+iplane) = YS + hdc_sigma(npl_off+iplane)*random(2)
	  if (iplane.eq.1 .or. iplane.eq.3 .or. iplane.eq.5) then
	    xdc(npl_off+iplane) = 0.   !y plane, no x information
	  else
	    ydc(npl_off+iplane) = 0.   !x-like plane, no y info
	  endif

	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_wire_thick
	  call project(PROJECT_xs,PROJECT_ys,drift)
          XS=PROJECT_XS
          YS=PROJECT_YS
	enddo
	radw = hdc_exit_thick/hdc_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

!rotate 45 degrees to compare to VDCs.  CHECK SIGN AND SIZE OF ROTATAION!!!
        xt=xs
        yt=ys
        call rotate_haxis(45.0d0,xt,yt)

	if (xt.gt.(hdc_2_bot-hdc_2x_offset) .or.
     >      xt.lt.(hdc_2_top-hdc_2x_offset) .or.
     >      yt.gt.(hdc_2_left-hdc_2y_offset) .or.
     >      yt.lt.(hdc_2_right-hdc_2y_offset) ) then
	  goto 500
	endif
	radw = hdc_cath_thick/hdc_cath_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C fit track to give new focal plane values, use LFIT from GENLIB

	do jchamber=1,hdc_nr_cham
	  npl_off = (jchamber-1)*hdc_nr_plan
	  do iplane=1,hdc_nr_plan
	    if (jchamber.eq.1) zdc(npl_off+iplane) = hdc_1_zpos +
     >          (iplane-0.5-0.5*hdc_nr_plan)*hdc_del_plane
	    if (jchamber.eq.2) zdc(npl_off+iplane) = hdc_2_zpos +
     >          (iplane-0.5-0.5*hdc_nr_plan)*hdc_del_plane
	  enddo
	enddo

	call lfit(zdc,xdc,12,0,dxfp4,xfp4,badf)
	call lfit(zdc,ydc,12,0,dyfp4,yfp4,badf)

	x_fp = dble(xfp4)
	y_fp = dble(yfp4)
	dx_fp = dble(dxfp4)
	dy_fp = dble(dyfp4)

C at last cathode foil of second drift chamber set, drift to hodoscopes

	drift = hscin_1x_zpos - hdc_2_zpos - 
     >          0.5*hdc_nr_plan*hdc_del_plane
	radw = drift/hair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(PROJECT_xs,PROJECT_ys,drift)
        XS=PROJECT_XS
        YS=PROJECT_YS
	if (ys.gt.(hscin_1x_left) .or.
     >      ys.lt.(hscin_1x_right)) then
	  goto 500
	endif
	radw = hscin_1x_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C finished first hodoscope, drift to cherenkov

	drift = hcer_zentrance - hscin_1x_zpos
	radw = drift/hair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(PROJECT_xs,PROJECT_ys,drift)
        XS=PROJECT_XS
        YS=PROJECT_YS

	radw = hcer_entr_thick/hcer_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = hcer_zmirror - hcer_zentrance
	radw = drift/hcer_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(PROJECT_xs,PROJECT_ys,drift)
        XS=PROJECT_XS
        YS=PROJECT_YS

	radw = hcer_mir_thick/hcer_mir_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = hcer_zexit - hcer_zmirror
	radw = drift/hcer_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(PROJECT_xs,PROJECT_ys,drift)
        XS=PROJECT_XS
        YS=PROJECT_YS

	radw = hcer_exit_thick/hcer_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C drift to second hodoscope

	drift = hscin_2x_zpos - hcer_zexit
	radw = drift/hair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(PROJECT_xs,PROJECT_ys,drift)
        XS=PROJECT_XS
        YS=PROJECT_YS

	if (ys.gt.(hscin_2x_left) .or.
     >      ys.lt.(hscin_2x_right)) then
	  goto 500
	endif
	radw = hscin_2x_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C Don't need to drift to calorimeter unless it's required in your trigger.
C Note that even with the standard PID trigger, the calorimeter is NOT
C required, since the trigger needs either the cerenkov OR the calorimeter.
C There is a seperate fiducial cut needed if you require the calorimeter
C in you analysis.  That cut is applied AFTER fitting the track (see below).

*	drift = hcal_4ta_zpos - hscin_2y_zpos
*	radw = drift/hair_radlen
*	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
*	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)
*	if (ys.gt.hcal_left .or. ys.lt.hcal_right .or.
*     >	   xs.gt.hcal_bottom .or. xs.lt.hcal_top) then
*	  rSTOP.cal = rSTOP.cal + 1
*	  stop_where=22.
*	  x_stop=xs
*	  y_stop=ys
*	  goto 500
*	endif

C If you use a calorimeter cut in your analysis, the engine applied a
C a fiducial cut at the calorimeter.  This is based on the position of the
C TRACK at the calorimeter, not the real position of the event.  Go to the
C back of the calorimeter since engine uses a FID cut at the back.
C The standard fiducial cut is 5 cm from the edges of the block.

	xcal = x_fp + dx_fp * hcal_4ta_zpos
	ycal = y_fp + dy_fp * hcal_4ta_zpos
*	if (ycal.gt.(hcal_left-5.0) .or. ycal.lt.(hcal_right+5.0) .or.
*     >	   xcal.gt.(hcal_bottom-5.0) .or. xcal.lt.(hcal_top+5.0)) then
*	  rSTOP.cal = rSTOP.cal + 1
*	  stop_where=23.
*	  x_stop=xs
*	  y_stop=ys
*	  goto 500
*	endif

	ok_hut = .true.

C We are done with this event, whether GOOD or BAD.

 500   continue

C ALL done!

	return
	end
