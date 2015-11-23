	subroutine transp(spectr,class)
!--------------------------------------------------------------------------------
!
! TRANSP - This subroutine transports a particle through the various
!   segments of a spectrometer. The passed variable CLASS determines
!   which transformation to use. Each transformation IS SEQUENTIAL,
!   and carries the particle from the last transformation to a particular
!   plane Z=const in the spectrometer.
!
!   NOTE: The coordinate system used here is the right handed "TRANSPORT"
!   coordinate system, in which +Z -> downstream, +X -> in median plane,
!   pointing in direction taken through bending magnet by high momentum rays,
!   +Y = transverse direction.
!
!   For an upward vertical bend spectrometer:
!
!   +X points towards the floor,
!   +Y points horizontally LEFT as one looks downstream.
!
! D. Potterveld - March 1993.
!
!   Modification History:
!
! August 1993.	(D. Potterveld) Modified to begin each transformation at
!		the pivot.
!
!  20-AUG-1993	(D. Potterveld) Modified to use COSY INFINITY forward maps.
!
!  03-MAR-1994	(D. Potterveld) Fixed bug to correctly compute cosy "A" and "B".
!
!  10-MAY-1995	(D. Potterveld) Switched to sequential transformations.
!
!  11-MAY-1995	(D. Potterveld) Added rejection of extremely small coeff.
!
!  30-OCT-1995	(D. Potterveld) Change to COSY-7 transport units.
!
!  16-SEP-1998  Check for decay of particle. dflag is true of the particle
!		has already decayed, so check for decay if dflag .eq. .false.
!		After pathlength is calculated, add extra check for decay
!		within the different path length.
!
!--------------------------------------------------------------------------------

	implicit 	none

	INCLUDE 'track.cmn'

! Arguments.

	integer*4       spectr			!HMS=1, SOS=2, HRSR=3, HRSL=4
	integer*4	class
	character*(*)	file			!Used in entry point below.

! Parameters.

	DOUBLE PRECISION coeff_min
	parameter	(coeff_min = 1.0d-14)	!coeff's smaller than this = 0

! Local declarations.

	integer*4	idummy,iterm

	DOUBLE PRECISION		delta_z		!pathlength difference.
	integer*4	i,j,k,kk,chan,n_classes

	character*132	str_line
	character*132	file_name

! Cosy reconstruction matrix elements.

	integer*4	max_elements
	parameter	(max_elements = 500)
	DOUBLE PRECISION coeff(nspectr,5,max_elements,max_class)
	integer*2	expon(nspectr,5,max_elements,max_class)
	integer*4	n_terms(nspectr,max_class)
	DOUBLE PRECISION		sum(5),ray(5),term,temp
	DOUBLE PRECISION		length(nspectr,max_class)
	integer*2	order
	integer*2	e1,e2,e3,e4                     !temp. exponants
	DOUBLE PRECISION	c1,c2,c3,c4,csum	!temp. coeffs.


! Function definitions.

	logical		locforunt
	real*8 grnd
        COMMON /DATDIR_I/ N_DAT_DIR
        COMMON /DATDIR_C/ DAT_DIR 
        INTEGER N_DAT_DIR 
	CHARACTER*100 DAT_DIR


! No amnesia between calls!!!

	save

! ================================ Executable Code =============================

! A word from the sponsor:
!  a) Using sequential matrix elements, variables with the ending 's' should be
!     used, otherwise they should end on '_transp'
!  b) the sequential matrix elements have been made with COSY 7, i.e. the units

! to be used are cm, mrad, not m, slopes

! Check that the path length passed to code (zd) is consistent with comments
! in MEs (if they exist)



! Pack local copy of input coordinates.

	  ray(1) = xs             	!cm.	( "X" )
	  ray(2) = dxdzs*1000.   	!mrad.	( "THETA" )
	  ray(3) = ys			!cm.	( "Y" )
	  ray(4) = dydzs*1000.          !mrad.	( "PHI" )
	  ray(5) = dpps			! Fractional "Delta P/P"

! Reset COSY sums.

	do i = 1,5
	   sum(i) = 0.
	enddo

! Compute COSY sums.

	k = class
	do i = 1,n_terms(spectr,k)
	  term = 1.0
	  do j = 1,5
	    temp = 1.0
	    if (expon(spectr,j,i,k).ne.0.) 
     >  	 temp = ray(j)**expon(spectr,j,i,k)
	    term = term*temp
	  enddo
	  sum(1) = sum(1) + term*coeff(spectr,1,i,k)		! NEW "X"
	  sum(2) = sum(2) + term*coeff(spectr,2,i,k)		! NEW "A"
	  sum(3) = sum(3) + term*coeff(spectr,3,i,k)		! NEW "Y"
	  sum(4) = sum(4) + term*coeff(spectr,4,i,k)		! NEW "B"
	  sum(5) = sum(5) + term*coeff(spectr,5,i,k)		! NEW "dL"
	enddo

! Unpack output coordinates. Note that DPPS is unchanged by transformation.
! Pathlength correction: real pathlength=nominal-sum(5), so delta_z=-sum(5)

	  xs    = sum(1)                                !cm
	  dxdzs = sum(2)/1000.         		!slope (mr)
	  ys    = sum(3)                         	!cm
	  dydzs = sum(4)/1000.                        !slope (mr)



	return

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

! Initialization entry points.

	entry transp_init_file(file,n_classes)

! Use passed filename.

	file_name = file
	goto 100

	entry transp_init(spectr,n_classes)

	  file_name= DAT_DIR(1:N_DAT_DIR)//'/hrs_forward_cosy.dat'

! Open input file.

100     if (.not.locforunt(chan))
     >  stop 'TRANSP_INIT: No I/O channels!'
	open (unit=chan,status='old',name=file_name)

! Strip away header.

	str_line = '!'
	do while (str_line(1:1).eq.'!')
	  read (chan,1001) str_line
	  if (str_line(1:8).eq.'!LENGTH:') then
	    read (str_line(10:),*) length(spectr,1)          !get length from comments
	    length(spectr,1)=100.*length(spectr,1)           !convert to cm.
	  endif
	enddo

! Read in the transformation tables.

	n_classes = 0
	do i = 1,max_class
	  n_terms(spectr,i) = 0
	  adrift(spectr,i) = .true.
	  driftdist(spectr,i) = 0.0
	enddo

	do while (.true.)
	  kk = n_classes + 1

! If too many transformations, complain!

	  if (kk.gt.max_class) 
     >      stop 'TRANSP_INIT: too many transformations!'

! Add data lines to table, looking for flag line.

	  do while (str_line(1:4).ne.' ---')
	    n_terms(spectr,kk) = n_terms(spectr,kk) + 1
	    if (n_terms(spectr,kk).gt.max_elements)
     >      stop 'TRANSP_INIT: too many COSY terms!'

! Read in MEs.  Coefficients are for calculating : X, XP, Y, YP, dZ
!		Exponents determine the powers of: X, XP, Y, YP, dZ, Delta

	    read (str_line,1200)
     >		(coeff(spectr,i,n_terms(spectr,kk),kk),i=1,5),
     >		(expon(spectr,j,n_terms(spectr,kk),kk),j=1,4),idummy,
     >		 expon(spectr,5,n_terms(spectr,kk),kk)

! Ignore time-of-flight term.

	    if (idummy.ne.0) then
	      if (coeff(spectr,1,n_terms(spectr,kk),kk).ne.0.0.or.
     >		  coeff(spectr,2,n_terms(spectr,kk),kk).ne.0.0.or.
     >		  coeff(spectr,3,n_terms(spectr,kk),kk).ne.0.0.or.
     >		  coeff(spectr,4,n_terms(spectr,kk),kk).ne.0.0)
     >			stop 'TRANSP_INIT: non-zero TOF terms!'
	      n_terms(spectr,kk) = n_terms(spectr,kk) - 1
	    endif

! Check to see if element is inconsistent with a field-free drift.
! 1st order terms require: diagonal terms have coeff=1
! 			   <x|xp> and <y|yp> terms have coeff=distance(m)/10
!			   other coeff=0
! Other order terms require all coeff=0

	    if (idummy.eq.0. .and. adrift(spectr,kk)) then
	      iterm = n_terms(spectr,kk)
	      e1 = expon(spectr,1,iterm,kk)
	      e2 = expon(spectr,2,iterm,kk)
	      e3 = expon(spectr,3,iterm,kk)
	      e4 = expon(spectr,4,iterm,kk)
	      c1 = coeff(spectr,1,iterm,kk)
	      c2 = coeff(spectr,2,iterm,kk)
	      c3 = coeff(spectr,3,iterm,kk)
	      c4 = coeff(spectr,4,iterm,kk)
	      csum = abs(c1)+abs(c2)+abs(c3)+abs(c4)
	      order = e1 + e2 + e3 + e4

	      if (order.eq.1) then
	        if (e1.eq.1) then
		  if (abs(c1-1.) .gt. coeff_min) 
     >    		adrift(spectr,kk)=.false.
		  if (abs(c2-0.) .gt. coeff_min) 
     >                  adrift(spectr,kk)=.false.
		  if (abs(c3-0.) .gt. coeff_min) 
     >                  adrift(spectr,kk)=.false.
		  if (abs(c4-0.) .gt. coeff_min) 
     >                 adrift(spectr,kk)=.false.
		else if (e2.eq.1) then
		  driftdist(spectr,kk)=1000.0*c1     !drift distance in cm.
		  if (abs(c2-1.) .gt. coeff_min) 
     >                adrift(spectr,kk)=.false.
		  if (abs(c3-0.) .gt. coeff_min) 
     >                adrift(spectr,kk)=.false.
		  if (abs(c4-0.) .gt. coeff_min) 
     >                adrift(spectr,kk)=.false.

		else if (e3.eq.1) then
		  if (abs(c1-0.) .gt. coeff_min) 
     >                adrift(spectr,kk)=.false.
		  if (abs(c2-0.) .gt. coeff_min) 
     >                adrift(spectr,kk)=.false.
		  if (abs(c3-1.) .gt. coeff_min) 
     >                adrift(spectr,kk)=.false.
		  if (abs(c4-0.) .gt. coeff_min) 
     >                adrift(spectr,kk)=.false.
		else if (e4.eq.1) then
		  if (abs(c1-0.) .gt. coeff_min) 
     >                adrift(spectr,kk)=.false.
		  if (abs(c2-0.) .gt. coeff_min) 
     >                adrift(spectr,kk)=.false.

		  if (abs(driftdist(spectr,kk)-1000.*c3).gt.coeff_min) 
     >					     adrift(spectr,kk)=.false. 
	  
		  if (abs(c4-1.) .gt. coeff_min) 
     >                adrift(spectr,kk)=.false.

		endif
	      else	!if order.ne.1
		if (abs(csum).gt.coeff_min) adrift(spectr,kk)=.false.
	      endif
	    endif

! Get next line from file.

	    read (chan,1001) str_line
	  enddo

! If flag line is seen, increment transformation counter.

	  n_classes = kk

	  if (adrift(spectr,kk) .and.
     >		abs(driftdist(spectr,kk)-length(spectr,kk))
     >                                  .gt.0.01) then
	    write(6,*) 
     >      'PROBLEM WITH TRANSFORMATION #',kk,' for spectrometer #',
     >      spectr
	    write(6,*) 
     >      'Appears to be a pure drift, driftdist=',
     >      driftdist(spectr,kk)
	    write(6,*) 
     >      'But the comments say that the length =',
     >      length(spectr,kk)
	  endif


! Read lines until a non-blank, non-comment non-terminal line is found.

150       read (chan,1001,end=200) str_line
	  if (str_line(1:8).eq.'!LENGTH:') then
	    read (str_line(10:),*) length(spectr,kk+1)                      !get length from comments
	    length(spectr,kk+1)=100.0*length(spectr,kk+1)                  !convert to cm.
	  endif
	  if (str_line(1:1).eq.'!'.or.str_line(1:4).eq.' ---'.or.
     >    str_line.eq.'    ') goto 150

	enddo


! Done with file.

200     close (unit=chan)



	return

! ============================== Format Statements =============================
1001    format(a)
1101    format(12x,6e11.3)
1102    format(10x,5f10.5)
1200    format(1x,5g14.7,1x,6i1)

	end






