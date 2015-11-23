	subroutine hrs_recon (delta_p,delta_t,delta_phi,y_tgt,fry)

!-----------------------------------------------------------------------------
! MC_HRSE_RECON : Reconstruct target quantities from tracks.
!		   This subroutine is part of the HRS program.
!
! Right-handed coordinates are assumed: X=down, Z=downstream, Y = (Z cross X)
!
! Author: D. Potterveld, ANL, 18-Mar-1993
!
! Modification History:
!
! 2-August-1995 (RMM) Hardcoded in an extra index on the coeff, expon, and
!                      n_terms variables to specify HRSE. (HRSE = 1)
!
!  19-AUG-1993	(DHP) Modified to use COSY INFINITY reconstruction coefficients.
!   January 2001 - Units are modified in the following manner
!	m --> cm
!	rad
!	frac variation --> %
!	
!----------------------------------------------------------------------------
     
	implicit none

	INCLUDE 'track.cmn'

	INTEGER*4 hrs
	parameter (hrs = 1)			!this is the HRS routine

! Argument definitions.

	DOUBLE PRECISION  delta_p,delta_t,delta_phi,y_tgt
	DOUBLE PRECISION  fry      !vertical position at the target (+y=down)

! Cosy reconstruction matrix elements.

	INTEGER*4	  max_elements
	parameter	  (max_elements = 1000)
	DOUBLE PRECISION  coeff(nspectr,4,max_elements)
	INTEGER*2	  expon(nspectr,5,max_elements)
	INTEGER*4	  n_terms(nspectr),max_order
	DOUBLE PRECISION  sum(4),hut(5),term

! Misc. variables.

	INTEGER*4	i,j
	INTEGER*4	chan
	logical*4	firsthrs	/.true./
	character*132	line
	character*132	file_name

! Functions.

	INTEGER*4	last_char
	logical*4	locforunt
	logical*4	firsttime	/.true./

! 
        COMMON /DATDIR_I/ N_DAT_DIR
        COMMON /DATDIR_C/ DAT_DIR 
        INTEGER N_DAT_DIR 
	CHARACTER*100 DAT_DIR
	save

! ============================= Executable Code ================================


! First time through, read in coefficients from data file.
! Use default filename.
	file_name = DAT_DIR(1:N_DAT_DIR)//'/hrs_recon_cosy.dat'

	if (firsttime) then
	   if (.not.locforunt(chan)) 
     >        stop 'MC_HRS_RECON: No I/O channels!'
           OPEN(UNIT=CHAN,FILE=file_name,STATUS='OLD',
     #     FORM='FORMATTED')

! Skip past header.

	   line = '!'
	   do while (line(1:1).eq.'!')
	      read (chan,1001) line
	   enddo

! Read in coefficients and exponents.
	   n_terms(hrs) = 0
	   max_order = 0
	   do while (line(1:4).ne.' ---')
	      n_terms(hrs) = n_terms(hrs) + 1
	      if (n_terms(hrs).gt.max_elements)
     >	      stop 'WCRECON: too many COSY terms!'
	      read (line,1200) (coeff(hrs,i,n_terms(hrs)),i=1,4),
     >			       (expon(hrs,j,n_terms(hrs)),j=1,5)
	      read (chan,1001) line
	      max_order = max(max_order, expon(hrs,1,n_terms(hrs)) +
     >				expon(hrs,2,n_terms(hrs)) +
     >				expon(hrs,3,n_terms(hrs)) +
     >				expon(hrs,4,n_terms(hrs)) +
     >				expon(hrs,5,n_terms(hrs)))
	   enddo
	   close (unit=chan)
	   firsttime = .false.
	endif

! Reset COSY sums.

	do i = 1,4
	   sum(i) = 0.
	enddo

! Convert hut quantities to right-handed coordinates, in meters and "radians".
! Make sure hut(5) is non-zero, to avoid taking 0.0**0 (which crashes)
	hut(1) = xs/100.d0		!cm --> m
	hut(2) = dxdzs			!slope ("radians")
	hut(3) = ys/100.d0		!cm --> m
	hut(4) = dydzs			!slope ("radians")
	hut(5) = fry/100.d0		!vert. position at target(cm-->m)
	if (abs(hut(5)).le.1.d-30) hut(5)=1.d-30

! Compute COSY sums.

	do i = 1,n_terms(hrs)
	   term = hut(1)**expon(hrs,1,i)*hut(2)**expon(hrs,2,i)
     >		* hut(3)**expon(hrs,3,i)*hut(4)**expon(hrs,4,i)
     >		* hut(5)**expon(hrs,5,i)
	   sum(1) = sum(1) + term*coeff(hrs,1,i)
	   sum(2) = sum(2) + term*coeff(hrs,2,i)
	   sum(3) = sum(3) + term*coeff(hrs,3,i)
	   sum(4) = sum(4) + term*coeff(hrs,4,i)
	enddo
     
! Load output values.

	delta_phi = sum(1)		! slope ("radians")
	y_tgt	  = sum(2)*100.d0	! m --> cm
	delta_t   = sum(3)		! slope ("radians")
	delta_p   = sum(4)*100.d0	! percent deviation
	return

! ============================ Format Statements ===============================

1001	format(a)
1200	format(1x,4g16.9,1x,5i1)

      END


















