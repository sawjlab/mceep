	logical*4 function locforunt(io_unit)

!------------------------------------------------------------------------------
!
! LOCFORUNT - Locate a free fortran I/O unit in the range
! [min_unit$,max_unit$]. If found, the free unit's number is returned
! in IO_UNIT, and locforunt returns .TRUE. Otherwise, IO_UNIT is
! unchanged, and we return .FALSE.
!
! DHP - May 1985
!
! DHP - March 1992, Modified for UNIX F77 compatibility. Now a LOGICAL*4 funct.
!----------------------------------------------------------------------------
	implicit integer*4 (a-z)

	parameter (min_unit = 10)
	parameter (max_unit = 29)
	logical*4 open
	save

!============================= Executable Code ================================

	locforunt = .false.			!assume failure
	do unit = min_unit,max_unit		!loop over range
	   inquire (unit=unit,opened=open)	!is it free?
	   if (.not.open) then			!yes...
	      io_unit = unit
	      locforunt = .true.
	      return
	   endif
	enddo
	return
	end
