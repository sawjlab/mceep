c
c ---------------------------------------------------------------------
c
c    Subroutine make_new_cwn
c
c    Purpose:   Read in HBOOK column-wise Ntuple and create new
c               Ntuple which is a copy of input one, plus the
c               values of the r-function ala J. LeRose.
c
c    Author:  P.E. Ulmer 
c             Based on routines of djh and mkj
c
c    Date:    27-SEP-2002
c
c    Method: 
c
c          See ~/mceep/utilities/README_CWN file before executing.
c
c          The routine get_include_cwn (within mceep_util)
c          must be executed FIRST in order to generate an include
c          containing the Ntuple structure.  Then the script
c          ~/mceep/utilities/build_cwn_include must be executed.
c          At this point, this routine can be run.
c
c          Reads column-wise HBOOK Ntuple and writes new Ntuple
c          which is a copy of the first plus the r-function
c          value(s).  Both electron and hadron r-function
c          values are written to the end of the given Ntuple.
c          If the Ntuple contains the coordinates of only one
c          spectrometer, then the other r-function is just
c          set to zero, but still written.
c          
c ---------------------------------------------------------------------
c
      subroutine make_new_cwn
c
      implicit none
c
      integer nwpawc
      integer ipaw
      parameter (nwpawc=2000000) 
      common /pawc/ ipaw(nwpawc)
c
      integer  ntuid
      integer  nchar,nevents,iev
      integer  irc,ierror,icycle,lrec
c
      real     rfunction,rfn_e,rfn_h
      real     y_l,dp_l,th_l,ph_l
      real     y_r,dp_r,th_r,ph_r
      logical  left,right,coll_e,coll_h
c
      character*1     ans
      character*60    filepre
      character*80    infile,outfile,tmpfile
      character*1000  block1
c
c ---------------------------------------------------------------------
c     Ntuple Variable Declarations
c ---------------------------------------------------------------------
c
      common /RFN/    rfn_e,rfn_h
c
      INCLUDE 'ntu.cmn'
c
      irc = 0
c
c ---------------------------------------------------------------------
c    Get input file name, Ntuple ID, etc.
c ---------------------------------------------------------------------
c
      write(6,*) '---------------------------------------'
      write(6,*) '  YOU MUST FIRST GENERATE INCLUDE FILE '
      write(6,*) '      using mceep_util option 5        '
      write(6,*) '   See ~/mceep/utilities/README_CWN    '
      write(6,*) '---------------------------------------'
      write(6,*) 'Enter input hbook file prefix > '
      read(5,*) filepre
      write(6,*) 'Enter input Ntuple ID > '
      read(5,*) ntuid
      tmpfile = filepre//'.hbook'
      call squeeze(tmpfile,tmpfile,nchar)
      infile = tmpfile(1:nchar)
      tmpfile = filepre//'_rfn.hbook'
      call squeeze(tmpfile,tmpfile,nchar)
      outfile = tmpfile(1:nchar)
c
c ---------------------------------------------------------------------
c    Determine which arms to compute the r-functions for.
c ---------------------------------------------------------------------
c
      left  = .false.        ! initialize
      right = .false.
c
      write(6,100) 
 100  format(/,' Left arm, Right arm or Both? (L,R,B) > ')
      read(5,*) ans
      if(ans.eq.'L' .or. ans.eq.'l') left  = .true.
      if(ans.eq.'R' .or. ans.eq.'r') right = .true.
      if(ans.eq.'B' .or. ans.eq.'b') then
         left  = .true.
         right = .true.
      endif
c
c ---------------------------------------------------------------------
c    Determine whether each arm has the 6 msr collimator or is "OPEN".
c ---------------------------------------------------------------------
c
      coll_e = .false.       ! initialize
      coll_h = .false.
c
      if(left) then
         write(6,110)
 110     format(/,' Left arm collimator? (Y,N) > ')
         read(5,*) ans
         if(ans.eq.'Y' .or. ans.eq.'y') coll_e = .true.
      endif
      if(right) then
         write(6,120)
 120     format(/,' Right arm collimator? (Y,N) > ')
         read(5,*) ans
         if(ans.eq.'Y' .or. ans.eq.'y') coll_h = .true.
      endif
c
      call hlimit(nwpawc)
c      
c ---------------------------------------------------------------------
c    Open input hbook file and determine number of events.
c ---------------------------------------------------------------------
c
      lrec = 4096
      call hropen(10,'HBKINPUT',infile,'p',lrec,irc)
      if (irc.ne.0) goto 888
      call hcdir('//HBKINPUT',' ')
      call hrin (ntuid, 9999, 0) 
c
      call hnoent(ntuid,nevents)  
      write(6,*) ' Number of events: ',nevents
c
      call hbname(ntuid,' ',0,'$CLEAR')
      call hbname(ntuid,'BLOCK1',FIRST_VAR,'$SET')
c
c ---------------------------------------------------------------------
c    Open the output hbook file and book the new Ntuple.
c ---------------------------------------------------------------------
c
      call hropen(20,'HBKOUTPUT',outfile,'NQ',lrec,irc)
      call hcdir('//HBKOUTPUT',' ')
      if (irc.ne.0) goto 889
      CALL hbnt(100,'ntu_out',' ')
      CALL hbname(100,'BLOCK1',FIRST_VAR,BLOCK1)
      CALL hbname(100,'RFN',RFN_E,'RFN_E,RFN_H')
c
c ---------------------------------------------------------------------
c    Call LeRose r-function initializers.
c ---------------------------------------------------------------------
c
      if(left)  call init_r_function(1,coll_e)
      if(right) call init_r_function(2,coll_h)
c
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c    Start the event loop.
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c
      do iev=1,nevents
c
c
c ---------------------------------------------------------------------
c    Retrieve the input variable values.
c ---------------------------------------------------------------------
c
         call hcdir ('//HBKINPUT', ' ')
         call hgnt(ntuid,iev,ierror)
c
c ---------------------------------------------------------------------
c    Call LeRose r-function routine and fill additional slots of
c    output Ntuple with values of r-function.
c    Syntax:     rfunction(iarm,ytg,delta,thtg,phtg)
c      with:     ytg=meters,delta=fraction,thtg=tan(rad),phtg=tan(rad)
c ---------------------------------------------------------------------
c
         INCLUDE 'ntu_user.cmn'
c       
         if(left) then                 ! electron
            rfn_e = rfunction(1,y_l,dp_l,th_l,ph_l)
         else
            rfn_e = 0.
         endif
c
         if(right) then                 ! hadron
            rfn_h = rfunction(2,y_r,dp_r,th_r,ph_r)
         else
            rfn_h = 0.
         endif
c
c ---------------------------------------------------------------------
c    Fill the output Ntuple.
c ---------------------------------------------------------------------
c
         call hcdir('//HBKOUTPUT',' ')
         call hfnt(100)
c
      enddo
c
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c    End of the event loop.
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c
c ---------------------------------------------------------------------
c     Close the output and input hbook files.
c ---------------------------------------------------------------------
c
      call hcdir('//HBKOUTPUT',' ')
      call hrout(100,icycle,'HBKOUTPUT')
      call hrend('HBKOUTPUT') 
c
      call hcdir('//HBKINPUT',' ')
      call hrend('HBKINPUT')
c
      goto 999
c
 888  write(*,*) ' Problem input file irc = ',irc
      return
 889  write(*,*) ' Problem output file irc = ',irc
      return
 999  write(*,*) ' end'
      end
