c
c ---------------------------------------------------------------------
c
c    Subroutine get_include_cwn
c
c    Purpose:   Read in hbook column-wise Ntuple and create include
c               file.  This is a preliminary step to adding r-function
c               values (done also by mceep_util).
c
c    Author:  P.E. Ulmer 
c
c    Date:    26-SEP-2002
c
c ---------------------------------------------------------------------
c
      subroutine get_include_cwn
c
      integer nwpawc
      integer ipaw
      parameter (nwpawc=2000000) 
      common /pawc/ ipaw(nwpawc)
c
      integer  ntuid,nchar
      integer  irc,lrec
c
      character*60    filepre
      character*80    infile,tmpfile
c
c ---------------------------------------------------------------------
c    Get input file name, Ntuple ID, etc.
c ---------------------------------------------------------------------
c
      write(6,*) 'Enter input hbook file prefix > '
      read(5,*) filepre
      write(6,*) 'Enter input Ntuple ID > '
      read(5,*) ntuid
      tmpfile = filepre//'.hbook'
      call squeeze(tmpfile,tmpfile,nchar)
      infile = tmpfile(1:nchar)
c      
c ---------------------------------------------------------------------
c    Open input hbook file.
c ---------------------------------------------------------------------
c
      call hlimit(nwpawc)
      lrec = 4096
      call hropen(10,'HBKINPUT',infile,'p',lrec,irc)
      if (irc.ne.0) goto 888
      call hcdir('//HBKINPUT',' ')
      call hrin (ntuid, 9999, 0) 
c   
c ---------------------------------------------------------------------
c    Write the include (actually a skeleton subroutine).
c    The procedure build_cwn should be executed after this code is run
c    in order to make the actual include file from this skeleton.
c ---------------------------------------------------------------------
c
      open(unit=11,name='skeleton',status='unknown',form='formatted')
      call huwfun(11,ntuid,'skeleton',0,'B')
      close(unit=11)
      write(6,*) 'NOW Execute ~/mceep/utilities/build_cwn_include'
c
c ---------------------------------------------------------------------
c     Close the input hbook file.
c ---------------------------------------------------------------------
c
      call hcdir('//HBKINPUT',' ')
      call hrend('HBKINPUT')
c
      goto 999
c
 888  write(*,*) ' Problem input file irc = ',irc
      return
 999  write(*,*) ' end'
      end
