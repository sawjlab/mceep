c
c ---------------------------------------------------------------------
c
c    Subroutine make_new_rwn
c    Purpose:   Read in hbook row-wise Ntuple and create new Ntuple
c               which is a copy of input one, plus the values
c               of the r-function ala J. LeRose.
c
c    Author:  M. Jones
c
c    Modifications:  
c          P. Ulmer 05-JUL-2001
c          Now reads arbitrary Ntuple structure and searches
c          within for set of variables needed by r-function.
c          The default is to write both electron and hadron
c          R-function values to the end of the given Ntuple.
c          If only one or the other set of variables is found 
c          in the Ntuple, then the other r-function is just
c          set to zero, but still written (it's easier on my
c          weary programming fingers).
c          
c ---------------------------------------------------------------------
c
      subroutine make_new_rwn
c
      implicit none
c
      integer nwpawc,iquest
      real blanc
      parameter (nwpawc=1000000) 
      common /pawc/  blanc(nwpawc)
      common /quest/ iquest(100)
c
      integer  irc,i,ierror,icycle,k,idn,imatch_e,imatch_h
      integer  nchar,nevents,iev,iconfig
      integer  id_e_th,id_e_y,id_e_ph,id_e_dp
      integer  id_h_th,id_h_y,id_h_ph,id_h_dp
      logical  coll_e,coll_h,match
      real     rfunction
c
      integer      ntuid,nvar
      real         data_in(100),data_out(100)
      real         rlow(100),rhigh(100)
      character*8  chtag(100),chtag_out(100)
      character*80 chtitl
c
      character*60 filepre
      character*80 infile,outfile,tmpfile
c
cpeu      call getenv('old_hbook', infile)
cpeu      call getenv('new_hbook', outfile)
c
c
c ---------------------------------------------------------------------
c    Get input file name, Ntuple ID, etc.
c ---------------------------------------------------------------------
c
      write(6,*) 'Enter input hbook file prefix > '
      read(5,*) filepre
      write(6,*) 'Enter input Ntuple ID > '
      read(5,*) ntuid
 10   write(6,100) 
 100  format(/,' L-arm coll.    & R-arm coll. ----------- 1 '
     &      ,/,' L-arm no coll. & R-arm coll. ----------- 2 '
     &      ,/,' L-arm coll.    & R-arm no coll. -------- 3 '
     &      ,/,' L-arm no coll. & R-arm no coll. -------- 4 '
     &     ,//,' Enter Option > ')
      read(5,*) iconfig
      if(iconfig .eq. 1) then 
         coll_e = .true.
         coll_h = .true.
      elseif(iconfig .eq. 2) then
         coll_e = .false.
         coll_h = .true.
      elseif(iconfig .eq. 3) then
         coll_e = .true.
         coll_h = .false.
      elseif(iconfig .eq. 4) then
         coll_e = .false.
         coll_h = .false.
      else
         goto 10
      endif
c 
      tmpfile = filepre//'.hbook'
      call squeeze(tmpfile,tmpfile,nchar)
      infile = tmpfile(1:nchar)
      tmpfile = filepre//'_rfn.hbook'
      call squeeze(tmpfile,tmpfile,nchar)
      outfile = tmpfile(1:nchar)
c
      call hlimit(nwpawc)
      iquest(10) = 65000
c      
c ---------------------------------------------------------------------
c    Open input hbook file and determine parameters of selected 
c    Ntuple (i.e. variable names, etc.).
c ---------------------------------------------------------------------
c
      call hropen(10,'HBKINPUT',infile,' ',1024,irc)
      if (irc.ne.0) goto 888
      call hcdir('//HBKINPUT',' ')
      call hgnpar(ntuid,'make_new_ntuple')
      call hnoent(ntuid,nevents)  
      call hgiven(ntuid,chtitl,nvar,chtag,rlow,rhigh)
      call hgiven(ntuid,chtitl,nvar,chtag,rlow,rhigh)  ! seem to need
                                                       ! another call
      write(6,*) ' Number of events: ',nevents
c
c ---------------------------------------------------------------------
c    Determine the indices for the variables needed to compute the
c    r-function.  Also, define the variable names for the output
c    Ntuple.
c ---------------------------------------------------------------------
c
      imatch_e = 0
      imatch_h = 0
      do i=1,nvar
         if(match(chtag(i),'E.TH_TG ',8)) then
            id_e_th = i 
            imatch_e = imatch_e + 1
         elseif(match(chtag(i),'E.Y_TG  ',8)) then
            id_e_y  = i
            imatch_e = imatch_e + 1
         elseif(match(chtag(i),'E.PH_TG ',8)) then
            id_e_ph = i
            imatch_e = imatch_e + 1
         elseif(match(chtag(i),'E.DP    ',8)) then
            id_e_dp = i
            imatch_e = imatch_e + 1
         endif
c                                       
         if(match(chtag(i),'H.TH_TG ',8)) then
            id_h_th = i
            imatch_h = imatch_h + 1
         elseif(match(chtag(i),'H.Y_TG  ',8)) then
            id_h_y  = i
            imatch_h = imatch_h + 1
         elseif(match(chtag(i),'H.PH_TG ',8)) then
            id_h_ph = i
            imatch_h = imatch_h + 1
         elseif(match(chtag(i),'H.DP    ',8)) then
            id_h_dp = i
            imatch_h = imatch_h + 1
         endif
c
         chtag_out(i) = chtag(i)
      enddo
c
      if(imatch_e .eq. 4) then
         write(6,*) ' Electron variables found in Ntuple '
      else
         write(6,*) ' Electron variables NOT   in Ntuple '
      endif
c
      if(imatch_h .eq. 4) then
         write(6,*) ' Hadron   variables found in Ntuple '
      else
         write(6,*) ' Hadron   variables NOT   in Ntuple '
      endif
c
      chtag_out(nvar+1) = 'RFN_E   '
      chtag_out(nvar+2) = 'RFN_H   '
c
c ---------------------------------------------------------------------
c    Open the output hbook file and book the new Ntuple.
c ---------------------------------------------------------------------
c
      call hropen(20,'HBKOUTPUT',outfile,'NQ',4096,irc)
      if (irc.ne.0) goto 889
      call hcdir('//HBKOUTPUT',' ')
      call hbookn(100,'ntu_out',nvar+2,'HBKOUTPUT',1024,chtag_out)
c
c ---------------------------------------------------------------------
c    Call LeRose r-function initializers.
c ---------------------------------------------------------------------
c
      if(imatch_e .eq. 4) call init_r_function(1,coll_e)
      if(imatch_h .eq. 4) call init_r_function(2,coll_h)
c
c ---------------------------------------------------------------------
c    Loop through the data and copy existing Ntuple into new output
c    Ntuple.
c ---------------------------------------------------------------------
c
      do iev=1,nevents
         call hcdir('//HBKINPUT',' ')
         call hgnf(ntuid,iev,data_in,ierror)
         do k=1,nvar
            data_out(k)=data_in(k)
         enddo
c
c ---------------------------------------------------------------------
c    Call LeRose r-function routine and fill last slots of output
c    Ntuple with values of r-function.
c    Syntax:     rfunction(iarm,ytg,delta,thtg,phtg)
c      with:     ytg=meters,delta=fraction,thtg=tan(rad),phtg=tan(rad)
c ---------------------------------------------------------------------
c
         if(imatch_e .eq. 4) then                 ! electron
            data_out(nvar+1)=rfunction(1,data_in(id_e_y),
     &                                   data_in(id_e_dp),
     &                                   data_in(id_e_th),
     &                                   data_in(id_e_ph))
         else
            data_out(nvar+1)=0.
         endif
c
         if(imatch_h .eq. 4) then                 ! hadron
            data_out(nvar+2)=rfunction(2,data_in(id_h_y),
     &                                   data_in(id_h_dp),
     &                                   data_in(id_h_th),
     &                                   data_in(id_h_ph))
         else
            data_out(nvar+2)=0.
         endif
c
         call hcdir('//HBKOUTPUT',' ')
         call hfn(100,data_out)
c         if(float(iev/100000)-float(iev)/100000. .eq. 0.) then
c            call hrout(100,icycle,'HBKOUTPUT')
c            call rzpurg(1)
c            write(6,*) ' calling hrout '
c         endif
      enddo
c
c ---------------------------------------------------------------------
c     Close the output hbook file.
c ---------------------------------------------------------------------
c
      call hcdir('//HBKOUTPUT',' ')
      call hrout(100,icycle,'HBKOUTPUT')
      call hrendc('HBKOUTPUT') 
      goto 999
c
 888  write(*,*) ' Problem input file irc = ',irc
      return
 889  write(*,*) ' Problem output file irc = ',irc
      return
 999  write(*,*) ' end'
      end
c
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c
c ---------------------------------------------------------------------
c     Logical Function MATCH
c ---------------------------------------------------------------------
c
      LOGICAL FUNCTION MATCH(S1,S2,LEN)
      CHARACTER*1 S1(LEN),S2(LEN)
      MATCH = .TRUE.
      DO I=1,LEN
        IF(S1(I).NE.S2(I))MATCH = .FALSE.
      ENDDO
      RETURN
      END

