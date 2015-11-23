       subroutine read_cosy_input_file(infile,mom,datdir,ndatdir)
*                                                                    *
**********************************************************************
c   input: infile from INPUT.F read
c
c   output: no direct output back but
c      reads file of spin matrix from COSY
c      which is a 3x3 tensor
c       gets:  
c              spimat(nline,3,3)
c              ip_spi(nline,4,3)
c
c
       IMPLICIT none
c
       character*300 infile
       character*100 datdir,tmpfile
       character*15 interfile(13) ! 13 == number of spin matrix files.
       
c
       integer icol1,icol2,inum
       integer il1,iskip
       integer nfile ! 13 == number of spin matrix files.
       integer ifind
       integer ndatdir
c
       double precision mom ! central momentum in MeV/c
       double precision momfile(13) ! 13 == number of spin matrix files.
       double precision sp_tmp11(3),sp_tmp12(3)
       double precision slope,momdiff,yint
c     
       INCLUDE 'spin.cmn'
c
       parameter (nfile=13)
       data momfile / 500.,750.,1000.,1250.
     >,1500.,1750.,2000.,2250.,2500.,2750.
     >,3000.,3250.,3500. /
       data interfile / 'stm_pc0500.spin'
     >,'stm_pc0750.spin'
     >,'stm_pc1000.spin'
     >,'stm_pc1250.spin'
     >,'stm_pc1500.spin'
     >,'stm_pc1750.spin'
     >,'stm_pc2000.spin'
     >,'stm_pc2250.spin'
     >,'stm_pc2500.spin'
     >,'stm_pc2750.spin'
     >,'stm_pc3000.spin'
     >,'stm_pc3250.spin'
     >,'stm_pc3500.spin' /
c
       if (infile .eq. 'internal') then
          write(*,*) ' Determining STM internally'
          do il1=2,nfile
             if (mom .ge. momfile(il1-1) .and.
     >          mom .lt. momfile(il1) ) ifind=il1
          enddo
          tmpfile=datdir(1:ndatdir)//'/'//interfile(ifind-1)
          write(*,*) ' opening file ', tmpfile
       open (unit=11,file=tmpfile,form='formatted',status='old')
          tmpfile=datdir(1:ndatdir)//'/'//interfile(ifind)
          write(*,*) ' opening file ', tmpfile
       open (unit=12,file=tmpfile,form='formatted',status='old')
c
       inum = 1
       il1 = 1 
       momdiff= momfile(ifind)-momfile(ifind-1)
 11     read(11,112,end=88,err=666)
     s        (sp_tmp11(icol1),icol1=1,3),
     s        (ip_spi(il1,icol2,inum),icol2=1,4), iskip,
     s        ip_spi(il1,5,inum)
 666   read(12,112,end=88,err=77)
     s        (sp_tmp12(icol1),icol1=1,3),
     s        (ip_spi(il1,icol2,inum),icol2=1,4), iskip,
     s        ip_spi(il1,5,inum)
       if (sp_tmp11(1).eq.0 .and.sp_tmp11(2).eq.0
     s   .and. sp_tmp11(3).eq.0) goto 77   
       do icol1=1,3
         slope= (sp_tmp12(icol1)-sp_tmp11(icol1))/momdiff
         yint= sp_tmp12(icol1) - slope*momfile(ifind)
         sp_ten(il1,icol1,inum)=slope*mom+yint
       enddo 
cxxx        write(3,112) 
cxxx     s        (sp_ten(il1,icol1,inum),icol1=1,3),
cxxx     s        (ip_spi(il1,icol2,inum),icol2=1,4), iskip,
cxxx     s        ip_spi(il1,5,inum)
       il1 = il1 + 1
       if (il1 .gt. 500) then
          write(*,*) ' Number of lines in COSY file:',
     s  infile,' is larger than 500'
          write(*,*) ' This is larger than array spimat'
          write(*,*) ' Need to change mceep'
          stop
       endif
       goto 11
 77    if (inum.eq.3) goto 88
       iligne_spi(inum)=il1-1
       il1=1
       inum = inum+1
       goto 11
 88    close(11)
       close(12)
c
c
       else ! read spin matrix directly from file
c
       open (unit=11,file=infile,form='formatted',
     s        status='old')
c
       inum = 1
       il1 = 1 
 1     read(11,112,end=888,err=777)
     s        (sp_ten(il1,icol1,inum),icol1=1,3),
     s        (ip_spi(il1,icol2,inum),icol2=1,4), iskip,
     s        ip_spi(il1,5,inum)
       if (sp_ten(il1,1,inum).eq.0 .and. sp_ten(il1,2,inum).eq.0
     s   .and. sp_ten(il1,3,inum).eq.0) goto 777   
cxxx       write(*,212) il1,
cxxx     s        (sp_ten(il1,icol1,inum),icol1=1,3),
cxxx     s        (ip_spi(il1,icol2,inum),icol2=1,4), iskip,
cxxx     s        ip_spi(il1,5,inum)
 112   format(1x,3g20.10,1x,6i1)
 212   format(i3,1x,3g20.10,1x,6i1)
       il1 = il1 + 1
       if (il1 .gt. 500) then
          write(*,*) ' Number of lines in COSY file:',
     s  infile,' is larger than 500'
          write(*,*) ' This is larger than array spimat'
          write(*,*) ' Need to change mceep'
          stop
       endif
       goto 1
 777   if (inum.eq.3) goto 888
       iligne_spi(inum)=il1-1
       il1=1
       inum = inum+1
       goto 1
 888   close(11)
c
      endif
c
        return
        end
