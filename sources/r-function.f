      subroutine init_r_function(iarm,collimator)
c      
c initializes arrays etc in the common block of the r-function package
c contains some diagnostic write statements you can comment them out if 
c you want.
c
c You must call init_r_function once before calling rfunction 
c (which you can then do over and over again without reinitializing).
c                                           -JJL 6/6/01
      parameter (narray=1400)
      common /lines/rline(2,12,narray,6),tmin(2,narray),tmax(2,narray),
     +              phimin(2,narray),phimax(2,narray),
     +              dy(2),dd(2),npts(2)
c
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILELOC
      INTEGER        N_DAT_DIR,IARM,NLOC      !NLOC is #characters in filepath
      CHARACTER*1    COLLIMATOR
c
      integer ndum
      character*24 daet
c
      IF(IARM.EQ.1 .AND. COLLIMATOR.EQ.'T') THEN
         FILELOC = DAT_DIR(1:N_DAT_DIR)//
     +        '/r-function/leftarm/collimator'
         NLOC = N_DAT_DIR + 30
      ELSEIF(IARM.EQ.1 .AND. COLLIMATOR.EQ.'F') THEN
         FILELOC = DAT_DIR(1:N_DAT_DIR)//
     +        '/r-function/leftarm/no-collimator'
         NLOC = N_DAT_DIR + 33
      ELSEIF(IARM.EQ.1 .AND. COLLIMATOR.EQ.'S') THEN
         FILELOC = DAT_DIR(1:N_DAT_DIR)//
     +        '/r-function/leftarm/septum'
         NLOC = N_DAT_DIR + 26
      ELSEIF(IARM.EQ.2 .AND. COLLIMATOR.EQ.'T') THEN
         FILELOC = DAT_DIR(1:N_DAT_DIR)//
     +        '/r-function/rightarm/collimator'
         NLOC = N_DAT_DIR + 31
      ELSEIF(IARM.EQ.2 .AND. COLLIMATOR.EQ.'F') THEN
         FILELOC = DAT_DIR(1:N_DAT_DIR)//
     +        '/r-function/rightarm/no-collimator'
         NLOC = N_DAT_DIR + 34
      ELSEIF(IARM.EQ.2 .AND. COLLIMATOR.EQ.'S') THEN
         FILELOC = DAT_DIR(1:N_DAT_DIR)//
     +        '/r-function/rightarm/septum'
         NLOC = N_DAT_DIR + 27
      ENDIF
c
c  for each of the 12 lines defining the theta, phi acceptance of the form:
c   0=A*theta+B*phi+C
c    rline(iarm,i,n,1)= yn (y value of nth set) i refers to the line number
c    rline(iarm,i,n,2)= deltan (delta value for nth set)
c    rline(iarm,i,n,3)=A (for the nth set)
c    rline(iarm,i,n,4)=B
c    rline(iarm,i,n,5)=C
c    rline(iarm,i,n,6)=R^2 or 2 for vertical lines 
c                          or 3 for lines through the origin
c
c  read in the files with the line parameters
c
      open(11,file=fileloc(1:nloc)//'/line1.dat',status='old')
      open(12,file=fileloc(1:nloc)//'/line2.dat',status='old')
      open(13,file=fileloc(1:nloc)//'/line3.dat',status='old')
      open(14,file=fileloc(1:nloc)//'/line4.dat',status='old')
      open(15,file=fileloc(1:nloc)//'/line5.dat',status='old')
      open(16,file=fileloc(1:nloc)//'/line6.dat',status='old')
      open(17,file=fileloc(1:nloc)//'/line7.dat',status='old')
      open(18,file=fileloc(1:nloc)//'/line8.dat',status='old')
      open(19,file=fileloc(1:nloc)//'/line9.dat',status='old')
      open(20,file=fileloc(1:nloc)//'/line10.dat',status='old')
      open(21,file=fileloc(1:nloc)//'/line11.dat',status='old')
      open(22,file=fileloc(1:nloc)//'/line12.dat',status='old')
      open(23,file=fileloc(1:nloc)//'/dy-dd.dat',status='old')
      open(30,file=fileloc(1:nloc)//'/min-max.dat',status='old')
      read(23,*)dy(iarm),dd(iarm)   !y and delta grid spacing
      close(23)
      write(6,*)'initializing r-function arrays'
      open(31,file=fileloc(1:nloc)//'/information.dat',status='old')
      read(31,'(a24)')daet
      write(6,*)'Files generated on:'
      write(6,'(a24)')daet
      read(31,'(a24)')daet
      write(6,'(a24)')daet
      read(31,'(a24)')daet
      write(6,'(a24)')daet
      close(31)
      do 10 i=1,narray
      
      do 11 n=1,12
      nin=n+10
      read(nin,*,end=12)ndum,(rline(iarm,n,i,j),j=1,6)
   11 continue
      read(30,*)ydum,ddum,tmin(iarm,i),tmax(iarm,i),
     +          phimin(iarm,i),phimax(iarm,i)
   10 continue
   12 npts(iarm)=i
      if(npts(iarm).ge.narray)then
      write(6,*)'ERROR in init-r-function: Need to increase narray'
      endif
c
      do n=11,30
         close(unit=n)
      enddo
c
      return
      end
      
      subroutine getindex(iarm,y,d,index)
      
c for given y and delta (y,d) find index for line parameters
      parameter (narray=1400)                  
      common /lines/rline(2,12,narray,6),tmin(2,narray),tmax(2,narray),
     +              phimin(2,narray),phimax(2,narray),
     +              dy(2),dd(2),npts(2)
      real y,d,ddummy,testa
      integer index,iarm
      ddummy=100000.
      index=0
c      write(6,*)'npts=',npts(iarm)
      do i=1,npts(iarm)
       a=((y-rline(iarm,1,i,1))**2)+((d-rline(iarm,1,i,2))**2)
       if(a.lt.ddummy)then
        index=i
        ddummy=a
       endif
      enddo
c      write(6,*)index
c test to see that a is not too big
      ddummy=sqrt(ddummy)
      testa=sqrt(dd(iarm)**2+dy(iarm)**2)/2
c      write(6,*)dd(iarm),dy(iarm),testa,ddummy
      if(ddummy.gt.testa) index=0   !outside the grid
      return
      end
      
      function rfunction(iarm,y,d,theta,phi)
c
c  sample rfunction call
c   y is in meters (y0)
c   d (delta) is a fraction
c   theta and phi are radians (theta_0 and phi_0)
c      
c  sample rfunction call:  answer=rfunction(iarm,y,d,theta,phi)
c      
c  answer is the distance in phi/theta space of the given  
c  theta and phi point from the edge of the acceptance for the given y 
c  and d negative numbers are outside the acceptance (-1000. 
c  is the default outside value)
c      
      
      parameter (narray=1400)                  
      common /lines/rline(2,12,narray,6),tmin(2,narray),tmax(2,narray),
     +              phimin(2,narray),phimax(2,narray),
     +              dy(2),dd(2),npts(2)
      real y,d,theta,phi,dline(16)
      integer index,iarm
      rfunction=100000.
      call getindex(iarm,y,d,index)
c      write(6,*)'index=',index,y,d 
      if(index.gt.0)then
      do 10 i=1,12   !loop through the 12 lines
c
c set up the signs  (see notes: 5/25-30/01)
c
      if(rline(iarm,i,index,4).eq.0.) go to 10  !ignore vertical lines, 
                                                !the box takes care of it      
      if(rline(iarm,i,index,6).eq.3.) go to 10  !ignore lines through 
                                                !the origin
      em=-rline(iarm,i,index,3)/rline(iarm,i,index,4) 
                                     !conventional slope of the line
      isign=+1
      if((tmax(iarm,index).gt.0.)
     +                 .and.(tmin(iarm,index).lt.0.))then  
                                      !origin is inside the acceptance
       go to 7   !all signs positive
       else
       rfunction=-1000.
       return
      endif
c
c  dline is the distance from the point to the line
c       
  7   dline(i)=isign*(
     +         (rline(iarm,i,index,3)*theta)
     +        +(rline(iarm,i,index,4)*phi)
     +        + rline(iarm,i,index,5))
     +        /sqrt(rline(iarm,i,index,3)**2+rline(iarm,i,index,4)**2)
c
c  then find the minimum
c      
      if(dline(i).lt.rfunction)rfunction=dline(i)
  10  continue
c
c make the simple box 
c     
      dline(13)=tmax(iarm,index)-theta
      if(dline(13).lt.rfunction)rfunction=dline(13)
      dline(14)=theta-tmin(iarm,index)
      if(dline(14).lt.rfunction)rfunction=dline(14)
      dline(15)=phimax(iarm,index)-phi
      if(dline(15).lt.rfunction)rfunction=dline(15)
      dline(16)=phi-phimin(iarm,index)
      if(dline(16).lt.rfunction)rfunction=dline(16)

      else
      rfunction=-1000.  !default no acceptance value
      endif
      return
      end
      
      
