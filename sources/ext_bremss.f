       subroutine ext_bremss(path,e0,erad,whichwall)
 
C------------------------------------------------------------------------------
c subroutine: ext_bremss
c     author: L.Todor
c       date: July 1999
c    purpose: generate an energy loss by bremsstrahlung
c             according to the known distribution of this quantity           
c  reference: Y.S.Tsai "Pair Production and Bremsstrahlung of Charged Leptons"
c             Review of Modern Physics 46, 1974 (815) & errata!)
c       used: similar method as in aeexb
c 
C------------------------------------------------------------------------------

       implicit none

       COMMON /extrad_wall/ btw 
       integer whichwall 
       double precision detrial,epsilon,iext,bt,gammabt,e0,erad
       double precision var1,thick,path,r1low
       double precision dgamma,r1,r2,btw(3)
       real rndm(2)
C
       INCLUDE 'eloss.cmn'
C
       thick=path*100.D0*rho 
       bt=b*thick
c       write(6,*) 'e0=',e0
       if (whichwall.eq.1) then
          bt=bt+btw(1)
       else if (whichwall.eq.2) then
          bt=bt+btw(2)
       else if (whichwall.eq.3) then
          bt=bt+btw(3)
       end if 
c       write(6,*) 'bt=', bt
       gammabt=dgamma(bt+1.D0)
c       write(6,*) 'gammabt',gammabt
 17    call RANECU(rndm,2) 
c       write(6,*) 'Randoms ',rndm(1),rndm(2)
       r1=dble(rndm(1))
       r2=dble(rndm(2))
       r1low=(0.001d0/e0)**(bt)
       if (r1.lt.r1low) then
             erad=0.D0
             return
       end if
       detrial=e0*(r1)**(1.D0/bt)
c       write(6,*) 'detrial=', detrial
       var1=detrial/e0
       iext=bt/(e0*gammabt)*(var1**(bt-1.D0))
     &         *(1-var1+0.75D0*var1**2)
       epsilon=r2*(1.25D0/e0)*(bt/gammabt)*var1**(bt-1.D0)
C------------------------------------------------------------------------------
c reject accept technique
C------------------------------------------------------------------------------
       if(epsilon.gt.iext) go to 17
       erad=detrial 
       return 
       end

        function dgamma(xx)
        double precision cof(6),stp,half,one,fpf,x,tmp,ser,xx
        double precision gammln,dgamma
        data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     +   -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
        data half,one,fpf/0.5d0,1.0d0,5.5d0/
        x=xx-one
        tmp=x+fpf
        tmp=(x+half)*log(tmp)-tmp
        ser=one
        do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11      continue
        gammln=tmp+log(stp*ser)
        dgamma=dexp(gammln)
        return
        end

