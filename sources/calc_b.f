       subroutine calc_b(thetae,model)
 
C------------------------------------------------------------------------------
c subroutine: calc_b
c     author: L.Todor
c       date: July 1999
c    purpose: calculate b (or 1/x0) parameter used in external bremsstrahlung 
c            distribution expression; x0=radiation length is 
c             the mean distance over which a high energy
c             electron losses all but 1/e of its energy 
c             by bresstrahlung
c  reference: Physical Review D, Particle and Fields, July 1996 
c         ! b=1/xo^2 unit is g/cm^2
C------------------------------------------------------------------------------

       implicit none
       COMMON /extrad_wall/ btw 
       double precision alpha,al,aconst
       double precision fz,lr,lrprim,btw(3),alw(3),rl(3),thetae
       integer model
C
       INCLUDE 'eloss.cmn'
C
C     order of materials below: mylar, air, Al
       data rl/28.7D0,30420.D0,8.9D0/

       aconst=1/716.408D0
       alpha=1.D0/137.036D0
       al=alpha*dble(ztg)
       al=al*al
       fz=1/(1+al)
       fz=fz+0.20206D0-0.0369D0*al
       fz=fz+0.0083D0*al*al-0.002D0*al*al*al
       fz=al*fz
       if(ztg.eq.1) then
            lr=5.31D0
            lrprim=6.144D0
       else if(ztg.eq.2) then
            lr=4.79D0
            lrprim=5.621D0
       else if(ztg.eq.3) then
            lr=4.74D0
            lrprim=5.805D0
       else if(ztg.eq.4) then
            lr=4.71D0
            lrprim=5.924D0
       else
            lr=log(184.15D0)-(log(dble(ztg)))/3.0D0
            lrprim=log(1194.0D0)-(log(dble(ztg)))/1.5D0
       end if
       b=aconst*(dble(ztg)*dble(ztg)*(lr-fz)+dble(ztg)*lrprim)/atg
c       write(6,*) 'b=',b
c approximative form
c       b=aconst*dble(ztg)*(dble(ztg)+1)*log(287.D0/sqrt(dble(ztg)))/atg
c       write(6,*) 'b=',b
c      setup call; calculate sum of bt for incoming and outgoing
c      through wall or endcap - except the effective path through target
c first Al wall for incoming beam
       alw(1)=0.007112D0 !cm
       btw(1)=alw(1)/rl(3)
   
cpeu - version 3.9
       if(model .eq. 1) btw(1) = 0.d0
cpeu - version 3.9

c outgoing  through lateral wall
       alw(2)=0.02032D0   ! cm
       if(thetae.ne.0.D0) then          
             btw(2)=alw(2)/(sin(thetae)*rl(3))+0.01778D0/rl(1)+
     &              58.D0/rl(2)+0.04064D0/rl(3)
       end if 

cpeu - version 3.9
       if(model .eq. 1)  btw(2)=0.01778D0/rl(1)+
     &              58.D0/rl(2)+0.04064D0/rl(3)
cpeu - version 3.9

c outgoing  through end cap
       alw(3)=0.011400D0  ! cm
       if(thetae.ne.(2*atan(1.D0))) then          
             btw(3)=alw(3)/(cos(thetae)*rl(3))+0.01778D0/rl(1)+
     &              58.D0/rl(2)+0.04064D0/rl(3)
       else 
          btw(3)=btw(2)
       end if 
       if(thetae.eq.0.D0) btw(2)=btw(3) 
            
cpeu - version 3.9
       if(model .eq. 1) btw(3)=btw(2)
cpeu - version 3.9

       return 
       end

