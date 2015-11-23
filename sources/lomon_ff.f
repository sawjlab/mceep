c
c Try to reproduce fits of Lomon in PRC64 035204
c
c N.B. !!!!! The factor N that is introduced in his paper, and included whenever
c            the dispersion relation fits to the rho are used, needs to be divided
c            by two compared to the values listed in the table in the paper.  I am
c            not sure if it is a typographical error in the equations or something
c            else.
c
c First set is labelled as GK(3) in his paper.  It includes the rho,omega,phi
c as traditional meson poles, plus the QCD constraint.
c
c Second set is labelled as DRN-GK(3) in his paper.  It includes the rho 
c from dispersion relations, and the omega and phi as traditional meson poles, plus
c the QCD constraint.
c 
c Third set is labelled at DRN-GK'(1) in his paper.  It includes the rho 
c from dispersion relations, the omega, phi, and rho' as meson pole terms, plus
c it includes the second form for the quark-nucleon form factors in the QCD 
c constraint.
c 
c Fourth set is his new fit, where he has included now also the omega'.  Other than
c this, it is the same form as DRN-GK'(1).  There was a typo in the numbers that he
c gave me for gphi/fphi for this fit, so I simply adjusted the value to give the best
c fit to the data. 
c
c
      subroutine lomonff(iset,q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
c
      implicit none

      integer*4 i,iset
      real*8 grhorat,kapparho,gomegarat,kappaomega,gomegaprat
      real*8 mrhop,kappaomegap
      real*8 gphirat,kappaphi,muphi,lambda1,lambdad,lambda2
      real*8 lambdaqcd,n,kappanu,kappas,mrho,momega,mphi,momegap
      real*8 q2,qt2
      real*8 f1rho,f2rho,f1d,f2d,f1iv,f2iv
      real*8 f1omega,f2omega,f1phi,f2phi
      real*8 f1is,f2is
      real*8 t1,t2,t3
      real*8 f1p,f2p,gep,gmp,tau,gd,mp,mu,ratio,gen,gmn,f1n,f2n,mun

c parameter fixed for all fits
      kappanu=3.706
      kappas=-0.12
      momega=0.784
      mphi=1.019
      mrho=0.776
      mrhop=1.45
      momegap=1.419
      mp=0.9382796
      mu=2.79
      mun=-1.91

c GK(3) in PRC
c
      if(iset .eq. 1) then
c
         grhorat=0.4466
         kapparho=4.3472
         gomegarat=0.4713
         kappaomega=21.762
         gphirat=-0.8461
         kappaphi=11.849
         muphi=1.1498
         lambda1=0.9006
         lambdad=1.7038
         lambda2=1.1336
         lambdaqcd=0.0312
c      
         tau=q2/(4.0*mp**2)
         qt2=q2*log((lambda2**2+q2)/lambdaqcd**2)/
     &        log(lambda2**2/lambdaqcd**2)
         f1rho=(lambda1**2/(lambda1**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))
         f1omega=f1rho
         f1phi=f1rho*(q2/(lambda1**2+q2))**1.5
         f1d=(lambdad**2/(lambdad**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))
         f2rho=(lambda1**2/(lambda1**2+qt2))**2*
     &        (lambda2**2/(lambda2**2+qt2))
         f2omega=f2rho
         f2phi=f2rho*((lambda1**2*(q2+muphi**2))/
     &        (muphi**2*(lambda1**2+q2)))**1.5
         f2d=(lambdad**2/(lambdad**2+qt2))**2*
     &        (lambda2**2/(lambda2**2+qt2))
         t2=grhorat*mrho**2/(mrho**2+q2)*f1rho
         t3=(1-grhorat)*f1d
         f1iv=t2+t3         
c         if(q2.eq.8.)write(*,*)q2,t2,t3
         t2=kapparho*grhorat*mrho**2/(mrho**2+q2)*f2rho
         t3=(kappanu-kapparho*grhorat)*f2d
         f2iv=t2+t3
c         if(q2.eq.8.)write(*,*)q2,t2,t3
         t1=gomegarat*(momega**2/(momega**2+q2))*f1omega
         t2=gphirat*(mphi**2/(mphi**2+q2))*f1phi
         t3=(1-gomegarat)*f1d
         f1is=t1+t2+t3
         t1=kappaomega*gomegarat*(momega**2/(momega**2+q2))*f2omega
         t2=kappaphi*gphirat*(mphi**2/(mphi**2+q2))*f2phi
         t3=(kappas-kappaomega*gomegarat-kappaphi*gphirat)*f2d
         f2is=t1+t2+t3
         f1p=(f1is+f1iv)/2.0
         f2p=(f2is+f2iv)/2.0
         f1n=(f1is-f1iv)/2.0
         f2n=(f2is-f2iv)/2.0
         gep=f1p-tau*f2p
         gmp=f1p+f2p
         gen=f1n-tau*f2n
         gmn=f1n+f2n
         gd=1.0/(1.0+q2/0.71)**2
         ratio=mu*gep/gmp

c
c set 2 : DRN-GK(3) in PRC 64 035204 
c

      elseif(iset .eq. 2) then
c
         grhorat=0.1013
         kapparho=-15.870
         gomegarat=0.6604
         kappaomega=8.847
         gphirat=-0.4054
         kappaphi=13.6415
         muphi=1.127
         lambda1=0.89361
         lambdad=1.0454
         lambda2=2.1614
         lambdaqcd=0.2452
         n=0.3919  ! divided by 2 compared to Table 1 in paper
c      
         tau=q2/(4.0*mp**2)
         qt2=q2*log((lambda2**2+q2)/lambdaqcd**2)/
     &        log(lambda2**2/lambdaqcd**2)
         f1rho=(lambda1**2/(lambda1**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))
         f1omega=f1rho
         f1phi=f1rho*(q2/(lambda1**2+q2))**1.5
         f1d=(lambdad**2/(lambdad**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))
         f2rho=(lambda1**2/(lambda1**2+qt2))**2*
     &        (lambda2**2/(lambda2**2+qt2))
         f2omega=f2rho
         f2phi=f2rho*((lambda1**2*(q2+muphi**2))/
     &        (muphi**2*(lambda1**2+q2)))**1.5
         f2d=(lambdad**2/(lambdad**2+qt2))**2*
     &        (lambda2**2/(lambda2**2+qt2))
         t1=n*((1.0317+0.0875/((1+q2/0.3176)**2))/(1+q2/0.5496))*f1rho
         t2=grhorat*mrhop**2/(mrhop**2+q2)*f1rho
         t3=(1-1.1192*n-grhorat)*f1d
c         if(q2.eq.8.)write(*,*)q2,t1,t2,t3
         f1iv=t1+t2+t3         
         t1=n*((5.7824+0.3907/((1+q2/0.1422)))/(1+q2/0.5362))*f2rho
         t2=kapparho*grhorat*mrhop**2/(mrhop**2+q2)*f2rho
c         if(q2.eq.8.)write(*,*)kappanu,6.1731*n,kapparho*grhorat,f2d
         t3=(kappanu-6.1731*n-kapparho*grhorat)*f2d
c         if(q2.eq.8.)write(*,*)q2,t1,t2,t3
         f2iv=t1+t2+t3
         t1=gomegarat*(momega**2/(momega**2+q2))*f1omega
         t2=gphirat*(mphi**2/(mphi**2+q2))*f1phi
         t3=(1-gomegarat)*f1d
         f1is=t1+t2+t3
         t1=kappaomega*gomegarat*(momega**2/(momega**2+q2))*f2omega
         t2=kappaphi*gphirat*(mphi**2/(mphi**2+q2))*f2phi
         t3=(kappas-kappaomega*gomegarat-kappaphi*gphirat)*f2d
         f2is=t1+t2+t3
         f1p=(f1is+f1iv)/2.0
         f2p=(f2is+f2iv)/2.0
         f1n=(f1is-f1iv)/2.0
         f2n=(f2is-f2iv)/2.0
         gep=f1p-tau*f2p
         gmp=f1p+f2p
         gen=f1n-tau*f2n
         gmn=f1n+f2n
         gd=1.0/(1.0+q2/0.71)**2
         ratio=mu*gep/gmp
c         write(*,100)q2,f1is,f1iv,f2is,f2iv,(gep/gd),(gmp/(mu*gd)),ratio
c
c DR-GK'(1)
c
      elseif(iset .eq. 3) then
c
         grhorat=0.0636
         kapparho=-0.4175
         gomegarat=0.7918
         kappaomega=5.1109
         gphirat=-0.3011
         kappaphi=13.4385
         muphi=1.1915
         lambda1=0.9660
         lambdad=1.3406
         lambda2=2.1382
         lambdaqcd=0.1163
         n=0.5000
c
         tau=q2/(4.0*mp**2)
         qt2=q2*log((lambdad**2+q2)/lambdaqcd**2)/
     &        log(lambdad**2/lambdaqcd**2)
         f1rho=(lambda1**2/(lambda1**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))
         f1omega=f1rho
         f1phi=f1rho*(q2/(lambda1**2+q2))**1.5
         f1d=(lambdad**2/(lambdad**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))
         f2rho=(lambda1**2/(lambda1**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))**2
         f2omega=f2rho
         f2phi=f2rho*((lambda1**2*(q2+muphi**2))/
     &        (muphi**2*(lambda1**2+q2)))**1.5
         f2d=(lambdad**2/(lambdad**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))**2
         t1=n*((1.0317+0.0875/((1+q2/0.3176)**2))/(1+q2/0.5496))*f1rho
         t2=grhorat*mrhop**2/(mrhop**2+q2)*f1rho
         t3=(1-1.1192*n-grhorat)*f1d
c         if(abs(q2-5.6).le.0.05)write(*,*)q2,t1,t2,t3,t1+t2+t3
         f1iv=t1+t2+t3         
         t1=n*((5.7824+0.3907/((1+q2/0.1422)))/(1+q2/0.5362))*f2rho
         t2=kapparho*grhorat*mrhop**2/(mrhop**2+q2)*f2rho
c         if(q2.ge.9.9)write(*,*)kappanu,6.1731*n,f2d
         t3=(kappanu-6.1731*n-kapparho*grhorat)*f2d
c         if(abs(q2-5.6).le.0.05)write(*,*)q2,t1,t2,t3,t1+t2+t3
         f2iv=t1+t2+t3
         t1=gomegarat*(momega**2/(momega**2+q2))*f1omega
         t2=gphirat*(mphi**2/(mphi**2+q2))*f1phi
         t3=(1-gomegarat)*f1d
c         if(abs(q2-5.6).le.0.05)write(*,*)q2,t1,t2,t3,t1+t2+t3
         f1is=t1+t2+t3
         t1=kappaomega*gomegarat*(momega**2/(momega**2+q2))*f2omega
         t2=kappaphi*gphirat*(mphi**2/(mphi**2+q2))*f2phi
         t3=(kappas-kappaomega*gomegarat-kappaphi*gphirat)*f2d
c         if(abs(q2-5.6).le.0.05)write(*,*)q2,t1,t2,t3,t1+t2+t3
         f2is=t1+t2+t3
         f1p=(f1is+f1iv)/2.0
         f2p=(f2is+f2iv)/2.0
         f1n=(f1is-f1iv)/2.0
         f2n=(f2is-f2iv)/2.0
         gep=f1p-tau*f2p
         gmp=f1p+f2p
         gen=f1n-tau*f2n
         gmn=f1n+f2n
         gd=1.0/(1.0+q2/0.71)**2
         ratio=mu*gep/gmp
c
c ---------------------------------------------------------------------------------
c  GKex(02S) in nucl-th/0203081 
c ---------------------------------------------------------------------------------
c
      elseif(iset .eq. 4) then
c
         grhorat=0.0401
         kapparho=6.8190
         gomegarat=0.6739
         kappaomega=0.8762
c      gphirat=-.1712   ! in article gphirat = 7.0172 must be typo
         gphirat=-.1676  
         kappaphi=7.0172
         muphi=0.8544
         lambda1=0.9407
         lambdad=1.2111
         lambda2=2.7891
         lambdaqcd=0.150
         n=0.5000
         gomegaprat=0.2552
         kappaomegap=1.4916
c
         tau=q2/(4.0*mp**2)
         qt2=q2*log((lambdad**2+q2)/lambdaqcd**2)/
     &        log(lambdad**2/lambdaqcd**2)
         f1rho=(lambda1**2/(lambda1**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))
         f1omega=f1rho
         f1phi=f1rho*(q2/(lambda1**2+q2))**1.5
         f1d=(lambdad**2/(lambdad**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))
         f2rho=(lambda1**2/(lambda1**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))**2
         f2omega=f2rho
         f2phi=f2rho*((lambda1**2*(q2+muphi**2))/
     &        (muphi**2*(lambda1**2+q2)))**1.5
         f2d=(lambdad**2/(lambdad**2+qt2))*
     &        (lambda2**2/(lambda2**2+qt2))**2
         t1=n*((1.0317+0.0875/((1+q2/0.3176)**2))/(1+q2/0.5496))*f1rho
         t2=grhorat*mrhop**2/(mrhop**2+q2)*f1rho
         t3=(1-1.1192*n-grhorat)*f1d
c         if(abs(q2-5.6).le.0.05)write(*,*)q2,t1,t2,t3,t1+t2+t3
         f1iv=t1+t2+t3         
         t1=n*((5.7824+0.3907/((1+q2/0.1422)))/(1+q2/0.5362))*f2rho
         t2=kapparho*grhorat*mrhop**2/(mrhop**2+q2)*f2rho
c         if(q2.ge.9.9)write(*,*)kappanu,6.1731*n,f2d
         t3=(kappanu-6.1731*n-kapparho*grhorat)*f2d
c         if(abs(q2-5.6).le.0.05)write(*,*)q2,t1,t2,t3,t1+t2+t3
         f2iv=t1+t2+t3
         t1=(gomegarat*(momega**2/(momega**2+q2))+
     &        gomegaprat*(momegap**2/(momegap**2+q2)))*f1omega
         t2=gphirat*(mphi**2/(mphi**2+q2))*f1phi
         t3=(1-gomegarat-gomegaprat)*f1d
c         if(abs(q2-5.6).le.0.05)write(*,*)q2,t1,t2,t3,t1+t2+t3
         f1is=t1+t2+t3
         t1=(kappaomega*gomegarat*(momega**2/(momega**2+q2))+
     &        kappaomegap*gomegaprat*(momegap**2/(momegap**2+
     &        q2)))*f2omega
         t2=kappaphi*gphirat*(mphi**2/(mphi**2+q2))*f2phi
         t3=(kappas-kappaomega*gomegarat-
     &        kappaomegap*gomegaprat-kappaphi*gphirat)*f2d
c         if(abs(q2-5.6).le.0.05)write(*,*)q2,t1,t2,t3,t1+t2+t3
         f2is=t1+t2+t3
         f1p=(f1is+f1iv)/2.0
         f2p=(f2is+f2iv)/2.0
         f1n=(f1is-f1iv)/2.0
         f2n=(f2is-f2iv)/2.0
         gep=f1p-tau*f2p
         gmp=f1p+f2p
         gen=f1n-tau*f2n
         gmn=f1n+f2n
         gd=1.0/(1.0+q2/0.71)**2
         ratio=mu*gep/gmp
c
      endif
c
      return
      end

