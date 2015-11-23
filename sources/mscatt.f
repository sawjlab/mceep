       subroutine multi_scang(thick,mom,varphi,varth,pmass,pz)
 
C------------------------------------------------------------------------------
c subroutine: multi_scang
c     author: L.Todor
c       date: July 1999
c    purpose: propose some quantities as angle variations due to multiple 
c             scattering effect
c  reference: William R. Leo "Techniques for Nuclear Particle" 
c             G.R.Lynch and O.I.Dahl "Approximation to Multiple 
c             Coulomb Scattering", Nuclear Instruments and Methods in
c             Physics Research B58 (1991) 6-10 
c
C------------------------------------------------------------------------------

       implicit none

       integer pz
       double precision mom,varphi,varth,mc2,en,beta,chic2,chia2,z,x,var
       double precision thick,pbe2,alpha,small,z23,zab,omega,v,f,th2 
       double precision pi,sigma,pmass,pzd
       real  rndm(2)
       parameter (pi=3.14159265358979324D0)
C
       INCLUDE 'eloss.cmn'
C
       mc2=pmass 
       en=sqrt(mom*mom+mc2*mc2)
       beta=mom/en
c       write(6,*) 'beta=',beta
       z=dble(ztg)
       pzd=dble(pz)
       x=thick*rho*100.D0
c       write(6,*) 'thick=',x
       alpha=1.D0/137.036D0
       pbe2=mom**2*beta**2
       chic2=0.157D0*pzd*(z*(z+1)/atg)*x/pbe2
c       write(6,*) 'chic2=',chic2,'pz=',pz,'atg=',atg 
       z23=exp((2.D0/3.D0)*dlog(z))
       small=2.007D0/100000.D0
       zab=pzd*z*alpha/beta
       chia2=small*z23*(1.D0+3.34D0*zab**2)/mom**2
c       write(6,*) 'chia2=',chia2 
C------------------------------------------------------------------------------
c not clear to how this f between 0.9 and 0.995 is chosen... 
C------------------------------------------------------------------------------
       f=0.99D0 
       omega=chic2/chia2
       v=omega/(2.D0*(1.D0-f))
       var=1.d0+v 
       th2=(2.D0*chic2/(1.D0+f**2))*((var/v)*log(var)-1.D0)       
       sigma=sqrt(th2/2.D0)
c       write(6,*) 'sigma=',sigma 
C------------------------------------------------------------------------------
c now I generate values according to a Gaussian distribution centered
c on 0 with width sigma; using Box-Muller method described in Numerical
c Recipes ch7.2 
C------------------------------------------------------------------------------
 11    call RANECU(rndm,2)
       if(rndm(1).eq.1.0.or.rndm(1).eq.0.0) go to 11
       if(rndm(2).eq.1.0.or.rndm(2).eq.0.0) go to 11
       varphi=sigma*(sqrt(-2.D0*log(dble(rndm(1)))))*
     &           cos(2.D0*pi*dble(rndm(2)))
       varth=sigma*(sqrt(-2.D0*log(dble(rndm(2)))))*
     *          cos(2.D0*pi*dble(rndm(1)))
c       write(6,*) 'varphi=',varphi       
c       write(6,*) 'varth=',varth       
       return 
       end
















