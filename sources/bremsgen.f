!!-------------------------------------------------------------------!!
!!  Full Angular Monte Carlo Simulation of Internal Bremsstrahlung   !!
!!  version 15.0 from 28 May 2004 (c) by FW                          !!
!!                                                                   !!
!!  26-aug-2004  P.E. Ulmer                                          !!
!!  I pass mceep angles with actual signs.  I then form thetae and   !!
!!  thetap, needed by this routine, as the absolute values of these  !!
!!  passed angles.  However, in forming the momentum vectors for the !!
!!  electron and proton, I do NOT explicity give the electron        !!
!!  y-component a negative value, but rather let the actual sign     !!
!!  of the angle determine the sign of this component.  This is      !!
!!  necessary since MCEEP can have the electron on either beam LEFT  !!
!!  or beam RIGHT, according to the input file setup.                !! 
!!                                                                   !!
!!-------------------------------------------------------------------!!

        subroutine bremsgen(QQ,E,thetae_mc,thetap_mc,cutoff,
     &                      Enew,Eprimenew,popnew)
        implicit none
        INCLUDE 'rada.cmn'
        double precision Emax,Lambda_flori,Lbpois,c_flori,Q,w,wtot
        double precision ran2,heaviside,energie
        external ran2,angles,heaviside
        integer brems_flag,interference,interf
        double precision E,Eprime,K,Kprime,QQ,po,pop,pp
        double precision thetae,thetap
        double precision wx,wy,wz
        double precision kx,ky,kz,kpx,kpy,kpz,ppx,ppy,ppz
        double precision Enew,Eprimenew,popnew,Emin
        double precision alpha,me,mp,me2,mp2,PI
        double precision w1,w2,w3
        double precision p1,p2,p3,p4,pp1,pp2,pp3,pp4,mu,nu

        double precision cutoff,sign_fix      ! (PEU)
        double precision thetae_mc,thetap_mc  ! (PEU)
c
        thetae = abs(thetae_mc)   ! bremsgen wants positive angles (PEU)
        thetap = abs(thetap_mc)   ! bremsgen wants positive angles (PEU)
c
cpeu        Emin = 0.001		! 1keV (0.001)
        Emin = cutoff             ! set to MCEEP cutoff value (PEU)
c
        alpha = 1.0/137.03599976d0
        me = 0.510998902d0	! electron mass in MeV
        mp = 938.271998d0	! proton mass in MeV
        me2 = 0.2611198778d0	! MeV^2
        mp2 = 880354.3422309d0	! MeV^2
        PI = 3.14159265d0
        Eprime = E - QQ/2./mp	! outgoing electron energy
        K = sqrt(E**2-me2)	! incoming electron 3-momentum
        Kprime = sqrt(Eprime**2-me2) ! outgoing electron 3-momentum
        po = mp			! incoming proton energy
        pop = mp+QQ/2./mp		! outgoing proton energy
        pp = sqrt(pop**2-mp2)

        p1 = alpha/2./PI**2*(E+K)/K*log((E+K)/(E-K))
        p2 = alpha/PI**2*log((Eprime+Kprime)/(Eprime-Kprime))
        p3 = 0.05*(1.-cos(thetae))
        p4 = alpha/2./PI**2*(pop+pp)/(pp)*log((pop+pp)/(pop-pp))

        PP1 = p1/(p1+p2+p3+p4)	! weight of e
        PP2 = p2/(p1+p2+p3+p4)	! weight of e'
        PP3 = p3/(p1+p2+p3+p4)	! weight of interference
        PP4 = p4/(p1+p2+p3+p4)	! weight of p'

 123    continue

        Q = 1
        w1 = 0.
        w2 = 0.
        w3 = 0.
        wtot = 0.
        Enew = 0.
        Eprimenew = 0.
        popnew = 0.

        kx = 0.
        ky = 0.
        kz = K  

        kpx = Kprime*sin(thetae_mc)
        kpy = 0.
        kpz = Kprime*cos(thetae_mc) 
        
        ppx = pp*sin(thetap_mc)
        ppy = 0.
        ppz = pp*cos(thetap_mc) 

 110    continue

        mu = ran2()

        if(mu.lt.pp1) then
           brems_flag = 1
           Emax = E
           call angles(wx,wy,wz,interf,1,QQ,E,thetae,thetap)
           interf = 0

        elseif(mu.gt.pp1.and.mu.lt.pp1+pp2) then
           brems_flag = 2
           Emax = Eprime
           call angles(wx,wy,wz,interf,2,QQ,E,thetae,thetap)
           interf = 0

        elseif(mu.gt.pp1+pp2.and.mu.lt.pp1+pp2+pp3) then
           interference = 1
           Emax = E
           call angles(wx,wy,wz,interf,4,QQ,E,thetae,thetap)
           if(interf.eq.1) brems_flag = 1
           if(interf.eq.2) brems_flag = 2

        else
           brems_flag = 3
           Emax = pop
           call angles(wx,wy,wz,interf,3,QQ,E,thetae,thetap)
           interf = 0

        endif

        Lambda_flori = 2.0*alpha/PI*(log(QQ/me2)+log(QQ/mp2) + 
     >       2.0*log(K/Kprime)-2.)

C
C generate a sequence n which is Poisson distributed according to 
C Lambda_flori
C

        Lbpois = Lambda_flori*log(Emax/Emin)
        c_flori = exp(-Lbpois)
        Q = Q*(ran2())
        if(Q.lt.c_flori) goto 100

C
C generate an energy according to pdf(w)=1/w between Emin and Emax:
C

        w = Emin*(Emax/Emin)**ran2()

        if (brems_flag.eq.1) then
           kx = kx - w*wy
           ky = ky + w*wx
           kz = kz - w*wz
           w1 = w1 + w
 
        else if (brems_flag.eq.2) then
           kpx = kpx - w*wy
           kpy = kpy + w*wx
           kpz = kpz - w*wz
           w2 = w2 + w
   
        else if (brems_flag.eq.3) then
           ppx = ppx - w*wy
           ppy = ppy + w*wx
           ppz = ppz - w*wz
           w3 = w3 + w
   
       end if

       wtot = wtot + w

        goto 110
 100    continue

        Enew = E - sqrt(kx**2+ky**2+kz**2+me2)
        Eprimenew = Eprime - sqrt(kpx**2+kpy**2+kpz**2+me2)
        popnew = pop - sqrt(ppx**2+ppy**2+ppz**2+mp2)

!        Enew = w1
!        Eprimenew = w2
!        popnew = w3

!       It's easier for MCEEP implementation if I make my own
!       version of the angle kicks (PEU).  The original versions
!       are commented out below.  I add winkel1-6 to my original
!       angles inside MCEEP.
!
!       Bremsgen coord     MCEEP coord     Physical def'n
!       --------------     -----------     --------------
!            z                  z           along beam
!           -x                  y           upward
!            y                  x           beam left (looking down)

        winkel1 = atan(kpy/sqrt(kpx**2+kpz**2))
        winkel2 = atan(kpx/kpz) - thetae_mc

        winkel3 = atan(ppy/sqrt(ppx**2+ppz**2))
        winkel4 = atan(ppx/ppz) - thetap_mc

        winkel5 = atan(ky/sqrt(kx**2+kz**2))
        winkel6 = atan(kx/kz)

        if(Enew.lt.-0.001) then 
           goto 123
        endif

        if(Eprimenew.lt.-0.001) then
           goto 123
        endif

        if(popnew.lt.-0.001) then 
           goto 123
        endif

        countall = countall + 1
        energie = Enew + Eprimenew + popnew

        if(energie.gt.225.) then
           countbad = countbad + 1
           goto 123
        endif

        end

!----------------------------------------------------------------------

C
C this subroutine generates events according to angular distribution
C
        subroutine angles(wx,wy,wz,interf,brems_flag,QQ,E,thetae,
     >                    thetap)
        implicit none
        INCLUDE 'rada.cmn'
        external ran2,heaviside
        double precision ran2,heaviside
        integer brems_flag,interf
        double precision E,Eprime,K,Kprime,QQ,po,pop,pp
        double precision thetae,thetap
        double precision x,phi,theta
        double precision faktor
        double precision A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,AAAA,MAYO
        double precision ok,okp,opp,kkp,kppp,kpp,wx,wy,wz
        double precision alpha,me,mp,me2,mp2,PI

        alpha=1.0/137.03599976d0
        me=0.510998902d0	! electron mass in MeV
        mp=938.271998d0		! proton mass in MeV
        me2=0.2611198778d0	!MeV^2
        mp2=880354.3422309d0	!MeV^2
        PI=3.14159265d0

        Eprime=E-QQ/2./mp	! outgoing electron energy
        K=sqrt(E**2-me2)	! incoming electron 3-momentum
        Kprime=sqrt(Eprime**2-me2) ! outgoing electron 3-momentum
        po=mp			! incoming proton energy
        pop=mp+QQ/2./mp		! outgoing proton energy
        pp=sqrt(pop**2-mp2)
        faktor=-alpha/4.0/PI**2

 299    continue
        x=0
        theta=0
        phi=0
C
C choose contribution to angular distribution according to
C the respective weights, get x, transfer x=cos theta to the
C correct coords...
C
        phi=2*PI*ran2()

        if(brems_flag.eq.1) then
           x=E/K-(E+K)/K*((E-K)/(E+K))**ran2()
           theta=acos(x)

           ok=k*cos(theta)
           okp=kprime*(-sin(thetae)*sin(theta)*sin(phi)
     >	       +cos(thetae)*cos(theta))
           opp=pp*(sin(thetap)*sin(theta)*sin(phi)
     >         +cos(thetap)*cos(theta))
           kkp=k*kprime*cos(thetae)
           kppp=kprime*pp*(-sin(thetae)*sin(thetap)
     >         +cos(thetae)*cos(thetap))
           kpp=k*pp*cos(thetap)

           A1=faktor*me**2/(Eprime-okp)**2
           A2=-2.*faktor*(E*Eprime-kkp)/(E-ok)/(Eprime-okp)
           A3=faktor*me**2/(E-ok)**2
           A4=-2.*faktor*(Eprime*pop-kppp)/(Eprime-okp)/(pop-opp)
           A5=2.*faktor*Eprime/(Eprime-okp)
           A6=2.*faktor*(E*pop-kpp)/(E-ok)/(pop-opp)
           A7=-2.*faktor*E/(E-ok)
           A8=faktor*mp**2/(pop-opp)**2
           A9=-2.*faktor*pop/(pop-opp)
           A10=faktor

           AAAA=A1+A2+A3+A4+A5+A6+A7+A8+A9+A10
           MAYO=alpha/2./PI**2*(E+K)/(E-ok)
     >         +alpha/PI**2*Kprime/(Eprime-okp)
     >         +0.05*heaviside(theta)*heaviside(thetae-theta)
     >         +alpha/2./PI**2*(pop+pp)/(pop-opp)

           if((ran2()*MAYO).gt.AAAA) goto 299

           wx=sin(theta)*cos(phi)
           wy=sin(theta)*sin(phi)
           wz=cos(theta)

        else if(brems_flag.eq.2) then
           x=Eprime/Kprime-(Eprime+Kprime)/Kprime*
     >     ((Eprime-Kprime)/(Eprime+Kprime))**ran2()
           theta=acos(x)

           ok=k*(sin(thetae)*sin(theta)*sin(phi)
     >     +cos(thetae)*cos(theta))
           okp=kprime*cos(theta)
           opp=pp*(sin(thetap+thetae)*sin(theta)*sin(phi)
     >     +cos(thetap+thetae)*cos(theta))
           kkp=k*kprime*cos(thetae)
           kppp=kprime*pp*cos(thetae+thetap)
           kpp=k*pp*(sin(thetae)*sin(thetae+thetap)
     >     +cos(thetae)*cos(thetae+thetap))

           A1=faktor*me**2/(Eprime-okp)**2
           A2=-2.*faktor*(E*Eprime-kkp)/(E-ok)/(Eprime-okp)
           A3=faktor*me**2/(E-ok)**2
           A4=-2.*faktor*(Eprime*pop-kppp)/(Eprime-okp)/(pop-opp)
           A5=2.*faktor*Eprime/(Eprime-okp)
           A6=2.*faktor*(E*pop-kpp)/(E-ok)/(pop-opp)
           A7=-2.*faktor*E/(E-ok)
           A8=faktor*mp**2/(pop-opp)**2
           A9=-2.*faktor*pop/(pop-opp)
           A10=faktor

           AAAA=A1+A2+A3+A4+A5+A6+A7+A8+A9+A10
           MAYO=alpha/2./PI**2*(E+K)/(E-ok)
     >     +alpha/PI**2*Kprime/(Eprime-okp)
     >     +0.05*heaviside(-theta)*heaviside(thetae+theta) 
     >     +alpha/2./PI**2*(pop+pp)/(pop-opp)

           if((ran2()*MAYO).gt.AAAA) goto 299

           wx=sin(theta)*cos(phi)
           wy=cos(thetae)*sin(theta)*sin(phi)-sin(thetae)*cos(theta)
           wz=sin(thetae)*sin(theta)*sin(phi)+cos(thetae)*cos(theta)

        elseif(brems_flag.eq.3) then
           x=pop/pp-(pop+pp)/pp*((pop-pp)/(pop+pp))**ran2()
           theta=acos(x)

           ok=k*(-sin(thetap)*sin(theta)*sin(phi)
     >     +cos(thetap)*cos(theta))
           okp=kprime*(-sin(thetap+thetae)*sin(theta)*sin(phi)
     >     +cos(thetap+thetae)*cos(theta))
           opp=pp*cos(theta)
           kkp=k*kprime*(sin(thetap)*sin(thetap+thetae)
     >     +cos(thetap)*cos(thetap+thetae))
           kppp=kprime*pp*cos(thetap+thetae)
           kpp=k*pp*cos(thetap)

           A1=faktor*me**2/(Eprime-okp)**2
           A2=-2.*faktor*(E*Eprime-kkp)/(E-ok)/(Eprime-okp)
           A3=faktor*me**2/(E-ok)**2
           A4=-2.*faktor*(Eprime*pop-kppp)/(Eprime-okp)/(pop-opp)
           A5=2.*faktor*Eprime/(Eprime-okp)
           A6=2.*faktor*(E*pop-kpp)/(E-ok)/(pop-opp)
           A7=-2.*faktor*E/(E-ok)
           A8=faktor*mp**2/(pop-opp)**2
           A9=-2.*faktor*pop/(pop-opp)
           A10=faktor

           AAAA=A1+A2+A3+A4+A5+A6+A7+A8+A9+A10
           MAYO=alpha/2./PI**2*(E+K)/(E-ok)
     >     +alpha/PI**2*Kprime/(Eprime-okp)
     >     +0.05*heaviside(theta-thetap)*heaviside(thetae+thetap-theta)
     >     +alpha/2./PI**2*(pop+pp)/(pop-opp)

           if((ran2()*MAYO).gt.AAAA) goto 299

           wx=sin(theta)*cos(phi)
           wy=cos(thetap)*sin(theta)*sin(phi)+sin(thetap)*cos(theta)
           wz=-sin(thetap)*sin(theta)*sin(phi)+cos(thetap)*cos(theta)

        elseif(brems_flag.eq.4) then
           theta = ran2()*thetae

!           if(theta.lt.0.5*thetae) then ! 50/50 partition
           if(ran2().lt.0.5) then ! coin toss
              interf = 1
              ok=k*cos(theta)
              okp=kprime*(-sin(thetae)*sin(theta)*sin(phi)
     >        +cos(thetae)*cos(theta))
              opp=pp*(sin(thetap)*sin(theta)*sin(phi)
     >        +cos(thetap)*cos(theta))
              kkp=k*kprime*cos(thetae)
              kppp=kprime*pp*(-sin(thetae)*sin(thetap)
     >        +cos(thetae)*cos(thetap))
              kpp=k*pp*cos(thetap)

              A1=faktor*me**2/(Eprime-okp)**2
              A2=-2.*faktor*(E*Eprime-kkp)/(E-ok)/(Eprime-okp)
              A3=faktor*me**2/(E-ok)**2
              A4=-2.*faktor*(Eprime*pop-kppp)/(Eprime-okp)/(pop-opp)
              A5=2.*faktor*Eprime/(Eprime-okp)
              A6=2.*faktor*(E*pop-kpp)/(E-ok)/(pop-opp)
              A7=-2.*faktor*E/(E-ok)
              A8=faktor*mp**2/(pop-opp)**2
              A9=-2.*faktor*pop/(pop-opp)
              A10=faktor
   
              AAAA=A1+A2+A3+A4+A5+A6+A7+A8+A9+A10
              MAYO=alpha/2./PI**2*(E+K)/(E-ok)
     >        +alpha/PI**2*Kprime/(Eprime-okp)
     >        +0.05*heaviside(theta)*heaviside(thetae-theta)
     >        +alpha/2./PI**2*(pop+pp)/(pop-opp)

              if((ran2()*MAYO).gt.AAAA) goto 299

              wx=sin(theta)*cos(phi)
              wy=sin(theta)*sin(phi)
              wz=cos(theta)

           else
              interf = 2

              ok=k*(sin(thetae)*sin(theta)*sin(phi)
     >        +cos(thetae)*cos(theta))
              okp=kprime*cos(theta)
              opp=pp*(sin(thetap+thetae)*sin(theta)*sin(phi)
     >        +cos(thetap+thetae)*cos(theta))
              kkp=k*kprime*cos(thetae)
              kppp=kprime*pp*cos(thetae+thetap)
              kpp=k*pp*(sin(thetae)*sin(thetae+thetap)
     >        +cos(thetae)*cos(thetae+thetap))

              A1=faktor*me**2/(Eprime-okp)**2
              A2=-2.*faktor*(E*Eprime-kkp)/(E-ok)/(Eprime-okp)
              A3=faktor*me**2/(E-ok)**2
              A4=-2.*faktor*(Eprime*pop-kppp)/(Eprime-okp)/(pop-opp)
              A5=2.*faktor*Eprime/(Eprime-okp)
              A6=2.*faktor*(E*pop-kpp)/(E-ok)/(pop-opp)
              A7=-2.*faktor*E/(E-ok)
              A8=faktor*mp**2/(pop-opp)**2
              A9=-2.*faktor*pop/(pop-opp)
              A10=faktor

              AAAA=A1+A2+A3+A4+A5+A6+A7+A8+A9+A10
              MAYO=alpha/2./PI**2*(E+K)/(E-ok)
     >        +alpha/PI**2*Kprime/(Eprime-okp)
     >        +0.05*heaviside(-theta)*heaviside(thetae+theta)
     >        +alpha/2./PI**2*(pop+pp)/(pop-opp)

              if((ran2()*MAYO).gt.AAAA) goto 299
      
              wx=sin(theta)*cos(phi)
              wy=cos(thetae)*sin(theta)*sin(phi)-sin(thetae)*cos(theta)
              wz=sin(thetae)*sin(theta)*sin(phi)+cos(thetae)*cos(theta)

           endif

        endif

        end

!----------------------------------------------------------------------

      double precision function heaviside(angle)
      implicit none
      external ran2,angles,bremsgen
      double precision angle

      if(angle.lt.0) then
         heaviside=0.
      else
         heaviside=1.
      endif

      return
      end

!----------------------------------------------------------------------

      double precision FUNCTION ran2()
      INTEGER IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2,idum
      integer idum 
      DATA idum/-123/
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)

      return
      END

