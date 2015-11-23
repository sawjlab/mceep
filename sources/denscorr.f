       subroutine denscorr(p,pmass,n)

       implicit none 
C------------------------------------------------------------------------------
c subroutine: denscorr
c     author: L.Todor
c       date: July 1999
c    purpose: Calculate the density effect for ionization loss
c             and shell correction coefficient
c  reference: R.M.Sternheimer, M.J.Berger & S.M. Seltzer
c             "Density Effect for Ionization Loss"
c             Atomic Data and Nuclear Data Tables 30, 261-271 (1984)
c             William R. Leo "Techniques for Nuclear Particle"
C-----------------------------------------------------------------------

       integer n
       double precision m 
       double precision beta,gamma,p,pmass,eta,oeta
       double precision x,x0,x1,c0,a,pow,i,rho0,delta0
       double precision cs1,cs2
       double precision hnu
C
       INCLUDE 'eloss.cmn'
C
        beta=p/sqrt(p**2+pmass**2)
c        write(6,*) 'Beta=',beta
        gamma=(sqrt(p**2+pmass**2))/pmass
c        write(6,*) 'Gamma=',gamma
c        write(6,*) 'prod=',beta*gamma
c       log 10 implemented below
        x=log(beta*gamma)/log(10.0D0)
C------------------------------------------------------------------------------
c  Substance identification made based on Z and A; 
c  however leave some room for A varying when there 
c  is a mixture of isotopes or in the first case
c  besides Hydrogen I accept Deuterium too
c and in the second case I accept He3 and He4
C------------------------------------------------------------------------------
        if(ztg.eq.1.and.(atg-1.D0).lt.1.1D0) then
C------------------------------------------------------------------------------
c  imep=mean excitation energy in MeV 
c  rho0=reference density in g/cm^3
C------------------------------------------------------------------------------
            imep=0.0000218
            rho0=0.06D0
            hnu=28.816d0*1.d-6*sqrt(rho*dble(ztg)/atg)
            c0=-2.d0*log(imep/hnu)-1.d0
            x0=0.4759D0
            x1=1.9215D0
            a=0.13483D0
            m=5.6249D0
            delta0=0.D0
        else if(ztg.eq.2.and.(atg-3.D0).lt.1.1D0) then
            imep=0.0000418D0
            rho0=0.00016632D0
            c0=-11.1393D0
            x0=2.2017D0
            x1=3.6122D0
            a=0.13443D0
            m=5.8347D0
            delta0=0.D0
C------------------------------------------------------------------------------
c    c0 was evaluated at rho0; make a correction to the actual density
C------------------------------------------------------------------------------
        else if(ztg.eq.6.and.(atg-12.D0).lt.1.1D0) then
            imep=0.000078D0
            if(rho.le.1.85D0) then
                rho0=1.7D0
                hnu=28.816d0*1.d-6*sqrt(rho*dble(ztg)/atg)
                c0=-2.d0*log(imep/hnu)-1.d0
                x0=0.048D0
                x1=2.5387D0
                a=0.20762D0
                m=2.9532D0
                delta0=0.14D0
           else if (rho.ge.2.1325D0) then
                rho0=2.2650D0
                hnu=28.816d0*1.d-6*sqrt(rho*dble(ztg)/atg)
                c0=-2.d0*log(imep/hnu)-1.d0
                x0=-0.0178D0
                x1=2.3415D0
                a=0.26142D0
                m=2.8697D0
                delta0=0.12D0
           else 
                rho0=2.0D0
                hnu=28.816d0*1.d-6*sqrt(rho*dble(ztg)/atg)
                c0=-2.d0*log(imep/hnu)-1.d0
                x0=-0.0351D0
                x1=2.486D0
                a=0.20240D0
                m=3.0036D0
                delta0=0.1D0
           end if
        else if(ztg.eq.13.and.(atg-27.D0).lt.1.1D0) then
            imep=0.000166D0
            rho0=2.6989D0
            hnu=28.816d0*1.d-6*sqrt(rho*dble(ztg)/atg)
            c0=-2.d0*log(imep/hnu)-1.d0
            x0=0.1708D0
            x1=3.0127D0
            a=0.08024D0
            m=3.6345D0
            delta0=0.12D0
        else if(ztg.eq.26.and.(atg-56.D0).lt.0.5D0) then
            imep=0.000286D0
            rho0=7.874D0
            hnu=28.816d0*1.d-6*sqrt(rho*dble(ztg)/atg)
            c0=-2.d0*log(imep/hnu)-1.d0
            x0=-0.0012D0
            x1=3.1531D0
            a=0.1468D0
            m=2.9632D0
            delta0=0.12D0
        else if(ztg.eq.82.and.(atg-207.D0).lt.0.5D0) then
            imep=0.000823D0
            rho0=11.35D0
            hnu=28.816d0*1.d-6*sqrt(rho*dble(ztg)/atg)
            c0=-2.d0*log(imep/hnu)-1.d0
            x0=0.3776D0
            x1=3.8073D0
            a=0.09359D0
            m=3.1608D0
            delta0=0.14D0
C------------------------------------------------------------------------------
c    For water I invent z=10 and a=18 
C------------------------------------------------------------------------------
        else if(ztg.eq.10.and.(atg-18.D0).lt.1.1D0) then
            imep=0.000075D0
            rho0=1.D0
            hnu=28.816d0*1.d-6*sqrt(rho*dble(ztg)/atg)
            c0=-2.d0*log(imep/hnu)-1.d0
            x0=0.24D0
            x1=2.8004D0
            a=0.09116D0
            m=3.4773D0
            delta0=0.D0
C------------------------------------------------------------------------------
c       if none of the above, return with delta and cshell 0 but still set
c       the imep with an empirical formula         
C------------------------------------------------------------------------------
        else  
            if (ztg.lt.13) then
               imep=12.D0*dble(ztg)+7.D0
            else
               imep=9.76D0*dble(ztg)+58.8*(dble(ztg))**(-0.19D0)
            end if
            imep=imep/1000000.D0
            cshell(n)=0.D0
            delta(n)=0.D0
            return  
        end if   
c        write(6,*) 'x=',x
        if (x.le.x0) then
           pow=exp(2.D0*(x-x0)*log(10.D0))
           delta(n)=delta0*pow
        end if
        if (x.lt.x1.and.x.gt.x0) then
           pow=exp(m*log(x1-x))
           delta(n)=4.6052D0*x+a*pow+c0
        end if
        if (x.ge.x1) then
           delta(n)=4.6052D0*x+c0
        end if
        eta=beta*gamma
        oeta=1.D0/(eta*eta)
        cs1=0.422377D0*oeta+0.0304043D0*oeta**2+0.00038106D0*oeta**3
        cs2=3.85019D0*oeta-0.1667989D0*oeta**2+0.00157955D0*oeta**3
        i=imep*1000.D0
        cshell(n)=cs1*i**2+cs2*i**3
        return
        end







