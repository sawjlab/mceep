C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       SUBROUTINE GET_SABJES_PAR
C
C       AUTHOR:  P.E. ULMER
C       DATE:    16-JUN-2000
C
C       PURPOSE: Read parameters from file deut_sabjes.dat.
C                See that file for format and description of parameters.
C                These parameters are for physics option 350 
C                (S. Jeschonnek d(e,e'p)n PWBA).
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
      SUBROUTINE get_sabjes_par
      IMPLICIT NONE
C
      COMMON /sabjes/ ipwba,iff,iwfmod,inonrel,iexfac,inos,inod
C
      INTEGER  ipwba,iff,iwfmod,inonrel,iexfac,inos,inod
C
      read(15,*) ipwba,iff,iwfmod,inonrel,iexfac,inos,inod
C
      return
      end
c
c
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c
c d(e,e'p)n code from Sabine Jeschonnek
c
c this code calculates the pwba results for d(e,e'p)n
c analytically
c-------------------------------------------------------
      subroutine wqpwba(elenin,q,omega,pmgev,protgev,
     #                  proth,protphi,pthe,phi,
     #                  ahelicity,wqfmmev)
c-------------------------------------------------------

      implicit real*8 (a-b,d-h,o-z)
      implicit complex*16 (c)
      implicit integer (i-n)

      COMMON /sabjes/ ipwba,iff,iwfmod,inonrel,iexfac,inos,inod
                         ! read from file deut_sabjes.dat

      dimension gch(-1:1,-1:1)
      dimension charnor(-1:1,-1:1,-1:1),cso(-1:1,-1:1,-1:1)
      dimension crho(-1:1,-1:1,-1:1),crhorot(-1:1,-1:1,-1:1)
      dimension cjprot(-1:1,-1:1,-1:1),cjmrot(-1:1,-1:1,-1:1)
      dimension cjmagy(-1:1,-1:1,-1:1), cjmagp(-1:1,-1:1,-1:1)
      dimension cjmagm(-1:1,-1:1,-1:1),cjp(-1:1,-1:1,-1:1)
      dimension cjconp(-1:1,-1:1,-1:1),cjconm(-1:1,-1:1,-1:1)
      dimension cjrelcop(-1:1,-1:1,-1:1),cjrelcom(-1:1,-1:1,-1:1)
      dimension cjm(-1:1,-1:1,-1:1)
      dimension prtlp(-1:1), prl(-1:1), prt(-1:1), prtt(-1:1)
      dimension prtl(-1:1), prtp(-1:1)
      dimension cjexrelp(-1:1,-1:1,-1:1),cjexrelm(-1:1,-1:1,-1:1) 
      dimension cjdiff(-1:1,-1:1,-1:1), cjsum(-1:1,-1:1,-1:1)

c dimension for the formfactors of Wally's subroutine
      dimension g(0:1,0:1)

c some constants
      pi = dacos(-1.d0)
      ci = (0.d0,1.d0)


      hbarc = 0.197327054d0
      hbarcmev = hbarc*1.d03

      alpha = 1.d0/137.0359896d0

      amp = 0.93827231d0
      amn = 0.93956563d0
      amd = 1.87561340d0


      ampsq = amp*amp
      amp4 = ampsq*ampsq
      amdsq = amd*amd


      amtarget = amd
      amejec = amp
      amres = amn

      ampfm = amp/hbarc

      cy00 = 1.d0/dsqrt(4.d0*pi)


cxxx      inos = 0     ! now read from file deut_sabjes.dat
cxxx      inod = 0     ! now read from file deut_sabjes.dat

      iconv = -1

c flag for the form factor: iff

cxxx      iff = 12      ! now read from file deut_sabjes.dat

c flag for the wave function model used

cxxx      iwfmod = 3      ! now read from file deut_sabjes.dat


c some switches
c if iqo=1, omega, q, and phi (of the missing momentum) are input
c and the misssing momentum is the running variable
c if iqo=0, the energy of the proton, the polar angle of the missin
c momentum and its azimuthal angle phi are input, and the missing
c momentum is the running variable


c fix the convention for the missing momentum
c ICONV = 1 means: pm = p - q (Rocco's convention)
c ICONV = -1 means: pm = q - p  
      if(iconv.eq.1) then
c         write(6,'(a9)')'#pm = p-q'
      elseif(iconv.eq.-1) then
c         write(6,'(a9)')'#pm = q-p'
      endif

c iff is the flag for the formfactor

      if(iwfmod.eq.1) then
c         write(6,'(a24)') '#Bonn wave function used'
      elseif(iwfmod.eq.0) then
c         write(6,'(a25)') '#Paris wave function used'
      elseif(iwfmod.eq.2) then
c         write(6,'(a26)') '#CDBonn wave function used'
      elseif(iwfmod.eq.3) then
c         write(6,'(a25)') '#A.V18 wave function used'
      elseif(iwfmod.eq.4) then
c         write(6,'(a26)') '#Wallys wave function used'
c      elseif(iwfmod.eq.5) then
c         write(6,'(a26)') '#Wallys wave f. table used'
c      else
         write(6,*)'no wave function specified',iwfmod
         stop
      endif

c      if(inos.eq.1) then
c         write(6,'(a18)')'no s wave included'
c      endif
c      if(inod.eq.1) then
c         write(6,'(a18)')'no d wave included'
c      endif
      


c improved relativistic operator, also needs the kinematics
c of the outgoing proton

c Kinetic energy of the outgoing proton in GeV

      tkinp = sqrt(protgev**2+amejec**2) - amejec

c flag for PWIA/PWBA

cxxx      ipwba = 1   ! now read from file deut_sabjes.dat


c version of current operator according to PRC57, 2438 (1998).

cxxx      inonrel = 0   ! now read from file deut_sabjes.dat
cxxx      iexfac = 1    ! now read from file deut_sabjes.dat


c option to use the exact relativistic factors
      if(inonrel.eq.1) then
         iexfac = 0
      endif
      
c      if(iexfac.eq.1) then
c         write(6,'(a30)')'#exact relativistic expressions used'
c      endif
c
c      write(6,'(a24,i5)')'#form factor param. used',iff

c
      if(ipwba.eq.1) then
c         write(6,'(a10)')'#PWBA used'
      elseif(ipwba.eq.0) then
c         write(6,'(a10)')'#PWIA used'
      else
         write(6,*) 'wrong input, IPWBA=',ipwba
         stop
      endif

      imgt = 1
      icvt = 1
      irelco = 1
      iexrel = 1
      ichar = 1
      iso = 1

      if (inonrel.eq.1) then
         irelco=0
         iso=0
         iexrel = 0
      endif
     
c      write(6,'(a30)')'#--------------------------------------------'
c      write(6,'(a39)') '#pieces of the current which contribute'
c      write(6,'(a14,2x,1i)') '#normal charge',ichar
c      write(6,'(a16,2x,1i)') '#spin orbit part',iso
c      write(6,'(a22,2x,1i)') '#magnetization current',imgt
c      write(6,'(a19,2x,1i)') '#convection current',icvt
c      write(6,'(a35,2x,1i)') 
c     ^              '#transverse relativistic correction',irelco
c      write(6,'(a30)')'#--------------------------------------------'


c read in the angles of the polarization axis for the
c initial state

c INPUT NEEDED for TARGET POLARIZATION
      rotthe = 0.d0
      rotphi = 0.d0

      rotcost = dcos(rotthe)
      rotsint = dsin(rotthe)
      rotfp = 0.5d0*(1.d0+rotcost)
      rotfm = 0.5d0*(1.d0-rotcost)
      rotfs = rotsint/dsqrt(2.d0)

      rotcosp = dcos(rotphi)
      rotsinp = dsin(rotphi)
      rotcos2p = dcos(2.d0*rotphi)
      rotsin2p = dsin(2.d0*rotphi)

      crotep = rotcosp+ci*rotsinp 
      crotem = rotcosp-ci*rotsinp
      crote2p = rotcos2p+ci*rotsin2p 
      crote2m = rotcos2p-ci*rotsin2p

c momentum of the outgoing proton in GeV, LAB system
cxxx      protgev = dsqrt(tkinp*tkinp+2.d0*amp*tkinp)  ! now passed directly
      protgevsq = protgev*protgev

      psint = dsin(pthe)
      pcost = dcos(pthe)

      
      if(omega.gt.q) then
         write(6,*)'timelike kinematics'
         stop
      endif
      
      ambda = omega/(2.d0*amp)
      appa = q/(2.d0*amp)
      tau = appa*appa-ambda*ambda

      pm = pmgev/hbarc
      pmgevsq = pmgev*pmgev

c to be consistent with the naming convention
c of the code seattle.f which i use here in part

         p = pm
         
c this is to make sure that the right values for pm
c and the wave functions are read in, even if the
c calculation does not start at pm = 5 MeV

c         if(iwfmod.eq.5) then
c            read(22,*) pmfile,um,wm,p1m,p2m
c         endif

         do mj = -1,1
            do msp = -1,1,2
               do msn = -1,1,2
                  charnor(mj,msp,msn) = 0.d0
                  cso(mj,msp,msn) = 0.d0
                  crho(mj,msp,msn) = 0.d0
                  crhorot(mj,msp,msn) = 0.d0
                  cjp(mj,msp,msn) = 0.d0
                  cjprot(mj,msp,msn) = 0.d0
                  cjm(mj,msp,msn) = 0.d0
                  cjmrot(mj,msp,msn) = 0.d0
                  cjconp(mj,msp,msn) = 0.d0
                  cjconm(mj,msp,msn) = 0.d0
                  cjmagy(mj,msp,msn) = 0.d0
                  cjmagp(mj,msp,msn) = 0.d0
                  cjmagm(mj,msp,msn) = 0.d0
                  cjrelcop(mj,msp,msn) = 0.d0
                  cjrelcom(mj,msp,msn) = 0.d0
                  enddo
               enddo
            enddo


c set up the kinematical factors


c put in the form factors

         bigqsqneg = q*q-omega*omega
         call sachs_sj(bigqsqneg,g,iff)

         gep = g(0,0)
         gmp = g(0,1)

         gen = g(1,0)
         gmn = g(1,1)



c for PWIA, switch off the neutron contribution

         if(ipwba.eq.0) then
            gen = 0.d0
            gmn = 0.d0
         endif

c Pauli Dirac form factors
         f1 = (ge+tau*gm)/(1.d0+tau)
         f2 = (gm-ge)/(1.d0+tau)


c other kinematic variables for Wally's single nucleon 
c response functions for positive energy states only

         fourq = dsqrt(q*q-omega*omega)

         epons = dsqrt(amp**2+pmgev**2)
         wdelta = amd*amd-2.d0*amd*epons
         wdeltasq = wdelta*wdelta


         sumc1 = -2.d0*tau*(f1+f2)**2+wdeltasq/(2.d0*
     ^        ampsq*amdsq)*(f1*f1+tau*f2*f2)-wdeltasq
     ^        *f2*f2/(8.d0*amp4)-wdeltasq*omega*f2*f2/
     ^        (4.d0*amp4*amd)+wdelta*omega*f1*(f1+f2)/
     ^        (ampsq*amd)

         p1 = wdeltasq/(2.d0*ampsq*amdsq)*(f1*f1+tau*f2*f2)
      
         p2 = -wdeltasq*f2*f2/(8.d0*amp4)

         p3 = -wdeltasq*omega*f2*f2/(4.d0*amp4*amd)

         p4 = wdelta*omega*f1*(f1+f2)/(ampsq*amd)


         sumc2 = 2.d0*(f1*f1+tau*f2*f2)/ampsq
         
         sumc3 = -wdelta*(f1*f1+tau*f2*f2)/(ampsq*amdsq)

         sumc4 = 0.d0

c on-shell values

         sumc1on = -2.d0*tau*(f1+f2)**2

         sumc2on = 2.d0*(f1*f1+tau*f2*f2)/ampsq

         sumc3on = 0.d0

         sumc4on = 0.d0


         pinx = dfloat(iconv)*pmgev*psint*dcos(phi)
         piny = dfloat(iconv)*pmgev*psint*dsin(phi)
         pinz = dfloat(iconv)*pmgev*pcost

         cpinp = -(pinx+ci*piny)/dsqrt(2.d0)
         cpinpst = conjg(cpinp)
         cpinm = (pinx-ci*piny)/dsqrt(2.d0)
         cpinmst = conjg(cpinm)

         pineps0 = (amd-epons)*q/fourq-pinz*omega/fourq

c         if(ionshell.eq.1) then
c            pineps0 = epons*q/fourq-pinz*omega/fourq
c         endif

         pineps0on = epons*q/fourq-pinz*omega/fourq

c         write(6,*) sumc1,sumc1on,sumc2,sumc3
c         write(6,*) pineps0,pineps0on

c calculate the single nucleon responses:

         w00 = sumc1+sumc2*pineps0**2+sumc3*2.d0*amd*q/fourq
     ^        *pineps0

         cwpp = sumc1*(-1.d0)+sumc2*cpinp*cpinpst

         cwmm = sumc1*(-1.d0)+sumc2*cpinm*cpinmst

         cwpm = sumc2*(cpinpst*cpinm)

         cwmp = sumc2*(cpinmst*cpinp)

         cw0p = sumc2*(-pineps0*cpinp)-sumc3*amd*q/fourq*
     ^        cpinp

         cwp0 = sumc2*(-cpinpst*pineps0)-sumc3*amd*q/fourq
     ^        *cpinpst

         cw0m = sumc2*(-pineps0*cpinm)-sumc3*amd*q/fourq
     ^        *cpinm

         cwm0 = sumc2*(-cpinmst*pineps0)-sumc3*amd*q/fourq
     ^        *cpinmst

         wrlsn = 0.5d0*w00*appa**2/tau

         cwrtsn = 0.5d0*(cwpp+cwmm)

         cwrttsn = 0.5d0*(cwpm+cwmp)

         cwrtlsn = -0.5d0*(-cw0p+cw0m-cwp0+cwm0)*appa/dsqrt(tau)

c calculate the on-shell single nucleon responses:
c for the factors -, appa, tau, see the relation between the
c Dmitrasinovic and Gross leptonic coeeficients and Bill vk,
c Table II of the Dmitrasinovic and Gross paper

         w00on = sumc1on+sumc2on*pineps0on**2

         cwppon = sumc1on*(-1.d0)+sumc2on*cpinp*cpinpst

         cwmmon = sumc1on*(-1.d0)+sumc2on*cpinm*cpinmst

         cwpmon = sumc2on*(cpinpst*cpinm)

         cwmpon = sumc2on*(cpinmst*cpinp)

         cw0pon = sumc2on*(-pineps0on*cpinp)

         cwp0on = sumc2on*(-cpinpst*pineps0on)

         cw0mon = sumc2on*(-pineps0on*cpinm)

         cwm0on = sumc2on*(-cpinmst*pineps0on)

         wrlsnon = 0.5d0*w00on*appa**2/tau

         cwrtsnon = 0.5d0*(cwppon+cwmmon)

         cwrttsnon = 0.5d0*(cwpmon+cwmpon)

         cwrtlsnon = -0.5d0*(-cw0pon+cw0mon-cwp0on+cwm0on)
     ^        *appa/dsqrt(tau)



c prepare the other dimensionless variables and 
c quantities needed for the relativistic expressions 
c for the current matrix elements

c in principle, eta for the proton goes to -eta for the neutron
c but it enters only quadratically in the equations, so f0 
c is the same for proton and neutron terms

c calculate f0, xi, xip:
               

         eta = pmgev/amp
         delta = eta*psint
         deltasq = delta*delta
         epsi = dsqrt(1.d0+eta*eta)

c prepare values for the neutron
         etan = protgev/amp
         deltan = etan*dsin(proth)
         deltasqn = deltan*deltan
         epsin = dsqrt(1.d0+etan*etan)

            if(iexfac.eq.1) then
               amu1 = appa*dsqrt(1.d0+tau)/(dsqrt(tau)*(epsi+ambda))
               amu2 = 2.d0*appa*dsqrt(1.d0+tau)/dsqrt(tau)
               amu2 = amu2/(1.d0+tau+epsi+ambda)

               amu1n = appa*dsqrt(1.d0+tau)/(dsqrt(tau)*(epsin+ambda))
               amu2n = 2.d0*appa*dsqrt(1.d0+tau)/dsqrt(tau)
               amu2n = amu2n/(1.d0+tau+epsin+ambda)

               auxx = delta*delta/(4.d0*(1.d0+tau))
               f0 = 1.d0/(amu1*dsqrt(1.d0+tau*amu2*amu2*auxx))

               auxxn = deltasqn/(4.d0*(1.d0+tau))
               f0n = 1.d0/(amu1n*dsqrt(1.d0+tau*amu2n*amu2n*auxxn))

               auxp = amu1*amu2*deltasq*tau*gmp*0.5d0
               xi0prot = appa/dsqrt(tau)*(gep+auxp/(1.d0+tau))
               auxn = amu1n*amu2n*deltasqn*tau*gmn*0.5d0
               xi0neu = appa/dsqrt(tau)*(gen+auxn/(1.d0+tau))

               xi0pprot = (amu1*gmp-0.5d0*amu2*gep)
     ^              /dsqrt(1.d0+tau)
               xi0pneu = (amu1n*gmn-0.5d0*amu2n*gen)
     ^              /dsqrt(1.d0+tau)

                  
               xi1prot = (amu1*gep+0.5d0*amu2*tau*gmp)
     ^              /dsqrt(1.d0+tau)
               xi1neu = (amu1n*gen+0.5d0*amu2n*tau*gmn)
     ^              /dsqrt(1.d0+tau)


               aux = amu1*amu2*deltasq*0.5d0/(1.d0+tau) 
               auxn = amu1n*amu2n*deltasqn*0.5d0/(1.d0+tau)

               xi1pprot = dsqrt(tau)/appa*gmp*(1.d0-aux)
               xi1pneu = dsqrt(tau)/appa*gmn*(1.d0-auxn)

               xi2pprot = 0.5d0*ambda*dsqrt(tau)*amu1*amu2*
     ^              gmp/appa**3
               xi2pneu = 0.5d0*ambda*dsqrt(tau)*amu1n*amu2n*
     ^              gmn/appa**3


               xi3pprot = dsqrt(tau)*amu1*amu2*(gep-gmp)*0.5d0
     ^              /(appa*(1.d0+tau))
               xi3pneu = dsqrt(tau)*amu1n*amu2n*(gen-gmn)*0.5d0
     ^              /(appa*(1.d0+tau))

            endif

c values for first order relativistic current
            if(iexfac.eq.0.and.inonrel.eq.0) then
               f0 = 1.d0
               f0n = 1.d0

               xi0prot = appa/dsqrt(tau)*gep
               xi0neu = appa/dsqrt(tau)*gen

               xi0pprot = (gmp-0.5d0*gep)/dsqrt(1.d0+tau)
               xi0pneu = (gmn-0.5d0*gen)/dsqrt(1.d0+tau)

               xi1prot = dsqrt(tau)/appa*(gep+0.5d0*tau*gmp)
               xi1neu = dsqrt(tau)/appa*(gen+0.5d0*tau*gmn)

               xi1pprot = dsqrt(tau)/appa*gmp
               xi1pneu = dsqrt(tau)/appa*gmn

c alternative form
               xi2pprot = dsqrt(tau)/(appa*2.d0*(1.d0+tau))*gmp
               xi2pneu = dsqrt(tau)/(appa*2.d0*(1.d0+tau))*gmn

               xi3pprot = 0.d0
               xi3pneu = 0.d0
            endif

c values for the nonrelativistic case
            if(inonrel.eq.1) then
               f0 = 1.d0
               f0n = 1.d0

               xi0prot = gep
               xi0neu = gen

               xi0pprot = 0.d0
               xi0pneu = 0.d0

               xi1prot = gep
               xi1neu = gen

               xi1pprot = gmp
               xi1pneu = gmn

               xi2pprot = 0.d0
               xi2pneu = 0.d0

               xi3pprot = 0.d0
               xi3pneu = 0.d0
            endif

c check if the current component is supposed to be calculated

            xi0prot = xi0prot*dfloat(ichar)
            xi0neu = xi0neu*dfloat(ichar)
            xi0pprot = xi0pprot*dfloat(iso)
            xi0pneu = xi0pneu*dfloat(iso)

            xi1prot = xi1prot*dfloat(icvt)
            xi1neu = xi1neu*dfloat(icvt)
            xi1pprot = xi1pprot*dfloat(imgt)
            xi1pneu = xi1pneu*dfloat(imgt)

            xi2pprot = xi2pprot*dfloat(irelco)
            xi2pneu = xi2pneu*dfloat(irelco)
            xi3pprot = xi3pprot*dfloat(iexrel)
            xi3pneu = xi3pneu*dfloat(iexrel)



         protmom = protgev/hbarc



         if(iwfmod.eq.1) then
            call momwave(pm,um,wm)
            call momwave(protmom,um1,wm1)
         elseif(iwfmod.eq.0) then
            call pamomwave(pm,um,wm)
            call pamomwave(protmom,um1,wm1)
         elseif(iwfmod.eq.2) then
            call cdbmomwave(pm,um,wm)
            call cdbmomwave(protmom,um1,wm1)
         elseif(iwfmod.eq.3) then
            call av18momwave(pm,um,wm)
            call av18momwave(protmom,um1,wm1)
         elseif(iwfmod.eq.4) then
            call wallymomwave(pm,um,wm)
            call wallymomwave(protmom,um1,wm1)
c         elseif(iwfmod.eq.5) then
c            read(22,*) pmfile,um,wm,p1m,p2m
c            um = um/(4.d0*pi)
c            wm = -wm/(4.d0*pi)
c            if(abs(pmfile-pm).ge.1.d-05) then
c               write(6,*) 'pm discrepancy in reading'
c               write(6,*) pm,pmfile,pmmev
c               stop
c            endif
         endif


c there must be a misprint in the paper,
c without the sign, no agreement between the numerical
c integration and the analytic result.
c also the other potentials (nijmegen) give a different sign
         wm = -wm
         wm1 = -wm1


         if (inod.eq.1) then
            wm = 0.d0
            wm1 = 0.d0
         endif

         if(inos.eq.1) then
            um = 0.d0
            um1 = 0.d0
         endif


c table of clebsch gordan coefficients for half integer values of m
c (1/2 arg1/2 1/2 arg2/2 | 1 (arg1+arg2)/2 )

         gch(1,1) = 1.d0
         gch(-1,-1) = 1.d0
         gch(-1,1) = dsqrt(0.5d0)
         gch(1,-1) = dsqrt(0.5d0)


c the index "2" refers to particle 2, i.e. the neutron,
c the index "1" refers to particle 1, i.e. the proton
c the PWBA matrix element is given by GEp*psi(2) + GEn*psi(1)


c spin-orbit matrix element

c here in this code, pm = iconv*pinitial
         etax = dfloat(iconv)*pmgev/amp*psint*dcos(phi)
         etay = dfloat(iconv)*pmgev/amp*psint*dsin(phi)

         do mj = -1,1
            do msp = -1,1,2
               do msn = -1,1,2

                  cfac = ci*4.d0*pi*appa

                  cfacprot = cfac*xi0pprot*gch(-msp,msn)
     ^                 *(ci*dfloat(msp)*etax+etay)*f0

                  cfacneu = cfac*xi0pneu*gch(msp,-msn)
     ^                 *(etay+ci*dfloat(msn)*etax)*f0n


                  msprot = (msn-msp)/2
                  mlprot = mj-msprot

                  msneu = (msp-msn)/2
                  mlneu = mj-msneu

c s wave 
                  cswav2 = -gc(0,mlprot,1,msprot,1,mj)*cy00*um

                  cswav1 = gc(0,mlneu,1,msneu,1,mj)*cy00*um1


c d wave contribution
                  if(mlprot.le.2.and.mlprot.ge.-2) then
                     cylm = cy2(mlprot,pthe,phi)
                  else
                     cylm = 0.d0
                  endif

                  if(mlneu.le.2.and.mlneu.ge.-2) then
                     cylm1 = cy2(mlneu,proth,protphi)
                  else
                     cylm1 = 0.d0
                  endif

                  cdwav2 = cylm*wm*gc(2,mlprot,1,msprot,1,mj)
                  cdwav1 = -cylm1*wm1*gc(2,mlneu,1,msneu,1,mj)


                  cso(mj,msp,msn) = cfacprot*(cswav2+cdwav2)+
     ^                 cfacneu*(cswav1+cdwav1)


               enddo
            enddo
         enddo


c normal charge operator


         facharp = xi0prot*4.d0*pi*f0
         facharn = xi0neu*4.d0*pi*f0n


         do mj = -1,1
            do msp = -1,1,2
               do msn = -1,1,2

                  ms = (msn+msp)/2
                  ml = mj-ms

                  cswav2 = gc(0,ml,1,ms,1,mj)*cy00*um
                  cswav1 = gc(0,ml,1,ms,1,mj)*cy00*um1


c d wave contribution
                  if(ml.le.2.and.ml.ge.-2) then
                     cylm = cy2(ml,pthe,phi)
                     cylm1 = cy2(ml,proth,protphi)
                  else
                     cylm = 0.d0
                     cylm1 = 0.d0
                  endif

                  cdwav2 = -cylm*wm*gc(2,ml,1,ms,1,mj)
                  cdwav1 = -cylm1*wm1*gc(2,ml,1,ms,1,mj)

                  charnor(mj,msp,msn) = 
     ^                 facharp*(cswav2+cdwav2)*gch(msp,msn)
     ^                 +facharn*(cswav1+cdwav1)*gch(msp,msn)



               enddo
            enddo
         enddo



c magnetization current


         facmagprot = -dsqrt(2.d0)*4.d0*pi*f0*appa*xi1pprot
         facmagneu = -dsqrt(2.d0)*4.d0*pi*f0n*appa*xi1pneu



         do mj = -1,1
            do msn = -1,1,2

               msp = 1
               ms = (msn-1)/2
               ml = mj-ms

               cswav = cy00*um*gc(0,ml,1,ms,1,mj)

               if(ml.le.2.and.ml.ge.-2) then
                     cylm = cy2(ml,pthe,phi)
                  else
                     cylm = 0.d0
                  endif

                  cdwav = cylm*wm*gc(2,ml,1,ms,1,mj)
        
                  cmjp = (cswav-cdwav)*gch(-1,msn)
                  
                  cjmagp(mj,msp,msn) = cmjp*facmagprot
                  cjmagp(mj,-msp,msn) = 0.d0
            enddo
         enddo

         do mj = -1,1
            do msp = -1,1,2

               msn = 1
               ms = (msp-1)/2
               ml = mj-ms

               cswav = cy00*um1*gc(0,ml,1,ms,1,mj)

               if(ml.le.2.and.ml.ge.-2) then
                     cylm = cy2(ml,proth,protphi)
                  else
                     cylm = 0.d0
                  endif

                  cdwav = cylm*wm1*gc(2,ml,1,ms,1,mj)
        
                  cmjp = (cswav-cdwav)*gch(msp,-1)
                  
                  cjmagp(mj,msp,msn) = 
     ^                 cjmagp(mj,msp,msn)+cmjp*facmagneu
            enddo
         enddo

         do mj = -1,1
            do msn = -1,1,2

               msp = -1
               ms = (msn+1)/2
               ml = mj-ms

               cswav = cy00*um*gc(0,ml,1,ms,1,mj)

               if(ml.le.2.and.ml.ge.-2) then
                     cylm = cy2(ml,pthe,phi)
                  else
                     cylm = 0.d0
                  endif

                  cdwav = cylm*wm*gc(2,ml,1,ms,1,mj)
        
                  cmjm = (cswav-cdwav)*gch(1,msn)
                  
                  cjmagm(mj,msp,msn) = cmjm*facmagprot
                  cjmagm(mj,-msp,msn) = 0.d0
            enddo
         enddo

         do mj = -1,1
            do msp = -1,1,2

               msn = -1
               ms = (msp+1)/2
               ml = mj-ms

               cswav = cy00*um1*gc(0,ml,1,ms,1,mj)

               if(ml.le.2.and.ml.ge.-2) then
                     cylm = cy2(ml,proth,protphi)
                  else
                     cylm = 0.d0
                  endif

                  cdwav = cylm*wm1*gc(2,ml,1,ms,1,mj)
        
                  cmjm = (cswav-cdwav)*gch(msp,1)
                  
                  cjmagm(mj,msp,msn) = 
     ^                 cjmagm(mj,msp,msn)+cmjm*facmagneu
            enddo
         enddo

c convection current

         facconprot = 4.d0*pi/dsqrt(2.d0)*f0*xi1prot
         facconneu = 4.d0*pi/dsqrt(2.d0)*f0n*xi1neu

         do mj = -1,1
            do msp = -1,1,2
               do msn = -1,1,2

                  ms = (msp+msn)/2
                  ml = mj-ms

                  cswav = um*cy00*gc(0,ml,1,ms,1,mj)
                  cswav1 = um1*cy00*gc(0,ml,1,ms,1,mj)
                              
                  if(ml.le.2.and.ml.ge.-2) then
                     cylm = cy2(ml,pthe,phi)
                     cylm1 = cy2(ml,proth,protphi)
                  else
                     cylm = 0.d0
                     cylm1 = 0.d0
                  endif

                  cdwav = cylm*wm*gc(2,ml,1,ms,1,mj)
                  cdwav1 = cylm1*wm1*gc(2,ml,1,ms,1,mj)

                  cmjcon = gch(msp,msn)*(cswav-cdwav)
                  cmjcon1 = gch(msp,msn)*(cswav1-cdwav1)

                  cjconp(mj,msp,msn) = 
     ^                 -(etax+ci*etay)*cmjcon*facconprot
     ^                 +(etax+ci*etay)*cmjcon1*facconneu
                  cjconm(mj,msp,msn) = 
     ^                 (etax-ci*etay)*cmjcon*facconprot
     ^                 -(etax-ci*etay)*cmjcon1*facconneu

               enddo
            enddo
         enddo


c first-order convective spin-orbit: (also known as Jrelcor)

         frelcor = -4.d0*pi*appa*appa/dsqrt(2.d0)
         frelcorprot = frelcor*xi2pprot*f0
         frelcorneu = frelcor*xi2pneu*f0n

         do mj = -1,1
            do msp = -1,1,2
               do msn = -1,1,2

                  ms = (msp+msn)/2
                  ml = mj-ms

                  cswav = um*cy00*gc(0,ml,1,ms,1,mj)
                  cswav1 = um1*cy00*gc(0,ml,1,ms,1,mj)
                              
                  if(ml.le.2.and.ml.ge.-2) then
                     cylm = cy2(ml,pthe,phi)
                     cylm1 = cy2(ml,proth,protphi)
                  else
                     cylm = 0.d0
                     cylm1 = 0.d0
                  endif

                  cdwav = cylm*wm*gc(2,ml,1,ms,1,mj)
                  cdwav1 = cylm1*wm1*gc(2,ml,1,ms,1,mj)

                  cmjrc = gch(msp,msn)*(cswav-cdwav)*dfloat(msp)
     ^                 *frelcorprot
                  cmjrc1 = gch(msp,msn)*(cswav1-cdwav1)*dfloat(msn)
     ^                 *frelcorneu


                  cjrelcop(mj,msp,msn) = 
     ^                 (etax+ci*etay)*cmjrc
     ^                 -(etax+ci*etay)*cmjrc1
                  cjrelcom(mj,msp,msn) = 
     ^                 (etax-ci*etay)*cmjrc
     ^                 -(etax-ci*etay)*cmjrc1

               enddo
            enddo
         enddo


c second-order convective spin-orbit: (also known as Jexrel)

         xi3pprot = dsqrt(tau)*amu1*amu2*(gep-gmp)/
     ^                 (2.d0*appa*(1.d0+tau))

         xi3pneu = dsqrt(tau)*amu1*amu2*(gen-gmn)/
     ^                 (2.d0*appa*(1.d0+tau))

         cfexrel = ci*appa*4.d0*pi/dsqrt(2.d0)

         if(iexrel.ne.1) cfexrel = 0.d0

         cfexrelprot = cfexrel*xi3pprot*f0
         cfexrelneu = cfexrel*xi3pneu*f0n


         do mj = -1,1
            do msp = -1,1,2
               do msn = -1,1,2

                  msprot = (-msp+msn)/2
                  mlprot = mj-msprot

                  msneu = (msp-msn)/2
                  mlneu = mj-msneu

                  cswav = um*cy00*gc(0,mlprot,1,msprot,1,mj)
                  cswav1 = um1*cy00*gc(0,mlneu,1,msneu,1,mj)
                              
                  if(mlprot.le.2.and.mlprot.ge.-2) then
                     cylm = cy2(mlprot,pthe,phi)
                  else
                     cylm = 0.d0
                  endif

                 if(mlneu.le.2.and.mlneu.ge.-2) then
                    cylm1 = cy2(mlneu,proth,protphi) 
                  else
                     cylm1 = 0.d0
                  endif


                 cdwav = cylm*wm*gc(2,mlprot,1,msprot,1,mj)
                 cdwav1 = cylm1*wm1*gc(2,mlneu,1,msneu,1,mj)



                 cmjxr = gch(-msp,msn)*cfexrelprot*(cswav-cdwav)
                 cmjxr1 = gch(msp,-msn)*cfexrelneu*(cswav1-cdwav1)


                  cjexrelp(mj,msp,msn) = (etax+ci*etay)*
     ^                 (-ci*dfloat(msp)*etax-etay)*cmjxr-
     ^                (etax+ci*etay)*(ci*dfloat(msn)*etax+etay)
     ^                *cmjxr1

                  cjexrelm(mj,msp,msn) = -(etax-ci*etay)*
     ^                 (-ci*dfloat(msp)*etax-etay)*cmjxr+
     ^                 (etax-ci*etay)*(ci*dfloat(msn)*etax+etay)
     ^                *cmjxr1

               enddo
            enddo
         enddo




c build up the charge operator matrix elements:
         do mj = -1,1
            do msp = -1,1,2
               do msn = -1,1,2
                  crho(mj,msp,msn) = charnor(mj,msp,msn)
     ^                  +cso(mj,msp,msn)

               enddo
            enddo
         enddo


c build up the current matrix elements:

         do mj = -1,1
            do msp = -1,1,2
               do msn = -1,1,2
                  cjp(mj,msp,msn) = cjmagp(mj,msp,msn)
     ^          +cjconp(mj,msp,msn)+cjrelcop(mj,msp,msn)
     ^          +cjexrelp(mj,msp,msn)
                  cjm(mj,msp,msn) = cjmagm(mj,msp,msn)
     ^          +cjconm(mj,msp,msn)+cjrelcom(mj,msp,msn)
     ^          +cjexrelm(mj,msp,msn)
               enddo
            enddo
         enddo

c build up an auxiliary field J+ - J- which will be used when calculating
c the transverse longitudinal response function
c and an auxiliary field J+ + J- which will be used when calculating
c the transverse longitudinal primed response function

         do mj = -1,1
            do msp = -1,1,2
               do msn = -1,1,2
                  cjdiff(mj,msp,msn) = 
     ^            cjp(mj,msp,msn) - cjm(mj,msp,msn)
                  cjsum(mj,msp,msn) =
     ^            cjp(mj,msp,msn) + cjm(mj,msp,msn)
               enddo
            enddo
         enddo




c+++++++++++++ROTATION OF THE POLARIZATION AXIS+++++++++++

c build up the matrix elements with respect to the rotated
c polarization axis for the initial state

         do msp = -1,1,2
            do msn = -1,1,2
               crhorot(1,msp,msn) = rotfp*crho(1,msp,msn)
     ^          -rotfs*crotep*crho(0,msp,msn)
     ^          +rotfm*crote2p*crho(-1,msp,msn)

               crhorot(0,msp,msn) = rotfs*crotem*crho(1,msp,msn)
     ^          +rotcost*crho(0,msp,msn)
     ^          -rotfs*crotep*crho(-1,msp,msn)

               crhorot(-1,msp,msn) = rotfm*crote2m*crho(1,msp,msn)
     ^           +rotfs*crotem*crho(0,msp,msn)
     ^           +rotfp*crho(-1,msp,msn)

               cjprot(1,msp,msn) = rotfp*cjp(1,msp,msn)
     ^          -rotfs*crotep*cjp(0,msp,msn)
     ^          +rotfm*crote2p*cjp(-1,msp,msn)

               cjprot(0,msp,msn) = rotfs*crotem*cjp(1,msp,msn)
     ^          +rotcost*cjp(0,msp,msn)
     ^          -rotfs*crotep*cjp(-1,msp,msn)

               cjprot(-1,msp,msn) = rotfm*crote2m*cjp(1,msp,msn)
     ^           +rotfs*crotem*cjp(0,msp,msn)
     ^           +rotfp*cjp(-1,msp,msn)

               cjmrot(1,msp,msn) = rotfp*cjm(1,msp,msn)
     ^          -rotfs*crotep*cjm(0,msp,msn)
     ^          +rotfm*crote2p*cjm(-1,msp,msn)

               cjmrot(0,msp,msn) = rotfs*crotem*cjm(1,msp,msn)
     ^          +rotcost*cjm(0,msp,msn)
     ^          -rotfs*crotep*cjm(-1,msp,msn)

               cjmrot(-1,msp,msn) = rotfm*crote2m*cjm(1,msp,msn)
     ^           +rotfs*crotem*cjm(0,msp,msn)
     ^           +rotfp*cjm(-1,msp,msn)
            enddo
         enddo

C++++++++++++++++++++++END OF ROTATION++++++++++++++++++++++++++++++

c build up responses

         rl = 0.d0
         do mj = -1,1
            prl(mj) = 0.d0
            do msp = -1,1,2
               do msn = -1,1,2
                  rl = rl
     ^              +crhorot(mj,msp,msn)*conjg(crhorot(mj,msp,msn))
                  prl(mj) = prl(mj)
     ^              +crhorot(mj,msp,msn)*conjg(crhorot(mj,msp,msn))
               enddo
            enddo
         enddo


         rt = 0.d0
         do mj = -1,1
            prt(mj) = 0.d0
            do msp = -1,1,2
               do msn = -1,1,2
                  rt = rt 
     ^                 + cjprot(mj,msp,msn)*conjg(cjprot(mj,msp,msn))
     ^                 + cjmrot(mj,msp,msn)*conjg(cjmrot(mj,msp,msn))
                  prt(mj) = prt(mj)
     ^                 + cjprot(mj,msp,msn)*conjg(cjprot(mj,msp,msn))
     ^                 + cjmrot(mj,msp,msn)*conjg(cjmrot(mj,msp,msn))
               enddo
            enddo
         enddo




         rtt = 0.d0
         do mj = -1,1
            prtt(mj) = 0.d0
            do msp = -1,1,2
               do msn = -1,1,2
                  rtt = rtt 
     ^            +2.d0*dreal(conjg(cjprot(mj,msp,msn))
     ^                       *cjmrot(mj,msp,msn))
                  prtt(mj) = prtt(mj)
     ^            +2.d0*dreal(conjg(cjprot(mj,msp,msn))
     ^                       *cjmrot(mj,msp,msn))
               enddo
            enddo
         enddo




         rtl = 0.d0
         do mj = -1,1
            prtl(mj) = 0.d0
            do msp = -1,1,2
               do msn = -1,1,2
                  rtl = rtl 
     ^             -2.d0*dreal(conjg(crhorot(mj,msp,msn))*
     ^             (cjprot(mj,msp,msn)-cjmrot(mj,msp,msn)))
                  prtl(mj) = prtl(mj)
     ^             -2.d0*dreal(conjg(crhorot(mj,msp,msn))*
     ^             (cjprot(mj,msp,msn)-cjmrot(mj,msp,msn)))
               enddo
            enddo
         enddo




         rtp = 0.d0
         do mj = -1,1
            prtp(mj) = 0.d0
            do msp = -1,1,2
               do msn = -1,1,2
                  rtp = rtp 
     ^                 + cjprot(mj,msp,msn)*conjg(cjprot(mj,msp,msn))
     ^                 - cjmrot(mj,msp,msn)*conjg(cjmrot(mj,msp,msn))
                  prtp(mj) = prtp(mj)
     ^                 + cjprot(mj,msp,msn)*conjg(cjprot(mj,msp,msn))
     ^                 - cjmrot(mj,msp,msn)*conjg(cjmrot(mj,msp,msn))
               enddo
            enddo
         enddo


         rtlp = 0.d0
         do mj = -1,1
            prtlp(mj) = 0.d0
            do msp = -1,1,2
               do msn = -1,1,2
                  rtlp = rtlp 
     ^             -2.d0*dreal(conjg(crhorot(mj,msp,msn))*
     ^             (cjprot(mj,msp,msn)+cjmrot(mj,msp,msn)))
                  prtlp(mj) = prtlp(mj)
     ^             -2.d0*dreal(conjg(crhorot(mj,msp,msn))*
     ^             (cjprot(mj,msp,msn)+cjmrot(mj,msp,msn)))
               enddo
            enddo
         enddo
 

         upanorm = 3.d0
         anorm = 1.d0

c unpolarized responses

      rl = rl/upanorm
      rt = rt/upanorm
      rtt = rtt/upanorm
      rtl = rtl/upanorm
      rtp = rtp/upanorm
      rtlp = rtlp/upanorm


c polarized responses

      do mj = -1,1
         prl(mj) = prl(mj)/anorm
         prt(mj) = prt(mj)/anorm
         prtt(mj) = prtt(mj)/anorm
         prtl(mj) = prtl(mj)/anorm
         prtp(mj) = prtp(mj)/anorm
         prtlp(mj) = prtlp(mj)/anorm
      enddo


c convert the dimension of the response functions from fm**3 to GeV**(-3)

      dimfac = hbarc**3

      rl = rl/dimfac
      rt = rt/dimfac
      rtt = rtt/dimfac
      rtl = rtl/dimfac
      rtp = rtp/dimfac
      rtlp = rtlp/dimfac

      specfun = specfun/dimfac
      specfunt = specfunt/dimfac
      specfuntt = specfuntt/dimfac
      specfuntl = specfuntl/dimfac
      specfuntlp = specfuntlp/dimfac



      do mj = -1,1
         prl(mj) = prl(mj)/dimfac
         prt(mj) = prt(mj)/dimfac
         prtt(mj) = prtt(mj)/dimfac
         prtl(mj) = prtl(mj)/dimfac
         prtp(mj) = prtp(mj)/dimfac
         prtlp(mj) = prtlp(mj)/dimfac
      enddo



c construction of the cross section

c calculation of the leptonic coefficients in the extreme
c relativistic limit

c calculate the electron angle from the beam energy and q,omega
c for ibeam=1, otherwise calculate the beam energy from the angle

      bigqsq = omega*omega-q*q

      efin = elenin-omega
      aux = dsqrt(-bigqsq*0.25d0/(elenin*efin))
      elthhalf = dasin(aux)
      elth = elthhalf*2.d0

      elthdeg = elth/pi*180.d0

      qrat = bigqsq/(q*q)
      tanelth = dtan(elth*0.5d0)


      vl = qrat*qrat
      vt = -0.5d0*qrat+tanelth**2
      vtt = 0.5d0*qrat
      vtl = qrat*dsqrt(-qrat+tanelth**2)/dsqrt(2.d0)
      vtp = dsqrt(-qrat+tanelth**2)*tanelth
      vtlp = qrat*tanelth/dsqrt(2.d0)


      product = vl*rl+vt*rt+vtt*rtt+vtl*rtl+
     ^     ahelicity*(vtp*rtp+vtlp*rtlp)

      product1 = vl*prl(1)+vt*prt(1)+vtt*prtt(1)+vtl*prtl(1)+
     ^     ahelicity*(vtp*prtp(1)+vtlp*prtlp(1))

      product0 = vl*prl(0)+vt*prt(0)+vtt*prtt(0)+vtl*prtl(0)+
     ^     ahelicity*(vtp*prtp(0)+vtlp*prtlp(0))

      productm1 = vl*prl(-1)+vt*prt(-1)
     ^     +vtt*prtt(-1)+vtl*prtl(-1)+
     ^     ahelicity*(vtp*prtp(-1)+vtlp*prtlp(-1))

c mott cross section (also in the extreme relativistic limit)

      angular = dcos(elth*0.5d0)/((dsin(elth*0.5d0))**2)
      sigmott = (angular*alpha/(2.d0*elenin))**2
      sigmottfm = sigmott*hbarc**2
      sigmottnb = 1.d07*sigmott*hbarc**2

c recoil factor

      enejec = tkinp+amejec
      recoil = abs(1.d0+(omega*protgev-enejec*q*dcos(proth))
     ^     /(amtarget*protgev))

c kinematic factor

      fakima = amejec*amres*protgev/(8.d0*pi*pi*pi*amtarget)

c units: Mott cross section GeV**-2  sr**-1
c        kinematic factor GeV**2
c        recoil factor is dimensionless
c so the dimensions stem from the response functions,
c which have GeV**-3 sr**-1
c cross section (units are GeV**-3)

      wq = fakima/recoil*sigmott*product


      wqpol1 = fakima/recoil*sigmott*product1
      wqpol0 = fakima/recoil*sigmott*product0
      wqpolm1 = fakima/recoil*sigmott*productm1

c in order to convert everything, it is multiplied by dimfac
c which takes everything from GeV**-3 to fm**3, and then it is
c multiplied by the appropriate conversion factor

c convert to fm**2 MeV**-1
      wqfmmev = wq/hbarcmev*dimfac

      wq1fmmev = wqpol1/hbarcmev*dimfac
      wq0fmmev = wqpol0/hbarcmev*dimfac
      wqm1fmmev = wqpolm1/hbarcmev*dimfac

c convert to nano barn MeV**-1

      wqnbmev = wq*1.d07/hbarcmev*dimfac

      wq1nbmev = wqpol1*1.d07/hbarcmev*dimfac
      wq0nbmev = wqpol0*1.d07/hbarcmev*dimfac
      wqm1nbmev = wqpolm1*1.d07/hbarcmev*dimfac

c tensor analysing power

      tenan = (wqpol1+wqpolm1-2.d0*wqpol0)/(wqpol1+wqpol0+wqpolm1)


c calculate the R_K according to Dmitrasinovic Gross:

      ros2l = wrlsn*ftwf*facws*8.d0*pi*pi*pi
      ros2t = dreal(cwrtsn)*ftwf*facws*8.d0*pi*pi*pi
      ros2tt = dreal(cwrttsn)*ftwf*facws*8.d0*pi*pi*pi
      ros2tl = dreal(cwrtlsn)*ftwf*facws*8.d0*pi*pi*pi

      ros2l = -ros2l/1.d03*bigqsq/(q*q)
      ros2t = ros2t/1.d03
      ros2tl = ros2tl/1.d03*dsqrt(-bigqsq)/q
      ros2tt = ros2tt/1.d03
      

      eprot = dsqrt(amp*amp+protgev*protgev)
      eneut = dsqrt(amn*amn+pmgev*pmgev)

c conversion to Wally's responses

      amass = 0.5d0*(amp+amn)
      winvmass = dsqrt((omega+amd)**2-q*q) 

      facws = amass**2/(2.d0*pi*pi)*amd/winvmass

c units of Rtilde: GeV**-1, as my responses have GeV**-3

c note that "bigqsq" is defined according to Bill, i.e.
c it is negative

      rltilde = -rl*facws*bigqsq/(q*q)
      rttilde = rt*facws
      rtttilde = rtt*facws/dcos(2.d0*protphi)
      rtltilde = rtl*facws*dsqrt(-bigqsq)/(q*dcos(protphi))

c conversion of the Rtilde from GeV**-1 to MeV**-1

      rltilde = rltilde/1.d03
      rttilde = rttilde/1.d03
      rtttilde = rtttilde/1.d03
      rtltilde = rtltilde/1.d03


c asymmetry aphi

      aphi = vtl*rtl/dcos(protphi)/(vl*rl+vt*rt+vtt*rtt)


c conversion to cm angles

      a = dsqrt(winvmass**2+q*q)/winvmass
      b = -q/winvmass

      eprotlab = dsqrt(amp*amp+protgevsq)

      protcmperp = protgev*dsin(proth)
      protcmz = eprotlab*b+protgev*dcos(proth)*a

      protcm = dsqrt(protcmperp**2+protcmz**2)


c thetapncm is really thetapcm

      thetapncm = dacos(protcmz/protcm)
      
      pthedeg = pthe/pi*180.d0

      thetapncmdeg = thetapncm/pi*180.d0


c check of kinematics:

      pmgevx = pmgev*dsin(pthe)*dcos(phi)
      pmgevy = pmgev*dsin(pthe)*dsin(phi)
      pmgevz = pmgev*dcos(pthe)

      protgevx = protgev*dsin(proth)*dcos(protphi)
      protgevy = protgev*dsin(proth)*dsin(protphi)
      protgevz = protgev*dcos(proth)

      prelx = 0.5d0*(protgevx-pmgevx)
      prely = 0.5d0*(protgevy-pmgevy)
      prelz = 0.5d0*(protgevz-pmgevz)

      prellab = dsqrt(prelx**2+prely**2+prelz**2)

c jacobian for the transformation from lab to c.m.
c from Jiri: d Omega_p^L / d Omega_kappa^cm

      appaj = dsqrt(.25d0*winvmass*winvmass-amp*amp)
      eflab = dsqrt(winvmass*winvmass+q*q)

      ajac = (appaj/protgev)**3*eflab/winvmass
      ajac = ajac*(1.d0+q*winvmass*dcos(thetapncm)/(2.d0*appaj*eflab))

      wqnbmevcm = wqnbmev*abs(ajac)
      wqos2nbmevcm = wqos2nbmev*abs(ajac)

c convert to microbarn instead of nanobarn
      wqmicbmevcm = wqnbmevcm*1.d-03


c conversion to Arenhoevel's f_xx in the c.m. system

      qcm = b*omega+a*q

      eneulab = dsqrt(amn*amn+pmgevsq)
      
      pmcmperp = pmgev*psint
      pmcmz = eneulab*b+pmgev*pcost*a

      pmcm = dsqrt(pmcmperp**2+pmcmz**2)

c protcm and pmcm are just the moduli
c actually appacm = appaj
      appacm = 0.5d0*(protcm+pmcm)


      facarrd = 3.d0*amass**2*alpha*appacm*hbarc/
     ^     (2.d0*pi*winvmass)

      root = dsqrt(qcm**2/q**2)


      farcml = rl*facarrd*root**2
      
      farcmt = rt*facarrd

      farcmtl = -rtl/dcos(protphi)*facarrd*root

      farcmtt = rtt/dcos(2.d0*protphi)*facarrd

c above, it is ok to divide out the LAB phi-anlge of the proton,
c as this is what is included in the R_K, which are the Raskin-Donnelly
c LAB system responses.
      


      return
      end
c-------------------------------------------------------
      function cy1(m,t,p)
c-------------------------------------------------------

      implicit real*8 (a-b,d-h,o-z)
      implicit complex*16 (c)
      implicit integer (i-n)

      ci = (0.d0,1.d0)
      pi = dacos(-1.d0)

      cexp = exp(ci*p*dfloat(m))


      if(m.eq.1) then
         arg = -0.5d0*dsqrt(1.5d0/pi)*dsin(t)
      elseif(m.eq.-1) then
         arg = 0.5d0*dsqrt(1.5d0/pi)*dsin(t)
      elseif(m.eq.0) then
         arg = 0.5d0*dsqrt(3.d0/pi)*dcos(t)
      else
         write(6,*) 'error in function cy1'
         write(6,*) 'm equals',m
         stop
      endif

      cy1 = arg*cexp

      return
      end
c-------------------------------------------------------
      function cy2(m,t,p)
c-------------------------------------------------------

      implicit real*8 (a-b,d-h,o-z)
      implicit complex*16 (c)
      implicit integer (i-n)

      ci = (0.d0,1.d0)
      pi = dacos(-1.d0)

      cexp = exp(ci*p*dfloat(m))


      if(m.eq.2.or.m.eq.-2) then 
         arg = 0.25d0*dsqrt(7.5d0/pi)*dsin(t)*dsin(t)
      elseif(m.eq.1) then
         arg = -0.5d0*dsqrt(7.5d0/pi)*dcos(t)*dsin(t)
      elseif(m.eq.-1) then
         arg = 0.5d0*dsqrt(7.5d0/pi)*dcos(t)*dsin(t)
      elseif(m.eq.0) then
         arg = 0.25d0*dsqrt(5.d0/pi)*
     ^         (2*dcos(t)*dcos(t)-dsin(t)*dsin(t))
      else
         write(6,*) 'error in function cy2'
         write(6,*) 'm equals',m
         stop
      endif

      cy2 = arg*cexp

      return
      end
c-------------------------------------------------------
      function cy3(m,t,p)
c-------------------------------------------------------

      implicit real*8 (a-b,d-h,o-z)
      implicit complex*16 (c)
      implicit integer (i-n)

      ci = (0.d0,1.d0)
      pi = dacos(-1.d0)

      cexp = exp(ci*p*dfloat(m))

      if(m.eq.3) then
         arg = -dsqrt(35.d0/pi)/8.d0*dsin(t)*dsin(t)*dsin(t)
      elseif(m.eq.-3) then
          arg = dsqrt(35.d0/pi)/8.d0*dsin(t)*dsin(t)*dsin(t)
      elseif(m.eq.2.or.m.eq.-2) then 
         arg = 0.25d0*dsqrt(52.5d0/pi)*dcos(t)*dsin(t)*dsin(t)
      elseif(m.eq.1) then
         arg = -dsqrt(21.d0/pi)/8.d0*
     ^   (4.d0*dcos(t)*dcos(t)*dsin(t)-dsin(t)*dsin(t)*dsin(t))
      elseif(m.eq.-1) then
         arg = dsqrt(21.d0/pi)/8.d0*
     ^   (4.d0*dcos(t)*dcos(t)*dsin(t)-dsin(t)*dsin(t)*dsin(t))
      elseif(m.eq.0) then
         arg = 0.25d0*dsqrt(7.d0/pi)*
     ^   (2*dcos(t)*dcos(t)*dcos(t)
     ^          -3.d0*dsin(t)*dsin(t)*dcos(t))
      else
         write(6,*) 'error in function cy3'
         write(6,*) 'm equals',m
         stop
      endif

      cy3 = arg*cexp

      return
      end
c-------------------------------------------------------
      function cy4(m,t,p)
c-------------------------------------------------------

      implicit real*8 (a-b,d-h,o-z)
      implicit complex*16 (c)
      implicit integer (i-n)

      ci = (0.d0,1.d0)
      pi = dacos(-1.d0)

      cexp = exp(ci*p*dfloat(m))

      ds = dsin(t)
      ds2 = ds*ds
      ds3 = ds2*ds
      ds4 = ds3*ds
      dc = dcos(t)
      dc2 = dc*dc
      dc3 = dc2*dc
      dc4 = dc3*dc

      if(m.eq.4.or.m.eq.-4) then
         arg = -dsqrt(630.d0/pi)/32.d0*ds4
      elseif(m.eq.3) then
          arg = -3.d0*dsqrt(35.d0/pi)/8.d0*ds3*dc
      elseif(m.eq.-3) then
          arg = 3.d0*dsqrt(35.d0/pi)/8.d0*ds3*dc
      elseif(m.eq.2.or.m.eq.-2) then 
         arg = 3.d0/8.d0*dsqrt(2.5d0/pi)*(6.d0*ds2*dc2-ds4)
      elseif(m.eq.1) then
         arg = 3.d0/8.d0*dsqrt(5.d0/pi)
     ^          *(3.d0*dc*ds3-4.d0*dc3*ds)
      elseif(m.eq.-1) then
         arg = -3.d0/8.d0*dsqrt(5.d0/pi)
     ^          *(3.d0*dc*ds3-4.d0*dc3*ds)
      elseif(m.eq.0) then
         arg = 3.d0/16.d0*dsqrt(1.d0/pi)
     ^          *(3.d0*ds4-24.d0*dc2*ds2+2.d0*dc4)
      else
         write(6,*) 'error in function cy4'
         write(6,*) 'm equals',m
         stop
      endif

      cy4 = arg*cexp

      return
      end
c------------------------------------------------------------
      function bes0(x)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      if(x.eq.0.d0) then
         bes0 = 1.d0
      else
         bes0 = dsin(x)/x
      endif

      return
      end
c------------------------------------------------------------
      function bes1(x)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      if(x.eq.0.d0) then
         bes1 = 0.d0
      else
         bes1 = dsin(x)/(x*x)-dcos(x)/x
      endif

      return
      end
c------------------------------------------------------------
      function bes2(x)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      x2 = x*x
      x3 = x2*x

      if(x.eq.0.d0) then
         bes2 = 0.d0
      else
         bes2 = dsin(x)*(3.d0/x3-1.d0/x)-dcos(x)*3.d0/x2
      endif

      return
      end
c------------------------------------------------------------
      function bes3(x)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      x2 = x*x
      x3 = x2*x
      x4 = x3*x

      if(x.eq.0.d0) then
         bes3 = 0.d0
      else
         bes3 = dsin(x)*(15.d0/x4-6.d0/x2)-
     ^        dcos(x)*(15.d0/x3-1.d0/x)
      endif

      return
      end
c------------------------------------------------------------
      function bes4(x)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      x2 = x*x
      x3 = x2*x
      x4 = x3*x
      x5 = x4*x

      if(x.eq.0.d0) then
         bes4 = 0.d0
      else
         bes4 = dsin(x)*(105.d0/x5-45.d0/x3+1.d0/x)+
     ^        dcos(x)*(10.d0/x2-105.d0/x4)
      endif

      return
      end
c------------------------------------------------------------
      subroutine momwave(p,um,wm)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      parameter(nmax=11)
      parameter(alpha=2.31609d-1,am0=9.d-1)
      dimension c(nmax),d(nmax),am(nmax)
 
      do j = 1,nmax
         am(j) = alpha+(j-1)*am0
      enddo
 
 
c coefficienti full model
      c(1) =  0.90457337d0
      c(2) = -0.35058661d0
      c(3) = -0.17635927d0
      c(4) = -0.10418261d2
      c(5) =  0.45089439d2
      c(6) = -0.14861947d3
      c(7) =  0.31779642d3
      c(8) = -0.37496518d3
      c(9) =  0.22560032d3
      c(10)= -0.54858290d2
 
      c(11) = 0.d0
      mx1 = nmax-1
      do k=1,mx1
         c(11) = c(11)-c(k)
      enddo
 
c coefficienti full model
      d(1) =  0.24133026d-1
      d(2) = -0.64430531d0
      d(3) =  0.51093352d0
      d(4) = -0.54419065d1
      d(5) =  0.15872034d2
      d(6) = -0.14742981d2
      d(7) =  0.44956539d1
      d(8) = -0.71152863d-1
 
      d(9) = 0.d0
      d(10) = 0.d0
      d(11) = 0.d0
      n3 = nmax-3
      sp2 = 0.d0
      sm = 0.d0
      sm2 = 0.d0
      do lp = 1,n3
         sp2 = sp2+d(lp)/am(lp)**2
         sm = sm+d(lp)
         sm2 = sm2+d(lp)*am(lp)**2
      enddo
 
      a = am(11)**2
      b = am(10)**2
      cc = am(9)**2
      d(9) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(9)**2
      b = am(11)**2
      cc = am(10)**2
      d(10) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(10)**2
      b = am(9)**2
      cc = am(11)**2
      d(11) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      ssp2 = 0.d0
      ssm = 0.d0
      ssm2 = 0.d0
      do lp = 1,nmax
         ssp2 = ssp2+d(lp)/am(lp)**2
         ssm = ssm+d(lp)
         ssm2 = ssm2+d(lp)*am(lp)**2
      enddo
c     print *,'*',ssp2,ssm,ssm2
 
 
      um = 0.d0
      wm = 0.d0

      do nn = 1,nmax
         um = um+c(nn)/(p*p+am(nn)*am(nn))
         wm = wm+d(nn)/(p*p+am(nn)*am(nn))
      enddo

      return
      end
c------------------------------------------------------------
      subroutine pamomwave(p,um,wm)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      parameter(nmax=13)
      parameter(alpha=2.3162461d-1,am0=1.d0)
      dimension c(nmax),d(nmax),am(nmax)
 
      do 111, j = 1,nmax
      am(j) = alpha+(j-1)*am0
111   continue
 

c coefficients paris potential (plb101 (1981) 139)
 
      c(1) =  0.88688076d0
      c(2) = -0.34717093d0
      c(3) = -0.30502380d1
      c(4) =  0.56207766d2
      c(5) = -0.74957334d3
      c(6) =  0.53365279d4
      c(7) = -0.22706863d5
      c(8) =  0.60434469d5
      c(9) = -0.10292058d6
      c(10)=  0.11223357d6
      c(11)= -0.75925226d5
      c(12)=  0.29059715d5
      c(13) = 0.d0
      mx1 = nmax-1
      do 112, k=1,mx1
      c(13) = c(13)-c(k)
112   continue
 



c coefficients paris potential (plb101 (1981) 139)
 
      d(1) =  0.23135193d-1
      d(2) = -0.85604572d0
      d(3) =  0.56068193d1
      d(4) = -0.69462922d2
      d(5) =  0.41631118d3
      d(6) = -0.12546621d4
      d(7) =  0.12387830d4
      d(8) =  0.33739172d4
      d(9) = -0.13041151d5
      d(10) = 0.19512524d5
      d(11) = 0.d0
      d(12) = 0.d0
      d(13) = 0.d0
 

      n3 = nmax-3
      sp2 = 0.d0
      sm = 0.d0
      sm2 = 0.d0
      do 113, lp = 1,n3
      sp2 = sp2+d(lp)/am(lp)**2
      sm = sm+d(lp)
      sm2 = sm2+d(lp)*am(lp)**2
113   continue
 
      a = am(13)**2
      b = am(12)**2
      cc = am(11)**2
      d(11) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(11)**2
      b = am(13)**2
      cc = am(12)**2
      d(12) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(12)**2
      b = am(11)**2
      cc = am(13)**2
      d(13) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 

 
      ssp2 = 0.d0
      ssm = 0.d0
      ssm2 = 0.d0
      do lp = 1,nmax
         ssp2 = ssp2+d(lp)/am(lp)**2
         ssm = ssm+d(lp)
         ssm2 = ssm2+d(lp)*am(lp)**2
      enddo
c     print *,'*',ssp2,ssm,ssm2
 
 
      um = 0.d0
      wm = 0.d0

      do nn = 1,nmax
         um = um+c(nn)/(p*p+am(nn)*am(nn))
         wm = wm+d(nn)/(p*p+am(nn)*am(nn))
      enddo

      return
      end

c------------------------------------------------------------
      subroutine av18momwave(p,um,wm)
c------------------------------------------------------------
c Rocco's parametrization of the Argonne V18 potential

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      parameter(nmax=12)
      dimension c(nmax),d(nmax),am(nmax)

      pi = dacos(-1.d0)
 
c mass parameters

      am(1) = 0.232500d+00
      am(2) = 0.500000d+00
      am(3) = 0.800000d+00
      am(4) = 0.120000d+01
      am(5) = 0.160000d+01
      am(6) = 0.200000d+01
      am(7) = 0.400000d+01
      am(8) = 0.600000d+01
      am(9) = 0.100000d+02
      am(10) = 0.140000d+02
      am(11) = 0.180000d+02
      am(12) = 0.220000d+02
 
c s wave coefficients

      c(1) =  0.105252223d+02
      c(2) =  0.124352529d+02
      c(3) = -0.687541641d+02
      c(4) =  0.239111042d+03
      c(5) = -0.441014422d+03
      c(6) =  0.300140328d+03
      c(7) = -0.230639939d+03
      c(8) =  0.409671540d+03
      c(9) = -0.733453611d+03
      c(10)=  0.123506081d+04
      c(11)= -0.120520606d+04
 
      c(12) = 0.d0
      mx1 = nmax-1
      do k=1,mx1
         c(12) = c(12)-c(k)
      enddo
 
c d wave coefficients
      d(1) =  0.280995496d+00
      d(2) =  0.334117629d-01
      d(3) = -0.727192237d+00
      d(4) = -0.302809607d+01
      d(5) = -0.903824982d+01
      d(6) =  0.496045967d+01
      d(7) = -0.271985613d+02
      d(8) =  0.125334598d+03
      d(9) = -0.346742235d+03


      d(10) = 0.d0
      d(11) = 0.d0
      d(12) = 0.d0

      n3 = nmax-3
      sp2 = 0.d0
      sm = 0.d0
      sm2 = 0.d0
      do lp = 1,n3
         sp2 = sp2+d(lp)/am(lp)**2
         sm = sm+d(lp)
         sm2 = sm2+d(lp)*am(lp)**2
      enddo
 
      a = am(12)**2
      b = am(11)**2
      cc = am(10)**2
      d(10) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(10)**2
      b = am(12)**2
      cc = am(11)**2
      d(11) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(11)**2
      b = am(10)**2
      cc = am(12)**2
      d(12) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
 
      um = 0.d0
      wm = 0.d0

      do nn = 1,nmax
         um = um+c(nn)/(p*p+am(nn)*am(nn))
         wm = wm+d(nn)/(p*p+am(nn)*am(nn))
      enddo

      um = um/(4.d0*pi)
      wm = wm/(4.d0*pi)


      return
      end


c------------------------------------------------------------
      subroutine cdbmomwave(p,um,wm)
c------------------------------------------------------------
c Rocco's parametrization of the CD Bonn potential (Machleidt)

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      parameter(nmax=12)
      dimension c(nmax),d(nmax),am(nmax)
 
      pi = dacos(-1.d0)

c mass parameters

      am(1) = 0.231600d+00
      am(2) = 0.500000d+00
      am(3) = 0.800000d+00
      am(4) = 0.120000d+01
      am(5) = 0.160000d+01
      am(6) = 0.200000d+01
      am(7) = 0.400000d+01
      am(8) = 0.600000d+01
      am(9) = 0.100000d+02
      am(10) = 0.140000d+02
      am(11) = 0.180000d+02
      am(12) = 0.220000d+02
 
c s wave coefficients

      c(1) = -0.112995238d+02
      c(2) =  0.342631587d+01
      c(3) = -0.185796560d+02
      c(4) =  0.697020450d+02
      c(5) = -0.126196725d+03
      c(6) =  0.109629157d+03
      c(7) = -0.216309234d+02
      c(8) = -0.505421595d+02
      c(9) =  0.188861755d+03
      c(10)= -0.417688387d+03
      c(11)=  0.496674109d+03
 
      c(12) = 0.d0

      do k=1,nmax-1
         c(12) = c(12)-c(k)
      enddo
 
c d wave coefficients
      d(1) =  -0.286623664d+00
      d(2) =   0.195863407d+00
      d(3) =  -0.196866168d+01
      d(4) =   0.208599899d+02
      d(5) =  -0.388561750d+02
      d(6) =   0.426753739d+02
      d(7) =  -0.457983138d+02
      d(8) =  -0.535081902d+00
      d(9) =   0.159713978d+03



      d(10) = 0.d0
      d(11) = 0.d0
      d(12) = 0.d0

      n3 = nmax-3
      sp2 = 0.d0
      sm = 0.d0
      sm2 = 0.d0
      do lp = 1,n3
         sp2 = sp2+d(lp)/am(lp)**2
         sm = sm+d(lp)
         sm2 = sm2+d(lp)*am(lp)**2
      enddo
 
      a = am(12)**2
      b = am(11)**2
      cc = am(10)**2
      d(10) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(10)**2
      b = am(12)**2
      cc = am(11)**2
      d(11) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(11)**2
      b = am(10)**2
      cc = am(12)**2
      d(12) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 

 
      um = 0.d0
      wm = 0.d0

      do nn = 1,nmax
         um = um+c(nn)/(p*p+am(nn)*am(nn))
         wm = wm+d(nn)/(p*p+am(nn)*am(nn))
      enddo

c take care of normalization and phase:

      um = -um
      wm = -wm

      um = um/(4.d0*pi)
      wm = wm/(4.d0*pi)

      return
      end
c------------------------------------------------------------
      subroutine wallymomwave(p,um,wm)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      parameter(nmax=12)
      dimension c(nmax),d(nmax),am(nmax)

      pi = dacos(-1.d0)
 
c mass parameters

      am(1) = 0.2316200d+00
      am(2) = 0.500000d+00
      am(3) = 0.800000d+00
      am(4) = 0.120000d+01
      am(5) = 0.160000d+01
      am(6) = 0.200000d+01
      am(7) = 0.400000d+01
      am(8) = 0.600000d+01
      am(9) = 0.100000d+02
      am(10) = 0.140000d+02
      am(11) = 0.180000d+02
      am(12) = 0.220000d+02
 
c s wave coefficients
      c(1) = 0.113835937E+02
      c(2) = -.818633452E+01
      c(3) = 0.596621395E+02
      c(4) = -.258482160E+03
      c(5) = 0.516421537E+03
      c(6) = -.402294428E+03
      c(7) = 0.153365492E+03
      c(8) = -.115551489E+03
      c(9) = 0.130772040E+03
      c(10)= -.218322352E+03
      c(11)= 0.222370915E+03


      c(12) = 0.d0
      mx1 = nmax-1
      do k=1,mx1
         c(12) = c(12)-c(k)
      enddo

      c(12)= -.937069211E+02


c d wave coefficients
      d(1) = -.304115815E+00
      d(2) = 0.109108134E+01
      d(3) = -.646883058E+01
      d(4) = 0.300888843E+02
      d(5) = -.488288588E+02
      d(6) = 0.484448585E+02
      d(7) = -.651042737E+02
      d(8) = 0.759630709E+02
      d(9) = -.112034936E+03

      d(10) = 0.d0
      d(11) = 0.d0
      d(12) = 0.d0

      n3 = nmax-3
      sp2 = 0.d0
      sm = 0.d0
      sm2 = 0.d0
      do lp = 1,n3
         sp2 = sp2+d(lp)/am(lp)**2
         sm = sm+d(lp)
         sm2 = sm2+d(lp)*am(lp)**2
      enddo
 
      a = am(12)**2
      b = am(11)**2
      cc = am(10)**2
      d(10) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(10)**2
      b = am(12)**2
      cc = am(11)**2
      d(11) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      a = am(11)**2
      b = am(10)**2
      cc = am(12)**2
      d(12) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
 
      ssp2 = 0.d0
      ssm = 0.d0
      ssm2 = 0.d0
      do lp = 1,nmax
         ssp2 = ssp2+d(lp)/am(lp)**2
         ssm = ssm+d(lp)
         ssm2 = ssm2+d(lp)*am(lp)**2
      enddo


      d(10) = 0.202201330E+03
      d(11) = -.215603795E+03
      d(12) = 0.931991144E+02



      um = 0.d0
      wm = 0.d0

      do nn = 1,nmax
         um = um+c(nn)/(p*p+am(nn)*am(nn))
         wm = wm+d(nn)/(p*p+am(nn)*am(nn))
      enddo


      anorm = 1.0047863577500d0

      um = um/(4.d0*pi*dsqrt(anorm))
      wm = wm/(4.d0*pi*dsqrt(anorm))


      return
      end
c------------------------------------------------------------
      function gc(j1,m1,j2,m2,j,m)
c------------------------------------------------------------
c this routine should be able to calculate Clebsch-Gordan Coefficients
c using formula (3.6.10) from Edmonds for angular momenta that
c are not too large


      implicit real*8 (a-b,d-h,o-z)
      implicit integer (i-n)
      implicit complex*16 (c)

      parameter (imax=20)

      dimension fac(0:imax)

c choice for 0! taken as otherwise the C-G coefficient for
c j=m would always vanish


      fac(0) = 1.d0
      fac(1) = 1.d0
      fac(2) = 2.d0
      fac(3) = 6.d0
      fac(4) = 24.d0
      fac(5) = 120.d0
      fac(6) = 720.d0

      if(imax.gt.6) then
         ndiff = imax-6
         aux = fac(6)
         do n=1,ndiff
            aux = aux*dfloat(6+n)
            fac(6+n) = aux
         enddo
      endif


      if(m1+m2.ne.m) then
         vc = 0.d0
      elseif(abs(m1).gt.j1) then
         vc = 0.d0
      elseif(abs(m2).gt.j2) then
         vc = 0.d0
      elseif(abs(m).gt.j) then
         vc = 0.d0
      elseif(j.gt.(j1+j2)) then
         vc = 0.d0
      elseif(j.lt.abs(j1-j2)) then
         vc = 0.d0
      elseif(j1.lt.0) then
         vc = 0.d0
      elseif(j2.lt.0) then
         vc = 0.d0
      elseif(j.lt.0) then
         vc = 0.d0
      else
         anum1 = dfloat(2*j+1)*fac(j1+j2-j)*fac(j1-m1)*fac(j2-m2)
         anum = anum1*fac(j+m)*fac(j-m)
         denom1 = fac(j1+j2+j+1)*fac(j1-j2+j)*fac(-j1+j2+j)
         denom = denom1*fac(j1+m1)*fac(j2+m2)         
         ffac = dsqrt(anum/denom)

         maxs = min(j1-m1,j-m)
         mins = max(0,j-j2-m1)
         sum = 0.d0
            do k=mins,maxs
               a = fac(j1+m1+k)*fac(j2+j-m1-k)
               bb = fac(k)*fac(j1-m1-k)*fac(j-m-k)
               b = bb*fac(j2-j+m1+k)
               sum = sum+dfloat((-1)**(k+j1-m1))*a/b
            enddo
            vc = ffac*sum

      endif

      gc = vc

      return
      end
c  This file contains:
c    sachs_sj,dircff_sj,ftog_sj,ftofp_sj,gtofp_sj,hohler_sj,
c    GariK_sj,ijl_sj,dipole_sj,GariKrumpelmann
c
c_______________________________________________________________________
c
c*****sachs_sj
c
      subroutine sachs_sj(q2,g,iff)
c
c     This subroutine returns Sachs form factors for a variety of
c     parameterizations.
c
c
c     INPUT
c             Q2.............negative of square of four-momentum transfer
c                            in Gev**2. Therefore Q2>0 for spacelike
c                            momentum transfers.
c             iff=1..........Hohler 8.2 parameterization with renormaliztion
c                            according to Friar.
c                =2..........Gari and Krumpelmann
c                =3..........Iachello, Jackson and Lande
c                =4..........dipole parameterization, Galster
c                =5..........dipole parameterization, F_1 neutron=0
c                =6..........dipole parameteriztion by Cioffi degli Atti
c                =7..........dipole G_E neutron=0.0
c                =8..........Gari and Krumpelmann, model 1
c                =9..........Gari and Krumpelmann, model 2
c                =10.........Gari and Krumpelmann, model 3
c                =11.........Mergell, Meissner and Drechsel
c                =12.........MMDA Alternate three-pole parameterization
c                =13.........Ulf-G Meissner, V. Mull, J. Speth and J. W. Van Orden,
c                            Phys. Lett. B 408, 381 (1997).
c                =14.........VO, Five pole parameterization with all parameters floating,
c                            best fit to old data.
c                =15.........Five pole parameterization that fits new preliminary CEBAF
c                            G_E^p data 
c
c     OUTPUT
c             g(iso,itype)...form factors
c
c             iso=0..........proton
c                 1..........neutron
c             itype=0........electric form factor
c                  =1........magnetic form factor
c
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c
      dimension g(0:1,0:1),f(2,2)
c
      if(iff.eq.1) then
        call hohler_sj(q2,f)
        call ftog_sj(q2,f,g)
      else if(iff.eq.2) then
        call garik_sj(q2,f)
        call ftog_sj(q2,f,g)
      else if(iff.eq.3) then
        call ijl_sj(q2,f)
        call ftog_sj(q2,f,g)
      else if(iff.ge.4.and.iff.le.7) then
        iffm3=iff-3
        call dipole_sj(q2,g,iffm3)
      else if(iff.ge.8.and.iff.le.10) then
        iffm7=iff-7
        call GariKrumpelmann(q2,iffm7,f)
        call ftog_sj(q2,f,g)
      else if(iff.ge.11.and.iff.le.15) then
	  iffm10=iff-10
        call mmd(q2,iffm10,f)
        call ftog_sj(q2,f,g)
      else
        print *,' unknown nucleon form factor type, iff=',iff
      endif
c
      return
      end
c
c_______________________________________________________________________
c
c*****dircff_sj
c
      subroutine dircff_sj(q2,fp,iff)
c
c     This subroutine returns Dirac and Pauli form factors for a variety of
c     parameterizations.
c
c
c     INPUT
c             Q2.............negative of square of four-momentum transfer
c                            in Gev**2. Therefore Q2>0 for spacelike
c                            momentum transfers.
c             iff=-2.........F_2p=1 all others =0
c                =-1.........F_1p=1 all others =0
c                =0..........all =0
c                =1..........Hohler 8.2 parameterization with renormaliztion
c                            according to Friar.
c                =2..........Gari and Krumpelmann
c                =3..........Iachello, Jackson and Lande
c                =4..........dipole parameterization
c                =5..........dipole parameterization, F_1 neutron=0
c                =6..........dipole parameteriztion by Cioffi degli Atti
c                =7..........dipole G_E neutron=0.0
c                =8..........Gari and Krumpelmann, model 1
c                =9..........Gari and Krumpelmann, model 2
c                =10.........Gari and Krumpelmann, model 3
c                =11.........Mergell, Meissner and Drechsel
c                =12.........MMDA Alternate three-pole parameterization
c                =13.........Ulf-G Meissner, V. Mull, J. Speth and J. W. Van Orden,
c                            Phys. Lett. B 408, 381 (1997).
c                =14.........VO, Five pole parameterization with all parameters floating,
c                            best fit to old data.
c                =15.........Five pole parameterization that fits new preliminary CEBAF
c                            G_E^p data 
c                
c
c     OUTPUT
c             fp(iso,itype)..form factors
c
c             iso=0..........proton
c                 1..........neutron
c             itype=1........Dirac form factor
c                  =2........Pauli form factor
c
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c
      dimension g(0:1,0:1),f(2,2),fp(0:1,2)
c
      if(iff.eq.-2) then
        fp(0,1)=0.
        fp(0,2)=1.
        fp(1,1)=0.
        fp(1,2)=0.
      else if(iff.eq.-1) then
        fp(0,1)=1.
        fp(0,2)=0.
        fp(1,1)=0.
        fp(1,2)=0.
      else if(iff.eq.0) then
        fp(0,1)=0.
        fp(0,2)=0.
        fp(1,1)=0.
        fp(1,2)=0.
      else if(iff.eq.1) then
        call hohler_sj(q2,f)
        call ftofp_sj(f,fp)
      else if(iff.eq.2) then
        call garik_sj(q2,f)
        call ftofp_sj(f,fp)
      else if(iff.eq.3) then
        call ijl_sj(q2,f)
        call ftofp_sj(f,fp)
      else if(iff.ge.4.and.iff.le.7) then
        iffm3=iff-3
        call dipole_sj(q2,g,iffm3)
        call gtofp_sj(q2,g,fp)
      else if(iff.ge.8.and.iff.le.10) then
        iffm7=iff-7
        call GariKrumpelmann(q2,iffm7,f)
        call ftofp_sj(f,fp)
      else if(iff.ge.11.and.iff.le.15) then
	  iffm10=iff-10
        call mmd(q2,iffm10,f)
        call ftofp_sj(f,fp)
      else
        print *,' unknown nucleon form factor type, iff=',iff
      endif
c
      return
      end

c
c_______________________________________________________________________
c
c*****ftog_sj
c
      subroutine ftog_sj(q2,f,g)
c
c     This subroutine translates the representation of Dirac
c     form factors to Sachs form factors.
c
      implicit real*8(a-h,o-z)
      implicit integer (i-n)
c
      real*8 mn
c
      dimension f(2,2),g(0:1,0:1),fp(0:1,2)
c
      data mn/.9389/
c
      tau=q2/(4.*mn**2)
c
      call ftofp_sj(f,fp)
c
      do 100 i=0,1
        g(i,0)=fp(i,1)-tau*fp(i,2)
        g(i,1)=fp(i,1)+fp(i,2)
 100  continue
c 
      return
      end
c
c_______________________________________________________________________
c
c*****ftofp_sj
c
      subroutine ftofp_sj(f,fp)
c
c     This subroutine translates isoscalar and isovector Dirac form 
c     factors to neutron and proton form factors.
c
      implicit real*8(a-h,o-z)
      implicit integer (i-n)
c
      dimension f(2,2),fp(0:1,2)
c
      do 100 i=1,2
        fp(0,i)=0.5*(f(2,i)+f(1,i))
        fp(1,i)=0.5*(f(2,i)-f(1,i))
 100  continue
c
      return
      end
c
c_______________________________________________________________________
c
c*****gtofp_sj
c
      subroutine gtofp_sj(q2,g,fp)
c
c     This subroutine translates Sachs form factors to Dirac form
c     factors.

      implicit real*8(a-h,o-z)
      implicit integer (i-n)
c
      real*8 mn
c
      dimension g(0:1,0:1),fp(0:1,2)
c
      data mn/.9389/
c
      tau=q2/(4.*mn**2)
c
      do 100 i=0,1
        fp(i,1)=(g(i,0)+tau*g(i,1))/(1.0+tau)
        fp(i,2)=(g(i,1)-g(i,0))/(1.+tau)
 100  continue
c
      return
      end
c
c_______________________________________________________________________
c
c*****hohler_sj
c
      subroutine hohler_sj(q2,f)
c
c     This subroutine calculates the electromagnetic form factors
c     using the Hohler parameterization 8.2.  Input is in GeV**2
c     and the square of the four-momentum transfer should be
c     positive.
c
c     REFERENCE: G. Hohler, et al., Nucl. Phys. B114, 505 (1976).
c
c     INPUT
c             Q2.............negative of square of four-momentum transfer
c                            in Gev**2. Therefore Q2>0 for spacelike
c                            momentum transfers.
c
c     OUTPUT
c             f(iso,itype)...form factors
c
c             iso=1..........isovector
c                 2..........isoscalar
c             itype=1........Dirac form factor F_1
c                  =2........Pauli form factor F_2
c
c     7-16-87
c            This version has been modified such that it takes into
c            account changes made by J. Friar to this parameterization
c            to give better behavior at the q2=0 point.
c
      implicit real*8(a-h,o-z)
      implicit integer (i-n)
      dimension srt(3,2),a(3,2,2),f(2,2),frov(2)
      data((srt(i,j),i=1,3),j=1,2)/1.21,2.45,2.95,0.783,1.02,1.8/
      data(((a(i,j,ij),i=1,3),j=1,2),ij=1,2)/0.05,-0.52,0.28,0.71
     #  ,-0.64,-0.13,-1.99,0.1956,0.1856,-0.11,0.13,-0.02/
c
      frov(1)=(0.955+0.09/(1.+q2/0.355)**2)/(1.+q2/0.536)/2.
      frov(2)=(5.335+0.962/(1.+q2/0.268))/(1.+q2/.603)/2.
c
      do 50 i1=1,2
        do 50 i2=1,2
          f(i1,i2)=frov(i2)*(1-(-1)**i1)/2.
          do 50 i3=1,3
50          f(i1,i2)=f(i1,i2)+a(i3,i1,i2)/(srt(i3,i1)**2+q2)
c
c                  Renormalization according to Friar
c
      f(1,1)=2.*f(1,1)*0.9954211
      f(2,1)=2.*f(2,1)*0.9930881
      f(1,2)=2.*f(1,2)*1.0052166
      f(2,2)=2.*f(2,2)*0.98891518
c
      return
      end
c
c_______________________________________________________________________
c
c*****GariK_sj
c
      subroutine garik_sj(q2,f)
c
c     This subroutine produces the form factors of Gari and Krumpelmann.
c
c
c     REFERENCE: M. Gari and W. Krumpelmann, Z. Phys. A322 (1985);
c                                            Phys. Lett. 173B, 10 (1986).
c
c     INPUT
c             Q2.............negative of square of four-momentum transfer
c                            in Gev**2. Therefore Q2>0 for spacelike
c                            momentum transfers.
c
c     OUTPUT
c             f(iso,itype)...form factors
c
c             iso=1..........isovector
c                 2..........isoscalar
c             itype=1........Dirac form factor F_1
c                  =2........Pauli form factor F_2
c
      implicit real*8(a-h,o-z)
c
      real*8 mrho,momega,kapv,kaps,kaprho,kapomg,lam1,lam2,lamqcd
c
      dimension f(2,2)
c
      data mrho/0.7661/,momega/0.784/,kapv/3.706/,kaps/-0.12/,
     &     kaprho/6.62/,kapomg/0.163/,crho/0.377/,comg/0.411/,
     &     lam1/0.795/,lam2/2.27/,lamqcd/0.29/
c
      q2tild=q2*log((lam2**2+q2)/lamqcd**2)/log((lam2/lamqcd)**2)
c
      d1=1./(1.+q2tild/lam1**2)
      d2=1./(1.+q2tild/lam2**2)
c
      f1=d1*d2
      f2=f1*d2
c
      drho=1./(1.+q2/mrho**2)
      domg=1./(1.+q2/momega**2)
c
      f(1,1)=(crho*drho+1.-crho)*f1
      f(2,1)=(comg*domg+1.-comg)*f1
      f(1,2)=(crho*kaprho*drho+kapv-crho*kaprho)*f2
      f(2,2)=(comg*kapomg*domg+kaps-comg*kapomg)*f2
c
      return
      end
c
c
c_______________________________________________________________________
c
c*****ijl_sj
c
      subroutine ijl_sj(q2,f)
c
c     This subroutine produces the form factors of Iachello, Jackson and
c     Lande.
c
c     REFERENCE: F. Iachello, A. D. Jackson and A. Lande, Phys. Lett. 43B,
c                191 (1973).
c
c     INPUT
c             Q2.............negative of square of four-momentum transfer
c                            in Gev**2. Therefore Q2>0 for spacelike
c                            momentum transfers.
c
c     OUTPUT
c             f(iso,itype)...form factors
c
c             iso=1..........isovector
c                 2..........isoscalar
c             itype=1........Dirac form factor F_1
c                  =2........Pauli form factor F_2
c
      implicit real*8(a-h,o-z)
c
      real*8 mpi,mrho,momega,mphi,kapv,kaps
c
      dimension f(2,2)
c
      data mpi/.13957/,mrho/0.765/,momega/0.784/,mphi/1.019/,
     &     gamrho/.112/,gam/0.25/,btarho/0.672/,btaomg/1.102/,
     &     btaphi/0.112/,alpphi/-0.052/,kaps/-0.120/,kapv/3.706/,
     &     pi/3.141592741012573/
c
      g=1./(1.+gam*q2)**2
c
      if(q2.gt.1.e-7) then
        alpha=2./pi*sqrt((q2+4.*mpi**2)/q2)*log((sqrt(q2+4.*mpi**2)+
     &        sqrt(q2))/(2.*mpi))
      else
        alpha=0.
      endif
c
      drho=(mrho**2+8.*gamrho*mpi/pi)/(mrho**2+q2+(4.*mpi**2+q2)*
     &gamrho*alpha/mpi)
      domg=1./(1.+q2/momega**2)
      dphi=1./(1.+q2/mphi**2)
c
      f(1,1)=g*(1.-btarho+btarho*drho)
      f(2,1)=g*(1-btaomg-btaphi+btaomg*domg+btaphi*dphi)
      f(1,2)=g*kapv*drho
      f(2,2)=g*((kaps-alpphi)*domg+alpphi*dphi)
c
      return
      end
c
c_______________________________________________________________________
c
c*****dipole_sj
c
      subroutine dipole_sj(q2,g,iff)
c
c     This subroutine produces several form factors which are expressed
c     as dipole fits to the Sachs form factors.
c
c
c     INPUT
c             Q2.............negative of square of four-momentum transfer
c                            in Gev**2. Therefore Q2>0 for spacelike
c                            momentum transfers.
c             iff=1..........Galster
c                =2..........F_1 neutron=0
c                =3..........Cioffi degli Atti
c                =4..........G_E neutron=0
c
c     OUTPUT
c             g(iso,itype)...form factors
c
c             iso=0..........proton
c                 1..........neutron
c             itype=0........electric form factor
c                  =1........magnetic form factor
c
      implicit real*8(a-h,o-z)
      implicit integer (i-n)
c
      real*8 mn,mup,mun
c
      dimension g(0:1,0:1)
c
      data mn/0.939/,mup/2.793/,mun/-1.913/
c
      tau=q2/(4.*mn**2)
c
c                          Galster 
c               S. Galster, et al., Nucl. Phys. B32, 221 (1971).
c
      if(iff.eq.1) then
        g(0,0)=1./(1.+q2/0.71)**2
        g(1,0)=-mun*tau*g(0,0)/(1.+5.6*tau)   !5.6
        g(0,1)=mup*g(0,0)
        g(1,1)=mun*g(0,0)
c
c           dipole with neutron Dirac form factor set to zero
c
      else if(iff.eq.2) then
        g(0,0)=1./(1.+q2/0.71)**2
        g(1,0)=-mun*tau*g(0,0)
        g(0,1)=mup*g(0,0)
        g(1,1)=mun*g(0,0)
c
c           dipole fit according to Cioffi degli Atti
c
      else if(iff.eq.3) then
        g(0,0)=1./(1.+q2/0.695)**2
        g(1,0)=-mun*tau*g(0,0)/(1.+4.*tau)
        g(0,1)=mup*g(0,0)
        g(1,1)=mun*g(0,0)
c
c                     G_E neutron = 0.
c
      else if(iff.eq.4) then
        g(0,0)=1./(1.+q2/0.71)**2
        g(1,0)=0.0
        g(0,1)=mup*g(0,0)
        g(1,1)=mun*g(0,0)
      endif
      return
      end
c
c_______________________________________________________________________
c
c*****GariKrumpelmann
c
      subroutine GariKrumpelmann(q2,i_version,f)
c
c     This subroutine produces the form factors of Gari and Krumpelmann.
c     This includes all versions in the second set of references
c
c
c     REFERENCE: M. Gari and W. Krumpelmann, Z. Phys. A322 (1985);
c                                            Phys. Lett. 173B, 10 (1986);
c
c                                            Phys. Lett. B274, 159 (1992);
c                                            Phys. Lett. B282, 483 (1992).
c
c     INPUT
c             Q2.............negative of square of four-momentum transfer
c                            in Gev**2. Therefore Q2>0 for spacelike
c                            momentum transfers.
c             i_version......1 to 3 corresponding to notation in 1992 
c                            papers.
c
c     OUTPUT
c             f(iso,itype)...form factors
c
c             iso=1..........isovector
c                 2..........isoscalar
c             itype=1........Dirac form factor F_1
c                  =2........Pauli form factor F_2
c
c      implicit real*8(a-h,o-z)
      implicit none
c
      real*8 m_rho,m_omega,m_phi,kappa_v,kappa_s,
     &       c_rho(3),kappa_rho(3),c_omega(3),kappa_omega(3),
     &       c_phi(3),kappa_phi(3),mu_phi(3),lambda_1(3),
     &       lambdaD_1(3),lambda_2(3),lambda_qcd(3),f(2,2)
c
      real*8 q2,q2_tilda,d1,d1d,d2,f1_rho,f1_omega,f1_D,
     &       f1_phi,f2_rho,f2_omega,f2_D,f2_phi,d_rho,
     &       d_omega,d_phi
c
      integer i_version
c
      data m_rho/0.766/,m_omega/0.784/,m_phi/1.019/,
     &     kappa_v/3.706/,kappa_s/-0.12/,
     &     c_rho/0.377,0.5927,0.5688/,
     &     kappa_rho/6.62,3.425,3.642/,
     &     c_omega/0.411,0.6212,0.5774/,
     &     kappa_omega/0.163,0.3941,0.4775/,
     &     c_phi/0.0,0.0,-0.666/,
     &     kappa_phi/0.0,0.0,-0.2378/,
     &     mu_phi/0.33,0.33,0.33/,
     &     lambda_1/0.795,0.867,0.823/,
     &     lambdaD_1/0.795,1.194,1.24/,
     &     lambda_2/2.270,2.063,1.95/,
     &     lambda_qcd/0.29,0.344,0.31/
c
      q2_tilda=q2*log((lambda_2(i_version)**2+q2)/
     &                lambda_qcd(i_version)**2)/
     &            log((lambda_2(i_version)/
     &                lambda_qcd(i_version))**2)
c
      d1=1./(1.+q2_tilda/lambda_1(i_version)**2)
      d1d=1./(1.+q2_tilda/lambdaD_1(i_version)**2)
      d2=1./(1.+q2_tilda/lambda_2(i_version)**2)
c
      f1_rho=d1*d2
      f1_omega=f1_rho
      f1_D=d1d*d2
      f1_phi=f1_rho*sqrt((q2/(lambda_1(i_version)**2+q2))**3)
c
      if(i_version.eq.1) then
        f2_rho=d2*f1_rho
        f2_D=d2*f1_D
      else 
        f2_rho=d1*f1_rho
        f2_D=d1d*f1_D
      end if
      f2_omega=f2_rho
c
      f2_phi=f2_rho*sqrt(((1.+q2/mu_phi(i_version)**2)/
     &       (1.+q2/lambda_1(i_version)**2))**3)
c
      d_rho=1./(1.+q2/m_rho**2)
      d_omega=1./(1.+q2/m_omega**2)
      d_phi=1./(1.+q2/m_phi**2)
c
      f(1,1)=c_rho(i_version)*d_rho*f1_rho
     &      +(1.0-c_rho(i_version))*f1_D
      f(2,1)=c_omega(i_version)*d_omega*f1_omega
     &      +c_phi(i_version)*d_phi*f1_phi
     &      +(1.0-c_omega(i_version))*f1_D
      f(1,2)=kappa_rho(i_version)*c_rho(i_version)*
     &        d_rho*f2_rho
     &      +(kappa_v-kappa_rho(i_version)*
     &        c_rho(i_version))*f2_D
      f(2,2)=kappa_omega(i_version)*c_omega(i_version)*
     &        d_omega*f2_omega
     &      +kappa_phi(i_version)*c_phi(i_version)*
     &        d_phi*f2_phi
     &      +(kappa_s-kappa_omega(i_version)*
     &        c_omega(i_version)
     &       -kappa_phi(i_version)*c_phi(i_version))*f2_D
c
      return
      end
c_______________________________________________________________________
c
c*****MMD
c
c      subroutine mmd(q2,f)
c
c     This subroutine produces the for factors of Mergell, Meissner
c     and Drechsel.
c
c     REFERENCE P. Mergell, Ulf-G. Meissner and D. Drechsel, Nucl. Phys.
c               A 596, 367 (1996).
c
c     INPUT
c             Q2.............negative of square of four-momentum transfer
c                            in Gev**2. Therefore Q2>0 for spacelike
c                            momentum transfers.
c
c     OUTPUT
c             f(iso,itype)...form factors
c
c             iso=1..........isovector
c                 2..........isoscalar
c             itype=1........Dirac form factor F_1
c                  =2........Pauli form factor F_2
c
c      implicit real*8(a-h,o-z)
c
c      dimension f(2,2)
c
c      fln=log((9.73+q2)/0.350)**(-2.148)
c
c      f(1,1)=2.*((13.16596+1.1245/(1+q2/0.318)**2)/(2.*(1+q2/0.550))
c     &       -38.868/(2.103+q2)+425.339/(2.723+q2)
c     &       -388.651/(2.856+q2))*fln
c      f(1,2)=2.*((73.7917+5.0212/(1+q2/0.142))/(2.*(1+q2/0.5360))
c     &       -73.495/(2.103+q2)+83.261/(2.723+q2)
c     &       -29.394/(2.856+q2))*fln
c      f(2,1)=2.*(9.457/(0.612+q2)-9.050/(1.038+q2)
c     &        -0.410/(2.560+q2))*fln
c      f(2,2)=2.*(-1.545/(0.612+q2)+1.987/(1.038+q2)
c     &        -0.436/(2.560+q2))*fln
c
c      return
c      end
c
c_______________________________________________________________________
c
c*****MMD
c
      subroutine mmd(q2,i_version,f)
c
c     This subroutine produces the for factors of Mergell, Meissner
c     and Drechsel.
c
c     REFERENCE P. Mergell, Ulf-G. Meissner and D. Drechsel, Nucl. Phys.
c               A 596, 367 (1996).
c
c     INPUT
c             Q2.............negative of square of four-momentum transfer
c                            in Gev**2. Therefore Q2>0 for spacelike
c                            momentum transfers.
c             i_version=1....MMD parameterization. This does not use the exact parameters
c                            from the reference above since they are do not reliably
c                            reproduce the fit and constraints.
c                      =2....Alternate three-pole parameterization
c                      =3....Ulf-G Meissner, V. Mull, J. Speth and J. W. Van Orden,
c                            Phys. Lett. B 408, 381 (1997).
c                      =4....Five pole parameterization with all parameters floating,
c                            best fit to old data.
c                      =5....Five pole parameterization that fits new preliminary CEBAF
c                            G_E^p data 
c
c     OUTPUT
c             f(iso,itype)...form factors
c
c             iso=1..........isovector
c                 2..........isoscalar
c             itype=1........Dirac form factor F_1
c                  =2........Pauli form factor F_2
c
      implicit none
c
      real*8 q2,fint(0:1,2),f(2,2),fln,factor          !,g(2)
c
      integer i,n,iso,iexp,i_version
c
      integer max_terms,n_terms(0:1,5)
c
      parameter(max_terms=5)
c
      real*8 c_s(max_terms,0:2,0:1,5),Q0Sqrd,LambdaSqrd(5),arho(2),
     &       brho(2),crho(2),drho(2),MaSqrd,MbSqrd,gamma
c 
      parameter(gamma=2.148)
c
      data MaSqrd/0.5/,MbSqrd/0.4/,Q0Sqrd/0.35/
      data arho/1.0317,5.7824/
      data brho/0.0875,0.3907/
      data crho/0.3176,0.1422/
      data drho/0.5496,0.5362/
	data LambdaSqrd(1)/9.733/
      data (n_terms(iso,1),iso=0,1)/3,3/
      data ((c_s(n,i,0,1),i=0,2),n=1,3)/
     &      6.115240E-01,   7.469750E-01,  -1.221656E-01,
     &      1.038361E+00,  -7.361827E-01,   1.609093E-01,
     &      2.585406E+00,  -3.998969E-02,  -3.983142E-02/
      data ((c_s(n,i,1,1),i=0,2),n=1,3)/
     &      2.102500E+00,  -3.414994E+00,  -6.544095E+00,
     &      2.734393E+00,   3.099864E+01,   7.134224E+00,
     &      2.868281E+00,  -2.809906E+01,  -2.081179E+00/
c
	data LambdaSqrd(2)/7.1549/
      data (n_terms(iso,2),iso=0,1)/3,3/
      data ((c_s(n,i,0,2),i=0,2),n=1,3)/
     &      6.115240E-01,   7.368497E-01,  -1.106676E-01,
     &      1.038361E+00,  -7.204214E-01,   1.310221E-01,
     &      4.422851E+00,  -1.104352E-01,  -2.986785E-02/
      data ((c_s(n,i,1,2),i=0,2),n=1,3)/
     &      2.102500E+00,  -4.904924E+00,  -5.681093E+00,
     &      2.734393E+00,   4.603449E+01,  -9.500462E-01,
     &      2.868281E+00,  -4.203599E+01,   5.292896E+00/
c
      data LambdaSqrd(3)/9.733/
      data (n_terms(iso,3),iso=0,1)/4,3/
      data ((c_s(n,i,0,3),i=0,2),n=1,4)/
     &      6.115240E-01,   6.773855E-01,  -9.657211E-02,
     &      1.038361E+00,  -1.930000E-02,  -4.000000E-03,
     &      1.254400E+00,  -7.987780E-01,   1.525659E-01,
     &      2.656775E+00,   1.215471E-01,  -5.340202E-02/
      data ((c_s(n,i,1,3),i=0,2),n=1,3)/
     &      2.102500E+00,  -3.607208E+00,  -6.534883E+00,
     &      2.722500E+00,   2.017498E+01,   5.946551E+00,
     &      2.958706E+00,  -1.710495E+01,  -9.027768E-01/
c
      data LambdaSqrd(4)/9.733/
      data (n_terms(iso,4),iso=0,1)/5,5/
	data ((c_s(n,i,0,4),i=0,2),n=1,5)/
     &      6.365155E-01,   7.816823E-01,   2.282614E-01,
     &      1.307334E+00,  -2.226998E-02,   9.921241E-02,
     &      1.225704E+00,  -5.698303E-01,   3.458286E+00,
     &      2.580973E+00,   1.073272E-01,  -2.756912E-01,
     &      1.095059E+00,  -3.173048E-01,  -3.514883E+00/
	data ((c_s(n,i,1,4),i=0,2),n=1,5)/
     &      2.432420E+00,   5.649976E+01,  -6.387368E+00,
     &      2.560801E+00,  -4.520188E+01,   7.035124E+00,
     &      2.072874E+00,  -1.159112E+01,  -4.658646E+00,
     &      3.038462E+00,  -8.905151E-02,   4.109003E+00,
     &      3.497939E+00,  -1.399468E-01,  -1.594393E+00/
c
      data LambdaSqrd(5)/9.733/
      data (n_terms(iso,5),iso=0,1)/5,5/
	data ((c_s(n,i,0,5),i=0,2),n=1,5)/
     &      6.780836E-01,   1.030995E+00,   2.596098E-01,
     &      1.138074E+00,  -2.295158E-01,   1.617209E+00,
     &      1.734596E+00,  -5.213922E-01,  -7.919872E-02,
     &      3.215552E+00,   3.100522E-01,  -8.486819E-02,
     &      9.577866E-01,  -5.884375E-01,  -1.716670E+00/
	data ((c_s(n,i,1,5),i=0,2),n=1,5)/
     &      1.597234E+00,  -1.573959E-01,  -7.707051E+00,
     &      1.916871E+00,   6.384136E-01,   6.882814E+00,
     &      8.032118E-01,  -8.592565E-02,   3.692524E-01,
     &      3.708225E+00,   7.088970E-02,  -6.981130E+00,
     &      4.195078E+00,  -9.613039E-01,   6.031879E+00/
c
      factor=fln(q2,LambdaSqrd(i_version),Q0Sqrd)
c
      do iso=0,1
c
        do i=1,2
c
          fint(iso,i)=0.0
c
          do n=1,n_terms(iso,i_version)
            fint(iso,i)=fint(iso,i)+c_s(n,i,iso,i_version)
     &             /fln(-c_s(n,0,iso,i_version),
     &         LambdaSqrd(i_version),Q0Sqrd)/(c_s(n,0,iso,i_version)+q2)
          end do
c
          fint(iso,i)=fint(iso,i)*factor
c
        end do
c
      end do
c
      do i=1,2
        iexp=3-i
        fint(1,i)=fint(1,i)+0.5d0*factor*(arho(i)/
     &         fln(-MaSqrd,LambdaSqrd(i_version),Q0Sqrd)
     &         +brho(i)/fln(-MbSqrd,LambdaSqrd(i_version),Q0Sqrd)/
     &         (1+q2/crho(i))**iexp)/(1+q2/drho(i))
      end do
c
      do iso=0,1
	  do i=1,2
	    f(2-iso,i)=2.d0*fint(iso,i)
	  end do
	end do
c
      return
      end
c
c_______________________________________________________________________
c
c*****fln.f
c
      function fln(Q2,LambdaSqrd,Q0Sqrd)
c
      implicit none
c
      real*8 fln,Q2
c
      integer max_terms
c
      parameter(max_terms=5)
c
      real*8 Q0Sqrd,LambdaSqrd,gamma
c
      parameter(gamma=2.148)
c      
      fln=log((LambdaSqrd+Q2)/Q0Sqrd)**(-gamma)
c
      return
      end


