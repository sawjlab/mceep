C PEEK.FOR
C AUTHOR: Pete Markowitz
C DATE:   11.02.97
C MODIFICATIONS:
C         07.07.99   P. Markowitz      Additional comments added.
C         20.03.01   P. Markowitz      Changed prescription to empirical fit.  
C         
C PURPOSE:
C This routine calculates p(e,e'K+)Lambda cross-sections from the empirical
C      fit of G.Niculescu et.al., Phys. Rev. Lett. 81, 4576 (1998)
C      and J.Cha (PhD theses, HALL C web page, 2001), in turn similar to 
C      P.Brauel et.al., Z.Phys. C, 3 p. 101 (1979).
C
C VARIABLES:
C   INPUT:
C       QMAG            !3-momentum of virtual photon           (MeV/c)
C       OMEGA           !energy of virtual photon               (MeV)
C       THETA_PQ        !angle between kaon and q               (rad)
C       THETA_CM        !angle between kaon and q in COM        (rad)
C       PHI_X           !angle between kaon and q planes        (rad)
C       pf_p_i          !momentum of kaon                       (MeV/c)
C       W               !Invariant mass                         (MeV)
C       qmu2_g          !4-momentum of virtual photon, squared  (GeV/c)^2
C       ef_p_i          !momentum of electron                   (MeV/c)
C   OUTPUT:
C       sigma_eep       !d3sigma/dEe'dOmegae'OmegaK            (fm^2/MeV/sr^2)
C ROUTINES CALLED: none
C 
      SUBROUTINE PEEK(SIGMA_EEP)

      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
      INCLUDE 'masses.cmn'
      INCLUDE 'var.cmn'
      DOUBLE PRECISION ME,MK,MKG,ML,MLG,MP,MPG,EK,ALPHA,PI,SG
      DOUBLE PRECISION PKCM,GTPR,FQ,GW,HT,IP,WG
      DOUBLE PRECISION SIGMA_EEP,DSIGMA,JACOB
      DOUBLE PRECISION PCM,BCM,GCM,EF_E_I,HC
c      INTEGER          
      HC   = 197.327053
      ME   = 0.51099906
      MK   = 493.677
      MKG  = .493677
      ML   = 1115.684
      MLG  = 1.115684
      MP   = 938.27231
      MPG  = MP/1000.D0
      WG   = W/1000.D0
      EK   = SQRT(MK*MK+PF_P_I*PF_P_I)
      ALPHA= 1.D0/137.0359895D0
      PI   = 3.141592654D0
      SG   = (W/1000.D0)*(W/1000.D0)
      EF_E_I = SQRT(PF_E_I**2 + ME**2)
C beta_cm and gamma_cm -- inputs all in MeV so outputs are dimensionless
C now determine momentum of kaon in K-Lambda cm frame, needed to properly add
C together the components of the cross-section
      PKCM = SQRT( (SG-(MKG+MLG)**2)*(SG-(MKG-MLG)**2)/4./SG )
      BCM  = QMAG/(OMEGA+MP)
      GCM  = 1.D0/SQRT(1.D0-BCM*BCM)
C CM TO LAB JACOBIAN
      JACOB =1.D0/((SIN(THETA_PQ)/SIN(THETA_CM))**2.D0*
     2 (COS(THETA_PQ)*COS(THETA_CM)+GCM*SIN(THETA_PQ)*SIN(THETA_CM)))
C virtual photon to electron beam flux conversion factor in MeV bins not GeV
      GTPR = ALPHA/2.D0/(PI**2.D0)*EF_E_I/E0_I*(SG-MPG**2.D0)/2.D0
     #                  /MPG/QMU2_G/(1.D0-EPSILON)/1000.D0
C calculate t and tmin
       TMIN=((-QMU2_G-MKG*MKG-MPG*MPG+MLG*MLG)/(2.0*WG))**2.0 
     #   -(sqrt( ((SG -QMU2_G -MPG*MPG)/(2.*WG))**2. +QMU2_G)
     #   -sqrt( ((SG+MKG*MKG -MLG*MLG)/(2.*WG))**2.-MKG*MKG))**2.
        PCM = (QMAG/1000)*MPG/WG
        T=TMIN-4.D0*(PCM/1000.D0)*PKCM/1000.*SIN(THETA_CM/2.D0)**2.D0
C
C---CALCULATE VIRTUAL PHOTOPRODUCTION CROSS SECTION IN fm^2/sr
C f(Q^2) gives dsig/dO = 3.832/(Q2+2.67)^2 in uB/sr at W=2.15, t=tmin
        FQ=3.832E-4/((QMU2_G+2.67)**2)
c g(W) gives i) phase space divided by phase space for W=2.15 
c (g(W)~pkcm/(s-mp^2)/W) plus ii) a resonance at 1.72 GeV
        GW=(PKCM/(SG-MPG*MPG)/(W/1000.))/8.4810E-2 
     #    + S2G/(SG*1.4369 +(SG-1.72)**2.D0)
c
c  h(tmin-t) falls exponentially
        HT=exp(-.2144*(TMIN-T))
c i(phi)
        IP=1-.01096*COS(PHI_X)+.0018*cos(2*PHI_X)
c
        DSIGMA=FQ*GW*HT*IP
        SIGMA_EEP=GTPR*JACOB*DSIGMA*(2+EPSILON)/3
C    include a debugging write statement but comment it out
c        write(*,*) PKCMP,PKCM, PKCM*1000.D0-PKCMP, P_CM
c        write(*,*) 'dsigma=',DSIGMA,PHI_X,IP,sigma_eep,HT
c        write(*,*) 't,tmin,ht=',T,TMIN,HT
c        write(*,*) 'g(w)=',GW,DSIGMA
c        write(*,*) '-----------------------------checked 3/22/01, pecm'
      RETURN
      END
