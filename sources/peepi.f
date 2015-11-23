C PEEPI.FOR

C AUTHOR: G.M. Huber
C DATE:   93.02.09

C MODIFICATIONS:
C         11-JAN-1995   P.E. Ulmer
C         Implementation in MCEEP now weights the cross section by
C         the pion survival probability rather than applying a test
C         for survival (i.e. a cut).
C
C         20-MAR-1996   G.M. Huber
C         Improved parameterization of cross section for large Q2, based
C         on Phys.Rev.D 17(1978)1693.
C
C         02-DEC-1996   Liming Qin
C         The cross section is transformed to the laboratory frame
C         from the center-of-mass system.
C
C PURPOSE:
C This routine calculates p(e,e'pi+)n cross-sections from a fit to the data of
C Brauel et al., Z.Phys.C. 3(1979)101.

C VARIABLES:
C   INPUT:
C       qmag            !3-momentum of virtual photon           (MeV/c)
C       omega           !energy of virtual photon               (MeV)
C       theta_pq        !angle between pi and q                 (rad)
C       pf_p_i          !momentum of pion                       (MeV/c)
C       qmu2_g          !4-momentum of virtual photon, squared  (GeV/c)^2
C       e0_i            !incident electron beam energy          (MeV)
C       tscat           !electron scattering angle              (rad)
C       pf_e_i          !momentum of scattered electron         (MeV/c)
C       epsilon         !epsilon                                (dimensionless)
C       phi_x           !angle between hadron & electron planes (rad)
C       eject_mass      !mass of the pion                       (MeV/c^2)
C   OUTPUT:
C       sigma_eep       !d3sigma/dEe'dOmegae'Omegapi            (fm^2/MeV/sr^2)

C ROUTINES CALLED:  None

      SUBROUTINE PEEPI(SIGMA_EEP)

      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z)
      INCLUDE 'masses.cmn'
      INCLUDE 'var.cmn'

      DOUBLE PRECISION MP,MPG,MN,ME,PI,ALPHA,EF_P_I,EF_E_I,T,S,SG
      DOUBLE PRECISION SIGL,SIGT,SIGLT,SIGTT,SIGMA_EEP,MRHO,FPI,FPI07
      DOUBLE PRECISION GTPR,A,PPICM,PGAMCM,GAMCM,BETCM

      ME   = 0.51099906
      MN   = 939.56563
      MP   = 938.27231
      MPG  = MP/1.E3

      ALPHA = 1./137.03598
      PI    = 3.141592654

C fits to p(e,e'pi+)n cross-sections

      EF_P_I = SQRT(PF_P_I**2 + EJECT_MASS**2)
      T = (QMAG-PF_P_I*COS(THETA_PQ))**2 + (PF_P_I*SIN(THETA_PQ))**2
     #          - (OMEGA-EF_P_I)**2
      T = T/1.E6                              !convert to (GeV//c)**2

      SIGL = 30.6*EXP(-14*ABS(T))+1.15*EXP(-0.6*ABS(T))
      SIGT = 2.82*EXP(-2*ABS(T))
      SIGLT= -0.14
      SIGTT= -2.23

C Now scale sigl by Q2*pion_formfactor

      MRHO = 0.712    !determined by empirical fit to formfactor (GeV/c2)
      IF (QMU2_G.LT.2.19) THEN
        FPI=1./(1.+QMU2_G/MRHO**2)
      ELSE
        FPI=1./(1.+2.19/MRHO**2)/QMU2_G
      ENDIF
      FPI2=FPI**2
      FPI07=1./(1.+0.70/MRHO**2)
      FPI072=FPI07**2

      SIGL=SIGL*(FPI2*QMU2_G)/(FPI072*0.7)

C Brauel scaled all his data to W=2.19 GeV, so rescale by 1/(W**2-mp**2)**2

      EF_E_I = SQRT(PF_E_I**2 + ME**2)

      S = (OMEGA+MP)**2 - QMAG**2
      SG= S/1.E6                                !convert to (GeV)**2

      SIGL = SIGL*(2.19**2-MPG**2)**2/(SG-MPG**2)**2
      SIGT = SIGT*(2.19**2-MPG**2)**2/(SG-MPG**2)**2

C Outside the range of validity of Brauel's Z.Phys. paper, use Bebek's Q2=10
C form factor paper (PRD 17 (1978) 1693).

      IF (QMU2_G.GE.2. .OR. ABS(T).GE.0.2) THEN
        R   = QMU2_G/4.2
        SIGL= (SIGL+SIGT)/(R+1)
        SIGT= SIGL*R
        SIGLT=0.
        SIGTT=0.
      ENDIF

C virtual photon to electron beam flux conversion factor

      GTPR = ALPHA/2./(PI**2)*EF_E_I/E0_I*(SG-MPG**2)/2./MPG/QMU2_G
     #                  /(1.-EPSILON)

C now determine momentum of pion in n-pi cm frame, needed to properly add
C together the components of the cross-section

      PPICM = SQRT( (S-(EJECT_MASS+MN)**2)*(S-(EJECT_MASS-MN)**2)/4./S )

      GAMCM = (OMEGA+MP)/SQRT(S)
      BETCM = QMAG/(OMEGA+MP)
      PGAMCM = (QMAG-BETCM*OMEGA)*GAMCM

      A = PGAMCM*PPICM/PI

      SIGMA_EEP = GTPR*( A*SIGT + EPSILON*A*SIGL
     #          + EPSILON*A*COS(2.*PHI_X)*SIGTT
     #          + SQRT(.5*EPSILON*(1.+EPSILON))*2.*A*COS(PHI_X)*SIGLT
     #          )/1.E13

C convert the cross section to lab frame from cm frame
C EXPLANATION:
C The cross section in the paper of Brauel et al (Z. Physik C3, 101(1979))
C is expressed in an invariant way and therefore is frame independent.
C If we convert the cross section to the expression in terms of solid
C angles by the factor q_cm*p_pi_cm/3.14 (Devenish and Lyth, Phys. Rev.
C D5, 47(1972)), the dOmega_e is measured in lab frame while the dOmega_pi
C is measured in the pion-neutron cm frame. (See Berends, Phys. Rev., D1,
C 2590 (1970) for example.)

      CMTOLAB = PF_P_I*PF_P_I*SQRT(S)/PPICM/
     #         ABS( (MP+OMEGA)*PF_P_I - EF_P_I*QMAG*COS(THETA_PQ) )
      SIGMA_EEP = SIGMA_EEP*CMTOLAB

      RETURN
      END
