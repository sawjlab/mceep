C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       SUBROUTINE:  ELECTRO_PROD
C
C       AUTHOR: R.W. Lourie
C       DATE:   1991
C       MODIFICATIONS:
C          1) P.E. Ulmer   8-OCT-1991
C             Convert routines to double precision for compatibility
C             with MCEEP.
C          2) G. A. Warren Aug-1996
C             Modified to match calculations of EPIPROD on a step-by-step
C             basis.  This process involved taking signs and factors out
C             in one place and putting them back in other places.  The
C             code for EPIPROD also originated from Robert Lourie, but it was
C             correct.  For instance, the born terms in this code were missing
C             a factor of i.  In addition, code did not boost polarization
C             observables from CM frame to Lab frame.  Also changed the proton
C             form factor to standard dipole rather than von Gehlen fit.
C          3) G. A. Warren July-1997
C             Removed rotation from CM to LAB frame.  Was not done
C             correctly.  Should be small effect.
C
C             N.B.  I trust this code only for neutral pion production.  I
C             cannot vouch for charge pion production!  -gaw
C
C
C       PURPOSE:
C          Computes cross section and polariation observables for pion
C          electroproduction on the proton. Uses helicity amplitude
C          formalism.  Nucleon pole (Born) terms included.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE electro_prod(pol_beam,sigma_eep,asym,
     #     PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

CXXX      DOUBLE PRECISION POL_CM(3),POL_REACT(3)

C
      CALL multipole_amplitudes               !multipoles
      CALL coefficients                       !legendre coefficients
      CALL resp_funcs                         !response functions
C
C --------------------------------------------------------------------
C     Now put pieces together for final cross sect. and polariations
C --------------------------------------------------------------------
C
      CALL xsec_and_pol(POL_BEAM,SIGMA_EEP,ASYM,
     #        PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD)

C GAW 97/9/22: There is an open question as to whether the polarizations
C              need to be rotated from the CM to LAB coordinates.  Such
C              a transformation was done in EPIPROD, but done incorrectly
C              (for versions earlier than 11/96) according to Jim Kelly.  
C              If a rotation is necessary, it is probably small.

      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
      SUBROUTINE multipole_amplitudes
C------------------------------------------------------------------------------
c
c       Computes multipole amplitudes for all major resonances up to
c       2 GeV invariant mass. Resonance parameters from Devenish+Lyth,
c       Nucl. Phys. B43 (1972), 228.
C------------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'var.cmn'
c
      DIMENSION sxqq2(30),axqq2(30),ap11(30),sp11(30)
      COMMON /p11_amplitudes/sxqq2,axqq2,ap11,sp11,amp_ap11,amp_sp11
      COMMON /p33_amplitudes/percent_c2,percent_e2
      COMMON /np11_elmts/ n_sp11,n_ap11
c
      LOGICAL   born_terms,use_p11_file
      COMMON /control_i/iopt_m1p,lmax
      COMMON /control_l/born_terms,use_p11_file
c
c       electric and magnetic multipole resonances
c
      DOUBLE COMPLEX s0p,s1p,s1m,m1m,m1p,e0p,e1p
      DOUBLE COMPLEX m2m,m2p,m3m,m3p,m4m,m4p,e2m,e2p,e3m,e3p,e4m,e4p
      COMMON /amplitudes/ s0p,s1p,s1m,m1m,m1p,e0p,e1p,
     #            m2m,m2p,m3m,m3p,m4m,m4p,e2m,e2p,e3m,e3p,e4m,e4p
      DOUBLE COMPLEX fact,factor
c
      COMMON /para/amp2,ampi2
c
      DOUBLE COMPLEX denom,e0p_0,em1_0,em2_0,em3_0,em4_0
      DOUBLE COMPLEX f_em1,f_em2,f_em3,f_em4
c
      PARAMETER (amp = 0.938279d0)
      PARAMETER (ampi = 0.1349642d0)

      amp2 = amp**2
      ampi2 = ampi**2
      root2 = sqrt(2.d0)

c
      qq2 = -qmu2_g                           !in GeV**2
      IF (QQ2.GT.0.d0) qq2 = - qq2
      gep = 1.d0/((1.d0-QQ2/0.71d0)**2)       !proton elastic form factor
      gpp = gep/(1.d0-QQ2/4.d0)               !addt'l factor for Delta
c
      w = w/1000.d0                           ! in GeV
      s = w*w                                 !CM energy^2
      ppic = f_ppic(w)                        !pion mom. in CM
      qc = f_qc(w,QQ2)
      qc0 = f_qc(w,0.d0)

c
c 0+ resonance  ----------------------------------------------------------------
c
      wr       = 1.505d0
      gammar   = 0.08d0
      xr       = 0.35d0
      sr       = wr**2
      ampr_e0p = 0.46d0
      ppicr    = f_ppic(wr)
      gamma    = gammar*(ppic/ppicr)
      denom    = dcmplx(sr-s,-wr*gamma)
      e0p_0    = wr*gamma*(ppicr/ppic)/denom
      e0p      = ampr_e0p*e0p_0*gep

      IF(iopt_m1p .LT. 0) e0p = dcmplx(0.d0,0.d0)

      s0p      = dcmplx(0.d0,0.d0)
c
c 1+ resonance (1232) ----------------------------------------------------------
c
      wr       = 1.232d0              !devinish/lyth
      gammar   = 0.114d0
      xr       = 0.167d0
      ampr_m1p = 3.52d0
      ampr_e1p = -percent_e2*ampr_m1p/100.d0
      em1_0    = f_em1(w,QQ2,ppic,qc,wr,gammar,xr)
      e1p      = ampr_e1p*gpp*em1_0
      m1p      = ampr_m1p*gpp*em1_0
      qcr      = f_qc(wr,QQ2)
      s1p      = (qc/qcr)*(-percent_c2/100.d0)*m1p

      IF (iopt_m1p.EQ.0) THEN
        m1p = dcmplx(0.d0,0.d0)
      ELSEIF (iopt_m1p.EQ.1) THEN
        m1p = -e1p
      ENDIF
      IF(iopt_m1p .EQ. -1)THEN
        s1p = dcmplx(0.d0,0.d0)
        e1p = dcmplx(0.d0,0.d0)
      ENDIF
c
c 1- resonance  ----------------------------------------------------------------
c
      wr     = 1.434d0
      gammar = 0.2d0
      xr     = 0.35d0
      IF (.NOT. use_p11_file) THEN
        ampr_m1m = amp_ap11
        em1_0    = f_em1(w,QQ2,ppic,qc,wr,gammar,xr)
        m1m      = ampr_m1m*gep*em1_0
        s1m      = qc/qcr*amp_sp11*gep*em1_0
      ELSE
        fact = factor(w,-qq2)
        qcr  = f_qc(wr,QQ2)
        m1m  = fact*rinterpq(axqq2,ap11,-qq2,n_ap11)*0.01973d0
        s1m  = fact*rinterpq(sxqq2,sp11,-qq2,n_sp11)*0.01973d0
        s1m  = -qcr/sqrt(-qq2)*s1m
      ENDIF

      IF (iopt_m1p .LT. 0) m1m = dcmplx(0.d0,0.d0)
      IF (iopt_m1p .LT. 0) s1m = dcmplx(0.d0,0.d0)
c
c 0+ resonance  ------------------------------------------------other resonances
c
      wr       = 1.63d0
      gammar   = 0.16d0
      xr       = 0.35d0
      sr       = wr**2
      ampr_e0p = -0.46d0
      ppicr    = f_ppic(wr)
      gamma    = gammar*(ppic/ppicr)
      denom    = dcmplx(sr-s,-wr*gamma)
      e0p_0    = wr*gamma*(ppicr/ppic)/denom

      IF (iopt_m1p .LT. 0) e0p_0 = dcmplx(0.d0,0.d0)
      e0p      = e0p + ampr_e0p*e0p_0*gep

      wr       = 1.7d0
      gammar   = 0.2d0
      xr       = 0.35d0
      sr       = wr**2
      ampr_e0p = 0.46d0
      ppicr    = f_ppic(wr)
      gamma    = gammar*(ppic/ppicr)
      denom    = dcmplx(sr-s,-wr*gamma)
      e0p_0    = wr*gamma*(ppicr/ppic)/denom

      IF (iopt_m1p .LT. 0) e0p_0 = dcmplx(0.d0,0.d0)
      e0p      =  e0p + ampr_e0p*e0p_0*gep
c
c 2+ resonance  ---------------------------------------------------------------
c
      wr       = 1.675d0
      gammar   = 0.134d0
      xr       = 0.35d0
      ampr_e2p = -0.03d0
      ampr_m2p =  0.07d0
      em2_0    = f_em2(w,QQ2,ppic,qc,wr,gammar,xr)
      e2p      = ampr_e2p*gep*em2_0
      m2p      = ampr_m2p*gep*em2_0
      IF (iopt_m1p .LT. 0) e2p = dcmplx(0.d0,0.d0)
      IF (iopt_m1p .LT. 0) m2p = dcmplx(0.d0,0.d0)
c
c 2- resonance  ----------------------------------------------------------------
c
      wr       = 1.519d0
      gammar   = 0.102d0
      xr       = 0.35d0
      ampr_e2m = 0.69d0
      ampr_m2m = 0.33d0
      em2_0    = f_em2(w,QQ2,ppic,qc,wr,gammar,xr)
      qcr      = f_qc(wr,QQ2)
      qcr0     = f_qc(wr,0.d0)
      r        = qcr/qcr0
      alpha    = ((ampr_e2m/ampr_m2m)**2)/3.d0
      g        = exp(-0.33d0*QQ2*QQ2)*sqrt((1.d0+alpha)/(alpha+r**4))
      e2m      = ampr_e2m*g*em2_0/((qc/qcr0)**2)
      m2m      = ampr_m2m*g*em2_0

      IF (iopt_m1p .LT. 0) e2m = dcmplx(0.d0,0.d0)
      IF (iopt_m1p .LT. 0) m2m = dcmplx(0.d0,0.d0)
c
c 3- resonance  ----------------------------------------------------------------
c
      wr       = 1.69d0
      gammar   = 0.104d0
      xr       = 0.35d0
      ampr_e3m = 0.28d0
      ampr_m3m = 0.14d0
      em3_0    = f_em3(w,QQ2,ppic,qc,wr,gammar,xr)
      qcr      = f_qc(wr,QQ2)
      qcr0     = f_qc(wr,0.d0)
      r        = qcr/qcr0
      alpha    = ((ampr_e3m/ampr_m3m)**2)/2.d0
      g        = exp(-0.33d0*QQ2*QQ2)*sqrt((1.d0+alpha)/(alpha+r**4))/r
      e3m      = ampr_e3m*g*em3_0/((qc/qcr0)**2)
      m3m      = ampr_m3m*g*em3_0

      IF (iopt_m1p .LT. 0) e3m = dcmplx(0.d0,0.d0)
      IF (iopt_m1p .LT. 0) m3m = dcmplx(0.d0,0.d0)
c
c 3+ resonance  ----------------------------------------------------------------
c
      wr       = 1.94d0
      gammar   = 0.2d0
      xr       = 0.35d0
      ampr_e3p = 0.004d0
      ampr_m3p = 0.13d0
      em3_0    = f_em3(w,QQ2,ppic,qc,wr,gammar,xr)
      e3p      = ampr_e3p*gep*em3_0
      m3p      = ampr_m3p*gep*em3_0

      IF (iopt_m1p .LT. 0) e3p = dcmplx(0.d0,0.d0)
      IF (iopt_m1p .LT. 0) m3p = dcmplx(0.d0,0.d0)
c
c 4- resonance  ----------------------------------------------------------------
c
      wr       = 2.19d0
      gammar   = 0.3d0
      xr       = 0.35d0
      ampr_e4m = 0.06d0
      ampr_m4m = 0.05d0
      em4_0    = f_em4(w,QQ2,ppic,qc,wr,gammar,xr)
      qcr0     = f_qc(wr,0.d0)
      e4m      = ampr_e4m*gep*em4_0/((qc/qcr0)**2)
      m4m      = ampr_m4m*gep*em4_0


      IF (iopt_m1p .LT. 0) e4m = dcmplx(0.d0,0.d0)
      IF (iopt_m1p .LT. 0) m4m = dcmplx(0.d0,0.d0)

      RETURN
      END

c
c-------------------------------------------------------------------------
c       Several useful little functions for manipulating complex variables
c-------------------------------------------------------------------------
c
      DOUBLE PRECISION FUNCTION preal(x,y)     !takes the Re part
      IMPLICIT NONE
      DOUBLE COMPLEX x,y

      preal = DREAL(dconjg(x)*y)

      RETURN
      END

c-------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION pimag(x,y)     !takes the Im part
      IMPLICIT NONE
      DOUBLE COMPLEX x,y

      pimag = dimag(dconjg(x)*y)

      RETURN
      END

c-------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION cnorm(x)       !takes |x|
      IMPLICIT NONE
      DOUBLE PRECISION pmag
      DOUBLE COMPLEX x

      cnorm = sqrt(pmag(x))

      RETURN
      END

c-------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION pmag(x)        !takes |x|^2
      IMPLICIT NONE
      DOUBLE COMPLEX x

      pmag= dconjg(x)*x

      RETURN
      END

c-------------------------------------------------------------------------

      DOUBLE COMPLEX FUNCTION prod(x,y)        !takes complex product
      IMPLICIT NONE
      DOUBLE COMPLEX x,y

      prod = dconjg(x)*y

      RETURN
      END
c
c-------------------------------------------------------------------------
c       Some kinematic functions
c-------------------------------------------------------------------------
c

      DOUBLE PRECISION FUNCTION f_qc(w,QQ2)      !photon momentum in COM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /para/amp2,ampi2

      qc2  = ((w*w-QQ2+amp2)/2.d0/w)**2 - amp2
      f_qc = sqrt(qc2)

      RETURN
      END

C-------------------------------------------------------------------------------
C  f_ppic
C  returns: momentum of pion in hadronic cm frame

      DOUBLE PRECISION FUNCTION f_ppic(w)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /para/amp2,ampi2

      epic   = (w*w+ampi2-amp2)/2.d0/w
      f_ppic = sqrt(epic**2-ampi2)

      RETURN
      END

c-------------------------------------------------------------------------
c       The following are the resonance functions (Breit-Wigners times
c       appropriate kinematic factors) for various mulitpolarities
c-------------------------------------------------------------------------
c
      DOUBLE COMPLEX FUNCTION f_em1(w,QQ2,ppic,qc,wr,gammar,xr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX denom

      ppicr = f_ppic(wr)
      sr    = wr*wr
      qc0   = f_qc(w,0.d0)
      qc0r  = f_qc(wr,0.d0)
      s     = w*w
      y     = (ppicr**2+xr*xr)/(ppic**2+xr*xr)
      gamma = gammar*((ppic/ppicr)**3)*y
      denom = dcmplx(sr-s,-wr*gamma)
      z = ppic/ppicr
      f_em1 = (qc/qc0)*wr*gamma/denom*z*(qc0/qc0r)/(z**3)
      RETURN
      END

C-------------------------------------------------------------------------------

      DOUBLE COMPLEX FUNCTION f_em2(w,QQ2,ppic,qc,wr,gammar,xr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX denom

      ppicr = f_ppic(wr)
      sr    = wr*wr
      qc0   = f_qc(w,0.d0)
      qc0r  = f_qc(wr,0.d0)
      s     = w*w
      y     = (ppicr**2+xr*xr)/(ppic**2+xr*xr)
      gamma = gammar*((ppic/ppicr)**5)*(y**2)
      denom = dcmplx(sr-s,-wr*gamma)
      z     = ppic/ppicr
      f_em2 = ((qc/qc0)**2)*wr*gamma/denom*(z**2)*((qc0/qc0r)**2)/(z**5)

      RETURN
      END

C-------------------------------------------------------------------------------

      DOUBLE COMPLEX FUNCTION f_em3(w,QQ2,ppic,qc,wr,gammar,xr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX denom

      ppicr = f_ppic(wr)
      sr    = wr*wr
      qc0   = f_qc(w,0.d0)
      qc0r  = f_qc(wr,0.d0)
      s     = w*w
      y     = (ppicr**2+xr*xr)/(ppic**2+xr*xr)
      gamma = gammar*((ppic/ppicr)**7)*(y**3)
      denom = dcmplx(sr-s,-wr*gamma)
      z     = ppic/ppicr
      f_em3 = ((qc/qc0)**3)*wr*gamma/denom*(z**3)*((qc0/qc0r)**3)/(z**7)

      RETURN
      END

C-------------------------------------------------------------------------------

      DOUBLE COMPLEX FUNCTION f_em4(w,QQ2,ppic,qc,wr,gammar,xr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX denom
      ppicr = f_ppic(wr)
      sr = wr*wr
      qc0 = f_qc(w,0.d0)
      qc0r = f_qc(wr,0.d0)
      s = w*w
      y = (ppicr**2+xr*xr)/(ppic**2+xr*xr)
      gamma = gammar*((ppic/ppicr)**9)*(y**4)
      denom = dcmplx(sr-s,-wr*gamma)
      z = ppic/ppicr
      f_em4 = ((qc/qc0)**4)*wr*gamma/denom*(z**4)*((qc0/qc0r)**4)/(z**9)

      RETURN
      END

c-----------------------------------------------------------------------
c       Here the 6 helicity amplitudes and then the 18 response functions
c       are computed.
c-----------------------------------------------------------------------

      SUBROUTINE resp_funcs
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'var.cmn'
      INCLUDE 'masses.cmn'

      LOGICAL   born_terms,use_p11_file
      COMMON /control_i/iopt_m1p,lmax
      COMMON /control_l/born_terms,use_p11_file
c
c       18 response functions
c

      COMMON wl_0,wl_n,wt_0,wt_n
      COMMON wtt_0,wtp_s,hborn(6)
      COMMON wtt_l,wtt_s,wtt_n,wtp_l
      COMMON wtl_0,wtlp_0,wtl_l,wtlp_l
      COMMON wtl_s,wtlp_s,wtl_n,wtlp_n
      COMMON h1,h2,h3,h4,h5,h6

      DOUBLE COMPLEX s0p,s1p,s1m,m1m,m1p,e0p,e1p
      DOUBLE COMPLEX m2m,m2p,m3m,m3p,m4m,m4p,e2m,e2p,e3m,e3p,e4m,e4p
      COMMON /amplitudes/ s0p,s1p,s1m,m1m,m1p,e0p,e1p,
     #         m2m,m2p,m3m,m3p,m4m,m4p,e2m,e2p,e3m,e3p,e4m,e4p

      DOUBLE COMPLEX H_Coeff
      DOUBLE COMPLEX H_Theta_Dep
C GAW 96/6/8: Change to 6x4 array.  Code previously added Coulomb part later
C      COMMON /factors1/H_Coeff(1:4,0:3)
C      COMMON /factors2/H_Theta_Dep(1:4,0:3)
      COMMON /factors1/H_Coeff(1:6,0:3)
      COMMON /factors2/H_Theta_Dep(1:6,0:3)

      DOUBLE COMPLEX h1,h2,h3,h4,h5,h6
      DOUBLE COMPLEX hp41,hm41,hp32,hm32
      DOUBLE COMPLEX cmplxi

      DATA pi/3.141592654/

      cmplxi = dcmplx(0.d0,1.d0)

      theta_cm_pi = pi - theta_cm
      s = sin(theta_cm_pi)
      c = cos(theta_cm_pi)
      s2 = sin(theta_cm_pi/2.d0)
      c2 = cos(theta_cm_pi/2.d0)
c
      root2 = sqrt(2.d0)
      qq2 = -qmu2_g                           !in GeV**2
      CALL legendre(theta_cm_pi)

      IF(born_terms) CALL bornt(-qq2,w,s,c,s2,c2,hborn,1)

      h1 = 0.D0                               !temporary abbreviations
      h2 = 0.D0
      h3 = 0.D0
      h4 = 0.D0
C GAW 96/6/8: add Coulomb part
      h5 = 0.D0
      h6 = 0.D0

      DO l = 0, lmax                          !sum up partial waves
        h1 = h1 + (H_Coeff(1,l) * H_Theta_Dep(1,l))
        h2 = h2 + (H_Coeff(2,l) * H_Theta_Dep(2,l))
        h3 = h3 + (H_Coeff(3,l) * H_Theta_Dep(3,l))
        h4 = h4 + (H_Coeff(4,l) * H_Theta_Dep(4,l))
C  GAW 96/6/8: add Coulomb part
        h5 = h5 + (H_Coeff(5,l) * H_Theta_Dep(5,l))
        h6 = h6 + (H_Coeff(6,l) * H_Theta_Dep(6,l))
      ENDDO

C GAW 96/6/10:  There is considerablee disagreement between EPIPROD and
c               MCEEP for these factors.  There was a know factor of i
c               that was originally missed in MCEEP.  Factors taken from
C               Lourie, Zeit Phys C50 (1991), 345, eq. A2.1
C               factor of root2 in h5,h6 differ from EPIPROD.  Signs taken from
C               EPIPROD.
C
c Original Code:
c
c       Actual helicity amplitudes
c
C      h1 = -h1*s*c2/root2 + dcmplx(hborn(1),0.d0)
C      h2 = +h2*root2*c2   + dcmplx(hborn(2),0.d0)
C      h3 = +h3*s*s2/root2 + dcmplx(hborn(3),0.d0)
C      h4 = -h4*root2*s2   + dcmplx(hborn(4),0.d0)
C      h5 = (-(s0p-s1m)+2.d0*s1p*(1.d0-3.d0*c))/root2*c2
C     #                     + dcmplx(hborn(5),0.d0)
C      h6 = ((s0p-s1m)+2.d0*s1p*(1.d0+3.d0*c))/root2*s2
C     #                     + dcmplx(hborn(6),0.d0)
c
c       switch Re<-->Im parts to account for overall factor of i
c
C       h1r = dreal(h1)
C       h1i = dimag(h1)
C       h1 = dcmplx(h1i,h1r)
C
C       h2r = dreal(h2)
C       h2i = dimag(h2)
C       h2 = dcmplx(-h2i,h2r)
C
C       h3r = dreal(h3)
C       h3i = dimag(h3)
C       h3 = dcmplx(-h3i,h3r)
C
C       h4r = dreal(h4)
C       h4i = dimag(h4)
C       h4 = dcmplx(h4i,h4r)
C
C       h5r = dreal(h5)
C       h5i = dimag(h5)
C       h5 = dcmplx(h5i,h5r)
C
C       h6r = dreal(h6)
C       h6i = dimag(h6)
C       h6 = dcmplx(-h6i,-h6r)


      h1 = +h1*s*c2/root2*cmplxi + dcmplx(0.d0,hborn(1))
      h2 = +h2*root2*c2*cmplxi   + dcmplx(0.d0,hborn(2))
      h3 = +h3*s*s2/root2*cmplxi + dcmplx(0.d0,hborn(3))
      h4 = +h4*root2*s2*cmplxi   + dcmplx(0.d0,hborn(4))
C GAW 96/6/13: Get rid of factors of root2
      h5 = -h5*c2*cmplxi   + dcmplx(0.d0,hborn(5))
      h6 = +h6*s2*cmplxi   + dcmplx(0.d0,hborn(6))

C GAW 96/6/14: Added ratio to correct longitudinal response functions
      ratio = qmag*eject_mass/w/sqrt(qmu2_g)/1.d6
      h5 = h5 * ratio
      h6 = h6 * ratio

C GAW 96/6/10: Sign of h(i) is opposite of EPIPROD.  They changed it to match
C               MW CGLN stuff.  Be Careful!!!


c
c       Combinations of helicity amplitudes that enter in response functions
c
      hp41 = h4+h1
      hm41 = h4-h1
      hp32 = h3+h2
      hm32 = h3-h2
c
c       Build the response functions -- Note that _s denotes the sideways
c       polariation component often denoted as _t for transverse.
c


C GAW 96/6/14: Multiple L components by 0.5
      wl_0   = (pmag(h5)+pmag(h6))*0.5
      wl_n   = (-2.*pimag(h6,h5))*0.5

      wt_0   = 0.5*(pmag(h1)+pmag(h2)+pmag(h3)+pmag(h4) )
      wt_n   =    pimag(h4,h2) + pimag(h3,h1)
      wtl_0  = -(preal(h5,hm41)  + preal(h6,hp32) )
      wtl_n  =   pimag(h6,hm41)  - pimag(h5,hp32)
      wtl_l   =   pimag(h5,hp41)  + pimag(h6,hm32)
      wtl_s   =   pimag(h6,hp41)  - pimag(h5,hm32)
      wtt_0  =   preal(h4,h1) - preal(h3,h2)
C GAW 96/6/13: Removed factors of root2
C Original Code:
C       wtt_n  =   (pimag(h1,h2) - pimag(h4,h3))/root2
C       wtt_l   = -(pimag(h3,h2) - pimag(h4,h1) )/root2
C       wtt_s   = -(pimag(h4,h3) + pimag(h1,h2) )/root2
C       wtlp_0 =   pimag(h5,hm41)  + pimag(h6,hp32)
C       wtlp_n =   root2*(preal(h6,hm41)  - preal(h5,hp32))
C       wtlp_l  = -(preal(h5,hp41)  + preal(h6,hm32) )*root2
C       wtlp_s  = -(preal(h6,hp41)  - preal(h5,hm32) )*root2
C       wtp_l   = 0.5*(pmag(h3)+pmag(h4)-pmag(h1)-pmag(h2) )
C       wtp_s   = preal(h2,h4) + preal(h1,h3)
      wtt_n  =   (pimag(h1,h2) - pimag(h4,h3))
      wtt_l   = -(pimag(h3,h2) - pimag(h4,h1) )
      wtt_s   = -(pimag(h4,h3) + pimag(h1,h2) )
      wtlp_0 =   pimag(h5,hm41)  + pimag(h6,hp32)
      wtlp_n =   (preal(h6,hm41)  - preal(h5,hp32))
      wtlp_l  = -(preal(h5,hp41)  + preal(h6,hm32) )
      wtlp_s  = -(preal(h6,hp41)  - preal(h5,hm32) )
      wtp_l   = 0.5*(pmag(h3)+pmag(h4)-pmag(h1)-pmag(h2) )
      wtp_s   = preal(h2,h4) + preal(h1,h3)

      RETURN
      END

c---------------------------------------------------------------------------
c       LEGENDRE calculates, for some angle ALPHA, the terms involving
c       the Legendre polynomials and their derivitives that appear in
c       the partial wave expansion of the helicity amplitudes
c---------------------------------------------------------------------------
c
c 96/6/8   GAW    Added Coulomb terms

      SUBROUTINE legendre(alpha)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      DOUBLE COMPLEX H_Theta_Dep
C GAW 96/6/8: Changed to 6x4 array
C      COMMON /factors2/H_Theta_Dep(1:4,0:3)
      COMMON /factors2/H_Theta_Dep(1:6,0:3)
c
      s = sin(alpha)
      c = cos(alpha)

c
c  Here the second index 0,1,2,3,4,5 labels the pi-N angular momentum
c
      H_Theta_Dep(1,0) = 0.d0
      H_Theta_Dep(1,1) = -3.d0
      H_Theta_Dep(1,2) = +3.d0 - 15.d0*c
      H_Theta_Dep(1,3) = 15.d0*c - 52.5d0*(c**2) + 7.5d0

      H_Theta_Dep(2,0) = -1.d0
      H_Theta_Dep(2,1) = +1.d0 - 3.d0*c
      H_Theta_Dep(2,2) = +3.d0*c - 7.5d0*(c**2) + 1.5d0
      H_Theta_Dep(2,3) = +7.5d0*(c**2) - 1.5d0 - 17.5d0*(c**3) + 7.5d0*c

      H_Theta_Dep(3,0) = 0.d0
      H_Theta_Dep(3,1) = +3.d0
      H_Theta_Dep(3,2) = +3.d0 + 15.d0*c
      H_Theta_Dep(3,3) = 15.d0*c + 52.5d0*(c**2) - 7.5d0

      H_Theta_Dep(4,0) = +1.d0
      H_Theta_Dep(4,1) = +1.d0 + 3.d0*c
      H_Theta_Dep(4,2) = +3.d0*c + 7.5d0*(c**2) -1.5d0
      H_Theta_Dep(4,3) = +7.5d0*(c**2) - 1.5d0 + 17.5d0*(c**3) - 7.5d0*c

C GAW 96/6/8: Added in Coulomb part

      DO i=0,3
        H_Theta_Dep(5,i) = H_Theta_Dep(2,i)
        H_Theta_Dep(6,i) = H_Theta_Dep(4,i)
      END DO

      RETURN
      END

c---------------------------------------------------------------------------
c       COEFFICIENTS calculates the linear combinations of the multipoles
c       that appear in the partial wave expansion of the helicity amplitudes.
c       These are the coefficients of the Legendre terms.
c---------------------------------------------------------------------------
c
C     96/6/8   GAW     Added Coulomb part to H_coeff after EPIPROD pattern

      SUBROUTINE coefficients
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'var.cmn'       !  GAW 96/6/8: need to bring in Qmu2_g & W
c
      DOUBLE COMPLEX H_Coeff
C GAW 96/6/8: Convert to 6x4 array
C      COMMON /factors1/H_Coeff(1:4,0:3)
      COMMON /factors1/H_Coeff(1:6,0:3)
      DOUBLE COMPLEX s0p,s1p,s1m,m1m,m1p,e0p,e1p
      DOUBLE COMPLEX m2m,m2p,m3m,m3p,m4m,m4p,e2m,e2p,e3m,e3p,e4m,e4p
      COMMON /amplitudes/ s0p,s1p,s1m,m1m,m1p,e0p,e1p,
     #            m2m,m2p,m3m,m3p,m4m,m4p,e2m,e2p,e3m,e3p,e4m,e4p
      DOUBLE COMPLEX c0p,c1m,c1p,c2m,c2p,c3m,c3p,c4m       ! GAW 96/6/8

C GAW 96/6/8: Convert Sl amplitudes to Cl amplitudes as in epiprod

      qvec = f_qc(W,-Qmu2_g)
      c0p = sqrt(Qmu2_g)/qvec*s0p
      c1p = sqrt(Qmu2_g)/qvec*2*s1p
      c2p = 0.d0
      c3p = 0.d0
      c1m = -sqrt(Qmu2_g)/qvec*s1m
      c2m = 0.d0
      c3m = 0.d0
      c4m = 0.d0

C  H_Coeff(1,l) = B_l+ - B_(l+1)-
      H_Coeff(1,0) = 0.d0
      H_Coeff(1,1) =           e1p -    m1p -    e2m -    m2m
      H_Coeff(1,2) =           e2p -    m2p -    e3m -    m3m
      H_Coeff(1,3) =           e3p -    m3p -    e4m -    m4m

C  H_Coeff(2,l) = A_l+ - A_(l+1)-
      H_Coeff(2,0) =               e0p                       -      m1m
      H_Coeff(2,1) = 0.5d0 * (3.d0*e1p +      m1p +      e2m - 3.d0*m2m)
      H_Coeff(2,2) =         (2.d0*e2p +      m2p +      e3m - 2.d0*m3m)
      H_Coeff(2,3) = 0.5d0 * (5.d0*e3p + 3.d0*m3p + 3.d0*e4m - 5.d0*m4m)

C  H_Coeff(3,l) = B_l+ + B_(l+1)-
      H_Coeff(3,0) = 0.d0
      H_Coeff(3,1) =       (   e1p -    m1p +    e2m +    m2m)
      H_Coeff(3,2) =       (   e2p -    m2p +    e3m +    m3m)
      H_Coeff(3,3) =       (   e3p -    m3p +    e4m +    m4m)

C  H_Coeff(4,l) = A_l+ + A_(l+1)-
      H_Coeff(4,0) =               e0p                       +      m1m
      H_Coeff(4,1) = 0.5d0 * (3.d0*e1p +      m1p -      e2m + 3.d0*m2m)
      H_Coeff(4,2) =         (2.d0*e2p +      m2p -      e3m + 2.d0*m3m)
      H_Coeff(4,3) = 0.5d0 * (5.d0*e3p + 3.d0*m3p - 3.d0*e4m + 5.d0*m4m)

C  H_Coeff(5,l) =  C_l+ - C_(l+1)-
      H_Coeff(5,0) = (  c0p -   c1m     )
      H_Coeff(5,1) = (  c1p -   c2m     )
      H_Coeff(5,2) = (  c2p -   c3m     )
      H_Coeff(5,3) = (  c3p -   c4m     )

C  H_Coeff(6,l) =  C_l+ + C_(l+1)-
      H_Coeff(6,0) = (  c0p +   c1m     )
      H_Coeff(6,1) = (  c1p +   c2m     )
      H_Coeff(6,2) = (  c2p +   c3m     )
      H_Coeff(6,3) = (  c3p +   c4m     )


      RETURN
      END
c
c---------------------------------------------------------------------------
      SUBROUTINE xsec_and_pol(POL_BEAM,SIGMA_EEP,ASYM,
     #     PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD)
c---------------------------------------------------------------------------
c
c       takes the response functions and computes the proton
c       polarization vector, electron analyzing power (ASYM) and
c       unpolarized cross section (fm^2/MeV/sr^2) in the LAB
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'var.cmn'
      INCLUDE 'masses.cmn'
c
      DOUBLE PRECISION k_gamma,jacobian
c
c       18 response functions
c
      COMMON wl_0,wl_n,wt_0,wt_n
      COMMON wtt_0,wtp_s,hborn(6)
      COMMON wtt_l,wtt_s,wtt_n,wtp_l
      COMMON wtl_0,wtlp_0,wtl_l,wtlp_l
      COMMON wtl_s,wtlp_s,wtl_n,wtlp_n
c
      DATA alpha/7.2974d-03/,pi/3.14159265d0/


      cth_pq = cos(theta_pq)
      q2 = -qmu2_g*1.d6        !in (MeV/c)^2
      w = w*1000.d0
      q_cm = eject_mass*qmag/w
      p_cm = sqrt(((w**2-eject_mass**2+miss_m**2)/2.d0/w)**2-miss_m**2)
c
c       Gamma is the so-called virtual photon flux
c

      gamma = -alpha/2.d0/pi**2*(e0_i-omega)/e0_i*(w**2-eject_mass**2)
     #          /2.d0/eject_mass/q2/(1.d0-epsilon)
c
      k_gamma = (w**2-eject_mass**2)/2.d0/eject_mass
      gamma = gamma*p_cm*w/k_gamma/eject_mass
      tt = tan(tscat/2.d0)
c
c       electron kinematic factors
c
      term = q2/qmag/qmag
      rho_l = -epsilon*term*(w/eject_mass)**2
      vl = 2.d0*rho_l
      vt = 1.0d0
      vtt = -epsilon
      vlt = -sqrt(rho_l*(1.d0+epsilon))
      vt_p = sqrt(1.d0-epsilon*epsilon)
      vlt_p = -sqrt(rho_l*(1.d0-epsilon))
c
      e_p = sqrt(pf_p_i**2 + eject_mass**2)
      cx = cos(phi_x)
      cx2 = cos(2.d0*phi_x)
      sx = sin(phi_x)
      sx2 = sin(2.d0*phi_x)
c
c      Build up cross section, electron analyzing power (ASYM) and
c      polarization functions (=cross section times polarization
c      component) from the response functions
c
      unpol_piece = vl*wl_0 + vt*wt_0 + vlt*wtl_0*cx + vtt*wtt_0*cx2

      asym_piece  = POL_BEAM*vlt_p*wtlp_0*sx

      pl_hi_piece = vtt*wtt_l*sx2 + vlt*wtl_l*sx

      pl_hd_piece = POL_BEAM*(vt_p*wtp_l + vlt_p*wtlp_l*cx)

      pt_hi_piece = vtt*wtt_s*sx2 + vlt*wtl_s*sx

      pt_hd_piece = POL_BEAM*(vt_p*wtp_s + vlt_p*wtlp_s*cx)

      pn_hi_piece = vl*wl_n + vt*wt_n + vlt*wtl_n*cx + vtt*wtt_n*cx2

      pn_hd_piece = POL_BEAM*vlt_p*wtlp_n*sx

      term1 = qmag*pf_p_i/q_cm/p_cm
      term2=abs(1.d0+omega/eject_mass-cth_pq*e_p*qmag/eject_mass/pf_p_i)

      jacobian = term1/term2                          !CM-->LAB

      SIGMA_EEP = gamma*jacobian*unpol_piece*1.0d-04          !in fm^2
      ASYM      = gamma*jacobian*asym_piece *1.0d-04
c
      PN_HI     = gamma*jacobian*pn_hi_piece   *1.0d-04
      PT_HI     = gamma*jacobian*pt_hi_piece   *1.0d-04
      PL_HI     = gamma*jacobian*pl_hi_piece   *1.0d-04
c
      PN_HD     = gamma*jacobian*pn_hd_piece   *1.0d-04
      PT_HD     = gamma*jacobian*pt_hd_piece   *1.0d-04
      PL_HD     = gamma*jacobian*pl_hd_piece   *1.0d-04
c
      RETURN
      END

c---------------------------------------------------------------------------
c       helicity amplitudes for Born terms
c---------------------------------------------------------------------------

      SUBROUTINE  bornt(q2,w,sinx,cosx,sinx2,cosx2,born,irea)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION born(6),hborn(6)
      DIMENSION
     *   xm(6,6),a(6),hi(6),iepsi(-1:+1),xn(-1:1,6),
     *   ai(-1:1,6),ieta(6),f1(-1:1),f2(-1:1),fp(-1:1)

      DATA iepsi/-1,1,1/,ieta/1,1,-1,1,-1,-1/
      DATA pi/3.1415927d0/,xmp/.9382796d0/,
     +  xmpip/0.139563d0/,xmpi0/0.1349630d0/,
     +  mag/1/
C GAW 96/6/12: Use pi0+- mass to calculate gr for proper reaction
C     DATA er/.302862d0/,gr/13.5656d0/,sq2/1.41421d0/
      DATA er/.302862d0/,sq2/1.41421d0/

      IF (irea.EQ.1.OR.irea.EQ.2) THEN
        xmpi = xmpi0
        gr   = 14.0280d0                 ! GAW 96/6/12
      ELSEIF (irea.EQ.3 .OR. irea.EQ.4) THEN
        xmpi = xmpip
        gr = 13.5656d0                   ! GAW 96/6/12
      ENDIF

      s      = w*w
      sigma  = s - xmp**2
      q0     = (sigma-q2)/(2.d0*w)
      q      = sqrt(q0**2 + q2)
      xk0pi  = (sigma+xmpi**2)/(2.d0*w)
      xkpi   = sqrt(xk0pi**2-xmpi**2)
      beta   = sigma + q2/2.d0
      xnue   = -(w + xmp)
      omega  =   w - xmp
      e1     =   w + xmp - q0
      e2     =   w + xmp - xk0pi
      g1     =  -w + xmp + q0
      g2     =  -w + xmp + xk0pi
      t      = 2.d0*xkpi*q*cosx-2.d0*xk0pi*q0+xmpi**2-q2
      u      = s - 2.d0*beta - t + xmpi**2
      x1s    = xmp*xmp - s
      x1u    = xmp*xmp - u
      p2     = t - xmpi*xmpi
      fp(-1) = er/(1.d0+q2/0.5d0)
      fp(0)  = 0.d0
      fp(1)  = 0.d0

C GAW 96/6/12: Replace original vG fit with standard dipole form
C      gep    = 1.d0/(1.d0+3.04d0*q2+1.54d0*q2*q2+.068d0*q2*q2*q2)
      gep    = 1.d0 / (1.d0 + q2/0.71)**2

      cap    = q2/(4.d0*xmp*xmp)
      f1p    = gep*(1.d0+2.793d0*cap)/(1.d0+cap)
      f1n    = gep*(  -1.913d0*cap)/(1.d0+cap)

      IF(mag.EQ.0)THEN
        f2p = 0.d0
        f2n = 0.d0
      ELSE

C GAW 96/6/12:  F2n Assumes that Gen are zero

        f2p = gep*1.793d0/(1.d0+cap)
        f2n = gep*(-1.913d0)/(1.d0+cap)
        f2p = f2p/(2.d0*xmp)
        f2n = f2n/(2.d0*xmp)
      END IF

      f1(+1)  = er*(f1p-f1n)
      f1(-1)  = f1(+1)
      f1(0 )  = er*(f1p+f1n)
      f2(+1)  = er*(f2p-f2n)
      f2(-1)  = f2(+1)
      f2(0 )  = er*(f2p+f2n)


      p1      = xk0pi-q0
      xm(1,1) = -4.d0*omega*e1
      xm(1,2) = beta*(q2+2.d0*xnue*p1-p2)
      xm(1,3) = -2.d0*e1*omega*omega+(e1-q0)*p2+q2*(-xnue+2.d0*p1)
      xm(1,4) = -2.d0*e1*omega*omega-(e1-q0)*p2-q2*(-xnue+2.d0*p1)
      xm(1,5) = -q2*(q2+2.d0*xnue*p1-p2)
      xm(1,6) = -4.d0*q2*e1

      p1      = -(xk0pi-q0)
      xm(2,1) = -4.d0*xnue*g1
      xm(2,2) = beta*(q2+2.d0*omega*p1-p2)
      xm(2,3) = -2.d0*g1*xnue*xnue+(g1+q0)*p2+q2*(-omega+2.d0*p1)
      xm(2,4) = -2.d0*g1*xnue*xnue-(g1+q0)*p2-q2*(-omega+2.d0*p1)
      xm(2,5) = -q2*(q2+2.d0*omega*p1-p2)
      xm(2,6) = -4.d0*q2*g1

      xm(3,1) = 0.d0
      xm(3,2) = -beta*2.d0*q*q
      xm(3,3) =  xnue*2.d0*q*q
      xm(3,4) = -xnue*2.d0*q*q
      xm(3,5) = q2*2.d0*q*q
      xm(3,6) = 0.d0

      xm(4,1) = 0.d0
      xm(4,2) = -beta*2.d0*q*q
      xm(4,3) =  omega*2.d0*q*q
      xm(4,4) = -omega*2.d0*q*q
      xm(4,5) = q2*2.d0*q*q
      xm(4,6) = 0.d0

      p1      =  2.d0*xk0pi-q0
      xm(5,1) =  2.d0*e1
      xm(5,2) =  beta*p1+(2.d0*w-1.5d0*q0)*p2
      xm(5,3) =  e1*omega-xnue*p1+p2
      xm(5,4) =  e1*omega+xnue*p1-p2
      xm(5,5) = -q2*p1+q0*p2
      xm(5,6) = -2.d0*omega*e1

      p1      = -2.d0*xk0pi+q0
      xm(6,1) =  2.d0*g1
      xm(6,2) =  beta*p1-(2.d0*w-1.5d0*q0)*p2
      xm(6,3) =  g1*xnue-omega*p1+p2
      xm(6,4) =  g1*xnue+omega*p1-p2
      xm(6,5) = -q2*p1-q0*p2
      xm(6,6) = -2.d0*xnue*g1

      DO i = 1,6
        DO j = 1,6
          xm(i,j) = xm(i,j)/(2.d0*q*q)
        ENDDO
      ENDDO

      hi(1) =  g1*sqrt(e1*e2)/(8.d0*pi*w)
      hi(2) = -e1*sqrt(g1*g2)/(8.d0*pi*w)
      hi(3) =  e2*sqrt(g1*g2)/(8.d0*pi*w)
      hi(4) =  g2*sqrt(e1*e2)/(8.d0*pi*w)
      hi(5) = -g1*sqrt(e1*e2)/(8*pi*w)
      hi(6) =  e1*sqrt(g1*g2)/(8.d0*pi*w)

      i1=0
      DO 71 i2=1,6
        xn(i1,i2) = gr*(1.d0/x1s + iepsi(i1)*ieta(i2)/x1u )
   71 CONTINUE

      ai(i1,1) = -f1(i1)*xn(i1,1)/2.d0
      ai(i1,2) =  f1(i1)*xn(i1,2)/p2
      ai(i1,3) =  f2(i1)*xn(i1,3)/2.d0
      ai(i1,4) =  f2(i1)*xn(i1,4)/2.d0

      IF (q2.NE.0.d0) ai(i1,5) = f1(i1)*xn(i1,5)/(2.d0*p2)
     #        + (1-iepsi(i1))*(fp(i1)-f1(i1))*gr/(q2*p2)

      IF (q2.EQ.0.d0) ai(i1,5) = 0.d0

      ai(i1,6) = 0.d0

      IF (irea.EQ.1.OR.irea.EQ.2) i1 =  1
      IF (irea.EQ.3.OR.irea.EQ.4) i1 = -1

      DO 72 i2=1,6
        xn(i1,i2) = gr*(1.d0/x1s + iepsi(i1)*ieta(i2)/x1u )
   72 CONTINUE

      ai(i1,1) = -f1(i1)*xn(i1,1)/2.d0
      ai(i1,2) =  f1(i1)*xn(i1,2)/p2
      ai(i1,3) =  f2(i1)*xn(i1,3)/2.d0
      ai(i1,4) =  f2(i1)*xn(i1,4)/2.d0

      IF (q2.NE.0.d0) ai(i1,5) = f1(i1)*xn(i1,5)/(2.d0*p2)
     #  + (1-iepsi(i1))*(fp(i1)-f1(i1))*gr/(q2*p2)

      IF (q2.EQ.0.d0) ai(i1,5)=0.d0

      ai(i1,6) = 0.d0


      IF (irea.EQ.1) THEN
        DO 73 i2=1,6
   73   a(i2) =  ai(0,i2) + ai(1,i2)
      END IF

      IF (irea.EQ.2) THEN
        DO 74 i2=1,6
   74   a(i2) = -ai(0,i2) + ai(1,i2)
      END IF

      IF (irea.EQ.3) THEN
        DO 75 i2=1,6
   75   a(i2) = -sq2*(ai(0,i2) + ai(-1,i2))
      END IF

      IF (irea.EQ.4) THEN
        DO 76 i2=1,6
   76   a(i2) =  sq2*(ai(0,i2) - ai(-1,i2))
      END IF


      DO 2 i=1,6
        hborn(i) = 0.d0
        DO 1 j=1,6
    1   hborn(i) = hborn(i)+a(j)*xm(i,j)
    2 hborn(i) = hborn(i)*hi(i)*19.732857d0


C GAW 96/6/11: From EPIPROD. Removed factors of +-1
C
C       born(1) =  cmplxi/root2*sthcm_pi*chth_pi* ( hborn(3) + hborn(4) )
C       born(2) = -cmplxi/root2*chth_pi         * ( hborn(1) + hborn(2) )
C       born(3) = -cmplxi/root2*sthcm_pi*shth_pi* ( hborn(3) - hborn(4) )
C       born(4) =  cmplxi/root2*shth_pi         * ( hborn(1) - hborn(2) )
C       born(5) =  cmplxi*sqrt(Q2)*chth_pi      * ( hborn(5) + hborn(6) )
C       born(6) = -cmplxi*sqrt(Q2)*shth_pi      * ( hborn(5) - hborn(6) )

      born(1) =  sinx*cosx2*(hborn(3)+hborn(4))/sq2
      born(2) = -cosx2*(hborn(1)+hborn(2))/sq2
      born(3) = -sinx*sinx2*(hborn(3)-hborn(4))/sq2
      born(4) =  sinx2   *  (hborn(1)-hborn(2))/sq2
C GAW 96/6/13: Replace q0 by sqrt(q2) to match EPIPROD & change sign of born(5)
      born(5) =  sqrt(q2)*cosx2*(hborn(5)+hborn(6))
      born(6) =  sqrt(q2)*sinx2*(hborn(5)-hborn(6))


      RETURN
      END
C
c-----------------------------------------------------------------------------
c       Resonance function used if external Roper amplitudes employed
c-----------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION factor(w,q2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX denom
      DATA amp/0.93828d0/,pi/3.14159265d0/,ampi/0.13496d0/
      DATA xpin/0.35d0/,rj/0.5d0/
      DATA wr/1.440d0/,gamma/0.2d0/
      alpha = 1.d0/137.d0

      sr    = wr*wr
      s     = w*w
      denom = wr/dcmplx(sr-s,-wr*gamma)

      q_pi_c2 = (wr**2-(ampi+amp)**2)*(wr**2-(ampi-amp)**2)/4.d0/wr**2
      q_kc2   = q2 + (wr**2-amp**2-q2)**2/4/wr/wr
      q_kw    = (wr**2-amp**2)/2.d0/wr

      factor1 = sqrt(  sqrt(q_kc2)*xpin*amp*gamma       /
     1               ( pi*wr*sqrt(q_pi_c2)*(2.d0*rj+1) ) )
      factor  = -denom*factor1

      RETURN
      END

