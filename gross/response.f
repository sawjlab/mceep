c_______________________________________________________________________
c
c*****respns
c
c     This subroutine constructs the response functions from the mixed
c     helicity-transversity amplitudes g(i).
c
c     INPUT
c          Value.........................x (if cm frame) or p (if lab frame)
c          Q2............................invariant four momentum contraction
c          theta.........................polar angle in cm or lab frame
c          iff...........................code to choose form factor type
c          ilab..........................code to choose lab frame evaluation
c
c     OUTPUT
c          resp(ideu,inuc,iresp).........array containing respons functions
c
c          ideu=0........................unpolarized deuteron
c          ideu=1........................sqrt(1.5)T_{10}
c          ideu=2........................T_{20)/sqrt(2)
c          ideu=3........................sqrt(3)Re(T_{22})
c          ideu=4........................sqrt(3)Im(T_{22})
c          ideu=5........................sqrt(1.5)Re(T_{11})
c          ideu=6........................sqrt(1.5)Im(T_{11})
c          ideu=7........................sqrt(1.5)Re(T_{21})
c          ideu=8........................sqrt(1.5)Im(T_{21})
c
c          inuc=0........................unpolarized nucleon
c          inuc=1........................P_n
c          inuc=2........................P_s
c          inuc=3........................P_l
c
c          iresp=1.......................R_L
c          iresp=2.......................R_T
c          iresp=3.......................R_{TT}^{(I)}
c          iresp=4.......................R_{LT}^{(I)}
c          iresp=5.......................R_{LT'}^{(I)}
c          iresp=6.......................R_{TT}^{(II)}
c          iresp=7.......................R_{LT}^{(II)}
c          iresp=8.......................R_{LT'}^{(II)}
c          iresp=9.......................R_{T'}
c
c
      SUBROUTINE respns(Value,Q2,theta,W,iff,ilab,resp)
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
      DIMENSION resp(0:8,0:3,9)
c
      COMPLEX g(18),curme(0:1,0:1,0:1,-1:1)
c
      DOUBLE PRECISION mn,mt,kappa2
      COMMON/dpara/mn,mt
c
      PARAMETER(pi=3.141592653589793)
c
      IF (ILab .EQ. 0) THEN
        x = Value
        OmegaL = Q2 / (2. * mn * x)
        qL = Sqrt (Q2 + OmegaL**2)
        W = Sqrt (MT**2 + 2. * MT * OmegaL - Q2)
        ThetaC = Theta
      ELSE
        OmegaL = (W**2 + Q2 - MT**2) / 2. / MT
        qL = Sqrt (Q2 + OmegaL**2)
        p1 = Value
        pcx = p1 * Sin (Theta)
        pcz = p1 * Cos (Theta) * Sqrt (W**2 + qL**2) / W -
     1           qL / W * Sqrt (p1**2 + mn**2)
        pc = Sqrt (pcx**2 + pcz**2)
cxxx        W = 2. * Sqrt (mn**2 + pc**2)
        ThetaC = ACos (pcz / pc)
        ThetaL = Theta
      END IF
      Q=sqrt(Q2)
c                      To calculate the covariant form of the
c                      cross section eta should be one
      eta=1.
c      eta=mt*ql/(w*q)    !version for mixed frame calculation
      phic=0.0
c
c              Note: The value of kappa**2 below is the kappa*2 of D and G
c                    multiplied by 2*M_T
c
      kappa2=mt*mn**2/(2.*pi**2*w)
c
c             initialize response function array
c
      DO 300 iresp=1,9
        DO 200 inuc=0,3
          DO 100 ideu=0,8
            resp(ideu,inuc,iresp)=0.
  100     CONTINUE
  200   CONTINUE
  300 CONTINUE
c
c          call routines to give matrix elements and convert to
c                            mixed amplitudes
c
      CALL pwia(ql,omegal,w,thetac,phic,iff,curme)
      CALL trnsvr(curme,g)
c      Do 965 iii = 1, 18
c  965    Write(6,*)'g(',iii,')',g(iii)
c
CXXX      WRITE(10,1) ql,omegal,w,thetac,phic
    1 FORMAT(/,6(2x,1pe11.4))
c
c            call subroutines to fill response function array
c
      CALL TABL10(g,eta,kappa2,resp)
      CALL TABL11(g,eta,kappa2,resp)
      CALL TABL12(g,eta,kappa2,resp)

c          Perform Wigner Rotation if response functions are to be
c          calculated for the lab frame

      IF (ILab .EQ. 1) THEN
        Gamma = Sqrt (W **2 + qL ** 2) / W
        GammaL = Sqrt (mn ** 2 + p1 ** 2) / mn
        Beta = qL / W / Gamma
        BetaL = p1 / mn / GammaL
        SinTht = Sin (ThetaL)
        CosTht = Cos (ThetaL)
        Denom = 1 + Gamma * GammaL * (1 - Beta * BetaL * CosTht)
        SinPhi = (Gamma * Beta * GammaL * BetaL -
     1             (GammaL - 1) * (Gamma - 1) * CosTht) * SinTht
     2             / Denom
        CosPhi = ((Gamma + GammaL) * SinTht ** 2 +
     1            (1 + Gamma * GammaL) * CosTht ** 2 -
     2            Gamma * GammaL * Beta * BetaL * CosTht) / Denom
        DO 400 IDeu = 0, 8
          DO 400 IResp = 1, 9
            RespPs = Resp (IDeu, 2, IResp)
            RespPl = Resp (IDeu, 3, IResp)
            Resp (IDeu, 2, IResp) = RespPs * CosPhi - RespPl * SinPhi
  400   Resp (IDeu, 3, IResp) = RespPs * SinPhi + RespPl * CosPhi
      END IF

      RETURN
      END
c
c_______________________________________________________________________
c
c*****TABL10
c
c     This subroutine reproduces the response functions as given in
c     table 10 of Gross and Dmitrasinovic.
c
      SUBROUTINE TABL10(g,eta,kappa2,resp)
c
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      DOUBLE PRECISION kappa2
c
      DIMENSION resp(0:8,0:3,9)
c
      COMPLEX g(18),a1,a2,a3,b1,b2,b3,c1,c2,d1,d2,e1,e2,e3,e4,
     &        f1,f2,f3,f4
c
      fac=4.*eta*kappa2/3.
c
      a1=conjg(g(2))*g(14)+conjg(g(4))*g(16)+conjg(g(6))*g(18)
      a2=conjg(g(2))*g(14)-conjg(g(6))*g(18)
      a3=conjg(g(2))*g(14)-2.*conjg(g(4))*g(16)+conjg(g(6))*g(18)
c
      b1=conjg(g(1))*g(13)+conjg(g(3))*g(15)+conjg(g(5))*g(17)
      b2=conjg(g(1))*g(13)-conjg(g(5))*g(17)
      b3=conjg(g(1))*g(13)-2.*conjg(g(3))*g(15)+conjg(g(5))*g(17)
c
      c1=conjg(g(2))*g(18)+conjg(g(6))*g(14)
      c2=conjg(g(2))*g(18)-conjg(g(6))*g(14)
c
      d1=conjg(g(1))*g(17)+conjg(g(5))*g(13)
      d2=conjg(g(1))*g(17)-conjg(g(5))*g(13)
c
      e1=conjg(g(1)+g(5))*g(16)+conjg(g(3))*(g(18)+g(14))
      e2=conjg(g(1)-g(5))*g(16)+conjg(g(3))*(g(18)-g(14))
      e3=conjg(g(1)-g(5))*g(16)-conjg(g(3))*(g(18)-g(14))
      e4=conjg(g(1)+g(5))*g(16)-conjg(g(3))*(g(18)+g(14))
c
      f1=conjg(g(2)+g(6))*g(15)+conjg(g(4))*(g(17)+g(13))
      f2=conjg(g(2)-g(6))*g(15)+conjg(g(4))*(g(17)-g(13))
      f3=conjg(g(2)-g(6))*g(15)-conjg(g(4))*(g(17)-g(13))
      f4=conjg(g(2)+g(6))*g(15)-conjg(g(4))*(g(17)+g(13))
c
      resp(0,0,4)=fac*real(a1+b1)
      resp(0,1,4)=fac*real(a1-b1)
      resp(1,0,4)=fac*real(a2+b2)
      resp(1,1,4)=fac*real(a2-b2)
      resp(2,0,4)=fac*real(a3+b3)
      resp(2,1,4)=fac*real(a3-b3)
      resp(3,0,4)=fac*real(c1+d1)
      resp(3,1,4)=fac*real(c1-d1)
      resp(4,0,4)=-fac*aimag(c2+d2)
      resp(4,1,4)=-fac*aimag(c2-d2)
c
      resp(5,2,4)=-fac*real(e1+f1)
      resp(5,3,4)=-fac*aimag(e1-f1)
      resp(6,2,4)=fac*aimag(e2+f2)
      resp(6,3,4)=-fac*real(e2-f2)
      resp(7,2,4)=-fac*real(e3+f3)
      resp(7,3,4)=-fac*aimag(e3-f3)
      resp(8,2,4)=fac*aimag(e4+f4)
      resp(8,3,4)=-fac*real(e4-f4)
c
      resp(0,0,5)=fac*aimag(a1+b1)
      resp(0,1,5)=fac*aimag(a1-b1)
      resp(1,0,5)=fac*aimag(a2+b2)
      resp(1,1,5)=fac*aimag(a2-b2)
      resp(2,0,5)=fac*aimag(a3+b3)
      resp(2,1,5)=fac*aimag(a3-b3)
      resp(3,0,5)=fac*aimag(c1+d1)
      resp(3,1,5)=fac*aimag(c1-d1)
      resp(4,0,5)=fac*real(c2+d2)
      resp(4,1,5)=fac*real(c2-d2)
c
      resp(5,2,5)=-fac*aimag(e1+f1)
      resp(5,3,5)=fac*real(e1-f1)
      resp(6,2,5)=-fac*real(e2+f2)
      resp(6,3,5)=-fac*aimag(e2-f2)
      resp(7,2,5)=-fac*aimag(e3+f3)
      resp(7,3,5)=fac*real(e3-f3)
      resp(8,2,5)=-fac*real(e4+f4)
      resp(8,3,5)=-fac*aimag(e4-f4)
c
      RETURN
      END
c
c_______________________________________________________________________
c
c*****TABL11
c
c     This subroutine reproduces the response functions as given in
c     table 11 of Gross and Dmitrasinovic.
c
      SUBROUTINE TABL11(g,eta,kappa2,resp)
c
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      DOUBLE PRECISION kappa2
c
      DIMENSION resp(0:8,0:3,9)
c
      COMPLEX g(18),c,d,alpha1,alpha2,beta1,beta2
c
      DO 100 iresp=1,3
c
        IF(iresp.GT.1) THEN
          fac=4.*kappa2/3.
        ELSE
          fac=kappa2*eta**2
        ENDIF
c
        ioff=18-6*iresp
        epsiln=-(-1)**iresp
c
        a1=abs(g(2+ioff))**2+abs(g(4+ioff))**2+abs(g(6+ioff))**2
        a2=abs(g(2+ioff))**2-abs(g(6+ioff))**2
        a3=abs(g(2+ioff))**2-2.*abs(g(4+ioff))**2+abs(g(6+ioff))**2
c
        b1=abs(g(1+ioff))**2+abs(g(3+ioff))**2+abs(g(5+ioff))**2
        b2=abs(g(1+ioff))**2-abs(g(5+ioff))**2
        b3=abs(g(1+ioff))**2-2.*abs(g(3+ioff))**2+abs(g(5+ioff))**2
c
        c=2.*conjg(g(2+ioff))*g(6+ioff)
c
        d=2.*conjg(g(1+ioff))*g(5+ioff)
c
        alpha1=2.*(conjg(g(1+ioff))*g(4+ioff)+conjg(g(3+ioff))*
     &         g(6+ioff))
        alpha2=2.*(conjg(g(1+ioff))*g(4+ioff)-conjg(g(3+ioff))*
     &         g(6+ioff))
c
        beta1=2.*(conjg(g(3+ioff))*g(2+ioff)+conjg(g(5+ioff))*g(4+ioff))
        beta2=2.*(conjg(g(3+ioff))*g(2+ioff)-conjg(g(5+ioff))*g(4+ioff))
c
        resp(0,0,iresp)=fac*(a1+b1)
        resp(0,1,iresp)=fac*epsiln*(a1-b1)
        resp(1,0,iresp)=fac*(a2+b2)
        resp(1,1,iresp)=fac*epsiln*(a2-b2)
        resp(2,0,iresp)=fac*(a3+b3)
        resp(2,1,iresp)=fac*epsiln*(a3-b3)
        resp(3,0,iresp)=fac*real(c+d)
        resp(3,1,iresp)=fac*epsiln*real(c-d)
        resp(4,0,iresp)=-fac*aimag(c+d)
        resp(4,1,iresp)=-fac*epsiln*aimag(c-d)
c
        resp(5,2,iresp)=-fac*real(alpha1+beta1)
        resp(5,3,iresp)=-fac*epsiln*aimag(alpha1+beta1)
        resp(6,2,iresp)=fac*aimag(alpha1-beta1)
        resp(6,3,iresp)=-fac*epsiln*real(alpha1-beta1)
        resp(7,2,iresp)=-fac*real(alpha2+beta2)
        resp(7,3,iresp)=-fac*epsiln*aimag(alpha2+beta2)
        resp(8,2,iresp)=fac*aimag(alpha2-beta2)
        resp(8,3,iresp)=-fac*epsiln*real(alpha2-beta2)
c
  100 CONTINUE
c
c                   replace rt+ and rt- with rt and rtt
c
      DO 300 inuc=0,3
        DO 200 ideu=0,8
          rt=0.5*(resp(ideu,inuc,2)+resp(ideu,inuc,3))
          rtt=0.5*(resp(ideu,inuc,2)-resp(ideu,inuc,3))
          resp(ideu,inuc,2)=rt
          resp(ideu,inuc,3)=rtt
  200   CONTINUE
  300 CONTINUE
c
      RETURN
      END
c
c_______________________________________________________________________
c
c*****TABL12
c
c     This subroutine reproduces the response functions as given in
c     table 12 of Gross and Dmitrasinovic.
c
      SUBROUTINE TABL12(g,eta,kappa2,resp)
c
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      DOUBLE PRECISION kappa2
c
      DIMENSION resp(0:8,0:3,9)
c
      COMPLEX g(18),a1,a2,a3,b1,b2,b3,c1,c2,d1,d2,e1,e2,e3,e4,
     &        f1,f2,f3,f4
c
      DO 100 ioff=0,12,12
        IF(ioff.EQ.0) THEN
          fac1=-4./3.*kappa2
          fac2=-fac1
          i1=6
          i2=9
        ELSE
          fac1=4./3.*kappa2*eta
          fac2=fac1
          i1=7
          i2=8
        ENDIF
c
        a1=conjg(g(8))*g(2+ioff)+conjg(g(10))*g(4+ioff)+conjg(g(12))*
     &     g(6+ioff)
        a2=conjg(g(8))*g(2+ioff)-conjg(g(12))*g(6+ioff)
        a3=conjg(g(8))*g(2+ioff)-2.*conjg(g(10))*g(4+ioff)+conjg(g(12))*
     &     g(6+ioff)
c
        b1=conjg(g(7))*g(1+ioff)+conjg(g(9))*g(3+ioff)+conjg(g(11))*
     &     g(5+ioff)
        b2=conjg(g(7))*g(1+ioff)-conjg(g(11))*g(5+ioff)
        b3=conjg(g(7))*g(1+ioff)-2.*conjg(g(9))*g(3+ioff)+conjg(g(11))*
     &    g(5+ioff)
c
        c1=conjg(g(8))*g(6+ioff)+conjg(g(12))*g(2+ioff)
        c2=conjg(g(8))*g(6+ioff)-conjg(g(12))*g(2+ioff)
c
        d1=conjg(g(7))*g(5+ioff)+conjg(g(11))*g(1+ioff)
        d2=conjg(g(7))*g(5+ioff)-conjg(g(11))*g(1+ioff)
c
        e1=conjg(g(7)+g(11))*g(4+ioff)+conjg(g(9))*(g(6+ioff)+g(2+ioff))
        e2=conjg(g(7)-g(11))*g(4+ioff)+conjg(g(9))*(g(6+ioff)-g(2+ioff))
        e3=conjg(g(7)-g(11))*g(4+ioff)-conjg(g(9))*(g(6+ioff)-g(2+ioff))
        e4=conjg(g(7)+g(11))*g(4+ioff)-conjg(g(9))*(g(6+ioff)+g(2+ioff))
c
        f1=conjg(g(8)+g(12))*g(3+ioff)+conjg(g(10))*(g(5+ioff)+
     &     g(1+ioff))
        f2=conjg(g(8)-g(12))*g(3+ioff)+conjg(g(10))*(g(5+ioff)-
     &     g(1+ioff))
        f3=conjg(g(8)-g(12))*g(3+ioff)-conjg(g(10))*(g(5+ioff)-
     &     g(1+ioff))
        f4=conjg(g(8)+g(12))*g(3+ioff)-conjg(g(10))*(g(5+ioff)+
     &     g(1+ioff))
c
        resp(0,2,i1)=fac1*aimag(a1+b1)
        resp(0,3,i1)=-fac1*real(a1-b1)
        resp(1,2,i1)=fac1*aimag(a2+b2)
        resp(1,3,i1)=-fac1*real(a2-b2)
        resp(2,2,i1)=fac1*aimag(a3+b3)
        resp(2,3,i1)=-fac1*real(a3-b3)
        resp(3,2,i1)=fac1*aimag(c1+d1)
        resp(3,3,i1)=-fac1*real(c1-d1)
        resp(4,2,i1)=fac1*real(c2+d2)
        resp(4,3,i1)=fac1*aimag(c2-d2)
c
        resp(5,0,i1)=-fac1*aimag(e1+f1)
        resp(5,1,i1)=-fac1*aimag(e1-f1)
        resp(6,0,i1)=-fac1*real(e2+f2)
        resp(6,1,i1)=-fac1*real(e2-f2)
        resp(7,0,i1)=-fac1*aimag(e3+f3)
        resp(7,1,i1)=-fac1*aimag(e3-f3)
        resp(8,0,i1)=-fac1*real(e4+f4)
        resp(8,1,i1)=-fac1*real(e4-f4)
c
        resp(0,2,i2)=fac2*real(a1+b1)
        resp(0,3,i2)=fac2*aimag(a1-b1)
        resp(1,2,i2)=fac2*real(a2+b2)
        resp(1,3,i2)=fac2*aimag(a2-b2)
        resp(2,2,i2)=fac2*real(a3+b3)
        resp(2,3,i2)=fac2*aimag(a3-b3)
        resp(3,2,i2)=fac2*real(c1+d1)
        resp(3,3,i2)=fac2*aimag(c1-d1)
        resp(4,2,i2)=-fac2*aimag(c2+d2)
        resp(4,3,i2)=fac2*real(c2-d2)
c
        resp(5,0,i2)=-fac2*real(e1+f1)
        resp(5,1,i2)=-fac2*real(e1-f1)
        resp(6,0,i2)=fac2*aimag(e2+f2)
        resp(6,1,i2)=fac2*aimag(e2-f2)
        resp(7,0,i2)=-fac2*real(e3+f3)
        resp(7,1,i2)=-fac2*real(e3-f3)
        resp(8,0,i2)=fac2*aimag(e4+f4)
        resp(8,1,i2)=fac2*aimag(e4-f4)
c
  100 CONTINUE
c
      RETURN
      END
c
c_______________________________________________________________________
c
c*****trnsvr
c
c     This subroutine constructs the mixed helity-transversity amplitudes
c     necessary to calculate the response functions from the helicity
c     amplitudes.
c
c
c     INPUT
c          curme(lnuc1,lnuc2,lgam,ldeu).........helicity matrix elements
c
c          lnuc1(lnuc2)=0................positive helicity for nucleon 1(2)
c                      =1................negative helicity for nucleon 1(2)
c          lgam=0,1......................photon helicity (particle 1)
c          ldeu=-1,0,1...................deuteron helicity (particle 2)
c
c      OUTPUT
c          g(i)..........................matrix elements in mixed basis
c          i=1,18........................index of matrix element
c
c
      SUBROUTINE trnsvr(curme,g)
c
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
c
      COMPLEX g(18),f(18),curme(0:1,0:1,0:1,-1:1),ip,im,rp,rm,
     &        op,om,ap,am,zz,gf1(12,12),gf2(13:18,13:18)
c
      PARAMETER(ip=(0.,1.),rp=(1.414213562373095,0.),op=(1.,0.),
     &          ap=ip*rp,im=-ip,rm=-rp,om=-op,am=-ap,zz=(0.,0.))
c
      DATA gf1/ip,ip,am,ap,im,im,om,op,rp,rp,op,om,
     &         ip,ip,ap,am,im,im,op,om,rp,rp,om,op,
     &         rp,rp,zz,zz,rp,rp,ap,am,zz,zz,ap,am,
     &         rp,rp,zz,zz,rp,rp,am,ap,zz,zz,am,ap,
     &         im,im,am,ap,ip,ip,op,om,rp,rp,om,op,
     &         im,im,ap,am,ip,ip,om,op,rp,rp,op,om,
     &         op,om,rp,rp,om,op,ip,ip,ap,am,im,im,
     &         om,op,rp,rp,op,om,ip,ip,am,ap,im,im,
     &         am,ap,zz,zz,am,ap,rp,rp,zz,zz,rp,rp,
     &         ap,am,zz,zz,ap,am,rp,rp,zz,zz,rp,rp,
     &         om,op,rp,rp,op,om,im,im,ap,am,ip,ip,
     &         op,om,rp,rp,om,op,im,im,am,ap,ip,ip/
c
      DATA gf2/ip,ip,am,ap,im,im,
     &         rp,rp,zz,zz,rp,rp,
     &         im,im,am,ap,ip,ip,
     &         om,op,rp,rp,op,om,
     &         ap,am,zz,zz,ap,am,
     &         op,om,rp,rp,om,op/
c
c
c               extract 18 necessary helicity matrix elements
c
      DO 200 ldeu=1,-1,-1
        f(14-ldeu)=curme(0,0,0,ldeu)
        f(17-ldeu)=curme(1,0,0,ldeu)
        DO 100 lnuc=0,1
          f(3+lnuc-2*ldeu)=curme(lnuc,lnuc,1,ldeu)
          f(9+lnuc-2*ldeu)=curme(lnuc,1-lnuc,1,ldeu)
  100   CONTINUE
  200 CONTINUE
c
c                              initialize g
c
      DO 250 i=1,18
        g(i)=(0.,0.)
  250 CONTINUE
c
c                 convert from helicity to mixed basis
c
      DO 400 i=1,12
        DO 300 j=1,12
          g(j)=g(j)+0.25*gf1(j,i)*f(i)
  300   CONTINUE
  400 CONTINUE
c
      DO 600 i=13,18
        DO 500 j=13,18
          g(j)=g(j)+0.5*gf2(j,i)*f(i)
  500   CONTINUE
  600 CONTINUE
c
CXXX      DO 700 i=1,18
CXXX        WRITE(10,*)i,f(i),g(i)
CXXX  700 CONTINUE
c
      RETURN
      END
