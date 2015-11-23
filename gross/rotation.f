c
c______________________________________________________________________
c
c*****rotone
c
c     This subroutine constructs the spin one rotation matrices
c
c
c       INPUT
c           theta.............polar angle
c
c       OUTPUT
c           d(lam,lamp).........spin rotation matrix
c
c                lam=-1,0,1.....spin projection
c                lamp=-1,0,1....spin projection
c
c
      SUBROUTINE rotone(theta,d1)
c
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
c
      DIMENSION d1(-1:1,-1:1)
c
      PARAMETER(root2=1.414213562373095,root2i=1./root2)
c
c          construct angular momentum rotation matrix d1 (lower case d)
c
      d1(0,0)=cos(theta)
      d1(1,1)=0.5*(1.+d1(0,0))
      d1(1,-1)=0.5*(1.-d1(0,0))
      d1(-1,0)=root2i*sin(theta)
      d1(1,0)=-d1(-1,0)
      d1(0,1)=d1(-1,0)
      d1(0,-1)=d1(1,0)
      d1(-1,-1)=d1(1,1)
      d1(-1,1)=d1(1,-1)
c
      RETURN
      END
c
c______________________________________________________________________
c
c*****rothlf
c
c     This subroutine constructs the spin one half rotation matrices
c
c
c       INPUT
c           theta.............polar angle
c
c       OUTPUT
c           d(lnuc,lnucp).........spin rotation matrix
c
c                lnuc=0.........positive nucleon helicity
c                    =1.........negative nucleon helicity
c                lnucp=0........positive nucleon helicity
c                     =1........negative nucleon helicity
c
c
      SUBROUTINE rothlf(theta,dhalf)
c
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
c
      DIMENSION dhalf(0:1,0:1)
c
c       construct angular momentum rotation matrix dhalf (lower case d)
c
      thetah=0.5*theta
c
      dhalf(0,0)=cos(thetah)
      dhalf(1,1)=dhalf(0,0)
      dhalf(1,0)=sin(thetah)
      dhalf(0,1)=-dhalf(1,0)
c
      RETURN
      END
