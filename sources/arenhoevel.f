C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE AREN_READ
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    06-Nov-2000
C       PURPOSE: Read in Arenhoevel response function files
C
C       MODIFICATIONS:
C            25-May-2001  PEU
C                Now get response function files from directory
C                associated with environment variable mceep_aren
C
C       The perl script aren_response.prl must first be used
C       to process each of Arenhoevel's files.
C
C       Unpolarized responses:  rl, rt, rlt, rtt
C
C       Each response is a five-index array:
C            rl(i,j,k,l,m), rt(i,j,k,l,m), ...
C
C
C          i:    Index for electron incident energy  (e0)
C          j:    Index for electron scattered energy (ef)
C          k:    Index for scattering angle          (th)
C          l:    Index for theta_CM                  (cm)
C          m:    Index for type of calculation (see below).
C
C          m=1:  PWBA          (NR)
C          m=2:  NORMAL        (NR)
C          m=3:  NORMAL+MEC    (NR)
C          m=4:  NORMAL+MEC+IC (NR)
C          m=5:  PWBA          (RC)
C          m=6:  NORMAL        (RC)
C          m=7:  NORMAL+MEC+IC (RC)
C
C                where NR = non-relativistic
C                      RC = relativisitic corrections included
C
C -----------------------------------------------------------------------
C
      SUBROUTINE AREN_READ
C
      IMPLICIT NONE
C
      INTEGER          ind(4),n,i,j,k,l,m,n_aren_dir
      DOUBLE PRECISION pi,kinvar(3),x,y
      CHARACTER*80     afile,dummy
      CHARACTER*200    afile_full
      CHARACTER*512    aren_dir_tmp
      CHARACTER*100    aren_dir
C
      PARAMETER (pi=3.14159265358979324D0) 
C
      INCLUDE 'aren_respon.cmn'
      INCLUDE 'aren_ranges.cmn'
C
C -----------------------------------------------------------------------
C     Open the file which lists all the Arenhoevel files.
C
C     The first four lines of this file must each contain
C     the low and high value and number of grid points for
C     a given variable in the order:
C
C         Line 1 -- e0 (MeV):   lo(1),hi(1),npts(1)
C         Line 2 -- ef (MeV):   lo(2),hi(2),npts(2)
C         Line 3 -- th (deg):   lo(3),hi(3),npts(3)
C         Line 4 -- cm (deg):   This line is special since
C                               it seems that Arenhoevel does not,
C                               by default, produce uniform grid
C                               in Theta_cm!
C
C         Line 4:  cm_1,cm_2,cm_3,cm_4,cm_d1,cm_d2,cm_d3,cm_npts
C
C                  For CM angle, it is assumed that there are
C                  three intervals each with uniform spacing
C                  (if you were clever enough to request uniform
C                  spacing then just set cm_d1=cm_d2=cm_d3
C                  and set cm_1/cm_4 to lo and hi values and
C                  set cm_2/cm_3 to be uniformly spaced):
C
C                      cm_1-cm_2, spacing of cm_d1
C                      cm_2-cm_3, spacing of cm_d2
C                      cm_3-cm_4, spacing of cm_d3
C                      cm_npts = total # of points in CM angle
C
C     Degree units are converted in this routine to radians for
C     use by the interpolator.
C -----------------------------------------------------------------------
C
      open(unit=11,name='arenhoevel_files.list',status='old',
     #     form='formatted')
C
      do n=1,3
         read(11,*) lo(n),hi(n),npts(n)  ! grid ranges, #pts
         del(n) = (hi(n)-lo(n))/dfloat(npts(n)-1)  ! get step sizes
      enddo
C
C -----------------------------------------------------------------------
C     Deal with possibly non-uniform grid in Theta_CM.
C -----------------------------------------------------------------------
C
      read(11,*) cm_1,cm_2,cm_3,cm_4,cm_d1,cm_d2,cm_d3,cm_npts
C
C -----------------------------------------------------------------------
C     Stuff the 4th element of lo, hi, npts and del arrays assuming
C     that the grid IS uniform.  I do this since others may want to
C     modify the code for uniform grid.  For non-uniform grid, I won't
C     use these values (or else there will clearly be trouble!).
C -----------------------------------------------------------------------
C
      lo(4) = cm_1
      hi(4) = cm_4
      npts(4) = cm_npts
      del(4)  = (hi(4)-lo(4))/dfloat(npts(4)-1)
C
C -----------------------------------------------------------------------
C     Convert units for grid ranges.
C -----------------------------------------------------------------------
C
      do n=3,4         ! Convert angles only (deg-->rad)
         lo(n)  = lo(n) * pi/180.d0
         hi(n)  = hi(n) * pi/180.d0
         del(n) = del(n)* pi/180.d0
      enddo
C
C -----------------------------------------------------------------------
C     Convert units for possibly non-uniform Theta_CM grid.
C -----------------------------------------------------------------------
C
      cm_1  = cm_1  * pi/180.d0
      cm_2  = cm_2  * pi/180.d0
      cm_3  = cm_3  * pi/180.d0
      cm_4  = cm_4  * pi/180.d0
      cm_d1 = cm_d1 * pi/180.d0
      cm_d2 = cm_d2 * pi/180.d0
      cm_d3 = cm_d3 * pi/180.d0
C
C ---------------------------------------------------------------------
C       Decode value of $mceep_aren.
C ---------------------------------------------------------------------
C
      CALL GETENV('mceep_aren',aren_dir_tmp)
      CALL SQUEEZE(aren_dir_tmp,aren_dir_tmp,n_aren_dir)
      aren_dir = aren_dir_tmp(1:n_aren_dir)
C
C -----------------------------------------------------------------------
C     Open the response function file.
C -----------------------------------------------------------------------
C
 1    read(11,*,end=99) afile
      afile_full = aren_dir(1:n_aren_dir)//'/'//afile
      open(unit=12,file=afile_full,status='old',form='formatted')
C
C -----------------------------------------------------------------------
C     Determine the type of calculation for this file and
C     associate with index m.  It is assumed that the filename
C     begins with an integer and that this integer corresponds with
C     the calculation type (see comment header in this file).
C     This should automatically be done by the perl script
C     aren_response.prl.
C -----------------------------------------------------------------------
C
      if(afile(1:1) .eq. '1') then
         m = 1
      elseif(afile(1:1) .eq. '2') then
         m = 2
      elseif(afile(1:1) .eq. '3') then
         m = 3
      elseif(afile(1:1) .eq. '4') then
         m = 4
      elseif(afile(1:1) .eq. '5') then
         m = 5
      elseif(afile(1:1) .eq. '6') then
         m = 6
      elseif(afile(1:1) .eq. '7') then
         m = 7
      else
         write(6,*) ' Arenhoevel file nomenclature error '
         write(6,*) ' Stopping execution of MCEEP '
         stop
      endif
C
C -----------------------------------------------------------------------
C     Determine the indices for e0, ef and theta for this file.
C -----------------------------------------------------------------------
C
      read(12,*) (kinvar(n),n=1,3)         ! e0, ef, th (respectively)
      kinvar(3) = kinvar(3) * pi/180.d0    ! convert th from deg-->rad
C
      do n=1,3
         ind(n) = nint((kinvar(n)-lo(n))/del(n)) + 1
         if(ind(n) .lt. 1 .or. ind(n) .gt. npts(n)) then
            write(6,*) ' Trying to index outside of Arenhoevel arrays '
            write(6,*) ' Stopping execution of MCEEP '
            stop
         endif
      enddo
C
      i = ind(1)    ! e0
      j = ind(2)    ! ef
      k = ind(3)    ! th
C
C -----------------------------------------------------------------------
C     Read the response functions.
C -----------------------------------------------------------------------
C
      read(12,*) dummy
      read(12,*) dummy
C
      do l=1,npts(4)       ! x and y are treated as dummies here
         read(12,*) x,y,  rl(i,j,k,l,m),  rt(i,j,k,l,m),
     #                   rlt(i,j,k,l,m), rtt(i,j,k,l,m)
      enddo
C
      close(unit=12)
C
      goto 1               ! Open and read next Arenhoevel file
C
 99   close(unit=11)       ! Done, so close the main file
C
 999  RETURN
      END
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C      Subroutine AREN_SIGMA
C
C      Author:  P.E. Ulmer
C      Date:    07-Nov-2000
C
C      Purpose:   Compute cross section given Arenhoevel's
C                 unpolarized response functions.  The kinematics
C                 are calculated by KINEM prior to this routine.
C
C      Units:
C                 rho_ij:        fm^-2          (density matrix)
C                 rij:           fm             (response funcs)
C                 qmu2:          (MeV/c)^2      (metric: >0)
C                 sigma          fm^2/MeV-sr^2  (cross section)
C
C      For the various expressions see Arenhovel's paper:
C              Phys. Rev. C43, 1022 (1991).
C
C      The Jacobian, dOmega_np^cm/dOmega_p, is calculated by the
C      routine DEUT_JACOB_CM.  This is required since MCEEP
C      wants all the cross section differentials in the LAB system
C      whereas Arenhoevel's cross section is differential in the
C      proton CoM solid angle.
C
C      The indices for the variables e0, ef, theta_scat and theta_cm,
C      must be passed to this routine.  A set of cross sections
C      (one for each type of calculation:  PWBA, Normal, ...) is
C      generated.
C
C ----------------------------------------------------------------------
C
      SUBROUTINE AREN_SIGMA(e0_ind,ef_ind,th_ind,cm_ind,e0,ef,th,cm,
     #           omega,w,qvec,qmu2,thpq,pr,phi,sigma)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION pi,hbarc,alpha,m_deut
      DOUBLE PRECISION e0,ef,th,cm,omega,w,qvec,qmu2,thpq,pr,phi
      DOUBLE PRECISION factor,qcm,beta,eta,qrat
      DOUBLE PRECISION rho_l,rho_t,rho_lt,rho_tt,cosphi,cos2phi
      DOUBLE PRECISION deut_jacob_cm,sigma(7)
C
      INTEGER e0_ind,ef_ind,th_ind,cm_ind,i,j,k,l,m
C
      PARAMETER (pi=3.14159265359D0)
      PARAMETER (hbarc=197.327D0)       ! in MeV-fm
      PARAMETER (alpha=7.297353D-3)
      PARAMETER (m_deut=1875.613D0)     ! deuteron mass in MeV
C
      INCLUDE 'aren_respon.cmn'
C
C
C ----------------------------------------------------------------------
C      Overall multiplicative factor for cross section.
C      The factor of hbarc is needed since the response functions
C      are in fm and I calculate the density matrix in MeV^2.
C      "factor" then produces a cross section in fm^2/MeV-sr^2.
C ----------------------------------------------------------------------
C
      factor = (alpha/(6.d0*pi*pi))*(ef/(e0*qmu2*qmu2))
     #         * hbarc
C
      qcm   = qvec*m_deut/w             ! 3-vector q in CoM system
      beta  = qvec/qcm                  ! boost from Lab --> CM
      eta   = (tan(th/2.d0))**2
      qrat  = qmu2/qvec**2
C
C ----------------------------------------------------------------------
C      Arenhoevel kinematic factors (in MeV^2).
C ----------------------------------------------------------------------
C
      rho_l  = beta**2 * qmu2 * qrat**2 / (2.d0*eta)
      rho_t  = 0.5d0 * qmu2 * (1.d0 + qrat/(2.d0*eta))
      rho_lt = beta * qmu2 * (qrat/eta) * sqrt((qrat+eta)/8.d0)
      rho_tt = -qmu2 * (qrat/(4.d0*eta))
C
      cosphi  = cos(phi)
      cos2phi = cos(2.d0*phi)
C
C ----------------------------------------------------------------------
C      Set indices for response functions.
C ----------------------------------------------------------------------
C
      i = e0_ind
      j = ef_ind
      k = th_ind
      l = cm_ind
C
C ----------------------------------------------------------------------
C      Assemble the cross sections for each model (m=1-7).
C      Response functions are in fm.
C ----------------------------------------------------------------------
C
      do m=1,7
         sigma(m) = factor * (  rho_l  *  rl(i,j,k,l,m)
     #                        + rho_t  *  rt(i,j,k,l,m)
     #                        + rho_lt * rlt(i,j,k,l,m) * cosphi
     #                        + rho_tt * rtt(i,j,k,l,m) * cos2phi )
C
         sigma(m) = sigma(m) * deut_jacob_cm(qvec,omega,thpq,cm,pr)
      enddo
C
      return
      end
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C      DOUBLE PRECISION FUNCTION DEUT_JACOB_CM
C
C      Purpose:
C
C        Calculate Jacobian from Center-of-Mass proton solid angle
C        to LAB for d(e,e'p)n.
C
C      This routine requires (calculated by KINEM):
C
C        QMAG, OMEGA, THETA_PQ, THETA_CM, PRMAG
C ----------------------------------------------------------------------
C
       DOUBLE PRECISION FUNCTION DEUT_JACOB_CM(QMAG,OMEGA,THETA_PQ,
     #                                         THETA_CM,PRMAG)
C
       IMPLICIT NONE
C
       DOUBLE PRECISION M_DEUT,M_PROT,M_NEUT
       DOUBLE PRECISION QMAG,OMEGA,THETA_PQ,THETA_CM,PRMAG
       DOUBLE PRECISION CTHPQ,STHPQ,E_NEUT,E_PROT,P_PROT
       DOUBLE PRECISION E_NP,M_NP,SQRT_ETA,SQRT_OP_ETA,BETA,GAMMA
       DOUBLE PRECISION P_PARA_CM,P_PERP_CM,P_PROT_CM
C
       PARAMETER (M_DEUT = 1875.613D0)       !Deuteron mass
       PARAMETER (M_PROT =  938.272D0)       !Proton mass
       PARAMETER (M_NEUT =  939.566D0)       !Neutron mass
C
       CTHPQ = COS(THETA_PQ)
       STHPQ = SIN(THETA_PQ)
C
       E_NEUT = SQRT(PRMAG**2   + M_NEUT**2) !Neutron energy in lab
       E_PROT = OMEGA + M_DEUT - E_NEUT
       P_PROT = SQRT(E_PROT**2 - M_PROT**2)  !Proton momentum in lab
C
       E_NP = E_PROT + E_NEUT                !Hadronic energy in lab
       M_NP = SQRT(E_NP**2-QMAG**2)          !Hadronic invariant mass
C
       SQRT_ETA = QMAG/M_NP                  !sqrt(eta)
       SQRT_OP_ETA = E_NP/M_NP               !sqrt(1+eta)
C
       BETA  = QMAG/(OMEGA+M_DEUT)           !Center of mass velocity
       GAMMA = 1.D0/SQRT(1.D0-BETA**2)
C
       P_PARA_CM = GAMMA*(P_PROT*CTHPQ-BETA*E_PROT)
       P_PERP_CM = P_PROT*STHPQ
       P_PROT_CM = SQRT(P_PARA_CM**2+P_PERP_CM**2) !ejectile COM momentum
C
       DEUT_JACOB_CM = ABS((P_PROT/P_PROT_CM)*(1.D0/SQRT_OP_ETA)
     #        /(1.D0-SQRT_ETA*E_PROT*CTHPQ/(SQRT_OP_ETA*P_PROT)))
C
       RETURN
       END
C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE AREN_INTERP
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    09-Nov-2000
C       PURPOSE: Interpolate Arenhoevel response functions
C
C       METHOD:  Since cross section variations can be extreme,
C                I first calculate the ratio of the Arenhoevel
C                cross section to PWIA at all of the vertices
C                of a cell surrounding the kinematic point for
C                a MCEEP event.  I then interpolate linearly
C                on this ratio and finally multiply the result
C                by the PWIA cross section at the correct
C                kinematics.  For speed purposes I use the
C                simplest possible PWIA cross section.
C
C       Unpolarized responses:  rl, rt, rlt, rtt
C
C       Each response is a five-index array:
C            rl(i,j,k,l,m), rt(i,j,k,l,m), ...
C
C       See the routine aren_read for the definition of the
C       various indices of these response functions.
C
C -----------------------------------------------------------------------
C
      SUBROUTINE AREN_INTERP(fail_cell,sigma_aren)
C
      IMPLICIT NONE
C
      COMMON /KINVAR/  KI,Q,Q2,QMU2,EP,PF,THETA_EP,CTH_PQ,PR
      DOUBLE PRECISION KI(3),Q(3),Q2,QMU2,EP,PF(3),THETA_EP,CTH_PQ,PR(3)
C
      DOUBLE PRECISION siglpt,siglt,sigtt,spectral,pwia_0,pwia_n
      DOUBLE PRECISION kinvar(4),val_lo(4),f(4),t1,t2,t3,t4
      DOUBLE PRECISION e0_n,ef_n,th_n,cm_n
      DOUBLE PRECISION omega_n,w_n,qvec_n,q2_n,qmu2_n
      DOUBLE PRECISION ep_n,pr_n,thpq_n
      DOUBLE PRECISION aren_n(7),aren_ratio(7),sigma_aren(7)
C
      INTEGER e0_ind,ef_ind,th_ind,cm_ind
      INTEGER ind_lo(4),tmp_ind,i,j,k,l,m,n
C
      LOGICAL fail_cell
C
      INCLUDE 'var.cmn'
      INCLUDE 'aren_ranges.cmn'
C
      fail_cell = .false.         ! initialize
C
C -----------------------------------------------------------------------
C     Set up array of kinematic variables for this event.
C -----------------------------------------------------------------------
C
      kinvar(1) = e0_i
      kinvar(2) = pf_e_i
      kinvar(3) = tscat
      kinvar(4) = theta_cm
C
C -----------------------------------------------------------------------
C     Determine N-dimensional cell for this event.
C
C     Define it by the lower index in each dimension.
C     The theta_cm index has to be treated specially since there
C     may be non-uniform grid in this variable.
C -----------------------------------------------------------------------
C
      do n=1,3
         ind_lo(n) = 1 + int( (kinvar(n)-lo(n)) / del(n) )
         val_lo(n) = dfloat(ind_lo(n)-1)*del(n) + lo(n)
      enddo
C
      if(kinvar(4).ge.cm_3) then
         tmp_ind = int( (kinvar(4)-cm_3) / cm_d3 )
         ind_lo(4) = 1 + nint((cm_2-cm_1)/cm_d1)
     #                 + nint((cm_3-cm_2)/cm_d2)
     #                 + tmp_ind
         del(4) = cm_d3
         val_lo(4) = dfloat(tmp_ind)*cm_d3 + cm_3
      elseif(kinvar(4).lt.cm_3 .and. kinvar(4).ge.cm_2) then
         tmp_ind = int( (kinvar(4)-cm_2) / cm_d2 )
         ind_lo(4) = 1 + nint((cm_2-cm_1)/cm_d1)
     #                 + tmp_ind
         del(4) = cm_d2
         val_lo(4) = dfloat(tmp_ind)*cm_d2 + cm_2
      elseif(kinvar(4).lt.cm_2) then
         tmp_ind = int( (kinvar(4)-cm_1) / cm_d1 )
         ind_lo(4) = 1 + tmp_ind
         del(4) = cm_d1
         val_lo(4) = dfloat(tmp_ind)*cm_d1 + cm_1
      endif
C
      do n=1,4
         if(ind_lo(n).lt.1 .or. ind_lo(n).ge.npts(n)) then
            fail_cell = .true.  ! can't find cell surrounding pt
            goto 999
         endif
         f(n) = (kinvar(n)-val_lo(n)) / del(n)  ! fraction of step
      enddo
C
      do m=1,7
         aren_ratio(m) = 0.d0         ! initialize
      enddo
C
C -----------------------------------------------------------------------
C     Got the cell, so proceed.
C
C     First, get PWIA cross section for the actual event kinematics.
C     For now, I will use model 100 with spectral function=10.
C     (ISPEC_OPT=10 is set in SUBROUTINE PHYSICS).
C     Here, SPEC_FAC is set to 1.d0 since this is overall scale
C     for the cross sections and irrelevant for this method.
C     Also, ignore recoil factor (REC_FAC) since it won't vary much
C     over the N-cell.
C -----------------------------------------------------------------------
C
      call off_shell_d(pf_e_i,q2,qmu2,ep,prec,theta_pq,
     #             tscat,siglpt,siglt,sigtt)
      pwia_0 = (siglpt+siglt*cos(phi_x)+sigtt*cos(2.d0*phi_x))
     #                  *spectral(.true.,1.d0,prmag,miss_m)
C
C -----------------------------------------------------------------------
C     Now, get PWIA and Arenhoevel cross sections for each vertex
C     of the N-dimensional cell.
C -----------------------------------------------------------------------
C
      do i=0,1
         do j=0,1
            do k=0,1
               do l=0,1
C
                     e0_ind = ind_lo(1) + i
                     ef_ind = ind_lo(2) + j
                     th_ind = ind_lo(3) + k
                     cm_ind = ind_lo(4) + l
C
                     e0_n = val_lo(1) + dfloat(i)*del(1)
                     ef_n = val_lo(2) + dfloat(j)*del(2)
                     th_n = val_lo(3) + dfloat(k)*del(3)
                     cm_n = val_lo(4) + dfloat(l)*del(4)
C
                     call deutkin(e0_n,ef_n,th_n,cm_n,
     #                            omega_n,w_n,qvec_n,q2_n,qmu2_n,
     #                            ep_n,pr_n,thpq_n)
C
                     call off_shell_d(ef_n,q2_n,qmu2_n,ep_n,pr_n,
     #                                thpq_n,th_n,siglpt,siglt,sigtt)
                     pwia_n = ( siglpt + siglt*cos(phi_x)
     #                        + sigtt*cos(2.d0*phi_x))
     #                        *spectral(.true.,1.d0,pr_n,miss_m)
C
                     call aren_sigma(e0_ind,ef_ind,th_ind,cm_ind,
     #                               e0_n,ef_n,th_n,cm_n,omega_n,w_n,
     #                               qvec_n,qmu2_n,thpq_n,pr_n,phi_x,
     #                               aren_n)
C
C -----------------------------------------------------------------------
C                    Finally, do the interpolation.
C -----------------------------------------------------------------------
C
                     t1 = 1.d0 - f(1) - dfloat(i)*(1.d0-2.d0*f(1))
                     t2 = 1.d0 - f(2) - dfloat(j)*(1.d0-2.d0*f(2))
                     t3 = 1.d0 - f(3) - dfloat(k)*(1.d0-2.d0*f(3))
                     t4 = 1.d0 - f(4) - dfloat(l)*(1.d0-2.d0*f(4))
C
                     do m=1,7
                        aren_ratio(m) =  aren_ratio(m)
     #                                + (aren_n(m)/pwia_n)*t1*t2*t3*t4
                     enddo
C
               enddo
            enddo
         enddo
      enddo
C
C -----------------------------------------------------------------------
C     Scale the final result by PWIA at the event kinematics.
C -----------------------------------------------------------------------
C
      do m=1,7
         sigma_aren(m) = aren_ratio(m) * pwia_0
      enddo
C
 999  RETURN
      END
C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE DEUTKIN
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    09-Nov-2000
C
C       PURPOSE: Simple kinematics routine for d(e,e'p)n for
C                use with the Arenhoevel interpolation routine.
C                (It's only needed to get PWIA cross sections
C                at the N-cell vertices).
C
C       INPUT:    e0     electron initial energy (MeV)
C                 ef     electron final   energy (MeV)
C                 th     scattering angle (rad)
C                 cm     np CoM angle     (rad)
C
C       OUTPUT:   q2     qvec^2 (MeV/c)^2
C                 qmu2   Q^2    (MeV/c)^2
C                 ep     proton energy   in Lab (MeV)
C                 pn     neutron (i.e. recoil) momentum in Lab (MeV/c)
C                 thpq   proton angle wrt q in Lab (rad)
C
C -----------------------------------------------------------------------
C
      SUBROUTINE DEUTKIN(e0,ef,th,cm,omega,w,qvec,q2,qmu2,ep,pn,thpq)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION e0,ef,th,cm
      DOUBLE PRECISION m_deut,m_prot,m_neut
      DOUBLE PRECISION qmu2,omega,q2,qvec,w,ep_cm,pp_cm,beta,gamma
      DOUBLE PRECISION ep,pp2,pp,en,pn2,pn,cthpq,thpq,test_cthpq
C
      PARAMETER (m_deut = 1875.613D0)     ! Deuteron mass
      PARAMETER (m_prot =  938.272D0)     ! Proton mass
      PARAMETER (m_neut =  939.566D0)     ! Neutron mass
C
      qmu2  = 4.d0*e0*ef*(sin(th/2.d0))**2
      omega = e0-ef
      q2    = qmu2 + omega**2             ! q_vec^2
      qvec  = sqrt(q2)
      w     = sqrt((omega+m_deut)**2-q2)  ! Invariant mass of final state
C
      ep_cm = (w**2 + m_prot**2 - m_neut**2)/(2.d0*w) ! prot energy in CoM
      pp_cm = sqrt(ep_cm**2 - m_prot**2)              ! prot mom    in CoM
      beta  = qvec/(omega+m_deut)                     ! -CoM velocity
      gamma = 1.d0/sqrt(1.d0-beta**2)
      ep    = gamma*ep_cm + gamma*beta*pp_cm*cos(cm)  ! prot energy in Lab
      pp2   = ep**2 - m_prot**2
      pp    = sqrt(pp2)                               ! prot mom    in Lab
      en    = omega + m_deut - ep                     ! neut energy in Lab
      pn2   = en**2 - m_neut**2
      pn    = sqrt(pn2)                               ! neut mom    in Lab
C
      cthpq = (q2 + pp2 - pn2) / (2.d0*pp*qvec)
C
C -----------------------------------------------------------------------
C     Protect against roundoff.  This is necessary, since the grid
C     kinematics can easily involve theta_cm = 0 or 180.
C -----------------------------------------------------------------------
C
      test_cthpq = 1.d0 - abs(cthpq)
      if(test_cthpq .lt. 0.D0 .and. test_cthpq .gt. -1.d-12) then
         cthpq = cthpq/abs(cthpq)
      endif
      thpq  = acos(cthpq)                             ! theta_pq    in Lab
C
      RETURN
      END










      







