C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C     Program   SYSTERR:
C
C               Performs analysis of systematic errors related to
C               kinematical uncertainties.  Requires pre-existing
C               Ntuple containing Transport vectors at target for
C               each particle (these are MCEEP variables 
C               62, 64, 66, 72, 74, 76).
C
C               Reads central kinematics and physics model choices.
C               Angles are converted to radians here.
C               Also reads the Ntuple and calls required standard
C               MCEEP routines.
C
C               Output:  An Ntuple identical to the input one, except:
C                        1) It's stripped of any events which lead
C                           to bogus kinematics for any of the 
C                           variations of the 9 kinematic variables.
C                        2) 10 additional entries of cross sections -
C                           1 for the nominal kinematics and the other
C                           9 for shifts of each of the 9 variables
C                           in turn (with all others at nominal values).
C
C
C     P.E. Ulmer   02-JAN-2002
C
C     Modifications:
C
C     INPUT:
C            hbook_pre:    Prefix of input HBOOK Ntuple file (w/o .hbook)
C
C            ntuid:        Ntuple ID within hbookfile
C
C            nev_max:      Maximum # of events to process
C
C            e1:           Energy of beam (MeV)
C            phi1:         In-plane angle of beam (deg)
C            theta1:       Out-of-plane angle of beam (deg)
C            e2:           Energy of scattered electron (MeV)
C            phi2:         In-plane angle of scattered electron (deg)
C            theta2:       Out-of-plane angle of scattered electron (deg)
C            px:           Momentum of ejectile (MeV/c)
C            phix:         In-plane angle of ejectile (deg)
C            thetax:       Out-of-plane angle of ejectile (deg)
C
C            atarg:        Target A
C            ztarg:        Target Z
C            mx:           Mass of ejectile (MeV)
C
C
C      OUTPUT:   This code produces the Ntuple required by
C                systerr.kumac.  The actual summing of errors is
C                done in toterr.
C
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C
      IMPLICIT NONE
c
      character*80     infile,chtitl
      character*60     hbook_pre
      character*80     hbook_in,hbook_out,tmpfile
      character*8      chtag(100),chtag_out(100)
      CHARACTER*512    DAT_DIR_TMP
      CHARACTER*100    DAT_DIR
c
      double precision e1,phi1,theta1,e2,phi2,theta2
      double precision px,phix,thetax,atarg,ztarg,deg_to_rad
      double precision eb,phib,thetab
      double precision spec_fac,pol_beam
      double precision cph_e,sph_e,cth_e,sth_e
      double precision cph_p,sph_p,cth_p,sth_p
      double precision theta_bp(2),cthbp_e,sthbp_e,cthbp_p,sthbp_p
      double precision rotx(3,3),roty(3,3),rotz(3,3),rotxy(3,3)
      double precision rot_e(3,3),rot_p(3,3),roti_e(3,3),roti_p(3,3)
      double precision trpt_t_e(6),trpt_t_p(6)
      double precision ki_vec(3),kf_vec(3)
      double precision sigma_eep,double_dummy(10),asymmetry,off
      double precision offset(0:9,0:9)
      double precision ki(3),q(3),q2,qmu2,ep,pf(3),theta_ep,cth_pq
      double precision pr(3)
      double precision sig(0:9)
C
C     Helicity dependent and independent polarizations.
C     Here, I include the variables but they are not really used
C     at present.
C
      DOUBLE PRECISION POL_N,POL_T,POL_L,POL_N_S,POL_T_S,POL_L_S
      DOUBLE PRECISION POL_X, POL_Y, POL_Z
      DOUBLE PRECISION POL_N_RHI,POL_T_RHI,POL_L_RHI
      DOUBLE PRECISION POL_N_RHD,POL_T_RHD,POL_L_RHD
      DOUBLE PRECISION POL_N_SHI,POL_T_SHI,POL_L_SHI
      DOUBLE PRECISION POL_N_SHD,POL_T_SHD,POL_L_SHD
      DOUBLE PRECISION POL_X_HI, POL_Y_HI, POL_Z_HI
      DOUBLE PRECISION POL_X_HD, POL_Y_HD, POL_Z_HD
c
      real             blanc
      real             data_in(100),data_out(100)
      real             rlow(100),rhigh(100)
c
      logical          bound,elastic,match,fail_kin,fail_any
      logical          accept_check,multiphoton,proton
      logical          radpeaking,radfull,accept_fcn
c
      integer          nfail,i,j
      integer          nwpawc,iquest
      integer          irc,ierror,icycle,k,idn,imatch_e,imatch_h
      integer          nchar,nevents,iev,nev_max
      integer          id_e_th,id_e_ph,id_e_dp
      integer          id_h_th,id_h_ph,id_h_dp
      integer          ntuid,nvar,lun
      integer          N_DAT_DIR
c
      INCLUDE 'labcoord.cmn'
      INCLUDE 'var.cmn'
      INCLUDE 'masses.cmn'
      INCLUDE 'masses_dat.cmn'
c
      parameter (deg_to_rad=0.017453d0)
      parameter (nwpawc=1000000)
c
      COMMON /PARTICLE/ PROTON
      COMMON /KINVAR/  KI,Q,Q2,QMU2,EP,PF,THETA_EP,CTH_PQ,PR
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
      common /pawc/  blanc(nwpawc)
      common /quest/ iquest(100)
c
      data theta_bp /-90.d0,-90.0/  ! assume vert. bending (upward) spect.
      data spec_fac /1.d0/          ! set to one, since irrelevant here
      data pol_beam /0.d0/          ! set to zero, since irrelevant here
c
      data off /1.d-3/              ! offset for momenta/angles
                                     ! units:  momenta - fractional error
                                     !         angles  - radians
C
C ---------------------------------------------------------------------
C       Decode value of $mceep_dat.
C ---------------------------------------------------------------------
C
      CALL GETENV('mceep_dat',DAT_DIR_TMP)
      CALL SQUEEZE(DAT_DIR_TMP,DAT_DIR_TMP,N_DAT_DIR)
      DAT_DIR = DAT_DIR_TMP(1:N_DAT_DIR)
c
c ------------------------------------------------------------------
c     Write out warnings since this is NOT full MCEEP analysis
c     here.  If you want full MCEEP analysis, then:  run MCEEP!
c ------------------------------------------------------------------
c
      write(6,*) ' NOTE: This code will not run properly for '
      write(6,*) '       the radiative case.  So enter 0 when '
      write(6,*) '       prompted below. '
      write(6,*) '  '
c
c ------------------------------------------------------------------
c     Read input file for central kinematics, masses and model
c     choices.
c ------------------------------------------------------------------
c
      write(6,*) ' Enter input file > '
      read(5,*) infile
      open(unit=1,file=infile,status='old',form='formatted')
c
      read(1,*) hbook_pre
      read(1,*) ntuid
      read(1,*) nev_max
      read(1,*) e1,phi1,theta1,e2,phi2,theta2,px,phix,thetax
      read(1,*) atarg,ztarg,eject_mass
c
 10   close(unit=1)
c
      call phys_choice(elastic,accept_check,bound,proton,
     #                 radpeaking,multiphoton,radfull,accept_fcn)
c
c ------------------------------------------------------------------
c     Convert angles to radians.
c ------------------------------------------------------------------
c
      phi1   = phi1   * deg_to_rad
      theta1 = theta1 * deg_to_rad
      phi2   = phi2   * deg_to_rad
      theta2 = theta2 * deg_to_rad
      phix   = phix   * deg_to_rad
      thetax = thetax * deg_to_rad
c
      do i=1,2
         theta_bp(i) = theta_bp(i) * deg_to_rad
      enddo
c
c ------------------------------------------------------------------
c     Get sines and cosines of angles.
c ------------------------------------------------------------------
c
      cph_e = cos(phi2)
      sph_e = sin(phi2)
      cth_e = cos(theta2)
      sth_e = sin(theta2)
c
      cph_p = cos(phix)
      sph_p = sin(phix)
      cth_p = cos(thetax)
      sth_p = sin(thetax)
c
      cthbp_e = cos(theta_bp(1))
      sthbp_e = sin(theta_bp(1))
c
      cthbp_p = cos(theta_bp(2))
      sthbp_p = sin(theta_bp(2))
c
c ------------------------------------------------------------------
c     Get matrices which rotate from LAB frame to spectrometer
c     Transport frame at the target.  Actually, here only the
c     inverse matrices are required since we directly read in the
c     Transport vector per event.
c ------------------------------------------------------------------
c
      call rotate_y(cph_e, sph_e,roty)
      call rotate_x(cth_e,-sth_e,rotx)
      call rotate_z(cthbp_e,sthbp_e,rotz)
      call mat_mult(3,3,3,rotx,roty,rotxy)
      call mat_mult(3,3,3,rotz,rotxy,rot_e)
c
      call rotate_y(cph_p, sph_p,roty)
      call rotate_x(cth_p,-sth_p,rotx)
      call rotate_z(cthbp_p,sthbp_p,rotz)
      call mat_mult(3,3,3,rotx,roty,rotxy)
      call mat_mult(3,3,3,rotz,rotxy,rot_p)
c
c ---------------------------------------------------------------------
c     Calculate inverse matrices:  Transport --> Lab
c ---------------------------------------------------------------------
c
      do i=1,3
        do j=1,3
          roti_e(i,j) = rot_e(j,i)
          roti_p(i,j) = rot_p(j,i)
        enddo
      enddo
c
c ---------------------------------------------------------------------
c     Calculate the target mass.
c ---------------------------------------------------------------------
c
      IF(ATARG .EQ. 1.) THEN
        MNUC = 938.279D0        !1H
      ELSEIF(ATARG .EQ. 2.) THEN
        MNUC = 1875.628D0       !2H
      ELSEIF(ATARG .EQ. 3. .AND. ZTARG .EQ. 1.) THEN
        MNUC = 2808.943D0       !3H
      ELSEIF(ATARG .EQ. 3. .AND. ZTARG .EQ. 2.) THEN
        MNUC = 2808.413D0       !3He
      ELSEIF(ATARG .EQ. 4. .AND. ZTARG .EQ. 2.) THEN
        MNUC = 3727.403D0       !4He
      ELSEIF(ATARG .EQ. 12 .AND. ZTARG .EQ. 6.) THEN
        MNUC = 11177.95D0       !12C
      ELSE
        MNUC = ATARG*(931.494D0-0.511D0)
      ENDIF
c
c ---------------------------------------------------------------------
c     Set the positional Transport coordinates to zero since
c     we only need the cross sections.  This is done here, since
c     it is not event dependent.
c ---------------------------------------------------------------------
c
      trpt_t_e(1) = 0.d0
      trpt_t_e(3) = 0.d0
      trpt_t_e(5) = 0.d0
c
      trpt_t_p(1) = 0.d0
      trpt_t_p(3) = 0.d0
      trpt_t_p(5) = 0.d0
c
c ---------------------------------------------------------------------
c    Form input and output HBOOK file names.
c ---------------------------------------------------------------------
c 
      tmpfile = hbook_pre//'.hbook'
      call squeeze2(tmpfile,tmpfile,nchar)
      hbook_in = tmpfile(1:nchar)
c
      tmpfile = hbook_pre//'_err.hbook'
      call squeeze2(tmpfile,tmpfile,nchar)
      hbook_out = tmpfile(1:nchar)
c
      call hlimit(nwpawc)
      iquest(10) = 65000
c      
c ---------------------------------------------------------------------
c    Open input hbook file and determine parameters of selected 
c    Ntuple (i.e. variable names, etc.).
c ---------------------------------------------------------------------
c
      call hropen(10,'HBKINPUT',hbook_in,' ',1024,irc)
      if (irc.ne.0) then
         write(6,*) ' Problem with input file irc = ',irc
         stop
      endif
      call hcdir('//HBKINPUT',' ')
      call hgnpar(ntuid,'make_new_ntuple')
      call hnoent(ntuid,nevents)  
      call hgiven(ntuid,chtitl,nvar,chtag,rlow,rhigh)
      call hgiven(ntuid,chtitl,nvar,chtag,rlow,rhigh)  ! seem to need
                                                       ! another call
      write(6,*) ' Number of events: ',nevents
c
c ---------------------------------------------------------------------
c    Determine the indices for the needed variables.
c    Also, define the additional variable names for the output Ntuple.
c ---------------------------------------------------------------------
c
      imatch_e = 0
      imatch_h = 0
      do i=1,nvar
         if(match(chtag(i),'TH_E_TG ',8)) then
            id_e_th = i 
            imatch_e = imatch_e + 1
         elseif(match(chtag(i),'PH_E_TG ',8)) then
            id_e_ph  = i
            imatch_e = imatch_e + 1
         elseif(match(chtag(i),'DP_E_TG ',8)) then
            id_e_dp = i
            imatch_e = imatch_e + 1
         endif
c                                       
         if(match(chtag(i),'TH_P_TG ',8)) then
            id_h_th = i
            imatch_h = imatch_h + 1
         elseif(match(chtag(i),'PH_P_TG ',8)) then
            id_h_ph = i
            imatch_h = imatch_h + 1
         elseif(match(chtag(i),'DP_P_TG ',8)) then
            id_h_dp = i
            imatch_h = imatch_h + 1
         endif
c
         chtag_out(i) = chtag(i)
      enddo
c
      if(imatch_e .eq. 3) then
         write(6,*) ' Electron variables found in Ntuple '
      else
         write(6,*) ' Electron variables NOT   in Ntuple '
         stop
      endif
c
      if(imatch_h .eq. 3) then
         write(6,*) ' Hadron   variables found in Ntuple '
      else
         write(6,*) ' Hadron   variables NOT   in Ntuple '
         stop
      endif
c
      chtag_out(nvar+1)  = 'SIG_0   '
      chtag_out(nvar+2)  = 'SIG_E1P '
      chtag_out(nvar+3)  = 'SIG_PH1P'
      chtag_out(nvar+4)  = 'SIG_TH1P'
      chtag_out(nvar+5)  = 'SIG_E2P '
      chtag_out(nvar+6)  = 'SIG_PH2P'
      chtag_out(nvar+7)  = 'SIG_TH2P'
      chtag_out(nvar+8)  = 'SIG_PXP '
      chtag_out(nvar+9)  = 'SIG_PHXP'
      chtag_out(nvar+10) = 'SIG_THXP'
c
c ---------------------------------------------------------------------
c    Open the output hbook file and book the new Ntuple.
c ---------------------------------------------------------------------
c
      call hropen(20,'HBKOUTPUT',hbook_out,'NQ',4096,irc)
      if (irc.ne.0) then
          write(6,*) ' Problem with output file irc = ',irc
          stop
      endif
      call hcdir('//HBKOUTPUT',' ')
      call hbookn(100,'ntu_out',nvar+10,'HBKOUTPUT',1024,chtag_out)
c
c
c ---------------------------------------------------------------------
c    Set up offset matrix.  This is a diagonal matrix so that
c    the offsets are applied successively to single variables
c    (all others being at their nominal values) in the order:
c         e1,phi1,theta1,e2,phi2,theta2,px,phix,thetax
c    Actually, for simplicity the offsets to the outgoing
c    particles are applied at the level of the target Transport
c    coordinates.
c ---------------------------------------------------------------------
c
      do i=0,9
         do j=0,9
            if(i.eq.j) then
               offset(i,j) = off    ! diagonal elements
            else
               offset(i,j) = 0.d0   ! off-diagonal elements
            endif
         enddo
      enddo
      offset(0,0) = 0.d0            ! this is the nominal case
                                     ! so no offset for first pass;
                                     ! actually this element is never
                                     ! used, but it makes me feel good 
c
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c    BEGIN EVENT LOOP
c ---------------------------------------------------------------------
c
      nfail = 0
      do iev=1,nevents
c
         if(iev .gt. nev_max) goto 9999   ! reached maximum number
c
         if(dfloat(iev/1000)-dfloat(iev)/1000. .eq. 0.)
     #          write(6,101) iev
  101    format(' Processing event: ',I9)
c
         fail_any = .false.               ! initialize
c
c ---------------------------------------------------------------------
c    Copy existing Ntuple into new output Ntuple.
c ---------------------------------------------------------------------
c
         call hcdir('//HBKINPUT',' ')
         call hgnf(ntuid,iev,data_in,ierror)
         do k=1,nvar
            data_out(k)=data_in(k)
         enddo
c
c ---------------------------------------------------------------------
c     Apply offsets.  Note that the Transport units are mrad and %
c     but the offsets were given units of radians and fractional
c     quantities.  Therefore we need conversion factors here.
c     We also need a conversion for the beam energy offset to go
c     from fractional error to error in MeV.
c
c     Get theta, phi and delta from the Ntuple Transport vectors.  
c     The x, y and L coordinates were set to zero above since we only 
c     need the cross sections.
c ---------------------------------------------------------------------
c
         do k=0,9                          ! offset loop
c
            eb          = e1               + offset(1,k)*e1
            phib        = phi1             + offset(2,k) 
            thetab      = theta1           + offset(3,k)
c
            e0_i        = eb               ! some physics codes need this
c
            trpt_t_e(6) = data_in(id_e_dp) + offset(4,k)*100.d0
            trpt_t_e(4) = data_in(id_e_ph) + offset(5,k)*1000.d0
            trpt_t_e(2) = data_in(id_e_th) + offset(6,k)*1000.d0
c
            trpt_t_p(6) = data_in(id_h_dp) + offset(7,k)*100.d0
            trpt_t_p(4) = data_in(id_h_ph) + offset(8,k)*1000.d0
            trpt_t_p(2) = data_in(id_h_th) + offset(9,k)*1000.d0
c
c ---------------------------------------------------------------------
c     Get Lab coordinates from the Transport vectors and make the
c     3-vectors for beam and scattered electron.
c ---------------------------------------------------------------------
c
            call trpt_to_lab(trpt_t_e,roti_e,e2,qi_e)
            call trpt_to_lab(trpt_t_p,roti_p,px,qi_p)
c
            call v3make(eb,phib,thetab,ki_vec)          !beam 3-vector
            call v3make(qi_e(6),qi_e(4),qi_e(2),kf_vec) !scatt e- 3-vector
c
c ---------------------------------------------------------------------
c     Call MCEEP kinematics routine.
c ---------------------------------------------------------------------
c
            call kinem(eb,ki_vec,qi_e(6),kf_vec,qi_p(6),
     #              qi_p(4),qi_p(2),elastic,.false.,.false.,fail_kin)
            if(fail_kin) then
               fail_any = .true.
               nfail = nfail + 1
               goto 999             ! get next event
            endif
c
c ---------------------------------------------------------------------
c     Recalculate "recoil factor".  This is necessary since kinem
c     is always called for the continuum case here, since we don't
c     want to force any correlations in this analysis induced by a
c     fixed missing mass.  However, the cross section still needs to
c     include the correct recoil factor.
c ---------------------------------------------------------------------
c
            call frec(bound,pr,pf,omega,ep,pf_p_i,recfac)
c
c ---------------------------------------------------------------------
c     Call MCEEP cross section routine.  Many of the variables are
c     currently irrelevant, in particular all the polarization
c     variables.
c ---------------------------------------------------------------------
c
            call physics(elastic,bound,spec_fac,pol_beam,
     #             sigma_eep,double_dummy,asymmetry,
     #          POL_N,POL_T,POL_L,POL_N_S,POL_T_S,POL_L_S,
     #          POL_X,POL_Y,POL_Z,
     #          POL_N_RHI,POL_T_RHI,POL_L_RHI,
     #          POL_N_RHD,POL_T_RHD,POL_L_RHD,
     #          POL_N_SHI,POL_T_SHI,POL_L_SHI,
     #          POL_N_SHD,POL_T_SHD,POL_L_SHD,
     #          POL_X_HI, POL_Y_HI, POL_Z_HI,
     #          POL_X_HD, POL_Y_HD, POL_Z_HD)
c
c ---------------------------------------------------------------------
c     If any cross section is zero, then get next event.
c ---------------------------------------------------------------------
c
            if(sigma_eep .eq. 0.d0) then
               fail_any = .true.
               nfail = nfail + 1
               goto 999             ! get next event
            endif
c
c ---------------------------------------------------------------------
c     Load the cross sections into data_out - for the output Ntuple.
c ---------------------------------------------------------------------
c
            sig(k) = sigma_eep
            data_out(nvar+k+1) = sig(k)
c
         enddo                               ! end of offset loop
c
c ---------------------------------------------------------------------
c     Fill the new Ntuple, but only if ALL cross sections for this 
c     event are valid!
c ---------------------------------------------------------------------
c
         if(.not. fail_any) then
           call hcdir('//HBKOUTPUT',' ')
           call hfn(100,data_out)
         endif
c
 999  enddo                     ! END OF EVENT LOOP
c
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c     Close the output hbook file.
c ---------------------------------------------------------------------
c
 9999 call hcdir('//HBKOUTPUT',' ')
      call hrout(100,icycle,'HBKOUTPUT')
      call hrendc('HBKOUTPUT') 
c
      write(6,1000) nfail
 1000 format(' # Events failing kinematics: ',i9)
c
c
      write(6,*) ' You should now execute systerr.kumac from PAW '
      write(6,*) ' After that, execute toterr (Fortran code) '
c
      stop
      end
c
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c
c ---------------------------------------------------------------------
c     Logical Function MATCH
c ---------------------------------------------------------------------
c
      LOGICAL FUNCTION MATCH(S1,S2,LEN)
      CHARACTER*1 S1(LEN),S2(LEN)
      MATCH = .TRUE.
      DO I=1,LEN
        IF(S1(I).NE.S2(I))MATCH = .FALSE.
      ENDDO
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       Extracts non-blank characters of any character string.
C------------------------------------------------------------------------------
C
      SUBROUTINE SQUEEZE2(CHARIN,CHAROUT,NCHAR)
      IMPLICIT NONE
C
      INTEGER NCHAR,N,J
      CHARACTER*80 CHARIN,CHAROUT
      N = 0
      DO J=1,80
        IF(CHARIN(J:J).NE.' ') THEN
          N = N+1
          CHAROUT(N:N) = CHARIN(J:J)
        ENDIF
      ENDDO
      NCHAR = N
      RETURN
      END
