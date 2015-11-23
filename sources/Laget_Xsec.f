C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C    

      SUBROUTINE LAGET_XSEC(XSEC,P_BEA,G_MAS,G_ENE,P_TGS,P_PHI,
     +                      SIG_L,SIG_T,SIG_LT,SIG_TT)

C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
C
C  Author: E. Voutier
C  Date: October 2005
C  Purpose: Determine D(e,e'p) cross section at current kinematics 
C  
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C  
C  Calculate the D(e,e'p) cross section for a selected kinematics from an 
C  input data  grid via
C  the interpolation  procedure specified by the user via the LAGET_INTP  
C  parameter (1 for linear, 2 for logarythmic in recoil momentum).  
C 
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
C
C  INPUT variables : 
C
C        P_BEA       3-momentum of the electron beam (MeV/c)         
C        G_MAS       quadrimomentum transfer (MeV^2)
C        G_ENE       energy transfer (MeV)
C        P_TGS       proton angle in the center of mass frame (rd)
C        P_PHI       out-of-plane angle of the proton (rd)
C        SIG_L       array of the longitudinal cross sections
C        SIG_T       array of the transverse cross sections
C        SIG_LT      array of the longitudinal-transverse cross sections
C        SIG_TT      array of the transverse-transverse cross sections
C
C  Passed via COMMON :
C
C        LAGET_INTP  interpolation index
C
C  OUTPUT variables :
C 
C        XSEC    actual D(e,e'p) cross section (µb.MeV-1.sr-2)
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

      IMPLICIT NONE
      
      INTEGER*4 NQ2,NOM,NTA
      INTEGER   LAGET_PS_FAIL,LAGET_GRID_FAIL
      INTEGER*4 LAGET_INTP,LAGET_PWIA,LAGET_FSI,LAGET_MEC

      REAL*8 INTRPOL,INTRLOG    
            
      REAL*8 XSEC,P_BEA,G_MAS,G_ENE,P_TGS,P_PHI
       
      REAL*8  SIG_L(NQ2,NOM,NTA),SIG_T(NQ2,NOM,NTA),
     +        SIG_LT(NQ2,NOM,NTA),
     +        SIG_TT(NQ2,NOM,NTA)
      
      REAL*8 P_CMS,E_CMS,P_ARA,P_ERP,P_MOM,P_ENE
      REAL*8 S_CMS,W_CMS,G_AMA,G_ABE
      REAL*8 E_SCA,P_SCA,C_TET
      REAL*8 E_BEA
      REAL*8 G_MOM
      
      REAL*8 XM_E,XM_P,XM_N,XM_D,XMD,SPN,DPN,XME2,XMN2,XMP2,XMD2,SPN2
      REAL*8 EPS,EPSP,FLUX,JCOB,SIGL,SIGT,SIGLT,SIGTT
      REAL*8 XB_MIN,XB_MAX,OM_MIN,OM_MAX
      REAL*8 CFIN,PI
      
      COMMON / SOME_MAS / DPN,XMD,SPN2,XMD2,XMN2
      COMMON / SIGMAS   / SIGL,SIGT,SIGLT,SIGTT
      COMMON / GRID_PAR / NQ2,NOM,NTA
      COMMON / LAGET_CNT/ LAGET_PS_FAIL,LAGET_GRID_FAIL
      COMMON / LAGET_MOD/ LAGET_INTP,LAGET_PWIA,LAGET_FSI,LAGET_MEC
     
      DATA  XM_E / 510.99890D-03 /,  XM_P / 938.27200D+00 /, 
     +      XM_N / 939.56533D+00 /,  XM_D / 187.56128D+01 /
     
      DATA  CFIN / 137.03599D+00 /
     
      PI = DACOS( -1.D+00 )

*---- Particle masses and relevant combinations                             

      XMD = XM_D
*-    Dimension 1 combinations
      SPN = XM_N + XM_P
      DPN = XM_N - XM_P
*-    Dimension 2 combinations                                     
      XME2 = XM_E * XM_E                                            
      XMD2 = XM_D * XM_D                                                       
      XMP2 = XM_P * XM_P
      XMN2 = XM_N * XM_N
      SPN2 =  SPN * SPN                                                      
      
*---- Initialisations

      XSEC  = 0.D+00
      SIGL  = 0.D+00
      SIGT  = 0.D+00
      SIGLT = 0.D+00
      SIGTT = 0.D+00 
      
*---- Beam energy           
                                             
      E_BEA = DSQRT( P_BEA * P_BEA + XME2 )
     
*---- Phase space restriction

      XB_MIN = E_BEA * DSQRT(G_MAS) + 
     +         P_BEA*DSQRT( G_MAS + 4.D+00*XME2 )
      XB_MIN = DSQRT(G_MAS) * XB_MIN / XM_P / 
     +         ( 4.D+00*P_BEA*P_BEA - G_MAS )
      XB_MAX = XM_D * G_MAS / XM_P / ( G_MAS + SPN2 - XMD2 )
      OM_MIN = 0.5D+00 * G_MAS / XM_P / XB_MAX
      OM_MAX = 0.5D+00 * G_MAS / XM_P / XB_MIN
      IF( (G_ENE.LT.OM_MIN).OR.(G_ENE.GT.OM_MAX) ) THEN
         LAGET_PS_FAIL = LAGET_PS_FAIL + 1
         GO TO 1      
      ENDIF

*---- Kinematics of the current event  
                                                                              
*-    Scattered electron 
      E_SCA = E_BEA - G_ENE                                 ! Energy
      P_SCA = DSQRT( E_SCA * E_SCA - XME2 )                 ! Momentum
      C_TET = E_BEA * E_SCA - XME2 - 0.5D+00 * G_MAS 
      C_TET = C_TET / P_BEA / P_SCA                         ! Polar angle 
                                                              ! cosinus
*-    Virtual photon
      G_MOM = DSQRT( G_MAS + G_ENE * G_ENE )                ! Momentum
*-    Center of mass frame
      S_CMS = XMD2 - G_MAS + 2.D+00*XM_D*G_ENE              ! Invariant squared
                                                              ! mass
      W_CMS = DSQRT( S_CMS )                                ! Invariant mass
      G_AMA = ( G_ENE + XM_D ) / W_CMS                      ! Gamma
      G_ABE = G_MOM / W_CMS                                 ! Gamma * Beta
*-    Knocked-out proton       
      P_CMS = ( S_CMS - SPN2 ) * ( S_CMS - DPN * DPN ) / S_CMS
      P_CMS = 0.5D+00 * DSQRT( P_CMS )                      ! Momentum in the 
                                                              ! cms frame
      E_CMS = DSQRT( P_CMS * P_CMS + XMP2 )                 ! Energy in the 
                                                              ! cms frame
      P_ARA = G_ABE * E_CMS + G_AMA * P_CMS * DCOS( P_TGS )
      P_ERP = P_CMS * DSIN( P_TGS )
      P_MOM = DSQRT( P_ARA*P_ARA + P_ERP*P_ERP )            ! Momentum in the 
                                                              ! lab frame
      P_ENE = G_AMA * E_CMS + G_ABE * P_CMS * DCOS( P_TGS ) ! Energy in the 
                                                              ! lab frame
      
*---- Virtual photon polarization

      EPS = G_MOM * G_MOM * (1.D+00 - C_TET) / (1.D+00 + C_TET) / G_MAS
      EPS = 1.D+00 / ( 1.D+00 + EPS + EPS ) 
      EPSP = DSQRT( 0.5D+00 * G_MAS * EPS * (1.D+00 + EPS) ) / G_ENE

*---- Virtual photon flux

      FLUX = 0.5D+00 * E_SCA * G_MOM / PI / PI / E_BEA / G_MAS / CFIN
      FLUX = FLUX / ( 1.D+00 - EPS )

*---- Calculation of the center of mass to the lab frame jacobian

      JCOB = G_ENE + XM_D - ( P_ENE * G_MOM * P_ARA / P_MOM / P_MOM )
      JCOB = W_CMS * P_MOM / P_CMS / DABS( JCOB )

*---- Interpolation of the individual response functions

      IF( LAGET_INTP.EQ.1 ) THEN
*-       Linear interpolation in angle     
         SIGL  = INTRPOL( SIG_L,G_MAS,G_ENE,P_TGS)
         SIGT  = INTRPOL( SIG_T,G_MAS,G_ENE,P_TGS)
         SIGLT = INTRPOL(SIG_LT,G_MAS,G_ENE,P_TGS)
         SIGTT = INTRPOL(SIG_TT,G_MAS,G_ENE,P_TGS)
      ELSE
*-       Logarythmic interpolation in recoil momentum
         SIGL  = INTRLOG( SIG_L,G_MAS,G_ENE,P_TGS)
         SIGT  = INTRLOG( SIG_T,G_MAS,G_ENE,P_TGS)
         SIGLT = INTRLOG(SIG_LT,G_MAS,G_ENE,P_TGS)
         SIGTT = INTRLOG(SIG_TT,G_MAS,G_ENE,P_TGS)
      ENDIF
*-    Safety check to protect against overflow (NaN)      
      IF((SIGL.EQ.0.D+00).AND.(SIGT.EQ.0.D+00).AND.(SIGLT.EQ.0.D+00)
     +                   .AND.(SIGTT.EQ.0.D+00)) GO TO 1

*---- Differential cross section (µb.MeV-1.sr-2)
 
      XSEC = FLUX * JCOB 
     +            * ( 
     +              SIGT 
     +            + EPS  * ( SIGL + SIGTT * DCOS(P_PHI+P_PHI) )
     +            - EPSP *   SIGLT * DCOS(P_PHI) 
     +		    )
      
    1 RETURN
      END

C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C    

      SUBROUTINE GRID_LOAD(SIG_L,SIG_T,SIG_LT,SIG_TT)
      
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C  Author: E. Voutier
C  Date: October 2005
C  Purpose: Load the partial cross section data grid
C  
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C  Load into the common area the  response  functions  of the D(e,e'p) 
C  reaction determined  in 
C  the framework of Jean-Marc Laget's formalism over a grid  sampled  
C  in  Q2, omega  and tta_p (the proton angle in the center of mass frame) 
C  in the range:
C
C            0.05 GeV2 <= Q2 <= 5.00 GeV2   step = 0.05 GeV2 
C           0.01 GeV <= omega <= 5.00 GeV   step = 0.01 GeV
C                 0 dg <= tta_p <= 180 dg   step = 2 dg 
C
C  The physics options  are  specified  according  to the convention 
C  encoded in the parameters Ipwia, Ifsi,and Imec that are part of the 
C  filename.
C
C     000  neutron contribution only : Pwia 
C     001  neutron contribution only : Pwia + Mec
C     010  neutron contribution only : Pwia + Fsi
C     011  neutron contribution only : Pwia + Fsi + Mec 
C     100  proton  contribution only : Pwia 
C     101  proton  contribution only : Pwia + Mec
C     110  proton  contribution only : Pwia + Fsi
C     111  proton  contribution only : Pwia + Fsi + Mec
C     200  neutron + proton contrib. : Pwia
C     201  neutron + proton contrib. : Pwia + Mec 
C     210  neutron + proton contrib. : Pwia + Fsi 
C     211  neutron + proton contrib. : Pwia + Fsi + Mec
C
C  The physics option parameters are passed via COMMON :
C
C     LAGET_PWIA
C     LAGET_FSI
C     LAGET_MEC
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - -  -  -  -  -  -  -  -  -
C
C  OUTPUT variables :
C 
C        SIG_L   array of the longitudinal cross sections
C        SIG_T   array of the transverse cross sections
C        SIG_LT  array of the longitudinal-transverse cross sections
C        SIG_TT  array of the transverse-transverse cross sections
C
C  These photoproduction like cross sections are connected to the usual
C  response functions via simple kinematic factors.    
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - -  -  -  -  -  -  -  -  -

      IMPLICIT NONE

      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
      COMMON / LAGET_MOD/ LAGET_INTP,LAGET_PWIA,LAGET_FSI,LAGET_MEC
      
      CHARACTER*4   Q2VAL,OMVAL
      CHARACTER*1   PW,FS,MC
      CHARACTER*300 FNAME,TMPNAM
      CHARACTER*100 DAT_DIR
      CHARACTER*50  SUBDIR
      
      INTEGER*4 NQ2,NOM,NTA
      INTEGER*4 N_Q2,N_OM,N_TA
      INTEGER*4 I_Q2,I_OM
      INTEGER*4 NCHAR,N_DAT_DIR
      INTEGER*4 LAGET_INTP,LAGET_PWIA,LAGET_FSI,LAGET_MEC

      PARAMETER ( N_Q2 = 100 )
      PARAMETER ( N_OM = 500 )
      PARAMETER ( N_TA =  91 )
           
      REAL*8  SIG_L(N_Q2,N_OM,N_TA),  SIG_T(N_Q2,N_OM,N_TA),
     +       SIG_LT(N_Q2,N_OM,N_TA), SIG_TT(N_Q2,N_OM,N_TA)
      REAL*8 Q2_STEP,OM_STEP,TA_STEP
      
      REAL*8 G_MAS,G_ENE
          
      COMMON / GRID_DEU / Q2_STEP,OM_STEP,TA_STEP
      COMMON / GRID_PAR / NQ2,NOM,NTA
     
      NQ2 = N_Q2
      NOM = N_OM
      NTA = N_TA
      
*-    Grid parameters

      Q2_STEP = 5.D+06 / DFLOAT(N_Q2)
      OM_STEP = 5.D+03 / DFLOAT(N_OM)
      TA_STEP = DACOS(-1.D+00) / DFLOAT(N_TA-1)
      
      WRITE(PW,'(I1)') LAGET_PWIA
      WRITE(FS,'(I1)') LAGET_FSI
      WRITE(MC,'(I1)') LAGET_MEC 
            
*---- User selection of the physics grid 
  
*-    No FSI nor MEC
      IF(LAGET_FSI.EQ.0) THEN
         SUBDIR = 'deut_laget/pwia/'
*-    With FSI and MEC	 
      ELSEIF(LAGET_MEC.EQ.1) THEN
         SUBDIR = 'deut_laget/pful/'
*-    With FSI only	 
      ELSE
         SUBDIR = 'deut_laget/pfsi/'
      ENDIF

*---- Load the physics selected grid	  
  
      DO I_Q2=1,N_Q2
         G_MAS = Q2_STEP * DFLOAT(I_Q2)
         WRITE(Q2VAL,'(F4.2)') 1.D-06 * G_MAS 
         DO I_OM=1,N_OM
            G_ENE = OM_STEP * DFLOAT(I_OM)
            WRITE(OMVAL,'(I4)') IDINT(G_ENE)
            IF(IDINT(G_ENE).LT.10) THEN
               OMVAL = '000'//OMVAL(4:4)
            ELSEIF(IDINT(G_ENE).LT.100) THEN 
               OMVAL = '00'//OMVAL(3:4)
            ELSEIF(IDINT(G_ENE).LT.1000) THEN 
               OMVAL = '0'//OMVAL(2:4)
            ENDIF
            TMPNAM = DAT_DIR(1:N_DAT_DIR)//'/'//SUBDIR//
     +               'q2_'//Q2VAL//'_om_'//OMVAL//
     +               '_'//PW//FS//MC//'.dat'
            CALL SQUEEZE(TMPNAM,FNAME,nchar)
            CALL READ_DATA(I_Q2,I_OM,SIG_L,SIG_T,SIG_LT,SIG_TT,FNAME) 
         ENDDO
      ENDDO
                                
      RETURN
      END
      
C
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
C    

      SUBROUTINE READ_DATA(IX,IY,XL,XT,XLT,XTT,FNAME)
      
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C  Author: E. Voutier
C  Date: October 2005
C  Purpose: Open and Read a basic data file with specified name 
C  
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
C  
C  Utility  tool for opening the cross section data  base  of Jean-Marc Laget 
C  and load partial cross sections into specified arrays. 
C  
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C  INPUT variables : 
C
C        IX      quadrimomentum transfer index         
C        IY      energy transfer index
C        FNAME   data file name
C
C  OUTPUT variables :
C 
C        XL      array of the longitudinal cross sections
C        XT      array of the transverse cross sections
C        XLT     array of the longitudinal-transverse cross sections
C        XTT     array of the transverse-transverse cross sections
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

      IMPLICIT NONE
      
      CHARACTER*(*) FNAME
      
      INTEGER*4 IX,IY,IZ,NQ2,NOM,NTA
        
      REAL*8 XLT(NQ2,NOM,NTA),XTT(NQ2,NOM,NTA)
      REAL*8 XL(NQ2,NOM,NTA),XT(NQ2,NOM,NTA)
      
      COMMON / GRID_PAR / NQ2,NOM,NTA

      OPEN(10,FILE=FNAME,STATUS='OLD')
      
      DO 100 IZ = 1,NTA
  100 READ(10,'(4(2X,E16.9))') 
     +      XL(IX,IY,IZ),XT(IX,IY,IZ),XLT(IX,IY,IZ),XTT(IX,IY,IZ) 
  
      CLOSE(10) 
  
      RETURN
      END

C
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
C    

      DOUBLE PRECISION FUNCTION INTRPOL(XDAT,X,Y,Z)
      
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
C
C  Author: E. Voutier
C  Date: October 2005
C  Purpose: 3-dimensional linear interpolation 
C  
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
C  
C  Determine the selected partial cross section (response function) via a 
C  3-dimensional linear interpolation. 
C  
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
C
C  INPUT variables : 
C
C        XDAT    3-dimensional array of cross section data         
C        X       quadrimomentum transfer
C        Y       energy transfer
C        Z       proton angle in the enter of mass frame
C
C  OUTPUT variable :
C
C        INTRPOL interpolated value of the cross section at (X,Y,Z)
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

      IMPLICIT NONE
      
      INTEGER*4 NQ2,NOM,NTA
        
      INTEGER*4 IX1,IY1,IZ1,IX2,IY2,IZ2

      REAL*8 Q2_STEP,OM_STEP,TA_STEP
      
      REAL*8 XDAT(NQ2,NOM,NTA),X,Y,Z
      
      REAL*8 S_111,S_121,S_112,S_122,S_211,S_221,S_212,S_222
      REAL*8 AX,AY,AZ,AXY,AXZ,AYZ,AXYZ
      
      COMMON / GRID_DEU / Q2_STEP,OM_STEP,TA_STEP
      COMMON / GRID_PAR / NQ2,NOM,NTA
      
      INTRPOL = 0.D+00
      
*---- Lower array index

      IX1 = IDINT( X / Q2_STEP )
      IY1 = IDINT( Y / OM_STEP )
      IZ1 = IDINT( Z / TA_STEP ) + 1
      
*---- Upper array index

      IX2 = IX1 + 1
      IY2 = IY1 + 1
      IZ2 = IZ1 + 1
*-    Correction for grid boundaries
      IF(IX1.EQ.NQ2) IX2 = IX1
      IF(IY1.EQ.NOM) IY2 = IY1
      IF(IZ1.EQ.NTA) IZ2 = IZ1
      
      IF( (IX1.LT.1).OR.(IX2.GT.NQ2).OR.(IY1.LT.1).OR.(IY2.GT.NOM) )THEN
         WRITE(6,*) ' Cross section grid out of range '
         INTRPOL = 0.D0
         RETURN
      ENDIF
      
      AX   = (  X - Q2_STEP * DFLOAT(IX1)   ) / Q2_STEP 
      AY   = (  Y - OM_STEP * DFLOAT(IY1)   ) / OM_STEP 
      AZ   = (  Z - TA_STEP * DFLOAT(IZ1-1) ) / TA_STEP
      AXY  =   AX * AY
      AXZ  =   AX * AZ
      AYZ  =   AY * AZ
      AXYZ =   AX * AY * AZ
      
      S_111 = XDAT(IX1,IY1,IZ1)
      S_121 = XDAT(IX1,IY2,IZ1)
      S_112 = XDAT(IX1,IY1,IZ2)
      S_122 = XDAT(IX1,IY2,IZ2)
      S_211 = XDAT(IX2,IY1,IZ1)
      S_221 = XDAT(IX2,IY2,IZ1)
      S_212 = XDAT(IX2,IY1,IZ2)
      S_222 = XDAT(IX2,IY2,IZ2)
      
      INTRPOL =   S_111 * 
     +           ( 1.D+00 - AX - AY - AZ + AXY + AXZ + AYZ - AXYZ )
     +	      +   S_121 *
     +	        ( AY - AXY - AYZ + AXYZ )
     +	      +   S_112 *
     +	        ( AZ - AXZ - AYZ + AXYZ )
     +	      +   S_122 *
     +	        ( AYZ - AXYZ )
     +	      +	  S_211 *
     + 	        ( AX - AXY - AXZ + AXYZ )
     +	      +	  S_221 *
     +	        ( AXY - AXYZ )
     +	      +	  S_212 *
     + 	        ( AXZ - AXYZ )	     	
     +	      +	  S_222 *
     +	        ( AXYZ )
     
      RETURN
      END

C
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
C
 
      DOUBLE PRECISION FUNCTION INTRLOG(XDAT,X,Y,Z)
       
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C  Author: E. Voutier
C  Date: October 2005
C  Purpose: 3-dimensional interpolation (2 linear + 1 logarythmic)
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C  Determine the selected partial cross section (response function) via a 
C  3-dimensional interpolation assuming a  linear  interpolation in X and Y, 
C  and a logarythmic  interpolation  in  f(Z) corresponding to the recoil
C  momentum. 
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
C
C  INPUT variables : 
C
C        XDAT    3-dimensional array of cross section data         
C        X       quadrimomentum transfer
C        Y       energy transfer
C        Z       proton angle in the center of mass frame
C
C  OUTPUT variable :
C
C        INTRLOG interpolated value of the cross section at (X,Y,Z)
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
 
      IMPLICIT NONE
       
      INTEGER*4 NQ2,NOM,NTA
         
      INTEGER*4 IX1,IY1,IZ1,IX2,IY2,IZ2
 
      REAL*8 Q2_STEP,OM_STEP,TA_STEP
       
      REAL*8 XDAT(NQ2,NOM,NTA),X,Y,Z
      
      REAL*8 INLGCOR
       
      REAL*8 S_111,S_121,S_112,S_122,S_211,S_221,S_212,S_222
      REAL*8 Q2_1,Q2_2,OM_1,OM_2,PT_1,PT_2
      REAL*8 S_11,S_12,S_21,S_22
      REAL*8 AX,AY,AXY
        
      COMMON / GRID_DEU / Q2_STEP,OM_STEP,TA_STEP
      COMMON / GRID_PAR / NQ2,NOM,NTA
       
      INTRLOG = 0.D+00
       
*---- Lower array index
 
      IX1 = IDINT( X / Q2_STEP )
      IY1 = IDINT( Y / OM_STEP )
      IZ1 = IDINT( Z / TA_STEP ) + 1
       
*---- Upper array index
 
      IX2 = IX1 + 1
      IY2 = IY1 + 1
      IZ2 = IZ1 + 1
*-    Correction for grid boundaries
      IF(IX1.EQ.NQ2) IX2 = IX1
      IF(IY1.EQ.NOM) IY2 = IY1
      IF(IZ1.EQ.NTA) IZ2 = IZ1
       
      IF( (IX1.LT.1).OR.(IX2.GT.NQ2).OR.(IY1.LT.1).OR.(IY2.GT.NOM) )THEN
         WRITE(6,*) ' Cross section grid out of range '
         INTRLOG = 0.D0
         RETURN
      ENDIF
       
      Q2_1 = Q2_STEP * DFLOAT(IX1)
      Q2_2 = Q2_STEP + DFLOAT(IX2)
      OM_1 = OM_STEP * DFLOAT(IY1)
      OM_2 = OM_STEP * DFLOAT(IY2)
      PT_1 = TA_STEP * DFLOAT(IZ1-1)
      PT_2 = TA_STEP * DFLOAT(IZ2-1)
      
      S_111 = XDAT(IX1,IY1,IZ1)
      S_121 = XDAT(IX1,IY2,IZ1)
      S_112 = XDAT(IX1,IY1,IZ2)
      S_122 = XDAT(IX1,IY2,IZ2)
      S_211 = XDAT(IX2,IY1,IZ1)
      S_221 = XDAT(IX2,IY2,IZ1)
      S_212 = XDAT(IX2,IY1,IZ2)
      S_222 = XDAT(IX2,IY2,IZ2)
           
      S_11 = INLGCOR(Q2_1,OM_1,PT_1,Z,PT_2,S_111,S_112)
      S_12 = INLGCOR(Q2_1,OM_2,PT_1,Z,PT_2,S_121,S_122)
      S_21 = INLGCOR(Q2_2,OM_1,PT_1,Z,PT_2,S_211,S_212)
      S_22 = INLGCOR(Q2_2,OM_2,PT_1,Z,PT_2,S_221,S_222)
              
      AX   = (  X - Q2_1 ) / Q2_STEP
      AY   = (  Y - OM_1 ) / OM_STEP
      AXY  =   AX * AY
       
      INTRLOG =   S_11 * ( 1.D+00 - AX - AY + AXY )
     +        +   S_12 * ( AY - AXY )
     +        +   S_21 * ( AX - AXY )
     +        +   S_22 * ( AXY )
      
      RETURN
      END
 
C
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
C
 
      DOUBLE PRECISION FUNCTION INLGCOR(X,Y,Z1,Z,Z2,F1,F2)
       
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
C
C  Author: E. Voutier
C  Date: October 2005
C  Purpose: 1-dimensional logarythmic interpolation
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C  Perform a logarythmic interpolation between the points Z1 and Z2 for 
C  non-singular F1 value.  Singular cases are reduced  to  a  linear 
C  interpolation  and correspond to  a nul F1 (phase space effects at the 
C  kinematic boundaries) and a sign change between F1 and F2.   
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C  INPUT variables : 
C
C        XDAT    3-dimensional array of cross section data         
C        X       quadrimomentum transfer
C        Y       energy transfer
C        Z1      lower bound of the proton angle
C        Z       proton angle in the center of mass frame
C        Z2      upper bound of the proton angle
C        F1      cross section value at (X,Y,Z1)
C        F2      cross section value at (X,Y,Z2)
C
C  OUTPUT variable :
C
C        INLGCOR interpolated value of the function at Z
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
 
      IMPLICIT NONE
      
      REAL*8 P_RECOIL
      
      REAL*8 X,Y,Z1,Z,Z2,F1,F2
      
      REAL*8 DZ,F1R,F1S,F2S
      REAL*8 R,R1,R2,DR

      INTEGER LAGET_PS_FAIL,LAGET_GRID_FAIL

      COMMON / LAGET_CNT/ LAGET_PS_FAIL,LAGET_GRID_FAIL
      
      INLGCOR = F1

      DZ = Z2 - Z1 

*---- Grid boundary
      
      IF( DZ.EQ.0.D+00 ) THEN
         LAGET_GRID_FAIL = LAGET_GRID_FAIL + 1
         GO TO 1
      ENDIF
      
*---- Sign extraction

      F1S = 1.D+00
      IF( F1.LT.0.D+00 ) F1S = -1.D+00
      F2S = 1.D+00
      IF( F2.LT.0.D+00 ) F2S = -1.D+00

*---- Recoil neutron momenta
      
      R  = P_RECOIL(X,Y, Z)
      R1 = P_RECOIL(X,Y,Z1)
      R2 = P_RECOIL(X,Y,Z2)
      DR = ( R - R1 ) / ( R2 - R1 )
      IF( DR.EQ.0.D+00 ) GO TO 1

*---- Quasi-threshold identification

      F1R = 1.D+06
      IF( F1.NE.0.D+00 ) F1R = DABS( F2/F1 )

*---- Linear interpolation in angle for quasi-threshold effects
      
      IF( F1R.GT.1.D+05 ) THEN
      
         INLGCOR = F1 + ( ( F2 - F1 ) * ( Z - Z1 ) / DZ ) 

      ELSEIF( DABS(F1S+F2S).LT.1.D+00 ) THEN
      
*---- Linear interpolation in p_r for oscillating sign

         INLGCOR = F1 + ( F2 - F1 ) * DR 
      
      ELSE
	  
*---- Logarythmic interpolation for constant sign
*     PEU:  Check for both signs negative.

         IF (F1S+F2S .LT. -1.5) THEN    ! Both signs negative
            INLGCOR = -( (-F1)**(1.D+00-DR) ) * ( (-F2)**DR )
         ELSE
            INLGCOR = ( F1**(1.D+00-DR) ) * ( F2**DR )
         ENDIF
      
      ENDIF
     
    1 RETURN
      END
       
C
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
C
 
      DOUBLE PRECISION FUNCTION P_RECOIL(G_MAS,G_ENE,P_TGS)
       
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
C
C  Author: E. Voutier
C  Date: October 2005
C  Purpose: calculate the 3-momentum of the recoil neutron
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
C
C  Determine the 3-momentum of the recoiling neutron in the  laboratory frame, 
C  for the current kinematics specified by the virtual photon and the center 
C  of mass angle of the proton.
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
C
C  INPUT variables : 
C
C        G_MAS    quadrimomentum transfer         
C        G_ENE    energy transfer
C        P_TGS    proton angle in the center of mass frame
C
C  OUTPUT variable :
C
C        P_RECOIL corresponding recoil momentum in the lab frame
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
 
      IMPLICIT NONE
      
      REAL*8 G_MAS,G_ENE,P_TGS
      
      REAL*8 DPN,XMD,SPN2,XMD2,XMN2
      
      REAL*8 S_CMS,W_CMS,R_CMS,E_CMS
      REAL*8 G_AMA,G_ABE,R_PAR,R_PER
      
      COMMON / SOME_MAS / DPN,XMD,SPN2,XMD2,XMN2
      
*---- Center of mass frame

      S_CMS = XMD2 - G_MAS + 2.D+00*XMD*G_ENE                 ! Invariant 
                                                                ! squared mass
      W_CMS = DSQRT( S_CMS )                                  ! Invariant mass
      G_AMA = ( G_ENE + XMD ) / W_CMS                         ! Gamma
      G_ABE = DSQRT( G_MAS + G_ENE * G_ENE ) / W_CMS          ! Gamma * Beta

*---- Recoil neutron in the center of mass frame

      R_CMS = ( S_CMS - SPN2 ) * ( S_CMS - DPN * DPN ) / S_CMS
      R_CMS = 0.5D+00 * DSQRT( R_CMS )                        ! Momentum in the
                                                                ! cms frame
      E_CMS = DSQRT( R_CMS * R_CMS + XMN2 )                   ! Energy in the 
                                                                ! cms frame
      
*---- Recoil neutron in the lab frame
      
      R_PAR =   G_ABE * E_CMS - G_AMA * R_CMS * DCOS( P_TGS ) ! Parallel mom.
                                                                ! component
      R_PER = - R_CMS * DSIN( P_TGS )                         ! Perp. mom.
                                                                ! component 
      
      P_RECOIL = DSQRT( R_PAR*R_PAR + R_PER*R_PER )           ! Recoil momentum
                                                                ! in lab frame
      
      RETURN
      END
      
C
C-----------------------------------------------------------------------------
C  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  
C-----------------------------------------------------------------------------
C
