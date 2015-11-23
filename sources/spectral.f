C
C -------------------------------------------------------------------------
C -------------------------------------------------------------------------
C       Subroutine SPECT_SETUP
C
C       AUTHORS:  A. Saha/P.E. Ulmer
C       DATE:     05-SEPT-1990
C
C       PURPOSE:
C               Select spectral function or momentum distribution
C               from a menu.  Set up all arrays or open and read
C               any necessary files.
C
C               All spectral functions are returned in fm^3 and
C               converted to MeV/c^-3 or (MeV/c^-3)-(MeV)^-1 (continuum)
C               within the routine SPECTRAL.
C               Momenta passed to each spectral function routine
C               are in MeV/c.
C -------------------------------------------------------------------------
C
      SUBROUTINE SPECT_SETUP(BOUND)
C
      IMPLICIT NONE
      COMMON /PHYSVAR/ IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
      INTEGER IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
      LOGICAL BOUND
C
    1 WRITE(6,800)
  800 FORMAT(
     # //' Unit Weighting:  ----------------------------------  1'
     # //' Enter Option or <RET> for more options>')
      READ(5,1001,ERR=1) ISPEC_OPT
      IF(ISPEC_OPT .NE. 0) GOTO 1111
C
    5 WRITE(6,900)
  900 FORMAT(
     # //' Deuteron:  Krautschneider ------------------------- 10'
     #  /'            Saha parametric fit to data ------------ 11'
     #  /'            from EPC code -------------------------- 12'
     #  /'            Bernheim data -------------------------- 13'
     #  /'            Van Orden parameterization ------------- 14'
     #  /'            Paris ---------------------------------- 15'
     #  /'            Bonn ----------------------------------- 16'
     #  /'            Argonne V18 ---------------------------- 17'
     #  /'            CD Bonn -------------------------------- 18'
     #  /'            Gross ---------------------------------- 19'
     # //' Enter Option or <RET> for more options>')
      READ(5,1001,ERR=5) ISPEC_OPT
 1001 FORMAT(I3)
      IF(ISPEC_OPT .NE. 0) GOTO 1111
C
   10 WRITE(6,910)
  910 FORMAT(
     # //' 3He:       Meier-Hajduk et al. ------------------- 20'
     #  /'            Argonne + VII ------------------------- 21'
     #  /'            Urbana + VII -------------------------- 22'
     #  /'            Jans/Marchand data -------------------- 23'
     #  /'            Salme --------------------------------- 24'
     # //' Enter Option or <RET> for more options>')
      READ(5,1001,ERR=10) ISPEC_OPT
      IF(ISPEC_OPT .NE. 0) GOTO 1111
C
   20 WRITE(6,920)
  920 FORMAT(
     #  /' 4He:       Ciofi degli Atti fit ------------------ 30'
     #  /'            Argonne + VII ------------------------- 31'
     #  /'            Urbana + VII -------------------------- 32'
     #  /'            v.d. Brand data (q=431 MeV/c) --------- 33'
     #  /'            v.d. Brand data (q=250 MeV/c) --------- 34'
     # //' Enter Option or <RET> for more options>')
      READ(5,1001,ERR=20) ISPEC_OPT
      IF(ISPEC_OPT .NE. 0) GOTO 1111
C
   30 WRITE(6,930)
  930 FORMAT(
     #  /' 16O:       1p1/2 Van Orden ----------------------- 40'
     #  /'            1p3/2 Van Orden ----------------------- 41'
     #  /'            1s1/2 Van Orden ----------------------- 42'
C
C BEGIN MODIFICATION 01
C -------------------------------------------------------------------------
C 010216 KF - add Udias calcs for 16O bound states
C calculations are Jose Udias, private communication
C RDWIA for E0 = 4050 MeV, Tp = 432 MeV, parallel kinematics
C file is /mceep/dat/16o_ju_dw.dat
C -------------------------------------------------------------------------
     #  /'            1p1/2 Udias RDWIA --------------------- 43'
     #  /'            1p3/2 Udias RDWIA --------------------- 44'
     #  /'            1s1/2 Udias RDWIA --------------------- 45'
C -------------------------------------------------------------------------
C END MODIFICATION 01
C
     # //' Enter Option or <RET> for more options>')
      READ(5,1001,ERR=30) ISPEC_OPT
      IF(ISPEC_OPT .NE. 0) GOTO 1111
C
C BEGIN MODIFICATION 02
C -------------------------------------------------------------------------
C 001010 KF - add Lapikas calcs for 208Pb bound states
C calculations are Louk Lapikas, private communication
C CDWIA for E0 = 674 MeV, Tp = 160 MeV, parallel kinematics
C file is /mceep/dat/208pb_ll_dw.dat
C -------------------------------------------------------------------------
   35 WRITE(6,935)
  935 FORMAT(
     #  /' 208Pb:     3s1/2 Lapikas DWIA -------------------- 80'
     #  /'            2d3/2 Lapikas DWIA -------------------- 81'
     #  /'           1h11/2 Lapikas DWIA -------------------- 82'
     #  /'            2d5/2 Lapikas DWIA -------------------- 83'
     #  /'            1g7/2 Lapikas DWIA -------------------- 84'
     #  /'            1g9/2 Lapikas DWIA -------------------- 85'
     # //' Enter Option or <RET> for more options>')
      READ(5,1001,ERR=35) ISPEC_OPT
      IF(ISPEC_OPT .NE. 0) GOTO 1111
C -------------------------------------------------------------------------
C END MODIFICATION 02
C
C BEGIN MODIFICATION 03
C -------------------------------------------------------------------------
C 010226 KF - added 100 to toggles 50 - 62 to free up some integers
C -------------------------------------------------------------------------
   40 WRITE(6,940)
  940 FORMAT(
     #  /' Generic:   Harmonic Oscillator  1s -------------- 150'
     #  /'            Harmonic Oscillator  1p -------------- 151'
     #  /'            Harmonic Oscillator  2s -------------- 152'
     #  /'            Harmonic Oscillator  1d -------------- 153'
     #  /'            Harmonic Oscillator  2p -------------- 154'
     #  /'            Harmonic Oscillator  1f -------------- 155'
     #  /'            Harmonic Oscillator  3s -------------- 156'
     #  /'            Harmonic Oscillator  2d -------------- 157' 
     #  /'            Harmonic Oscillator  1g -------------- 158'
     # //' Enter Option or <RET> for more options>')
      READ(5,1001,ERR=40) ISPEC_OPT
      IF(ISPEC_OPT .NE. 0) GOTO 1111
C
   50 WRITE(6,950)
  950 FORMAT(
     #  /' From file: S(p,E) (continuum) ------------------ 160'
     #  /'            Phi(p)^2 (bound state) -------------- 161'
     #  /'            Phi(p)^2 and f(E) (continuum) ------- 162'
     # //' Enter Option or <RET> to see all options again>')
C -------------------------------------------------------------------------
C END MODIFICATION 03
C
      READ(5,1001,ERR=50) ISPEC_OPT
C
 1111 IF(ISPEC_OPT .EQ. 1) THEN
        IF(.NOT. BOUND) CALL GET_EMDIST         !Unit weighting
C
      ELSEIF(ISPEC_OPT .EQ. 10) THEN
        IF(.NOT. BOUND) CALL GET_EMDIST         !Deuteron-Krautschneider
C
      ELSEIF(ISPEC_OPT .EQ. 11) THEN
        IF(.NOT. BOUND) CALL GET_EMDIST         !Deuteron - Saha
C
      ELSEIF(ISPEC_OPT .EQ. 12) THEN
        IF(.NOT. BOUND) CALL GET_EMDIST         !Deuteron - EPC
C
      ELSEIF(ISPEC_OPT .EQ. 13) THEN
        CALL GET_DEUT_BERN                      !Deuteron - Bernheim
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 14) THEN
        CALL GET_DEUT_VO                        !Deuteron - Van Orden
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 15) THEN
        CALL GET_DEUT_PARIS                     !Deuteron - Paris
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 16) THEN
        CALL GET_DEUT_BONN                      !Deuteron - Bonn
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 17) THEN
        CALL GET_DEUT_AV18                      !Deuteron - AV18
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 18) THEN
        CALL GET_DEUT_CDBONN                    !Deuteron - CD Bonn
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 19) THEN
        CALL GET_DEUT_GROSS                     !Deuteron - Gross
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 20) THEN
        CALL HE3READ(BOUND)                     !3He - Meier Hajduk
C
      ELSEIF(ISPEC_OPT .EQ. 21) THEN
        CALL GET_HE3_ARGONNE                    !3He - Argonne
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 22) THEN
        CALL GET_HE3_URBANA                     !3He - Urbana
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 23) THEN
        CALL GET_HE3_EXP                        !3He - Experimental fit
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 24) THEN
        CALL HE3READ_SALME(BOUND)               !3He - Salme
C
      ELSEIF(ISPEC_OPT .EQ. 30) THEN
        IF(.NOT. BOUND) CALL GET_EMDIST         !4He - Ciofi degli Atti
C
      ELSEIF(ISPEC_OPT .EQ. 31) THEN
        CALL GET_HE4_ARGONNE                    !4He - Argonne
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 32) THEN
        CALL GET_HE4_URBANA                     !4He Urbana
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 33) THEN
        CALL GET_HE4_EXP_Q1                     !4He - data (q=431 MeV/c)
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 34) THEN
        CALL GET_HE4_EXP_Q2                     !4He - data (q=250 MeV/c)
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 40) THEN
        CALL GET_O16P12_VO                      !16O p1/2 - Van Orden
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 41) THEN
        CALL GET_O16P32_VO                      !16O p3/2 - Van Orden
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 42) THEN
        CALL GET_O16S12_VO                      !16O s1/2 - Van Orden
        IF(.NOT. BOUND) CALL GET_EMDIST
C
C BEGIN MODIFICATION 04
C -------------------------------------------------------------------------
C 010216 KF - add Udias calcs for 16O bound states
C RDWIA for E0 = 4050 MeV, Tp = 432 MeV, parallel kinematics
C -------------------------------------------------------------------------
      ELSEIF(ISPEC_OPT .EQ. 43) THEN
        CALL GET_O16P12_UD                      !16O p1/2 - Udias
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 44) THEN
        CALL GET_O16P32_UD                      !16O p3/2 - Udias
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 45) THEN
        CALL GET_O16S12_UD                      !16O s1/2 - Udias
        IF(.NOT. BOUND) CALL GET_EMDIST
C -------------------------------------------------------------------------
C END MODIFICATION 04
C
C BEGIN MODIFICATION 05
C -------------------------------------------------------------------------
C 001010 KF - add Lapikas calcs for 208Pb bound states
C CDWIA for E0 = 674 MeV, Tp = 160 MeV, parallel kinematics
C -------------------------------------------------------------------------
      ELSEIF(ISPEC_OPT .EQ. 80) THEN
        CALL GET_PB2S12_LL                      !208Pb  3s1/2 - Lapikas
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 81) THEN
        CALL GET_PB2D32_LL                      !208Pb  2d3/2 - Lapikas
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 82) THEN
        CALL GET_PB2H12_LL                      !208Pb 1h11/2 - Lapikas
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 83) THEN
        CALL GET_PB2D52_LL                      !208Pb  2d5/2 - Lapikas
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 84) THEN
        CALL GET_PB1G72_LL                      !208Pb  1g7/2 - Lapikas
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 85) THEN
        CALL GET_PB1G92_LL                      !208Pb  1g9/2 - Lapikas
        IF(.NOT. BOUND) CALL GET_EMDIST
C -------------------------------------------------------------------------
C END MODIFICATION 05
C
C BEGIN MODIFICATION 06
C -------------------------------------------------------------------------
C 010226 KF - added 100 to toggles 50 - 62 to free up some integers
C -------------------------------------------------------------------------
      ELSEIF(ISPEC_OPT .EQ. 150 .OR.
     #       ISPEC_OPT .EQ. 151 .OR.
     #       ISPEC_OPT .EQ. 152 .OR.
     #       ISPEC_OPT .EQ. 153 .OR.
     #       ISPEC_OPT .EQ. 154 .OR.
     #       ISPEC_OPT .EQ. 155 .OR.
     #       ISPEC_OPT .EQ. 156 .OR.                
     #       ISPEC_OPT .EQ. 157 .OR.
     #       ISPEC_OPT .EQ. 158)     THEN      
        CALL GET_OSC_PARAM                            !Harm. osc.
        IF(.NOT. BOUND) CALL GET_EMDIST
C
      ELSEIF(ISPEC_OPT .EQ. 160) THEN
        IF(BOUND) WRITE(6,2010)                   !From file: S(p,E)
 2010   FORMAT(/' Override:  CONTINUUM case assumed ')
        BOUND = .FALSE. !Continuum only
        CALL GET_SPECTRAL
C
      ELSEIF(ISPEC_OPT .EQ. 161) THEN
        IF(.NOT. BOUND) WRITE(6,2020)            !From file: phi(p)^2
 2020   FORMAT(/' Override:  BOUND state case assumed ')
        BOUND = .TRUE.  !Bound state only
        CALL GET_PHISQ
C
      ELSEIF(ISPEC_OPT .EQ. 162) THEN
        IF(BOUND) WRITE(6,2010)                   !From file: phi(p)^2*f(E)
        BOUND = .FALSE. !Continuum only
        CALL GET_PHISQ
        CALL GET_EMDIST
C -------------------------------------------------------------------------
C END MODIFICATION 06
C
      ELSE
        GOTO 5          !Nonexistent option selected
C
      ENDIF
C
      RETURN
      END
C
C -------------------------------------------------------------------------
C -------------------------------------------------------------------------
C       Function SPECTRAL
C
C       AUTHORS:  A. Saha/P.E. Ulmer
C       DATE:     05-SEPT-1990
C
C       PURPOSE:
C               Get the value of the momentum distribution or
C               spectral function for a particular event (kinematics)
C               by evaluating a formula or by interpolating on
C               arrays formed in subroutine SPECT_SETUP.
C
C               Momenta passed to each routine are in MeV/c.
C               All momentum distributions are returned in fm^3.
C
C               The UNITS are changed so that the function
C               SPECTRAL is in MeV/c^-3 (bound state) or
C               (MeV/c^-3)-(MeV)^-1 (continuum).
C -------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION SPECTRAL(BOUND,SPEC_FAC,PRMAG,MISS_M)
C
      IMPLICIT NONE
      COMMON /PHYSVAR/ IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
C
      DOUBLE PRECISION EVAL_EMDIST,EVAL_PHISQ,EVAL_SPECTRAL
      DOUBLE PRECISION DEUTDIST,DEUTSAHA,DEUT_EPC,DEUT_BERN,DEUT_VO
      DOUBLE PRECISION DEUT_PARIS
      DOUBLE PRECISION DEUT_BONN,DEUT_AV18,DEUT_CDBONN,DEUT_GROSS
      DOUBLE PRECISION HE3SPEC,HE3_ARGONNE,HE3_URBANA,HE3_EXP
      DOUBLE PRECISION HE3SPEC_SALME
      DOUBLE PRECISION HE4_CIOFI,HE4_ARGONNE,HE4_URBANA
      DOUBLE PRECISION HE4_EXP_Q1,HE4_EXP_Q2
C
C BEGIN MODIFICATION 07
C -------------------------------------------------------------------------
C 010216 KF - add Udias calcs for 16O bound states
C RDWIA for E0 = 4050 MeV, Tp = 432 MeV, parallel kinematics
C -------------------------------------------------------------------------
      DOUBLE PRECISION O16P12_UD,O16P32_UD,O16S12_UD
C -------------------------------------------------------------------------
C END MODIFICATION 07
C
      DOUBLE PRECISION O16P12_VO,O16P32_VO,O16S12_VO
C
C BEGIN MODIFICATION 08
C -------------------------------------------------------------------------
C 001010 KF - add Lapikas calcs for 208Pb bound states
C CDWIA for E0 = 674 MeV, Tp = 160 MeV, parallel kinematics
C -------------------------------------------------------------------------
      DOUBLE PRECISION PB2S12_LL,PB2D32_LL,PB2H12_LL
      DOUBLE PRECISION PB2D52_LL,PB1G72_LL,PB1G92_LL
C -------------------------------------------------------------------------
C END MODIFICATION 08
C
      DOUBLE PRECISION HARM_OSC
      DOUBLE PRECISION PRMAG,MISS_M,F_EM,SPEC_FAC,HBARC,HBARC3
C
      INTEGER IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
      LOGICAL BOUND
C
      PARAMETER (HBARC  = 197.327)    ! Hbar*C    in  MeV-fm
      PARAMETER (HBARC3 = 7.68369D6)  !(Hbar*C)^3 in (MeV-fm)^3
C
C ---------------------------------------------------------------------
C       Unit weighting
C ---------------------------------------------------------------------
C
      IF(ISPEC_OPT .EQ. 1) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       Deuteron - Krautschneider
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 10) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = DEUTDIST(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       Deuteron - Saha
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 11) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = DEUTSAHA(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       Deuteron - EPC
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 12) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = DEUT_EPC(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       Deuteron - Bernheim
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 13) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = DEUT_BERN(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       Deuteron - Van Orden
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 14) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = DEUT_VO(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       Deuteron - Paris
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 15) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = DEUT_PARIS(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       Deuteron - Bonn
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 16) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = DEUT_BONN(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       Deuteron - Argonne V18
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 17) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = DEUT_AV18(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       Deuteron - CD Bonn
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 18) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = DEUT_CDBONN(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       Deuteron - Gross
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 19) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = DEUT_GROSS(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       3He - Meier Hajduk
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 20) THEN
        IF(BOUND) THEN
          SPECTRAL = HE3SPEC(BOUND,MISS_M,PRMAG)/HBARC3
        ELSE
          SPECTRAL = HE3SPEC(BOUND,MISS_M,PRMAG)/HBARC3/HBARC
        ENDIF
C
C ---------------------------------------------------------------------
C       3He - Argonne
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 21) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = HE3_ARGONNE(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       3He - Urbana
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 22) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = HE3_URBANA(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       3He - Experimental fit
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 23) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = HE3_EXP(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       3He - Salme
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 24) THEN
        IF(BOUND) THEN
          SPECTRAL = HE3SPEC_SALME(BOUND,MISS_M,PRMAG)/HBARC3
        ELSE
          SPECTRAL = HE3SPEC_SALME(BOUND,MISS_M,PRMAG)/HBARC3/HBARC
        ENDIF
C
C ---------------------------------------------------------------------
C       4He - Ciofi degli Atti
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 30) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = HE4_CIOFI(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       4He - Argonne
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 31) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = HE4_ARGONNE(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       4He - Urbana
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 32) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = HE4_URBANA(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       4He - data (q=431 MeV/c)
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 33) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = HE4_EXP_Q1(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       4He - data (q=250 MeV/c)
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 34) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = HE4_EXP_Q2(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       16O p1/2 - Van Orden
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 40) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = O16P12_VO(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       16O p3/2 - Van Orden
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 41) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = O16P32_VO(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       16O s1/2 - Van Orden
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 42) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = O16S12_VO(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C 
C BEGIN MODIFICATION 09
C ---------------------------------------------------------------------
C 010216 KF - add Udias calcs for 16O bound states
C RDWIA for E0 = 4050 MeV, Tp = 432 MeV, parallel kinematics
C ---------------------------------------------------------------------
C       16O p1/2 - Udias
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 43) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = O16P12_UD(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       16O p3/2 - Udias
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 44) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = O16P32_UD(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       16O s1/2 - Udias
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 45) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = O16S12_UD(PRMAG)*F_EM/HBARC3
C ---------------------------------------------------------------------
C END MODIFICATION 09
C 
C BEGIN MODIFICATION 10
C ---------------------------------------------------------------------
C 001010 KF - add Lapikas for 208Pb bound states
C CDWIA for E0 = 674 MeV, Tp = 160 MeV, parallel kinematics
C ---------------------------------------------------------------------
C       208Pb 3s1/2 - Lapikas
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 80) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = PB2S12_LL(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       208Pb 2d3/2 - Lapikas
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 81) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = PB2D32_LL(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       208Pb 1h11/2 - Lapikas
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 82) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = PB2H12_LL(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       208Pb 2d5/2 - Lapikas
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 83) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = PB2D52_LL(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       208Pb 1g7/2 - Lapikas
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 84) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = PB1G72_LL(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       208Pb 1g9/2 - Lapikas
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 85) THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL = PB1G92_LL(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C END MODIFICATION 10
C
C BEGIN MODIFICATION 11
C -------------------------------------------------------------------------
C 010226 KF - added 100 to toggles 50 - 62 to free up some integers
C -------------------------------------------------------------------------
C       Harmonic oscillator  
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 150 .OR.
     #       ISPEC_OPT .EQ. 151 .OR.
     #       ISPEC_OPT .EQ. 152 .OR.
     #       ISPEC_OPT .EQ. 153 .OR.
     #       ISPEC_OPT .EQ. 154 .OR.
     #       ISPEC_OPT .EQ. 155 .OR.
     #       ISPEC_OPT .EQ. 156 .OR.
     #       ISPEC_OPT .EQ. 157 .OR.
     #       ISPEC_OPT .EQ. 158)    THEN
        IF(BOUND) THEN
          F_EM = 1.D0
        ELSE
          F_EM = EVAL_EMDIST(MISS_M)  !MeV^-1
        ENDIF
        SPECTRAL= HARM_OSC(PRMAG)*F_EM/HBARC3
C
C ---------------------------------------------------------------------
C       From file:  S(p,E)
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 160) THEN
        SPECTRAL = EVAL_SPECTRAL(PRMAG,MISS_M)/HBARC3
C
C ---------------------------------------------------------------------
C       From file: phi(p)^2
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 161) THEN
        SPECTRAL = EVAL_PHISQ(PRMAG)/HBARC3
C
C ---------------------------------------------------------------------
C       From file: phi(p)^2*f(E)
C ---------------------------------------------------------------------
C
      ELSEIF(ISPEC_OPT .EQ. 162) THEN
        SPECTRAL = EVAL_PHISQ(PRMAG)*EVAL_EMDIST(MISS_M)/HBARC3
C
      ENDIF
C
C ---------------------------------------------------------------------
C END MODIFICATION 11
C
C ---------------------------------------------------------------------
C       Multiply spectral function by spectroscopic factor.
C ---------------------------------------------------------------------
C
      SPECTRAL = SPECTRAL*SPEC_FAC
C
      RETURN
      END
C
C -------------------------------------------------------------------
C -------------------------------------------------------------------
C       DEUTERON:  KRAUTSCHNEIDER MOMENTUM DISTRIBUTION.
C
C         Reference:    F. Krautschneider, Ph.D. thesis,
C                       Bonn University, BONN-IR-76-37 (1976);
C                       L. Hulthen and M. Sagawara,
C                       Handb. d. Phys. Bd. 39, S.1.
C
C       Momenta are in MeV/c.
C       Constants are from Krautschneider's Ph.D. thesis
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1
C -------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DEUTDIST(PR)
      IMPLICIT NONE
C
      DOUBLE PRECISION PI,HBARC,ANORM,PR,PR2,T1,T2
C
      PARAMETER (PI = 3.1415926D0)
      PARAMETER (HBARC = 197.3286D0)  !Hbar*c (MeV-fm)
      ANORM = 0.638D0                 !Normalize to one proton (units: MeV)
      PR2 = PR**2
      T1 = PR2 + 2088.D0
      T2 = PR2 + 67600.D0
      DEUTDIST = HBARC**3*4.D0*ANORM*PI*(1.D0/T1-1.D0/T2)**2  !in fm^3
      RETURN
      END
C
C -------------------------------------------------------------------
C -------------------------------------------------------------------
C       DEUTERON:  A. SAHA FIT TO d(e,e'p) DATA
C
C         Reference:    A. Saha, private communication.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.810
C -------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DEUTSAHA(PR)
      IMPLICIT NONE
C
      DOUBLE PRECISION PM,PR,A,B,FPM
C
      PM = PR*0.001D0                 !convert to GeV/c
C
  101 IF (PM.GT.0.015) GO TO 102
      A = -100.*(ALOG10(0.1209)-ALOG10(0.1133))
      B = -100.*(.005*(-5.+ALOG10(0.1133))-.015*(-5.+ALOG10(0.1209)))
      GO TO 201
  102 IF (PM.GT.0.025) GO TO 103
      A = -100.*(ALOG10(0.1133)-ALOG10(0.0875))
      B = -100.*(.015*(-5.+ALOG10(0.0875))-.025*(-5.+ALOG10(0.1133)))
      GO TO 201
  103 IF (PM.GT.0.035) GO TO 104
      A = -100.*(ALOG10(0.875)-ALOG10(0.604))
      B = -100.*(.025*(-6.+ALOG10(0.604))-.035*(-6.+ALOG10(0.875)))
      GO TO 201
  104 IF (PM.GT.0.045) GO TO 105
      A = -100.*(ALOG10(0.604)-ALOG10(0.375))
      B = -100.*(.035*(-6.+ALOG10(0.375))-.045*(-6.+ALOG10(0.604)))
      GO TO 201
  105 IF (PM.GT.0.055) GO TO 106
      A = -100.*(ALOG10(0.375)-ALOG10(0.188))
      B = -100.*(.045*(-6.+ALOG10(0.188))-.055*(-6.+ALOG10(0.375)))
      GO TO 201
  106 IF (PM.GT.0.065) GO TO 107
      A = -100.*(ALOG10(0.188)-ALOG10(0.129))
      B = -100.*(.055*(-6.+ALOG10(0.129))-.065*(-6.+ALOG10(0.188)))
      GO TO 201
  107 IF (PM.GT.0.175) GO TO 108
      A = -9.091*(ALOG10(0.129)-ALOG10(0.002))
      B = -9.091*(.065*(-6.+ALOG10(0.002))-.175*(-6.+ALOG10(0.129)))
      GO TO 201
  108 A = -6.25*(ALOG10(0.2)-ALOG10(0.0048))
      B = -6.25*(.175*(-8.+ALOG10(0.0048))-.335*(-8.+ALOG10(0.2)))
C
  201 FPM = (A*PM)+B
      DEUTSAHA = ((197.3286)**3.)*(10.0**FPM)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       DEUTERON:  FROM EPC PROGRAM OF LIGHTBODY AND O'CONNELL
C
C         Reference:    J.W. Lightbody, Jr. and J.S. O'Connell,
C                       Computers in Physics, May/June 1988, p. 57.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.936
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DEUT_EPC(PR)
      IMPLICIT NONE
C
      DOUBLE PRECISION PP,PR,SGS
C
      PP=PR/197.3286                  !Convert to fm^-1
      SGS=3.694-7.428*PP-2.257*PP**2
      SGS=SGS+3.618*PP**3-1.377*PP**4+.221*PP**5-.103*PP**6
C        IF(SGS.LT.-293.)GO TO 1
      SGS=EXP(SGS)
      DEUT_EPC = SGS/.18825/4./3.141592654
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       DEUTERON:  SACLAY EXPT. BERNHEIM ET AL. (READS THE DATA FILE)
C
C         Reference:    M. Bernheim, et al., Nucl. Phys. A365, 349 (1981).
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.776
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_DEUT_BERN
      IMPLICIT NONE
      COMMON /BERNHEIM/ P_BERN,S_BERN
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_BERN(34),S_BERN(34),DS
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/d_saclay.exp'
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(/)
      DO I=1,34
        READ (3,*) P_BERN(I),S_BERN(I),DS
        S_BERN(I) = S_BERN(I)*((197.3286)**3)         !convert to fm^3
      ENDDO
      CLOSE(UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       DEUTERON:  SACLAY EXPT. BERNHEIM ET AL. (INTERPOLATES ON DATA)
C
C         Reference:    M. Bernheim, et al., Nucl. Phys. A365, 349 (1981).
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.776
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DEUT_BERN(PR)
      IMPLICIT NONE
      COMMON /BERNHEIM/ P_BERN,S_BERN
C
      DOUBLE PRECISION P_BERN(34),S_BERN(34),RINTERPEXP,PR
C
      DEUT_BERN = RINTERPEXP(P_BERN,S_BERN,PR,34)       !in fm^3
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       DEUTERON:  VAN ORDEN DISTRIBUTION (READS THE DISTRIBUTION FILE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_DEUT_VO
      IMPLICIT NONE
      COMMON /DEUTVO/ P_VO,S_VO
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_VO(101),S_VO(101)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/d_vo.dat'
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(/)
      DO I = 1,101
        READ (3,*) P_VO(I),S_VO(I)
        P_VO(I) = P_VO(I) * 1000.D0          !convert to MeV/c
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       DEUTERON:  VAN ORDEN DISTRIBUTION (INTERPOLATES TO GET VALUE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DEUT_VO(PR)
      IMPLICIT NONE
      COMMON /DEUTVO/ P_VO,S_VO
C
      DOUBLE PRECISION P_VO(101),S_VO(101),PR,RINTERPEXP
C
      DEUT_VO = RINTERPEXP(P_VO,S_VO,PR,101)
C
      RETURN
      END
c
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c       DEUTERON:  PARIS WAVEFUNCTION (GET THE COEFFICIENTS)
c
c       Reference:  M. Lacombe et al., Phys. Lett. 101B, 139 (1981).
c
c       Normalization:  4pi x Integral(dp p^2 momdist) = 1
c
c       Coded by:  P.E. Ulmer on 1/14/99
c------------------------------------------------------------------------------
c
      SUBROUTINE GET_DEUT_PARIS
c
      implicit none
c
      COMMON /DEUTPARIS/ C,D,M
c
      double precision c(13),d(13),m(13),alph,m0
      double precision sum_c,sum_d1,sum_d2,sum_d3
      integer nc,n1(3),n2(3),n3(3),j
c
c *** DATA FOR PARIS ***
      DATA c(1),c(2),c(3),c(4),c(5),c(6),
     &     c(7),c(8),c(9),c(10),c(11),c(12)
     &  / 0.88688076D+00, -0.34717093D+00, -0.30502380D+01,
     &    0.56207766D+02, -0.74957334D+03,  0.53365279D+04,
     &   -0.22706863D+05,  0.60434469D+05, -0.10292058D+06,
     &    0.11223357D+06, -0.75925226D+05,  0.29059715D+05/
      DATA d(1),d(2),d(3),d(4),d(5),d(6),
     &     d(7),d(8),d(9),d(10)
     &  / 0.23135193D-01, -0.85604572D+00,  0.56068193D+01,
     &   -0.69462922D+02,  0.41631118D+03, -0.12546621D+04,
     &    0.12387830D+04,  0.33739172D+04, -0.13041151D+05,
     &    0.19512524D+05/
      DATA nc /13/
      DATA alph,m0 /0.23162461D+00 , 1.0D+00/
c
      do j=1,nc
         m(j) = alph + m0*dfloat(j-1)
      enddo
c
      sum_c = 0.d0
      do j=1,nc-1
         sum_c = sum_c + c(j)
      enddo
      c(nc) = -sum_c
c
      sum_d1 = 0.d0
      sum_d2 = 0.d0
      sum_d3 = 0.d0
c
      do j=1,nc-3
         sum_d1 = sum_d1 + d(j)/m(j)**2
         sum_d2 = sum_d2 + d(j)
         sum_d3 = sum_d3 + d(j)*m(j)**2
      enddo
c
      n1(1) = nc-2
      n2(1) = nc-1
      n3(1) = nc
c
      n1(2) = nc-1
      n2(2) = nc
      n3(2) = nc-2
c
      n1(3) = nc
      n2(3) = nc-2
      n3(3) = nc-1
c
      do j=1,3
         d(n1(j)) = (m(n1(j))**2/
     &        ((m(n3(j))**2-m(n1(j))**2)*(m(n2(j))**2-m(n1(j))**2)))
     &        * (-m(n2(j))**2*m(n3(j))**2*sum_d1 
     &        + (m(n2(j))**2+m(n3(j))**2)*sum_d2
     &        - sum_d3)
      enddo
c
      return
      end
c
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c       DEUTERON:  PARIS WAVEFUNCTION (Compute Momdist)
c
c       Reference:  M. Lacombe et al., Phys. Lett. 101B, 139 (1981).
c
c       Normalization:  4pi x Integral(dp p^2 momdist) = 1
c
c       Coded by:  P.E. Ulmer on 1/14/99
c------------------------------------------------------------------------------
c
      DOUBLE PRECISION FUNCTION DEUT_PARIS(PR)
      IMPLICIT NONE
      COMMON /DEUTPARIS/ C,D,M
C
      double precision c(13),d(13),m(13),sum_uc,sum_wc,u,w
      double precision pr,p_invfm,pi,sq2opi,hbarc
      integer nc,j
c
      parameter (pi=3.141592654d0)
      parameter (sq2opi=0.79788456d0) !sqrt(2/pi)
      parameter (hbarc = 197.3286d0)  !hbar*c (MeV-fm)
c
      data nc /13/
c
      p_invfm = pr/hbarc    ! convert to fm^-1
c
      sum_uc = 0.d0
      sum_wc = 0.d0
      do j=1,nc
         sum_uc = sum_uc + c(j)/(p_invfm**2+m(j)**2)
         sum_wc = sum_wc + d(j)/(p_invfm**2+m(j)**2)
      enddo
      u = sq2opi*sum_uc
      w = sq2opi*sum_wc
c
      DEUT_PARIS = (u**2 + w**2)/(4.d0*pi)       ! in fm^3
c
      return
      end
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       The following deuteron momentum distributions were kindly
C       given to me by Sabine Jeschonnek - 3/15/00.
C       I have slightly modified them for compatibility with MCEEP.
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
c
c------------------------------------------------------------
      subroutine get_deut_bonn
c------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c
      common/deutbonn/c,d,am
c
      double precision c(14),d(14),am(14)
c
      parameter(nmax=11)
      parameter(alpha=2.31609d-1,am0=9.d-1)
c
      do j = 1,nmax
         am(j) = alpha+(j-1)*am0
      enddo
c 
c coefficients, full model
c
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
c 
      c(11) = 0.d0
c
      mx1 = nmax-1
c
      do k=1,mx1
         c(11) = c(11)-c(k)
      enddo
c 
c coefficients, full model
c
      d(1) =  0.24133026d-1
      d(2) = -0.64430531d0
      d(3) =  0.51093352d0
      d(4) = -0.54419065d1
      d(5) =  0.15872034d2
      d(6) = -0.14742981d2
      d(7) =  0.44956539d1
      d(8) = -0.71152863d-1
c 
      d(9) = 0.d0
      d(10) = 0.d0
      d(11) = 0.d0
c
      n3 = nmax-3
      sp2 = 0.d0
      sm = 0.d0
      sm2 = 0.d0
c
      do lp = 1,n3
         sp2 = sp2+d(lp)/am(lp)**2
         sm = sm+d(lp)
         sm2 = sm2+d(lp)*am(lp)**2
      enddo
c 
      a = am(11)**2
      b = am(10)**2
      cc = am(9)**2
      d(9) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c 
      a = am(9)**2
      b = am(11)**2
      cc = am(10)**2
      d(10) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c 
      a = am(10)**2
      b = am(9)**2
      cc = am(11)**2
      d(11) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c
      return
      end
c
c------------------------------------------------------------
      DOUBLE PRECISION FUNCTION deut_bonn(pr)
c------------------------------------------------------------
c
      implicit none
c
      common/deutbonn/c,d,am
c
      double precision c(14),d(14),am(14)
      double precision pr,p,um,wm,pi,sq2opi,hbarc
      integer nn,nmax
c
      parameter (nmax=11)
      parameter (pi=3.141592654d0)
      parameter (sq2opi=0.79788456d0) !sqrt(2/pi)
      parameter (hbarc = 197.3286d0)  !hbar*c (MeV-fm)
c
      p = pr/hbarc    ! convert to fm^-1
c
      um = 0.d0
      wm = 0.d0
c
      do nn = 1,nmax
         um = um+c(nn)/(p*p+am(nn)*am(nn))
         wm = wm+d(nn)/(p*p+am(nn)*am(nn))
      enddo
c
      um = sq2opi*um
      wm = sq2opi*wm
c
      deut_bonn = (um**2 + wm**2)/(4.d0*pi)       ! in fm^3
c
      return
      end
c
c------------------------------------------------------------
      subroutine get_deut_av18
c------------------------------------------------------------
c
c Rocco's parametrization of the Argonne V18 potential
c
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c
      common/deutav18/c,d,am
c
      double precision c(14),d(14),am(14)
c
      parameter(nmax=12)
c 
c mass parameters
c
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
c 
c s wave coefficients
c
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
c 
      c(12) = 0.d0
c
      mx1 = nmax-1
c
      do k=1,mx1
         c(12) = c(12)-c(k)
      enddo
c 
c d wave coefficients
c
      d(1) =  0.280995496d+00
      d(2) =  0.334117629d-01
      d(3) = -0.727192237d+00
      d(4) = -0.302809607d+01
      d(5) = -0.903824982d+01
      d(6) =  0.496045967d+01
      d(7) = -0.271985613d+02
      d(8) =  0.125334598d+03
      d(9) = -0.346742235d+03
c
      d(10) = 0.d0
      d(11) = 0.d0
      d(12) = 0.d0
c
      n3 = nmax-3
      sp2 = 0.d0
      sm = 0.d0
      sm2 = 0.d0
c
      do lp = 1,n3
         sp2 = sp2+d(lp)/am(lp)**2
         sm = sm+d(lp)
         sm2 = sm2+d(lp)*am(lp)**2
      enddo
c 
      a = am(12)**2
      b = am(11)**2
      cc = am(10)**2
      d(10) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c 
      a = am(10)**2
      b = am(12)**2
      cc = am(11)**2
      d(11) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c 
      a = am(11)**2
      b = am(10)**2
      cc = am(12)**2
      d(12) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c 
      return
      end
c
c------------------------------------------------------------
      DOUBLE PRECISION FUNCTION deut_av18(pr)
c------------------------------------------------------------
c
      implicit none
c
      common/deutav18/c,d,am
c
      double precision c(14),d(14),am(14)
      double precision pr,p,um,wm,pi,sq2opi,hbarc
      integer nn,nmax
c
      parameter(nmax=12)
      parameter (pi=3.141592654d0)
      parameter (sq2opi=0.79788456d0) !sqrt(2/pi)
      parameter (hbarc = 197.3286d0)  !hbar*c (MeV-fm)
c
      p = pr/hbarc    ! convert to fm^-1
c
      um = 0.d0
      wm = 0.d0
c
      do nn = 1,nmax
         um = um+c(nn)/(p*p+am(nn)*am(nn))
         wm = wm+d(nn)/(p*p+am(nn)*am(nn))
      enddo
c
      um = sq2opi*um/(4.d0*pi)
      wm = sq2opi*wm/(4.d0*pi)
c
      deut_av18 = (um**2 + wm**2)/(4.d0*pi)       ! in fm^3
c
      return
      end
c------------------------------------------------------------
      subroutine get_deut_cdbonn
c------------------------------------------------------------
c
c Rocco's parametrization of the CD Bonn potential (Machleidt)
c
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c
      common/deutcdbonn/c,d,am
c
      double precision c(14),d(14),am(14)
c
      parameter(nmax=12)
c
c mass parameters
c
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
c 
c s wave coefficients
c
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
c 
      c(12) = 0.d0
c
      do k=1,nmax-1
         c(12) = c(12)-c(k)
      enddo
c 
c d wave coefficients
c
      d(1) =  -0.286623664d+00
      d(2) =   0.195863407d+00
      d(3) =  -0.196866168d+01
      d(4) =   0.208599899d+02
      d(5) =  -0.388561750d+02
      d(6) =   0.426753739d+02
      d(7) =  -0.457983138d+02
      d(8) =  -0.535081902d+00
      d(9) =   0.159713978d+03
c
      d(10) = 0.d0
      d(11) = 0.d0
      d(12) = 0.d0
c
      n3 = nmax-3
      sp2 = 0.d0
      sm = 0.d0
      sm2 = 0.d0
c
      do lp = 1,n3
         sp2 = sp2+d(lp)/am(lp)**2
         sm = sm+d(lp)
         sm2 = sm2+d(lp)*am(lp)**2
      enddo
c 
      a = am(12)**2
      b = am(11)**2
      cc = am(10)**2
      d(10) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c 
      a = am(10)**2
      b = am(12)**2
      cc = am(11)**2
      d(11) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c 
      a = am(11)**2
      b = am(10)**2
      cc = am(12)**2
      d(12) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c 
      return
      end
c
c------------------------------------------------------------
      DOUBLE PRECISION FUNCTION deut_cdbonn(pr)
c------------------------------------------------------------
c
      implicit none
c
      common/deutcdbonn/c,d,am
c
      double precision c(14),d(14),am(14)
      double precision pr,p,um,wm,pi,sq2opi,hbarc
      integer nn,nmax
c
      parameter(nmax=12)
      parameter (pi=3.141592654d0)
      parameter (sq2opi=0.79788456d0) !sqrt(2/pi)
      parameter (hbarc = 197.3286d0)  !hbar*c (MeV-fm)
c
      p = pr/hbarc    ! convert to fm^-1
c
      um = 0.d0
      wm = 0.d0
c
      do nn = 1,nmax
         um = um+c(nn)/(p*p+am(nn)*am(nn))
         wm = wm+d(nn)/(p*p+am(nn)*am(nn))
      enddo
c
      um = sq2opi*um/(4.d0*pi)
      wm = sq2opi*wm/(4.d0*pi)
c
c take care of phase:
c
      um = -um
      wm = -wm
c
      deut_cdbonn = (um**2 + wm**2)/(4.d0*pi)       ! in fm^3
c
      return
      end
c
c------------------------------------------------------------
      subroutine get_deut_gross
c------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c
      common/deutgross/c,d,am
c
      double precision c(14),d(14),am(14)
c
      parameter(nmax=12)
c
c mass parameters
c
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
c 
c s wave coefficients
c
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
c
      c(12) = 0.d0
c
      mx1 = nmax-1
c
      do k=1,mx1
         c(12) = c(12)-c(k)
      enddo
c
c      c(12)= -.937069211E+02
c
c d wave coefficients
c
      d(1) = -.304115815E+00
      d(2) = 0.109108134E+01
      d(3) = -.646883058E+01
      d(4) = 0.300888843E+02
      d(5) = -.488288588E+02
      d(6) = 0.484448585E+02
      d(7) = -.651042737E+02
      d(8) = 0.759630709E+02
      d(9) = -.112034936E+03
c
      d(10) = 0.d0
      d(11) = 0.d0
      d(12) = 0.d0
c
      n3 = nmax-3
      sp2 = 0.d0
      sm = 0.d0
      sm2 = 0.d0
c
      do lp = 1,n3
         sp2 = sp2+d(lp)/am(lp)**2
         sm = sm+d(lp)
         sm2 = sm2+d(lp)*am(lp)**2
      enddo
c 
      a = am(12)**2
      b = am(11)**2
      cc = am(10)**2
      d(10) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c 
      a = am(10)**2
      b = am(12)**2
      cc = am(11)**2
      d(11) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c 
      a = am(11)**2
      b = am(10)**2
      cc = am(12)**2
      d(12) = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2)
c
c      d(10) = 0.202201330E+03
c      d(11) = -.215603795E+03
c      d(12) = 0.931991144E+02
c
      return
      end
c
c------------------------------------------------------------
      DOUBLE PRECISION FUNCTION deut_gross(pr)
c------------------------------------------------------------
c
      implicit none
c
      common/deutgross/c,d,am
c
      double precision c(14),d(14),am(14)
      double precision pr,p,um,wm,pi,sq2opi,hbarc,anorm
      integer nn,nmax
c
      parameter(nmax=12)
      parameter (pi=3.141592654d0)
      parameter (sq2opi=0.79788456d0) !sqrt(2/pi)
      parameter (hbarc = 197.3286d0)  !hbar*c (MeV-fm)
c
      p = pr/hbarc    ! convert to fm^-1
c
      um = 0.d0
      wm = 0.d0
c
      do nn = 1,nmax
         um = um+c(nn)/(p*p+am(nn)*am(nn))
         wm = wm+d(nn)/(p*p+am(nn)*am(nn))
      enddo
c
      anorm = 1.0047863577500d0
c
      um = sq2opi*um/(4.d0*pi*sqrt(anorm))
      wm = sq2opi*wm/(4.d0*pi*sqrt(anorm))
c
      deut_gross = (um**2 + wm**2)/(4.d0*pi)       ! in fm^3
c
      return
      end
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       3HE:  MEIER-HAJDUK SPECTRAL FCN (READS THE FILE)
C
C         Reference:    H. Meier-Hajduk et al., Nucl. Phys.
C                       A395, 332 (1983) and private communication.
C
C         Note the units are fm^3 for the bound state (2-body breakup)
C         and fm^4 for the continuum (3-body breakup).  The 2-body
C         breakup channel accounts for roughly 92% of the T=0
C         strength.
C
C       Normalization:  4pi x Integral(dp   p^2 momdist) = 1.62 (2-body)
C                       4pi x Integral(dpdE p^2 momdist) = 1.06 (3-body)
C------------------------------------------------------------------------------
C
      SUBROUTINE HE3READ(BOUND)
      IMPLICIT NONE
      COMMON /HE3/ E,P,S0,S1
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION E(15),P(25),S0(15,25),S1(15,25)
      INTEGER I,J,JI,JF,JMIN,JMAX,K,N_DAT_DIR
      LOGICAL BOUND
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/he3sf.dat'
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ(3,*) (E(I),I=1,5)
      READ(3,*) (E(I),I=6,10)
      READ(3,*) (E(I),I=11,15)
C
      DO I=1,5
        JI=(I-1)*5+1
        JF=(I-1)*5+5
        READ(3,*)(P(K),K=JI,JF)
      ENDDO
C
      DO J=1,15
        DO I=1,5
          JI=(I-1)*5+1
          IF (I.EQ.5) JF=(I-1)*5+4
          IF (I.NE.5) JF=(I-1)*5+5
          READ(3,*)(S0(J,K),K=JI,JF)
        ENDDO
      ENDDO
C
      DO J=1,15
        DO I=1,5
          JI=(I-1)*5+1
          IF (I.EQ.5) JF=(I-1)*5+4
          IF (I.NE.5) JF=(I-1)*5+5
          READ(3,*)(S1(J,K),K=JI,JF)
        ENDDO
      ENDDO
C
C       Zero out the appropriate piece of each spectral function
C       for bound state or continuum.
C
      IF(BOUND) THEN
        JMIN = 2
        JMAX = 15
      ELSE
        JMIN = 1
        JMAX = 1
      ENDIF
C
      DO J=JMIN,JMAX
        DO K=1,24
          S0(J,K) = 0.
          S1(J,K) = 0.
        ENDDO
      ENDDO
C
      CLOSE(UNIT=3)
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       3HE:  MEIER-HAJDUK SPECTRAL FCN (INTERPOLATES TO GET VALUE)
C
C         Reference:    H. Meier-Hajduk et al., Nucl. Phys.
C                       A395, 332 (1983) and private communication.
C
C         Note the units are fm^3 for the bound state (2-body breakup)
C         and fm^4 for the continuum (3-body breakup).  The 2-body
C         breakup channel accounts for roughly 92% of the T=0
C         strength.
C
C       Normalization:  4pi x Integral(dp   p^2 momdist) = 1.62 (2-body)
C                       4pi x Integral(dpdE p^2 momdist) = 1.06 (3-body)
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HE3SPEC(BOUND,EM,PM)
      IMPLICIT NONE
      COMMON /HE3/ E,P,S0,S1
C
      DOUBLE PRECISION E(15),P(25),S0(15,25),S1(15,25)
      DOUBLE PRECISION EM,PM,E_M,P_M,CON,S_0,S_1,S0E1,S0E2,S1E1,S1E2
      INTEGER I,J
      LOGICAL BOUND
C
C------------------------------------------------------------------------------
C       Convert quantities.
C------------------------------------------------------------------------------
C
      PARAMETER (CON = 1.0)
      E_M = EM/197.33                 !MeV   to fm^-1
C      P_M = PM/197.33                 !MeV/c to fm^-1
C------------------------------------------------------------------------------
C  Take absolute value to be safe
C------------------------------------------------------------------------------
C
      P_M = DABS(PM) / 197.33                 !MeV/c to fm^-1
C
C------------------------------------------------------------------------------
C  Determine the missing energy bin, id est the value of index I: 
C------------------------------------------------------------------------------
C
      IF (.NOT.BOUND) THEN
        DO I = 1, 14
          IF (E_M.LE.E(I)) GO TO 301
        ENDDO
        HE3SPEC = 0.   ! NOT TABULATED UP TO THIS ENERGY
        RETURN
      ENDIF
 301  CONTINUE
C
C------------------------------------------------------------------------------
C Determine the missing momentum bin, id est the value of index J:
C------------------------------------------------------------------------------
C
      DO J = 1, 24
        IF (P_M.LE.P(J)) GO TO 302
      ENDDO
      HE3SPEC = 0.   ! NOT TABULATED UP TO THIS MOMENTUM
      RETURN
 302  CONTINUE
C
C------------------------------------------------------------------------------
C Treat BOUND and CONTINUUM state separately:
C------------------------------------------------------------------------------
C
C BOUND STATE:
C
      IF (BOUND) THEN    ! BOUND = NOT CONTINUUM
        IF (E_M.GT.E(1)) THEN 
          WRITE(6,*)' PROBLEM BINDING ENERGY ERROR'
          WRITE(6,*)'(file spectral.f, func he3spec)'
          STOP
      ENDIF
        IF (J.EQ.1) THEN
          S_0 = S0(1,1) / 197.33
        ELSE
          S_0 = (((P_M-P(J-1))/(P(J)-P(J-1)))*(S0(1,J)-S0(1,J-1))
     #          + S0(1,J-1))/197.33
        ENDIF	  
        S_1 = 0. ! ALWAYS 0. FOR BOUND STATES
C
C CONTINUUM STATE:
C
      ELSE               ! NOT BOUND = CONTINUUM
      
        IF (I.LE.2) THEN    ! The continuum spectral fcn
          HE3SPEC = 0.      ! is zero below E_M=E(2)
          RETURN            ! so don't interpolate
C                             unless E_M .gt. E(2)	  
        ELSE     ! I.gt.2
          IF (J.GT.1) THEN  ! interpolation in EM and PM
            S0E1 = ((E_M-E(I-1))/(E(I)-E(I-1)))*(S0(I,J-1)-S0(I-1,J-1))
     #           + S0(I-1,J-1)
            S0E2 = ((E_M-E(I-1))/(E(I)-E(I-1)))*(S0(I,J)-S0(I-1,J))
     #           + S0(I-1,J)
            S1E1 = ((E_M-E(I-1))/(E(I)-E(I-1)))*(S1(I,J-1)-S1(I-1,J-1))
     #           + S1(I-1,J-1)
            S1E2 = ((E_M-E(I-1))/(E(I)-E(I-1)))*(S1(I,J)-S1(I-1,J))
     #           + S1(I-1,J)
            S_0 = ((P_M-P(J-1))/(P(J)-P(J-1)))*(S0E2-S0E1) + S0E1
            S_1 = ((P_M-P(J-1))/(P(J)-P(J-1)))*(S1E2-S1E1) + S1E1
          ELSE   ! J = 1, interpolation only in EM, not in PM
C
C------------------------------------------------------------------------------
C The value of the spectral function is taken equal to its first tabulated
C value in he3sf.dat (at PM = P(J=1)= 5.05 MeV/c)
C and not extrapolated towards PM=0 using both P(J=1) and P(J=2) because its
C derivative versus PM becomes much smaller with a decreasing PM.
C We keep here the interpolation versus EM.	  
C------------------------------------------------------------------------------
C
            S_0 = ((E_M-E(I-1))/(E(I)-E(I-1)))*(S0(I,1)-S0(I-1,1))
     #           + S0(I-1,1)
            S_1 = ((E_M-E(I-1))/(E(I)-E(I-1)))*(S1(I,1)-S1(I-1,1))
     #           + S1(I-1,1)
          ENDIF
        ENDIF

      ENDIF                ! BOUND OR CONTINUUM
C
C      HE3SPEC = 3.*(S_0+S_1/3.)*CON           !in fm^3 or fm^4
      HE3SPEC = 3.*S_0+S_1          !in fm^3 or fm^4
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       3HE:  ARGONNE + MODEL VII (SET UP THE ARRAYS)
C
C         These parameters are from a variational calculation using
C         the Argonne potential for the 2-body breakup of 3He.
C
C         Reference:    R. Schiavilla, V.R. Pandharipande and
C                       R.B. Wiringa, Nucl. Phys. A449, 219 (1986)
C                       and private communication.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C                       where the integral extends to 500 MeV/c
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_HE3_ARGONNE
      IMPLICIT NONE
      COMMON /HE3ARG/ P,Y
C
      DOUBLE PRECISION P(16),A0(16),A2(16),Y(16)
      INTEGER I
C
      DATA A0 / 93.08,82.54,59.55,38.3,24.0,14.9,9.17,5.61,3.33,
     #          1.96,1.12,0.58,0.26,0.06,-0.04,-0.09/
      DATA A2 / 0.0,-0.375,-1.08,-1.53,-1.63,-1.54,-1.37,-1.17,
     #          -0.98,-0.8,-0.65,-0.53,-0.44,-0.36,-0.3,-0.25/
C
      DO I = 1,16
        P(I) = (I-1)*0.16*0.19733
        Y(I) = 1.5*((A0(I)**2)+(A2(I)**2))/((2.*3.141592654)**4)/1.38
      ENDDO
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       3HE:  ARGONNE + MODEL VII (INTERPOLATES TO GET VALUE)
C
C         These parameters are from a variational calculation using
C         the Argonne potential for the 2-body breakup of 3He.
C
C         Reference:    R. Schiavilla, V.R. Pandharipande and
C                       R.B. Wiringa, Nucl. Phys. A449, 219 (1986)
C                       and private communication.
C
C         Note:  this parameterization becomes unreasonable beyond
C                Pr=500 MeV/c; thus, the distribution is set to 0
C                beyond this.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C                       where the integral extends to 500 MeV/c.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HE3_ARGONNE(PR)
      IMPLICIT NONE
C
      INCLUDE 'summary.cmn'
C
      COMMON /HE3ARG/ P,Y
C
      DOUBLE PRECISION P(16),Y(16),PR,PM,RINTERPEXP
C
      IF(PR .GT. 500.D0) THEN
CXXX    TYPE *,' Pr > 500 MeV/c --> setting distribution to 0 '
        I_SPECT_OOR = I_SPECT_OOR + 1
        HE3_ARGONNE = 0.D0
      ELSE
        PM = PR*0.001D0                 !convert to GeV/c
        HE3_ARGONNE = RINTERPEXP(P,Y,PM,16)
      ENDIF
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       3HE:  URBANA + MODEL VII (SET UP THE ARRAYS)
C
C         These parameters are from a variational calculation using
C         the Urbana potential for the 2-body breakup of 3He.
C
C         Reference:    R. Schiavilla, V.R. Pandharipande and
C                       R.B. Wiringa, Nucl. Phys. A449, 219 (1986)
C                       and private communication.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C                       where the integral extends to 500 MeV/c.
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_HE3_URBANA
      IMPLICIT NONE
      COMMON /HE3URB/ P,Y
C
      DOUBLE PRECISION P(16),A0(16),A2(16),Y(16)
      INTEGER I
C
      DATA A0 / 97.8,85.99,60.7,38.2,23.5,14.4,8.75,5.37,3.25,
     #          1.96,1.17,0.64,0.31,0.1,-0.02,-0.07/
      DATA A2 / 0.0,-0.461,-1.31,-1.79,-1.86,-1.74,-1.53,-1.3,
     #          -1.07,-0.87,-0.7,-0.56,-0.45,-0.37,-0.29,-0.23/
C
      DO I = 1,16
        P(I) = (I-1)*0.16*0.19733
        Y(I) = 1.5*((A0(I)**2)+(A2(I)**2))/((2.*3.141592654)**4)/1.38
      ENDDO
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       3HE:  URBANA + MODEL VII (INTERPOLATES TO GET VALUE)
C
C         These parameters are from a variational calculation using
C         the Urbana potential for the 2-body breakup of 3He.
C
C         Reference:    R. Schiavilla, V.R. Pandharipande and
C                       R.B. Wiringa, Nucl. Phys. A449, 219 (1986)
C                       and private communication.
C
C         Note:  this parameterization becomes unreasonable beyond
C                Pr=500 MeV/c; thus, the distribution is set to 0
C                beyond this.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C                       where the integral extends to 500 MeV/c.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HE3_URBANA(PR)
      IMPLICIT NONE
C
      INCLUDE 'summary.cmn'
C
      COMMON /HE3URB/ P,Y
C
      DOUBLE PRECISION P(16),Y(16),PR,PM,RINTERPEXP
C
      IF(PR .GT. 500.D0) THEN
CXXX    TYPE *,' Pr > 500 MeV/c --> setting distribution to 0 '
        I_SPECT_OOR = I_SPECT_OOR + 1
        HE3_URBANA = 0.D0
      ELSE
        PM = PR*0.001D0                 !convert to GeV/c
        HE3_URBANA = RINTERPEXP(P,Y,PM,16)
      ENDIF
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       3HE:   SACLAY EXPT. JANS + MARCHAND (READS THE DATA FILE)
C
C         References:  Pr=  5-305 MeV/c:  E. Jans et al., Phys. Rev.
C                                         Lett. 49, 974 (1982).
C                      Pr=318-360 MeV/c:  C. Marchand et al., Phys. Rev.
C                                         Lett. 60, 1703 (1988).
C
C         Note these data are for the 2-body break-up channel only.
C         Also, the data for 155-195 MeV/c are averages of the
C         Kin.I and Kin.II data sets.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.07
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_HE3_EXP
      IMPLICIT NONE
      COMMON /HE3EXP/ P_HE3,S_HE3
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_HE3(36),S_HE3(36),DS
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/3he_saclay.exp'
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(///)
      DO I=1,36
        READ (3,*) P_HE3(I),S_HE3(I),DS
        S_HE3(I) = S_HE3(I)*((0.19733)**3)   !convert to fm^3
      ENDDO
      CLOSE(UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       3HE:  SACLAY EXPT. JANS + MARCHAND (INTERPOLATES TO GET VALUE)
C
C         References:  Pr=  5-305 MeV/c:  E. Jans et al., Phys. Rev.
C                                         Lett. 49, 974 (1982).
C                      Pr=318-360 MeV/c:  C. Marchand et al., Phys. Rev.
C                                         Lett. 60, 1703 (1988).
C
C         Note these data are for the 2-body break-up channel only.
C         Also, the data for 155-195 MeV/c are averages of the
C         Kin.I and Kin.II data sets.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.07
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HE3_EXP(PR)
      IMPLICIT NONE
      COMMON /HE3EXP/ P_HE3,S_HE3
C
      DOUBLE PRECISION P_HE3(36),S_HE3(36),PR,RINTERPEXP
C
      HE3_EXP = RINTERPEXP(P_HE3,S_HE3,PR,36)   !in fm^3
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       3HE:  SALME (READS THE FILE)
C
C         Reference:    A. Kievsky, E. Pace, G. Salme and
C                       M. Viviani, Phys. C56, 64 (1997) and
C                       J. Templon, priv. comm. (1999).
C
C         Note the units are fm^3 for the bound state (2-body breakup)
C         and fm^4 for the continuum (3-body breakup).
C
C       Normalization:  4pi x Integral(dp   p^2 momdist) = 1.326 (2-body)
C                       4pi x Integral(dpdE p^2 momdist) = 0.638 (3-body)
C------------------------------------------------------------------------------
C
      SUBROUTINE HE3READ_SALME(BOUND)
C
      IMPLICIT NONE
C
      COMMON /HE3_SALME/ N_PM,N_EM,E,P,SPECT_2BBU,SPECT_CONT
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION E(139),P(81)
      DOUBLE PRECISION SPECT_2BBU(81),SPECT_CONT(81,139)
      INTEGER I,J,N_PM,N_EM,N_DAT_DIR
      LOGICAL BOUND
C
      IF(BOUND) THEN
         FILENAME = DAT_DIR(1:N_DAT_DIR)//'/he3-salme-2bbu.dat'
      ELSE
         FILENAME = DAT_DIR(1:N_DAT_DIR)//'/he3-salme-cont.dat'
      ENDIF
C
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
C
      IF(BOUND) THEN
         READ(3,*)N_PM
         DO I=1,N_PM
            READ(3,*) P(I),SPECT_2BBU(I)
         ENDDO
      ELSE
         READ(3,*)N_PM,N_EM
         DO I=1,N_PM
            READ(3,*) P(I)
            DO J=1,N_EM
               READ(3,*) E(J),SPECT_CONT(I,J)
            ENDDO
         ENDDO
      ENDIF
C
      CLOSE(UNIT=3)
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       3HE:  SALME SPECTRAL FCN (INTERPOLATES TO GET VALUE)
C
C       Bound state:
C           A linear interpolation on the log of S(P) is performed.
C       Continuum:
C           First a linear interpolation in energy is performed
C           for each momentum value straddling the kinematic point.
C           Then a linear interpolation on the logs of these
C           two values is performed.  Thus, the energy variable is
C           interpolated linearly, whereas the momentum variable is
C           interpolated logarithmically.
C
C         Reference:    A. Kievsky, E. Pace, G. Salme and
C                       M. Viviani, Phys. C56, 64 (1997) and
C                       J. Templon, priv. comm. (1999).
C
C         Note the units are fm^3 for the bound state (2-body breakup)
C         and fm^4 for the continuum (3-body breakup).
C
C       Normalization:  4pi x Integral(dp   p^2 momdist) = 1.326 (2-body)
C                       4pi x Integral(dpdE p^2 momdist) = 0.638 (3-body)
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HE3SPEC_SALME(BOUND,EM,PM)
C
      IMPLICIT NONE
C
      COMMON /HE3_SALME/ N_PM,N_EM,E,P,SPECT_2BBU,SPECT_CONT
C
      DOUBLE PRECISION E(139),P(81)
      DOUBLE PRECISION SPECT_2BBU(81),SPECT_CONT(81,139)
      DOUBLE PRECISION EM,PM,FRAC_E,FRAC_P
      DOUBLE PRECISION S1,S2,SE1,SE2,S11,S12,S21,S22,HBARC,HBARC3
      INTEGER I,J,N_PM,N_EM
      LOGICAL BOUND
C
      PARAMETER (HBARC  = 197.327D0)  ! Hbar*C    in  MeV-fm
      PARAMETER (HBARC3 = 7.68369D6)  !(Hbar*C)^3 in (MeV-fm)^3
C
      IF(BOUND) THEN
        DO I=1,N_PM
          IF(PM.LE.P(I)) GO TO 301
        ENDDO
        HE3SPEC_SALME = 0.D0   ! NOT TABULATED UP TO THIS MOMENTUM
        RETURN
      ELSE
        DO I=1,N_PM
          IF(PM.LE.P(I)) GOTO 250
        ENDDO
        HE3SPEC_SALME = 0.D0   ! NOT TABULATED UP TO THIS MOMENTUM
        RETURN
 250    DO J=1,N_EM
          IF(EM.LE.E(J)) GOTO 301
        ENDDO
        HE3SPEC_SALME = 0.D0   ! NOT TABULATED UP TO THIS ENERGY
        RETURN
      ENDIF
C
 301  CONTINUE        ! Successful:  found indices I and J
C
      IF(BOUND) THEN
         IF(I.EQ.1) THEN
          HE3SPEC_SALME = SPECT_2BBU(1)
	ELSE
          FRAC_P = (PM-P(I-1)) / (P(I)-P(I-1))
          S1 = LOG(SPECT_2BBU(I-1))
          S2 = LOG(SPECT_2BBU(I))
          HE3SPEC_SALME = FRAC_P * (S2-S1) + S1
          HE3SPEC_SALME = EXP(HE3SPEC_SALME)
        ENDIF	  
      ELSE
          IF (J.EQ.1) THEN
             HE3SPEC_SALME = 0.D0   ! Do not interpolate for EM < E(1)
             RETURN
          ENDIF
          IF (I.EQ.1) THEN
             FRAC_E = (EM-E(J-1)) / (E(J)-E(J-1))
             S21 = SPECT_CONT(I,J-1)
             S22 = SPECT_CONT(I,J)
             HE3SPEC_SALME = FRAC_E * (S22-S21) + S21
          ENDIF
	  IF (I.GT.1 .AND. J.GT.1) THEN  ! interpolation in EM and PM
            FRAC_E = (EM-E(J-1)) / (E(J)-E(J-1))
            FRAC_P = (PM-P(I-1)) / (P(I)-P(I-1))
            S11 = SPECT_CONT(I-1,J-1)
            S12 = SPECT_CONT(I-1,J)
            S21 = SPECT_CONT(I,J-1)
            S22 = SPECT_CONT(I,J)
            SE1 = LOG(FRAC_E * (S12-S11) + S11)
            SE2 = LOG(FRAC_E * (S22-S21) + S21)
            HE3SPEC_SALME = EXP(FRAC_P *(SE2-SE1) + SE1)
          ENDIF
      ENDIF
C
      IF(.NOT. BOUND) THEN
        HE3SPEC_SALME = HE3SPEC_SALME * HBARC
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       4HE:  CIOFI DEGLI ATTI, ET AL.
C
C       Note for Pr > 225 MeV/c this parameterization yields
C       unphysical results.  Therefore, the momentum distribution
C       is set to zero above this value.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C                       where the integral extends from 0 to 225 MeV/c
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HE4_CIOFI(PR)
      IMPLICIT NONE
C
      INCLUDE 'summary.cmn'
C
      DOUBLE PRECISION PR,PM2,SPM
C
      IF(PR .GT. 225.d0) THEN
CXXX    TYPE *,' Pr > 225 MeV/c --> setting mom-dist to zero '
        I_SPECT_OOR = I_SPECT_OOR + 1
        HE4_CIOFI = 0.D0
      ELSE
        PM2 = (PR/197.3286)**2                  !in fm^2
        SPM = 1.0354-(1.3664*PM2)+(0.271*PM2*PM2)
        HE4_CIOFI = (10.**SPM)/(2.*3.141592654)
      ENDIF
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       4HE:  ARGONNE + MODEL VII (SET UP THE ARRAYS)
C
C         These parameters are for the t+p breakup channel
C         using the Argonne potential.
C
C         Reference:    R. Schiavilla, V.R. Pandharipande and
C                       R.B. Wiringa, Nucl. Phys. A449, 219 (1986)
C                       and private communication.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.92
C                       where the integral extends to 430 MeV/c.
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_HE4_ARGONNE
      IMPLICIT NONE
      COMMON /HE4ARG/ P,Y
C
      DOUBLE PRECISION P(17),S(17),Y(17)
      INTEGER I
C
      DATA S / 54.79,51.68,43.68,33.77,24.49,17.03,11.51,7.61,4.94,
     #          3.14,1.94,1.15,0.64,0.31,0.11,-0.007,-0.08/
C
      DO I = 1,17
        P(I) = (I-1)*0.16*0.19733
        Y(I) = 2.*(S(I)**2)/((2.*3.141592654)**4)/1.65
      ENDDO
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       4HE:  ARGONNE + MODEL VII (INTERPOLATES TO GET VALUE)
C
C         These parameters are for the t+p breakup channel
C         using the Argonne potential.
C
C         Reference:    R. Schiavilla, V.R. Pandharipande and
C                       R.B. Wiringa, Nucl. Phys. A449, 219 (1986)
C                       and private communication.
C
C         Note:  this parameterization becomes unreasonable beyond
C                Pr=430 MeV/c; thus, the distribution is set to 0
C                beyond this.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.92
C                       where the integral extends to 430 MeV/c.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HE4_ARGONNE(PR)
      IMPLICIT NONE
C
      INCLUDE 'summary.cmn'
C
      COMMON /HE4ARG/ P,Y
C
      DOUBLE PRECISION P(17),Y(17),PR,PM,RINTERPEXP
C
      IF(PR .GT. 430.D0) THEN
CXXX    TYPE *,' Pr > 430 MeV/c --> setting distribution to 0 '
        I_SPECT_OOR = I_SPECT_OOR + 1
        HE4_ARGONNE = 0.D0
      ELSE
        PM = PR*0.001D0                         !convert to GeV/c
        HE4_ARGONNE = RINTERPEXP(P,Y,PM,17)
      ENDIF
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       4HE:  URBANA + MODEL VII (SET UP THE ARRAYS)
C
C         These parameters are for the t+p breakup channel
C         using the Urbana potential.
C
C         Reference:    R. Schiavilla, V.R. Pandharipande and
C                       R.B. Wiringa, Nucl. Phys. A449, 219 (1986)
C                       and private communication.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C                       where the integral extends to 430 MeV/c.
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_HE4_URBANA
      IMPLICIT NONE
      COMMON /HE4URB/ P,Y
C
      DOUBLE PRECISION P(17),S(17),Y(17)
      INTEGER I
C
      DATA S / 46.32,44.27,38.78,31.46,24.0,17.50,12.31,8.40,5.56,
     #          3.56,2.19,1.29,0.71,0.36,0.14,0.007,-0.07 /
C
      DO I = 1,17
        P(I) = (I-1)*0.16*0.19733
        Y(I) = 2.*(S(I)**2)/((2.*3.141592654)**4)/1.54
      ENDDO
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       4HE:  URBANA + MODEL VII (INTERPOLATES TO GET VALUE)
C
C         These parameters are for the t+p breakup channel
C         using the Urbana potential.
C
C         Reference:    R. Schiavilla, V.R. Pandharipande and
C                       R.B. Wiringa, Nucl. Phys. A449, 219 (1986)
C                       and private communication.
C
C         Note:  this parameterization becomes unreasonable beyond
C                Pr=430 MeV/c; thus, the distribution is set to 0
C                beyond this.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C                       where the integral extends to 430 MeV/c.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HE4_URBANA(PR)
      IMPLICIT NONE
C
      INCLUDE 'summary.cmn'
C
      COMMON /HE4URB/ P,Y
C
      DOUBLE PRECISION P(17),Y(17),PR,PM,RINTERPEXP
C
      IF(PR .GT. 430.D0) THEN
CXXX    TYPE *,' Pr > 430 MeV/c --> setting distribution to 0 '
        I_SPECT_OOR = I_SPECT_OOR + 1
        HE4_URBANA = 0.D0
      ELSE
        PM = PR*0.001D0                         !convert to GeV/c
        HE4_URBANA = RINTERPEXP(P,Y,PM,17)
      ENDIF
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       4HE:  NIKHEF EXPT. J.v.d. BRAND ET AL. (READS THE DATA FILE)
C                    q=431 MeV/c
C
C         These data correspond to the t+p breakup channel only.
C
C         Reference:    J.F.J. van den Brand et al., Phys. Rev.
C                       Lett. 60, 2006 (1988) and private communication.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.139
C                       where the integral extends from 0 to 207.5 MeV/c
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_HE4_EXP_Q1
      IMPLICIT NONE
      COMMON /HE4EXP_Q1/ P_HE4_Q1,S_HE4_Q1
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_HE4_Q1(40),S_HE4_Q1(40)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/4he_nik_q1.exp'
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(///)
      DO I=1,40
        READ (3,*) P_HE4_Q1(I),S_HE4_Q1(I)
        S_HE4_Q1(I) = S_HE4_Q1(I)*((0.19733)**3)   !convert to fm^3
      ENDDO
      CLOSE(UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       4HE:  NIKHEF EXPT. J.v.d. BRAND ET AL. (INTERPOLATES TO GET VALUE)
C                    q=431 MeV/c
C
C         These data correspond to the t+p breakup channel only.
C
C         Reference:    J.F.J. van den Brand et al., Phys. Rev.
C                       Lett. 60, 2006 (1988) and private communication.
C
C       Note that this data set extends only to Pr = 207.5 MeV/c
C       Therefore, the momentum distribution is set to zero above
C       this value.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.139
C                       where the integral extends from 0 to 207.5 MeV/c
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HE4_EXP_Q1(PR)
      IMPLICIT NONE
C
      INCLUDE 'summary.cmn'
C
      COMMON /HE4EXP_Q1/ P_HE4_Q1,S_HE4_Q1
C
      DOUBLE PRECISION P_HE4_Q1(40),S_HE4_Q1(40),PR,RINTERPEXP
C
      IF(PR .GT. 207.5D0) THEN
CXXX    TYPE *,' Pr > 207.5 MeV/c --> setting mom-dist to zero '
        I_SPECT_OOR = I_SPECT_OOR + 1
        HE4_EXP_Q1 = 0.D0
      ELSE
        HE4_EXP_Q1 = RINTERPEXP(P_HE4_Q1,S_HE4_Q1,PR,40)  !in fm^3
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       4HE:  NIKHEF EXPT. J.v.d. BRAND ET AL. (READS THE DATA FILE)
C                    q=250 MeV/c
C
C         These data correspond to the t+p breakup channel only.
C
C         Reference:    J.F.J. van den Brand et al., Phys. Rev.
C                       Lett. 60, 2006 (1988) and private communication.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.405 (0.998)
C                       where integral extends from 112.5 (0.0) to 362.5 MeV/c
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_HE4_EXP_Q2
      IMPLICIT NONE
      COMMON /HE4EXP_Q2/ P_HE4_Q2,S_HE4_Q2
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_HE4_Q2(51),S_HE4_Q2(51)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/4he_nik_q2.exp'
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(///)
      DO I=1,51
        READ (3,*) P_HE4_Q2(I),S_HE4_Q2(I)
        S_HE4_Q2(I) = S_HE4_Q2(I)*((0.19733)**3)   !convert to fm^3
      ENDDO
      CLOSE(UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       4HE:  NIKHEF EXPT. J.v.d. BRAND ET AL. (INTERPOLATES TO GET VALUE)
C                    q=250 MeV/c
C
C         These data correspond to the t+p breakup channel only.
C
C         Reference:    J.F.J. van den Brand et al., Phys. Rev.
C                       Lett. 60, 2006 (1988) and private communication.
C
C       Note that this data set extends from Pr = 112.5 to 362.5 MeV/c.
C       Extrapolations give fairly reasonable values for Pr < 112.5 MeV/c
C       (within a factor of two of the q=431 MeV/c data set).  However,
C       extrapolation beyond 362.5 MeV/c yields unreasonable values.
C       Therefore, the momentum distribution is set to zero above
C       this range.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.405 (0.998)
C                       where integral extends from 112.5 (0.0) to 362.5 MeV/c
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HE4_EXP_Q2(PR)
      IMPLICIT NONE
C
      INCLUDE 'summary.cmn'
C
      COMMON /HE4EXP_Q2/ P_HE4_Q2,S_HE4_Q2
C
      DOUBLE PRECISION P_HE4_Q2(51),S_HE4_Q2(51),PR,RINTERPEXP
C
      IF(PR .GT. 362.5D0) THEN
CXXX    TYPE *,' Pr > 362.5 MeV/c --> setting mom-dist to zero '
        I_SPECT_OOR = I_SPECT_OOR + 1
        HE4_EXP_Q2 = 0.D0
      ELSE
        HE4_EXP_Q2 = RINTERPEXP(P_HE4_Q2,S_HE4_Q2,PR,51)  !in fm^3
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C BEGIN MODIFICATION 12
C------------------------------------------------------------------------------
C       010216 KF - add subroutine
C
C       16O:   P1/2 - UDIAS  (READS THE DISTRIBUTION FILE)
C       integration range is [0, 895] MeV/c
C       reference:  private communciaton
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.532
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_O16P12_UD
      IMPLICIT NONE
      COMMON /O16P12_DW/ P_O16P12_UD,S_O16P12_UD
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_O16P12_UD(180),S_O16P12_UD(180)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/16o_ju_dw.dat'
      WRITE(6,*) ''
      WRITE(6,*) 'PLEASE NOTE:  THIS RDWIA BOUND-STATE MOMENTUM'
      WRITE(6,*) 'DISTRIBUTION WAS CALCULATED SPECIFICALLY FOR'
      WRITE(6,*) 'PARALLEL KINEMATICS, E0 = 4045 MEV, TP = 432 MEV.'
      WRITE(6,*) ''
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(//)
      DO I = 1,180
        READ (3,*) P_O16P12_UD(I),S_O16P12_UD(I)
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C       010216 KF - add subroutine
C
C       16O:   P1/2 - UDIAS  (INTERPOLATES TO GET VALUE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.532
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION O16P12_UD(PR)
      IMPLICIT NONE
      COMMON /O16P12_DW/ P_O16P12_UD,S_O16P12_UD
C
      DOUBLE PRECISION P_O16P12_UD(180),S_O16P12_UD(180),PR,RINTERPEXP
C
      O16P12_UD = RINTERPEXP(P_O16P12_UD,S_O16P12_UD,PR,180)
C
      RETURN
      END
C
C---------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C       010216 KF - add subroutine
C
C       16O:   P3/2 - UDIAS  (READS THE DISTRIBUTION FILE)
C       integration range is [0, 895] MeV/c
C       reference:  private communciaton
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.205
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_O16P32_UD
      IMPLICIT NONE
      COMMON /O16P32_DW/ P_O16P32_UD,S_O16P32_UD
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_O16P32_UD(180),S_O16P32_UD(180)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/16o_ju_dw.dat'
      WRITE(6,*) ''
      WRITE(6,*) 'PLEASE NOTE:  THIS RDWIA BOUND-STATE MOMENTUM'
      WRITE(6,*) 'DISTRIBUTION WAS CALCULATED SPECIFICALLY FOR'
      WRITE(6,*) 'PARALLEL KINEMATICS, E0 = 4045 MEV, TP = 432 MEV.'
      WRITE(6,*) ''
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(182(/)/)
      DO I = 1,180
        READ (3,*) P_O16P32_UD(I),S_O16P32_UD(I)
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C       010216 KF - add subroutine
C
C       16O:   P3/2 - UDIAS  (INTERPOLATES TO GET VALUE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.205
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION O16P32_UD(PR)
      IMPLICIT NONE
      COMMON /O16P32_DW/ P_O16P32_UD,S_O16P32_UD
C
      DOUBLE PRECISION P_O16P32_UD(180),S_O16P32_UD(180),PR,RINTERPEXP
C
      O16P32_UD = RINTERPEXP(P_O16P32_UD,S_O16P32_UD,PR,180)
C
      RETURN
      END
C
C---------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C       010216 KF - add subroutine
C
C       16O:   S1/2 - UDIAS  (READS THE DISTRIBUTION FILE)
C       integration range is [0, 895] MeV/c
C       reference:  private communciaton
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.464
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_O16S12_UD
      IMPLICIT NONE
      COMMON /O16S12_DW/ P_O16S12_UD,S_O16S12_UD
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_O16S12_UD(180),S_O16S12_UD(180)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/16o_ju_dw.dat'
      WRITE(6,*) ''
      WRITE(6,*) 'PLEASE NOTE:  THIS RDWIA BOUND-STATE MOMENTUM'
      WRITE(6,*) 'DISTRIBUTION WAS CALCULATED SPECIFICALLY FOR'
      WRITE(6,*) 'PARALLEL KINEMATICS, E0 = 4045 MEV, TP = 432 MEV.'
      WRITE(6,*) ''
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(363(/)/)
      DO I = 1,180
        READ (3,*) P_O16S12_UD(I),S_O16S12_UD(I)
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C       010216 KF - add subroutine
C
C       16O:   S1/2 - UDIAS  (INTERPOLATES TO GET VALUE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.464
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION O16S12_UD(PR)
      IMPLICIT NONE
      COMMON /O16S12_DW/ P_O16S12_UD,S_O16S12_UD
C
      DOUBLE PRECISION P_O16S12_UD(180),S_O16S12_UD(180),PR,RINTERPEXP
C
      O16S12_UD = RINTERPEXP(P_O16S12_UD,S_O16S12_UD,PR,180)
C
      RETURN
      END
C
C---------------------------------------------------------------------------
C END MODIFICATION 12
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       16O:   P1/2 - VAN ORDEN  (READS THE DISTRIBUTION FILE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_O16P12_VO
      IMPLICIT NONE
      COMMON /O16P12/ P_O16P12_VO,S_O16P12_VO
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_O16P12_VO(93),S_O16P12_VO(93)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/16o_vo.dat'
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(//)
      DO I = 1,93
        READ (3,*) P_O16P12_VO(I),S_O16P12_VO(I)
        P_O16P12_VO(I) = P_O16P12_VO(I)*1000.D0   !convert to MeV/c
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       16O:   P1/2 - VAN ORDEN  (INTERPOLATES TO GET VALUE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION O16P12_VO(PR)
      IMPLICIT NONE
      COMMON /O16P12/ P_O16P12_VO,S_O16P12_VO
C
      DOUBLE PRECISION P_O16P12_VO(93),S_O16P12_VO(93),PR,RINTERPEXP
C
      O16P12_VO = RINTERPEXP(P_O16P12_VO,S_O16P12_VO,PR,93)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       16O:  P3/2 - VAN ORDEN  (READS THE DISTRIBUTION FILE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 4
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_O16P32_VO
      IMPLICIT NONE
      COMMON /O16P32/ P_O16P32_VO,S_O16P32_VO
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
C
      DOUBLE PRECISION P_O16P32_VO(93),S_O16P32_VO(93)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/16o_vo.dat'
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(96(/))
      DO I = 1,93
        READ (3,*) P_O16P32_VO(I),S_O16P32_VO(I)
        P_O16P32_VO(I) = P_O16P32_VO(I)*1000.D0   !convert to MeV/c
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       16O:  P3/2 - VAN ORDEN  (INTERPOLATES TO GET VALUE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 4
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION O16P32_VO(PR)
      IMPLICIT NONE
      COMMON /O16P32/ P_O16P32_VO,S_O16P32_VO
C
      DOUBLE PRECISION P_O16P32_VO(93),S_O16P32_VO(93),PR,RINTERPEXP
C
      O16P32_VO = RINTERPEXP(P_O16P32_VO,S_O16P32_VO,PR,93)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       16O:  S1/2 - VAN ORDEN  (READS THE DISTRIBUTION FILE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_O16S12_VO
      IMPLICIT NONE
      COMMON /O16S12/ P_O16S12_VO,S_O16S12_VO
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_O16S12_VO(93),S_O16S12_VO(93)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/16o_vo.dat'
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(190(/))
      DO I = 1,93
        READ (3,*) P_O16S12_VO(I),S_O16S12_VO(I)
        P_O16S12_VO(I) = P_O16S12_VO(I)*1000.D0   !convert to MeV/c
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       16O:  S1/2 - VAN ORDEN  (INTERPOLATES TO GET VALUE)
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION O16S12_VO(PR)
      IMPLICIT NONE
      COMMON /O16S12/ P_O16S12_VO,S_O16S12_VO
C
      DOUBLE PRECISION P_O16S12_VO(93),S_O16S12_VO(93),PR,RINTERPEXP
C
      O16S12_VO = RINTERPEXP(P_O16S12_VO,S_O16S12_VO,PR,93)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C BEGIN MODIFICATION 13
C---------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       001010 KF - add 208Pb from Lapikas
C       CDWIA for E0 = 674 MeV, Tp = 160 MeV, parallel kinematics
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  3S1/2 - LAPIKAS  (READS THE DISTRIBUTION FILE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.780
C
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_PB2S12_LL
      IMPLICIT NONE
      COMMON /PB2S12/ P_PB2S12_LL,S_PB2S12_LL
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_PB2S12_LL(33),S_PB2S12_LL(33)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/208pb_ll_dw.dat'
      WRITE(6,*) ''
      WRITE(6,*) 'PLEASE NOTE:  THIS NRDWIA BOUND-STATE MOMENTUM'
      WRITE(6,*) 'DISTRIBUTION WAS CALCULATED SPECIFICALLY FOR'
      WRITE(6,*) 'PARALLEL KINEMATICS, E0 = 674 MEV, TP = 160 MEV.'
      WRITE(6,*) ''
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(//)
      DO I = 1,33
        READ (3,*) P_PB2S12_LL(I),S_PB2S12_LL(I)
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  3S1/2 - LAPIKAS  (INTERPOLATES TO GET VALUE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 0.780
C
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION PB2S12_LL(PR)
      IMPLICIT NONE
      COMMON /PB2S12/ P_PB2S12_LL,S_PB2S12_LL
C
      DOUBLE PRECISION P_PB2S12_LL(33),S_PB2S12_LL(33),PR,RINTERPEXP
C
      IF (PR .LE. 300.) THEN
        PB2S12_LL = RINTERPEXP(P_PB2S12_LL,S_PB2S12_LL,PR,33)
      ELSE
        PB2S12_LL = 0.
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  2D3/2 - LAPIKAS  (READS THE DISTRIBUTION FILE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.232
C
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_PB2D32_LL
      IMPLICIT NONE
      COMMON /PB2D32/ P_PB2D32_LL,S_PB2D32_LL
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_PB2D32_LL(33),S_PB2D32_LL(33)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/208pb_ll_dw.dat'
      WRITE(6,*) ''
      WRITE(6,*) 'PLEASE NOTE:  THIS NRDWIA BOUND-STATE MOMENTUM'
      WRITE(6,*) 'DISTRIBUTION WAS CALCULATED SPECIFICALLY FOR'
      WRITE(6,*) 'PARALLEL KINEMATICS, E0 = 674 MEV, TP = 160 MEV.'
      WRITE(6,*) ''
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(36(/))
      DO I = 1,33
        READ (3,*) P_PB2D32_LL(I),S_PB2D32_LL(I)
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  2D3/2 - LAPIKAS  (INTERPOLATES TO GET VALUE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.232
C
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION PB2D32_LL(PR)
      IMPLICIT NONE
      COMMON /PB2D32/ P_PB2D32_LL,S_PB2D32_LL
C
      DOUBLE PRECISION P_PB2D32_LL(33),S_PB2D32_LL(33),PR,RINTERPEXP
C
      IF (PR .LE. 300.) THEN
        PB2D32_LL = RINTERPEXP(P_PB2D32_LL,S_PB2D32_LL,PR,33)
      ELSE
        PB2D32_LL = 0.
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  1H11/2 - LAPIKAS  (READS THE DISTRIBUTION FILE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2.932
C
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_PB2H12_LL
      IMPLICIT NONE
      COMMON /PB2H12/ P_PB2H12_LL,S_PB2H12_LL
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_PB2H12_LL(33),S_PB2H12_LL(33)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/208pb_ll_dw.dat'
      WRITE(6,*) ''
      WRITE(6,*) 'PLEASE NOTE:  THIS NRDWIA BOUND-STATE MOMENTUM'
      WRITE(6,*) 'DISTRIBUTION WAS CALCULATED SPECIFICALLY FOR'
      WRITE(6,*) 'PARALLEL KINEMATICS, E0 = 674 MEV, TP = 160 MEV.'
      WRITE(6,*) ''
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(70(/))
      DO I = 1,33
        READ (3,*) P_PB2H12_LL(I),S_PB2H12_LL(I)
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  1H11/2 - LAPIKAS  (INTERPOLATES TO GET VALUE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2.932
C
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION PB2H12_LL(PR)
      IMPLICIT NONE
      COMMON /PB2H12/ P_PB2H12_LL,S_PB2H12_LL
C
      DOUBLE PRECISION P_PB2H12_LL(33),S_PB2H12_LL(33),PR,RINTERPEXP
C
      IF (PR .LE. 300.) THEN
        PB2H12_LL = RINTERPEXP(P_PB2H12_LL,S_PB2H12_LL,PR,33)
      ELSE
        PB2H12_LL = 0.
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  2D5/2 - LAPIKAS  (READS THE DISTRIBUTION FILE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.954
C
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_PB2D52_LL
      IMPLICIT NONE
      COMMON /PB2D52/ P_PB2D52_LL,S_PB2D52_LL
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_PB2D52_LL(33),S_PB2D52_LL(33)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/208pb_ll_dw.dat'
      WRITE(6,*) ''
      WRITE(6,*) 'PLEASE NOTE:  THIS NRDWIA BOUND-STATE MOMENTUM'
      WRITE(6,*) 'DISTRIBUTION WAS CALCULATED SPECIFICALLY FOR'
      WRITE(6,*) 'PARALLEL KINEMATICS, E0 = 674 MEV, TP = 160 MEV.'
      WRITE(6,*) ''
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(104(/))
      DO I = 1,33
        READ (3,*) P_PB2D52_LL(I),S_PB2D52_LL(I)
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  2D5/2 - LAPIKAS  (INTERPOLATES TO GET VALUE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1.954
C
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION PB2D52_LL(PR)
      IMPLICIT NONE
      COMMON /PB2D52/ P_PB2D52_LL,S_PB2D52_LL
C
      DOUBLE PRECISION P_PB2D52_LL(33),S_PB2D52_LL(33),PR,RINTERPEXP
C
      IF (PR .LE. 300.) THEN
        PB2D52_LL = RINTERPEXP(P_PB2D52_LL,S_PB2D52_LL,PR,33)
      ELSE
        PB2D52_LL = 0.
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  1G7/2 - LAPIKAS  (READS THE DISTRIBUTION FILE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2.669
C
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_PB1G72_LL
      IMPLICIT NONE
      COMMON /PB1G72/ P_PB1G72_LL,S_PB1G72_LL
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_PB1G72_LL(33),S_PB1G72_LL(33)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/208pb_ll_dw.dat'
      WRITE(6,*) ''
      WRITE(6,*) 'PLEASE NOTE:  THIS NRDWIA BOUND-STATE MOMENTUM'
      WRITE(6,*) 'DISTRIBUTION WAS CALCULATED SPECIFICALLY FOR'
      WRITE(6,*) 'PARALLEL KINEMATICS, E0 = 674 MEV, TP = 160 MEV.'
      WRITE(6,*) ''
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(138(/))
      DO I = 1,33
        READ (3,*) P_PB1G72_LL(I),S_PB1G72_LL(I)
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  1G7/2 - LAPIKAS  (INTERPOLATES TO GET VALUE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2.669
C
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION PB1G72_LL(PR)
      IMPLICIT NONE
      COMMON /PB1G72/ P_PB1G72_LL,S_PB1G72_LL
C
      DOUBLE PRECISION P_PB1G72_LL(33),S_PB1G72_LL(33),PR,RINTERPEXP
C
      IF (PR .LE. 300.) THEN
        PB1G72_LL = RINTERPEXP(P_PB1G72_LL,S_PB1G72_LL,PR,33)
      ELSE
        PB1G72_LL = 0.
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  1G9/2 - LAPIKAS  (READS THE DISTRIBUTION FILE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2.843
C
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_PB1G92_LL
      IMPLICIT NONE
      COMMON /PB1G92/ P_PB1G92_LL,S_PB1G92_LL
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100  DAT_DIR
      CHARACTER*200  FILENAME
      DOUBLE PRECISION P_PB1G92_LL(33),S_PB1G92_LL(33)
      INTEGER I,N_DAT_DIR
C
      FILENAME = DAT_DIR(1:N_DAT_DIR)//'/208pb_ll_dw.dat'
      WRITE(6,*) ''
      WRITE(6,*) 'PLEASE NOTE:  THIS NRDWIA BOUND-STATE MOMENTUM'
      WRITE(6,*) 'DISTRIBUTION WAS CALCULATED SPECIFICALLY FOR'
      WRITE(6,*) 'PARALLEL KINEMATICS, E0 = 674 MEV, TP = 160 MEV.'
      WRITE(6,*) ''
      OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD')
      READ (3,101)
  101 FORMAT(172(/))
      DO I = 1,33
        READ (3,*) P_PB1G92_LL(I),S_PB1G92_LL(I)
      ENDDO
      CLOSE (UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
C
C       001010 KF - add subroutine
C
C       208Pb:  1G9/2 - LAPIKAS  (INTERPOLATES TO GET VALUE)
C
C       Louk Lapikas, private communication
C       integration range is [0, 300] MeV/c
C       Normalization:  4pi x Integral(dp p^2 momdist) = 2.843
C
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION PB1G92_LL(PR)
      IMPLICIT NONE
      COMMON /PB1G92/ P_PB1G92_LL,S_PB1G92_LL
C
      DOUBLE PRECISION P_PB1G92_LL(33),S_PB1G92_LL(33),PR,RINTERPEXP
C
      IF (PR .LE. 300.) THEN
        PB1G92_LL = RINTERPEXP(P_PB1G92_LL,S_PB1G92_LL,PR,33)
      ELSE
        PB1G92_LL = 0.
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C---------------------------------------------------------------------------
C END MODIFICATION 13
C
C------------------------------------------------------------------------------
C       SUBROUTINE GET_OSC_PARAM
C
C       AUTHOR:  P.E. ULMER
C       DATE:    AUG-30-1990
C
C       PURPOSE: Get oscillator parameter for Harmonic Oscillator
C                momentum distribution.
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_OSC_PARAM
      IMPLICIT NONE
      COMMON /OSC/B
C
      DOUBLE PRECISION B
C
      WRITE(6,10)' Enter Osc. param. (fm) (def.=1.67)>'
   10 FORMAT(A)
      READ(5,20) B
   20 FORMAT(F12.5)
      IF(B .EQ. 0.) B = 1.67D0
C
      RETURN
      END
C
C -------------------------------------------------------------------
C -------------------------------------------------------------------
C       Harmonic oscillator momentum distribution for 1S, 1P, 2S, 1d,
C       2P, 1f, 3S, 2d and 1g particles in fm^3.
C
C       Normalization:  4pi x Integral(dp p^2 momdist) = 1
C -------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION HARM_OSC(PR)
      IMPLICIT NONE
C
      COMMON /PHYSVAR/ IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
      COMMON /OSC/B
C
      DOUBLE PRECISION PR,B,SQRT_PI,PI,X,ARG,MPLR
      INTEGER   IELAST_OPT,IPHYS_OPT,ISPEC_OPT,IMODEL
C
      PARAMETER (PI = 3.1415926D0, SQRT_PI = 1.7724539D0)
C
      X    = PR/197.3286D0         !convert to fm^-1
      ARG  = X*B
      ARG  = ARG * ARG
      MPLR = B**3 * EXP(-ARG)/SQRT_PI/PI
C
      IF (ISPEC_OPT .EQ. 150)  THEN                               !1s-state   
           HARM_OSC = MPLR                                             
C
      ELSEIF(ISPEC_OPT .EQ. 151) THEN                             !1p-state
           HARM_OSC = MPLR * (2.D0/3.D0) * ARG
C
      ELSEIF(ISPEC_OPT .EQ. 152) THEN                             !2s-state
          HARM_OSC = MPLR * (2.D0/3.D0) * (3.D0/2.D0 - ARG)**2
C
      ELSEIF(ISPEC_OPT .EQ. 153) THEN                             !1d-state
          HARM_OSC = MPLR * (4.D0/15.D0) * ARG**2
C
      ELSEIF(ISPEC_OPT .EQ. 154) THEN                             !2p-state
          HARM_OSC = MPLR * (4.D0/15.D0) * ARG *
     #            (5.D0/2.D0 - ARG)**2
C
      ELSEIF(ISPEC_OPT .EQ. 155) THEN                             !1f-state 
          HARM_OSC = MPLR * (8.D0/105.D0) * ARG**3
C
      ELSEIF(ISPEC_OPT .EQ. 156) THEN                             !3s-state  
          HARM_OSC = MPLR * (8.D0/15.D0)*(15.D0/8.D0 -
     #            5.D0* ARG/2.D0 + ARG**2/2)**2
C
      ELSEIF(ISPEC_OPT .EQ. 157) THEN                             !2d-state 
          HARM_OSC = MPLR * (8.D0/105.D0) * (ARG *
     #              (7.D0/2.D0 - ARG))**2
C
      ELSE                                                       !1g-state
          HARM_OSC = MPLR * (16.D0 * ARG**4 / 945.D0)
      ENDIF
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE GET_SPECTRAL
C
C       AUTHOR:  P.E. ULMER
C       DATE:    AUG-30-1990
C
C       PURPOSE: Get spectral function, S(p,E), from user-supplied file.
C
C       The spectral function should be NORMALIZED so that:
C               the integral of S(p,E)p^2 dp dE x 4pi = 1
C
C       UNITS:  fm^3 MeV^-1
C
C       Further, it is assumed that each value of p has
C       the same number of E elements (rectangular array)
C       and that the spacing between successive p values
C       and between successive E values is constant.
C
C       The first line in the external file should specify
C       the number of p elements followed by the number of E elements.
C       For the format of the remainder of the file, see the code.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_SPECTRAL
      IMPLICIT NONE
      COMMON /SPECVAR_D/ E_ARR,P_ARR,SPECT_FCN,DELTA_P,DELTA_E
      COMMON /SPECVAR_I/ NBINS_P,NBINS_E
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100 DAT_DIR
      DOUBLE PRECISION E_ARR(200),P_ARR(200),SPECT_FCN(200,200)
      DOUBLE PRECISION DELTA_P,DELTA_E
      INTEGER NBINS_P,NBINS_E,I,J,N_DAT_DIR
      CHARACTER*200 SPECTRAL_FILE
C
      WRITE(6,10) ' Enter file for Spectral Distribution >'
   10 FORMAT(A)
      READ(5,20)SPECTRAL_FILE
   20 FORMAT(A20)
      SPECTRAL_FILE = DAT_DIR(1:N_DAT_DIR)//'/'//SPECTRAL_FILE
C
      OPEN(UNIT=3,FILE=SPECTRAL_FILE,STATUS='OLD',
     #          FORM='FORMATTED')
      READ(3,*) NBINS_P,NBINS_E
C
      DO I=1,NBINS_P
        READ(3,*) P_ARR(I)
        DO J=1,NBINS_E
          READ(3,*) E_ARR(J),SPECT_FCN(I,J)
        ENDDO
      ENDDO
      DELTA_P = P_ARR(2) - P_ARR(1)
      DELTA_E = E_ARR(2) - E_ARR(1)
C
      CLOSE(UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       FUNCTION EVAL_SPECTRAL
C
C       AUTHOR:  P.E. ULMER
C       DATE:    SEP-05-1990
C
C       PURPOSE: Evaluate spectral distribution by interpolation
C                on 2-dimensional array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION EVAL_SPECTRAL(PRMAG,MISS_M)
      IMPLICIT NONE
      COMMON /SPECVAR_D/ E_ARR,P_ARR,SPECT_FCN,DELTA_P,DELTA_E
      COMMON /SPECVAR_I/ NBINS_P,NBINS_E
C
      DOUBLE PRECISION E_ARR(200),P_ARR(200),SPECT_FCN(200,200)
      DOUBLE PRECISION PRMAG,MISS_M,DELTA_P,DELTA_E,RINTERP2D
      INTEGER NBINS_P,NBINS_E
C
      EVAL_SPECTRAL = RINTERP2D(P_ARR,E_ARR,SPECT_FCN,PRMAG,MISS_M,
     #          NBINS_P,NBINS_E,DELTA_P,DELTA_E)        !fm^3 MeV^-1
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE GET_PHISQ
C
C       AUTHOR:  P.E. ULMER
C       DATE:    AUG-30-1990
C
C       PURPOSE: Get momentum distribution, phi(p)^2, from
C                user-supplied file.
C
C       The momentum distribution should be NORMALIZED so that:
C               the integral of phi(p)^2 p^2 dp x 4pi = 1
C
C       UNITS:  fm^3
C
C       Further, it is assumed that each value of p has
C       the same number of E elements (rectangular array).
C
C       The first line in the external file should specify
C       the number of p elements followed by the number of E elements.
C       For the format of the remainder of the file, see the code.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_PHISQ
      IMPLICIT NONE
      COMMON /PHISQ_D/ PR_ARR,PHI_DAT
      COMMON /PHISQ_I/ NBINS_PR
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100 DAT_DIR
      DOUBLE PRECISION PR_ARR(200),PHI_DAT(200)
      INTEGER NBINS_PR,I,N_DAT_DIR
      CHARACTER*200 PHISQ_FILE
C
      WRITE(6,10) ' Enter file for phi^2 distribution >'
   10 FORMAT(A)
      READ(5,20)PHISQ_FILE
   20 FORMAT(A20)
      PHISQ_FILE = DAT_DIR(1:N_DAT_DIR)//'/'//PHISQ_FILE
C
      OPEN(UNIT=3,FILE=PHISQ_FILE,STATUS='OLD',
     #          FORM='FORMATTED')
C
      DO I=1,200
        READ(3,*,END=100)PR_ARR(I),PHI_DAT(I)
        NBINS_PR = I
      ENDDO
  100 CLOSE(UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       FUNCTION EVAL_PHISQ
C
C       AUTHOR:  P.E. ULMER
C       DATE:    SEP-05-1990
C
C       PURPOSE: Evaluate momentum distribution by interpolation on array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION EVAL_PHISQ(PR)
      IMPLICIT NONE
      COMMON /PHISQ_D/ PR_ARR,PHI_DAT
      COMMON /PHISQ_I/ NBINS_PR
C
      DOUBLE PRECISION PR_ARR(200),PHI_DAT(200),PR,RINTERPEXP
      INTEGER NBINS_PR
C
      EVAL_PHISQ = RINTERPEXP(PR_ARR,PHI_DAT,PR,NBINS_PR) !fm^3
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE GET_EMDIST
C
C       AUTHOR:  P.E. ULMER
C       DATE:    AUG-30-1990
C
C       PURPOSE: Get missing energy distribution for continuum.
C
C       It is assumed that the F(Em) distribution is NORMALIZED:
C               integral of f(Em) dEm = 1
C
C       UNITS:  MeV^-1
C------------------------------------------------------------------------------
C
      SUBROUTINE GET_EMDIST
      IMPLICIT NONE
      COMMON /EMDIST_D/ EM_ARR,EMDIST
      COMMON /EMDIST_I/ NBINS_EM
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      CHARACTER*100 DAT_DIR
      DOUBLE PRECISION EM_ARR(1000),EMDIST(1000)
      INTEGER NBINS_EM,I,N_DAT_DIR
      CHARACTER*200 EMFILE
C
      WRITE(6,10) ' Enter file for Em distribution >'
   10 FORMAT(A)
      READ(5,20)EMFILE
   20 FORMAT(A20)
      EMFILE = DAT_DIR(1:N_DAT_DIR)//'/'//EMFILE
      OPEN(UNIT=3,FILE=EMFILE,STATUS='OLD',FORM='FORMATTED')
      DO I=1,1000
        READ(3,*,END=88)EM_ARR(I),EMDIST(I)
        NBINS_EM = I
      ENDDO
   88 CLOSE(UNIT=3)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       FUNCTION EVAL_EMDIST
C
C       AUTHOR:  P.E. ULMER
C       DATE:    SEP-05-1990
C
C       PURPOSE: Evaluate missing mass distribution by interpolation on array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION EVAL_EMDIST(MISS_M)
      IMPLICIT NONE
      COMMON /EMDIST_D/ EM_ARR,EMDIST
      COMMON /EMDIST_I/ NBINS_EM
C
      DOUBLE PRECISION EM_ARR(1000),EMDIST(1000),MISS_M,RINTERPEXP
      INTEGER NBINS_EM
C
      EVAL_EMDIST = RINTERPEXP(EM_ARR,EMDIST,MISS_M,NBINS_EM)   !MeV^-1
C
      RETURN
      END
