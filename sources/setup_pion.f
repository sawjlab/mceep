C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       SUBROUTINE:  SETUP_PION
C
C       AUTHOR: R.W. Lourie
C       DATE:   1991
C       MODIFICATIONS:
C
C       PURPOSE:
C          Reads the file 'RES_OPT_DAT'
C          to get parameters/options for resonance electroproduction
C          subroutine ELECTRO_PROD.  (RES_OPT_DAT is the name of the
C          file given as the answer to the question asked in PHYS_CHOICE
C          when option 500 selected.)
C
C       Structure of RES_OPT_DAT:
C
C       LMAX    : Maximum partial wave in pi-N system (can be 1,2 or 3)
C       IOPT_M1P: =0 to turn off M1+ of Delta, =1 makes E1+=M1+ (pQCD),
C                 =-1 to turn off Delta quadrupole and ALL other resonances,
C                 =anything <-1 to turn off all other resonances.
C
C       PERCENT_E2: Ratio of E2/M1 for Delta.
C       PERCENT_C2: Ratio of C2/M1 for Delta.
C
C       BORN_TERMS  : =1 to use Born terms, 0 to set them to 0.
C       USE_P11_FILE: =1 if Roper amplitudes are contained in external file.
C
C       P11_FILE: Name of ext. file containing Roper A1/2 and S1/2
C                 amplitudes (if USE_P11_FILE=1). Note that two files are
C                 required but only one name (files should be
C                 named P11_FILE_A.DAT or _S.DAT).
C       AMP_AP11: Numerical value for P11 transverse amplitude (if
C                 USE_P11_FILE=0).
C       AMP_SP11: Numerical value for P11 longitudinal amplitude (if
C                 USE_P11_FILE=0).
C
C------------------------------------------------------------------------------
C
      SUBROUTINE setup_pion(ilun)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER p11_file*300,sp11_file*300,ap11_file*300
      LOGICAL   born_terms,use_p11_file
      DIMENSION sxqq2(30),axqq2(30),ap11(30),sp11(30)
      INTEGER   nchar
C
      COMMON /p11_amplitudes/sxqq2,axqq2,ap11,sp11,amp_ap11,amp_sp11
      COMMON /p33_amplitudes/percent_c2,percent_e2
      COMMON /np11_elmts/ n_sp11,n_ap11
      COMMON /control_i/iopt_m1p,lmax
      COMMON /control_l/born_terms,use_p11_file
C
      READ(ilun,*)lmax,iopt_m1p
      READ(ilun,*)percent_e2,percent_c2
      READ(ilun,*)born_terms,use_p11_file
C
      IF(use_p11_file)THEN
        READ(ilun,100) p11_file
100     FORMAT(a)
        CALL SQUEEZE(p11_file,p11_file,nchar)
      ELSE
        READ(ilun,*)amp_ap11,amp_sp11
      ENDIF
C
      IF(use_p11_file)THEN
        sp11_file = p11_file(1:nchar)//'_s.dat'
        ap11_file = p11_file(1:nchar)//'_a.dat'
        OPEN(unit=8,file=sp11_file,status='old')
        OPEN(unit=9,file=ap11_file,status='old')
        DO j=1,30
          READ(8,*,end=777)sxqq2(j),sp11(j)
          n_sp11 = j
        ENDDO
  777   CLOSE(8)
        DO j=1,30
          READ(9,*,end=778)axqq2(j),ap11(j)
          n_ap11 = j
        ENDDO
  778   CLOSE(9)
      ENDIF
C
      RETURN
      END
