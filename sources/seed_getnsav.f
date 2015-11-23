C
C ----------------------------------------------------------------------
C       Subroutine GETSEED
C
C       PURPOSE:
C           Get the initial random seeds from the file seeds_start.dat.
C           If the file doesn't exist, create it and load it
C           with default seeds.
C
C       MODIFICATIONS:
C           6/20/85   Change from two seeds to one.
C           1/23/91   Change back to two seeds with eye for future
C                     expansion.
C          10/16/92   PEU:  Now uses different file name than for
C                     ending seeds allowing replay with the same seeds.
C ----------------------------------------------------------------------
C
      SUBROUTINE GETSEED
C
      CHARACTER*100 DAT_DIR
      CHARACTER*200 INFILE
      INTEGER ISEEDS0(2)
      COMMON /RANCOM/ ISEEDS0
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
C     Read the seeds from seeds_start.dat, call RECUIN to initialize
C     the seeds and save them in RANCOM for future reference.
C     If seeds_start.dat does not exist, pick default seeds and create
C     a new file.
C
      INFILE = DAT_DIR(1:N_DAT_DIR)//'/seeds_start.dat'
      OPEN(UNIT=10,STATUS='OLD',FILE=INFILE,ERR=200)
      READ(10,*) iseeds0(1),iseeds0(2)
      CLOSE(UNIT=10)
      CALL RECUIN(iseeds0(1),iseeds0(2))      ! For RANECU
      RETURN
C
  200 CONTINUE
      OPEN(UNIT=10,STATUS='NEW',FILE=INFILE)
      iseeds0(1) = 12345
      iseeds0(2) = 67890
      WRITE(10,300) iseeds0
  300 FORMAT(i15,',',i15)
      CLOSE(unit=10)
      CALL recuin(iseeds0(1),iseeds0(2))      ! For RANECU
      RETURN
C
      END
C
C ----------------------------------------------------------------------
C       Subroutine SAVESEED
C
C       PURPOSE:
C           Write out the current seeds to the file seeds_end.dat.
C
C       MODIFICATIONS:
C          10/16/92   PEU:  Now uses different file name than for
C                     starting seeds allowing replay with the same
C                     seeds.
C ----------------------------------------------------------------------
C
      SUBROUTINE SAVESEED
C
      CHARACTER*100 DAT_DIR
      CHARACTER*200 OUTFILE
      INTEGER ISEEDS0(2),ISEEDS(2)
      COMMON /RANCOM/ ISEEDS0
      COMMON /DATDIR_C/ DAT_DIR
      COMMON /DATDIR_I/ N_DAT_DIR
C
      OUTFILE = DAT_DIR(1:N_DAT_DIR)//'/seeds_end.dat'
      OPEN(UNIT=10,STATUS='UNKNOWN',FILE=OUTFILE)
      CALL RECUUT(iseeds(1),iseeds(2))
      WRITE(10,300) iseeds
  300 FORMAT(i15,',',i15)
      CLOSE(UNIT=10)
      RETURN
C
      END
