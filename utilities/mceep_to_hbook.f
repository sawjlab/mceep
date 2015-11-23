C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C     SUBROUTINE: MCEEP_TO_HBOOK
C
C     AUTHOR: P.E. Ulmer
C     DATE:   16-JUL-1997
C     MODIFICATIONS:
C     
C     9-JUL-1997 (PEU)
C     Now generates a single hbook file from multiple
C     top_drawer files and MCEEP-ntuple files.  This program can
C     be used with the script, make_hbook.com, to scan the MCEEP
C     input file for 1-D and 2-D plots and Ntuples to generate a 
C     single hbook file.
C
C     PURPOSE:
C     
C     Produce HBOOK file from MCEEP output files.
C
C     ROUTINES CALLED:
C               READ_MCEEP_1D  - reads in 1-D MCEEP Topdrawer file.
C               READ_MCEEP_2D  - reads in 2-D MCEEP Topdrawer file.
C               MCEEP_TO_HBOOK_NTU  - processes MCEEP 'NTU' file.
C               MCEEP_TO_HBOOK_NTM  - processes MCEEP 'NTM' file.
C
C
C     Routines READ_MCEEP_1D and READ_MCEEP_2D call subroutines
C     MAKE_HBOOK_1D and MAKE_HBOOK_2D respectively which produce
C     the HBOOK histograms.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE MCEEP_TO_HBOOK
C
      IMPLICIT NONE
C     
      COMMON /PAWC/HMEMOR
C
      INCLUDE 'var.cmn'            !Defines variable names
      INCLUDE 'spectrometer.cmn'   !Defines variable names
C     
      INTEGER HMEMOR(4000000)
      CHARACTER*80 OUTFILE
      CHARACTER*1 PLOT_TYPE
      INTEGER ID,ISTAT,LRECL
C
      INCLUDE 'var_dat.cmn'
      INCLUDE 'spectrometer_dat.cmn'
C
      DATA LRECL /1024/                  !Record length
C
C------------------------------------------------------------------------------
C     Open HBOOK file.
C------------------------------------------------------------------------------
C     
      WRITE(6,101) '$Enter HBOOK file > '
 101  FORMAT(A)
      READ(5,102) OUTFILE
 102  FORMAT(A80)
      CALL HLIMIT(4000000)
      CALL HROPEN(7,'FILE',OUTFILE,'N',LRECL,ISTAT)
      ID = 0                                    !histogram ID
C
C------------------------------------------------------------------------------
C     Determine type of plot
C------------------------------------------------------------------------------
C
    1 WRITE(6,110)
 110  FORMAT(/' 1D, 2D, NTU, NTM  (1, 2, N, M or <RET> to Exit) > ')
      READ(5,112) PLOT_TYPE
  112 FORMAT(A1)
      IF(PLOT_TYPE .EQ. ' ') GOTO 999
      ID = ID + 1                               !Update histogram ID
C
      IF(PLOT_TYPE.EQ.'1') THEN
         CALL READ_MCEEP_1D(ID)
      ELSEIF(PLOT_TYPE.EQ.'2') THEN
         CALL READ_MCEEP_2D(ID)
      ELSEIF(PLOT_TYPE.EQ.'N' .OR. PLOT_TYPE.EQ.'n') THEN
         CALL MCEEP_TO_HBOOK_NTU(ID)
      ELSEIF(PLOT_TYPE.EQ.'M' .OR. PLOT_TYPE.EQ.'m') THEN
         CALL MCEEP_TO_HBOOK_NTM(ID)
      ENDIF
C
      GOTO 1
 999  CALL HREND('FILE')
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       SUBROUTINE READ_MCEEP_1D
C
C       PURPOSE:
C               Produces HBOOK file from a 1-D MCEEP-generated
C               TOPDRAWER file.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE READ_MCEEP_1D(ID)
C
      IMPLICIT NONE
C
      CHARACTER*80 INFILE
      CHARACTER*12 DUMMY
      REAL X(500),Y(500)
      REAL XMIN,XMAX,DELTAX
      INTEGER ID,NBINS,NCHAR,I
C
      WRITE(6,10) '$Enter input file > '
   10 FORMAT(A)
      READ(5,15) INFILE
   15 FORMAT(A80)
      CALL SQUEEZE(INFILE,INFILE,NCHAR)
C
      OPEN(UNIT=1,FILE=INFILE,STATUS='OLD',FORM='FORMATTED')
      DO I=1,7
        READ(1,20) DUMMY
   20   FORMAT(A12)
      ENDDO
C
      DO I=1,500
        READ(1,*,ERR=199) X(I),Y(I)
        NBINS = I
      ENDDO
  199 CLOSE(UNIT=1)
C
      DELTAX = X(2) - X(1)
      XMIN = X(1)     - DELTAX/2.
      XMAX = X(NBINS) + DELTAX/2.
C
      CALL MAKE_HBOOK_1D(ID,XMIN,XMAX,NBINS,Y,INFILE(1:NCHAR))
C
      RETURN
      END
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C     SUBROUTINE MAKE_HBOOK_1D
C
C       Produces a 1D HBOOK file from a passed array.
C ------------------------------------------------------------------
C
      SUBROUTINE MAKE_HBOOK_1D(ID,XMIN,XMAX,NBINS,Y,INFILE)
C
      IMPLICIT NONE
C
      COMMON /PAWC/HMEMOR
      REAL XMIN,XMAX,XVAL,YVAL,Y(500)
      INTEGER HMEMOR(4000000)
      INTEGER NBINS,N,I,ID
      CHARACTER*(*) INFILE
C
C ------------------------------------------------------------------
C       Set up histogram:  Call HBOOK1
C ------------------------------------------------------------------
C
      CALL HBOOK1(ID,INFILE,NBINS,XMIN,XMAX,0.)
C
C ------------------------------------------------------------------
C       Fill histogram:  Call HFILL
C ------------------------------------------------------------------
C
      DO I=1,NBINS
        XVAL = (DFLOAT(I)-0.5)*(XMAX-XMIN)/DFLOAT(NBINS)+XMIN
        YVAL = Y(I)
        CALL HFILL(ID,XVAL,0.,YVAL)
      ENDDO
C
C ------------------------------------------------------------------
C     Write into file.
C ------------------------------------------------------------------
C     
      CALL HROUT(ID,N,' ')
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C
C       SUBROUTINE READ_MCEEP_2D
C
C       PURPOSE:
C               Produces HBOOK file from a 2-D MCEEP-generated
C               TOPDRAWER file.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE READ_MCEEP_2D(ID)
C
      IMPLICIT NONE
C
      CHARACTER*80 INFILE
      CHARACTER*12 DUMMY
      REAL X(500),Y(500)
      REAL XMIN,XMAX,YMIN,YMAX
      REAL ARR2D(500,500)
      INTEGER ID,NBINSX,NBINSY,NCHAR,I,J
C
      WRITE(6,10) '$Enter input file >'
   10 FORMAT(A)
      READ(5,15) INFILE
   15 FORMAT(A80)
      CALL SQUEEZE(INFILE,INFILE,NCHAR)
C
      OPEN(UNIT=1,FILE=INFILE,STATUS='OLD',FORM='FORMATTED')
      DO I=1,7
        READ(1,20) DUMMY
   20   FORMAT(A12)
      ENDDO
C
      DO I=1,500
        READ(1,25,ERR=199) X(I)
   25   FORMAT(T3,F10.3)
        NBINSX = I
      ENDDO
  199 NBINSX = NBINSX-1
      BACKSPACE 1
C
      DO J=1,500
        READ(1,31,ERR=299) Y(J)
   31   FORMAT(T12,F12.5)
        NBINSY = J
        DO I=1,NBINSX
          READ(1,35) ARR2D(I,J)
   35     FORMAT(T3,E14.5)
        ENDDO
        READ(1,20) DUMMY
        READ(1,20) DUMMY
      ENDDO
  299 CLOSE(UNIT=1)
C
      XMIN = X(1)
      XMAX = X(NBINSX)
      YMIN = Y(1)
      YMAX = Y(NBINSY)
C
      CALL MAKE_HBOOK_2D(ID,XMIN,XMAX,YMIN,YMAX,NBINSX,NBINSY,ARR2D,
     #     INFILE(1:NCHAR))
C
      RETURN
      END
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C     SUBROUTINE MAKE_HBOOK_2D
C
C       Produces a 2D HBOOK file from a passed array.
C ------------------------------------------------------------------
C
      SUBROUTINE MAKE_HBOOK_2D(ID,XMIN,XMAX,YMIN,YMAX,NX,NY,ARR2D,
     #     INFILE)
C
      IMPLICIT NONE
C
      COMMON /PAWC/HMEMOR
      REAL ARR2D(500,500)
      REAL XMIN,XMAX,YMIN,YMAX,XVAL,YVAL
      INTEGER HMEMOR(4000000)
      INTEGER NX,NY,N,I,J,ID
      CHARACTER*(*) INFILE
C
C ------------------------------------------------------------------
C       Set up histogram:  Call HBOOK2
C ------------------------------------------------------------------
C
      CALL HBOOK2(ID,INFILE,NX+1,XMIN,XMAX,NY+1,YMIN,YMAX,0.)
C
C ------------------------------------------------------------------
C       Fill histogram:  Call HFILL
C ------------------------------------------------------------------
C
      DO I=1,NX
        DO J=1,NY
          XVAL = DFLOAT(I)*(XMAX-XMIN)/DFLOAT(NX)+XMIN
          YVAL = DFLOAT(J)*(YMAX-YMIN)/DFLOAT(NY)+YMIN
          CALL HFILL(ID,XVAL,YVAL,ARR2D(I,J))
        ENDDO
      ENDDO
C
C ------------------------------------------------------------------
C     Write to file.
C ------------------------------------------------------------------
C
      CALL HROUT(ID,N,' ')
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C      MCEEP_TO_HBOOK_NTU
C
C      Reads in a MCEEP N-Tuple ('NTU') file and converts it
C      to HBOOK format.
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
C
      SUBROUTINE MCEEP_TO_HBOOK_NTU(ID)
C
      IMPLICIT NONE
C
      COMMON /PAWC/ HMEMOR
C
      CHARACTER*80     INFILE
      CHARACTER*5      NTU_TYPE
      CHARACTER*8      NTU_VAR_NAME(102)
      REAL             NTU_VAR(102)
      INTEGER          HMEMOR(4000000)
      INTEGER          IWEIGHT,N_WTS,N_VAR,IND_VAR(100)
      INTEGER          LUN_INP
      INTEGER          N_ELMNTS,ID
      INTEGER          NPRIME
      INTEGER          I,J,N,NCHAR
C
      INCLUDE 'var.cmn'            !Defines variable names
      INCLUDE 'spectrometer.cmn'   !Defines variable names
C
      DATA LUN_INP /1/                   !Logical units for infile
      DATA NPRIME /5000/                 !primary alloc. (# words)
C
C ---------------------------------------------------------------------
C     Read in the MCEEP N-Tuple file header.
C ---------------------------------------------------------------------
C
      WRITE(6,10) '$Enter MCEEP NTU N-Tuple filename >'
10    FORMAT(A)
      READ(5,20) INFILE
20    FORMAT(A80)
      CALL SQUEEZE(INFILE,INFILE,NCHAR)
      OPEN(UNIT=LUN_INP,FILE=INFILE,STATUS='OLD',FORM='FORMATTED')
      READ(LUN_INP,30) NTU_TYPE     !Kinematics or Transport?
30    FORMAT(1X,A5)
      READ(LUN_INP,40) IWEIGHT,N_VAR !PID, # indep. variables in N-Tuple
40    FORMAT(1X,I5,1X,I5)
      READ(LUN_INP,50) (IND_VAR(J),J=1,N_VAR) !Indices of N-Tuple
                                              ! independent variables
50    FORMAT(1X,100(I3,1X))
C
C ---------------------------------------------------------------------
C     Determine the number of weights to read in.
C
C        If IWEIGHT .LT. 0  then we need two weights to define
C            bin-averaged cross sections, polarizations, etc.
C            For example, the avg. cross section is defined by
C            dividing the yield for a given bin by the
C            phase space volume of that bin.  Therefore, we need
C            both the yield and the phase space weights.
C
C        If IWEIGHT .GE. 0  then we only need one weight which
C            defines the yield or partial yield for a given
C            MCEEP "event".  (Except when IWEIGHT=9999; see below.)
C
C        If IWEIGHT .EQ. 9999 then we are dealing with the 
C            acceptance Ntuple.  In this case, there are no weights
C            though the Jacobians are included as independent
C            variables.
C ---------------------------------------------------------------------
C
      N_WTS = 1
      IF(IWEIGHT .LT. 0) N_WTS = 2
      IF(IWEIGHT .EQ. 9999) N_WTS = 0
      N_ELMNTS = N_VAR + N_WTS        !Total number of elements in
                                      ! HBOOK array.
C
C ---------------------------------------------------------------------
C     Set up labels for Kinematics or Transport N-Tuple
C     independent variables.
C ---------------------------------------------------------------------
C
      IF(NTU_TYPE .EQ. 'TRANS') THEN
         DO I=1,N_VAR
            NTU_VAR_NAME(I+N_WTS) = TR_LABEL(IND_VAR(I))//'   '
         ENDDO
      ELSEIF(NTU_TYPE .EQ. 'KINEM') THEN
         DO I=1,N_VAR
            NTU_VAR_NAME(I+N_WTS) = VAR_NAME(IND_VAR(I))
         ENDDO
      ENDIF
C
C ---------------------------------------------------------------------
C     Set up labels for Kinematics or Transport N-Tuple
C     dependent variables (i.e. weights).  These labels as well as
C     the values are packed into the appropriate arrays starting
C     with index N_VAR+1.
C ---------------------------------------------------------------------
C
      IF(N_WTS .NE. 0) THEN
         NTU_VAR_NAME(1) = 'Numer_wt'
         IF(IWEIGHT .LT. 0) THEN
            NTU_VAR_NAME(2) = 'Denom_wt'
         ENDIF
      ENDIF
C
CXXX  NTU_VAR_NAME(1) = REACT_NAME(IWEIGHT)
CXXX  IF(IWEIGHT .LT. 0) THEN
CXXX     CALL GET_PID_DEN(IWEIGHT,IWEIGHT_DEN)   
CXXX     NTU_VAR_NAME(2) = REACT_NAME(IWEIGHT_DEN)
CXXX  ENDIF
C
C ---------------------------------------------------------------------
C     Book the N-Tuple.
C ---------------------------------------------------------------------
C
      CALL HBOOKN(ID,INFILE,N_ELMNTS,'FILE',NPRIME,NTU_VAR_NAME)
C
C ---------------------------------------------------------------------
C     Fill the N-tuple.
C ---------------------------------------------------------------------
C
200   IF(N_WTS .NE. 0) THEN
         READ(LUN_INP,210,END=999) (NTU_VAR(I),I=1,N_WTS)
      ENDIF
210   FORMAT(100(1X,1PE12.5))
      READ(LUN_INP,210,END=999) (NTU_VAR(I+N_WTS),I=1,N_VAR)
      CALL HFN(ID,NTU_VAR)
      GOTO 200
C
999   CALL HROUT(ID,N,' ')
      CLOSE(UNIT=LUN_INP)
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C      MCEEP_TO_HBOOK_NTM
C
C      Reads in a MCEEP N-Tuple ('NTM') file and converts it
C      to HBOOK format.
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
C
      SUBROUTINE MCEEP_TO_HBOOK_NTM(ID)
C
      IMPLICIT NONE
C
      COMMON /PAWC/ HMEMOR
C
      CHARACTER*80     INFILE
      CHARACTER*8      NTU_VAR_NAME(200)
      REAL             NTU_VAR(200)
      INTEGER          HMEMOR(4000000)
      INTEGER          N_WT,N_VAR,IND_WT(100),IND_VAR(100)
      INTEGER          LUN_INP
      INTEGER          N_ELMNTS,ID
      INTEGER          NPRIME
      INTEGER          I,J,N,NCHAR
C
      INCLUDE 'var.cmn'            !Defines variable names
C
      DATA LUN_INP /1/                   !Logical units for infile
      DATA NPRIME /5000/                 !primary alloc. (# words)
C
C ---------------------------------------------------------------------
C     Read in the MCEEP N-Tuple file header.
C ---------------------------------------------------------------------
C
      WRITE(6,10) '$Enter MCEEP NTM N-Tuple filename >'
10    FORMAT(A)
      READ(5,20) INFILE
20    FORMAT(A80)
      CALL SQUEEZE(INFILE,INFILE,NCHAR)
      OPEN(UNIT=LUN_INP,FILE=INFILE,STATUS='OLD',FORM='FORMATTED')
      READ(LUN_INP,40) N_WT,N_VAR !# weights, # indep. variables
40    FORMAT(1X,I5,1X,I5)
      READ(LUN_INP,50) (IND_WT(J),J=1,N_WT)   !Indices of N-Tuple
                                              ! weights
      READ(LUN_INP,50) (IND_VAR(J),J=1,N_VAR) !Indices of N-Tuple
                                              ! independent variables
50    FORMAT(1X,100(I3,1X))
C
C ---------------------------------------------------------------------
C     Set up labels for Kinematics or Transport N-Tuple
C     dependent variables (i.e. weights).  These labels as well as
C     the values are packed into the appropriate arrays starting
C     with index N_VAR+1.
C ---------------------------------------------------------------------
C
      N_ELMNTS = N_WT + N_VAR
      DO I=1,N_WT
         NTU_VAR_NAME(I) = REACT_NAME(IND_WT(I))
      ENDDO
      DO I=1,N_VAR
         NTU_VAR_NAME(I+N_WT) = VAR_NAME(IND_VAR(I))
      ENDDO
C
C ---------------------------------------------------------------------
C     Book the N-Tuple.
C ---------------------------------------------------------------------
C
      CALL HBOOKN(ID,INFILE,N_ELMNTS,'FILE',NPRIME,NTU_VAR_NAME)
C
C ---------------------------------------------------------------------
C     Fill the N-tuple.
C ---------------------------------------------------------------------
C
200   READ(LUN_INP,210,END=999) (NTU_VAR(I),I=1,N_WT)
210   FORMAT(100(1X,1PE12.5))
      READ(LUN_INP,210) (NTU_VAR(I+N_WT),I=1,N_VAR)
      CALL HFN(ID,NTU_VAR)
      GOTO 200
C
999   CALL HROUT(ID,N,' ')
      CLOSE(UNIT=LUN_INP)
C
      RETURN
      END
