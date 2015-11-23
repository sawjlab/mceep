C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE TOP_1D
C
C       Writes a Topdrawer file of passed array (ARR) with error DARR.
C       XMIN and XMAX are the limits for the x-variable and
C       NBINS is the number of bins in the array (from 1 to NBINS).
C       Logical flag ERRWRITE = .true. will cause points with error bars
C       to be plotted.  Otherwise a histogram is plotted.
C------------------------------------------------------------------------------
C
      SUBROUTINE TOP_1D(XMIN,XMAX,NBINS,IOOB,TMP_FILE,LTIT,BTIT,
     #                  ARR,DARR,ERRWRITE)
      IMPLICIT NONE
C
      DOUBLE PRECISION ARR(500),DARR(500)
      DOUBLE PRECISION XMIN,XMAX,SUM_Y,SUM_XY,X,CENT,DELTAX
      INTEGER NBINS,IOOB,NCHAR,IL,I
      CHARACTER*300 OUTFILE,TMP_FILE
      CHARACTER*(*) LTIT,BTIT
      LOGICAL ERRWRITE
C
      CALL SQUEEZE(TMP_FILE,TMP_FILE,NCHAR)
      IL = 25
      OUTFILE = TMP_FILE(1:NCHAR)
C
      OPEN(UNIT=IL,FILE=OUTFILE,STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE(IL,'(''   set size 7.88 by 10.5 '')')
      WRITE(IL,'(''   set font duplex '')')
      WRITE(IL,'(''   set intensity 1 '')')
      WRITE(IL,'(''   title top '''' '',A,'' '''' '')')
     #                  TMP_FILE(1:NCHAR)
      WRITE(IL,'(''   title bottom '''' '',A,'' '''' '')') BTIT
      WRITE(IL,'(''   title left '''' '',A,'' '''' '')') LTIT
      IF(ERRWRITE)THEN
        WRITE(IL,'(''   set order x y dy '')')
      ELSE
        WRITE(IL,'(''   set order x y '')')
      ENDIF
C
      SUM_Y = 0.              !Initialize sum
      SUM_XY = 0.
C
      DELTAX = (XMAX-XMIN)/DFLOAT(NBINS)      !bin width
C
      DO I=1,NBINS
        X = (DFLOAT(I)-0.5D0)*DELTAX + XMIN
        IF(ERRWRITE)THEN
          WRITE(IL,10)X,ARR(I),DARR(I)
        ELSE
          WRITE(IL,11)X,ARR(I)
        ENDIF
   10   FORMAT(1X,F12.5,1X,E14.5,1X,E14.5)
   11   FORMAT(1X,F12.5,1X,E14.5)
        SUM_XY = SUM_XY + X*ARR(I)
        SUM_Y = SUM_Y + ARR(I)
      ENDDO
C
      IF(SUM_Y .NE. 0.) THEN
        CENT = SUM_XY/SUM_Y           !Centroid
      ELSE
        CENT = 0.
      ENDIF
C
      IF(ERRWRITE)THEN
        WRITE(IL,'(''   plot points '')')
      ELSE
        WRITE(IL,'(''   histogram solid '')')
      ENDIF
C
      WRITE(IL,
     # '(''   title 1.8 9.1 '''' Centroid: '',T40,E12.4,'' '''' '')')
     #         CENT
      WRITE(IL,
     # '(''   title 1.8 8.8 '''' Sum: '',T40,E12.4,'' '''' '')')
     #         SUM_Y
      WRITE(IL,
     # '(''   title 1.8 8.2 '''' Out of bounds: '',T40,I12,'' '''' '')')
     #         IOOB
C
      CLOSE(UNIT=IL)
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C       SUBROUTINE TOP_2D
C
C       Writes a 2-dimensional Topdrawer file of passed array (ARR2D).
C       XMIN (YMIN) and XMAX (YMAX) are the limits for the
C       x-variable (y-variable) and NBINSX and NBINSY are the number
C       of bins in each variable respectively (from 1 to NBINSX(Y)).
C------------------------------------------------------------------------------
C
      SUBROUTINE TOP_2D(XMIN,XMAX,YMIN,YMAX,NBINSX,NBINSY,
     #                  TMP_FILE,ARR2D)
      IMPLICIT NONE
C
      CHARACTER*300 OUTFILE,TMP_FILE
      DOUBLE PRECISION ARR2D(50,50),X(50),Y(50)
      DOUBLE PRECISION XMIN,XMAX,YMIN,YMAX,DELTAX,DELTAY
      INTEGER NBINSX,NBINSY,NCHAR,IL,I,J
C
      CALL SQUEEZE(TMP_FILE,TMP_FILE,NCHAR)
      IL = 25
      OUTFILE = TMP_FILE(1:NCHAR)
C
      DELTAX = (XMAX-XMIN)/DFLOAT(NBINSX)    !bin width for X
      DO I=1,NBINSX
        X(I) = (DFLOAT(I)-0.5D0)*DELTAX + XMIN
      ENDDO
C
      DELTAY = (YMAX-YMIN)/DFLOAT(NBINSY)    !bin width for Y
      DO J=1,NBINSY
        Y(J) = (DFLOAT(J)-0.5D0)*DELTAY + YMIN
      ENDDO
C
      OPEN(UNIT=IL,FILE=OUTFILE,STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE(IL,'(''   set size 7.88 by 10.5 '')')
      WRITE(IL,'(''   set font duplex '')')
      WRITE(IL,'(''   set intensity 1 '')')
      WRITE(IL,'(''   set outline off '')')
      WRITE(IL,'(''   title top '''' '',A,'' '''' '')')
     #                  TMP_FILE(1:NCHAR)
      WRITE(IL,'(''   set three '')')
      WRITE(IL,'(''   read mesh for x = '')')
      DO I=1,NBINSX
        WRITE(IL,20) X(I)
   20   FORMAT(T3,F10.3)
      ENDDO
      DO J=1,NBINSY
        WRITE(IL,30) Y(J),(ARR2D(I,J),I=1,NBINSX)
   30   FORMAT(/,T3,' for y = ',F12.5,2X,' z = ',100(/,T3,E14.5))
      ENDDO
      WRITE(IL,'(''   histogram '')')
      WRITE(IL,'(''   plot axes '')')
C
      CLOSE(UNIT=IL)
      RETURN
      END
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C       SUBROUTINE TOP_SCAT
C
C       Opens and sets up 2-dimensional scatter plots.  The titles
C       and limits are set up here.  The data is written in event by
C       event in the main routine.
C----------------------------------------------------------------------
C
      SUBROUTINE TOP_SCAT(LUNI,TMP_FILE,BTIT2D,LTIT2D,SQUEEZE_AXES,
     #          X2MIN,X2MAX,Y2MIN,Y2MAX)
      IMPLICIT NONE
C
      DOUBLE PRECISION X2MIN,X2MAX,Y2MIN,Y2MAX
      INTEGER LUNI,NCHAR_F,NCHAR_B,NCHAR_L
      CHARACTER*300 FILE2D,TMP_FILE
      CHARACTER*(*) BTIT2D,LTIT2D
      LOGICAL SQUEEZE_AXES                    !Squeeze axis titles?
C
      CALL SQUEEZE(TMP_FILE,TMP_FILE,NCHAR_F)
      IF(SQUEEZE_AXES) CALL SQUEEZE(BTIT2D,BTIT2D,NCHAR_B)
      IF(SQUEEZE_AXES) CALL SQUEEZE(LTIT2D,LTIT2D,NCHAR_L)
      FILE2D = TMP_FILE(1:NCHAR_F)
C
      OPEN(UNIT=LUNI,FILE=FILE2D,STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE(LUNI,'(''   set size 7.88 by 10.5 '')')
      WRITE(LUNI,'(''   set font duplex '')')
      WRITE(LUNI,'(''   set intensity 1 '')')
      WRITE(LUNI,
     #  '(''   title top '''' '',A,'' '''' '')') TMP_FILE(1:NCHAR_F)
      IF(SQUEEZE_AXES) THEN
        WRITE(LUNI,
     #    '(''   title bottom '''' '',A,'' '''' '')') BTIT2D(1:NCHAR_B)
        WRITE(LUNI,
     #    '(''   title left '''' '',A,'' '''' '')') LTIT2D(1:NCHAR_L)
      ELSE
        WRITE(LUNI,
     #    '(''   title bottom '''' '',A,'' '''' '')') BTIT2D
        WRITE(LUNI,
     #    '(''   title left '''' '',A,'' '''' '')') LTIT2D
      ENDIF
      WRITE(LUNI,
     #  '(''   set limits X from '',E12.4,'' to '',E12.4,'' '')')
     #     X2MIN,X2MAX
      WRITE(LUNI,
     #  '(''   set limits Y from '',E12.4,'' to '',E12.4,'' '')')
     #     Y2MIN,Y2MAX
      RETURN
      END
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C       SUBROUTINE NTU_SETUP
C
C       Opens and sets up N-Tuple files.
C
C       The particle ID tag and variable indices are written
C       here.  The data is written in event by event in the main
C       routine.
C
C       For IWEIGHT .GE. 0, a temporary file is opened since all
C       sigmas must be divided by the number of events into the
C       acceptance, SUM_ACC; SUM_ACC is not known until the end of the
C       MCEEP run.
C----------------------------------------------------------------------
C
      SUBROUTINE NTU_SETUP(LUNI,TMP_FILE,IWEIGHT,NTU_TRANS,
     #                     NVAR,IND_VAR)
      IMPLICIT NONE
C
      INTEGER LUNI,NCHAR_F,IWEIGHT,NVAR,IND_VAR(NVAR),I
      CHARACTER*300 FILENTU,FILENTU_TMP,TMP_FILE
      LOGICAL NTU_TRANS
C
      CALL SQUEEZE(TMP_FILE,TMP_FILE,NCHAR_F)
      FILENTU     = TMP_FILE(1:NCHAR_F)
      FILENTU_TMP = 'TEMP_'//TMP_FILE(1:NCHAR_F)
C
      OPEN(UNIT=LUNI,FILE=FILENTU,STATUS='UNKNOWN',FORM='FORMATTED')
      IF(IWEIGHT .GE. 0)
     #      OPEN(UNIT=LUNI+10,FILE=FILENTU_TMP,STATUS='UNKNOWN',
     #           FORM='FORMATTED')
C
      IF (NTU_TRANS) THEN
         WRITE(LUNI,10) ' TRANS '
         IF(IWEIGHT .GE. 0) WRITE(LUNI+10,10) ' TRANS '
10       FORMAT(A)
      ELSE
         WRITE(LUNI,10) ' KINEM '
         IF(IWEIGHT .GE. 0) WRITE(LUNI+10,10) ' KINEM '
      ENDIF
      WRITE(LUNI,20) IWEIGHT,NVAR
      IF(IWEIGHT .GE. 0) WRITE(LUNI+10,20) IWEIGHT,NVAR
20    FORMAT(1X,I5,1X,I5)
      WRITE(LUNI,30) (IND_VAR(I),I=1,NVAR)
      IF(IWEIGHT .GE. 0) WRITE(LUNI+10,30) (IND_VAR(I),I=1,NVAR)
30    FORMAT(1X,100(I3,1X))
C
      RETURN
      END
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C       SUBROUTINE NTM_SETUP
C
C       Opens and sets up Multiple Weight N-Tuple files.
C
C       The particle ID tags and variable indices are written
C       here.  The data is written in event by event in the main
C       routine.
C
C       For this type of N-tuple, it is assumed that all the weight
C       indices are >0.   A temporary file is opened since all
C       sigmas must be divided by the number of events into the
C       acceptance, SUM_ACC; SUM_ACC is not known until the end of the
C       MCEEP run.
C----------------------------------------------------------------------
C
      SUBROUTINE NTM_SETUP(LUNI,TMP_FILE,NWT,IND_WT,NVAR,IND_VAR)
      IMPLICIT NONE
C
      INTEGER LUNI,NCHAR_F,NWT,IND_WT(NWT),NVAR,IND_VAR(NVAR),I
      CHARACTER*300 FILENTU,FILENTU_TMP,TMP_FILE
C
      CALL SQUEEZE(TMP_FILE,TMP_FILE,NCHAR_F)
      FILENTU     = TMP_FILE(1:NCHAR_F)
      FILENTU_TMP = 'TEMP_'//TMP_FILE(1:NCHAR_F)
C
      OPEN(UNIT=LUNI,FILE=FILENTU,STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(UNIT=LUNI+10,FILE=FILENTU_TMP,STATUS='UNKNOWN',
     #           FORM='FORMATTED')
C
      WRITE(LUNI,20) NWT,NVAR
      WRITE(LUNI+10,20) NWT,NVAR
20    FORMAT(1X,I5,1X,I5)
      WRITE(LUNI,30)    (IND_WT(I),I=1,NWT)
      WRITE(LUNI,30)    (IND_VAR(I),I=1,NVAR)
      WRITE(LUNI+10,30) (IND_WT(I),I=1,NWT)
      WRITE(LUNI+10,30) (IND_VAR(I),I=1,NVAR)
30    FORMAT(1X,100(I3,1X))
C
      RETURN
      END

