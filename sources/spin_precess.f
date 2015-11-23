C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C     SUBROUTINE SPIN_PRECESS
C
C     AUTHOR:    M. NOZAR
C     DATE:      10-JUL-1991
C     PURPOSE:
C                Calling routine for spin precession calculation.
C                Calls appropriate routines which return 3X3
C                spin rotation matrices.
C
C                For each subsequent call to this routine (each
C                occurrence of a "POL" element in the input file) the
C                next spin precession calculation subroutine is called.
C                Up to 10 different spin subroutines are currently
C                allowed which accomodates a spectrometer consisting
C                of up to 10 different magnetic elements.  At present,
C                a simple 1st-order calculation is employed
C                (SPIN_1ST_ORDER) where the precession angle for a
C                given element is assumed to be proportional to the
C                net bend angle for that element.
C
C                The user may incorporate more sophistocated
C                prescriptions by commenting out the calls to
C                SPIN_1ST_ORDER and enabling the commented lines
C                (CUSER lines) in this routine and by modifying
C                ../sources/spin_dat.for
C                (containing subroutines SPIN_1 to SPIN_10 which
C                currently merely return 3X3 identity matrices).
C
C
C                NOTE1:
C                      Any user routine replacing SPIN_1ST_ORDER
C                      will need to include SPECTROMETER.CMN
C                      so that both TVEC (the current Transport
C                      vector) and PF_P (the spectrometer central
C                      momentum) can be accessed.
C
C                NOTE2:
C                      The matrices returned by SPIN_1ST_ORDER
C                      are functions of the current Transport
C                      vector.  SPIN_1ST_ORDER implicitly assumes
C                      that this vector corresponds to the ray at
C                      the ENTRANCE of the element.  Thus, any "POL"
C                      elements in the input file should precede
C                      the corresponding "MAT" element.  (To compute
C                      the net bend angle for a given magnetic element,
C                      SPIN_1ST_ORDER compares the current (i.e.
C                      entrance) Transport vector with the same
C                      vector after multiplication by the appropriate
C                      Tranport matrix (i.e. exit vector).  This
C                      multiplication is done explicitly in 
C                      SPIN_1ST_ORDER.)
C
C                      If the user wishes to replace
C                      the first order routines with routines
C                      wherein the matrix for a given spectrometer
C                      element is based on the Transport coordinates
C                      at the exit of the element, then the "POL"
C                      element in the input file should follow any
C                      Transport 6X6 matrix for that element and vice
C                      versa.
C
C ----------------------------------------------------------------------
C
      SUBROUTINE SPIN_PRECESS (N,BEND_0,PMAT_EL)
      IMPLICIT NONE
      DOUBLE PRECISION BEND_0,PMAT_EL(3,3)
      INTEGER N
C      
C ----------------------------------------------------------------------
C     Get the Polarization matrix for each element
C ----------------------------------------------------------------------
C
      GOTO (1,2,3,4,5,6,7,8,9,10)N
      WRITE(6,100)
100   FORMAT (' ILLEGAL ENTRY FOR SPIN TRANSPORT')  
      STOP
C
1     CALL SPIN_1ST_ORDER(BEND_0,PMAT_EL)
CUSER CALL SPIN_1(BEND_0,PMAT_EL)
      RETURN
C
2     CALL SPIN_1ST_ORDER(BEND_0,PMAT_EL)
CUSER CALL SPIN_2(BEND_0,PMAT_EL)
      RETURN
C
3     CALL SPIN_1ST_ORDER(BEND_0,PMAT_EL)
CUSER CALL SPIN_3(BEND_0,PMAT_EL)
      RETURN
C
4     CALL SPIN_1ST_ORDER(BEND_0,PMAT_EL)
CUSER CALL SPIN_4(BEND_0,PMAT_EL)
      RETURN
C
5     CALL SPIN_1ST_ORDER(BEND_0,PMAT_EL)
CUSER CALL SPIN_5(BEND_0,PMAT_EL)
      RETURN
C
6     CALL SPIN_1ST_ORDER(BEND_0,PMAT_EL)
CUSER CALL SPIN_6(BEND_0,PMAT_EL)
      RETURN
C
7     CALL SPIN_1ST_ORDER(BEND_0,PMAT_EL)
CUSER CALL SPIN_7(BEND_0,PMAT_EL)
      RETURN
C
8     CALL SPIN_1ST_ORDER(BEND_0,PMAT_EL)
CUSER CALL SPIN_8(BEND_0,PMAT_EL)
      RETURN
C
9     CALL SPIN_1ST_ORDER(BEND_0,PMAT_EL) 
CUSER CALL SPIN_9(BEND_0,PMAT_EL)
      RETURN
C
10    CALL SPIN_1ST_ORDER(BEND_0,PMAT_EL) 
CUSER CALL SPIN_10(BEND_0,PMAT_EL)
      RETURN
C
      END
