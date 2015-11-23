C
C -----------------------------------------------------------------------
C -----------------------------------------------------------------------
C
C       SUBROUTINE INDEX_ACCEPT_FCN
C
C       AUTHOR:  Paul E. Ulmer
C       DATE:    11-DEC-2004
C       PURPOSE: Setup the variable indices for the acceptance
C                function Ntuple.  Obviously, this implies that
C                certain indices in var.cmn/var_dat.cmn should
C                never be changed.
C
C -----------------------------------------------------------------------
C
      SUBROUTINE INDEX_ACCEPT_FCN(NVAR,IND_ACCEPT_FCN)
C
      IMPLICIT NONE
C
      INTEGER NVAR,IND_ACCEPT_FCN(NVAR),I,J
C
      IND_ACCEPT_FCN(1) = 164  ! Electron solid angle Jacobian
      IND_ACCEPT_FCN(2) = 165  ! Hadron   solid angle Jacobian
C
      J = 2
C
      DO I=61,66               ! Electron initial Transport vector
         J = J + 1
         IND_ACCEPT_FCN(J) = I
      ENDDO
C
      DO I=71,76               ! Hadron   initial Transport vector
         J = J + 1
         IND_ACCEPT_FCN(J) = I
      ENDDO
C
      DO I=101,106             ! Electron final Transport vector
         J = J + 1
         IND_ACCEPT_FCN(J) = I
      ENDDO
C
      DO I=111,116             ! Hadron   final Transport vector
         J = J + 1
         IND_ACCEPT_FCN(J) = I
      ENDDO
C
      RETURN
      END
